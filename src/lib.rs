//! Differential genotype-missingness test, a clean-room Rust port of PLINK 1.9
//! `--test-missing`.
//!
//! For each variant we form the 2×2 table of (missing vs non-missing genotype) ×
//! (case vs control) and report the two-sided Fisher's exact p-value with the
//! per-group missing-call fractions, in PLINK's `.missing` column order:
//! `CHR SNP F_MISS_A F_MISS_U P`. Only nondegenerate variants — those with at
//! least one missing call across cases and controls — are written, matching
//! PLINK.

mod bed;
mod fisher;
mod fmt;
mod pheno;
mod testmissing;

use anyhow::Result;
use bed::Fileset;
use rayon::prelude::*;
use std::collections::HashMap;
use std::io::Write;
use std::path::Path;
use testmissing::{ChromClass, StatusMasks, VariantResult, classify_chrom};

pub use pheno::{MissingCode, PhenoOptions};

/// Variants processed per streamed `.bed` block — bounds resident genotype
/// memory to `BLOCK × bytes_per_variant`, while staying large enough to amortise
/// rayon's per-block fan-out.
const BLOCK: usize = 4096;

/// Open a fileset for a `--test-missing` scan.
pub fn open(prefix: &Path) -> Result<Fileset> {
    Fileset::open(prefix)
}

/// PLINK's numeric chromosome code (X→23, Y→24, XY→25, MT→26; others pass
/// through), used for the CHR column.
fn chrom_code_into(out: &mut String, chrom: &str) {
    out.push_str(match chrom {
        "X" | "x" => "23",
        "Y" | "y" => "24",
        "XY" | "xy" => "25",
        "MT" | "mt" | "M" | "m" => "26",
        other => other,
    });
}

/// Stream the test over the fileset and write the `.missing` report. Genotypes
/// are read in blocks; each block's tables are computed in parallel and its
/// nondegenerate rows written before the next block is read, so peak memory
/// stays near one block regardless of fileset size.
pub fn run(fs: &mut Fileset, opts: &PhenoOptions, out: &mut impl Write) -> Result<()> {
    let status = pheno::load_status(&fs.samples, opts)?;
    let sex: Vec<u8> = fs.samples.iter().map(|s| s.sex).collect();
    let masks = StatusMasks::build(&status, &sex, fs.bytes_per_variant);
    let bpv = fs.bytes_per_variant;
    let classes: Vec<ChromClass> = fs
        .variants
        .iter()
        .map(|v| classify_chrom(&v.chrom))
        .collect();

    let snpw = snp_width(&fs.variants);
    write_header(out, snpw)?;

    // F_MISS_A/U are k/den for one of a few fixed denominators, so their %.4g
    // forms are tabulated once and looked up rather than re-formatted per row.
    let mut fmiss = FmissCache::default();
    // A variant's P (and its %.4g text) is fixed by its 2×2 counts; dense data
    // repeats only a handful of count combinations, so memoise on them.
    let mut pcache: HashMap<(u64, u64, u64, u64), String> = HashMap::new();
    let lf = fisher::LnFactTable::new(fs.n_samples());

    // The single-thread path is the perfgate baseline and avoids rayon's
    // per-block fan-out tax; multiple threads decode a block in parallel.
    let parallel = rayon::current_num_threads() > 1;
    let mut block = Vec::with_capacity(BLOCK * bpv);
    let mut results: Vec<VariantResult> = Vec::with_capacity(BLOCK);
    let mut scratch = RowScratch::default();
    let mut first = 0usize;
    loop {
        let count = fs.read_block(BLOCK, &mut block)?;
        if count == 0 {
            break;
        }
        if parallel {
            results.clear();
            (0..count)
                .into_par_iter()
                .map(|i| masks.result(&block[i * bpv..i * bpv + bpv], classes[first + i]))
                .collect_into_vec(&mut results);
            for (i, r) in results.iter().enumerate() {
                if r.nondegenerate() {
                    write_variant(
                        out,
                        &mut scratch,
                        &mut fmiss,
                        &mut pcache,
                        &lf,
                        &fs.variants[first + i],
                        r,
                        snpw,
                    )?;
                }
            }
        } else {
            for i in 0..count {
                let r = masks.result(&block[i * bpv..i * bpv + bpv], classes[first + i]);
                if r.nondegenerate() {
                    write_variant(
                        out,
                        &mut scratch,
                        &mut fmiss,
                        &mut pcache,
                        &lf,
                        &fs.variants[first + i],
                        &r,
                        snpw,
                    )?;
                }
            }
        }
        first += count;
    }
    Ok(())
}

/// PLINK's `.missing` SNP column width: 4 for ids up to 4 chars, else the
/// longest id plus 2.
fn snp_width(variants: &[bed::Variant]) -> usize {
    let max_id = variants.iter().map(|v| v.id.len()).max().unwrap_or(0);
    if max_id <= 4 { 4 } else { max_id + 2 }
}

fn write_header(out: &mut impl Write, snpw: usize) -> std::io::Result<()> {
    writeln!(
        out,
        "{:>4} {:>snpw$} {:>12} {:>12} {:>12} ",
        "CHR", "SNP", "F_MISS_A", "F_MISS_U", "P"
    )
}

/// Reusable per-row scratch so the hot writer makes no allocations.
#[derive(Default)]
struct RowScratch {
    line: String,
    chrom: String,
}

/// Memoised `%.4g(num/den)` for the handful of denominators a run uses (the
/// case/control counts), so F_MISS_A/U are looked up instead of re-formatted.
#[derive(Default)]
struct FmissCache {
    tables: Vec<(u64, Vec<String>)>,
}

impl FmissCache {
    fn get(&mut self, num: u64, den: u64) -> &str {
        let idx = match self.tables.iter().position(|(d, _)| *d == den) {
            Some(i) => i,
            None => {
                let mut scratch = String::new();
                let table: Vec<String> = (0..=den)
                    .map(|k| {
                        let mut s = String::new();
                        #[allow(clippy::cast_precision_loss)]
                        fmt::g_into(&mut s, &mut scratch, k as f64 / den as f64);
                        s
                    })
                    .collect();
                self.tables.push((den, table));
                self.tables.len() - 1
            }
        };
        &self.tables[idx].1[num as usize]
    }
}

#[allow(clippy::too_many_arguments)]
fn write_variant(
    out: &mut impl Write,
    s: &mut RowScratch,
    fmiss: &mut FmissCache,
    pcache: &mut HashMap<(u64, u64, u64, u64), String>,
    lf: &fisher::LnFactTable,
    var: &bed::Variant,
    r: &VariantResult,
    snpw: usize,
) -> std::io::Result<()> {
    s.line.clear();
    s.chrom.clear();
    chrom_code_into(&mut s.chrom, &var.chrom);
    pad_right(&mut s.line, &s.chrom, 4);
    s.line.push(' ');
    pad_right(&mut s.line, &var.id, snpw);

    s.line.push(' ');
    let fa = fmiss.get(r.miss_case, r.n_case);
    pad_right(&mut s.line, fa, 12);
    s.line.push(' ');
    let fu = fmiss.get(r.miss_ctrl, r.n_ctrl);
    pad_right(&mut s.line, fu, 12);

    let key = (r.miss_case, r.n_case, r.miss_ctrl, r.n_ctrl);
    let p_str = pcache.entry(key).or_insert_with(|| {
        let mut s = String::new();
        fmt::g_into(&mut s, &mut String::new(), r.p(lf));
        s
    });
    s.line.push(' ');
    pad_right(&mut s.line, p_str, 12);

    s.line.push('\n');
    out.write_all(s.line.as_bytes())
}

/// Append `cell` to `out` right-justified in a field of `width` (spaces left).
fn pad_right(out: &mut String, cell: &str, width: usize) {
    for _ in 0..width.saturating_sub(cell.len()) {
        out.push(' ');
    }
    out.push_str(cell);
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::path::PathBuf;

    fn pack(codes: &[u8]) -> Vec<u8> {
        let mut bytes = vec![0u8; codes.len().div_ceil(4)];
        for (s, &c) in codes.iter().enumerate() {
            bytes[s / 4] |= (c & 0b11) << ((s % 4) * 2);
        }
        bytes
    }

    /// `rows` are raw per-sample 2-bit codes (one inner vec per variant).
    fn write_fileset(dir: &Path, rows: &[Vec<u8>], aff: &[u8]) -> PathBuf {
        let prefix = dir.join("toy");
        let mut bim = File::create(prefix.with_extension("bim")).unwrap();
        for (v, _) in rows.iter().enumerate() {
            writeln!(bim, "1\trs{v}\t0\t{}\tA\tG", 100 + v).unwrap();
        }
        let mut fam = File::create(prefix.with_extension("fam")).unwrap();
        for (s, a) in aff.iter().enumerate() {
            writeln!(fam, "F{s}\tI{s}\t0\t0\t1\t{a}").unwrap();
        }
        let mut bed = File::create(prefix.with_extension("bed")).unwrap();
        bed.write_all(&[0x6c, 0x1b, 0x01]).unwrap();
        for r in rows {
            bed.write_all(&pack(r)).unwrap();
        }
        prefix
    }

    fn scan(prefix: &Path) -> Vec<Vec<String>> {
        let mut fs = open(prefix).unwrap();
        let opts = PhenoOptions {
            pheno_file: None,
            missing: MissingCode::default(),
            allow_no_sex: true,
            case_control_01: false,
        };
        let mut out = Vec::new();
        run(&mut fs, &opts, &mut out).unwrap();
        String::from_utf8(out)
            .unwrap()
            .lines()
            .map(|l| l.split_whitespace().map(str::to_string).collect())
            .collect()
    }

    #[test]
    fn degenerate_variants_are_skipped() {
        let dir = tempfile::tempdir().unwrap();
        // 4 samples: case, ctrl, case, ctrl.
        let aff = [2u8, 1, 2, 1];
        // rs0 fully called (skipped), rs1 has a missing case.
        let rows = vec![vec![0b00, 0b10, 0b11, 0b00], vec![0b01, 0b00, 0b00, 0b00]];
        let prefix = write_fileset(dir.path(), &rows, &aff);
        let lines = scan(&prefix);
        // header + only rs1.
        assert_eq!(lines.len(), 2);
        assert_eq!(lines[0], vec!["CHR", "SNP", "F_MISS_A", "F_MISS_U", "P"]);
        assert_eq!(lines[1][1], "rs1");
    }

    #[test]
    fn reports_fractions_and_p() {
        let dir = tempfile::tempdir().unwrap();
        // 8 samples: even = case, odd = control.
        let aff = [2u8, 1, 2, 1, 2, 1, 2, 1];
        // cases 0,2 missing (2/4); control 3 missing (1/4).
        let rows = vec![vec![0b01, 0b00, 0b01, 0b01, 0b11, 0b10, 0b00, 0b11]];
        let prefix = write_fileset(dir.path(), &rows, &aff);
        let lines = scan(&prefix);
        assert_eq!(lines.len(), 2);
        let r = &lines[1];
        assert_eq!(r[2], "0.5"); // F_MISS_A = 2/4
        assert_eq!(r[3], "0.25"); // F_MISS_U = 1/4
        let p: f64 = r[4].parse().unwrap();
        assert!((0.0..=1.0).contains(&p));
    }

    #[test]
    fn header_has_trailing_space_data_does_not() {
        let dir = tempfile::tempdir().unwrap();
        let aff = [2u8, 1, 2, 1];
        let rows = vec![vec![0b01, 0b00, 0b00, 0b00]];
        let prefix = write_fileset(dir.path(), &rows, &aff);
        let mut fs = open(&prefix).unwrap();
        let opts = PhenoOptions {
            pheno_file: None,
            missing: MissingCode::default(),
            allow_no_sex: true,
            case_control_01: false,
        };
        let mut out = Vec::new();
        run(&mut fs, &opts, &mut out).unwrap();
        let text = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert!(lines[0].ends_with(" P "), "header tail: {:?}", lines[0]);
        assert!(!lines[1].ends_with(' '), "data tail: {:?}", lines[1]);
    }

    #[test]
    fn chrom_codes() {
        let code = |c: &str| {
            let mut s = String::new();
            chrom_code_into(&mut s, c);
            s
        };
        assert_eq!(code("X"), "23");
        assert_eq!(code("MT"), "26");
        assert_eq!(code("7"), "7");
    }

    #[test]
    fn snp_width_rule() {
        let v = |id: &str| bed::Variant {
            chrom: "1".into(),
            id: id.into(),
        };
        assert_eq!(snp_width(&[v("r")]), 4);
        assert_eq!(snp_width(&[v("rsab")]), 4); // len 4
        assert_eq!(snp_width(&[v("rsabc")]), 7); // len 5 → +2
        assert_eq!(snp_width(&[v("rsLONGNAME123")]), 15); // len 13 → +2
    }
}
