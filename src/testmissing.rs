//! Differential genotype-missingness test (`plink --test-missing`).
//!
//! For each variant we build the 2×2 table
//!
//! ```text
//!            missing   non-missing
//!   cases       a           b
//!   controls    c           d
//! ```
//!
//! and report the two-sided Fisher's exact p-value, alongside the missing-call
//! fractions in cases (`F_MISS_A = a/(a+b)`) and controls (`F_MISS_U = c/(c+d)`).
//! Only nondegenerate variants — those with at least one missing call across
//! cases and controls — are written, matching PLINK.
//!
//! Missing genotypes are the packed 2-bit code `0b01` (lo-bit set, hi-bit clear).
//! PLINK additionally treats *het-haploid* calls as missing: a heterozygous
//! (`0b10`) male call on chromosome X, and any het male call on Y; on Y females
//! are dropped entirely. Per variant we mask the row's missing (and, on the
//! sex chromosomes, het) lanes against precomputed per-group lane masks and
//! popcount, so each variant costs a popcount sweep rather than a per-sample
//! branch.

use crate::fisher;
use crate::pheno::Status;

/// Chromosome class determining haploid handling.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum ChromClass {
    /// Autosomes, MT, and the XY pseudo-autosomal region: het is a valid call.
    Diploid,
    /// X: heterozygous male calls are treated as missing.
    X,
    /// Y: females are excluded; heterozygous male calls are treated as missing.
    Y,
}

pub fn classify_chrom(chrom: &str) -> ChromClass {
    match chrom {
        "23" | "X" | "x" => ChromClass::X,
        "24" | "Y" | "y" => ChromClass::Y,
        _ => ChromClass::Diploid,
    }
}

/// One reported variant's table summary.
pub struct VariantResult {
    pub miss_case: u64,
    pub n_case: u64,
    pub miss_ctrl: u64,
    pub n_ctrl: u64,
}

impl VariantResult {
    pub fn p(&self, lf: &fisher::LnFactTable) -> f64 {
        let nonmiss_case = self.n_case - self.miss_case;
        let nonmiss_ctrl = self.n_ctrl - self.miss_ctrl;
        fisher::two_sided(
            lf,
            self.miss_case,
            nonmiss_case,
            self.miss_ctrl,
            nonmiss_ctrl,
        )
    }
    /// PLINK writes only nondegenerate variants: those with at least one missing
    /// call in cases or controls.
    pub fn nondegenerate(&self) -> bool {
        self.miss_case + self.miss_ctrl > 0
    }

    #[cfg(test)]
    fn f_miss_a(&self) -> f64 {
        frac(self.miss_case, self.n_case)
    }
    #[cfg(test)]
    fn f_miss_u(&self) -> f64 {
        frac(self.miss_ctrl, self.n_ctrl)
    }
}

#[cfg(test)]
fn frac(num: u64, den: u64) -> f64 {
    if den == 0 {
        f64::NAN
    } else {
        #[allow(clippy::cast_precision_loss)]
        {
            num as f64 / den as f64
        }
    }
}

/// Lo-bit of each 2-bit lane across a 64-bit word.
const LANE_LO: u64 = 0x5555_5555_5555_5555;

/// Per-group lane masks: bit `2*s` set iff sample `s` is in that group. Laid out
/// like a packed genotype row so a row word ANDs directly against a mask word.
/// The `male` masks restrict each group to males for X/Y haploid handling.
pub struct StatusMasks {
    case: Vec<u64>,
    ctrl: Vec<u64>,
    case_male: Vec<u64>,
    ctrl_male: Vec<u64>,
    n_case: u64,
    n_ctrl: u64,
    n_case_male: u64,
    n_ctrl_male: u64,
}

impl StatusMasks {
    /// `sex` is PLINK's `.fam` code (1 = male, 2 = female, 0 = unknown).
    pub fn build(status: &[Status], sex: &[u8], bytes_per_variant: usize) -> Self {
        let words_per_row = bytes_per_variant.div_ceil(8);
        let mut case = vec![0u64; words_per_row];
        let mut ctrl = vec![0u64; words_per_row];
        let mut case_male = vec![0u64; words_per_row];
        let mut ctrl_male = vec![0u64; words_per_row];
        let mut n_case = 0u64;
        let mut n_ctrl = 0u64;
        let mut n_case_male = 0u64;
        let mut n_ctrl_male = 0u64;
        for (s, st) in status.iter().enumerate() {
            let (word, bit) = (s >> 5, (s & 31) * 2);
            let male = sex[s] == 1;
            match st {
                Status::Case => {
                    case[word] |= 1u64 << bit;
                    n_case += 1;
                    if male {
                        case_male[word] |= 1u64 << bit;
                        n_case_male += 1;
                    }
                }
                Status::Control => {
                    ctrl[word] |= 1u64 << bit;
                    n_ctrl += 1;
                    if male {
                        ctrl_male[word] |= 1u64 << bit;
                        n_ctrl_male += 1;
                    }
                }
                Status::Missing => {}
            }
        }
        Self {
            case,
            ctrl,
            case_male,
            ctrl_male,
            n_case,
            n_ctrl,
            n_case_male,
            n_ctrl_male,
        }
    }

    pub fn result(&self, row: &[u8], class: ChromClass) -> VariantResult {
        match class {
            ChromClass::Diploid => {
                let (miss_case, miss_ctrl) = scan_diploid(row, &self.case, &self.ctrl);
                VariantResult {
                    miss_case,
                    n_case: self.n_case,
                    miss_ctrl,
                    n_ctrl: self.n_ctrl,
                }
            }
            // X: all samples count; het males also count as missing.
            ChromClass::X => {
                let (miss_case, miss_ctrl) = scan_haploid(
                    row,
                    &self.case,
                    &self.ctrl,
                    &self.case_male,
                    &self.ctrl_male,
                );
                VariantResult {
                    miss_case,
                    n_case: self.n_case,
                    miss_ctrl,
                    n_ctrl: self.n_ctrl,
                }
            }
            // Y: females are dropped; only males count, het males are missing.
            ChromClass::Y => {
                let (miss_case, miss_ctrl) = scan_haploid(
                    row,
                    &self.case_male,
                    &self.ctrl_male,
                    &self.case_male,
                    &self.ctrl_male,
                );
                VariantResult {
                    miss_case,
                    n_case: self.n_case_male,
                    miss_ctrl,
                    n_ctrl: self.n_ctrl_male,
                }
            }
        }
    }
}

/// Pad the trailing partial row word with zero lanes (decode HomA1, never bad).
#[inline]
fn last_word(row: &[u8], full: usize) -> u64 {
    let mut bytes = [0u8; 8];
    bytes[..row.len() - full * 8].copy_from_slice(&row[full * 8..]);
    u64::from_le_bytes(bytes)
}

/// Autosomal/MT/XY missing tally: only the `0b01` no-call lane is bad. The
/// branch-free zip over full words lets the loop vectorise; the trailing partial
/// word is handled once.
#[inline]
fn scan_diploid(row: &[u8], case: &[u64], ctrl: &[u64]) -> (u64, u64) {
    let full = row.len() / 8;
    let mut mc = 0u32;
    let mut uc = 0u32;
    for (chunk, (&cm, &um)) in row[..full * 8]
        .chunks_exact(8)
        .zip(case.iter().zip(ctrl.iter()))
    {
        let word = u64::from_le_bytes(chunk.try_into().unwrap());
        let miss = word & !(word >> 1) & LANE_LO;
        mc += (miss & cm).count_ones();
        uc += (miss & um).count_ones();
    }
    if full < case.len() {
        let word = last_word(row, full);
        let miss = word & !(word >> 1) & LANE_LO;
        mc += (miss & case[full]).count_ones();
        uc += (miss & ctrl[full]).count_ones();
    }
    (u64::from(mc), u64::from(uc))
}

/// Sex-chromosome tally: the `0b01` no-call lane is bad for everyone in
/// `case`/`ctrl`, and the het `0b10` lane is additionally bad for the (male)
/// `case_h`/`ctrl_h` groups.
#[inline]
fn scan_haploid(
    row: &[u8],
    case: &[u64],
    ctrl: &[u64],
    case_h: &[u64],
    ctrl_h: &[u64],
) -> (u64, u64) {
    let full = row.len() / 8;
    let mut mc = 0u32;
    let mut uc = 0u32;
    let tally = |word: u64, w: usize, mc: &mut u32, uc: &mut u32| {
        let miss = word & !(word >> 1) & LANE_LO;
        let het = (word >> 1) & !word & LANE_LO;
        *mc += ((miss & case[w]) | (het & case_h[w])).count_ones();
        *uc += ((miss & ctrl[w]) | (het & ctrl_h[w])).count_ones();
    };
    for (w, chunk) in row[..full * 8].chunks_exact(8).enumerate() {
        tally(
            u64::from_le_bytes(chunk.try_into().unwrap()),
            w,
            &mut mc,
            &mut uc,
        );
    }
    if full < case.len() {
        tally(last_word(row, full), full, &mut mc, &mut uc);
    }
    (u64::from(mc), u64::from(uc))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pack(codes: &[u8]) -> Vec<u8> {
        let mut bytes = vec![0u8; codes.len().div_ceil(4)];
        for (s, &c) in codes.iter().enumerate() {
            bytes[s / 4] |= (c & 0b11) << ((s % 4) * 2);
        }
        bytes
    }

    // Genotype 2-bit codes: 0b00 HomA1, 0b10 Het, 0b11 HomA2, 0b01 missing.
    #[test]
    fn splits_missing_by_status() {
        let status = [
            Status::Case,
            Status::Control,
            Status::Case,
            Status::Control,
            Status::Case,
            Status::Control,
            Status::Case,
            Status::Control,
        ];
        let sex = [1u8; 8];
        let masks = StatusMasks::build(&status, &sex, 8usize.div_ceil(4));
        assert_eq!(masks.n_case, 4);
        assert_eq!(masks.n_ctrl, 4);

        // Missing at samples 0,2 (cases) and 3 (control).
        let row = pack(&[0b01, 0b00, 0b01, 0b01, 0b11, 0b10, 0b00, 0b11]);
        let r = masks.result(&row, ChromClass::Diploid);
        assert_eq!(r.miss_case, 2);
        assert_eq!(r.miss_ctrl, 1);
        assert!((r.f_miss_a() - 0.5).abs() < 1e-12);
        assert!((r.f_miss_u() - 0.25).abs() < 1e-12);
        assert!(r.nondegenerate());
    }

    #[test]
    fn fully_called_variant_is_degenerate() {
        let status = [Status::Case, Status::Control, Status::Case, Status::Control];
        let sex = [1u8; 4];
        let masks = StatusMasks::build(&status, &sex, 1);
        let row = pack(&[0b00, 0b10, 0b11, 0b00]);
        let r = masks.result(&row, ChromClass::Diploid);
        assert!(!r.nondegenerate());
        assert_eq!(r.f_miss_a(), 0.0);
        assert_eq!(r.p(&fisher::LnFactTable::new(4)), 1.0);
    }

    // On X, a het male is missing; a het female is a valid call.
    #[test]
    fn x_chrom_het_male_is_missing() {
        // 4 samples: case-male, case-female, ctrl-male, ctrl-female.
        let status = [Status::Case, Status::Case, Status::Control, Status::Control];
        let sex = [1u8, 2, 1, 2];
        let masks = StatusMasks::build(&status, &sex, 1);
        // All het.
        let row = pack(&[0b10, 0b10, 0b10, 0b10]);
        let dip = masks.result(&row, ChromClass::Diploid);
        assert_eq!(dip.miss_case, 0); // het valid when diploid
        let x = masks.result(&row, ChromClass::X);
        assert_eq!(x.miss_case, 1); // case-male het → missing; case-female het valid
        assert_eq!(x.miss_ctrl, 1); // ctrl-male het → missing; ctrl-female het valid
        assert_eq!(x.n_case, 2); // females still counted
        assert_eq!(x.n_ctrl, 2);
        assert!((x.f_miss_a() - 0.5).abs() < 1e-12);
    }

    // On X, a no-call (0b01) for a female still counts as missing.
    #[test]
    fn x_chrom_female_nocall_counts() {
        let status = [Status::Case, Status::Case, Status::Control, Status::Control];
        let sex = [1u8, 2, 1, 2]; // case-male, case-female, ctrl-male, ctrl-female
        let masks = StatusMasks::build(&status, &sex, 1);
        // case-female explicit missing; others homozygous.
        let row = pack(&[0b00, 0b01, 0b00, 0b00]);
        let x = masks.result(&row, ChromClass::X);
        assert_eq!(x.miss_case, 1);
        assert_eq!(x.miss_ctrl, 0);
        assert_eq!(x.n_case, 2);
    }

    // On Y, females are excluded; het males are missing.
    #[test]
    fn y_chrom_excludes_females() {
        let status = [Status::Case, Status::Case, Status::Control, Status::Control];
        let sex = [1u8, 2, 1, 2]; // case-male, case-female, ctrl-male, ctrl-female
        let masks = StatusMasks::build(&status, &sex, 1);
        // case-male homozygous (valid), ctrl-male missing, females homozygous.
        let row = pack(&[0b00, 0b00, 0b01, 0b00]);
        let y = masks.result(&row, ChromClass::Y);
        assert_eq!(y.n_case, 1); // only the case male
        assert_eq!(y.n_ctrl, 1); // only the control male
        assert_eq!(y.miss_case, 0);
        assert_eq!(y.miss_ctrl, 1);
        assert!((y.f_miss_u() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn spans_word_boundary() {
        let mut status = vec![Status::Control; 40];
        status[33] = Status::Case;
        status[35] = Status::Case;
        let sex = vec![1u8; 40];
        let masks = StatusMasks::build(&status, &sex, 40usize.div_ceil(4));
        assert_eq!(masks.n_case, 2);
        assert_eq!(masks.n_ctrl, 38);

        let mut codes = vec![0b11u8; 40];
        codes[33] = 0b01; // case missing
        codes[5] = 0b01; // control missing
        let row = pack(&codes);
        let r = masks.result(&row, ChromClass::Diploid);
        assert_eq!(r.miss_case, 1);
        assert_eq!(r.miss_ctrl, 1);
    }
}
