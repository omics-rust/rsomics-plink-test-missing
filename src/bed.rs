//! Streaming PLINK1 `.bed`/`.bim`/`.fam` reader.
//!
//! Genotypes are read from the `.bed` in variant-major blocks rather than mapped
//! or slurped, so resident memory stays near one block of rows regardless of
//! fileset size. The tiny `.bim`/`.fam` are parsed up front.

use anyhow::{Context, Result, bail};
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

pub struct Variant {
    pub chrom: String,
    pub id: String,
}

pub struct Sample {
    pub fid: String,
    pub iid: String,
    pub sex: u8,
    pub phen: String,
}

pub struct Fileset {
    pub variants: Vec<Variant>,
    pub samples: Vec<Sample>,
    pub bytes_per_variant: usize,
    bed: BufReader<File>,
}

const MAGIC: [u8; 3] = [0x6c, 0x1b, 0x01];

impl Fileset {
    pub fn open(prefix: &Path) -> Result<Self> {
        let variants = parse_bim(&prefix.with_extension("bim"))?;
        let samples = parse_fam(&prefix.with_extension("fam"))?;
        let bytes_per_variant = samples.len().div_ceil(4);

        let path = prefix.with_extension("bed");
        let mut file = File::open(&path).with_context(|| format!("open {}", path.display()))?;
        let mut magic = [0u8; 3];
        file.read_exact(&mut magic)
            .with_context(|| format!("read .bed header from {}", path.display()))?;
        if magic[..2] != MAGIC[..2] {
            bail!("{} is not a PLINK .bed (bad magic)", path.display());
        }
        if magic[2] != 0x01 {
            bail!(
                "{} is not variant-major (only SNP-major .bed supported)",
                path.display()
            );
        }

        let expected = bytes_per_variant as u64 * variants.len() as u64;
        let have = file.metadata()?.len() - 3;
        if have != expected {
            bail!(
                "{}: genotype payload is {have} bytes, expected {expected} for {} variants × {} samples",
                path.display(),
                variants.len(),
                samples.len()
            );
        }

        Ok(Self {
            variants,
            samples,
            bytes_per_variant,
            bed: BufReader::with_capacity(1 << 20, file),
        })
    }

    pub fn n_variants(&self) -> usize {
        self.variants.len()
    }

    pub fn n_samples(&self) -> usize {
        self.samples.len()
    }

    /// Read the next `max` variants' packed rows into `buf`, returning how many
    /// were read (0 at end of file). `buf` is resized to `count × bpv`.
    pub fn read_block(&mut self, max: usize, buf: &mut Vec<u8>) -> Result<usize> {
        let bpv = self.bytes_per_variant;
        buf.resize(max * bpv, 0);
        let mut filled = 0;
        while filled < max * bpv {
            let n = self.bed.read(&mut buf[filled..])?;
            if n == 0 {
                break;
            }
            filled += n;
        }
        if filled % bpv != 0 {
            bail!("truncated .bed: partial variant row");
        }
        buf.truncate(filled);
        Ok(filled / bpv)
    }
}

fn parse_bim(path: &Path) -> Result<Vec<Variant>> {
    let f = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let mut out = Vec::new();
    for (lineno, line) in BufReader::new(f).lines().enumerate() {
        let line = line?;
        let t: Vec<&str> = line.split_whitespace().collect();
        if t.is_empty() {
            continue;
        }
        if t.len() < 6 {
            bail!(
                "{}:{}: expected 6 columns, got {}",
                path.display(),
                lineno + 1,
                t.len()
            );
        }
        out.push(Variant {
            chrom: t[0].to_string(),
            id: t[1].to_string(),
        });
    }
    Ok(out)
}

fn parse_fam(path: &Path) -> Result<Vec<Sample>> {
    let f = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let mut out = Vec::new();
    for (lineno, line) in BufReader::new(f).lines().enumerate() {
        let line = line?;
        let t: Vec<&str> = line.split_whitespace().collect();
        if t.is_empty() {
            continue;
        }
        if t.len() < 6 {
            bail!(
                "{}:{}: expected 6 columns, got {}",
                path.display(),
                lineno + 1,
                t.len()
            );
        }
        out.push(Sample {
            fid: t[0].to_string(),
            iid: t[1].to_string(),
            // plink treats any sex code other than 1/2 as unknown, including
            // the standard missing sentinel -9 and 0.
            sex: match t[4].trim() {
                "1" => 1,
                "2" => 2,
                _ => 0,
            },
            phen: t[5].to_string(),
        });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn reads_blocks_matching_layout() {
        let dir = tempfile::tempdir().unwrap();
        let prefix = dir.path().join("t");
        let mut bim = File::create(prefix.with_extension("bim")).unwrap();
        for v in 0..5 {
            writeln!(bim, "1 snp{v} 0 {} A G", 100 + v).unwrap();
        }
        let mut fam = File::create(prefix.with_extension("fam")).unwrap();
        for s in 0..6 {
            writeln!(fam, "F{s} I{s} 0 0 1 1").unwrap();
        }
        let bpv = 6usize.div_ceil(4); // 2
        let mut bed = File::create(prefix.with_extension("bed")).unwrap();
        bed.write_all(&MAGIC).unwrap();
        for v in 0..5u8 {
            bed.write_all(&[v, v.wrapping_add(1)]).unwrap();
        }
        drop(bed);

        let mut fs = Fileset::open(&prefix).unwrap();
        assert_eq!(fs.n_variants(), 5);
        assert_eq!(fs.n_samples(), 6);
        assert_eq!(fs.bytes_per_variant, bpv);

        let mut buf = Vec::new();
        assert_eq!(fs.read_block(3, &mut buf).unwrap(), 3);
        assert_eq!(buf, vec![0, 1, 1, 2, 2, 3]);
        assert_eq!(fs.read_block(3, &mut buf).unwrap(), 2);
        assert_eq!(buf, vec![3, 4, 4, 5]);
        assert_eq!(fs.read_block(3, &mut buf).unwrap(), 0);
    }
}
