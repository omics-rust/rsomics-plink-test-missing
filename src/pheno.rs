//! Case/control affection loading with PLINK's sample-matching semantics.
//!
//! Affection comes from `.fam` column 6 by default, or a `--pheno` file
//! (`FID IID value`, optional `FID IID` header). PLINK's default coding is
//! `1 = control`, `2 = case`; `0` and the missing code (`-9`) are unaffected by
//! the test and the sample is dropped. `--1` shifts the coding to `0 = control`,
//! `1 = case`. A sample whose phenotype is missing, or — without `--allow-no-sex`
//! — whose sex is unknown, is excluded.

use crate::bed::Sample;
use anyhow::{Context, Result, bail};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// The integer code PLINK treats as a missing phenotype (`-9` default).
#[derive(Clone, Copy)]
pub struct MissingCode(pub i64);

impl Default for MissingCode {
    fn default() -> Self {
        MissingCode(-9)
    }
}

/// Case/control status for one sample.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Status {
    Case,
    Control,
    /// Missing/excluded — does not enter either count.
    Missing,
}

pub struct PhenoOptions<'a> {
    pub pheno_file: Option<&'a Path>,
    pub missing: MissingCode,
    pub allow_no_sex: bool,
    /// `--1`: control = 0, case = 1 (default is control = 1, case = 2).
    pub case_control_01: bool,
}

/// Per-sample case/control status in fileset order. Errors when no sample is a
/// case or no sample is a control (the test would be undefined for every
/// variant).
pub fn load_status(samples: &[Sample], opts: &PhenoOptions) -> Result<Vec<Status>> {
    let raw: Vec<String> = match opts.pheno_file {
        Some(path) => from_value_file(samples, path)?,
        None => samples.iter().map(|s| s.phen.clone()).collect(),
    };

    let mut status: Vec<Status> = raw
        .iter()
        .map(|v| classify(v, opts.missing, opts.case_control_01))
        .collect();

    if !opts.allow_no_sex {
        for (s, st) in samples.iter().zip(status.iter_mut()) {
            if s.sex == 0 {
                *st = Status::Missing;
            }
        }
    }

    let n_case = status.iter().filter(|s| **s == Status::Case).count();
    let n_ctrl = status.iter().filter(|s| **s == Status::Control).count();
    if n_case == 0 || n_ctrl == 0 {
        bail!(
            "--test-missing needs both cases and controls (found {n_case} cases, {n_ctrl} controls)"
        );
    }
    Ok(status)
}

fn classify(raw: &str, missing: MissingCode, case_control_01: bool) -> Status {
    if raw.eq_ignore_ascii_case("na") || raw.parse::<i64>() == Ok(missing.0) {
        return Status::Missing;
    }
    let (case, control) = if case_control_01 {
        ("1", "0")
    } else {
        ("2", "1")
    };
    if raw == case {
        Status::Case
    } else if raw == control {
        Status::Control
    } else {
        Status::Missing
    }
}

/// Read the first value column of a keyed phenotype file, aligned to fileset
/// order; samples absent from the file get the missing code so they drop out.
fn from_value_file(samples: &[Sample], path: &Path) -> Result<Vec<String>> {
    let f = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let mut rows: HashMap<(String, String), String> = HashMap::new();
    for (lineno, line) in BufReader::new(f).lines().enumerate() {
        let line = line?;
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.is_empty() {
            continue;
        }
        if lineno == 0
            && fields.len() >= 2
            && fields[0].eq_ignore_ascii_case("fid")
            && fields[1].eq_ignore_ascii_case("iid")
        {
            continue;
        }
        if fields.len() < 3 {
            bail!(
                "{}:{}: expected FID IID value, got {} fields",
                path.display(),
                lineno + 1,
                fields.len()
            );
        }
        rows.insert(
            (fields[0].to_string(), fields[1].to_string()),
            fields[2].to_string(),
        );
    }
    Ok(samples
        .iter()
        .map(|s| {
            rows.get(&(s.fid.clone(), s.iid.clone()))
                .cloned()
                .unwrap_or_else(|| "-9".to_string())
        })
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn sample(fid: &str, iid: &str, sex: u8, phen: &str) -> Sample {
        Sample {
            fid: fid.into(),
            iid: iid.into(),
            sex,
            phen: phen.into(),
        }
    }

    fn opts<'a>(pheno: Option<&'a Path>, allow_no_sex: bool) -> PhenoOptions<'a> {
        PhenoOptions {
            pheno_file: pheno,
            missing: MissingCode::default(),
            allow_no_sex,
            case_control_01: false,
        }
    }

    #[test]
    fn fam_affection_default_coding() {
        let s = [
            sample("F", "A", 1, "2"),
            sample("F", "B", 1, "1"),
            sample("F", "C", 1, "0"),
            sample("F", "D", 1, "-9"),
        ];
        let st = load_status(&s, &opts(None, true)).unwrap();
        assert_eq!(
            st,
            vec![
                Status::Case,
                Status::Control,
                Status::Missing,
                Status::Missing
            ]
        );
    }

    #[test]
    fn one_flag_shifts_coding() {
        let s = [sample("F", "A", 1, "1"), sample("F", "B", 1, "0")];
        let mut o = opts(None, true);
        o.case_control_01 = true;
        let st = load_status(&s, &o).unwrap();
        assert_eq!(st, vec![Status::Case, Status::Control]);
    }

    #[test]
    fn no_sex_dropped_without_flag() {
        // The sex-0 case drops; a sexed case and control keep the test defined.
        let s = [
            sample("F", "A", 0, "2"),
            sample("F", "B", 1, "1"),
            sample("F", "C", 2, "2"),
        ];
        let st = load_status(&s, &opts(None, false)).unwrap();
        assert_eq!(st, vec![Status::Missing, Status::Control, Status::Case]);
    }

    #[test]
    fn pheno_file_overrides_fam() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("ph.txt");
        let mut f = File::create(&p).unwrap();
        writeln!(f, "FID IID PHE").unwrap();
        writeln!(f, "F A 2").unwrap();
        writeln!(f, "F B 1").unwrap();
        let s = [sample("F", "A", 1, "-9"), sample("F", "B", 1, "-9")];
        let st = load_status(&s, &opts(Some(&p), true)).unwrap();
        assert_eq!(st, vec![Status::Case, Status::Control]);
    }

    #[test]
    fn errors_when_no_controls() {
        let s = [sample("F", "A", 1, "2"), sample("F", "B", 1, "2")];
        assert!(load_status(&s, &opts(None, true)).is_err());
    }
}
