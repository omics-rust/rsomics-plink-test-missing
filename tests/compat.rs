//! Differential compatibility against PLINK 1.9 `--test-missing`.
//!
//! The checked-in golden `small.missing` was produced once by running PLINK 1.9
//! on the committed `small.{bed,bim,fam}` fileset, so this test needs no PLINK at
//! run time. We run our binary on the same fileset and compare the `.missing`
//! report. Output is byte-identical to PLINK's (verified at authoring time on a
//! 500×400 fixture); the byte check asserts that, and a field-level check
//! pinpoints any future divergence (CHR/SNP exact; F_MISS_A/F_MISS_U exact at
//! PLINK's `%.4g`; P within one unit of PLINK's last printed digit).

use std::path::{Path, PathBuf};
use std::process::Command;

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-plink-test-missing"))
}

fn golden_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden")
}

fn run_ours(prefix: &Path, out_prefix: &Path) {
    let status = Command::new(ours())
        .arg(prefix)
        .args([
            "--allow-no-sex".as_ref(),
            "--out".as_ref(),
            out_prefix.as_os_str(),
            "-t1".as_ref(),
        ])
        .status()
        .expect("run rsomics-plink-test-missing");
    assert!(status.success(), "rsomics-plink-test-missing failed");
}

fn rows(text: &str) -> Vec<Vec<String>> {
    text.lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| l.split_whitespace().map(str::to_string).collect())
        .collect()
}

/// Significant figures PLINK prints in `.missing` (its `%.4g` width).
const PLINK_SIG_FIGS: i32 = 4;

/// Numeric agreement at PLINK's display precision: our value sits within one
/// unit of PLINK's last printed significant figure.
fn num_match(ours: &str, golden: &str) -> bool {
    let (a, b) = (ours.parse::<f64>().unwrap(), golden.parse::<f64>().unwrap());
    if b == 0.0 {
        return a.abs() < 1e-12;
    }
    let ulp = 10f64.powf(b.abs().log10().floor() - f64::from(PLINK_SIG_FIGS - 1));
    (a - b).abs() <= ulp + b.abs() * 1e-12
}

#[test]
fn matches_plink_golden() {
    let dir = golden_dir();
    let scratch = tempfile::tempdir().expect("scratch dir");
    let out = scratch.path().join("ours");
    run_ours(&dir.join("small"), &out);

    let got = std::fs::read_to_string(out.with_extension("missing")).expect("read ours");
    let want = std::fs::read_to_string(dir.join("small.missing")).expect("read golden");

    assert_eq!(got, want, "byte-level mismatch vs PLINK golden");

    let g = rows(&got);
    let w = rows(&want);
    assert_eq!(
        g.len(),
        w.len(),
        "row count: ours {} golden {}",
        g.len(),
        w.len()
    );
    assert_eq!(g[0], w[0], "header mismatch");

    for (gi, wi) in g.iter().zip(&w).skip(1) {
        assert_eq!(gi.len(), 5, "ours row has {} fields", gi.len());
        assert_eq!(wi.len(), 5, "golden row has {} fields", wi.len());
        // CHR SNP — exact.
        assert_eq!(gi[0], wi[0], "CHR: ours {gi:?} golden {wi:?}");
        assert_eq!(gi[1], wi[1], "SNP: ours {gi:?} golden {wi:?}");
        // F_MISS_A F_MISS_U P — at PLINK's display precision.
        for (name, c) in [("F_MISS_A", 2), ("F_MISS_U", 3), ("P", 4)] {
            assert!(
                num_match(&gi[c], &wi[c]),
                "{name}: ours {} golden {} (row {gi:?})",
                gi[c],
                wi[c]
            );
        }
    }
}
