use criterion::{Criterion, criterion_group, criterion_main};
use rsomics_plink_test_missing::{MissingCode, PhenoOptions, open, run};
use std::hint::black_box;
use std::path::PathBuf;

fn bench_test_missing(c: &mut Criterion) {
    // Set RSOMICS_TESTMISSING_BFILE to a representative-large fileset prefix and
    // RSOMICS_TESTMISSING_PHENO to its case/control phenotype file; otherwise the
    // in-repo golden is used.
    let prefix = std::env::var("RSOMICS_TESTMISSING_BFILE")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/golden/small"));
    let pheno = std::env::var("RSOMICS_TESTMISSING_PHENO")
        .map(PathBuf::from)
        .ok();

    let opts = PhenoOptions {
        pheno_file: pheno.as_deref(),
        missing: MissingCode::default(),
        allow_no_sex: true,
        case_control_01: false,
    };

    c.bench_function("test_missing", |b| {
        b.iter(|| {
            let mut fs = open(&prefix).expect("open fileset");
            let mut sink = Vec::new();
            run(&mut fs, black_box(&opts), &mut sink).unwrap();
            black_box(sink);
        });
    });
}

criterion_group!(benches, bench_test_missing);
criterion_main!(benches);
