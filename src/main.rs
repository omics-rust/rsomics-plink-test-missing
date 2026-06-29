use clap::Parser;
use rsomics_plink_test_missing::{MissingCode, PhenoOptions, open, run as run_scan};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::process::ExitCode;

#[derive(Parser)]
#[command(
    name = "rsomics-plink-test-missing",
    about = "Differential genotype-missingness test between cases and controls (plink --test-missing)",
    version
)]
struct Cli {
    /// Path prefix for the .bed/.bim/.fam fileset (without extension).
    bfile: PathBuf,

    /// Output prefix; writes <OUT>.missing (plink --out).
    #[arg(long, default_value = "plink")]
    out: PathBuf,

    /// Case/control phenotype file (FID IID value, optional FID/IID header);
    /// without it the affection is read from .fam column 6.
    #[arg(long)]
    pheno: Option<PathBuf>,

    /// Phenotype is 0 = control, 1 = case (default is 1 = control, 2 = case).
    #[arg(long = "1")]
    case_control_01: bool,

    /// Integer code marking a missing phenotype value.
    #[arg(long, default_value_t = -9)]
    missing_phenotype: i64,

    /// Keep samples with unknown sex (sex code 0) instead of dropping them.
    #[arg(long)]
    allow_no_sex: bool,

    /// Number of rayon worker threads (0 = rayon default).
    #[arg(short = 't', long, default_value_t = 0)]
    threads: usize,
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e:#}");
            ExitCode::FAILURE
        }
    }
}

fn run(cli: Cli) -> anyhow::Result<()> {
    if cli.threads != 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(cli.threads)
            .build_global()
            .ok();
    }

    let mut fs = open(&cli.bfile)?;
    let opts = PhenoOptions {
        pheno_file: cli.pheno.as_deref(),
        missing: MissingCode(cli.missing_phenotype),
        allow_no_sex: cli.allow_no_sex,
        case_control_01: cli.case_control_01,
    };

    let path = cli.out.with_extension("missing");
    let mut w = BufWriter::with_capacity(1 << 20, File::create(&path)?);
    run_scan(&mut fs, &opts, &mut w)?;
    w.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cli_definition_is_valid() {
        <Cli as clap::CommandFactory>::command().debug_assert();
    }
}
