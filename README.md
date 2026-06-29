# rsomics-plink-test-missing

Differential genotype-missingness test between cases and controls — a
single-binary Rust port of `plink --test-missing`.

For each variant a 2×2 table of (missing vs non-missing genotype) × (case vs
control) is tested by a two-sided Fisher's exact test, surfacing
platform/batch artefacts where missingness differs by affection status. The
report follows PLINK's `.missing` column order.

## Usage

```
rsomics-plink-test-missing <bfile> --out result
rsomics-plink-test-missing <bfile> --pheno pheno.txt --out result
```

- `<bfile>` is the path prefix for a PLINK1 binary fileset (`.bed`/`.bim`/`.fam`),
  the same convention as `plink --bfile`.
- Affection status is read from `.fam` column 6 by default (`1 = control`,
  `2 = case`; `0`/`-9`/unmatched samples are dropped), or from a `--pheno FILE`
  (`FID IID value`, optional `FID IID` header).
- `--1` switches the coding to `0 = control`, `1 = case` (matching `plink --1`).
- `--allow-no-sex` keeps samples with unknown sex (matching `plink --allow-no-sex`).
- `--missing-phenotype N` sets the missing-value integer code (default `-9`).
- `-t/--threads N` sets the worker-thread count (per-variant work is parallel).

Output is `<out>.missing` with columns `CHR SNP F_MISS_A F_MISS_U P`, one line
per nondegenerate variant (a variant with at least one missing call across
cases and controls). `F_MISS_A` is the missing-call fraction in cases,
`F_MISS_U` in controls, and `P` is the two-sided Fisher's exact p-value.

## Method

The `.bed` is read in variant-major blocks and each block is processed in
parallel, so resident memory stays near one block regardless of fileset size.
Case- and control-membership are precomputed as packed lane masks laid out like
a genotype row; per variant the row's missing lanes (the 2-bit code `0b01`) are
masked against each group's lane mask and popcounted, giving the missing counts
in cases and controls in one sweep.

The 2×2 table's two-sided Fisher's exact p-value sums every hypergeometric table
no more likely than the observed one. A log-factorial table built once (up to the
sample count) makes each test transcendental-free; the mode's probability is
exponentiated once and the rest of the tail is reached by the exact
multiplicative recurrence `P(x±1)/P(x)`, stopping when a tail term is too small to
change the running sum. F_MISS fractions and per-table p-values are memoised on
their counts, so the repeated tables a real fileset produces are formatted and
tested once.

## Origin

This crate is an independent Rust reimplementation of PLINK 1.9 `--test-missing`
based on:

- Purcell et al. 2007, "PLINK: a tool set for whole-genome association and
  population-based linkage analyses" (Am J Hum Genet, doi:10.1086/519795)
- Chang et al. 2015, "Second-generation PLINK" (GigaScience,
  doi:10.1186/s13742-015-0047-8)
- The PLINK 1.9 `--test-missing` / `.missing` documentation
  (<https://www.cog-genomics.org/plink/1.9/assoc>,
  <https://www.cog-genomics.org/plink/1.9/formats#missing>)
- The public PLINK 1.9 binary fileset specification
  (<https://www.cog-genomics.org/plink/1.9/formats>)
- Black-box behaviour testing against the `plink` 1.9 binary

No source code from the GPL PLINK upstream was used as reference during
implementation. The two-sided Fisher's exact test is the standard
probability-convention p-value (the path SciPy's `fisher_exact` and R's
`fisher.test` reduce to). Test fixtures are independently generated; the
committed golden expectations were produced once by running PLINK 1.9 on the
fixture so the compatibility test runs without PLINK installed.

License: MIT OR Apache-2.0.
Upstream credit: PLINK 1.9 (Christopher Chang et al., GPLv3) —
<https://www.cog-genomics.org/plink/>.
