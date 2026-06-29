//! Two-sided Fisher's exact test on a 2×2 table.
//!
//! PLINK `--test-missing` reports the two-sided Fisher p-value for the table
//! [case, control] × [missing, non-missing]. The two-sided p-value is the sum
//! of the hypergeometric probabilities of every table (with the same margins) no
//! more likely than the observed one — the standard "probability" convention
//! shared by R's `fisher.test` and SciPy's `fisher_exact(alternative="two-sided")`.
//!
//! Log-factorials are precomputed once (`LnFactTable`), so the mode's
//! probability is a handful of table lookups and an `exp`; every other table is
//! reached by the exact multiplicative recurrence `P(x±1)/P(x)` (one multiply and
//! divide, no transcendental). Probabilities are summed with a relative cutoff
//! (`1 + 1e-7`) on the ratio to the observed table — matching PLINK's own
//! tolerance so borderline tables fall on the same side of the boundary — and
//! the monotone tails stop once a term is too small to move the running sum.

const EPSILON: f64 = 1.0 + 1e-7;

/// Precomputed `ln(k!)` for `k` up to the sample count, so a per-variant test
/// needs no `lgamma` call at all — the five log-factorials of its margins are
/// table lookups.
pub struct LnFactTable(Vec<f64>);

impl LnFactTable {
    pub fn new(max: usize) -> Self {
        let mut t = Vec::with_capacity(max + 1);
        let mut acc = 0.0;
        t.push(0.0); // ln(0!) = 0
        for k in 1..=max {
            #[allow(clippy::cast_precision_loss)]
            {
                acc += (k as f64).ln();
            }
            t.push(acc);
        }
        Self(t)
    }

    #[inline]
    fn get(&self, k: u64) -> f64 {
        self.0[k as usize]
    }
}

/// Two-sided Fisher's exact p-value for the 2×2 table `[[a, b], [c, d]]`,
/// with the margins held fixed, using a precomputed log-factorial table.
/// Returns 1.0 for a degenerate table (an empty row or column).
#[must_use]
pub fn two_sided(lf: &LnFactTable, a: u64, b: u64, c: u64, d: u64) -> f64 {
    let row1 = a + b;
    let row2 = c + d;
    let col1 = a + c;
    let col2 = b + d;
    let n = row1 + row2;
    if row1 == 0 || row2 == 0 || col1 == 0 || col2 == 0 {
        return 1.0;
    }

    // Free corner `x` (the top-left cell) ranges over [lo, hi]; the other cells
    // are determined: top-right = row1-x, bottom-left = col1-x,
    // bottom-right = row2-(col1-x).
    let lo = col1.saturating_sub(row2);
    let hi = col1.min(row1);

    // Hypergeometric mode, where the unnormalised probability peaks.
    let mode = ((row1 + 1) * (col1 + 1) / (n + 2)).clamp(lo, hi);
    let ln_const = lf.get(row1) + lf.get(row2) + lf.get(col1) + lf.get(col2) - lf.get(n);
    let p_mode = (ln_const
        - lf.get(mode)
        - lf.get(col1 - mode)
        - lf.get(row1 - mode)
        - lf.get(row2 - (col1 - mode)))
    .exp();

    let p_obs = p_mode * ratio_to_mode(a, mode, row1, row2, col1);
    let cutoff = p_obs * EPSILON;

    let mut sum = if p_mode <= cutoff { p_mode } else { 0.0 };
    // Both walks leave the mode monotonically decreasing, so once a term is below
    // the cutoff it stays there, and once it is too small to change the running
    // sum in f64 the rest of the tail is negligible — stop there.
    //
    // Walk up: P(x+1) = P(x) · (row1-x)(col1-x) / ((x+1)(row2-col1+x+1)).
    let mut p = p_mode;
    for x in mode..hi {
        let d_cell = row2 - (col1 - x); // current bottom-right
        p *= f(row1 - x) * f(col1 - x) / (f(x + 1) * f(d_cell + 1));
        if p <= cutoff {
            let next = sum + p;
            if next == sum {
                break;
            }
            sum = next;
        }
    }
    // Walk down: P(x-1) = P(x) · x(row2-col1+x) / ((row1-x+1)(col1-x+1)).
    p = p_mode;
    for x in (lo + 1..=mode).rev() {
        let d_cell = row2 - (col1 - x);
        p *= f(x) * f(d_cell) / (f(row1 - x + 1) * f(col1 - x + 1));
        if p <= cutoff {
            let next = sum + p;
            if next == sum {
                break;
            }
            sum = next;
        }
    }
    sum.min(1.0)
}

#[inline]
#[allow(clippy::cast_precision_loss)]
fn f(x: u64) -> f64 {
    x as f64
}

/// `P(x) / P(mode)` for the hypergeometric, by the exact multiplicative walk —
/// used to place the observed table without a second `lgamma`.
fn ratio_to_mode(x: u64, mode: u64, row1: u64, row2: u64, col1: u64) -> f64 {
    let mut r = 1.0;
    if x > mode {
        for k in mode..x {
            let d_cell = row2 - (col1 - k);
            r *= f(row1 - k) * f(col1 - k) / (f(k + 1) * f(d_cell + 1));
        }
    } else {
        for k in (x + 1..=mode).rev() {
            let d_cell = row2 - (col1 - k);
            r *= f(k) * f(d_cell) / (f(row1 - k + 1) * f(col1 - k + 1));
        }
    }
    r
}

#[cfg(test)]
mod tests {
    use super::{LnFactTable, two_sided};

    fn close(got: f64, want: f64) {
        let rel = (got - want).abs() / want.abs().max(f64::MIN_POSITIVE);
        assert!(rel <= 1e-12, "got {got:e} want {want:e} rel {rel:e}");
    }

    fn ts(a: u64, b: u64, c: u64, d: u64) -> f64 {
        let lf = LnFactTable::new((a + b + c + d) as usize);
        two_sided(&lf, a, b, c, d)
    }

    // Reference values from scipy.stats.fisher_exact(..., alternative="two-sided").
    #[test]
    fn matches_scipy_two_sided() {
        close(ts(1, 9, 11, 3), 0.002_759_456_185_220_083_6);
        close(ts(8, 2, 1, 5), 0.034_965_034_965_034_975);
        close(ts(10, 10, 10, 10), 1.0);
        close(ts(5, 0, 0, 5), 0.007_936_507_936_507_938);
        close(ts(2, 8, 7, 3), 0.069_778_518_694_927_36);
        close(ts(100, 50, 40, 110), 4.467_383_093_506_087e-12);
        close(ts(3, 1, 1, 3), 0.485_714_285_714_285_65);
    }

    #[test]
    fn degenerate_tables_are_one() {
        assert_eq!(ts(0, 0, 5, 5), 1.0);
        assert_eq!(ts(5, 5, 0, 0), 1.0);
        assert_eq!(ts(0, 5, 0, 5), 1.0);
        assert_eq!(ts(5, 0, 5, 0), 1.0);
    }

    // No missing genotypes at all in either group → P = 1 (PLINK also skips the
    // row, but the test itself must be well-defined).
    #[test]
    fn no_missing_is_one() {
        assert_eq!(ts(0, 100, 0, 80), 1.0);
    }
}
