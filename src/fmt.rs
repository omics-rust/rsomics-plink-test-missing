//! Faithful C `printf("%.4g")` — PLINK's number format for the `.missing`
//! F_MISS_A / F_MISS_U / P columns.
//!
//! Keeps 4 significant figures, switches to scientific notation when the decimal
//! exponent is `< -4` or `>= 4`, strips trailing zeros, and pads the exponent to
//! two signed digits (`5e-05`, `1.036e-45`).

const P: i32 = 4;

/// Append the `%.4g` form of `x` to `out`, using `scratch` as reusable working
/// space so the hot path allocates nothing.
pub fn g_into(out: &mut String, scratch: &mut String, x: f64) {
    use std::fmt::Write as _;
    if x.is_nan() {
        out.push_str("nan");
        return;
    }
    if x == 0.0 {
        out.push('0');
        return;
    }
    if x.is_infinite() {
        out.push_str(if x < 0.0 { "-inf" } else { "inf" });
        return;
    }

    scratch.clear();
    write!(scratch, "{:.*e}", (P - 1) as usize, x).unwrap();
    let (mant, exp_str) = scratch.split_once('e').expect("scientific form has e");
    let exp: i32 = exp_str.parse().expect("scientific form has exponent");

    if !(-4..P).contains(&exp) {
        let mant = strip_trailing(mant);
        write!(
            out,
            "{mant}e{}{:02}",
            if exp < 0 { '-' } else { '+' },
            exp.abs()
        )
        .unwrap();
    } else {
        let frac = (P - 1 - exp).max(0) as usize;
        let mant_start = out.len();
        write!(out, "{x:.frac$}").unwrap();
        strip_trailing_in_place(out, mant_start);
    }
}

/// Trim trailing zeros (and a dangling point) from a fractional decimal.
fn strip_trailing(s: &str) -> &str {
    if s.contains('.') {
        s.trim_end_matches('0').trim_end_matches('.')
    } else {
        s
    }
}

/// Same trim, applied only to the decimal written into `out` from `start`.
fn strip_trailing_in_place(out: &mut String, start: usize) {
    if out[start..].contains('.') {
        let kept = start
            + out[start..]
                .trim_end_matches('0')
                .trim_end_matches('.')
                .len();
        out.truncate(kept);
    }
}

#[cfg(test)]
mod tests {
    use super::g_into;

    fn g(x: f64) -> String {
        let mut s = String::new();
        let mut scratch = String::new();
        g_into(&mut s, &mut scratch, x);
        s
    }

    // Reference values from C printf("%.4g", x).
    #[test]
    fn matches_c_printf_4g() {
        assert_eq!(g(0.0), "0");
        assert_eq!(g(1.0), "1");
        assert_eq!(g(0.533_333_333_3), "0.5333");
        assert_eq!(g(0.133_333_333_3), "0.1333");
        assert_eq!(g(0.050_174_912_54), "0.05017");
        assert_eq!(g(0.651_340_996_2), "0.6513");
        assert_eq!(g(0.245_077_461_3), "0.2451");
        assert_eq!(g(0.006_754_674_575), "0.006755");
        assert_eq!(g(0.2), "0.2");
        assert_eq!(g(5e-5), "5e-05");
        assert_eq!(g(1.035_856_563_790_026_7e-45), "1.036e-45");
    }
}
