//! Special functions needed to reproduce limma's moderated t-test.
//!
//! `diffpeak()` delegates its statistics to limma, so matching it numerically
//! means implementing the handful of special functions limma's empirical-Bayes
//! step calls. Each one is written here rather than pulled from a crate so the
//! algorithms match limma's own (which reach `statmod` for `trigammaInverse`),
//! and each is pinned against an R oracle in `testdata/`.
//!
//! The functions, and where limma uses them:
//!
//! - [`digamma`] / [`trigamma`] / [`tetragamma`]: `fitFDist()` estimates the
//!   prior degrees of freedom from `mean(trigamma(df1/2))`, and needs the
//!   second derivative for its Newton step.
//! - [`logmdigamma`]: `log(x) - digamma(x)`, used to put the variance estimates
//!   on the log scale before their mean and spread are taken.
//! - [`trigamma_inverse`]: inverts trigamma to turn a variance-of-logs into
//!   `df2`, limma's prior degrees of freedom.
//! - [`t_two_sided_p`]: the moderated p-value.
//! - [`regularized_incomplete_beta`] / [`betacf`]: what the t tail reduces to.
//!
//! All are continuous-argument; limma only ever calls them with positive `x`.

/// Lanczos coefficients for `ln(gamma(x))`, `g = 7`, `n = 9`.
///
/// The same set Numerical Recipes uses; it is accurate to about 15 significant
/// digits for `x > 0`, which is the range limma works in. Quoted verbatim, so
/// clippy's precision lint is allowed here: shortening them would obscure the
/// provenance and invite a transcription error.
#[allow(clippy::excessive_precision)]
const LANCZOS_COEFFICIENTS: [f64; 9] = [
    0.999_999_999_999_809_93,
    676.520_368_121_885_1,
    -1_259.139_216_722_402_8,
    771.323_428_777_653_13,
    -176.615_029_162_140_6,
    12.507_343_278_686_905,
    -0.138_571_095_265_720_12,
    9.984_369_578_019_572e-6,
    1.505_632_735_149_311_6e-7,
];

/// `ln(gamma(x))` for `x > 0` via the Lanczos approximation.
pub fn ln_gamma(x: f64) -> f64 {
    if x < 0.5 {
        // Reflection keeps the series in its accurate range.
        let pi = std::f64::consts::PI;
        return (pi / (pi * x).sin()).ln() - ln_gamma(1.0 - x);
    }
    let shifted = x - 1.0;
    let mut series = LANCZOS_COEFFICIENTS[0];
    let t = shifted + 7.5;
    for (index, coefficient) in LANCZOS_COEFFICIENTS.iter().enumerate().skip(1) {
        series += coefficient / (shifted + index as f64);
    }
    0.5 * (2.0 * std::f64::consts::PI).ln() + (shifted + 0.5) * t.ln() - t + series.ln()
}

/// Argument shift for the digamma/trigamma recurrences.
///
/// Chosen by measuring the asymptotic series against the exact values
/// `digamma(1) = -gamma` and `trigamma(1) = pi^2/6`: 6 gives ~1e-10, 10 gives
/// ~7e-13, and 20 reaches ~3e-16, i.e. machine precision.
const SERIES_SHIFT: f64 = 20.0;

/// Digamma function `psi(x) = d/dx ln(gamma(x))` for `x > 0`.
///
/// The recurrence pushes the argument up to `SERIES_SHIFT` and the asymptotic
/// series finishes the job. The shift is not cosmetic: at 6 the series is only
/// accurate to ~1e-10, which is far short of the double precision limma's
/// variance fit needs (`evar` is a difference of quantities that nearly cancel,
/// so a 1e-10 error there visibly moves `df.prior`).
pub fn digamma(x: f64) -> f64 {
    let mut value = 0.0;
    let mut argument = x;
    // Recurrence: shift upward until the asymptotic series converges quickly.
    while argument < SERIES_SHIFT {
        value -= 1.0 / argument;
        argument += 1.0;
    }
    // Asymptotic expansion in 1/x.
    let inverse = 1.0 / argument;
    let inverse_squared = inverse * inverse;
    value += argument.ln()
        - 0.5 * inverse
        - inverse_squared
            * (1.0 / 12.0
                - inverse_squared
                    * (1.0 / 120.0
                        - inverse_squared * (1.0 / 252.0 - inverse_squared * (1.0 / 240.0))));
    value
}

/// Trigamma function `psi'(x)` for `x > 0`.
pub fn trigamma(x: f64) -> f64 {
    let mut value = 0.0;
    let mut argument = x;
    while argument < SERIES_SHIFT {
        value += 1.0 / (argument * argument);
        argument += 1.0;
    }
    let inverse = 1.0 / argument;
    let inverse_squared = inverse * inverse;
    value
        + inverse
            * (1.0
                + 0.5 * inverse
                + inverse_squared
                    * (1.0 / 6.0
                        - inverse_squared
                            * (1.0 / 30.0
                                - inverse_squared
                                    * (1.0 / 42.0 - inverse_squared * (1.0 / 30.0)))))
}

/// Tetragamma function `psi''(x)` for `x > 0`, the derivative Newton's method
/// needs to invert [`trigamma`].
pub fn tetragamma(x: f64) -> f64 {
    let mut value = 0.0;
    let mut argument = x;
    while argument < SERIES_SHIFT {
        value -= 2.0 / (argument * argument * argument);
        argument += 1.0;
    }
    let inverse = 1.0 / argument;
    let inverse_squared = inverse * inverse;
    value
        - inverse_squared
            * (1.0
                + inverse
                + inverse_squared
                    * (1.0 / 2.0
                        - inverse_squared
                            * (1.0 / 6.0
                                - inverse_squared
                                    * (1.0 / 6.0 - inverse_squared * (3.0 / 20.0)))))
}

/// `logmdigamma(x) = log(x) - digamma(x)`.
///
/// limma reaches `statmod` for this; the identity is what that function
/// computes, confirmed against it to 2.6e-13 in the oracle.
pub fn logmdigamma(x: f64) -> f64 {
    x.ln() - digamma(x)
}

/// Inverse of [`trigamma`], matching `statmod::trigammaInverse`.
///
/// The two asymptotic branches matter: `fitFDist()` calls this with very small
/// arguments when the residual variances are tightly clustered, and the Newton
/// iteration started at `0.5 + 1/x` would otherwise converge slowly or not at
/// all.
pub fn trigamma_inverse(x: f64) -> f64 {
    if x <= 0.0 {
        return f64::NAN;
    }
    // Asymptotes, exactly as statmod switches them.
    if x > 1e7 {
        return 1.0 / x.sqrt();
    }
    if x < 1e-6 {
        return 1.0 / x;
    }

    let mut y = 0.5 + 1.0 / x;
    for _ in 0..100 {
        let first = trigamma(y);
        let second = tetragamma(y);
        if second == 0.0 {
            break;
        }
        // Newton on f(y) = trigamma(y) - x, written the way statmod does it.
        let step = first * (1.0 - first / x) / second;
        y += step;
        if step.abs() < 1e-12 * y.abs() {
            break;
        }
    }
    y
}

/// Continued fraction for the incomplete beta function (Lentz's algorithm).
///
/// Ported from Numerical Recipes' `betacf`, which is what R's `pbeta` uses in
/// spirit; the number of iterations and the `EPS`/`FPMIN` guards match.
fn betacf(a: f64, b: f64, x: f64) -> f64 {
    const MAX_ITERATIONS: usize = 300;
    const EPSILON: f64 = 3.0e-16;
    const TINY: f64 = 1.0e-300;

    let qab = a + b;
    let qap = a + 1.0;
    let qam = a - 1.0;
    let mut c = 1.0;
    let mut d = 1.0 - qab * x / qap;
    if d.abs() < TINY {
        d = TINY;
    }
    d = 1.0 / d;
    let mut h = d;

    for m in 1..=MAX_ITERATIONS {
        let m = m as f64;
        let m2 = 2.0 * m;
        // Even step.
        let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if d.abs() < TINY {
            d = TINY;
        }
        c = 1.0 + aa / c;
        if c.abs() < TINY {
            c = TINY;
        }
        d = 1.0 / d;
        h *= d * c;
        // Odd step.
        let aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if d.abs() < TINY {
            d = TINY;
        }
        c = 1.0 + aa / c;
        if c.abs() < TINY {
            c = TINY;
        }
        d = 1.0 / d;
        let delta = d * c;
        h *= delta;
        if (delta - 1.0).abs() < EPSILON {
            break;
        }
    }
    h
}

/// Regularised incomplete beta `I_x(a, b)` for `a, b > 0` and `x` in `[0, 1]`.
pub fn regularized_incomplete_beta(a: f64, b: f64, x: f64) -> f64 {
    if !(0.0..=1.0).contains(&x) {
        return f64::NAN;
    }
    if x == 0.0 {
        return 0.0;
    }
    if x == 1.0 {
        return 1.0;
    }

    // The continued fraction converges quickly only for x < (a+1)/(a+b+2), so
    // the symmetry `I_x(a,b) = 1 - I_{1-x}(b,a)` covers the other half.
    let front = (ln_gamma(a + b) - ln_gamma(a) - ln_gamma(b) + a * x.ln() + b * (1.0 - x).ln()).exp();
    if x < (a + 1.0) / (a + b + 2.0) {
        front * betacf(a, b, x) / a
    } else {
        1.0 - front * betacf(b, a, 1.0 - x) / b
    }
}

/// Two-sided t tail `2 * P(T > |t|)` with `df` degrees of freedom.
///
/// R computes `2 * pt(-abs(t), df)`. Rather than implement a t distribution,
/// this uses the identity that makes the tail a regularised incomplete beta:
///
/// ```text
/// 2 * P(T > |t|) = I_x(df/2, 1/2),   x = df / (df + t^2)
/// ```
///
/// which is exact for real `df` and avoids needing a separate t-CDF.
pub fn t_two_sided_p(t: f64, df: f64) -> f64 {
    if df <= 0.0 || !df.is_finite() || !t.is_finite() {
        return f64::NAN;
    }
    let x = df / (df + t * t);
    regularized_incomplete_beta(df / 2.0, 0.5, x)
}

#[cfg(test)]
// The R-derived expected values below are quoted at full precision so their
// provenance stays visible; see the same note in the `diffpeak` test module.
#[allow(clippy::excessive_precision)]
mod tests {
    use super::*;

    /// Euler-Mascheroni constant, used to check `digamma` against its exact
    /// value at 1 and 1/2.
    const EULER_MASCHERONI: f64 = 0.577_215_664_901_532_9;

    /// Relative comparison, since these functions span many magnitudes.
    ///
    /// The tolerances used below are set from `SERIES_SHIFT`'s measured
    /// accuracy, not the other way around.
    fn close(actual: f64, expected: f64, tolerance: f64) -> bool {
        if actual == expected {
            return true;
        }
        let scale = actual.abs().max(expected.abs()).max(1e-300);
        (actual - expected).abs() <= tolerance * scale
    }

    #[test]
    fn digamma_matches_known_values() {
        // psi(1) = -gamma, psi(0.5) = -gamma - 2ln2.
        assert!(close(digamma(1.0), -EULER_MASCHERONI, 1e-15));
        let expected_half = -EULER_MASCHERONI - 2.0 * 2.0f64.ln();
        assert!(close(digamma(0.5), expected_half, 1e-15));
    }

    #[test]
    fn trigamma_matches_known_values() {
        // psi'(1) = pi^2/6, psi'(0.5) = pi^2/2.
        assert!(close(trigamma(1.0), std::f64::consts::PI.powi(2) / 6.0, 1e-15));
        assert!(close(trigamma(0.5), std::f64::consts::PI.powi(2) / 2.0, 1e-15));
    }

    #[test]
    fn tetragamma_is_the_derivative_of_trigamma() {
        // Central difference, which is the cheapest independent check.
        for x in [0.7f64, 1.5, 3.0, 8.0, 40.0] {
            let h = 1e-5 * x.max(1.0);
            let numeric = (trigamma(x + h) - trigamma(x - h)) / (2.0 * h);
            // The central difference is only second-order accurate, so the
            // tolerance reflects the scheme rather than the function.
            assert!(
                close(numeric, tetragamma(x), 1e-5),
                "tetragamma({x}): numeric {numeric} vs {:.17}",
                tetragamma(x)
            );
        }
    }

    #[test]
    fn logmdigamma_is_log_minus_digamma() {
        for x in [0.25, 1.0, 2.5, 10.0] {
            assert!(close(logmdigamma(x), x.ln() - digamma(x), 1e-15));
        }
    }

    #[test]
    fn trigamma_inverse_inverts_trigamma() {
        // The round trip is what limma relies on: fitFDist turns an
        // "evar" into a df via this function.
        for x in [1e-3, 0.05, 0.2, 1.0, 5.0, 50.0] {
            let y = trigamma_inverse(x);
            assert!(
                close(trigamma(y), x, 1e-11),
                "trigamma_inverse({x}) = {y}, trigamma = {}",
                trigamma(y)
            );
        }
    }

    #[test]
    fn trigamma_inverse_matches_its_asymptotes() {
        // Both branches are switch points in statmod, so check either side.
        assert!(close(trigamma_inverse(1e-8), 1e8, 1e-12));
        assert!(close(trigamma_inverse(1e10), 1e-5, 1e-12));
    }

    #[test]
    fn incomplete_beta_satisfies_its_bounds() {
        assert_eq!(regularized_incomplete_beta(2.0, 3.0, 0.0), 0.0);
        assert_eq!(regularized_incomplete_beta(2.0, 3.0, 1.0), 1.0);
        // I_x(a,b) is increasing in x.
        let mut previous = 0.0;
        for step in 0..=20 {
            let x = step as f64 / 20.0;
            let value = regularized_incomplete_beta(2.0, 3.0, x);
            assert!(value >= previous - 1e-15, "not increasing at x = {x}");
            previous = value;
        }
    }

    #[test]
    fn incomplete_beta_uses_both_symmetry_branches() {
        // x near 1 exercises the `1 - I_{1-x}(b,a)` branch, and the two must
        // agree with the direct series at the switch point.
        for a in [0.5, 1.0, 2.0, 7.0] {
            for b in [0.5, 1.0, 3.0, 12.5] {
                let switch = (a + 1.0) / (a + b + 2.0);
                let below = regularized_incomplete_beta(a, b, switch - 1e-9);
                let above = regularized_incomplete_beta(a, b, switch + 1e-9);
                assert!(
                    close(below, above, 1e-6),
                    "discontinuity at the switch for a={a}, b={b}: {below} vs {above}"
                );
            }
        }
    }

    #[test]
    fn t_tail_is_symmetric_and_bounded() {
        assert!(close(t_two_sided_p(0.0, 5.0), 1.0, 1e-14));
        assert!(close(t_two_sided_p(2.0, 5.0), t_two_sided_p(-2.0, 5.0), 1e-15));
        // Large t drives the tail to zero.
        assert!(t_two_sided_p(1e6, 10.0) < 1e-20);
        // Known value: 2*P(T>1) with df=1 is 0.5.
        assert!(close(t_two_sided_p(1.0, 1.0), 0.5, 1e-14));
        // A hand-checkable pair with df = 2, where the tail has the closed form
        // 1 - t/sqrt(t^2 + 2). R's `2 * pt(-2, 2)` is 0.1835034190722739.
        assert!(close(t_two_sided_p(2.0, 2.0), 0.18350341907227391, 1e-14));
        // And df = 4, where 2*pt(-2, 4) = 0.1161165235168153.
        assert!(close(t_two_sided_p(2.0, 4.0), 0.1161165235168153, 1e-14));
    }

    #[test]
    fn degenerate_inputs_do_not_panic() {
        assert!(t_two_sided_p(f64::NAN, 5.0).is_nan());
        assert!(t_two_sided_p(1.0, 0.0).is_nan());
        assert!(t_two_sided_p(1.0, f64::INFINITY).is_nan());
        assert!(trigamma_inverse(0.0).is_nan());
        assert!(regularized_incomplete_beta(1.0, 1.0, 2.0).is_nan());
    }
}
