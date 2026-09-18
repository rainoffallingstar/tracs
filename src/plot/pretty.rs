//! R's `pretty()` algorithm, ported from R's `src/appl/pretty.c` (`R_pretty`) and
//! `base::pretty.default`.
//!
//! `track_plot()` places axis ticks and gene-direction arrows with `pretty()`, and
//! `ggplot-rs` uses the extended-Wilkinson algorithm instead. The two disagree, so
//! ticks would drift; this module reproduces R's behaviour exactly and the results
//! are fed to `ScaleContinuous::with_breaks()` to override the engine's own choice.
//!
//! The port follows the C source closely, including its floating-point subtleties
//! (`rounding_eps`, the `i_small` branch, and the `high.u.bias`/`u5.bias` unit
//! selection), because those determine where the ticks land.

/// `rounding_eps` from `pretty.c` (was `1e-5` before R 0.65).
const ROUNDING_EPS: f64 = 1e-10;

/// Defaults from `base::pretty.default`.
const DEFAULT_N: i32 = 5;
const DEFAULT_SHRINK_SML: f64 = 0.75;
const DEFAULT_HIGH_U_BIAS: f64 = 1.5;
/// `u5.bias = 0.5 + 1.5 * high.u.bias`
const DEFAULT_U5_BIAS: f64 = 0.5 + 1.5 * DEFAULT_HIGH_U_BIAS;
/// `f.min = 2^-20`
const DEFAULT_F_MIN: f64 = 9.5367431640625e-7;

/// Result of the interval computation: the tick sequence plus the step size.
#[derive(Clone, Debug, PartialEq)]
pub struct Pretty {
    /// Tick positions, equivalent to R's returned numeric vector.
    pub values: Vec<f64>,
    /// The chosen unit (tick spacing). R returns this as `z$unit`.
    pub unit: f64,
    /// The number of intervals, R's `z$n`.
    pub intervals: i32,
}

impl Pretty {
    /// Tick positions as a `Vec<f64>`, for direct use as scale breaks.
    pub fn breaks(&self) -> Vec<f64> {
        self.values.clone()
    }

    /// First tick, i.e. R's `z$l` (the lower bound after `bounds = TRUE`).
    pub fn low(&self) -> f64 {
        self.values.first().copied().unwrap_or(0.0)
    }

    /// Last tick, i.e. R's `z$u`.
    pub fn high(&self) -> f64 {
        self.values.last().copied().unwrap_or(0.0)
    }
}

/// Port of `pretty.default(x, n = 5L, ...)`; the interval is `[min(x), max(x)]`.
///
/// Non-finite inputs are dropped, matching R. An empty (or all-non-finite) input
/// yields no ticks, and a degenerate `lo == up` is handled the way R does.
pub fn pretty(x: &[f64]) -> Pretty {
    pretty_with(x, DEFAULT_N)
}

/// `pretty.default()` with an explicit `n`.
pub fn pretty_with(x: &[f64], n: i32) -> Pretty {
    let finite: Vec<f64> = x.iter().copied().filter(|v| v.is_finite()).collect();
    if finite.is_empty() {
        return Pretty {
            values: Vec::new(),
            unit: 0.0,
            intervals: 0,
        };
    }
    let lo = finite.iter().copied().fold(f64::INFINITY, f64::min);
    let up = finite.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    pretty_bounds(lo, up, n)
}

/// `pretty(c(lo, up), n)`.
pub fn pretty_range(lo: f64, up: f64, n: i32) -> Pretty {
    pretty_bounds(lo, up, n)
}

/// Core of `pretty.default`: `R_pretty()` followed by `seq.int(l, u, length.out = n + 1)`.
fn pretty_bounds(lo_in: f64, up_in: f64, n: i32) -> Pretty {
    // `min.n = n %/% 3L`
    let min_n = n.div_euclid(3);

    let (mut lo, mut up, mut ndiv) = (lo_in, up_in, n);
    let unit = r_pretty(
        &mut lo,
        &mut up,
        &mut ndiv,
        min_n,
        DEFAULT_SHRINK_SML,
        [DEFAULT_HIGH_U_BIAS, DEFAULT_U5_BIAS, DEFAULT_F_MIN],
        // `eps.correct = 0L` in pretty.default, and `bounds = TRUE`.
        0,
        true,
    );

    // R builds the result as `seq.int(z$l, z$u, length.out = n + 1L)`; with
    // `bounds = TRUE` the bounds are already ns*unit / nu*unit, so an evenly
    // spaced sequence over ndiv intervals reproduces it.
    let count = if ndiv > 0 { ndiv + 1 } else { 0 };
    let mut values: Vec<f64> = if count > 0 {
        (0..count).map(|i| lo + (up - lo) * (i as f64) / (ndiv as f64)).collect()
    } else {
        Vec::new()
    };

    // pretty.default snaps tiny values to exact zero:
    //   delta <- diff(range(l, u) / n)
    //   if (any(abs(s) < 1e-14 * delta)) s[...] <- 0
    if ndiv > 0 {
        let delta = (up - lo) / (ndiv as f64);
        let threshold = 1e-14 * delta;
        for value in values.iter_mut() {
            if value.abs() < threshold {
                *value = 0.0;
            }
        }
    }

    Pretty {
        values,
        unit,
        intervals: ndiv,
    }
}

/// Direct port of `R_pretty()` from `src/appl/pretty.c`.
///
/// With `return_bounds = true` (R's `bounds = TRUE`) `*lo`/`*up` are widened to
/// the tick bounds; with `false` they become the integer tick counts `ns`/`nu`.
/// `ndiv` is always overwritten with the interval count.
#[allow(clippy::too_many_arguments)]
fn r_pretty(
    lo: &mut f64,
    up: &mut f64,
    ndiv: &mut i32,
    min_n: i32,
    shrink_sml: f64,
    high_u_fact: [f64; 3],
    eps_correction: i32,
    return_bounds: bool,
) -> f64 {
    let high_u_bias = high_u_fact[0];
    let u5_bias = high_u_fact[1];
    let f_min = high_u_fact[2];

    let lo_input = *lo;
    let up_input = *up;
    let dx = up_input - lo_input;

    // cell := "scale"
    let cell_initial;
    let i_small;
    if dx == 0.0 && up_input == 0.0 {
        // up == lo == 0
        cell_initial = 1.0;
        i_small = true;
    } else {
        let cell = lo_input.abs().max(up_input.abs());
        // U = upper bound on cell/unit
        let mut u = 1.0
            + if u5_bias >= 1.5 * high_u_bias + 0.5 {
                1.0 / (1.0 + high_u_bias)
            } else {
                1.5 / (1.0 + u5_bias)
            };
        u *= (std::cmp::max(*ndiv, 1) as f64) * f64::EPSILON;
        // "added times 3, as several calculations here"
        i_small = dx < cell * u * 3.0;
        cell_initial = cell;
    }

    let cell = if i_small {
        let mut cell = cell_initial;
        if cell > 10.0 {
            cell = 9.0 + cell / 10.0;
        }
        cell *= shrink_sml;
        if min_n > 1 {
            cell /= min_n as f64;
        }
        cell
    } else {
        let mut cell = dx;
        if cell.is_finite() {
            if *ndiv > 1 {
                cell /= *ndiv as f64;
            }
        } else if *ndiv >= 2 {
            cell = up_input / (*ndiv as f64) - lo_input / (*ndiv as f64);
        }
        cell
    };

    // Clamp against subnormal underflow and overflow, as the C code does.
    const MAX_F: f64 = 1.25;
    let subsmall = if f_min * f64::MIN_POSITIVE == 0.0 {
        f64::MIN_POSITIVE
    } else {
        f_min * f64::MIN_POSITIVE
    };
    let cell = if cell < subsmall {
        subsmall
    } else if cell > f64::MAX / MAX_F {
        f64::MAX / MAX_F
    } else {
        cell
    };

    // base <= cell < 10*base
    let base = 10f64.powf(cell.log10().floor());

    // unit: one of {1,2,5,10} * base, favouring 5 over 2 when h5 > h.
    let mut unit = base;
    let mut u = 2.0 * base - cell;
    if u < high_u_bias * (cell - unit) {
        unit = 2.0 * base;
        u = 5.0 * base - cell;
        if u < u5_bias * (cell - unit) {
            unit = 5.0 * base;
            u = 10.0 * base - cell;
            if u < high_u_bias * (cell - unit) {
                unit = 10.0 * base;
            }
        }
    }

    let mut ns = (lo_input / unit + ROUNDING_EPS).floor();
    let mut nu = (up_input / unit - ROUNDING_EPS).ceil();

    if eps_correction != 0 && (eps_correction > 1 || !i_small) {
        let e = f64::EPSILON;
        let d_max = f64::MAX * (1.0 - (e * 0.5));
        // Move *lo to the left.
        if lo_input < 0.0 {
            *lo *= 1.0 + e;
        } else if lo_input > 0.0 {
            *lo *= 1.0 - e;
        } else {
            *lo = -unit.min(f64::MIN_POSITIVE);
        }
        // And *up to the right.
        if up_input < 0.0 {
            *up *= 1.0 - e;
        } else if up_input > 0.0 {
            if up_input < d_max {
                *up *= 1.0 + e;
            }
        } else {
            *up = unit.min(f64::MIN_POSITIVE);
        }
    }

    while ns * unit > *lo + ROUNDING_EPS * unit {
        ns -= 1.0;
    }
    while !(ns * unit).is_finite() {
        ns += 1.0;
    }
    while nu * unit < *up - ROUNDING_EPS * unit {
        nu += 1.0;
    }
    while !(nu * unit).is_finite() {
        nu -= 1.0;
    }

    let mut k = (0.5 + nu - ns) as i32;
    if k < min_n {
        // Ensure nu - ns == min_n.
        let k_diff = min_n - k;
        if lo_input == 0.0 && ns == 0.0 && up_input != 0.0 {
            nu += k_diff as f64;
        } else if up_input == 0.0 && nu == 0.0 && lo_input != 0.0 {
            ns -= k_diff as f64;
        } else if ns >= 0.0 {
            nu += (k_diff / 2) as f64;
            ns -= (k_diff / 2 + k_diff % 2) as f64;
        } else {
            ns -= (k_diff / 2) as f64;
            nu += (k_diff / 2 + k_diff % 2) as f64;
        }
        *ndiv = min_n;
    } else {
        *ndiv = k;
        k = 0;
    }
    let _ = k;

    if return_bounds {
        if ns * unit < *lo {
            *lo = ns * unit;
        }
        if nu * unit > *up {
            *up = nu * unit;
        }
    } else {
        *lo = ns;
        *up = nu;
    }

    unit
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Values captured from R 4.6.0's `pretty()` for the intervals `track_plot()`
    /// actually produces: genomic loci and bigWig max-signal ranges.
    #[test]
    fn matches_r_reference_intervals() {
        let cases: &[(f64, f64, &[f64])] = &[
            (
                158_145_820.0,
                158_156_686.0,
                &[
                    158_144_000.0,
                    158_146_000.0,
                    158_148_000.0,
                    158_150_000.0,
                    158_152_000.0,
                    158_154_000.0,
                    158_156_000.0,
                    158_158_000.0,
                ],
            ),
            (0.0, 51.93, &[0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0]),
            (0.0, 8.8, &[0.0, 2.0, 4.0, 6.0, 8.0, 10.0]),
            (1.0, 10.0, &[0.0, 2.0, 4.0, 6.0, 8.0, 10.0]),
            (-5.0, 5.0, &[-6.0, -4.0, -2.0, 0.0, 2.0, 4.0, 6.0]),
            (0.0, 1.0, &[0.0, 0.2, 0.4, 0.6, 0.8, 1.0]),
            (
                100_000.0,
                1_000_000.0,
                &[0.0, 200_000.0, 400_000.0, 600_000.0, 800_000.0, 1_000_000.0],
            ),
            (
                0.0,
                1_000_000.0,
                &[0.0, 200_000.0, 400_000.0, 600_000.0, 800_000.0, 1_000_000.0],
            ),
            (
                1_000_000.0,
                2_000_000.0,
                &[
                    1_000_000.0,
                    1_200_000.0,
                    1_400_000.0,
                    1_600_000.0,
                    1_800_000.0,
                    2_000_000.0,
                ],
            ),
        ];

        for (lo, up, expected) in cases {
            let got = pretty_range(*lo, *up, DEFAULT_N);
            assert_eq!(
                got.values.len(),
                expected.len(),
                "tick count for [{lo}, {up}]: got {:?}",
                got.values
            );
            for (index, (actual, want)) in got.values.iter().zip(expected.iter()).enumerate() {
                assert!(
                    (actual - want).abs() <= 1e-9 * want.abs().max(1.0),
                    "tick {index} for [{lo}, {up}]: got {actual}, want {want} (all: {:?})",
                    got.values
                );
            }
        }
    }

    #[test]
    fn handles_degenerate_and_empty_input() {
        assert!(pretty(&[]).values.is_empty());
        assert!(pretty(&[f64::NAN, f64::INFINITY]).values.is_empty());

        // Degenerate intervals are handled by R's min.n padding, which widens the
        // range rather than returning the point itself. Expectations below are the
        // literal R 4.6.0 outputs.
        let cases: &[(f64, &[f64])] = &[
            (0.0, &[-1.0, 0.0]),
            (2.0, &[0.0, 2.0]),
            (-1.0, &[-1.0, 0.0]),
            (5.0, &[0.0, 5.0]),
        ];
        for (value, expected) in cases {
            let got = pretty_range(*value, *value, DEFAULT_N);
            assert_eq!(
                got.values, *expected,
                "degenerate pretty({value}) should match R"
            );
        }
    }

    #[test]
    fn bounds_cover_the_input_range() {
        // The whole point of `bounds = TRUE`: ticks must span the requested range.
        for (lo, up) in [(0.3, 0.7), (-3.2, 7.4), (1234.5, 98765.4), (0.0, 51.93)] {
            let got = pretty_range(lo, up, DEFAULT_N);
            assert!(
                got.low() <= lo && got.high() >= up,
                "ticks {:?} do not cover [{lo}, {up}]",
                got.values
            );
        }
    }
}
