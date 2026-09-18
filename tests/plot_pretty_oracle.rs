//! Validates the [`tracs::plot::pretty`] port against a captured R oracle.
//!
//! `track_plot()` positions axis ticks with R's `pretty()`. `ggplot-rs` uses the
//! extended-Wilkinson algorithm instead, so the port exists to reproduce R's
//! behaviour exactly; this test is what keeps the port honest.
//!
//! The oracle in `testdata/pretty_r_oracle.tsv` was produced by R 4.6.0:
//!
//! ```r
//! set.seed(11)
//! fm <- function(v) paste(format(v, scientific = FALSE, trim = TRUE, digits = 15), collapse = ",")
//! cases <- list(c(0, 0), c(1, 1), c(-3, -3), c(0.5, 0.5),
//!   c(158145820, 158156686), c(0, 51.93), c(0, 8.8), c(1, 10), c(-5, 5), c(0, 1),
//!   c(1e5, 1e6), c(0, 1e6), c(1e6, 2e6), c(-1, -0.001), c(0.001, 0.002),
//!   c(123.456, 789.012), c(-98765.4, -1234.5), c(0, 1e-8), c(1e9, 2e9),
//!   c(46913486, 46964325), c(7565097, 7590856), c(128747680, 128753674))
//! for (c_ in cases) cat(sprintf("%.17g\t%.17g\t%s\n", c_[1], c_[2], fm(pretty(c_))))
//! for (i in 1:40) {
//!   a <- runif(1, -1e7, 1e7); b <- a + 10^runif(1, -6, 7)
//!   cat(sprintf("%.17g\t%.17g\t%s\n", a, b, fm(pretty(c(a, b)))))
//! }
//! ```
//!
//! Regenerate with the script above when the oracle needs extending.

use std::fs;
use std::path::PathBuf;

use anyhow::{anyhow, Context, Result};

/// The `pretty` module lives in the binary crate, so the test includes it directly
/// rather than through a library path.
#[path = "../src/plot/pretty.rs"]
mod pretty;

fn oracle_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("testdata")
        .join("pretty_r_oracle.tsv")
}

/// Relative tolerance used when comparing tick values.
///
/// R stores these as IEEE doubles, and the port follows the same operations, so
/// the only expected divergence is float formatting in the oracle file (which is
/// written with 15 significant digits).
const REL_TOLERANCE: f64 = 1e-12;

fn approx_eq(a: f64, b: f64) -> bool {
    if a == b {
        return true;
    }
    let scale = a.abs().max(b.abs()).max(1.0);
    (a - b).abs() <= REL_TOLERANCE * scale
}

#[test]
fn port_matches_r_oracle() -> Result<()> {
    let path = oracle_path();
    let text = fs::read_to_string(&path).with_context(|| format!("read oracle: {path:?}"))?;

    let mut checked = 0usize;
    for (line_number, line) in text.lines().enumerate() {
        let line = line.trim_end();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 3 {
            return Err(anyhow!(
                "malformed oracle row {} (expected 3 tab-separated fields): {line:?}",
                line_number + 1
            ));
        }

        let lo: f64 = fields[0].parse().with_context(|| {
            format!("parse lo at line {}: {:?}", line_number + 1, fields[0])
        })?;
        let up: f64 = fields[1].parse().with_context(|| {
            format!("parse up at line {}: {:?}", line_number + 1, fields[1])
        })?;
        let expected: Vec<f64> = if fields[2].is_empty() {
            Vec::new()
        } else {
            fields[2]
                .split(',')
                .map(|value| {
                    value.trim().parse::<f64>().with_context(|| {
                        format!("parse tick {value:?} at line {}", line_number + 1)
                    })
                })
                .collect::<Result<Vec<f64>>>()?
        };

        let got = pretty::pretty_range(lo, up, 5);

        if got.values.len() != expected.len() {
            return Err(anyhow!(
                "tick count mismatch for [{lo}, {up}]: got {} ({:?}), R has {} ({:?})",
                got.values.len(),
                got.values,
                expected.len(),
                expected
            ));
        }
        for (index, (actual, want)) in got.values.iter().zip(expected.iter()).enumerate() {
            if !approx_eq(*actual, *want) {
                return Err(anyhow!(
                    "tick {index} mismatch for [{lo}, {up}]: got {actual}, R has {want}\n  got:  {:?}\n  R:    {:?}",
                    got.values,
                    expected
                ));
            }
        }
        checked += 1;
    }

    if checked == 0 {
        return Err(anyhow!("oracle contained no cases: {path:?}"));
    }
    // Guard against the oracle being silently truncated.
    if checked < 50 {
        return Err(anyhow!(
            "oracle only had {checked} cases; expected at least 50"
        ));
    }
    eprintln!("validated {checked} intervals against the R oracle");

    Ok(())
}
