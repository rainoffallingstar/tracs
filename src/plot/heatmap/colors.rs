//! Sequential colour ramps for heatmaps.
//!
//! `profile_heatmap()` accepts any of grDevices' ~110 `hcl.pals(type =
//! "sequential")` palettes and interpolates them in HCL space. Reproducing that
//! whole catalogue (and the HCL interpolation maths) is out of proportion for a
//! heatmap, so this module ships the palettes that are actually useful here as
//! pre-sampled anchor colours, interpolated linearly in sRGB between them.
//!
//! Each ramp is stored **light to dark**, matching R's orientation after it
//! reverses `hcl.colors()`: low signal is pale, high signal is saturated.

/// A named sequential palette.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Palette {
    /// Name as accepted by `--col-pal`.
    pub name: &'static str,
    /// Anchor colours from lightest to darkest.
    pub anchors: &'static [&'static str],
}

/// Palettes available to `--col-pal`.
///
/// `Blues` is first because it is `profile_heatmap()`'s default.
pub const PALETTES: &[Palette] = &[
    // Sampled from R's hcl.colors(9, "Blues") then reversed.
    Palette {
        name: "Blues",
        anchors: &[
            "#F4FAFE", "#DEEEF7", "#C1DBEC", "#A1C4E0", "#7FABD3", "#5C90C6", "#3573B9",
            "#305596", "#273871",
        ],
    },
    // A perceptually smoother blue ramp, also light-to-dark.
    Palette {
        name: "Viridis",
        anchors: &[
            "#FDE725", "#DCE319", "#B8DE29", "#95D840", "#73D055", "#55C667", "#3CBB75",
            "#2D708E", "#440154",
        ],
    },
    Palette {
        name: "Greys",
        anchors: &[
            "#FFFFFF", "#F0F0F0", "#D9D9D9", "#BDBDBD", "#969696", "#737373", "#525252",
            "#252525", "#000000",
        ],
    },
    // R's single-hue reds; useful for marking enrichment specifically.
    Palette {
        name: "Reds",
        anchors: &[
            "#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D",
            "#A50F15", "#67000D",
        ],
    },
];

/// Default palette name, matching `profile_heatmap()`'s `col_pal`.
#[allow(dead_code)] // referenced by tests; the CLI default lives in clap's arg spec
pub const DEFAULT_PALETTE: &str = "Blues";

/// Resolves a palette by name, case-insensitively.
///
/// Returns `None` for an unknown name so callers can report the valid options
/// rather than silently falling back.
pub fn find_palette(name: &str) -> Option<&'static Palette> {
    let wanted = name.trim();
    PALETTES
        .iter()
        .find(|palette| palette.name.eq_ignore_ascii_case(wanted))
}

/// Names of every available palette, for error messages and `--help`.
pub fn palette_names() -> Vec<&'static str> {
    PALETTES.iter().map(|palette| palette.name).collect()
}

/// Parses a `#rrggbb` colour into 8-bit components.
fn parse_hex(color: &str) -> Option<(u8, u8, u8)> {
    let hex = color.strip_prefix('#')?;
    if hex.len() != 6 {
        return None;
    }
    let channel = |range: std::ops::Range<usize>| u8::from_str_radix(&hex[range], 16).ok();
    Some((channel(0..2)?, channel(2..4)?, channel(4..6)?))
}

/// Linearly interpolates between two hex colours.
fn mix(from: &str, to: &str, fraction: f64) -> String {
    match (parse_hex(from), parse_hex(to)) {
        (Some((r0, g0, b0)), Some((r1, g1, b1))) => {
            let blend = |a: u8, b: u8| -> u8 {
                let value = a as f64 + (b as f64 - a as f64) * fraction.clamp(0.0, 1.0);
                value.round().clamp(0.0, 255.0) as u8
            };
            format!(
                "#{:02X}{:02X}{:02X}",
                blend(r0, r1),
                blend(g0, g1),
                blend(b0, b1)
            )
        }
        // A malformed anchor should not abort rendering; fall back to the start.
        _ => from.to_string(),
    }
}

/// Expands a palette into `steps` discrete colours, light to dark.
///
/// `profile_heatmap()` builds a 255-entry ramp. The count is caller-controlled
/// here so tests can use a handful of colours and the renderer can match R.
pub fn resolve_ramp(palette: &Palette, steps: usize) -> Vec<String> {
    let anchors = palette.anchors;
    // A zero-step ramp means "no colour scale at all", which is distinct from
    // "one colour"; checking this first keeps the two cases from collapsing.
    if steps == 0 || anchors.is_empty() {
        return Vec::new();
    }
    if anchors.len() == 1 || steps == 1 {
        return vec![anchors[0].to_string()];
    }

    (0..steps)
        .map(|step| {
            let position = step as f64 / (steps - 1) as f64 * (anchors.len() - 1) as f64;
            let lower = position.floor() as usize;
            let upper = (lower + 1).min(anchors.len() - 1);
            let fraction = position - lower as f64;
            if lower == upper {
                anchors[lower].to_string()
            } else {
                mix(anchors[lower], anchors[upper], fraction)
            }
        })
        .collect()
}

/// Reverses a ramp, matching `profile_heatmap()`'s `revpal = TRUE`.
pub fn reverse_ramp(ramp: &mut [String]) {
    ramp.reverse();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_palette_exists_and_is_blues() {
        let palette = find_palette(DEFAULT_PALETTE).expect("default palette must exist");
        assert_eq!(palette.name, "Blues");
    }

    #[test]
    fn palette_lookup_is_case_insensitive() {
        assert!(find_palette("blues").is_some());
        assert!(find_palette("  BLUES ").is_some());
        assert!(find_palette("nope").is_none());
    }

    #[test]
    fn ramp_spans_light_to_dark_and_matches_the_anchor_ends() {
        let palette = find_palette("Blues").unwrap();
        let ramp = resolve_ramp(palette, 255);
        assert_eq!(ramp.len(), 255);
        // Ends must be exactly the anchor colours so the extremes are not drifted
        // by interpolation.
        assert_eq!(ramp.first().unwrap(), "#F4FAFE", "low end should be pale");
        assert_eq!(ramp.last().unwrap(), "#273871", "high end should be saturated");

        // Monotonically darkening: compare summed channels across the ramp.
        let brightness = |hex: &str| -> u32 {
            let (r, g, b) = parse_hex(hex).unwrap();
            r as u32 + g as u32 + b as u32
        };
        let first_quarter = brightness(&ramp[0]);
        let last_quarter = brightness(&ramp[254]);
        assert!(
            last_quarter < first_quarter,
            "ramp should darken: {first_quarter} -> {last_quarter}"
        );
    }

    #[test]
    fn ramp_interpolates_between_anchors() {
        let palette = Palette {
            name: "test",
            anchors: &["#000000", "#FFFFFF"],
        };
        let ramp = resolve_ramp(&palette, 3);
        assert_eq!(ramp, vec!["#000000", "#808080", "#FFFFFF"]);
    }

    #[test]
    fn degenerate_step_counts_are_handled() {
        let palette = find_palette("Blues").unwrap();
        // Asking for one colour yields the light end rather than dividing by zero.
        assert_eq!(resolve_ramp(palette, 1), vec!["#F4FAFE".to_string()]);
        assert!(resolve_ramp(palette, 0).is_empty());
    }

    #[test]
    fn reversing_flips_the_orientation() {
        let palette = find_palette("Blues").unwrap();
        let mut ramp = resolve_ramp(palette, 5);
        let light_first = ramp.clone();
        reverse_ramp(&mut ramp);
        assert_eq!(ramp[0], light_first[4]);
        assert_eq!(ramp[4], light_first[0]);
    }

    #[test]
    fn all_shipped_ramps_are_light_to_dark() {
        // Guards against a palette being added with its anchors reversed, which
        // would silently invert the heatmap's meaning.
        for palette in PALETTES {
            let ramp = resolve_ramp(palette, 32);
            let (r0, g0, b0) = parse_hex(&ramp[0]).unwrap();
            let (r1, g1, b1) = parse_hex(ramp.last().unwrap()).unwrap();
            let start = r0 as u32 + g0 as u32 + b0 as u32;
            let end = r1 as u32 + g1 as u32 + b1 as u32;
            assert!(
                end < start,
                "palette {} is not light-to-dark: {} -> {}",
                palette.name,
                ramp[0],
                ramp.last().unwrap()
            );
        }
    }
}
