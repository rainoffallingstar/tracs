//! Validates the [`tracs::plot::layout`] port against a captured R oracle.
//!
//! `track_plot()` decides panel order and relative heights in `.make_layout()`,
//! whose `layout_ord` handling is easy to get subtly wrong: `s` (scale) is
//! injected with `c(lord, "s")` *before* intersecting with the available keys,
//! so an empty `layout_ord` puts the scale panel first rather than last.
//!
//! The oracle in `testdata/layout_r_oracle.tsv` was produced by R 4.6.0:
//!
//! ```r
//! library(data.table)
//! make_layout <- function(ntracks, ntracks_h = 3, cytoband = TRUE, cytoband_h = 1,
//!                         genemodel = TRUE, genemodel_h = 1, chrHMM = TRUE, chrHMM_h = 1,
//!                         loci = TRUE, loci_h = 2, scale_track_height = 1, lord = NULL) {
//!   lo_h_ord <- list("p" = loci_h, "b" = rep(ntracks_h, ntracks), "h" = chrHMM_h,
//!                    "g" = genemodel_h, "c" = cytoband_h, "s" = scale_track_height)
//!   keep <- c(if (loci) "p", "b", if (chrHMM) "h", if (genemodel) "g", "s", if (cytoband) "c")
//!   keep <- keep[!is.na(keep)]
//!   lo_h_ord <- lo_h_ord[keep]
//!   lord <- c(lord, "s")
//!   lord <- c(intersect(lord, names(lo_h_ord)), setdiff(names(lo_h_ord), lord))
//!   lo_heights <- unlist(lo_h_ord[lord], use.names = FALSE)
//!   # ... (the data.table block computing `data`, then:)
//!   lord_expanded <- character(0)
//!   for (k in lord) lord_expanded <- c(lord_expanded, rep(k, if (k == "b") ntracks else 1))
//!   paste(paste(lord_expanded[data], lo_heights, sep = ":", collapse = ","))
//! }
//! ```
//!
//! Regenerate with that script when new cases are needed.

use std::fs;
use std::path::PathBuf;

use anyhow::{anyhow, Context, Result};

/// `layout.rs` lives in the binary crate, so it is included directly.
#[path = "../src/plot/layout.rs"]
mod layout;

use layout::{make_layout, LayoutRequest, TrackKind};

fn oracle_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("testdata")
        .join("layout_r_oracle.tsv")
}

fn parse_flag(field: &str) -> Result<bool> {
    match field.trim() {
        "0" => Ok(false),
        "1" => Ok(true),
        other => Err(anyhow!("expected 0/1 flag, got {other:?}")),
    }
}

#[test]
fn layout_matches_r_oracle() -> Result<()> {
    let path = oracle_path();
    let text = fs::read_to_string(&path).with_context(|| format!("read oracle: {path:?}"))?;

    let mut checked = 0usize;
    for (line_number, line) in text.lines().enumerate() {
        let line = line.trim_end();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 13 {
            return Err(anyhow!(
                "malformed oracle row {} (expected 13 tab-separated fields, got {}): {line:?}",
                line_number + 1,
                fields.len()
            ));
        }
        let context = || format!("oracle line {}", line_number + 1);

        let bigwig_count: usize = fields[0].trim().parse().with_context(context)?;
        let layout_ord: Vec<TrackKind> = fields[1]
            .chars()
            .filter_map(TrackKind::from_key)
            .collect();
        let has_peaks = parse_flag(fields[2])?;
        let has_chromhmm = parse_flag(fields[3])?;
        let has_gene = parse_flag(fields[4])?;
        let has_cytoband = parse_flag(fields[5])?;
        let bigwig_height: f64 = fields[6].trim().parse().with_context(context)?;
        let peaks_height: f64 = fields[7].trim().parse().with_context(context)?;
        let gene_height: f64 = fields[8].trim().parse().with_context(context)?;
        let scale_height: f64 = fields[9].trim().parse().with_context(context)?;
        let chromhmm_height: f64 = fields[10].trim().parse().with_context(context)?;
        let cytoband_height: f64 = fields[11].trim().parse().with_context(context)?;

        let request = LayoutRequest {
            bigwig_height,
            peaks_height,
            gene_height,
            scale_height,
            chromhmm_height,
            cytoband_height,
            bigwig_count,
            has_peaks,
            has_chromhmm,
            has_gene,
            has_cytoband,
            layout_ord,
        };
        let resolved = make_layout(&request);
        let observed: Vec<String> = resolved
            .panels
            .iter()
            .map(|panel| {
                // Compare heights as integers: the oracle uses whole numbers and
                // R's printing would otherwise hide genuine differences.
                format!("{}:{}", panel.kind.key(), panel.height.round() as i64)
            })
            .collect();
        let observed = observed.join(",");

        if observed != fields[12] {
            return Err(anyhow!(
                "layout mismatch at {}\n  got:  {observed}\n  R:    {}",
                context(),
                fields[12]
            ));
        }
        checked += 1;
    }

    if checked == 0 {
        return Err(anyhow!("oracle contained no cases: {path:?}"));
    }
    eprintln!("validated {checked} layout configurations against the R oracle");

    Ok(())
}
