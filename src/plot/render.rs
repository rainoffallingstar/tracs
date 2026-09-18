//! Top-level track renderer: layout, panels, and PDF/SVG output.
//!
//! This is the replacement for `Rscript trackplot.R`. It reads the intermediates
//! written by `track-extract`, resolves the panel layout the same way
//! `.make_layout()` does, draws each track with its own y-scale, and converts the
//! composed SVG to PDF.

use std::path::Path;

use anyhow::{anyhow, Context, Result};

use crate::plot::fonts;
use crate::plot::io;
use crate::plot::layout::{make_layout, LayoutRequest, TrackKind};
use crate::plot::svg::{
    self, draw_gene_panel, draw_ideogram_panel, draw_peaks_panel, draw_scale_panel,
    draw_signal_panel, SvgWriter, XAxis, YAxis,
};

/// Output format for the renderer.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OutputFormat {
    Pdf,
    Svg,
}

impl OutputFormat {
    /// Infers the format from a path's extension, defaulting to PDF.
    pub fn from_path(path: &Path) -> Self {
        match path
            .extension()
            .and_then(|extension| extension.to_str())
            .map(|extension| extension.to_ascii_lowercase())
            .as_deref()
        {
            Some("svg") => OutputFormat::Svg,
            _ => OutputFormat::Pdf,
        }
    }
}

/// Everything the renderer needs, mirroring the `track_plot()` arguments that
/// affect drawing.
#[derive(Clone, Debug)]
pub struct RenderOptions {
    /// Total canvas width in points.
    pub width: f64,
    /// Total canvas height in points.
    pub height: f64,
    /// Track colours, cycled across samples when fewer are given.
    pub colors: Vec<String>,
    /// Show y-axis ticks instead of `[min-max]` annotations.
    pub show_axis: bool,
    /// Draw the ideogram panel.
    pub show_ideogram: bool,
    /// Draw the gene model panel.
    pub draw_gene_track: bool,
    /// Place track names to the left instead of as centred titles.
    pub track_names_to_left: bool,
    /// Explicit display names, overriding sample names.
    pub track_names: Option<Vec<String>>,
    /// Base font size in points.
    pub font_size: f64,
    /// Panel heights (relative), forwarded to the layout engine.
    pub bigwig_height: f64,
    pub peaks_height: f64,
    pub gene_height: f64,
    pub scale_height: f64,
    pub chromhmm_height: f64,
    pub cytoband_height: f64,
    /// Explicit y maxima, one per sample (or cycled).
    pub y_max: Option<Vec<f64>>,
    /// Auto-scale all bigWig tracks to a shared maximum.
    pub group_auto_scale: bool,
    /// User panel order (`layout_ord`).
    pub layout_ord: Vec<TrackKind>,
    /// Width reserved for the left margin (axis labels / track names).
    pub left_margin: f64,
    /// Right margin.
    pub right_margin: f64,
    /// Space above each panel for its centred title.
    pub panel_title_space: f64,
    /// Peaks sets, already loaded as `(name, intervals)`.
    pub peaks: Vec<(String, Vec<(u64, u64)>)>,
}

impl Default for RenderOptions {
    fn default() -> Self {
        Self {
            width: 864.0,
            height: 432.0,
            colors: Vec::new(),
            show_axis: false,
            show_ideogram: true,
            draw_gene_track: true,
            track_names_to_left: false,
            track_names: None,
            font_size: 10.0,
            bigwig_height: 3.0,
            peaks_height: 2.0,
            gene_height: 2.0,
            scale_height: 2.0,
            chromhmm_height: 1.0,
            cytoband_height: 2.0,
            y_max: None,
            group_auto_scale: false,
            layout_ord: Vec::new(),
            left_margin: 60.0,
            right_margin: 12.0,
            panel_title_space: 14.0,
            peaks: Vec::new(),
        }
    }
}

/// Inputs resolved from a `--work-dir`.
#[derive(Clone, Debug)]
pub struct RenderInputs {
    pub tracks: Vec<io::SampleTrack>,
    pub region: io::Region,
    pub transcripts: Vec<io::Transcript>,
    pub cytobands: Vec<io::Cytoband>,
}

/// Loads every intermediate file present in `work_dir`.
///
/// Missing optional files (gene models, cytobands) are treated as empty rather
/// than errors, matching how the extraction step omits them when the caller
/// passes `--no-gene-models` or `--no-cytoband`.
pub fn load_inputs(work_dir: &Path) -> Result<RenderInputs> {
    let tracks_path = work_dir.join("tracks.tsv");
    let meta_path = work_dir.join("meta.tsv");

    let tracks = io::read_tracks(&tracks_path)
        .with_context(|| format!("read tracks from {work_dir:?}"))?;
    let region = io::read_region(&meta_path)
        .with_context(|| format!("read region from {work_dir:?}"))?;

    let gene_models_path = work_dir.join("gene_models.tsv");
    let transcripts = if gene_models_path.exists() {
        io::read_transcripts(&gene_models_path)?
    } else {
        Vec::new()
    };

    let cytoband_path = work_dir.join("cytoband.tsv");
    let cytobands = if cytoband_path.exists() {
        io::read_cytobands(&cytoband_path, &region.chromosome)?
    } else {
        Vec::new()
    };

    Ok(RenderInputs {
        tracks,
        region,
        transcripts,
        cytobands,
    })
}

/// Resolves the y maximum for each bigWig track.
///
/// Reproduces `track_plot()`'s logic: explicit `y_max` wins, then
/// `groupAutoScale` uses one shared maximum, otherwise each track scales to its
/// own maximum. Values are rounded to two decimals as R does.
fn resolve_y_maxima(tracks: &[io::SampleTrack], options: &RenderOptions) -> Vec<f64> {
    if let Some(explicit) = options.y_max.as_ref().filter(|values| !values.is_empty()) {
        return (0..tracks.len())
            .map(|index| explicit[index % explicit.len()])
            .collect();
    }

    if options.group_auto_scale {
        let shared = tracks
            .iter()
            .map(io::SampleTrack::max_signal)
            .filter(|value| value.is_finite())
            .fold(f64::NEG_INFINITY, f64::max);
        // A track with no signal at all would otherwise produce a zero range.
        let shared = if shared.is_finite() && shared > 0.0 {
            shared
        } else {
            1.0
        };
        return vec![round_two(shared); tracks.len()];
    }

    tracks
        .iter()
        .map(|track| {
            let max = track.max_signal();
            if max.is_finite() && max > 0.0 {
                round_two(max)
            } else {
                1.0
            }
        })
        .collect()
}

/// Rounds to two decimals, matching `round(plot_height, digits = 2)`.
fn round_two(value: f64) -> f64 {
    (value * 100.0).round() / 100.0
}

/// Composes the SVG document for the given inputs.
pub fn render_svg(inputs: &RenderInputs, options: &RenderOptions) -> Result<String> {
    if inputs.tracks.is_empty() {
        return Err(anyhow!("no tracks to render"));
    }

    let layout_request = LayoutRequest {
        bigwig_height: options.bigwig_height,
        peaks_height: options.peaks_height,
        gene_height: options.gene_height,
        scale_height: options.scale_height,
        chromhmm_height: options.chromhmm_height,
        cytoband_height: options.cytoband_height,
        bigwig_count: inputs.tracks.len(),
        has_peaks: !options.peaks.is_empty(),
        // chromHMM panels are not produced by track-extract yet.
        has_chromhmm: false,
        has_gene: options.draw_gene_track && !inputs.transcripts.is_empty(),
        has_cytoband: options.show_ideogram && !inputs.cytobands.is_empty(),
        layout_ord: options.layout_ord.clone(),
    };
    let layout = make_layout(&layout_request);

    // Panels share the same padding so their data areas line up vertically.
    let plot_left = options.left_margin;
    let plot_right = options.width - options.right_margin;
    if plot_right <= plot_left {
        return Err(anyhow!(
            "canvas too narrow: left margin {} leaves no room for data (width {})",
            options.left_margin,
            options.width
        ));
    }

    // Reserve title space on every panel so heights stay proportional.
    let units_per_point = if layout.total_height > 0.0 {
        options.height / layout.total_height
    } else {
        0.0
    };

    let colors = svg::resolve_colors(&options.colors, inputs.tracks.len());
    let y_maxima = resolve_y_maxima(&inputs.tracks, options);

    let mut writer = SvgWriter::new(options.width, options.height);
    let mut cursor_y = 0.0f64;

    for panel in &layout.panels {
        let panel_height = (panel.height * units_per_point).max(1.0);
        let panel_top = cursor_y;
        cursor_y += panel_height;

        let axis = XAxis {
            plot_left,
            plot_right,
            data_start: inputs.region.start as f64,
            data_end: inputs.region.end as f64,
        };

        writer.panel(0.0, panel_top, options.width, panel_height, |draw| {
            // Every panel reserves the same title space, then uses what remains
            // for its data area. This keeps baselines aligned across tracks.
            let title_space = if options.track_names_to_left {
                0.0
            } else {
                options.panel_title_space
            };
            let plot_top = title_space;
            let plot_bottom = (panel_height - 2.0).max(plot_top + 1.0);

            match panel.kind {
                TrackKind::BigWig => {
                    let index = panel.bigwig_index.unwrap_or(0);
                    let Some(track) = inputs.tracks.get(index) else {
                        return;
                    };
                    let y_max = y_maxima.get(index).copied().unwrap_or(1.0);
                    let name = options
                        .track_names
                        .as_ref()
                        .and_then(|names| names.get(index))
                        .cloned()
                        .unwrap_or_else(|| track.sample.clone());
                    draw_signal_panel(
                        draw,
                        track,
                        axis,
                        YAxis {
                            plot_top,
                            plot_bottom,
                            y_max,
                        },
                        colors.get(index).map(String::as_str).unwrap_or(svg::DEFAULT_TRACK_COLOR),
                        options.show_axis,
                        &name,
                        options.track_names_to_left,
                        options.font_size,
                    );
                }
                TrackKind::Gene => {
                    draw_gene_panel(draw, &inputs.transcripts, axis, options.font_size);
                }
                TrackKind::Scale => {
                    draw_scale_panel(
                        draw,
                        axis,
                        &inputs.region.chromosome,
                        inputs.region.start,
                        inputs.region.end,
                        options.font_size,
                    );
                }
                TrackKind::Cytoband => {
                    draw_ideogram_panel(
                        draw,
                        &inputs.cytobands,
                        &inputs.region.chromosome,
                        inputs.region.start,
                        inputs.region.end,
                        options.font_size,
                    );
                }
                TrackKind::Peaks => {
                    draw_peaks_panel(draw, &options.peaks, axis, options.font_size);
                }
                TrackKind::ChromHmm => {
                    // chromHMM tracks are not produced by track-extract yet; the
                    // layout reserves no space for them (`has_chromhmm` is false),
                    // so this arm is unreachable in practice.
                }
            }
        });
    }

    Ok(writer.finish())
}

/// Converts an SVG document to PDF bytes.
///
/// The svg is laid out at 72 dpi so one SVG user unit equals one point, which
/// makes the resulting page dimensions match `options.width`/`height` exactly.
pub fn svg_to_pdf(svg_document: &str) -> Result<Vec<u8>> {
    let options = fonts::usvg_options();
    let tree = svg2pdf::usvg::Tree::from_str(svg_document, &options)
        .map_err(|error| anyhow!("parse composed SVG: {error:?}"))?;
    let pdf = svg2pdf::to_pdf(
        &tree,
        svg2pdf::ConversionOptions::default(),
        svg2pdf::PageOptions::default(),
    )
    .map_err(|error| anyhow!("convert SVG to PDF: {error:?}"))?;

    if pdf.len() < 5 || &pdf[..5] != b"%PDF-" {
        return Err(anyhow!(
            "SVG conversion produced an invalid PDF ({} bytes)",
            pdf.len()
        ));
    }
    Ok(pdf)
}

/// Renders `inputs` to `out_path`, choosing the format from the extension.
pub fn render_to_file(inputs: &RenderInputs, options: &RenderOptions, out_path: &Path) -> Result<()> {
    let svg_document = render_svg(inputs, options)?;

    match OutputFormat::from_path(out_path) {
        OutputFormat::Svg => {
            std::fs::write(out_path, svg_document)
                .with_context(|| format!("write SVG to {out_path:?}"))?;
        }
        OutputFormat::Pdf => {
            let pdf = svg_to_pdf(&svg_document)?;
            std::fs::write(out_path, pdf).with_context(|| format!("write PDF to {out_path:?}"))?;
        }
    }
    Ok(())
}

/// Renders directly from a work directory.
pub fn render_work_dir(
    work_dir: &Path,
    options: &RenderOptions,
    out_path: &Path,
) -> Result<RenderInputs> {
    let inputs = load_inputs(work_dir)?;
    render_to_file(&inputs, options, out_path)?;
    Ok(inputs)
}

/// Convenience for tests: renders to PDF bytes without touching the filesystem.
pub fn render_pdf_bytes(inputs: &RenderInputs, options: &RenderOptions) -> Result<Vec<u8>> {
    let svg_document = render_svg(inputs, options)?;
    svg_to_pdf(&svg_document)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::plot::io::{Cytoband, Region, SampleTrack, SignalBin, Transcript};

    fn sample_track(name: &str, peak: f64) -> SampleTrack {
        let bins = (0..10)
            .map(|index| SignalBin {
                chromosome: "chr1".to_string(),
                start: 1000 + index * 100,
                end: 1100 + index * 100,
                size: 100,
                // Ramp up to `peak` so bars have varying heights.
                max: peak * (index as f64 + 1.0) / 10.0,
            })
            .collect();
        SampleTrack {
            sample: name.to_string(),
            bins,
        }
    }

    fn inputs(track_count: usize) -> RenderInputs {
        RenderInputs {
            tracks: (0..track_count)
                .map(|index| sample_track(&format!("s{index}"), 10.0 * (index as f64 + 1.0)))
                .collect(),
            region: Region {
                chromosome: "chr1".to_string(),
                start: 1000,
                end: 2000,
                binsize: 100,
                loci: "chr1:1000-2000".to_string(),
            },
            transcripts: vec![Transcript {
                chromosome: "chr1".to_string(),
                strand: "+".to_string(),
                transcript: "T1".to_string(),
                gene: "G1".to_string(),
                start: 1100,
                end: 1800,
                exons: vec![(1100, 1200), (1700, 1800)],
            }],
            cytobands: vec![Cytoband {
                start: 0,
                end: 2000,
                stain: "gneg".to_string(),
                color: "#FFFFFF".to_string(),
            }],
        }
    }

    #[test]
    fn infers_format_from_extension() {
        assert_eq!(
            OutputFormat::from_path(Path::new("out.pdf")),
            OutputFormat::Pdf
        );
        assert_eq!(
            OutputFormat::from_path(Path::new("out.svg")),
            OutputFormat::Svg
        );
        assert_eq!(
            OutputFormat::from_path(Path::new("out.SVG")),
            OutputFormat::Svg
        );
        // Unknown or missing extensions default to PDF.
        assert_eq!(
            OutputFormat::from_path(Path::new("out")),
            OutputFormat::Pdf
        );
    }

    #[test]
    fn per_track_scaling_differs_when_not_grouped() {
        let inputs = inputs(2);
        let options = RenderOptions::default();
        let maxima = resolve_y_maxima(&inputs.tracks, &options);
        // Each track scales to its own peak.
        assert_eq!(maxima.len(), 2);
        assert!(
            (maxima[0] - 10.0).abs() < 0.01,
            "first track should peak at 10, got {}",
            maxima[0]
        );
        assert!(
            (maxima[1] - 20.0).abs() < 0.01,
            "second track should peak at 20, got {}",
            maxima[1]
        );
    }

    #[test]
    fn group_auto_scale_shares_one_maximum() {
        let inputs = inputs(2);
        let options = RenderOptions {
            group_auto_scale: true,
            ..Default::default()
        };
        let maxima = resolve_y_maxima(&inputs.tracks, &options);
        assert_eq!(maxima[0], maxima[1], "shared scale must match");
        assert!((maxima[0] - 20.0).abs() < 0.01, "shared max is the largest peak");
    }

    #[test]
    fn explicit_y_max_wins_and_cycles() {
        let inputs = inputs(3);
        let options = RenderOptions {
            y_max: Some(vec![99.0]),
            group_auto_scale: true,
            ..Default::default()
        };
        let maxima = resolve_y_maxima(&inputs.tracks, &options);
        assert_eq!(maxima, vec![99.0, 99.0, 99.0]);
    }

    #[test]
    fn silent_track_gets_a_usable_range() {
        // A track whose bins are all zero would map every bar to the baseline;
        // the renderer must fall back to a non-zero range instead of dividing by 0.
        let mut track = sample_track("silent", 0.0);
        for bin in &mut track.bins {
            bin.max = 0.0;
        }
        let options = RenderOptions::default();
        let maxima = resolve_y_maxima(&[track], &options);
        assert_eq!(maxima, vec![1.0]);
    }

    #[test]
    fn svg_contains_every_requested_panel() {
        let inputs = inputs(2);
        let svg_document = render_svg(&inputs, &RenderOptions::default()).expect("render svg");

        // Outer document + one nested panel per track in the layout.
        assert!(svg_document.starts_with("<svg"));
        assert!(svg_document.ends_with("</svg>"));
        // Two bigWigs + gene + scale + ideogram.
        assert!(
            svg_document.matches("<svg ").count() >= 6,
            "expected at least 6 nested svg elements, found {}",
            svg_document.matches("<svg ").count()
        );
        assert!(svg_document.contains("s0") && svg_document.contains("s1"));
        assert!(svg_document.contains("T1 [G1]"), "gene label missing");
        assert!(svg_document.contains("chr1:1000-2000"), "loci label missing");
    }

    #[test]
    fn renders_a_valid_pdf() {
        let inputs = inputs(2);
        let pdf = render_pdf_bytes(&inputs, &RenderOptions::default()).expect("render pdf");
        assert!(pdf.len() > 1000, "pdf suspiciously small: {} bytes", pdf.len());
        assert_eq!(&pdf[..5], b"%PDF-");
    }

    #[test]
    fn optional_panels_can_be_disabled() {
        let inputs = inputs(1);
        let options = RenderOptions {
            show_ideogram: false,
            draw_gene_track: false,
            ..Default::default()
        };
        let svg_document = render_svg(&inputs, &options).expect("render");
        assert!(!svg_document.contains("T1 [G1]"), "gene panel should be absent");
        // The scale panel is always drawn, regardless of the optional flags.
        assert!(svg_document.contains("chr1:1000-2000"));
    }

    #[test]
    fn rejects_canvas_with_no_room_for_data() {
        let inputs = inputs(1);
        let options = RenderOptions {
            width: 50.0,
            left_margin: 60.0,
            ..Default::default()
        };
        assert!(render_svg(&inputs, &options).is_err());
    }

    #[test]
    fn every_panel_title_uses_sans_serif_stack() {
        let inputs = inputs(2);
        let svg_document = render_svg(&inputs, &RenderOptions::default()).expect("render");
        // Text must be resolvable by the bundled fallback font.
        assert!(svg_document.contains("sans-serif"));
        assert!(!svg_document.contains("font-family=\"\""));
    }
}
