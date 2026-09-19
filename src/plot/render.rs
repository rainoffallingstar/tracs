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
    self, draw_gene_panel, draw_ideogram_panel, draw_overlay_panel, draw_peaks_panel,
    draw_scale_panel, draw_signal_panel, SvgWriter, XAxis, YAxis,
};

/// One `par(mar=)` unit in points.
///
/// R's default `cin` is `c(0.15, 0.2)` inches, so one margin line is 0.2in =
/// 14.4pt at the 72dpi the renderer works in.
const R_LINE_HEIGHT: f64 = 14.4;

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
    /// Explicit y minima, one per sample (or cycled).
    pub y_min: Option<Vec<f64>>,
    /// Auto-scale all bigWig tracks to a shared maximum.
    pub group_auto_scale: bool,
    /// Draw every bigWig in one panel as a line plot (`track_overlay`).
    pub track_overlay: bool,
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
            y_min: None,
            group_auto_scale: false,
            track_overlay: false,
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

/// Resolved y-axis range for one bigWig panel.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct YRange {
    pub min: f64,
    pub max: f64,
}

impl YRange {
    /// Guards against a degenerate range, e.g. a track whose bins are all zero.
    fn sanitised(self) -> Self {
        let max = if self.max.is_finite() && self.max > self.min {
            self.max
        } else {
            self.min + 1.0
        };
        let min = if self.min.is_finite() { self.min } else { 0.0 };
        YRange { min, max }
    }
}

/// Resolves the y range for each bigWig panel, reproducing `track_plot()`.
///
/// `track_plot()` uses `c(min(x$max), max(x$max))` per track rather than
/// `0 ..= max`, so a track that never reaches zero still fills its panel. The
/// precedence is: explicit `y_max`/`y_min` wins, then `groupAutoScale` shares one
/// range across tracks, otherwise each track gets its own. Values are rounded to
/// two decimals as R does.
fn resolve_y_ranges(tracks: &[io::SampleTrack], options: &RenderOptions) -> Vec<YRange> {
    let per_track: Vec<(f64, f64)> = tracks
        .iter()
        .map(|track| {
            let max = track.max_signal();
            let min = track.min_signal();
            (min, max)
        })
        .collect();

    // Explicit `y_max` / `y_min` override everything, cycling when shorter.
    let explicit_max = options.y_max.as_ref().filter(|values| !values.is_empty());
    let explicit_min = options
        .y_min
        .as_ref()
        .filter(|values| !values.is_empty());
    if explicit_max.is_some() || explicit_min.is_some() {
        return (0..tracks.len())
            .map(|index| {
                let (min, max) = per_track[index];
                YRange {
                    min: explicit_min
                        .map(|values| values[index % values.len()])
                        .unwrap_or(min),
                    max: explicit_max
                        .map(|values| values[index % values.len()])
                        .unwrap_or(max),
                }
                .sanitised()
            })
            .collect();
    }

    if options.group_auto_scale {
        let shared_min = per_track
            .iter()
            .map(|(min, _)| *min)
            .filter(|value| value.is_finite())
            .fold(f64::INFINITY, f64::min);
        let shared_max = per_track
            .iter()
            .map(|(_, max)| *max)
            .filter(|value| value.is_finite())
            .fold(f64::NEG_INFINITY, f64::max);
        let shared = YRange {
            min: if shared_min.is_finite() { shared_min } else { 0.0 },
            max: shared_max,
        }
        .sanitised();
        return vec![YRange {
            min: round_two(shared.min),
            max: round_two(shared.max),
        }; tracks.len()];
    }

    per_track
        .into_iter()
        .map(|(min, max)| {
            let range = YRange { min, max }.sanitised();
            YRange {
                min: round_two(range.min),
                max: round_two(range.max),
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
        // R sets `ntracks = 1` for the overlay so all samples share one panel.
        bigwig_count: if options.track_overlay {
            1
        } else {
            inputs.tracks.len()
        },
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
    let y_ranges = resolve_y_ranges(&inputs.tracks, options);

    // Overlay collapses every bigWig into a single panel, exactly as R does with
    // `ntracks = 1`. Its range spans the min/max across *all* samples.
    let overlay_range = if options.track_overlay {
        resolve_y_ranges(
            &inputs.tracks,
            &RenderOptions {
                group_auto_scale: true,
                ..options.clone()
            },
        )
        .first()
        .copied()
    } else {
        None
    };

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

        // `track_plot()` sets per-track `par(mar=)`; reproducing those keeps the
        // data area proportions close to R's rather than using one shared inset.
        let margins = panel.kind.margins_lines();
        let top_inset = margins[0] * R_LINE_HEIGHT;
        let bottom_inset = margins[2] * R_LINE_HEIGHT;
        let _ = margins[1];

        writer.panel(0.0, panel_top, options.width, panel_height, |draw| {
            // Track titles sit in the top inset; the data area uses the rest.
            let title_space = if options.track_names_to_left {
                top_inset
            } else {
                top_inset.max(options.panel_title_space)
            };
            let plot_top = title_space;
            let plot_bottom = (panel_height - bottom_inset).max(plot_top + 1.0);

            match panel.kind {
                TrackKind::BigWig => {
                    let index = panel.bigwig_index.unwrap_or(0);

                    // Overlay mode funnels every sample into this one panel,
                    // matching R's `ntracks = 1`.
                    if let Some(range) = overlay_range {
                        let names: Vec<String> = inputs
                            .tracks
                            .iter()
                            .enumerate()
                            .map(|(sample_index, track)| {
                                options
                                    .track_names
                                    .as_ref()
                                    .and_then(|names| names.get(sample_index))
                                    .cloned()
                                    .unwrap_or_else(|| track.sample.clone())
                            })
                            .collect();
                        draw_overlay_panel(
                            draw,
                            &inputs.tracks,
                            axis,
                            YAxis::new(plot_top, plot_bottom, range.min, range.max),
                            &colors,
                            options.show_axis,
                            &names,
                            options.font_size,
                        );
                        return;
                    }

                    let Some(track) = inputs.tracks.get(index) else {
                        return;
                    };
                    let range = y_ranges.get(index).copied().unwrap_or(YRange {
                        min: 0.0,
                        max: 1.0,
                    });
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
                        YAxis::new(plot_top, plot_bottom, range.min, range.max),
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
#[allow(dead_code)]
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
    fn per_track_ranges_are_independent_without_grouping() {
        let inputs = inputs(2);
        let options = RenderOptions::default();
        let ranges = resolve_y_ranges(&inputs.tracks, &options);
        assert_eq!(ranges.len(), 2);

        // Each track spans its own min..max, matching track_plot()'s
        // `c(min(x$max), max(x$max))` rather than a shared 0..max.
        assert!((ranges[0].max - 10.0).abs() < 0.01, "track 0 max: {}", ranges[0].max);
        assert!((ranges[1].max - 20.0).abs() < 0.01, "track 1 max: {}", ranges[1].max);
        // The sample ramp starts at 1/10 of the peak, so the minimum is non-zero.
        assert!((ranges[0].min - 1.0).abs() < 0.01, "track 0 min: {}", ranges[0].min);
        assert!((ranges[1].min - 2.0).abs() < 0.01, "track 1 min: {}", ranges[1].min);
    }

    #[test]
    fn group_auto_scale_shares_one_range() {
        let inputs = inputs(2);
        let options = RenderOptions {
            group_auto_scale: true,
            ..Default::default()
        };
        let ranges = resolve_y_ranges(&inputs.tracks, &options);
        assert_eq!(ranges[0], ranges[1], "shared scale must match exactly");
        // Spans the union of both tracks: min over both, max over both.
        assert!((ranges[0].min - 1.0).abs() < 0.01, "shared min: {}", ranges[0].min);
        assert!((ranges[0].max - 20.0).abs() < 0.01, "shared max: {}", ranges[0].max);
    }

    #[test]
    fn explicit_y_max_and_min_override_and_cycle() {
        let inputs = inputs(3);
        let options = RenderOptions {
            y_max: Some(vec![99.0]),
            y_min: Some(vec![-5.0]),
            group_auto_scale: true,
            ..Default::default()
        };
        let ranges = resolve_y_ranges(&inputs.tracks, &options);
        assert_eq!(ranges.len(), 3);
        for range in &ranges {
            assert_eq!(range.max, 99.0);
            assert_eq!(range.min, -5.0);
        }
    }

    #[test]
    fn explicit_y_max_alone_keeps_the_track_minimum() {
        let inputs = inputs(1);
        let options = RenderOptions {
            y_max: Some(vec![50.0]),
            ..Default::default()
        };
        let ranges = resolve_y_ranges(&inputs.tracks, &options);
        assert_eq!(ranges[0].max, 50.0);
        // y_min was not given, so the track's own minimum is kept.
        assert!((ranges[0].min - 1.0).abs() < 0.01, "min: {}", ranges[0].min);
    }

    #[test]
    fn silent_track_gets_a_usable_range() {
        // A track whose bins are all zero would collapse the axis; the renderer
        // must fall back to a non-zero span instead of dividing by zero.
        let mut track = sample_track("silent", 0.0);
        for bin in &mut track.bins {
            bin.max = 0.0;
        }
        let options = RenderOptions::default();
        let ranges = resolve_y_ranges(&[track], &options);
        assert!(
            ranges[0].max > ranges[0].min,
            "degenerate range: {:?}",
            ranges[0]
        );
    }

    #[test]
    fn overlay_collapses_bigwig_panels_and_spans_all_samples() {
        let inputs = inputs(3);
        let options = RenderOptions {
            track_overlay: true,
            show_ideogram: false,
            draw_gene_track: false,
            ..Default::default()
        };

        // R sets `ntracks = 1`, so the layout holds exactly one bigWig panel.
        let layout = make_layout(&LayoutRequest {
            bigwig_height: options.bigwig_height,
            peaks_height: options.peaks_height,
            gene_height: options.gene_height,
            scale_height: options.scale_height,
            chromhmm_height: options.chromhmm_height,
            cytoband_height: options.cytoband_height,
            bigwig_count: 1,
            has_peaks: false,
            has_chromhmm: false,
            has_gene: false,
            has_cytoband: false,
            layout_ord: Vec::new(),
        });
        let bigwig_panels = layout
            .panels
            .iter()
            .filter(|panel| panel.kind == TrackKind::BigWig)
            .count();
        assert_eq!(bigwig_panels, 1, "overlay must use a single bigWig panel");

        let svg_document = render_svg(&inputs, &options).expect("render overlay");
        // Every sample is drawn as a line, and all three names appear in the key.
        let paths = svg_document.matches("<path").count();
        assert_eq!(paths, 3, "expected one polyline per sample, got {paths}");
        for name in ["s0", "s1", "s2"] {
            assert!(svg_document.contains(name), "overlay legend missing {name}");
        }
    }

    #[test]
    fn overlay_is_off_by_default() {
        let inputs = inputs(2);
        let options = RenderOptions::default();
        assert!(!options.track_overlay);
        let svg_document = render_svg(&inputs, &options).expect("render");
        // Bars use <rect>, so a default render must contain no sample polylines.
        assert!(
            !svg_document.contains("<path d=\"M "),
            "bar mode should not emit overlay polylines"
        );
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
