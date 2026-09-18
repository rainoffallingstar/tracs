//! SVG track panels.
//!
//! Each track is rendered as a standalone SVG fragment with its own coordinate
//! system, so every panel keeps an independent y-scale (the core requirement of
//! `track_plot()`). A small SVG writer is used rather than a plotting crate
//! because the layout needs direct control over panel geometry: the charting
//! libraries evaluated for this project either expose no panel rectangle, or
//! corrupt margins/ids when many panels are composed into one document.
//!
//! Coordinates are emitted in SVG user units (1 unit == 1 pt at 72 dpi), so the
//! final page size is simply the canvas size converted to points.

use std::fmt::Write as _;

use crate::plot::io::{Cytoband, SampleTrack, Transcript};
use crate::plot::pretty;

/// Colour used for gene model exons and introns (`track_plot()`'s `exon_col`).
pub const GENE_COLOR: &str = "#192a56";
/// Colour for the ideogram's highlighted region.
pub const IDEOGRAM_HIGHLIGHT: &str = "#d35400";
/// Border colour used by peaks and cytoband panels.
pub const BORDER_COLOR: &str = "#34495e";
/// Default track colour (`track_plot()`'s documented `col` default).
pub const DEFAULT_TRACK_COLOR: &str = "#2f3640";

/// Escape the five characters that matter inside XML text/attribute values.
fn escape(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    for ch in text.chars() {
        match ch {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&apos;"),
            _ => out.push(ch),
        }
    }
    out
}

/// Incrementally builds an SVG document.
///
/// Panels are written as nested `<svg>` elements positioned with `x`/`y`, which
/// gives each one an independent coordinate system without touching global
/// state.
pub struct SvgWriter {
    body: String,
    width: f64,
    height: f64,
}

impl SvgWriter {
    pub fn new(width: f64, height: f64) -> Self {
        Self {
            body: String::new(),
            width,
            height,
        }
    }

    /// Appends a nested panel occupying `(x, y, width, height)`.
    ///
    /// The closure receives a writer whose coordinates start at the panel's
    /// own origin.
    pub fn panel<F>(&mut self, x: f64, y: f64, width: f64, height: f64, draw: F)
    where
        F: FnOnce(&mut PanelWriter),
    {
        let mut panel = PanelWriter {
            body: String::new(),
            width,
            height,
        };
        draw(&mut panel);
        let _ = write!(
            self.body,
            "<svg x=\"{x:.3}\" y=\"{y:.3}\" width=\"{width:.3}\" height=\"{height:.3}\" \
             viewBox=\"0 0 {width:.3} {height:.3}\" overflow=\"visible\">{}</svg>",
            panel.body
        );
    }

    /// Renders the document, ensuring the panel canvas is at least as large as
    /// the content so labels are not clipped away.
    pub fn finish(self) -> String {
        format!(
            "<svg width=\"{w:.3}\" height=\"{h:.3}\" viewBox=\"0 0 {w:.3} {h:.3}\" \
             xmlns=\"http://www.w3.org/2000/svg\">{}</svg>",
            self.body,
            w = self.width,
            h = self.height
        )
    }
}

/// Drawing surface for a single panel.
///
/// `plot_left`/`plot_right` delineate the data area; labels are placed in the
/// margins. Panels are emitted with `<rect>` for bars and `<line>`/`<path>` for
/// rules and introns, which is what makes the output independent of any
/// plotting engine's internal state.
pub struct PanelWriter {
    body: String,
    width: f64,
    height: f64,
}

impl PanelWriter {
    /// Total panel width.
    pub fn width(&self) -> f64 {
        self.width
    }

    /// Total panel height.
    pub fn height(&self) -> f64 {
        self.height
    }

    /// Filled rectangle.
    pub fn rect(&mut self, x: f64, y: f64, width: f64, height: f64, fill: &str) {
        let _ = write!(
            self.body,
            "<rect x=\"{x:.3}\" y=\"{y:.3}\" width=\"{width:.3}\" height=\"{height:.3}\" fill=\"{}\"/>",
            escape(fill)
        );
    }

    /// Rectangle with an outline.
    pub fn rect_stroked(
        &mut self,
        x: f64,
        y: f64,
        width: f64,
        height: f64,
        fill: &str,
        stroke: &str,
        stroke_width: f64,
    ) {
        let _ = write!(
            self.body,
            "<rect x=\"{x:.3}\" y=\"{y:.3}\" width=\"{width:.3}\" height=\"{height:.3}\" \
             fill=\"{}\" stroke=\"{}\" stroke-width=\"{stroke_width:.3}\"/>",
            escape(fill),
            escape(stroke)
        );
    }

    /// Straight line.
    pub fn line(&mut self, x1: f64, y1: f64, x2: f64, y2: f64, stroke: &str, stroke_width: f64) {
        let _ = write!(
            self.body,
            "<line x1=\"{x1:.3}\" y1=\"{y1:.3}\" x2=\"{x2:.3}\" y2=\"{y2:.3}\" \
             stroke=\"{}\" stroke-width=\"{stroke_width:.3}\"/>",
            escape(stroke)
        );
    }

    /// Dashed horizontal line, used for the scale bar.
    pub fn dashed_line(&mut self, x1: f64, y: f64, x2: f64, stroke: &str) {
        let _ = write!(
            self.body,
            "<line x1=\"{x1:.3}\" y1=\"{y:.3}\" x2=\"{x2:.3}\" y2=\"{y:.3}\" \
             stroke=\"{}\" stroke-width=\"1\" stroke-dasharray=\"4 3\"/>",
            escape(stroke)
        );
    }

    /// Text label.
    ///
    /// `anchor` is `"start"`, `"middle"` or `"end"`. The font stack ends in
    /// `sans-serif` so the bundled fallback in [`crate::plot::fonts`] applies when
    /// the host has no fonts.
    pub fn text(
        &mut self,
        x: f64,
        y: f64,
        label: &str,
        font_size: f64,
        anchor: &str,
        fill: &str,
    ) {
        let _ = write!(
            self.body,
            "<text x=\"{x:.3}\" y=\"{y:.3}\" font-size=\"{font_size:.1}\" \
             font-family=\"Inter, Helvetica, Arial, sans-serif\" \
             fill=\"{}\" text-anchor=\"{}\">{}</text>",
            escape(fill),
            escape(anchor),
            escape(label)
        );
    }

    /// Dashed vertical line, used for axis ticks.
    pub fn tick(&mut self, x: f64, y1: f64, y2: f64, stroke: &str) {
        let _ = write!(
            self.body,
            "<line x1=\"{x:.3}\" y1=\"{y1:.3}\" x2=\"{x:.3}\" y2=\"{y2:.3}\" \
             stroke=\"{}\" stroke-width=\"1\"/>",
            escape(stroke)
        );
    }
}

/// Horizontal placement inside a panel: where the data area starts and ends.
#[derive(Clone, Copy, Debug)]
pub struct XAxis {
    /// Left edge of the data area, in panel units.
    pub plot_left: f64,
    /// Right edge of the data area.
    pub plot_right: f64,
    /// Data range mapped onto the area.
    pub data_start: f64,
    pub data_end: f64,
}

impl XAxis {
    /// Maps a genomic coordinate to a panel x position.
    pub fn map(&self, value: f64) -> f64 {
        let span = self.data_end - self.data_start;
        if span <= 0.0 {
            return self.plot_left;
        }
        let fraction = (value - self.data_start) / span;
        self.plot_left + fraction * (self.plot_right - self.plot_left)
    }

    /// Width of the data area.
    pub fn plot_width(&self) -> f64 {
        (self.plot_right - self.plot_left).max(0.0)
    }
}

/// Vertical placement for a track whose y axis runs `0 ..= y_max`.
#[derive(Clone, Copy, Debug)]
pub struct YAxis {
    /// Top of the data area (y = `y_max`).
    pub plot_top: f64,
    /// Baseline (y = 0).
    pub plot_bottom: f64,
    pub y_max: f64,
}

impl YAxis {
    /// Maps a signal value to a panel y position.
    pub fn map(&self, value: f64) -> f64 {
        if self.y_max <= 0.0 {
            return self.plot_bottom;
        }
        let fraction = (value / self.y_max).clamp(0.0, 1.0);
        self.plot_bottom - fraction * (self.plot_bottom - self.plot_top)
    }

    /// Height of the data area.
    pub fn plot_height(&self) -> f64 {
        (self.plot_bottom - self.plot_top).max(0.0)
    }
}

/// Formats a genomic coordinate the way `track_plot()` does: megabases and
/// hundred-kilobases are abbreviated, everything else prints as-is.
///
/// The branch order matters (`> 1e6` before `> 1e5`), and the K divisor is
/// `1e5`, not `1e3` -- both mirror `track_plot()` line ~894.
pub fn format_coordinate(value: f64) -> String {
    if value > 1e6 {
        format!("{}M", trim_number(value / 1e6))
    } else if value > 100_000.0 {
        format!("{}K", trim_number(value / 1e5))
    } else {
        trim_number(value)
    }
}

/// Formats a number without a trailing `.0`, matching R's default printing.
fn trim_number(value: f64) -> String {
    if value.fract() == 0.0 {
        format!("{}", value as i64)
    } else {
        let text = format!("{value}");
        text.trim_end_matches('0').trim_end_matches('.').to_string()
    }
}

/// Formats a y-axis tick value, matching `track_plot()`'s rounding.
pub fn format_tick(value: f64) -> String {
    trim_number(value)
}

/// Draws a bigWig signal panel: one filled bar per bin, scaled to `y_max`.
///
/// When `show_axis` is false the y range is annotated as `[min-max]` in the top
/// left corner, which is what `track_plot()` does instead of axis ticks.
pub fn draw_signal_panel(
    panel: &mut PanelWriter,
    track: &SampleTrack,
    x_axis: XAxis,
    y_axis: YAxis,
    color: &str,
    show_axis: bool,
    track_name: &str,
    track_name_left: bool,
    font_size: f64,
) {
    let baseline = y_axis.plot_bottom;
    for bin in &track.bins {
        if !bin.max.is_finite() || bin.max <= 0.0 {
            continue;
        }
        let left = x_axis.map(bin.start as f64);
        let right = x_axis.map(bin.end as f64);
        let width = (right - left).max(0.5);
        let top = y_axis.map(bin.max);
        let height = (baseline - top).max(0.0);
        if height > 0.0 {
            panel.rect(left, top, width, height, color);
        }
    }

    // Baseline.
    panel.line(
        x_axis.plot_left,
        baseline,
        x_axis.plot_right,
        baseline,
        "black",
        1.0,
    );

    if show_axis {
        let ticks = pretty::pretty_range(0.0, y_axis.y_max, 5);
        for tick in &ticks.values {
            if *tick < 0.0 || *tick > y_axis.y_max {
                continue;
            }
            let y = y_axis.map(*tick);
            panel.line(x_axis.plot_left - 4.0, y, x_axis.plot_left, y, "black", 1.0);
            panel.text(
                x_axis.plot_left - 6.0,
                y + 3.0,
                &format_tick(*tick),
                font_size,
                "end",
                "black",
            );
        }
    } else {
        // R prints "[0-60]" style range annotations when the axis is hidden.
        let label = format!("[0-{}]", format_tick(y_axis.y_max));
        panel.text(
            x_axis.plot_left,
            y_axis.plot_top + font_size,
            &label,
            font_size,
            "start",
            "black",
        );
    }

    draw_track_name(
        panel,
        x_axis,
        y_axis,
        track_name,
        track_name_left,
        font_size,
    );
}

/// Places the track label either inside the panel (left) or as a centred title.
fn draw_track_name(
    panel: &mut PanelWriter,
    x_axis: XAxis,
    y_axis: YAxis,
    name: &str,
    to_left: bool,
    font_size: f64,
) {
    if name.is_empty() {
        return;
    }
    if to_left {
        // R anchors these at the left edge of the data area, right-aligned.
        panel.text(
            x_axis.plot_left - 6.0,
            y_axis.plot_top + font_size,
            name,
            font_size,
            "end",
            "black",
        );
    } else {
        panel.text(
            (x_axis.plot_left + x_axis.plot_right) / 2.0,
            y_axis.plot_top - 2.0,
            name,
            font_size,
            "middle",
            "black",
        );
    }
}

/// Draws the gene model panel: intron line, exon boxes and strand arrows.
pub fn draw_gene_panel(
    panel: &mut PanelWriter,
    transcripts: &[Transcript],
    x_axis: XAxis,
    font_size: f64,
) {
    let count = transcripts.len();
    if count == 0 {
        return;
    }

    // One row per transcript, spread evenly over the panel height.
    let row_height = panel.height() / count as f64;
    for (index, transcript) in transcripts.iter().enumerate() {
        let centre = row_height * (index as f64 + 0.5);

        // Intron line spanning the transcript.
        panel.line(
            x_axis.map(transcript.start as f64),
            centre,
            x_axis.map(transcript.end as f64),
            centre,
            GENE_COLOR,
            1.0,
        );

        // Exon boxes, drawn taller than the intron line.
        let exon_half_height = (row_height * 0.25).min(6.0);
        for (exon_start, exon_end) in &transcript.exons {
            let left = x_axis.map(*exon_start as f64);
            let right = x_axis.map(*exon_end as f64);
            panel.rect(
                left,
                centre - exon_half_height,
                (right - left).max(0.5),
                exon_half_height * 2.0,
                GENE_COLOR,
            );
        }

        // Strand arrows at R's pretty() tick positions along the transcript.
        if transcript.end > transcript.start {
            let arrows = pretty::pretty_range(
                transcript.start as f64,
                transcript.end as f64,
                // Fewer arrows than the default 5 intervals, since each arrow is
                // a glyph rather than a tick.
                3,
            );
            let glyph = if transcript.strand == "+" { ">" } else { "<" };
            for position in &arrows.values {
                let clamped = position.clamp(transcript.start as f64, transcript.end as f64);
                panel.text(
                    x_axis.map(clamped),
                    centre + font_size * 0.35,
                    glyph,
                    font_size,
                    "middle",
                    GENE_COLOR,
                );
            }
        }

        // R labels genes to the right of the model; transcripts also show their id.
        let label = if transcript.gene.is_empty() {
            transcript.transcript.clone()
        } else {
            format!("{} [{}]", transcript.transcript, transcript.gene)
        };
        panel.text(
            x_axis.plot_left,
            centre - font_size * 0.5,
            &label,
            font_size,
            "start",
            "black",
        );
    }
}

/// Draws the coordinate scale panel: a dashed ruler with tick labels.
pub fn draw_scale_panel(
    panel: &mut PanelWriter,
    x_axis: XAxis,
    chromosome: &str,
    start: u64,
    end: u64,
    font_size: f64,
) {
    let ruler_y = panel.height() * 0.5;
    panel.dashed_line(x_axis.plot_left, ruler_y, x_axis.plot_right, "black");

    let ticks = pretty::pretty_range(start as f64, end as f64, 5);
    for position in &ticks.values {
        let x = x_axis.map(*position);
        // R draws the tick as a short vertical mark above the ruler.
        panel.tick(x, ruler_y - 4.0, ruler_y, "black");
        panel.text(
            x,
            ruler_y - 6.0,
            &format_coordinate(*position),
            font_size,
            "middle",
            "black",
        );
    }

    // R annotates the full region at the right edge.
    panel.text(
        x_axis.plot_right,
        ruler_y - 6.0 - font_size - 2.0,
        &format!("{chromosome}:{start}-{end}"),
        font_size,
        "end",
        "black",
    );
}

/// Draws the ideogram panel: cytobands for the whole chromosome with the
/// plotted region highlighted.
pub fn draw_ideogram_panel(
    panel: &mut PanelWriter,
    bands: &[Cytoband],
    chromosome: &str,
    region_start: u64,
    region_end: u64,
    font_size: f64,
) {
    if bands.is_empty() {
        return;
    }
    // The ideogram spans the entire chromosome, not just the plotted region, so
    // it uses its own axis covering `0 ..= chromosome end`.
    let chromosome_end = bands.iter().map(|band| band.end).max().unwrap_or(0);
    if chromosome_end == 0 {
        return;
    }
    let ideogram_axis = XAxis {
        plot_left: 0.0,
        plot_right: panel.width(),
        data_start: 0.0,
        data_end: chromosome_end as f64,
    };

    let band_top = panel.height() * 0.25;
    let band_height = panel.height() * 0.5;
    for band in bands {
        let left = ideogram_axis.map(band.start as f64);
        let right = ideogram_axis.map(band.end as f64);
        panel.rect_stroked(
            left,
            band_top,
            (right - left).max(0.5),
            band_height,
            &band.color,
            BORDER_COLOR,
            0.5,
        );
    }

    // Highlighted region, drawn taller so it stands out over the bands.
    let highlight_left = ideogram_axis.map(region_start as f64);
    let highlight_right = ideogram_axis.map(region_end as f64);
    panel.rect(
        highlight_left,
        band_top - 2.0,
        (highlight_right - highlight_left).max(1.0),
        band_height + 4.0,
        IDEOGRAM_HIGHLIGHT,
    );

    panel.text(
        0.0,
        band_top + band_height / 2.0,
        chromosome,
        font_size,
        "end",
        "black",
    );
}

/// Draws the peaks panel: one row per region set, with overlapping intervals
/// filled.
pub fn draw_peaks_panel(
    panel: &mut PanelWriter,
    region_sets: &[(String, Vec<(u64, u64)>)],
    x_axis: XAxis,
    font_size: f64,
) {
    let count = region_sets.len();
    if count == 0 {
        return;
    }
    let row_height = panel.height() / count as f64;

    for (index, (name, regions)) in region_sets.iter().enumerate() {
        let centre = row_height * (index as f64 + 0.5);
        // Light band showing the full plotted span for this row.
        panel.rect(
            x_axis.plot_left,
            centre - 0.5,
            x_axis.plot_width(),
            1.0,
            "gray90",
        );
        for (start, end) in regions {
            let left = x_axis.map(*start as f64);
            let right = x_axis.map(*end as f64);
            let height = (row_height * 0.8).max(2.0);
            panel.rect(
                left,
                centre - height / 2.0,
                (right - left).max(0.5),
                height,
                BORDER_COLOR,
            );
        }
        panel.text(
            x_axis.plot_left - 6.0,
            centre + font_size * 0.35,
            name,
            font_size,
            "end",
            "black",
        );
    }
}

/// Assigns each track a colour, cycling the palette when there are more tracks
/// than colours.
///
/// Mirrors `track_plot()`'s `col = rep(col, length(summary_list))`.
pub fn resolve_colors(requested: &[String], count: usize) -> Vec<String> {
    if count == 0 {
        return Vec::new();
    }
    if requested.is_empty() {
        return vec![DEFAULT_TRACK_COLOR.to_string(); count];
    }
    (0..count)
        .map(|index| requested[index % requested.len()].clone())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn escapes_xml_metacharacters() {
        assert_eq!(escape("a&b"), "a&amp;b");
        assert_eq!(escape("<tag>"), "&lt;tag&gt;");
        assert_eq!(escape("\"q\""), "&quot;q&quot;");
    }

    #[test]
    fn x_axis_maps_and_clamps() {
        let axis = XAxis {
            plot_left: 10.0,
            plot_right: 110.0,
            data_start: 100.0,
            data_end: 200.0,
        };
        assert_eq!(axis.map(100.0), 10.0);
        assert_eq!(axis.map(200.0), 110.0);
        assert_eq!(axis.map(150.0), 60.0);
        assert_eq!(axis.plot_width(), 100.0);
    }

    #[test]
    fn y_axis_maps_from_baseline_upward() {
        let axis = YAxis {
            plot_top: 0.0,
            plot_bottom: 100.0,
            y_max: 50.0,
        };
        assert_eq!(axis.map(0.0), 100.0);
        assert_eq!(axis.map(50.0), 0.0);
        assert_eq!(axis.map(25.0), 50.0);
        // Values beyond the axis are clamped, never drawn outside the panel.
        assert_eq!(axis.map(200.0), 0.0);
        assert_eq!(axis.plot_height(), 100.0);
    }

    #[test]
    fn coordinate_formatting_matches_trackplot_rules() {
        // > 1e6 => M (using 1e6 divisor)
        assert_eq!(format_coordinate(158_156_686.0), "158.156686M");
        // > 1e5 => K, and the divisor is 1e5, not 1e3
        assert_eq!(format_coordinate(158_000.0), "1.58K");
        // Plain below the thresholds
        assert_eq!(format_coordinate(12_345.0), "12345");
        assert_eq!(format_coordinate(0.0), "0");
    }

    #[test]
    fn resolves_colors_by_cycling() {
        let colors = resolve_colors(&["#aaa".to_string(), "#bbb".to_string()], 5);
        assert_eq!(colors, vec!["#aaa", "#bbb", "#aaa", "#bbb", "#aaa"]);

        let defaults = resolve_colors(&[], 2);
        assert_eq!(defaults, vec![DEFAULT_TRACK_COLOR, DEFAULT_TRACK_COLOR]);

        assert!(resolve_colors(&[], 0).is_empty());
    }

    #[test]
    fn writer_nests_panels_with_own_coordinates() {
        let mut writer = SvgWriter::new(200.0, 100.0);
        writer.panel(0.0, 0.0, 200.0, 50.0, |panel| {
            panel.rect(0.0, 0.0, 10.0, 10.0, "red");
        });
        writer.panel(0.0, 50.0, 200.0, 50.0, |panel| {
            panel.text(5.0, 5.0, "label", 10.0, "start", "black");
        });
        let svg = writer.finish();

        // Two nested panels, each with its own viewBox.
        assert_eq!(svg.matches("<svg").count(), 3, "outer + two panels");
        assert!(svg.contains("viewBox=\"0 0 200.000 50.000\""));
        assert!(svg.contains("y=\"50.000\""));
        // The label must stay inside its own panel, not the outer document.
        assert!(svg.contains(">label</text>"));
    }

    #[test]
    fn writer_emits_sans_serif_fallback_in_font_stack() {
        let mut writer = SvgWriter::new(100.0, 50.0);
        writer.panel(0.0, 0.0, 100.0, 50.0, |panel| {
            panel.text(1.0, 1.0, "x", 10.0, "start", "black");
        });
        let svg = writer.finish();
        assert!(
            svg.contains("sans-serif"),
            "font stack must end in sans-serif so the bundled fallback applies"
        );
    }
}
