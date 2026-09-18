//! Panel layout, ported from `trackplot.R`'s `.make_layout()`.
//!
//! `track_plot()` stacks one panel per track and lets the caller reorder them
//! with `layout_ord`. Two details matter for matching R's output:
//!
//! 1. `s` (scale) is injected into the requested order via `c(lord, "s")` before
//!    `intersect(lord, names(...))`, so it comes **first** when `layout_ord` is
//!    empty and **last** when the caller named any keys. Any keys the caller
//!    omitted are then appended in canonical order.
//! 2. Panels are *drawn* in a fixed order (`p, b, h, g, s, c`) while *heights*
//!    follow the resolved order. R reconciles the two by computing a draw index
//!    per row, which is what this module reproduces.
//!
//! Both behaviours are pinned by `tests/plot_layout_oracle.rs` against output
//! captured from R 4.6.0.

/// The six panel kinds `track_plot()` knows about.
///
/// The single-letter keys match R's (`layout_ord` is user-facing and documented
/// as `c("p", "b", "h", "g", "c")`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum TrackKind {
    /// Peaks / regions track (`p`).
    Peaks,
    /// bigWig signal tracks (`b`); one row per sample.
    BigWig,
    /// chromHMM track (`h`).
    ChromHmm,
    /// Gene model track (`g`).
    Gene,
    /// Coordinate scale bar (`s`); always present.
    Scale,
    /// Cytoband / ideogram track (`c`).
    Cytoband,
}

impl TrackKind {
    /// The single-letter key used by `layout_ord`.
    pub fn key(self) -> char {
        match self {
            TrackKind::Peaks => 'p',
            TrackKind::BigWig => 'b',
            TrackKind::ChromHmm => 'h',
            TrackKind::Gene => 'g',
            TrackKind::Scale => 's',
            TrackKind::Cytoband => 'c',
        }
    }

    /// Parses a `layout_ord` letter. Unknown letters are ignored, as in R.
    pub fn from_key(key: char) -> Option<Self> {
        match key {
            'p' => Some(TrackKind::Peaks),
            'b' => Some(TrackKind::BigWig),
            'h' => Some(TrackKind::ChromHmm),
            'g' => Some(TrackKind::Gene),
            's' => Some(TrackKind::Scale),
            'c' => Some(TrackKind::Cytoband),
            _ => None,
        }
    }

    /// Draw order used by `track_plot()`: `p, b, h, g, s, c` maps to 1..=6.
    fn draw_order(self) -> u8 {
        match self {
            TrackKind::Peaks => 1,
            TrackKind::BigWig => 2,
            TrackKind::ChromHmm => 3,
            TrackKind::Gene => 4,
            TrackKind::Scale => 5,
            TrackKind::Cytoband => 6,
        }
    }
}

/// Which panels and heights to include, mirroring `.make_layout()`'s arguments.
#[derive(Clone, Debug)]
pub struct LayoutRequest {
    /// One height per bigWig sample (`bw_track_height`).
    pub bigwig_height: f64,
    /// `peaks_track_height`.
    pub peaks_height: f64,
    /// `gene_track_height`.
    pub gene_height: f64,
    /// `scale_track_height`.
    pub scale_height: f64,
    /// `chromHMM_track_height`.
    pub chromhmm_height: f64,
    /// `cytoband_track_height`.
    pub cytoband_height: f64,
    /// Number of bigWig panels (one per sample).
    pub bigwig_count: usize,
    /// Whether a peaks track is being drawn.
    pub has_peaks: bool,
    /// Whether a chromHMM track is being drawn.
    pub has_chromhmm: bool,
    /// Whether a gene track is being drawn (`draw_gene_track`).
    pub has_gene: bool,
    /// Whether the ideogram is drawn (`show_ideogram`).
    pub has_cytoband: bool,
    /// User order from `layout_ord`.
    pub layout_ord: Vec<TrackKind>,
}

/// One resolved panel: what to draw, how tall, and in which draw slot.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Panel {
    pub kind: TrackKind,
    /// Index into the bigWig sample list; only meaningful for `TrackKind::BigWig`.
    pub bigwig_index: Option<usize>,
    /// Relative height, straight from the `*_track_height` arguments.
    pub height: f64,
    /// 1-based slot in the draw sequence (R's `ord_req4`).
    pub draw_index: usize,
}

/// A fully resolved layout: panels in draw order plus the height total.
#[derive(Clone, Debug, PartialEq)]
pub struct Layout {
    /// Panels in the order they should be emitted top to bottom.
    pub panels: Vec<Panel>,
    /// Sum of all panel heights, used to scale to a pixel canvas.
    pub total_height: f64,
}

impl Layout {
    /// Panel heights as a fraction of the total, for pixel allocation.
    pub fn fractional_heights(&self) -> Vec<f64> {
        if self.total_height <= 0.0 {
            return vec![0.0; self.panels.len()];
        }
        self.panels
            .iter()
            .map(|panel| panel.height / self.total_height)
            .collect()
    }
}

/// Resolves which panels to draw, their heights, and their draw order.
///
/// This is a port of `.make_layout()`. The draw-order reconciliation is the
/// subtle part: R computes, for each user-ordered row, the position it occupies
/// once all rows are sorted by the fixed draw order, then sorts the result back
/// into user order. Net effect: `draw_index` values are a permutation of
/// `1..=n`, but heights stay keyed to the user's sequence.
pub fn make_layout(request: &LayoutRequest) -> Layout {
    // R builds `lo_h_ord` over the keys that are actually kept; `b` expands to
    // one entry per bigWig. The `s` key is added to any user order.
    let mut ordered: Vec<TrackKind> = Vec::new();
    for kind in request.layout_ord.iter().copied() {
        if is_kept(kind, request) && !ordered.contains(&kind) {
            ordered.push(kind);
        }
    }
    // `lord = c(lord, "s")` then any missing kept keys are appended in the
    // canonical order of `lo_h_ord`.
    if is_kept(TrackKind::Scale, request) && !ordered.contains(&TrackKind::Scale) {
        ordered.push(TrackKind::Scale);
    }
    for kind in all_kinds() {
        if is_kept(kind, request) && !ordered.contains(&kind) {
            ordered.push(kind);
        }
    }

    // Expand each key into its rows, in the user order, recording heights.
    let mut rows: Vec<(TrackKind, Option<usize>, f64)> = Vec::new();
    for kind in ordered {
        let height = height_for(kind, request);
        match kind {
            TrackKind::BigWig => {
                for index in 0..request.bigwig_count {
                    rows.push((kind, Some(index), height));
                }
            }
            _ => rows.push((kind, None, height)),
        }
    }

    // R sorts the rows by the fixed draw order and records each row's position;
    // the original (user) sequence is then the order panels appear top to bottom.
    let mut by_draw: Vec<usize> = (0..rows.len()).collect();
    by_draw.sort_by_key(|&index| rows[index].0.draw_order());
    let mut draw_index_of_row = vec![0usize; rows.len()];
    for (slot, &row_index) in by_draw.iter().enumerate() {
        draw_index_of_row[row_index] = slot + 1;
    }

    let panels: Vec<Panel> = rows
        .into_iter()
        .enumerate()
        .map(|(index, (kind, bigwig_index, height))| Panel {
            kind,
            bigwig_index,
            height,
            draw_index: draw_index_of_row[index],
        })
        .collect();

    let total_height = panels.iter().map(|panel| panel.height).sum();
    Layout {
        panels,
        total_height,
    }
}

/// The canonical key order R uses when appending missing panels.
fn all_kinds() -> [TrackKind; 6] {
    [
        TrackKind::Peaks,
        TrackKind::BigWig,
        TrackKind::ChromHmm,
        TrackKind::Gene,
        TrackKind::Scale,
        TrackKind::Cytoband,
    ]
}

/// Mirrors `.make_layout()`'s `keep` vector: which keys are present at all.
fn is_kept(kind: TrackKind, request: &LayoutRequest) -> bool {
    match kind {
        // `"b"` is always kept, and peaks/chromHMM are data-driven.
        TrackKind::BigWig => true,
        TrackKind::Peaks => request.has_peaks,
        TrackKind::ChromHmm => request.has_chromhmm,
        TrackKind::Gene => request.has_gene,
        TrackKind::Scale => true,
        TrackKind::Cytoband => request.has_cytoband,
    }
}

fn height_for(kind: TrackKind, request: &LayoutRequest) -> f64 {
    match kind {
        TrackKind::Peaks => request.peaks_height,
        TrackKind::BigWig => request.bigwig_height,
        TrackKind::ChromHmm => request.chromhmm_height,
        TrackKind::Gene => request.gene_height,
        TrackKind::Scale => request.scale_height,
        TrackKind::Cytoband => request.cytoband_height,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn base_request() -> LayoutRequest {
        LayoutRequest {
            bigwig_height: 3.0,
            peaks_height: 2.0,
            gene_height: 2.0,
            scale_height: 2.0,
            chromhmm_height: 1.0,
            cytoband_height: 2.0,
            bigwig_count: 2,
            has_peaks: false,
            has_chromhmm: false,
            has_gene: true,
            has_cytoband: false,
            layout_ord: Vec::new(),
        }
    }

    #[test]
    fn default_order_matches_r_when_layout_ord_is_empty() {
        // With `lord = NULL`, R does `c(lord, "s")` first, so `s` lands at the
        // front and the remaining keys follow. Verified against the oracle row
        // with an empty layout_ord: "s:2,b:3,b:3,g:2".
        let layout = make_layout(&base_request());
        let observed: Vec<String> = layout
            .panels
            .iter()
            .map(|p| format!("{}:{}", p.kind.key(), p.height as i64))
            .collect();
        assert_eq!(observed.join(","), "s:2,b:3,b:3,g:2");

        assert_eq!(layout.panels[0].kind, TrackKind::Scale);
        assert_eq!(layout.panels[1].bigwig_index, Some(0));
        assert_eq!(layout.panels[2].bigwig_index, Some(1));
        assert_eq!(layout.panels[3].kind, TrackKind::Gene);
    }

    #[test]
    fn explicit_layout_ord_matches_r() {
        // Oracle row: n=2 lord=gb => "g:2,b:3,b:3,b:3,s:2"
        let mut request = base_request();
        request.bigwig_count = 3;
        request.layout_ord = vec![TrackKind::Gene, TrackKind::BigWig];
        let layout = make_layout(&request);
        let observed: Vec<String> = layout
            .panels
            .iter()
            .map(|p| format!("{}:{}", p.kind.key(), p.height as i64))
            .collect();
        assert_eq!(observed.join(","), "g:2,b:3,b:3,b:3,s:2");
    }

    #[test]
    fn full_layout_ord_matches_r() {
        // Oracle row: n=2 lord=cpbhg with every optional panel on
        // => "c:2,p:2,b:3,b:3,h:1,g:2,s:2"
        let mut request = base_request();
        request.has_peaks = true;
        request.has_chromhmm = true;
        request.has_cytoband = true;
        request.layout_ord = vec![
            TrackKind::Cytoband,
            TrackKind::Peaks,
            TrackKind::BigWig,
            TrackKind::ChromHmm,
            TrackKind::Gene,
        ];
        let layout = make_layout(&request);
        let observed: Vec<String> = layout
            .panels
            .iter()
            .map(|p| format!("{}:{}", p.kind.key(), p.height as i64))
            .collect();
        assert_eq!(observed.join(","), "c:2,p:2,b:3,b:3,h:1,g:2,s:2");
    }

    #[test]
    fn bigwig_only_request_matches_r() {
        // Oracle row: n=1 lord=b, no gene/ideogram => "b:3,s:2"
        let mut request = base_request();
        request.bigwig_count = 1;
        request.has_gene = false;
        request.layout_ord = vec![TrackKind::BigWig];
        let layout = make_layout(&request);
        let observed: Vec<String> = layout
            .panels
            .iter()
            .map(|p| format!("{}:{}", p.kind.key(), p.height as i64))
            .collect();
        assert_eq!(observed.join(","), "b:3,s:2");
    }

    #[test]
    fn heights_follow_arguments() {
        let layout = make_layout(&base_request());
        // 2 bigWigs at 3 + gene 2 + scale 2
        assert_eq!(layout.total_height, 3.0 + 3.0 + 2.0 + 2.0);

        // Look panels up by kind rather than index, since the scale panel comes
        // first when `layout_ord` is empty.
        let bigwig_heights: Vec<f64> = layout
            .panels
            .iter()
            .filter(|p| p.kind == TrackKind::BigWig)
            .map(|p| p.height)
            .collect();
        assert_eq!(bigwig_heights, vec![3.0, 3.0]);

        let gene_height = layout
            .panels
            .iter()
            .find(|p| p.kind == TrackKind::Gene)
            .map(|p| p.height);
        assert_eq!(gene_height, Some(2.0));

        let scale_height = layout
            .panels
            .iter()
            .find(|p| p.kind == TrackKind::Scale)
            .map(|p| p.height);
        assert_eq!(scale_height, Some(2.0));
    }

    #[test]
    fn scale_is_always_appended_even_when_omitted() {
        let mut request = base_request();
        request.layout_ord = vec![TrackKind::Gene, TrackKind::BigWig];
        let layout = make_layout(&request);
        let kinds: Vec<TrackKind> = layout.panels.iter().map(|p| p.kind).collect();
        assert!(
            kinds.contains(&TrackKind::Scale),
            "scale panel must always be present: {kinds:?}"
        );
        // User order is honoured for the keys they named.
        assert_eq!(kinds[0], TrackKind::Gene);
        assert_eq!(kinds[1], TrackKind::BigWig);
    }

    #[test]
    fn draw_indices_are_a_permutation() {
        let mut request = base_request();
        request.has_peaks = true;
        request.has_chromhmm = true;
        request.has_cytoband = true;
        request.layout_ord = vec![
            TrackKind::Cytoband,
            TrackKind::ChromHmm,
            TrackKind::Peaks,
            TrackKind::BigWig,
            TrackKind::Gene,
        ];
        let layout = make_layout(&request);

        let mut indices: Vec<usize> = layout.panels.iter().map(|p| p.draw_index).collect();
        indices.sort_unstable();
        let expected: Vec<usize> = (1..=layout.panels.len()).collect();
        assert_eq!(indices, expected, "draw indices must cover 1..=n exactly");
    }

    #[test]
    fn draw_order_matches_r_canonical_sequence() {
        // With an explicit user order, draw indices must still follow the fixed
        // p, b, h, g, s, c sequence regardless of the user's sequence.
        let mut request = base_request();
        request.has_peaks = true;
        request.layout_ord = vec![TrackKind::Peaks, TrackKind::BigWig, TrackKind::Gene];
        let layout = make_layout(&request);

        let peaks = layout
            .panels
            .iter()
            .find(|p| p.kind == TrackKind::Peaks)
            .expect("peaks panel");
        let gene = layout
            .panels
            .iter()
            .find(|p| p.kind == TrackKind::Gene)
            .expect("gene panel");
        let scale = layout
            .panels
            .iter()
            .find(|p| p.kind == TrackKind::Scale)
            .expect("scale panel");
        assert!(peaks.draw_index < gene.draw_index);
        assert!(gene.draw_index < scale.draw_index);
    }

    #[test]
    fn optional_panels_are_excluded_when_disabled() {
        let mut request = base_request();
        request.has_gene = false;
        request.has_cytoband = true;
        request.has_chromhmm = false;
        let layout = make_layout(&request);
        let kinds: Vec<TrackKind> = layout.panels.iter().map(|p| p.kind).collect();
        assert!(!kinds.contains(&TrackKind::Gene));
        assert!(!kinds.contains(&TrackKind::ChromHmm));
        assert!(kinds.contains(&TrackKind::Cytoband));
    }

    #[test]
    fn unknown_layout_letters_are_ignored() {
        assert_eq!(TrackKind::from_key('x'), None);
        assert_eq!(TrackKind::from_key('p'), Some(TrackKind::Peaks));
        assert_eq!(TrackKind::from_key('b'), Some(TrackKind::BigWig));
    }

    #[test]
    fn fractional_heights_sum_to_one() {
        let layout = make_layout(&base_request());
        let sum: f64 = layout.fractional_heights().iter().sum();
        assert!((sum - 1.0).abs() < 1e-12, "fractions summed to {sum}");
    }
}
