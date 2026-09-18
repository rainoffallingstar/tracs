//! Native track renderer.
//!
//! Replaces the previous `Rscript trackplot.R` step: the extraction side already
//! produced `tracks.tsv`/`meta.tsv` in Rust, and this module turns those into a
//! PDF (or SVG) without an R dependency.
//!
//! The layout follows `track_plot()`: a vertical stack of independently
//! y-scaled panels (bigWig signal, gene models, scale bar, peaks, chromHMM,
//! ideogram) whose order and relative heights come from `layout_ord`.

pub mod fonts;
pub mod io;
pub mod layout;
pub mod pretty;
pub mod svg;
