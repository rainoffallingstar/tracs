use std::collections::HashMap;
use std::ffi::OsString;
use std::fs;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::Command as ProcCommand;
use std::process::Stdio;
use std::time::Duration;

use anyhow::{anyhow, Context, Result};
use bigtools::bed::bedparser::parse_bed;
use bigtools::utils::misc::stats_for_bed_item;
use bigtools::BigWigRead;
use clap::{Parser, Subcommand};
use mysql::params;
use mysql::prelude::Queryable;
use serde_json::Value as JsonValue;

mod plot;

/// Base font size for track labels, in points.
///
/// R's `track_plot()` uses base graphics defaults with `cex` scaling; this is
/// the equivalent flat size for the native renderer.
const DEFAULT_FONT_SIZE: f64 = 10.0;

/// Default left margin for panels, leaving room for y-axis labels and track
/// names. `track_plot()` uses `4` or `2` lines depending on `show_axis`; at the
/// default font size one line is roughly 12pt, giving ~48 or ~24pt.
const DEFAULT_LEFT_MARGIN: f64 = 48.0;

#[derive(Parser, Debug)]
#[command(
    name = "tracs",
    about = "Rust replacements for bwtool commands used by trackplot.R (summary/matrix), plus optional higher-level helpers",
    version
)]
struct Cli {
    #[command(subcommand)]
    cmd: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Summarize values in a bigWig over intervals from a BED file.
    Summary(SummaryArgs),
    /// Generate a matrix of tiled averages around interval anchors.
    Matrix(MatrixArgs),
    /// Extract binned max signal tracks for a locus/gene across multiple bigWigs.
    TrackExtract(TrackExtractArgs),
    /// Run `track-extract` then render tracks natively into a PDF.
    #[command(alias = "plot")]
    PlotTrack(PlotTrackArgs),
    /// Draw a profile plot: mean/median signal around a focal point per sample.
    Profile(ProfileArgs),
    /// Draw a heatmap of matrices around a focal point, one panel per sample.
    Heatmap(HeatmapArgs),
    /// Draw a PCA scatter plot of samples from summary tables.
    Pca(PcaArgs),
    /// Draw a volcano plot from a differential-binding results table.
    Volcano(VolcanoArgs),
    /// Summarize HOMER `annotatePeaks.pl` output as stacked annotation bars.
    HomerAnnots(HomerAnnotsArgs),
}

#[derive(Parser, Debug)]
struct HomerAnnotsArgs {
    /// One or more HOMER `annotatePeaks.pl` output files (repeatable), one per
    /// sample. The `Annotation` column is selected by name, so the usual default
    /// column set works as-is.
    #[arg(long = "anno")]
    anno: Vec<PathBuf>,

    /// Sample names, comma-separated and matching `--anno` order. Defaults to
    /// each file's name up to its first dot, as `summarize_homer_annots()` does.
    #[arg(long = "sample")]
    samples: Vec<String>,

    /// Keep HOMER's literal `NA` annotation as a real category.
    ///
    /// Off by default because R drops those peaks: its palette names a colour
    /// for `NA`, but `fread()` reads the text as a missing value and `%in%`
    /// never matches it, so unannotated peaks vanish while still shrinking every
    /// other segment. Enabling this draws them in the palette's `gray70` and
    /// makes each bar sum to 1.
    #[arg(long = "keep-unannotated", default_value_t = false)]
    keep_unannotated: bool,

    /// Legend font size multiplier, matching `summarize_homer_annots()`'s
    /// `legend_font_size`.
    #[arg(long = "legend-font-size", default_value_t = 1.0)]
    legend_font_size: f64,

    /// Output path; `.svg` or `.pdf` decides the format.
    #[arg(long = "out")]
    out: PathBuf,

    /// Output working directory for intermediate files (optional).
    #[arg(long = "work-dir")]
    work_dir: Option<PathBuf>,

    /// Figure width in inches.
    #[arg(long = "width", default_value_t = 6.0)]
    width: f64,

    /// Figure height in inches.
    #[arg(long = "height", default_value_t = 4.0)]
    height: f64,

    /// Draw axis ticks and labels.
    #[arg(long = "show-axis", default_value_t = true, action = clap::ArgAction::Set)]
    show_axis: bool,
}

#[derive(Parser, Debug)]
struct VolcanoArgs {
    /// Differential-binding results table, tab- or comma-separated, with a
    /// header row. Must carry a log fold change, a p-value and an adjusted
    /// p-value column.
    ///
    /// Works with `limma::topTable()`, `DESeq2::results()`, `edgeR::topTags()`,
    /// or any table with those three columns.
    #[arg(long = "results")]
    results: PathBuf,

    /// Column holding the log fold change (aliases: `logFC`, `log2FoldChange`).
    #[arg(long = "logfc-col")]
    log_fc_col: Option<String>,

    /// Column holding the raw p-value (aliases: `P.Value`, `pvalue`).
    #[arg(long = "p-col")]
    p_col: Option<String>,

    /// Column holding the adjusted p-value (aliases: `adj.P.Val`, `padj`).
    #[arg(long = "padj-col")]
    padj_col: Option<String>,

    /// Significance threshold on the adjusted p-value, matching
    /// `volcano_plot(fdr = ...)`. A peak is significant when `padj < fdr`.
    #[arg(long = "fdr", default_value_t = 0.1)]
    fdr: f64,

    /// Colour for significantly up peaks.
    #[arg(long = "upcol", default_value = plot::volcano::DEFAULT_UP_COLOR)]
    upcol: String,

    /// Colour for significantly down peaks.
    #[arg(long = "downcol", default_value = plot::volcano::DEFAULT_DOWN_COLOR)]
    downcol: String,

    /// Point opacity, matching `volcano_plot()`'s `alpha` (R's `alpha.f`).
    #[arg(long = "alpha", default_value_t = 0.6)]
    alpha: f64,

    /// Point size multiplier, matching `volcano_plot()`'s `size`.
    #[arg(long = "point-size", default_value_t = 0.8)]
    point_size: f64,

    /// Plot title. Overrides the `contrast` attribute read from the table's
    /// comment header, if any.
    #[arg(long = "title")]
    title: Option<String>,

    /// Output path; `.svg` or `.pdf` decides the format.
    #[arg(long = "out")]
    out: PathBuf,

    /// Output working directory for intermediate files (optional).
    #[arg(long = "work-dir")]
    work_dir: Option<PathBuf>,

    /// Figure width in inches.
    #[arg(long = "width", default_value_t = 6.0)]
    width: f64,

    /// Figure height in inches.
    #[arg(long = "height", default_value_t = 6.0)]
    height: f64,

    /// Draw axis ticks and labels.
    #[arg(long = "show-axis", default_value_t = true, action = clap::ArgAction::Set)]
    show_axis: bool,
}

#[derive(Parser, Debug)]
struct HeatmapArgs {
    /// One or more matrix files from `tracs matrix` (repeatable).
    #[arg(long = "matrix")]
    matrices: Vec<PathBuf>,

    /// Sample names, comma-separated and matching `--matrix` order.
    /// Defaults to each matrix's file stem.
    #[arg(long = "sample")]
    samples: Vec<String>,

    /// Bases upstream of the focal point (must match how the matrices were built).
    #[arg(long, default_value_t = 2500)]
    up: u32,

    /// Bases downstream of the focal point (must match how the matrices were built).
    #[arg(long, default_value_t = 2500)]
    down: u32,

    /// Row ordering within each panel: `mean` or `median`.
    #[arg(long = "sort-by", default_value = "mean")]
    sort_by: String,

    /// Sequential palette name (e.g. `Blues`, `Viridis`, `Greys`, `Reds`).
    #[arg(long = "col-pal", default_value = "Blues")]
    col_pal: String,

    /// Reverse the palette (light where dark was).
    #[arg(long = "revpal", default_value_t = false)]
    revpal: bool,

    /// Lower colour limit, comma-separated per sample. Defaults to the matrix min.
    #[arg(long = "zmin")]
    zmin: Option<String>,

    /// Upper colour limit, comma-separated per sample. Defaults to the max row mean,
    /// matching `profile_heatmap()`.
    #[arg(long = "zmax")]
    zmax: Option<String>,

    /// Output path; `.svg` or `.pdf` decides the format.
    #[arg(long = "out")]
    out: PathBuf,

    /// Output working directory for intermediate files (optional).
    #[arg(long = "work-dir")]
    work_dir: Option<PathBuf>,

    /// Figure width in inches.
    #[arg(long = "width", default_value_t = 6.0)]
    width: f64,

    /// Height of each heatmap panel in inches.
    #[arg(long = "height", default_value_t = 3.0)]
    height: f64,

    /// Draw colour bars and axis labels.
    #[arg(long = "show-axis", default_value_t = true, action = clap::ArgAction::Set)]
    show_axis: bool,
}

#[derive(Parser, Debug)]
struct PcaArgs {
    /// One or more summary tables, each in `extract_summary()` orientation:
    /// rows are regions, columns are samples. Repeatable.
    ///
    /// A `tracs matrix` file is not a summary table, because its columns are
    /// bins rather than samples. Build summary tables with
    /// `tracs summary -with-sum`, one file per sample, then pass them here.
    #[arg(long = "summary")]
    summaries: Vec<PathBuf>,

    /// Sample names, comma-separated and matching the summary table's column order.
    /// Defaults to each table's file stem.
    #[arg(long = "sample")]
    samples: Vec<String>,

    /// Group label per sample, comma-separated and matching `--sample` order.
    /// Points are coloured by group, as `pca_plot(color_by = ...)` does.
    #[arg(long = "color-by")]
    color_by: Option<String>,

    /// Point colours, comma-separated. Defaults to `pca_plot()`'s palette.
    #[arg(long = "col")]
    col: Option<String>,

    /// Number of most-variable regions to keep, matching `pca_plot()`'s `top`.
    #[arg(long = "top", default_value_t = 1000)]
    top: usize,

    /// Apply `log2(x + offset)` before the PCA, matching `pca_plot(log2 = TRUE)`.
    #[arg(long = "log2", default_value_t = false)]
    log2: bool,

    /// Offset used by `--log2`, matching `profile_plot()`'s default.
    #[arg(long = "log2-offset", default_value_t = 0.1)]
    log2_offset: f64,

    /// Component drawn on the x axis, 1-based (`pca_plot()`'s `xpc`).
    #[arg(long = "xpc", default_value_t = 1)]
    xpc: usize,

    /// Component drawn on the y axis, 1-based (`pca_plot()`'s `ypc`).
    #[arg(long = "ypc", default_value_t = 2)]
    ypc: usize,

    /// Flip the sign of the x component. R's component signs are arbitrary and
    /// can differ between builds, so this matches a specific reference figure.
    #[arg(long = "flip-x", default_value_t = false)]
    flip_x: bool,

    /// Flip the sign of the y component.
    #[arg(long = "flip-y", default_value_t = false)]
    flip_y: bool,

    /// Draw the variance-explained scree panel beside the scatter plot
    /// (`pca_plot(show_cree = TRUE)`).
    #[arg(long = "show-cree", default_value_t = true, action = clap::ArgAction::Set)]
    show_cree: bool,

    /// Sample label size multiplier, matching `pca_plot()`'s `lab_size`.
    /// Use 0 to hide the labels.
    #[arg(long = "lab-size", default_value_t = 1.0)]
    lab_size: f64,

    /// Point size multiplier, matching `pca_plot()`'s `size`.
    #[arg(long = "point-size", default_value_t = 1.0)]
    point_size: f64,

    /// Output path; `.svg` or `.pdf` decides the format.
    #[arg(long = "out")]
    out: PathBuf,

    /// Output working directory for intermediate files (optional).
    #[arg(long = "work-dir")]
    work_dir: Option<PathBuf>,

    /// Figure width in inches.
    #[arg(long = "width", default_value_t = 6.0)]
    width: f64,

    /// Figure height in inches.
    #[arg(long = "height", default_value_t = 5.0)]
    height: f64,

    /// Draw axis ticks and labels.
    #[arg(long = "show-axis", default_value_t = true, action = clap::ArgAction::Set)]
    show_axis: bool,
}

#[derive(Parser, Debug)]
struct ProfileArgs {
    /// One or more matrix files from `tracs matrix` (repeatable).
    #[arg(long = "matrix")]
    matrices: Vec<PathBuf>,

    /// Sample names, comma-separated and matching `--matrix` order.
    /// Defaults to each matrix's file stem.
    #[arg(long = "sample")]
    samples: Vec<String>,

    /// Bases upstream of the focal point (must match how the matrices were built).
    #[arg(long, default_value_t = 2500)]
    up: u32,

    /// Bases downstream of the focal point (must match how the matrices were built).
    #[arg(long, default_value_t = 2500)]
    down: u32,

    /// How replicates are collapsed into one line: `mean` or `median`.
    #[arg(long = "stat", default_value = "mean")]
    stat: String,

    /// Group labels, comma-separated and matching `--matrix` order. When given,
    /// samples sharing a label are pooled into one line (R's `condition`).
    #[arg(long = "condition")]
    condition: Option<String>,

    /// Line colours, comma-separated. Defaults to `profile_plot()`'s palette.
    #[arg(long = "col")]
    col: Option<String>,

    /// Output path; `.svg` or `.pdf` decides the format.
    #[arg(long = "out")]
    out: PathBuf,

    /// Output working directory for intermediate files (optional).
    #[arg(long = "work-dir")]
    work_dir: Option<PathBuf>,

    /// Figure width in inches.
    #[arg(long = "width", default_value_t = 6.0)]
    width: f64,

    /// Figure height in inches.
    #[arg(long = "height", default_value_t = 4.0)]
    height: f64,

    /// x axis label.
    #[arg(long)]
    xlab: Option<String>,

    /// y axis label.
    #[arg(long)]
    ylab: Option<String>,

    /// Draw axis ticks and labels.
    #[arg(long = "show-axis", default_value_t = true, action = clap::ArgAction::Set)]
    show_axis: bool,
}

#[derive(Parser, Debug)]
struct SummaryArgs {
    /// Include a sum column (kept for bwtool CLI compatibility).
    #[arg(long = "with-sum", default_value_t = false)]
    with_sum: bool,

    /// Keep bed columns in output (we always output chr/start/end/size).
    #[arg(long = "keep-bed", default_value_t = false)]
    keep_bed: bool,

    /// Print header.
    #[arg(long = "header", default_value_t = false)]
    header: bool,

    /// Input BED path (3+ columns). Uses 0-based, half-open coordinates.
    bed: PathBuf,

    /// Input bigWig path.
    bigwig: PathBuf,

    /// Output path.
    out: PathBuf,
}

#[derive(Parser, Debug)]
struct MatrixArgs {
    /// Anchor at interval start.
    #[arg(long = "starts", default_value_t = false)]
    starts: bool,

    /// Anchor at interval end.
    #[arg(long = "ends", default_value_t = false)]
    ends: bool,

    /// Bin size for tiled averages (required by trackplot.R usage).
    #[arg(long = "tiled-averages")]
    tiled_averages: u32,

    /// Region size, formatted as UP:DOWN (e.g. 2500:2500).
    size: String,

    /// Input BED path (3+ columns).
    bed: PathBuf,

    /// Input bigWig path.
    bigwig: PathBuf,

    /// Output path.
    out: PathBuf,
}

#[derive(Parser, Debug)]
struct TrackExtractArgs {
    /// Output directory to write results into.
    #[arg(long = "out-dir")]
    out_dir: PathBuf,

    /// Bin size (bp) to compute max signal over.
    #[arg(long, default_value_t = 10)]
    binsize: u32,

    /// Padding (bp) to extend both sides of the locus/gene.
    #[arg(long, default_value_t = 0)]
    padding: i64,

    /// Target region, formatted as chr:start-end (commas allowed).
    #[arg(long)]
    loci: Option<String>,

    /// Gene query to extract. Accepted forms:
    /// - gene symbol (e.g. `CD1D`)
    /// - Entrez ID (e.g. `912`)
    /// - Ensembl gene id (e.g. `ENSG00000158473`)
    #[arg(long)]
    gene: Option<String>,

    /// Reference genome build used for UCSC refGene lookup when using --gene (if --gtf isn't provided or doesn't match).
    #[arg(long, default_value = "hg19")]
    build: String,

    /// Cytoband table name for ideogram (UCSC), e.g. cytoBand.
    #[arg(long = "ideo-tbl", default_value = "cytoBand")]
    ideo_tbl: String,

    /// Skip UCSC cytoband fetch/output.
    #[arg(long = "no-cytoband", default_value_t = false)]
    no_cytoband: bool,

    /// Skip gene model extraction/output (gene_models.tsv).
    #[arg(long = "no-gene-models", default_value_t = false)]
    no_gene_models: bool,

    /// GTF file for gene model lookup when using --gene (optional).
    #[arg(long)]
    gtf: Option<PathBuf>,

    /// Input bigWig paths (repeatable).
    #[arg(long = "bigwig")]
    bigwigs: Vec<PathBuf>,

    /// Optional sample names (repeatable, same length as --bigwig). Defaults to file basenames.
    #[arg(long = "sample")]
    samples: Vec<String>,
}

#[derive(Parser, Debug)]
struct PlotTrackArgs {
    /// Output PDF path.
    #[arg(long = "out")]
    out_pdf: PathBuf,

    /// Output working directory for intermediate TSVs (optional). If not provided, a temp dir is created.
    #[arg(long = "work-dir")]
    work_dir: Option<PathBuf>,

    /// PDF width in inches.
    #[arg(long = "pdf-width", default_value_t = 12.0)]
    pdf_width: f64,

    /// PDF height in inches.
    #[arg(long = "pdf-height", default_value_t = 6.0)]
    pdf_height: f64,

    /// Bin size (bp) to compute max signal over.
    #[arg(long, default_value_t = 10)]
    binsize: u32,

    /// Padding (bp) to extend both sides of the locus/gene.
    #[arg(long, default_value_t = 0)]
    padding: i64,

    /// Target region, formatted as chr:start-end (commas allowed).
    #[arg(long)]
    loci: Option<String>,

    /// Gene query to extract. Accepted forms:
    /// - gene symbol (e.g. `CD1D`)
    /// - Entrez ID (e.g. `912`)
    /// - Ensembl gene id (e.g. `ENSG00000158473`)
    #[arg(long)]
    gene: Option<String>,

    /// Reference genome build used for UCSC refGene lookup when using --gene/--loci.
    #[arg(long, default_value = "hg19")]
    build: String,

    /// Cytoband table name for ideogram (UCSC), e.g. cytoBand.
    #[arg(long = "ideo-tbl", default_value = "cytoBand")]
    ideo_tbl: String,

    /// Show ideogram track.
    #[arg(long = "show-ideogram", default_value_t = true, action = clap::ArgAction::Set)]
    show_ideogram: bool,

    /// Draw gene track.
    #[arg(long = "draw-gene-track", default_value_t = true, action = clap::ArgAction::Set)]
    draw_gene_track: bool,

    /// Draw all bigWigs in a single overlay track as line plot.
    #[arg(long = "track-overlay", default_value_t = false)]
    track_overlay: bool,

    /// Collapse transcripts by gene.
    #[arg(long = "collapse-txs", default_value_t = true, action = clap::ArgAction::Set)]
    collapse_txs: bool,

    /// Optional GTF for offline gene lookup (best with Ensembl IDs).
    #[arg(long)]
    gtf: Option<PathBuf>,

    /// Optional coldata TSV with columns `bw_files` and `bw_sample_names` (same shape as `read_coldata()` output).
    /// If provided, `--bigwig/--sample` are ignored.
    #[arg(long = "coldata")]
    coldata: Option<PathBuf>,

    /// Input bigWig paths (repeatable). Ignored when `--coldata` is provided.
    #[arg(long = "bigwig")]
    bigwigs: Vec<PathBuf>,

    /// Optional sample names (repeatable, same length as --bigwig). Defaults to file basenames. Ignored when `--coldata` is provided.
    #[arg(long = "sample")]
    samples: Vec<String>,

    /// Colors for tracks, comma-separated (e.g. "#d34,#2980b9").
    /// Use `auto` to assign a distinct palette automatically.
    #[arg(long = "col", default_value = "auto")]
    col: String,

    /// Auto-scale y-range across all tracks (track_plot `groupAutoScale`).
    #[arg(long = "group-auto-scale", default_value_t = false, action = clap::ArgAction::Set)]
    group_auto_scale: bool,

    /// Custom y-max, comma-separated (length 1 or N tracks).
    #[arg(long = "y-max")]
    y_max: Option<String>,

    /// Custom y-min, comma-separated (length 1 or N tracks).
    #[arg(long = "y-min")]
    y_min: Option<String>,

    /// Only show these transcript IDs in the gene track (comma-separated; track_plot `txname`).
    #[arg(long = "txname")]
    txname: Option<String>,

    /// Only show these gene names in the gene track (comma-separated; track_plot `genename`).
    #[arg(long = "genename")]
    genename: Option<String>,

    /// Show y-axis scale.
    #[arg(long = "show-axis", default_value_t = false, action = clap::ArgAction::Set)]
    show_axis: bool,

    /// Track display names, comma-separated (track_plot `track_names`).
    #[arg(long = "track-names")]
    track_names: Option<String>,

    /// Track name x-position (track_plot `track_names_pos`).
    #[arg(long = "track-names-pos", default_value_t = 0.0)]
    track_names_pos: f64,

    /// Place track names on the left (track_plot `track_names_to_left`).
    #[arg(
        long = "track-names-to-left",
        default_value_t = false,
        action = clap::ArgAction::Set
    )]
    track_names_to_left: bool,

    /// Gene track font scale (track_plot `gene_fsize`).
    #[arg(long = "gene-fsize", default_value_t = 1.0)]
    gene_fsize: f64,

    /// Order bigWig tracks by sample name (comma-separated; track_plot `bw_ord`).
    #[arg(long = "bw-ord")]
    bw_ord: Option<String>,

    /// Layout order (comma-separated; track_plot `layout_ord`, default "p,b,h,g,c").
    #[arg(long = "layout-ord", default_value = "p,b,h,g,c")]
    layout_ord: String,

    /// Regions BED/TSV to mark (first 3 cols: chr, start, end).
    #[arg(long = "regions-bed")]
    regions_bed: Option<PathBuf>,

    /// bw track height (track_plot `bw_track_height`).
    #[arg(long = "bw-track-height", default_value_t = 3.0)]
    bw_track_height: f64,

    /// peaks track height (track_plot `peaks_track_height`).
    #[arg(long = "peaks-track-height", default_value_t = 2.0)]
    peaks_track_height: f64,

    /// gene track height (track_plot `gene_track_height`).
    #[arg(long = "gene-track-height", default_value_t = 2.0)]
    gene_track_height: f64,

    /// scale track height (track_plot `scale_track_height`).
    #[arg(long = "scale-track-height", default_value_t = 2.0)]
    scale_track_height: f64,

    /// chromHMM track height (track_plot `chromHMM_track_height`).
    #[arg(long = "chromhmm-track-height", default_value_t = 1.0)]
    chromhmm_track_height: f64,

    /// cytoband/ideogram track height (track_plot `cytoband_track_height`).
    #[arg(long = "cytoband-track-height", default_value_t = 2.0)]
    cytoband_track_height: f64,

    /// Left margin (track_plot `left_mar`). Default NULL (trackplot decides based on `show_axis`).
    #[arg(long = "left-mar")]
    left_mar: Option<f64>,

    /// Peaks BED files to draw as top peak tracks (repeatable; track_plot `peaks`).
    #[arg(long = "peaks")]
    peaks: Vec<PathBuf>,

    /// Peaks track names (comma-separated; track_plot `peaks_track_names`).
    #[arg(long = "peaks-track-names")]
    peaks_track_names: Option<String>,

    /// chromHMM BED files with 4 columns (chr,start,end,name; repeatable; track_plot `chromHMM`).
    #[arg(long = "chromhmm")]
    chromhmm: Vec<PathBuf>,

    /// chromHMM track names (comma-separated; track_plot `chromHMM_names`).
    #[arg(long = "chromhmm-names")]
    chromhmm_names: Option<String>,

    /// chromHMM color mapping as `state=color` pairs (comma-separated; track_plot `chromHMM_cols`).
    /// Example: `1=red,2=orange,3=purple`.
    #[arg(long = "chromhmm-cols")]
    chromhmm_cols: Option<String>,

    /// Fetch chromHMM tracks from UCSC tables (repeatable table names). Implemented in Rust (no `mysql` binary).
    #[arg(long = "ucsc-chromhmm")]
    ucsc_chromhmm: Vec<String>,

    /// Highlight box color (track_plot `boxcol`).
    #[arg(long = "boxcol", default_value = "#ffc41a")]
    boxcol: String,

    /// Highlight box alpha (track_plot `boxcolalpha`).
    #[arg(long = "boxcolalpha", default_value_t = 0.4)]
    boxcolalpha: f64,
}

fn main() -> Result<()> {
    let cli = Cli::parse_from(normalize_bwtoolish_args(std::env::args_os()));
    match cli.cmd {
        Command::Summary(args) => cmd_summary(args),
        Command::Matrix(args) => cmd_matrix(args),
        Command::TrackExtract(args) => cmd_track_extract(args),
        Command::PlotTrack(args) => cmd_plot_track(args),
        Command::Profile(args) => cmd_profile(args),
        Command::Heatmap(args) => cmd_heatmap(args),
        Command::Pca(args) => cmd_pca(args),
        Command::Volcano(args) => cmd_volcano(args),
        Command::HomerAnnots(args) => cmd_homer_annots(args),
    }
}

fn normalize_bwtoolish_args<I>(args: I) -> Vec<OsString>
where
    I: IntoIterator<Item = OsString>,
{
    // bwtool uses single-dash long options (e.g. `-with-sum`, `-tiled-averages=50`).
    // clap expects `--with-sum`. Normalize for compatibility.
    let mut out: Vec<OsString> = Vec::new();
    for a in args.into_iter() {
        let s = a.to_string_lossy();
        let replaced = match s.as_ref() {
            "-with-sum" => "--with-sum".into(),
            "-keep-bed" => "--keep-bed".into(),
            "-header" => "--header".into(),
            "-starts" => "--starts".into(),
            "-ends" => "--ends".into(),
            _ => {
                if let Some(rest) = s.strip_prefix("-tiled-averages=") {
                    format!("--tiled-averages={rest}").into()
                } else {
                    a
                }
            }
        };
        out.push(replaced);
    }
    out
}

fn cmd_summary(args: SummaryArgs) -> Result<()> {
    // Keep flag compatibility even if we don't use them beyond validation.
    let _ = args.keep_bed;
    let _ = args.with_sum;

    let bed_file = File::open(&args.bed).with_context(|| format!("open bed: {:?}", args.bed))?;
    let bigwig_file =
        BigWigRead::open_file(&args.bigwig).with_context(|| format!("open bigwig: {:?}", args.bigwig))?;
    let mut bigwig = bigwig_file.cached();

    let chrom_lens = chrom_lengths(&bigwig);

    let out_file = File::create(&args.out).with_context(|| format!("create out: {:?}", args.out))?;
    let mut w = BufWriter::new(out_file);

    if args.header {
        // Keep the specific column names trackplot.R expects (`start`, `end`, `size`, `max`).
        // Extra columns are fine; trackplot.R selects a subset.
        writeln!(&mut w, "chromosome\tstart\tend\tsize\tsum\tmin\tmax\tmean")?;
    }

    let reader = BufReader::new(bed_file);
    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let Some(parsed) = parse_bed(&line) else {
            continue;
        };
        let (chrom, entry) = parsed?;

        let start = entry.start;
        let end = entry.end;
        if end <= start {
            continue;
        }

        let chr_len = chrom_lens.get(chrom).copied().unwrap_or(0);
        let query_start = start.min(chr_len);
        let query_end = end.min(chr_len);

        let mut stats = if chr_len == 0 || query_end <= query_start {
            None
        } else {
            Some(stats_for_bed_item(
                chrom,
                bigtools::BedEntry {
                    start: query_start,
                    end: query_end,
                    rest: entry.rest.clone(),
                },
                &mut bigwig,
            )?)
        };

        // Convert no-signal cases to zeros (trackplot expects numeric).
        let (sum, min, max) = match stats.as_mut() {
            None => (0.0, 0.0, 0.0),
            Some(s) if s.bases == 0 => (0.0, 0.0, 0.0),
            Some(s) => (s.sum, s.min, s.max),
        };

        let size = end - start;
        let mean = if size == 0 { 0.0 } else { sum / (size as f64) };
        writeln!(
            &mut w,
            "{chrom}\t{start}\t{end}\t{size}\t{sum}\t{min}\t{max}\t{mean}"
        )?;
    }

    w.flush()?;
    Ok(())
}

/// Reads a `tracs matrix` output file into rows of numbers.
///
/// Matrix files are headerless and whitespace/tab separated, one row per region.
/// Non-numeric tokens yield `NaN` so a malformed cell degrades to "missing"
/// rather than aborting the whole plot.
fn read_matrix_file(path: &Path) -> Result<Vec<Vec<f64>>> {
    let text = fs::read_to_string(path).with_context(|| format!("read matrix: {path:?}"))?;
    let mut rows = Vec::new();
    for line in text.lines() {
        if line.trim().is_empty() {
            continue;
        }
        let row: Vec<f64> = line
            .split_whitespace()
            .map(|token| token.parse::<f64>().unwrap_or(f64::NAN))
            .collect();
        if !row.is_empty() {
            rows.push(row);
        }
    }
    if rows.is_empty() {
        return Err(anyhow!("matrix file has no numeric rows: {path:?}"));
    }
    Ok(rows)
}

fn cmd_profile(args: ProfileArgs) -> Result<()> {
    if args.matrices.is_empty() {
        return Err(anyhow!("profile: at least one --matrix is required"));
    }
    if args.up == 0 && args.down == 0 {
        return Err(anyhow!("profile: --up and --down cannot both be 0"));
    }

    let stat = plot::profile::SummaryStat::from_name(&args.stat)
        .ok_or_else(|| anyhow!("profile: --stat must be 'mean' or 'median' (got {:?})", args.stat))?;

    // Sample names come from --sample when given, else each matrix's file stem.
    let names: Vec<String> = split_csv(&args.samples.join(","));
    let samples: Vec<(String, Vec<Vec<f64>>)> = args
        .matrices
        .iter()
        .enumerate()
        .map(|(index, path)| {
            let name = names.get(index).cloned().unwrap_or_else(|| {
                path.file_stem()
                    .and_then(|value| value.to_str())
                    .unwrap_or("sample")
                    .to_string()
            });
            Ok((name, read_matrix_file(path)?))
        })
        .collect::<Result<Vec<_>>>()?;

    let conditions = args
        .condition
        .as_ref()
        .map(|raw| split_csv(raw))
        .filter(|values| !values.is_empty());
    if let Some(conditions) = conditions.as_ref() {
        if conditions.len() != samples.len() {
            return Err(anyhow!(
                "profile: --condition count ({}) must match --matrix count ({})",
                conditions.len(),
                samples.len()
            ));
        }
    }

    // With `--condition`, replicates are pooled per group, matching
    // profile_summarize(condition = ...).
    let summarized = match conditions.as_ref() {
        Some(conditions) => plot::profile::summarize_by_condition(&samples, conditions, stat),
        None => plot::profile::summarize_samples(&samples, stat),
    };
    let series: Vec<plot::profile::ProfileSeries> = summarized
        .into_iter()
        .map(|(name, values)| plot::profile::ProfileSeries { name, values })
        .collect();

    let requested_colors = args.col.as_ref().map(|raw| split_csv(raw)).unwrap_or_default();
    let colors = plot::profile::resolve_profile_colors(&requested_colors, series.len());

    let width = args.width * 72.0;
    let height = args.height * 72.0;
    let mut writer = plot::svg::SvgWriter::new(width, height);
    writer.panel(0.0, 0.0, width, height, |panel| {
        plot::profile::draw_profile_panel(
            panel,
            &series,
            &colors,
            args.up,
            args.down,
            args.show_axis,
            DEFAULT_FONT_SIZE,
            args.xlab.as_deref().unwrap_or(""),
            args.ylab.as_deref().unwrap_or(""),
        );
    });
    let svg_document = writer.finish();

    match plot::render::OutputFormat::from_path(&args.out) {
        plot::render::OutputFormat::Svg => {
            fs::write(&args.out, svg_document)
                .with_context(|| format!("write SVG: {:?}", args.out))?;
        }
        plot::render::OutputFormat::Pdf => {
            let pdf = plot::render::svg_to_pdf(&svg_document)?;
            fs::write(&args.out, pdf).with_context(|| format!("write PDF: {:?}", args.out))?;
        }
    }

    // Keep the summarized profile alongside the figure when a work dir is given,
    // so callers can plot from the numbers instead of re-deriving them.
    if let Some(work_dir) = args.work_dir.as_ref() {
        fs::create_dir_all(work_dir)
            .with_context(|| format!("create work dir: {work_dir:?}"))?;
        let mut out = BufWriter::new(File::create(work_dir.join("profile_summary.tsv"))?);
        write!(&mut out, "bin")?;
        for s in &series {
            write!(&mut out, "\t{}", sanitize_tsv_value(&s.name))?;
        }
        writeln!(&mut out)?;
        let nbins = series.iter().map(|s| s.values.len()).max().unwrap_or(0);
        for bin in 0..nbins {
            write!(&mut out, "{bin}")?;
            for s in &series {
                match s.values.get(bin) {
                    Some(value) if value.is_finite() => write!(&mut out, "\t{value}")?,
                    _ => write!(&mut out, "\tNA")?,
                }
            }
            writeln!(&mut out)?;
        }
        out.flush()?;
    }

    Ok(())
}

/// Parses a comma-separated list of per-sample numbers.
fn parse_csv_numbers(raw: Option<&String>) -> Option<Vec<f64>> {
    raw.map(|text| {
        split_csv(text)
            .iter()
            .filter_map(|value| value.trim().parse::<f64>().ok())
            .collect()
    })
    .filter(|values: &Vec<f64>| !values.is_empty())
}

fn cmd_heatmap(args: HeatmapArgs) -> Result<()> {
    if args.matrices.is_empty() {
        return Err(anyhow!("heatmap: at least one --matrix is required"));
    }
    if args.up == 0 && args.down == 0 {
        return Err(anyhow!("heatmap: --up and --down cannot both be 0"));
    }

    let sort_by = plot::heatmap::SortBy::from_name(&args.sort_by).ok_or_else(|| {
        anyhow!(
            "heatmap: --sort-by must be 'mean' or 'median' (got {:?})",
            args.sort_by
        )
    })?;

    let palette = plot::heatmap::colors::find_palette(&args.col_pal).ok_or_else(|| {
        anyhow!(
            "heatmap: unknown --col-pal {:?}; available: {}",
            args.col_pal,
            plot::heatmap::colors::palette_names().join(", ")
        )
    })?;

    let names = split_csv(&args.samples.join(","));
    let samples: Vec<(String, Vec<Vec<f64>>)> = args
        .matrices
        .iter()
        .enumerate()
        .map(|(index, path)| {
            let name = names.get(index).cloned().unwrap_or_else(|| {
                path.file_stem()
                    .and_then(|value| value.to_str())
                    .unwrap_or("sample")
                    .to_string()
            });
            Ok((name, read_matrix_file(path)?))
        })
        .collect::<Result<Vec<_>>>()?;

    let z_mins = parse_csv_numbers(args.zmin.as_ref());
    let z_maxs = parse_csv_numbers(args.zmax.as_ref());
    let panels = plot::heatmap::build_panels(
        &samples,
        sort_by,
        z_mins.as_deref(),
        z_maxs.as_deref(),
    );

    // R builds a 255-entry ramp; match that so the gradient resolves identically.
    let mut ramp = plot::heatmap::colors::resolve_ramp(palette, 255);
    if args.revpal {
        plot::heatmap::colors::reverse_ramp(&mut ramp);
    }

    // One panel per sample, stacked vertically like R's `layout()`.
    let panel_height = args.height * 72.0;
    let total_height = panel_height * panels.len() as f64;
    let width = args.width * 72.0;
    let mut writer = plot::svg::SvgWriter::new(width, total_height);
    for (index, panel) in panels.iter().enumerate() {
        writer.panel(0.0, index as f64 * panel_height, width, panel_height, |draw| {
            plot::heatmap::draw_heatmap_panel(
                draw,
                panel,
                &ramp,
                args.up,
                args.down,
                args.show_axis,
                DEFAULT_FONT_SIZE,
            );
        });
    }
    let svg_document = writer.finish();

    match plot::render::OutputFormat::from_path(&args.out) {
        plot::render::OutputFormat::Svg => {
            fs::write(&args.out, svg_document)
                .with_context(|| format!("write SVG: {:?}", args.out))?;
        }
        plot::render::OutputFormat::Pdf => {
            let pdf = plot::render::svg_to_pdf(&svg_document)?;
            fs::write(&args.out, pdf).with_context(|| format!("write PDF: {:?}", args.out))?;
        }
    }

    // Keep the resolved colour limits so a caller can reproduce or adjust them.
    if let Some(work_dir) = args.work_dir.as_ref() {
        fs::create_dir_all(work_dir)
            .with_context(|| format!("create work dir: {work_dir:?}"))?;
        let mut out = BufWriter::new(File::create(work_dir.join("heatmap_limits.tsv"))?);
        writeln!(&mut out, "sample\tz_min\tz_max\trows\tcolumns")?;
        for panel in &panels {
            let columns = panel.matrix.iter().map(Vec::len).max().unwrap_or(0);
            writeln!(
                &mut out,
                "{}\t{}\t{}\t{}\t{}",
                sanitize_tsv_value(&panel.name),
                panel.z_min,
                panel.z_max,
                panel.matrix.len(),
                columns
            )?;
        }
        out.flush()?;
    }

    Ok(())
}

/// Reads a `tracs summary` output file into one value column.
///
/// `extract_summary()` keeps only the `sum` column from each `bwtool summary`
/// output, so that is what is read here. The file has a header row in the
/// `bwtool summary -header -with-sum` form (`chromosome start end size sum ...`),
/// which is detected by its non-numeric second field rather than by position, so
/// headerless files work too.
///
/// Values are taken from the column named `sum` when a header is present, and
/// from the fifth column otherwise, matching `x[,.(chromosome, start, end, size, sum)]`.
fn read_summary_column(path: &Path) -> Result<Vec<f64>> {
    let text = fs::read_to_string(path).with_context(|| format!("read summary: {path:?}"))?;
    let mut values = Vec::new();
    let mut sum_index: Option<usize> = None;

    for (line_number, line) in text.lines().enumerate() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = trimmed.split_whitespace().collect();

        // A header row is any row whose `sum` field is not a number.
        if sum_index.is_none() && line_number < 4 {
            if let Some(position) = fields.iter().position(|field| *field == "sum") {
                sum_index = Some(position);
                continue;
            }
        }

        let position = sum_index.unwrap_or(4);
        let value = fields
            .get(position)
            .map(|field| field.parse::<f64>().unwrap_or(f64::NAN))
            .unwrap_or(f64::NAN);
        values.push(value);
    }

    if values.is_empty() {
        return Err(anyhow!("summary file has no data rows: {path:?}"));
    }
    Ok(values)
}

fn cmd_pca(args: PcaArgs) -> Result<()> {
    if args.summaries.is_empty() {
        return Err(anyhow!("pca: at least one --summary is required"));
    }
    if args.xpc == 0 || args.ypc == 0 {
        return Err(anyhow!("pca: --xpc/--ypc are 1-based and must be >= 1"));
    }
    if args.top == 0 {
        return Err(anyhow!("pca: --top must be >= 1"));
    }

    // Each summary file supplies one sample column. The sample order is the
    // order the files were given, matching how `extract_summary()` cbind's the
    // per-bigWig columns.
    let names = split_csv(&args.samples.join(","));
    let columns: Vec<(String, Vec<f64>)> = args
        .summaries
        .iter()
        .enumerate()
        .map(|(index, path)| {
            let name = names.get(index).cloned().unwrap_or_else(|| {
                path.file_stem()
                    .and_then(|value| value.to_str())
                    .unwrap_or("sample")
                    .to_string()
            });
            Ok((name, read_summary_column(path)?))
        })
        .collect::<Result<Vec<_>>>()?;

    let table = plot::pca::SummaryTable::from_sample_columns(&columns);
    if table.n_regions() == 0 {
        return Err(anyhow!("pca: summary tables contained no regions"));
    }
    if table.n_samples() < 2 {
        return Err(anyhow!(
            "pca: need at least two samples, got {}",
            table.n_samples()
        ));
    }

    let (pca, ranking) = plot::pca::fit_summary_table(
        &table,
        args.top,
        args.log2,
        args.log2_offset,
    )
    .ok_or_else(|| {
        anyhow!(
            "pca: cannot fit with {} regions across {} samples",
            table.n_regions(),
            table.n_samples()
        )
    })?;

    // `pca_plot()` indexes components by name, so a request beyond the fitted
    // count is a user error worth reporting rather than silently plotting zero.
    let x_component = args.xpc - 1;
    let y_component = args.ypc - 1;
    if x_component >= pca.components.len() || y_component >= pca.components.len() {
        return Err(anyhow!(
            "pca: asked for PC{} vs PC{} but only {} components exist \
             (prcomp returns min(samples, regions); with {} samples and {} regions used)",
            args.xpc,
            args.ypc,
            pca.components.len(),
            pca.n_samples,
            pca.n_features
        ));
    }

    // Colours: without --color-by everything is black, matching `pca_plot()`.
    // With a group label per sample, the distinct labels take palette colours in
    // first-seen order, which is what `condition_colors[...][.N, condition]` does.
    let group_labels = args
        .color_by
        .as_ref()
        .map(|raw| split_csv(raw))
        .unwrap_or_default();
    let legend_entries: Vec<(String, String)> = if group_labels.is_empty() {
        Vec::new()
    } else {
        let mut distinct: Vec<String> = Vec::new();
        for label in &group_labels {
            if !distinct.contains(label) {
                distinct.push(label.clone());
            }
        }
        let requested = split_csv(&args.col.clone().unwrap_or_default());
        let palette = plot::profile::resolve_profile_colors(&requested, distinct.len());
        distinct
            .into_iter()
            .zip(palette)
            .collect()
    };

    let point_colors: Vec<String> = if group_labels.is_empty() {
        vec![plot::pca::panel::DEFAULT_POINT_COLOR.to_string(); pca.n_samples]
    } else {
        (0..pca.n_samples)
            .map(|sample| {
                group_labels
                    .get(sample)
                    .and_then(|label| {
                        legend_entries
                            .iter()
                            .find(|(name, _)| name == label)
                            .map(|(_, color)| color.clone())
                    })
                    .unwrap_or_else(|| plot::pca::panel::DEFAULT_POINT_COLOR.to_string())
            })
            .collect()
    };

    let mut points = plot::pca::panel::scatter_points(
        &pca,
        x_component,
        y_component,
        &table.sample_names,
        &point_colors,
    );
    if args.flip_x {
        for point in points.iter_mut() {
            point.x = -point.x;
        }
    }
    if args.flip_y {
        for point in points.iter_mut() {
            point.y = -point.y;
        }
    }

    let x_title = plot::pca::panel::axis_title(
        &format!("PC{}", args.xpc),
        pca.variance_explained(x_component),
    );
    let y_title = plot::pca::panel::axis_title(
        &format!("PC{}", args.ypc),
        pca.variance_explained(y_component),
    );

    // `show_cree` splits the figure into a scatter panel plus a scree panel,
    // mirroring R's `layout(matrix(c(1, 2), ncol = 2))`.
    let width = args.width * 72.0;
    let height = args.height * 72.0;
    let (scatter_width, scree_width) = if args.show_cree {
        (width * 0.72, width * 0.28)
    } else {
        (width, 0.0)
    };

    let mut writer = plot::svg::SvgWriter::new(width, height);
    writer.panel(0.0, 0.0, scatter_width, height, |draw| {
        plot::pca::panel::draw_scatter_panel(
            draw,
            &points,
            &x_title,
            &y_title,
            &legend_entries,
            args.show_axis,
            args.lab_size,
            args.point_size,
            DEFAULT_FONT_SIZE,
        );
    });
    if args.show_cree {
        writer.panel(scatter_width, 0.0, scree_width, height, |draw| {
            plot::pca::panel::draw_scree_panel(
                draw,
                &scree_names(&pca),
                &scree_values(&pca),
                args.show_axis,
                DEFAULT_FONT_SIZE,
            );
        });
    }
    let svg_document = writer.finish();

    match plot::render::OutputFormat::from_path(&args.out) {
        plot::render::OutputFormat::Svg => {
            fs::write(&args.out, svg_document)
                .with_context(|| format!("write SVG: {:?}", args.out))?;
        }
        plot::render::OutputFormat::Pdf => {
            let pdf = plot::render::svg_to_pdf(&svg_document)?;
            fs::write(&args.out, pdf).with_context(|| format!("write PDF: {:?}", args.out))?;
        }
    }

    // Record the component table so a caller can read the variance shares and
    // reuse the scores without re-fitting.
    if let Some(work_dir) = args.work_dir.as_ref() {
        fs::create_dir_all(work_dir)
            .with_context(|| format!("create work dir: {work_dir:?}"))?;

        let mut out = BufWriter::new(File::create(work_dir.join("pca_components.tsv"))?);
        writeln!(&mut out, "component\tsdev\tvariance_explained")?;
        for index in 0..pca.components.len() {
            writeln!(
                &mut out,
                "PC{}\t{}\t{}",
                index + 1,
                pca.sdev(index),
                pca.variance_explained(index)
            )?;
        }
        out.flush()?;

        let mut out = BufWriter::new(File::create(work_dir.join("pca_scores.tsv"))?);
        write!(&mut out, "sample")?;
        for index in 0..pca.components.len() {
            write!(&mut out, "\tPC{}", index + 1)?;
        }
        writeln!(&mut out)?;
        for sample in 0..pca.n_samples {
            write!(
                &mut out,
                "{}",
                sanitize_tsv_value(&table.sample_names[sample])
            )?;
            for component in &pca.components {
                write!(
                    &mut out,
                    "\t{}",
                    component.scores.get(sample).copied().unwrap_or(f64::NAN)
                )?;
            }
            writeln!(&mut out)?;
        }
        out.flush()?;

        // Which regions actually entered the fit; `--top` can drop most of them.
        let mut out = BufWriter::new(File::create(work_dir.join("pca_regions.tsv"))?);
        writeln!(&mut out, "rank\tregion_index\tstandard_deviation\tused")?;
        let used = table.n_regions().min(args.top);
        for (rank, region_index) in ranking.iter().enumerate() {
            writeln!(
                &mut out,
                "{}\t{}\t{}\t{}",
                rank + 1,
                region_index,
                plot::pca::region_sd(&table.regions[*region_index]),
                if rank < used { 1 } else { 0 }
            )?;
        }
        out.flush()?;
    }

    Ok(())
}

/// Component names for the scree panel, e.g. `PC1`, `PC2`, ...
fn scree_names(pca: &plot::pca::Pca) -> Vec<String> {
    (1..=pca.components.len())
        .map(|index| format!("PC{index}"))
        .collect()
}

/// Variance shares for the scree panel, matching `pca_var_explained`.
fn scree_values(pca: &plot::pca::Pca) -> Vec<f64> {
    pca.components
        .iter()
        .map(|component| component.variance_explained)
        .collect()
}

/// Splits one table row on tabs or commas and trims the fields.
///
/// The differential-binding tables users have on hand come from several tools,
/// and `topTable()` writes a TSV while `write.csv()`-style exports are commas.
/// Accepting both avoids forcing a format conversion first.
fn split_table_row(line: &str) -> Vec<String> {
    if line.contains('\t') {
        line.split('\t').map(|field| field.trim().to_string()).collect()
    } else {
        line.split(',').map(|field| field.trim().to_string()).collect()
    }
}

/// Reads a differential-binding results table into [`plot::volcano::Peak`]s.
///
/// Returns the peaks and the contrast label found in the table's comment header,
/// if any. `topTable()` writes its `contrast` attribute as a `# contrast: ...`
/// comment when the table is saved with `write.table()`, which is where the
/// volcano title comes from.
///
/// Column lookup is by name, with the aliases the common tools emit. A header row
/// is required: guessing which unnamed column is the fold change would be worse
/// than saying so.
fn read_results_table(
    path: &Path,
    log_fc_col: Option<&str>,
    p_col: Option<&str>,
    padj_col: Option<&str>,
) -> Result<(Vec<plot::volcano::Peak>, Option<String>)> {
    let text = fs::read_to_string(path).with_context(|| format!("read results: {path:?}"))?;

    // `limma`'s attributes survive as comments; look for a contrast label.
    let mut contrast: Option<String> = None;
    let mut lines = text.lines().peekable();
    while let Some(line) = lines.peek() {
        let trimmed = line.trim();
        if !trimmed.starts_with('#') {
            break;
        }
        if let Some(rest) = trimmed.trim_start_matches('#').trim().strip_prefix("contrast") {
            let value = rest.trim_start_matches([':', '=', ' ']).trim();
            if !value.is_empty() {
                contrast = Some(value.to_string());
            }
        }
        lines.next();
    }

    let header_line = lines
        .find(|line| !line.trim().is_empty())
        .ok_or_else(|| anyhow!("results table is empty: {path:?}"))?;
    let header = split_table_row(header_line);

    let log_fc_aliases = ["logfc", "log2foldchange", "log2fc"];
    let p_aliases = ["p.value", "pvalue", "p_value", "pval"];
    let padj_aliases = ["adj.p.val", "adj.pval", "padj", "p.adjust", "fdr"];

    let find_column = |explicit: Option<&str>,
                       aliases: &[&str],
                       what: &str,
                       flag: &str|
     -> Result<usize> {
        if let Some(name) = explicit {
            return header
                .iter()
                .position(|field| field.eq_ignore_ascii_case(name))
                .ok_or_else(|| {
                    anyhow!(
                        "volcano: {what} column {name:?} not found; table has: {}",
                        header.join(", ")
                    )
                });
        }
        header
            .iter()
            .position(|field| aliases.iter().any(|alias| field.eq_ignore_ascii_case(alias)))
            .ok_or_else(|| {
                anyhow!(
                    "volcano: could not find a {what} column (looked for {}); \
                     table has: {}. Use --{flag} to name it explicitly.",
                    aliases.join(", "),
                    header.join(", ")
                )
            })
    };

    let log_fc_index = find_column(
        log_fc_col,
        &log_fc_aliases,
        "log fold change",
        "logfc-col",
    )?;
    let p_index = find_column(p_col, &p_aliases, "p-value", "p-col")?;
    let padj_index = find_column(
        padj_col,
        &padj_aliases,
        "adjusted p-value",
        "padj-col",
    )?;

    let mut peaks = Vec::new();
    for line in lines {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let fields = split_table_row(trimmed);
        // A row may legitimately carry NA, which R reads as a missing value
        // rather than an error, so parsing failures become NaN instead of
        // aborting the whole table.
        let number_at = |index: usize| -> f64 {
            fields
                .get(index)
                .map(|field| field.trim().parse::<f64>().unwrap_or(f64::NAN))
                .unwrap_or(f64::NAN)
        };
        peaks.push(plot::volcano::Peak {
            log_fold_change: number_at(log_fc_index),
            p_value: number_at(p_index),
            adjusted_p_value: number_at(padj_index),
        });
    }

    if peaks.is_empty() {
        return Err(anyhow!("results table has no data rows: {path:?}"));
    }
    Ok((peaks, contrast))
}

fn cmd_volcano(args: VolcanoArgs) -> Result<()> {
    if !(args.fdr > 0.0) || !args.fdr.is_finite() {
        return Err(anyhow!("volcano: --fdr must be a positive number"));
    }

    let (peaks, contrast) =
        read_results_table(
            &args.results,
            args.log_fc_col.as_deref(),
            args.p_col.as_deref(),
            args.padj_col.as_deref(),
        )?;
    let title = args
        .title
        .clone()
        .or(contrast)
        .unwrap_or_default();
    let volcano = plot::volcano::Volcano::new(peaks, title);

    if volcano.drawable().count() == 0 {
        return Err(anyhow!(
            "volcano: no peaks have finite logFC and p-value; nothing to plot"
        ));
    }
    // R's `range(res$logFC)` returns NA when any fold change is missing, so it
    // aborts and draws nothing. Plotting the usable peaks is more useful, but
    // the difference from R must not be silent.
    if volcano.skipped() > 0 {
        eprintln!(
            "volcano: skipping {} of {} peaks that cannot be plotted \
             (a non-finite logFC, or a p-value of 0)",
            volcano.skipped(),
            volcano.peaks.len()
        );
    }
    // R builds `ylim = c(0, max(-log10(P.Value)))`, which is infinite when any
    // p-value is exactly 0, and then `plot()` aborts. Report it the way R's
    // error does rather than emitting a figure with no y scale.
    if !volcano.y_max().is_finite() {
        return Err(anyhow!(
            "volcano: a p-value of 0 makes -log10(p) infinite; R fails here too \
             (\"need finite 'ylim' values\"). Drop or floor those rows."
        ));
    }

    let (down_count, up_count) = volcano.counts(args.fdr);
    let width = args.width * 72.0;
    let height = args.height * 72.0;
    let mut writer = plot::svg::SvgWriter::new(width, height);
    writer.panel(0.0, 0.0, width, height, |draw| {
        plot::volcano::draw_volcano_panel(
            draw,
            &volcano,
            args.fdr,
            &args.upcol,
            &args.downcol,
            args.alpha,
            args.point_size,
            args.show_axis,
            DEFAULT_FONT_SIZE,
        );
    });
    let svg_document = writer.finish();

    match plot::render::OutputFormat::from_path(&args.out) {
        plot::render::OutputFormat::Svg => {
            fs::write(&args.out, svg_document)
                .with_context(|| format!("write SVG: {:?}", args.out))?;
        }
        plot::render::OutputFormat::Pdf => {
            let pdf = plot::render::svg_to_pdf(&svg_document)?;
            fs::write(&args.out, pdf).with_context(|| format!("write PDF: {:?}", args.out))?;
        }
    }

    // Record the classification so a caller can reuse it without re-deriving it.
    if let Some(work_dir) = args.work_dir.as_ref() {
        fs::create_dir_all(work_dir)
            .with_context(|| format!("create work dir: {work_dir:?}"))?;
        let mut out = BufWriter::new(File::create(work_dir.join("volcano_summary.tsv"))?);
        writeln!(
            &mut out,
            "total\tdrawable\tskipped\tsignificant_down\tsignificant_up\tfdr\tx_min\tx_max\ty_max"
        )?;
        let (x_min, x_max) = volcano.x_range().unwrap_or((f64::NAN, f64::NAN));
        writeln!(
            &mut out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            volcano.peaks.len(),
            volcano.drawable().count(),
            volcano.skipped(),
            down_count,
            up_count,
            args.fdr,
            x_min,
            x_max,
            volcano.y_max()
        )?;
        out.flush()?;
    }

    Ok(())
}

fn cmd_homer_annots(args: HomerAnnotsArgs) -> Result<()> {
    if args.anno.is_empty() {
        return Err(anyhow!("homer-annots: at least one --anno is required"));
    }
    if !(args.legend_font_size > 0.0) || !args.legend_font_size.is_finite() {
        return Err(anyhow!(
            "homer-annots: --legend-font-size must be a positive number"
        ));
    }

    let names = split_csv(&args.samples.join(","));
    let samples: Vec<plot::homer::Sample> = args
        .anno
        .iter()
        .enumerate()
        .map(|(index, path)| {
            let name = names.get(index).cloned().unwrap_or_default();
            plot::homer::parse_sample(path, &name)
        })
        .collect::<Result<Vec<_>>>()?;

    let summary = plot::homer::build_summary(&samples, args.keep_unannotated);    if summary.categories.is_empty() {
        return Err(anyhow!(
            "homer-annots: no peaks fell into a plottable category. R draws only \
             the categories in its fixed palette ({}); this input had none of them.",
            plot::homer::ANNOTATION_PALETTE
                .iter()
                .map(|(name, _)| *name)
                .collect::<Vec<_>>()
                .join(", ")
        ));
    }

    // R drops anything outside the palette silently, which is why its bars stop
    // short of 1. Surface that rather than reproducing the surprise, and call
    // out the `NA` case specifically because it follows from `fread()`'s NA
    // coercion rather than from a decision about that category.
    for dropped in &summary.dropped {
        if dropped.is_unmapped_na() {
            eprintln!(
                "homer-annots: the `NA` category ({} unannotated peak{}) is not drawn. \
                 R names a palette colour for it, but `fread()` reads the literal NA text \
                 as a missing value and `%in%` never matches NA, so the entry is \
                 unreachable. Pass --keep-unannotated to draw it.",
                dropped.total(),
                if dropped.total() == 1 { "" } else { "s" }
            );
        } else {
            eprintln!(
                "homer-annots: category {:?} is not in the palette, so it is not drawn \
                 ({} peak{} across all samples). Fractions still divide by the full peak \
                 count, so the bars sum to less than 1.",
                dropped.name,
                dropped.total(),
                if dropped.total() == 1 { "" } else { "s" }
            );
        }
    }
    for (index, name) in summary.samples.iter().enumerate() {
        let drawn = summary.column_sum(index);
        if drawn < 1.0 - 1e-9 {
            eprintln!(
                "homer-annots: {name}: drawn segments cover {drawn:.4} of the bar; \
                 the rest is in categories outside the palette."
            );
        }
    }

    let width = args.width * 72.0;
    let height = args.height * 72.0;
    let mut writer = plot::svg::SvgWriter::new(width, height);
    writer.panel(0.0, 0.0, width, height, |draw| {
        plot::homer::draw_homer_panel(
            draw,
            &summary,
            args.legend_font_size,
            args.show_axis,
            DEFAULT_FONT_SIZE,
        );
    });
    let svg_document = writer.finish();

    match plot::render::OutputFormat::from_path(&args.out) {
        plot::render::OutputFormat::Svg => {
            fs::write(&args.out, svg_document)
                .with_context(|| format!("write SVG: {:?}", args.out))?;
        }
        plot::render::OutputFormat::Pdf => {
            let pdf = plot::render::svg_to_pdf(&svg_document)?;
            fs::write(&args.out, pdf).with_context(|| format!("write PDF: {:?}", args.out))?;
        }
    }

    // Record the matrix so a caller can reuse the numbers without re-parsing.
    if let Some(work_dir) = args.work_dir.as_ref() {
        fs::create_dir_all(work_dir)
            .with_context(|| format!("create work dir: {work_dir:?}"))?;
        let mut out = BufWriter::new(File::create(work_dir.join("homer_annotations.tsv"))?);
        write!(&mut out, "annotation")?;
        for name in &summary.samples {
            write!(&mut out, "\t{}", sanitize_tsv_value(name))?;
        }
        writeln!(&mut out)?;
        for (category_index, category) in summary.categories.iter().enumerate() {
            write!(&mut out, "{}", sanitize_tsv_value(category))?;
            for sample_index in 0..summary.samples.len() {
                write!(
                    &mut out,
                    "\t{}",
                    summary.fractions[category_index][sample_index]
                )?;
            }
            writeln!(&mut out)?;
        }
        // A trailing row makes the dropped share explicit rather than leaving
        // the reader to notice the columns do not sum to 1.
        write!(&mut out, "__dropped__")?;
        for sample_index in 0..summary.samples.len() {
            write!(&mut out, "\t{}", 1.0 - summary.column_sum(sample_index))?;
        }
        writeln!(&mut out)?;
        out.flush()?;

        let mut out = BufWriter::new(File::create(work_dir.join("homer_counts.tsv"))?);
        writeln!(&mut out, "sample\tn_peaks\tunannotated")?;
        for (index, name) in summary.samples.iter().enumerate() {
            writeln!(
                &mut out,
                "{}\t{}\t{}",
                sanitize_tsv_value(name),
                summary.n_peaks[index],
                samples[index].unannotated()
            )?;
        }
        out.flush()?;

        // R's `leg` column (`Annotation [N]`) is what its commented-out pie
        // branch would have labelled, and it is the per-sample count anyone
        // checking the picture by hand wants, so record it.
        let mut out = BufWriter::new(File::create(work_dir.join("homer_legend.tsv"))?);
        write!(&mut out, "annotation")?;
        for name in &summary.samples {
            write!(&mut out, "\t{}", sanitize_tsv_value(name))?;
        }
        writeln!(&mut out)?;
        for (category_index, labels) in summary
            .legend_labels(&samples)
            .into_iter()
            .enumerate()
        {
            write!(
                &mut out,
                "{}",
                sanitize_tsv_value(&summary.categories[category_index])
            )?;
            for label in labels {
                write!(&mut out, "\t{}", sanitize_tsv_value(&label))?;
            }
            writeln!(&mut out)?;
        }
        out.flush()?;
    }

    Ok(())
}

fn cmd_matrix(args: MatrixArgs) -> Result<()> {
    if args.starts && args.ends {
        return Err(anyhow!("matrix: only one of --starts/--ends may be set"));
    }
    let (up, down) = parse_up_down(&args.size)?;
    if up == 0 && down == 0 {
        return Err(anyhow!("matrix: size must be non-zero (got {})", args.size));
    }
    let binsize = args.tiled_averages;
    if binsize == 0 {
        return Err(anyhow!("matrix: --tiled-averages must be > 0"));
    }

    let total = up
        .checked_add(down)
        .ok_or_else(|| anyhow!("matrix: up+down overflow"))?;
    if total % binsize != 0 {
        return Err(anyhow!(
            "matrix: (up+down) must be divisible by binsize ({} % {} != 0)",
            total,
            binsize
        ));
    }
    let nbins = (total / binsize) as usize;

    let bed_file = File::open(&args.bed).with_context(|| format!("open bed: {:?}", args.bed))?;
    let bigwig_file =
        BigWigRead::open_file(&args.bigwig).with_context(|| format!("open bigwig: {:?}", args.bigwig))?;
    let mut bigwig = bigwig_file.cached();
    let chrom_lens = chrom_lengths(&bigwig);

    let out_file = File::create(&args.out).with_context(|| format!("create out: {:?}", args.out))?;
    let mut w = BufWriter::new(out_file);

    let reader = BufReader::new(bed_file);
    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let Some(parsed) = parse_bed(&line) else {
            continue;
        };
        let (chrom, entry) = parsed?;
        let start = entry.start;
        let end = entry.end;
        if end <= start {
            // Still emit a row of zeros to preserve row count.
            write_zero_row(&mut w, nbins)?;
            continue;
        }

        let anchor = if args.starts {
            start
        } else if args.ends {
            end
        } else {
            start + (end - start) / 2
        };

        let region_start = anchor as i64 - up as i64;
        let region_end = anchor as i64 + down as i64;

        let chr_len = chrom_lens.get(chrom).copied().unwrap_or(0);
        if chr_len == 0 || region_end <= 0 || region_start as u64 >= chr_len as u64 {
            write_zero_row(&mut w, nbins)?;
            continue;
        }

        let query_start = region_start.max(0) as u32;
        let query_end = (region_end.min(chr_len as i64)) as u32;

        let mut sums = vec![0.0f64; nbins];
        if query_end > query_start {
            let interval = bigwig
                .get_interval(chrom, query_start, query_end)?
                .collect::<Result<Vec<_>, _>>()
                .map_err(|e| anyhow!("bigwig interval error: {e}"))?;
            for v in interval {
                let seg_start = v.start.max(query_start) as i64;
                let seg_end = v.end.min(query_end) as i64;
                if seg_end <= seg_start {
                    continue;
                }
                let mut pos = seg_start;
                while pos < seg_end {
                    let bin_idx = ((pos - region_start) / binsize as i64) as isize;
                    if bin_idx < 0 || (bin_idx as usize) >= nbins {
                        break;
                    }
                    let bin_idx = bin_idx as usize;
                    let bin_start = region_start + bin_idx as i64 * binsize as i64;
                    let bin_end = bin_start + binsize as i64;
                    let overlap_end = seg_end.min(bin_end);
                    let overlap_len = overlap_end - pos;
                    sums[bin_idx] += (overlap_len as f64) * (v.value as f64);
                    pos = overlap_end;
                }
            }
        }

        // Emit averages as a single numeric row.
        for (i, s) in sums.iter().enumerate() {
            let avg = s / (binsize as f64);
            if i + 1 == nbins {
                writeln!(&mut w, "{avg}")?;
            } else {
                write!(&mut w, "{avg}\t")?;
            }
        }
    }

    w.flush()?;
    Ok(())
}

fn cmd_track_extract(args: TrackExtractArgs) -> Result<()> {
    if args.bigwigs.is_empty() {
        return Err(anyhow!("track-extract: at least one --bigwig is required"));
    }
    let have_loci = args.loci.is_some();
    let have_gene = args.gene.is_some();
    if have_loci == have_gene {
        return Err(anyhow!(
            "track-extract: provide exactly one of --loci or --gene"
        ));
    }
    if !args.samples.is_empty() && args.samples.len() != args.bigwigs.len() {
        return Err(anyhow!(
            "track-extract: --sample count ({}) must match --bigwig count ({})",
            args.samples.len(),
            args.bigwigs.len()
        ));
    }
    if args.binsize == 0 {
        return Err(anyhow!("track-extract: --binsize must be > 0"));
    }

    fs::create_dir_all(&args.out_dir)
        .with_context(|| format!("create out dir: {:?}", args.out_dir))?;

    let (chr, start, end, gene_models_gene) = if let Some(loci) = args.loci.as_deref() {
        let (chr, start, end) = parse_loci(loci)?;
        (chr, start, end, None)
    } else {
        let gene_raw = args.gene.as_deref().unwrap();
        let gm = match args.gtf.as_ref() {
            None => {
                let gene_sym = normalize_gene_to_symbol(gene_raw, &args.build)?;
                query_ucsc_refgene_by_symbol(&args.build, &gene_sym)?
            }
            Some(gtf) => {
                let gene_ensg = normalize_gene_to_ensg(gene_raw, &args.build)?;
                match parse_gtf_for_gene(gtf, &gene_ensg) {
                    Ok(gm) => gm,
                    Err(_) => {
                        // Fallbacks:
                        // 1) Some GTFs carry the gene symbol in `gene_name`/`gene`.
                        // 2) If that still fails, fallback to UCSC (symbol) when reachable.
                        match parse_gtf_for_gene(gtf, gene_raw) {
                            Ok(gm) => gm,
                            Err(_) => {
                                let gene_sym = normalize_gene_to_symbol(gene_raw, &args.build)?;
                                query_ucsc_refgene_by_symbol(&args.build, &gene_sym)?
                            }
                        }
                    }
                }
            }
        };
        let chr = gm.chr.clone();
        (chr, gm.start, gm.end, Some(gm))
    };

    let (region_start, region_end) = apply_padding(start, end, args.padding)?;
    let bins = gen_bins(region_start, region_end, args.binsize);
    if bins.is_empty() {
        return Err(anyhow!("track-extract: no bins produced for region"));
    }

    // Write meta for R to reconstruct loci.
    let meta_path = args.out_dir.join("meta.tsv");
    {
        let mut w = BufWriter::new(File::create(&meta_path)?);
        writeln!(&mut w, "key\tvalue")?;
        writeln!(&mut w, "chr\t{chr}")?;
        writeln!(&mut w, "start\t{region_start}")?;
        writeln!(&mut w, "end\t{region_end}")?;
        writeln!(&mut w, "binsize\t{}", args.binsize)?;
        writeln!(&mut w, "loci\t{}:{}-{}", chr, region_start, region_end)?;
    }

    let gene_models_region = if args.no_gene_models {
        None
    } else if gene_models_gene.is_some() {
        // For `--gene`, keep the gene-specific models.
        gene_models_gene.clone()
    } else {
        // For `--loci`, fetch gene models overlapping the (padded) region.
        match query_ucsc_refgene_by_region(&args.build, &chr, region_start, region_end) {
            Ok(gm) if !gm.transcripts.is_empty() => Some(gm),
            _ => None,
        }
    };

    if !args.no_cytoband {
        let cytoband = query_ucsc_cytoband(&args.build, &chr, &args.ideo_tbl)?;
        let cytoband_path = args.out_dir.join("cytoband.tsv");
        let mut w = BufWriter::new(File::create(&cytoband_path)?);
        writeln!(&mut w, "chr\tstart\tend\tband\tstain\tcolor")?;
        for c in cytoband {
            writeln!(
                &mut w,
                "{}\t{}\t{}\t{}\t{}\t{}",
                c.chr, c.start, c.end, c.band, c.stain, c.color
            )?;
        }
    }

    // Optionally write gene models if we looked them up.
    if let Some(gm) = gene_models_region.as_ref() {
        let gene_models_path = args.out_dir.join("gene_models.tsv");
        let mut w = BufWriter::new(File::create(&gene_models_path)?);
        writeln!(
            &mut w,
            "chr\tstart\tend\tstrand\ttx\tgene\texon_start\texon_end"
        )?;
        for tx in &gm.transcripts {
            for (exon_start, exon_end) in &tx.exons {
                writeln!(
                    &mut w,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    tx.chr,
                    tx.start,
                    tx.end,
                    tx.strand,
                    tx.tx,
                    tx.gene,
                    exon_start,
                    exon_end
                )?;
            }
        }
    }

    let sample_names: Vec<String> = if args.samples.is_empty() {
        args.bigwigs
            .iter()
            .map(|p| {
                p.file_name()
                    .and_then(|s| s.to_str())
                    .unwrap_or("sample")
                    .split('.')
                    .next()
                    .unwrap_or("sample")
                    .to_string()
            })
            .collect()
    } else {
        args.samples.clone()
    };

    let tracks_path = args.out_dir.join("tracks.tsv");
    let out_file = File::create(&tracks_path).with_context(|| format!("create out: {:?}", tracks_path))?;
    let mut w = BufWriter::new(out_file);
    writeln!(&mut w, "sample\tchromosome\tstart\tend\tsize\tmax")?;

    for (idx, bw_path) in args.bigwigs.iter().enumerate() {
        let sample = &sample_names[idx];
        let bw =
            BigWigRead::open_file(bw_path).with_context(|| format!("open bigwig: {:?}", bw_path))?;
        let mut bw = bw.cached();
        let chrom_lens = chrom_lengths(&bw);

        let chr_len = chrom_lens.get(&chr).copied().unwrap_or(0);
        let query_start = region_start.min(chr_len);
        let query_end = bins
            .last()
            .map(|b| b.end)
            .unwrap_or(region_end)
            .min(chr_len);

        let mut maxes: Vec<f32> = vec![f32::NAN; bins.len()];
        if chr_len != 0 && query_end > query_start {
            let interval = bw
                .get_interval(&chr, query_start, query_end)?
                .collect::<Result<Vec<_>, _>>()
                .map_err(|e| anyhow!("bigwig interval error: {e}"))?;

            for v in interval {
                let ov_start = v.start.max(query_start);
                let ov_end = v.end.min(query_end);
                if ov_end <= ov_start {
                    continue;
                }
                let first = ((ov_start - region_start) / args.binsize) as usize;
                let last = ((ov_end - 1 - region_start) / args.binsize) as usize;
                let last = last.min(maxes.len().saturating_sub(1));
                for b in first..=last {
                    let cur = maxes[b];
                    if cur.is_nan() || v.value > cur {
                        maxes[b] = v.value;
                    }
                }
            }
        }

        for (bidx, bin) in bins.iter().enumerate() {
            let max = maxes[bidx];
            let max = if max.is_nan() { 0.0 } else { max };
            writeln!(
                &mut w,
                "{sample}\t{chr}\t{}\t{}\t{}\t{}",
                bin.start,
                bin.end,
                bin.size,
                max
            )?;
        }
    }

    w.flush()?;
    Ok(())
}

fn cmd_plot_track(args: PlotTrackArgs) -> Result<()> {
    let have_loci = args.loci.is_some();
    let have_gene = args.gene.is_some();
    if have_loci == have_gene {
        return Err(anyhow!("plot-track: provide exactly one of --loci or --gene"));
    }

    let work_dir = match args.work_dir.clone() {
        Some(d) => {
            fs::create_dir_all(&d).with_context(|| format!("create work dir: {:?}", d))?;
            d
        }
        None => create_temp_workdir("tracktools_plot_")?,
    };

    let (bigwigs, sample_names) = if let Some(coldata_path) = args.coldata.as_ref() {
        read_coldata_tsv(coldata_path).with_context(|| format!("read coldata: {:?}", coldata_path))?
    } else {
        if args.bigwigs.is_empty() {
            return Err(anyhow!(
                "plot-track: provide either --coldata <tsv> or at least one --bigwig"
            ));
        }
        if !args.samples.is_empty() && args.samples.len() != args.bigwigs.len() {
            return Err(anyhow!(
                "plot-track: --sample count ({}) must match --bigwig count ({})",
                args.samples.len(),
                args.bigwigs.len()
            ));
        }
        let sample_names: Vec<String> = if args.samples.is_empty() {
            args.bigwigs
                .iter()
                .map(|p| {
                    p.file_name()
                        .and_then(|s| s.to_str())
                        .unwrap_or("sample")
                        .split('.')
                        .next()
                        .unwrap_or("sample")
                        .to_string()
                })
                .collect()
        } else {
            args.samples.clone()
        };
        (args.bigwigs.clone(), sample_names)
    };
    let resolved_cols = resolve_track_colors(&args.col, bigwigs.len());

    // 1) Run extraction
    let extract_args = TrackExtractArgs {
        out_dir: work_dir.clone(),
        binsize: args.binsize,
        padding: args.padding,
        loci: args.loci.clone(),
        gene: args.gene.clone(),
        build: args.build.clone(),
        ideo_tbl: args.ideo_tbl.clone(),
        no_cytoband: !args.show_ideogram,
        no_gene_models: !args.draw_gene_track,
        gtf: args.gtf.clone(),
        bigwigs: bigwigs.clone(),
        samples: sample_names.clone(),
    };
    cmd_track_extract(extract_args)?;

    // Read meta to get final plotted region (padded).
    let meta_path = work_dir.join("meta.tsv");
    let (plot_chr, _plot_start, _plot_end) =
        read_meta_region(&meta_path).context("read plotted region from meta.tsv")?;

    // If requested, fetch UCSC chromHMM tables into local files (Rust-native; no mysql binary).
    let mut chromhmm_paths: Vec<PathBuf> = Vec::new();
    for p in &args.chromhmm {
        chromhmm_paths.push(
            fs::canonicalize(p).with_context(|| format!("chromhmm file not found: {:?}", p))?,
        );
    }
    if !args.ucsc_chromhmm.is_empty() {
        for tbl in &args.ucsc_chromhmm {
            let rows = query_ucsc_chromhmm_by_chr(&args.build, &plot_chr, tbl)?;
            let out_path = work_dir.join(format!("chromHMM_{tbl}.tsv"));
            write_ucsc_chromhmm_tsv(&out_path, &rows)?;
            chromhmm_paths.push(out_path);
        }
    }
    let chromhmm_paths: Vec<PathBuf> = chromhmm_paths
        .into_iter()
        .map(|p| fs::canonicalize(&p).unwrap_or(p))
        .collect();

    let peaks_abs: Vec<PathBuf> = args
        .peaks
        .iter()
        .map(|p| fs::canonicalize(p).with_context(|| format!("peaks file not found: {:?}", p)))
        .collect::<Result<Vec<_>>>()?;

    // Validate optional names vectors.
    if let Some(names) = args.peaks_track_names.as_ref() {
        let n = split_csv(names).len();
        if n != peaks_abs.len() {
            return Err(anyhow!(
                "plot-track: --peaks-track-names count ({n}) must match --peaks count ({})",
                peaks_abs.len()
            ));
        }
    }
    if let Some(names) = args.chromhmm_names.as_ref() {
        let n = split_csv(names).len();
        if n != chromhmm_paths.len() {
            return Err(anyhow!(
                "plot-track: --chromhmm-names count ({n}) must match total chromHMM track count ({})",
                chromhmm_paths.len()
            ));
        }
    }

    // 2) Write coldata mapping for R
    let coldata_path = work_dir.join("coldata.tsv");
    {
        let mut w = BufWriter::new(File::create(&coldata_path)?);
        writeln!(&mut w, "bw_files\tbw_sample_names")?;
        for (bw, sn) in bigwigs.iter().zip(sample_names.iter()) {
            writeln!(&mut w, "{}\t{}", bw.display(), sn)?;
        }
    }

    // 3) Write plotting params for R (key/value TSV).
    let plot_params_path = work_dir.join("plot_params.tsv");
    let mut kv: Vec<(String, String)> = Vec::new();
    kv.push(("ref_build".to_string(), args.build.clone()));
    kv.push(("show_ideogram".to_string(), bool_to_r(args.show_ideogram)));
    kv.push(("draw_gene_track".to_string(), bool_to_r(args.draw_gene_track)));
    kv.push(("track_overlay".to_string(), bool_to_r(args.track_overlay)));
    kv.push(("collapse_txs".to_string(), bool_to_r(args.collapse_txs)));
    kv.push(("col".to_string(), resolved_cols.clone()));
    kv.push((
        "group_auto_scale".to_string(),
        bool_to_r(args.group_auto_scale),
    ));
    kv.push(("y_max".to_string(), args.y_max.clone().unwrap_or_default()));
    kv.push(("y_min".to_string(), args.y_min.clone().unwrap_or_default()));
    kv.push(("txname".to_string(), args.txname.clone().unwrap_or_default()));
    kv.push(("genename".to_string(), args.genename.clone().unwrap_or_default()));
    kv.push(("show_axis".to_string(), bool_to_r(args.show_axis)));
    kv.push((
        "track_names".to_string(),
        args.track_names.clone().unwrap_or_default(),
    ));
    kv.push(("track_names_pos".to_string(), args.track_names_pos.to_string()));
    kv.push((
        "track_names_to_left".to_string(),
        bool_to_r(args.track_names_to_left),
    ));
    kv.push(("gene_fsize".to_string(), args.gene_fsize.to_string()));
    kv.push(("bw_ord".to_string(), args.bw_ord.clone().unwrap_or_default()));
    kv.push(("layout_ord".to_string(), args.layout_ord.clone()));
    kv.push((
        "regions_bed".to_string(),
        args.regions_bed
            .as_ref()
            .map(|p| p.display().to_string())
            .unwrap_or_default(),
    ));
    kv.push(("bw_track_height".to_string(), args.bw_track_height.to_string()));
    kv.push((
        "peaks_track_height".to_string(),
        args.peaks_track_height.to_string(),
    ));
    kv.push((
        "gene_track_height".to_string(),
        args.gene_track_height.to_string(),
    ));
    kv.push((
        "scale_track_height".to_string(),
        args.scale_track_height.to_string(),
    ));
    kv.push((
        "chromhmm_track_height".to_string(),
        args.chromhmm_track_height.to_string(),
    ));
    kv.push((
        "cytoband_track_height".to_string(),
        args.cytoband_track_height.to_string(),
    ));
    kv.push((
        "left_mar".to_string(),
        args.left_mar.map(|v| v.to_string()).unwrap_or_default(),
    ));
    kv.push(("boxcol".to_string(), args.boxcol.clone()));
    kv.push(("boxcolalpha".to_string(), args.boxcolalpha.to_string()));
    kv.push(("peaks".to_string(), join_paths_csv(&peaks_abs)));
    kv.push((
        "peaks_track_names".to_string(),
        args.peaks_track_names.clone().unwrap_or_default(),
    ));
    kv.push(("chromhmm".to_string(), join_paths_csv(&chromhmm_paths)));
    kv.push((
        "chromhmm_names".to_string(),
        args.chromhmm_names.clone().unwrap_or_default(),
    ));
    kv.push((
        "chromhmm_cols".to_string(),
        args.chromhmm_cols.clone().unwrap_or_default(),
    ));
    write_kv_tsv(&plot_params_path, &kv)?;

    // 4) Render natively (no Rscript): compose SVG then convert to PDF/SVG.
    let render_options = plot::render::RenderOptions {
        width: args.pdf_width * 72.0,
        height: args.pdf_height * 72.0,
        colors: split_csv(&resolved_cols),
        show_axis: args.show_axis,
        show_ideogram: args.show_ideogram,
        draw_gene_track: args.draw_gene_track,
        track_names_to_left: args.track_names_to_left,
        track_names: args.track_names.as_ref().map(|names| split_csv(names)),
        font_size: DEFAULT_FONT_SIZE,
        bigwig_height: args.bw_track_height,
        peaks_height: args.peaks_track_height,
        gene_height: args.gene_track_height,
        scale_height: args.scale_track_height,
        chromhmm_height: args.chromhmm_track_height,
        cytoband_height: args.cytoband_track_height,
        y_max: args.y_max.as_ref().map(|values| {
            split_csv(values)
                .iter()
                .filter_map(|value| value.trim().parse::<f64>().ok())
                .collect()
        }),
        y_min: args.y_min.as_ref().map(|values| {
            split_csv(values)
                .iter()
                .filter_map(|value| value.trim().parse::<f64>().ok())
                .collect()
        }),
        group_auto_scale: args.group_auto_scale,
        track_overlay: args.track_overlay,
        layout_ord: split_csv(&args.layout_ord)
            .iter()
            .filter_map(|key| key.trim().chars().next())
            .filter_map(plot::layout::TrackKind::from_key)
            .collect(),
        left_margin: args.left_mar.unwrap_or(DEFAULT_LEFT_MARGIN),
        peaks: peaks_abs
            .iter()
            .enumerate()
            .map(|(index, path)| {
                let name = args
                    .peaks_track_names
                    .as_ref()
                    .map(|names| split_csv(names))
                    .and_then(|names| names.get(index).cloned())
                    .unwrap_or_else(|| {
                        path.file_name()
                            .and_then(|value| value.to_str())
                            .unwrap_or("peaks")
                            .to_string()
                    });
                let regions = plot::io::read_regions(path).unwrap_or_default();
                (name, regions)
            })
            .collect(),
        chromhmm: {
            let names: Vec<String> = args
                .chromhmm_names
                .as_ref()
                .map(|raw| split_csv(raw))
                .unwrap_or_default();
            let mut tracks = Vec::new();
            for (index, path) in chromhmm_paths.iter().enumerate() {
                // R strips the wgEncodeBroadHmm/HMM affixes from the track name,
                // so default to the file stem and let the renderer clean it.
                let name = names.get(index).cloned().unwrap_or_else(|| {
                    path.file_stem()
                        .and_then(|value| value.to_str())
                        .unwrap_or("chromHMM")
                        .to_string()
                });
                tracks.push(
                    plot::io::read_chromhmm(path, &name, &plot_chr)
                        .with_context(|| format!("read chromHMM track: {path:?}"))?,
                );
            }
            tracks
        },
        chromhmm_cols: args
            .chromhmm_cols
            .as_ref()
            .map(|raw| {
                split_csv(raw)
                    .iter()
                    .filter_map(|pair| pair.split_once('='))
                    .map(|(key, value)| (key.trim().to_string(), value.trim().to_string()))
                    .collect()
            })
            .unwrap_or_default(),
        ..plot::render::RenderOptions::default()
    };

    plot::render::render_work_dir(&work_dir, &render_options, &args.out_pdf)
        .with_context(|| format!("render tracks to {:?}", args.out_pdf))?;

    Ok(())
}

fn read_coldata_tsv(path: &PathBuf) -> Result<(Vec<PathBuf>, Vec<String>)> {
    let f = File::open(path).with_context(|| format!("open coldata: {:?}", path))?;
    let r = BufReader::new(f);
    let mut lines = r.lines();

    let first = match lines.next() {
        None => return Err(anyhow!("coldata is empty: {:?}", path)),
        Some(l) => l?,
    };

    let mut bw_idx = 0usize;
    let mut sn_idx = 1usize;
    let mut have_header = false;

    let cols: Vec<&str> = first.split('\t').collect();
    if cols.iter().any(|c| *c == "bw_files") || cols.iter().any(|c| *c == "bw_sample_names") {
        have_header = true;
        bw_idx = cols
            .iter()
            .position(|c| *c == "bw_files")
            .ok_or_else(|| anyhow!("coldata missing column bw_files"))?;
        sn_idx = cols
            .iter()
            .position(|c| *c == "bw_sample_names")
            .ok_or_else(|| anyhow!("coldata missing column bw_sample_names"))?;
    }

    let mut bigwigs = Vec::new();
    let mut samples = Vec::new();

    if !have_header {
        let cols: Vec<&str> = first.split('\t').collect();
        if cols.len() < 2 {
            return Err(anyhow!(
                "coldata must have at least 2 tab-separated columns (bw_files, bw_sample_names)"
            ));
        }
        bigwigs.push(PathBuf::from(cols[0]));
        samples.push(cols[1].to_string());
    }

    for line in lines {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        let max_idx = bw_idx.max(sn_idx);
        if cols.len() <= max_idx {
            continue;
        }
        bigwigs.push(PathBuf::from(cols[bw_idx]));
        samples.push(cols[sn_idx].to_string());
    }

    if bigwigs.is_empty() {
        return Err(anyhow!("coldata contains no rows: {:?}", path));
    }
    if samples.iter().any(|s| s.trim().is_empty()) {
        return Err(anyhow!("coldata has empty bw_sample_names"));
    }
    Ok((bigwigs, samples))
}

fn sanitize_tsv_value(s: &str) -> String {
    s.replace('\t', " ").replace('\n', " ").replace('\r', " ")
}

fn write_kv_tsv(path: &PathBuf, kv: &[(String, String)]) -> Result<()> {
    let mut w = BufWriter::new(File::create(path).with_context(|| format!("create: {:?}", path))?);
    writeln!(&mut w, "key\tvalue")?;
    for (k, v) in kv {
        writeln!(
            &mut w,
            "{}\t{}",
            sanitize_tsv_value(k),
            sanitize_tsv_value(v)
        )?;
    }
    w.flush()?;
    Ok(())
}

fn bool_to_r(b: bool) -> String {
    if b { "TRUE" } else { "FALSE" }.to_string()
}

fn split_csv(s: &str) -> Vec<String> {
    s.split(',')
        .map(|x| x.trim())
        .filter(|x| !x.is_empty())
        .map(|x| x.to_string())
        .collect()
}

fn resolve_track_colors(col_arg: &str, n_tracks: usize) -> String {
    let s = col_arg.trim();
    if !s.eq_ignore_ascii_case("auto") {
        return s.to_string();
    }
    if n_tracks == 0 {
        return String::new();
    }
    auto_palette(n_tracks).join(",")
}

fn auto_palette(n: usize) -> Vec<String> {
    // Tableau 10 (distinct, color-blind friendly-ish).
    const TABLEAU10: [&str; 10] = [
        "#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7",
        "#9C755F", "#BAB0AC",
    ];
    if n <= TABLEAU10.len() {
        return TABLEAU10[..n].iter().map(|s| s.to_string()).collect();
    }
    let mut out: Vec<String> = TABLEAU10.iter().map(|s| s.to_string()).collect();
    for i in TABLEAU10.len()..n {
        // Evenly spaced hues for the remainder.
        let h = (i as f64) * (360.0 / (n as f64));
        let (r, g, b) = hsl_to_rgb(h, 0.65, 0.50);
        out.push(format!("#{:02X}{:02X}{:02X}", r, g, b));
    }
    out
}

fn hsl_to_rgb(h_deg: f64, s: f64, l: f64) -> (u8, u8, u8) {
    let h = (h_deg % 360.0) / 360.0;
    let s = s.clamp(0.0, 1.0);
    let l = l.clamp(0.0, 1.0);
    if s == 0.0 {
        let v = (l * 255.0).round() as u8;
        return (v, v, v);
    }
    let q = if l < 0.5 { l * (1.0 + s) } else { l + s - l * s };
    let p = 2.0 * l - q;
    let r = hue_to_rgb(p, q, h + 1.0 / 3.0);
    let g = hue_to_rgb(p, q, h);
    let b = hue_to_rgb(p, q, h - 1.0 / 3.0);
    (
        (r * 255.0).round() as u8,
        (g * 255.0).round() as u8,
        (b * 255.0).round() as u8,
    )
}

fn hue_to_rgb(p: f64, q: f64, mut t: f64) -> f64 {
    if t < 0.0 {
        t += 1.0;
    }
    if t > 1.0 {
        t -= 1.0;
    }
    if t < 1.0 / 6.0 {
        return p + (q - p) * 6.0 * t;
    }
    if t < 1.0 / 2.0 {
        return q;
    }
    if t < 2.0 / 3.0 {
        return p + (q - p) * (2.0 / 3.0 - t) * 6.0;
    }
    p
}

fn join_paths_csv(paths: &[PathBuf]) -> String {
    paths
        .iter()
        .map(|p| p.display().to_string())
        .collect::<Vec<_>>()
        .join(",")
}

fn read_meta_region(path: &PathBuf) -> Result<(String, u32, u32)> {
    let f = File::open(path).with_context(|| format!("open meta: {:?}", path))?;
    let r = BufReader::new(f);
    let mut chr: Option<String> = None;
    let mut start: Option<u32> = None;
    let mut end: Option<u32> = None;

    for (idx, line) in r.lines().enumerate() {
        let line = line?;
        if idx == 0 && line.starts_with("key") {
            continue;
        }
        let Some((k, v)) = line.split_once('\t') else {
            continue;
        };
        match k.trim() {
            "chr" => chr = Some(v.trim().to_string()),
            "start" => start = Some(v.trim().parse::<u32>().context("parse meta start")?),
            "end" => end = Some(v.trim().parse::<u32>().context("parse meta end")?),
            _ => {}
        }
    }

    Ok((
        chr.ok_or_else(|| anyhow!("meta.tsv missing chr"))?,
        start.ok_or_else(|| anyhow!("meta.tsv missing start"))?,
        end.ok_or_else(|| anyhow!("meta.tsv missing end"))?,
    ))
}

#[derive(Clone, Debug)]
struct ChromHmmRecord {
    chr: String,
    start: u32,
    end: u32,
    name: String,
}

fn query_ucsc_chromhmm_by_chr(build: &str, chr: &str, table: &str) -> Result<Vec<ChromHmmRecord>> {
    if !is_safe_ucsc_value(build) {
        return Err(anyhow!(
            "build contains unsupported characters for UCSC query: {build}"
        ));
    }
    if !is_safe_ucsc_value(table) {
        return Err(anyhow!(
            "chromHMM table contains unsupported characters for UCSC query: {table}"
        ));
    }
    let chr = if chr.starts_with("chr") {
        chr.to_string()
    } else {
        format!("chr{chr}")
    };
    if !is_safe_ucsc_value(&chr) {
        return Err(anyhow!(
            "chromosome contains unsupported characters for UCSC query: {chr}"
        ));
    }

    let mut conn = connect_ucsc_mysql(build)?;
    let sql = format!(
        "select chrom, chromStart, chromEnd, name from {table} WHERE chrom = :chr"
    );
    let rows: Vec<(String, u32, u32, String)> = conn
        .exec(sql, params! { "chr" => chr.as_str() })
        .context("UCSC chromHMM query")?;

    Ok(rows
        .into_iter()
        .map(|(chr, start, end, name)| ChromHmmRecord {
            chr,
            start,
            end,
            name,
        })
        .collect())
}

fn write_ucsc_chromhmm_tsv(path: &PathBuf, rows: &[ChromHmmRecord]) -> Result<()> {
    let mut w = BufWriter::new(File::create(path).with_context(|| format!("create: {:?}", path))?);
    for r in rows {
        writeln!(&mut w, "{}\t{}\t{}\t{}", r.chr, r.start, r.end, r.name)?;
    }
    w.flush()?;
    Ok(())
}

fn chrom_lengths<R>(bw: &BigWigRead<R>) -> HashMap<String, u32> {
    let mut m = HashMap::new();
    for c in bw.chroms() {
        m.insert(c.name.clone(), c.length);
    }
    m
}

fn parse_up_down(s: &str) -> Result<(u32, u32)> {
    let (a, b) = s
        .split_once(':')
        .ok_or_else(|| anyhow!("size must be formatted as UP:DOWN (got {s})"))?;
    let up: u32 = a.parse().with_context(|| format!("invalid UP in size: {a}"))?;
    let down: u32 = b.parse().with_context(|| format!("invalid DOWN in size: {b}"))?;
    Ok((up, down))
}

fn write_zero_row(w: &mut BufWriter<File>, nbins: usize) -> io::Result<()> {
    for i in 0..nbins {
        if i + 1 == nbins {
            writeln!(w, "0")?;
        } else {
            write!(w, "0\t")?;
        }
    }
    Ok(())
}

fn create_temp_workdir(prefix: &str) -> Result<PathBuf> {
    let base = std::env::temp_dir();
    let pid = std::process::id();
    for i in 0..100 {
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap_or_default()
            .as_nanos();
        let d = base.join(format!("{prefix}{pid}_{nanos}_{i}"));
        if fs::create_dir(&d).is_ok() {
            return Ok(d);
        }
    }
    Err(anyhow!("failed to create a temporary work directory"))
}

#[derive(Clone, Debug)]
struct Bin {
    start: u32,
    end: u32,
    size: u32,
}

fn gen_bins(start: u32, end: u32, binsize: u32) -> Vec<Bin> {
    let mut bins = Vec::new();
    let mut s = start;
    while s <= end {
        let e = s.saturating_add(binsize);
        bins.push(Bin {
            start: s,
            end: e,
            size: e.saturating_sub(s),
        });
        s = s.saturating_add(binsize);
        if bins.len() > 50_000_000 {
            break;
        }
    }
    bins
}

fn apply_padding(start: u32, end: u32, padding: i64) -> Result<(u32, u32)> {
    if start >= end {
        return Err(anyhow!("invalid region: end must be > start"));
    }
    let s = start as i64 - padding;
    let e = end as i64 + padding;
    let s = s.max(0);
    let e = e.max(0);
    let s: u32 = s
        .try_into()
        .map_err(|_| anyhow!("region start out of range"))?;
    let e: u32 = e
        .try_into()
        .map_err(|_| anyhow!("region end out of range"))?;
    if s >= e {
        return Err(anyhow!("invalid padded region: end must be > start"));
    }
    Ok((s, e))
}

fn parse_loci(loci: &str) -> Result<(String, u32, u32)> {
    let (chr, rest) = loci
        .split_once(':')
        .ok_or_else(|| anyhow!("invalid loci (missing ':'): {loci}"))?;
    let (start_s, end_s) = rest
        .split_once('-')
        .ok_or_else(|| anyhow!("invalid loci (missing '-'): {loci}"))?;
    let start_s = start_s.replace(',', "");
    let end_s = end_s.replace(',', "");
    let start: u32 = start_s
        .parse()
        .with_context(|| format!("invalid loci start: {start_s}"))?;
    let end: u32 = end_s
        .parse()
        .with_context(|| format!("invalid loci end: {end_s}"))?;
    if start >= end {
        return Err(anyhow!("invalid loci: end must be > start ({loci})"));
    }
    Ok((chr.to_string(), start, end))
}

#[derive(Clone, Debug)]
struct TranscriptModel {
    chr: String,
    strand: String,
    tx: String,
    gene: String,
    start: u32,
    end: u32,
    exons: Vec<(u32, u32)>,
}

#[derive(Clone, Debug)]
struct GeneModels {
    chr: String,
    start: u32,
    end: u32,
    transcripts: Vec<TranscriptModel>,
}

fn parse_gtf_for_gene(gtf: &PathBuf, gene_query: &str) -> Result<GeneModels> {
    let f = File::open(gtf).with_context(|| format!("open gtf: {:?}", gtf))?;
    let r = BufReader::new(f);

    let mut gene_chr: Option<String> = None;
    let mut gene_start: Option<u32> = None;
    let mut gene_end: Option<u32> = None;

    let mut txs: HashMap<String, TranscriptModel> = HashMap::new();

    for line in r.lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 9 {
            continue;
        }
        let chr = cols[0];
        let feature = cols[2];
        let start: u32 = match cols[3].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let end: u32 = match cols[4].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let strand = cols[6];
        let info = cols[8];

        // `gene_name` and `gene_id` are both parsed out of `info`, so a line can
        // only match when `info` contains the query. Checking that first avoids
        // allocating an attribute map for every line of a full-genome GTF.
        if !info.contains(gene_query) {
            continue;
        }

        let attrs = parse_gtf_attrs(info);
        let gene_name = attrs
            .get("gene_name")
            .or_else(|| attrs.get("gene"))
            .or_else(|| attrs.get("Name"))
            .or_else(|| attrs.get("gene_id"))
            .map(|s| s.as_str())
            .unwrap_or("");
        let gene_id = attrs.get("gene_id").map(|s| s.as_str()).unwrap_or("");

        let is_match = gene_name == gene_query
            || gene_id == gene_query
            || info.contains(gene_query);
        if !is_match {
            continue;
        }

        match gene_chr.as_deref() {
            None => gene_chr = Some(chr.to_string()),
            Some(c) if c != chr => {
                // Prefer the first chromosome, matching trackplot.R behavior.
                continue;
            }
            _ => {}
        }

        gene_start = Some(gene_start.map_or(start, |s| s.min(start)));
        gene_end = Some(gene_end.map_or(end, |e| e.max(end)));

        if feature != "exon" {
            continue;
        }

        let tx = match attrs.get("transcript_id") {
            Some(t) => t.to_string(),
            None => continue,
        };
        let gene_out = if !gene_name.is_empty() {
            gene_name.to_string()
        } else if !gene_id.is_empty() {
            gene_id.to_string()
        } else {
            gene_query.to_string()
        };

        let entry = txs.entry(tx.clone()).or_insert_with(|| TranscriptModel {
            chr: chr.to_string(),
            strand: strand.to_string(),
            tx,
            gene: gene_out,
            start,
            end,
            exons: Vec::new(),
        });
        entry.start = entry.start.min(start);
        entry.end = entry.end.max(end);
        entry.exons.push((start, end));
    }

    let chr = gene_chr.ok_or_else(|| anyhow!("gene not found in gtf: {gene_query}"))?;
    let start = gene_start.ok_or_else(|| anyhow!("gene not found in gtf: {gene_query}"))?;
    let end = gene_end.ok_or_else(|| anyhow!("gene not found in gtf: {gene_query}"))?;

    let mut transcripts: Vec<TranscriptModel> = txs.into_values().collect();
    for tx in &mut transcripts {
        tx.exons.sort_by_key(|(s, _)| *s);
        tx.exons.dedup();
    }
    transcripts.sort_by(|a, b| a.tx.cmp(&b.tx));

    Ok(GeneModels {
        chr,
        start,
        end,
        transcripts,
    })
}

fn parse_gtf_attrs(info: &str) -> HashMap<String, String> {
    // Very small, permissive parser for `key "value";` entries.
    let mut m = HashMap::new();
    for part in info.split(';') {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }
        let (k, v) = match part.split_once(' ') {
            Some((k, v)) => (k.trim(), v.trim()),
            None => continue,
        };
        let v = v.trim_matches('"').trim();
        if !k.is_empty() && !v.is_empty() {
            m.insert(k.to_string(), v.to_string());
        }
    }
    m
}

fn query_ucsc_refgene_by_symbol(build: &str, gene: &str) -> Result<GeneModels> {
    // Match trackplot.R:
    //   select chrom, txStart, txEnd, strand, name, name2, exonStarts, exonEnds
    //   from refGene where name2="<gene>"
    if !is_safe_ucsc_value(gene) {
        return Err(anyhow!(
            "gene contains unsupported characters for UCSC query: {gene}"
        ));
    }
    if !is_safe_ucsc_value(build) {
        return Err(anyhow!(
            "build contains unsupported characters for UCSC query: {build}"
        ));
    }

    let mut conn = connect_ucsc_mysql(build)?;
    let sql = "select chrom, txStart, txEnd, strand, name, name2, exonStarts, exonEnds from refGene WHERE name2 = :gene";
    let rows: Vec<(String, u32, u32, String, String, String, String, String)> = conn
        .exec(sql, mysql::params! { "gene" => gene })
        .context("UCSC refGene query")?;
    let mut gene_models: Vec<TranscriptModel> = Vec::new();

    for (chr, start, end, strand, tx, gene_name, exon_starts, exon_ends) in rows {
        let starts: Vec<u32> = exon_starts
            .split(',')
            .filter(|s| !s.is_empty())
            .map(|s| s.parse::<u32>())
            .collect::<std::result::Result<Vec<_>, _>>()
            .with_context(|| format!("invalid exonStarts: {exon_starts}"))?;
        let ends: Vec<u32> = exon_ends
            .split(',')
            .filter(|s| !s.is_empty())
            .map(|s| s.parse::<u32>())
            .collect::<std::result::Result<Vec<_>, _>>()
            .with_context(|| format!("invalid exonEnds: {exon_ends}"))?;

        let mut exons = Vec::new();
        for (s, e) in starts.into_iter().zip(ends.into_iter()) {
            if e > s {
                exons.push((s, e));
            }
        }
        exons.sort_by_key(|(s, _)| *s);
        exons.dedup();

        gene_models.push(TranscriptModel {
            chr,
            strand,
            tx,
            gene: gene_name,
            start,
            end,
            exons,
        });
    }

    if gene_models.is_empty() {
        return Err(anyhow!("no transcripts found for gene: {gene}"));
    }

    let chr = gene_models[0].chr.clone();
    let mut start = gene_models[0].start;
    let mut end = gene_models[0].end;
    gene_models.retain(|t| t.chr == chr);
    for t in &gene_models {
        start = start.min(t.start);
        end = end.max(t.end);
    }

    Ok(GeneModels {
        chr,
        start,
        end,
        transcripts: gene_models,
    })
}

fn query_ucsc_refgene_by_region(build: &str, chr: &str, start: u32, end: u32) -> Result<GeneModels> {
    if !is_safe_ucsc_value(build) {
        return Err(anyhow!(
            "build contains unsupported characters for UCSC query: {build}"
        ));
    }
    let chr = if chr.starts_with("chr") {
        chr.to_string()
    } else {
        format!("chr{chr}")
    };
    if !is_safe_ucsc_value(&chr) {
        return Err(anyhow!(
            "chromosome contains unsupported characters for UCSC query: {chr}"
        ));
    }
    if start >= end {
        return Err(anyhow!("invalid region for refGene query"));
    }

    let mut conn = connect_ucsc_mysql(build)?;
    let sql = "select chrom, txStart, txEnd, strand, name, name2, exonStarts, exonEnds \
               from refGene \
               WHERE chrom = :chr AND txStart < :end AND txEnd > :start";
    let rows: Vec<(String, u32, u32, String, String, String, String, String)> = conn
        .exec(
            sql,
            params! {
                "chr" => chr.as_str(),
                "start" => start,
                "end" => end,
            },
        )
        .context("UCSC refGene region query")?;

    let mut gene_models: Vec<TranscriptModel> = Vec::new();
    for (chr, tx_start, tx_end, strand, tx, gene_name, exon_starts, exon_ends) in rows {
        let starts: Vec<u32> = exon_starts
            .split(',')
            .filter(|s| !s.is_empty())
            .map(|s| s.parse::<u32>())
            .collect::<std::result::Result<Vec<_>, _>>()
            .with_context(|| format!("invalid exonStarts: {exon_starts}"))?;
        let ends: Vec<u32> = exon_ends
            .split(',')
            .filter(|s| !s.is_empty())
            .map(|s| s.parse::<u32>())
            .collect::<std::result::Result<Vec<_>, _>>()
            .with_context(|| format!("invalid exonEnds: {exon_ends}"))?;

        let mut exons = Vec::new();
        for (s, e) in starts.into_iter().zip(ends.into_iter()) {
            if e > s {
                exons.push((s, e));
            }
        }
        exons.sort_by_key(|(s, _)| *s);
        exons.dedup();

        gene_models.push(TranscriptModel {
            chr,
            strand,
            tx,
            gene: gene_name,
            start: tx_start,
            end: tx_end,
            exons,
        });
    }

    if gene_models.is_empty() {
        return Ok(GeneModels {
            chr,
            start,
            end,
            transcripts: Vec::new(),
        });
    }

    let chr_out = gene_models[0].chr.clone();
    gene_models.retain(|t| t.chr == chr_out);
    Ok(GeneModels {
        chr: chr_out,
        start,
        end,
        transcripts: gene_models,
    })
}

fn is_safe_ucsc_value(s: &str) -> bool {
    !s.is_empty()
        && s.chars()
            .all(|c| c.is_ascii_alphanumeric() || c == '_' || c == '-' || c == '.')
}

#[derive(Clone, Debug)]
struct CytobandRecord {
    chr: String,
    start: u32,
    end: u32,
    band: String,
    stain: String,
    color: String,
}

fn query_ucsc_cytoband(build: &str, chr: &str, tbl: &str) -> Result<Vec<CytobandRecord>> {
    if !is_safe_ucsc_value(build) {
        return Err(anyhow!(
            "build contains unsupported characters for UCSC query: {build}"
        ));
    }
    if !is_safe_ucsc_value(tbl) {
        return Err(anyhow!(
            "ideo table contains unsupported characters for UCSC query: {tbl}"
        ));
    }
    let chr = if chr.starts_with("chr") {
        chr.to_string()
    } else {
        format!("chr{chr}")
    };
    if !is_safe_ucsc_value(&chr) {
        return Err(anyhow!(
            "chromosome contains unsupported characters for UCSC query: {chr}"
        ));
    }

    let mut conn = connect_ucsc_mysql(build)?;
    let sql = format!(
        "select chrom, chromStart, chromEnd, name, gieStain from {tbl} WHERE chrom = :chr"
    );
    let rows: Vec<(String, u32, u32, String, String)> = conn
        .exec(sql, mysql::params! { "chr" => chr.as_str() })
        .context("UCSC cytoband query")?;

    let mut out: Vec<CytobandRecord> = Vec::new();
    for (chr, start, end, band, stain) in rows {
        let color = stain_to_color(&stain);
        out.push(CytobandRecord {
            chr,
            start,
            end,
            band,
            stain,
            color,
        });
    }

    if out.is_empty() {
        return Err(anyhow!("no cytoband records returned for {chr} ({build})"));
    }
    Ok(out)
}

fn connect_ucsc_mysql(build: &str) -> Result<mysql::Conn> {
    let hosts = [
        "genome-euro-mysql.soe.ucsc.edu",
        "genome-mysql.soe.ucsc.edu",
        "genome-mysql.cse.ucsc.edu",
    ];
    let mut last_err: Option<String> = None;

    for host in hosts {
        let builder = mysql::OptsBuilder::new()
            .ip_or_hostname(Some(host))
            .user(Some("genome"))
            .db_name(Some(build))
            .tcp_port(3306)
            .tcp_connect_timeout(Some(std::time::Duration::from_secs(5)));

        match mysql::Conn::new(builder) {
            Ok(conn) => return Ok(conn),
            Err(e) => last_err = Some(format!("{e}")),
        }
    }

    Err(anyhow!(
        "UCSC mysql connect failed for build '{build}' on all hosts: {}\nHint: if UCSC is unreachable, use `--gtf` for `--gene` lookups and/or disable ideogram (`--no-cytoband` or `--show-ideogram false`).",
        last_err.unwrap_or_else(|| "unknown error".to_string())
    ))
}

fn stain_to_color(stain: &str) -> String {
    match stain {
        "gneg" => "#FFFFFF".to_string(),
        "acen" => "#660033".to_string(),
        "gvar" => "#660099".to_string(),
        "stalk" => "#6600CC".to_string(),
        _ => {
            if let Some(num) = stain.strip_prefix("gpos") {
                if let Ok(pct) = num.parse::<f64>() {
                    // Match trackplot.R: i <- round(256 - i*2.56)
                    let v = (256.0 - pct * 2.56).round().clamp(0.0, 255.0) as u8;
                    return format!("#{:02X}{:02X}{:02X}", v, v, v);
                }
            }
            "#FFFFFF".to_string()
        }
    }
}

#[derive(Clone, Debug)]
enum GeneInputKind {
    Entrez(String),
    Ensg(String),
    Symbol(String),
}

fn classify_gene_input(s: &str) -> GeneInputKind {
    let s = s.trim();
    if !s.is_empty() && s.chars().all(|c| c.is_ascii_digit()) {
        return GeneInputKind::Entrez(s.to_string());
    }
    let up = s.to_ascii_uppercase();
    if up.starts_with("ENSG") && up[4..].chars().all(|c| c.is_ascii_digit() || c == '.') {
        // Strip version suffix if present.
        return GeneInputKind::Ensg(up.split('.').next().unwrap_or(&up).to_string());
    }
    GeneInputKind::Symbol(s.to_string())
}

fn is_human_build(build: &str) -> bool {
    build.to_ascii_lowercase().starts_with("hg")
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum GeneLookupMode {
    Auto,
    Local,
    Online,
    None,
}

fn gene_lookup_mode() -> GeneLookupMode {
    match std::env::var("TRACKTOOLS_GENE_LOOKUP")
        .unwrap_or_else(|_| "auto".to_string())
        .to_ascii_lowercase()
        .as_str()
    {
        "local" => GeneLookupMode::Local,
        "online" => GeneLookupMode::Online,
        "none" => GeneLookupMode::None,
        _ => GeneLookupMode::Auto,
    }
}

fn have_cmd(cmd: &str) -> bool {
    ProcCommand::new(cmd)
        .arg("--help")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok()
}

fn orgdb_sqlite_path(build: &str) -> Option<PathBuf> {
    let p = std::env::var("TRACKTOOLS_ORGDB_SQLITE").ok();
    if let Some(p) = p.as_deref() {
        let p = p.trim();
        if !p.is_empty() {
            let pb = PathBuf::from(p);
            if pb.exists() {
                return Some(pb);
            }
        }
    }

    if !is_human_build(build) {
        return None;
    }
    let candidates = [
        "/usr/lib/R/library/org.Hs.eg.db/extdata/org.Hs.eg.sqlite",
        "/usr/local/lib/R/site-library/org.Hs.eg.db/extdata/org.Hs.eg.sqlite",
        "/usr/local/lib/R/library/org.Hs.eg.db/extdata/org.Hs.eg.sqlite",
    ];
    for c in candidates {
        let p = PathBuf::from(c);
        if p.exists() {
            return Some(p);
        }
    }
    None
}

fn sqlite3_query_one(db: &Path, sql: &str) -> Result<Option<String>> {
    let out = ProcCommand::new("sqlite3")
        .args(["-batch", "-noheader"])
        .arg(db)
        .arg(sql)
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .output()
        .with_context(|| format!("run sqlite3 query: {sql}"))?;
    if !out.status.success() {
        return Err(anyhow!(
            "sqlite3 failed (status {}): {}",
            out.status,
            String::from_utf8_lossy(&out.stderr).trim()
        ));
    }
    let s = String::from_utf8_lossy(&out.stdout).trim().to_string();
    if s.is_empty() {
        Ok(None)
    } else {
        Ok(Some(s.lines().next().unwrap_or("").trim().to_string()))
    }
}

#[derive(Clone, Debug)]
struct OrgDbCli {
    sqlite: PathBuf,
}

impl OrgDbCli {
    fn open_for_build(build: &str) -> Option<Self> {
        if !have_cmd("sqlite3") {
            return None;
        }
        let sqlite = orgdb_sqlite_path(build)?;
        Some(Self { sqlite })
    }

    fn internal_id_by_symbol_or_alias(&self, sym: &str) -> Result<Option<String>> {
        let sql = format!(
            "select _id from gene_info where symbol = '{sym}' limit 1;"
        );
        if let Some(id) = sqlite3_query_one(&self.sqlite, &sql)? {
            return Ok(Some(id));
        }
        let sql = format!(
            "select _id from alias where alias_symbol = '{sym}' limit 1;"
        );
        sqlite3_query_one(&self.sqlite, &sql)
    }

    fn internal_id_by_entrez(&self, entrez: &str) -> Result<Option<String>> {
        let sql = format!("select _id from genes where gene_id = '{entrez}' limit 1;");
        sqlite3_query_one(&self.sqlite, &sql)
    }

    fn internal_id_by_ensg(&self, ensg: &str) -> Result<Option<String>> {
        let sql = format!("select _id from ensembl where ensembl_id = '{ensg}' limit 1;");
        sqlite3_query_one(&self.sqlite, &sql)
    }

    fn symbol_by_internal_id(&self, id: &str) -> Result<Option<String>> {
        let sql = format!("select symbol from gene_info where _id = {id} limit 1;");
        sqlite3_query_one(&self.sqlite, &sql)
    }

    fn ensg_by_internal_id(&self, id: &str) -> Result<Option<String>> {
        let sql = format!("select ensembl_id from ensembl where _id = {id} limit 1;");
        sqlite3_query_one(&self.sqlite, &sql)
    }
}

fn mygene_lookup(
    gene_raw: &str,
) -> Result<Option<(Option<String>, Option<String>, Option<String>)>> {
    let gene = gene_raw.trim();
    if gene.is_empty() {
        return Ok(None);
    }
    let client = match reqwest::blocking::Client::builder()
        .user_agent("tracs/0.1")
        .timeout(Duration::from_secs(8))
        .build()
    {
        Ok(c) => c,
        Err(_) => return Ok(None),
    };

    let fields = "symbol,entrezgene,ensembl.gene";
    let hit: JsonValue = if gene.chars().all(|c| c.is_ascii_digit()) {
        let url = format!("https://mygene.info/v3/gene/{gene}");
        match client.get(url).query(&[("fields", fields)]).send() {
            Ok(resp) => match resp.error_for_status() {
                Ok(r) => r.json().unwrap_or(JsonValue::Null),
                Err(_) => return Ok(None),
            },
            Err(_) => return Ok(None),
        }
    } else {
        let url = "https://mygene.info/v3/query";
        let resp = match client
            .get(url)
            .query(&[
                ("q", gene),
                ("species", "human"),
                ("size", "1"),
                ("fields", fields),
            ])
            .send()
        {
            Ok(r) => r,
            Err(_) => return Ok(None),
        };
        let resp = match resp.error_for_status() {
            Ok(r) => r,
            Err(_) => return Ok(None),
        };
        let j: JsonValue = match resp.json() {
            Ok(v) => v,
            Err(_) => return Ok(None),
        };
        j.get("hits")
            .and_then(|h| h.as_array())
            .and_then(|a| a.first())
            .cloned()
            .unwrap_or(JsonValue::Null)
    };

    if hit.is_null() {
        return Ok(None);
    }

    let sym = hit.get("symbol").and_then(|v| v.as_str()).map(|s| s.to_string());
    let entrez = hit.get("entrezgene").and_then(|v| {
        if let Some(n) = v.as_i64() {
            Some(n.to_string())
        } else if let Some(s) = v.as_str() {
            Some(s.to_string())
        } else {
            None
        }
    });
    let ensg = hit.get("ensembl").and_then(|ens| {
        if let Some(obj) = ens.as_object() {
            obj.get("gene").and_then(|g| g.as_str()).map(|s| s.to_string())
        } else if let Some(arr) = ens.as_array() {
            for e in arr {
                if let Some(obj) = e.as_object() {
                    if let Some(g) = obj.get("gene").and_then(|g| g.as_str()) {
                        return Some(g.to_string());
                    }
                }
            }
            None
        } else {
            None
        }
    });

    if sym.is_none() && entrez.is_none() && ensg.is_none() {
        Ok(None)
    } else {
        Ok(Some((sym, entrez, ensg)))
    }
}

fn normalize_gene_to_symbol(gene_raw: &str, build: &str) -> Result<String> {
    normalize_gene_to_symbol_with_mode(gene_raw, build, gene_lookup_mode())
}

fn normalize_gene_to_symbol_with_mode(
    gene_raw: &str,
    build: &str,
    mode: GeneLookupMode,
) -> Result<String> {
    let kind = classify_gene_input(gene_raw);
    if !is_human_build(build) {
        return match kind {
            GeneInputKind::Symbol(sym) => {
                if !is_safe_ucsc_value(&sym) {
                    Err(anyhow!(
                        "gene contains unsupported characters for UCSC query: {sym} (use a simple symbol/name)"
                    ))
                } else {
                    Ok(sym)
                }
            }
            GeneInputKind::Entrez(entrez) => Err(anyhow!(
                "gene looks like an Entrez ID ({entrez}) but automatic Entrez→symbol mapping is only supported for hg* (human) builds; pass a gene symbol/name"
            )),
            GeneInputKind::Ensg(ensg) => Err(anyhow!(
                "gene looks like an Ensembl ID ({ensg}) but automatic Ensembl→symbol mapping is only supported for hg* (human) builds; pass a gene symbol/name"
            )),
        };
    }

    // Without a GTF we query UCSC `refGene.name2`, which expects a gene symbol/name.
    let db = OrgDbCli::open_for_build(build);

    match kind {
        GeneInputKind::Symbol(sym) => {
            if !is_safe_ucsc_value(&sym) {
                return Err(anyhow!(
                    "gene contains unsupported characters for UCSC query: {sym} (use --gtf for local lookup)"
                ));
            }
            if let Some(db) = db.as_ref() {
                if let Some(id) = db.internal_id_by_symbol_or_alias(&sym)? {
                    if let Some(canon) = db.symbol_by_internal_id(&id)? {
                        return Ok(canon);
                    }
                }
            }
            Ok(sym)
        }
        GeneInputKind::Entrez(entrez) => {
            if let Some(db) = db.as_ref() {
                let Some(id) = db.internal_id_by_entrez(&entrez)? else {
                    return Err(anyhow!("unknown Entrez ID: {entrez}"));
                };
                let Some(sym) = db.symbol_by_internal_id(&id)? else {
                    return Err(anyhow!("Entrez ID {entrez} has no symbol mapping"));
                };
                return Ok(sym);
            }
            if mode == GeneLookupMode::None || mode == GeneLookupMode::Local {
                return Err(anyhow!(
                    "gene looks like an Entrez ID ({entrez}) but no local orgdb sqlite was found for build {build}; set TRACKTOOLS_ORGDB_SQLITE or pass a gene symbol"
                ));
            }
            if let Some((sym, _entrez2, _ensg)) = mygene_lookup(&entrez)? {
                if let Some(sym) = sym {
                    return Ok(sym);
                }
            }
            Err(anyhow!("failed to map Entrez ID to gene symbol: {entrez}"))
        }
        GeneInputKind::Ensg(ensg) => {
            if let Some(db) = db.as_ref() {
                let Some(id) = db.internal_id_by_ensg(&ensg)? else {
                    return Err(anyhow!("unknown Ensembl gene id: {ensg}"));
                };
                let Some(sym) = db.symbol_by_internal_id(&id)? else {
                    return Err(anyhow!("Ensembl gene id {ensg} has no symbol mapping"));
                };
                return Ok(sym);
            }
            if mode == GeneLookupMode::None || mode == GeneLookupMode::Local {
                return Err(anyhow!(
                    "gene looks like an Ensembl ID ({ensg}) but no local orgdb sqlite was found for build {build}; set TRACKTOOLS_ORGDB_SQLITE or pass a gene symbol"
                ));
            }
            if let Some((sym, _entrez2, _ensg2)) = mygene_lookup(&ensg)? {
                if let Some(sym) = sym {
                    return Ok(sym);
                }
            }
            Err(anyhow!("failed to map Ensembl gene id to gene symbol: {ensg}"))
        }
    }
}

fn normalize_gene_to_ensg(gene_raw: &str, build: &str) -> Result<String> {
    normalize_gene_to_ensg_with_mode(gene_raw, build, gene_lookup_mode())
}

fn normalize_gene_to_ensg_with_mode(
    gene_raw: &str,
    build: &str,
    mode: GeneLookupMode,
) -> Result<String> {
    let kind = classify_gene_input(gene_raw);
    if !is_human_build(build) {
        return match kind {
            GeneInputKind::Ensg(ensg) => Ok(ensg),
            _ => Err(anyhow!(
                "automatic symbol/Entrez→Ensembl conversion is only supported for hg* (human) builds; pass an Ensembl gene id (ENSG...) or a gene identifier that matches your GTF"
            )),
        };
    }

    // When a GTF is supplied, this repo's GTFs use Ensembl gene ids as gene_name.
    let db = OrgDbCli::open_for_build(build);

    match kind {
        GeneInputKind::Ensg(ensg) => Ok(ensg),
        GeneInputKind::Entrez(entrez) => {
            if let Some(db) = db.as_ref() {
                let Some(id) = db.internal_id_by_entrez(&entrez)? else {
                    return Err(anyhow!("unknown Entrez ID: {entrez}"));
                };
                let Some(ensg) = db.ensg_by_internal_id(&id)? else {
                    return Err(anyhow!("Entrez ID {entrez} has no Ensembl mapping"));
                };
                return Ok(ensg);
            }
            if mode == GeneLookupMode::None || mode == GeneLookupMode::Local {
                return Err(anyhow!(
                    "gene looks like an Entrez ID ({entrez}) but no local orgdb sqlite was found for build {build}; set TRACKTOOLS_ORGDB_SQLITE or pass an Ensembl gene id"
                ));
            }
            if let Some((_sym, _entrez2, ensg)) = mygene_lookup(&entrez)? {
                if let Some(ensg) = ensg {
                    return Ok(ensg);
                }
            }
            Err(anyhow!("failed to map Entrez ID to Ensembl gene id: {entrez}"))
        }
        GeneInputKind::Symbol(sym) => {
            if !is_safe_ucsc_value(&sym) {
                return Err(anyhow!("gene contains unsupported characters: {sym}"));
            }
            if let Some(db) = db.as_ref() {
                let Some(id) = db.internal_id_by_symbol_or_alias(&sym)? else {
                    return Err(anyhow!("unknown gene symbol: {sym}"));
                };
                let Some(ensg) = db.ensg_by_internal_id(&id)? else {
                    return Err(anyhow!("gene symbol {sym} has no Ensembl mapping"));
                };
                return Ok(ensg);
            }
            if mode == GeneLookupMode::None || mode == GeneLookupMode::Local {
                return Err(anyhow!(
                    "gene looks like a symbol ({sym}) but no local orgdb sqlite was found for build {build}; set TRACKTOOLS_ORGDB_SQLITE or pass an Ensembl gene id"
                ));
            }
            if let Some((_sym2, _entrez2, ensg)) = mygene_lookup(&sym)? {
                if let Some(ensg) = ensg {
                    return Ok(ensg);
                }
            }
            Err(anyhow!("failed to map gene symbol to Ensembl gene id: {sym}"))
        }
    }
}

#[cfg(test)]
mod gene_norm_tests {
    use super::*;

    #[test]
    fn classify_gene_inputs() {
        match classify_gene_input("912") {
            GeneInputKind::Entrez(x) => assert_eq!(x, "912"),
            _ => panic!("expected entrez"),
        }
        match classify_gene_input("ENSG00000158473.5") {
            GeneInputKind::Ensg(x) => assert_eq!(x, "ENSG00000158473"),
            _ => panic!("expected ensg"),
        }
        match classify_gene_input("CD1D") {
            GeneInputKind::Symbol(x) => assert_eq!(x, "CD1D"),
            _ => panic!("expected symbol"),
        }
    }

    #[test]
    fn normalize_cd1d_hg19_with_local_orgdb() -> Result<()> {
        if OrgDbCli::open_for_build("hg19").is_none() {
            eprintln!("orgdb sqlite/sqlite3 not available; skipping");
            return Ok(());
        }
        assert_eq!(
            normalize_gene_to_symbol_with_mode("912", "hg19", GeneLookupMode::Local)?,
            "CD1D"
        );
        assert_eq!(
            normalize_gene_to_symbol_with_mode(
                "ENSG00000158473",
                "hg19",
                GeneLookupMode::Local
            )?,
            "CD1D"
        );
        assert_eq!(
            normalize_gene_to_ensg_with_mode("CD1D", "hg19", GeneLookupMode::Local)?,
            "ENSG00000158473"
        );
        Ok(())
    }
}
