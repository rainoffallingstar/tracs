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
    /// Run `track-extract` then call R to plot tracks into a PDF.
    #[command(alias = "plot")]
    PlotTrack(PlotTrackArgs),
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
    /// Path to `trackplot.R` (this repository file).
    #[arg(long = "trackplot-r", default_value = "trackplot.R")]
    trackplot_r: PathBuf,

    /// Output PDF path.
    #[arg(long = "out")]
    out_pdf: PathBuf,

    /// Path to `Rscript` (defaults to `Rscript` in PATH).
    #[arg(long = "rscript", default_value = "Rscript")]
    rscript: String,

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

    let trackplot_r = fs::canonicalize(&args.trackplot_r)
        .with_context(|| format!("trackplot-r not found: {:?}", args.trackplot_r))?;

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
    kv.push(("col".to_string(), resolved_cols));
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

    // 4) Write a small R script that constructs summary_list and calls track_plot()
    let r_script_path = work_dir.join("plot_track.R");
    fs::write(
        &r_script_path,
        r##"
args <- commandArgs(trailingOnly = TRUE)
trackplot_r <- args[[1]]
work_dir <- args[[2]]
out_pdf <- args[[3]]
pdf_w <- as.numeric(args[[4]])
pdf_h <- as.numeric(args[[5]])

suppressPackageStartupMessages(library(data.table))
source(trackplot_r)

params <- fread(file.path(work_dir, "plot_params.tsv"))
if (!all(c("key","value") %in% colnames(params)) && ncol(params) >= 2) {
  colnames(params)[1:2] <- c("key","value")
}
get_param <- function(k, default = "") {
  v <- params[key %in% k, value]
  if (length(v) < 1 || is.na(v[1]) || !nzchar(v[1])) return(default)
  v[1]
}
get_bool <- function(k, default = FALSE) {
  v <- tolower(get_param(k, ifelse(default, "true", "false")))
  v %in% c("true", "t", "1", "yes", "y")
}
get_num <- function(k, default = NA_real_) {
  v <- get_param(k, "")
  if (!nzchar(v)) return(default)
  as.numeric(v)
}
get_csv <- function(k) {
  v <- get_param(k, "")
  if (!nzchar(v)) return(NULL)
  strsplit(v, ",", fixed = TRUE)[[1]]
}
get_csv_num <- function(k) {
  v <- get_csv(k)
  if (is.null(v)) return(NULL)
  as.numeric(v)
}
get_named_vec <- function(k) {
  v <- get_csv(k)
  if (is.null(v)) return(NULL)
  kv <- strsplit(v, "=", fixed = TRUE)
  keys <- vapply(kv, function(x) ifelse(length(x) >= 1, x[[1]], ""), character(1))
  vals <- vapply(kv, function(x) ifelse(length(x) >= 2, x[[2]], ""), character(1))
  ok <- nzchar(keys) & nzchar(vals)
  vals <- vals[ok]
  keys <- keys[ok]
  if (length(keys) == 0) return(NULL)
  names(vals) <- keys
  vals
}

ref_build <- get_param("ref_build", "hg19")
show_ideogram <- get_bool("show_ideogram", TRUE)
draw_gene_track <- get_bool("draw_gene_track", TRUE)
track_overlay <- get_bool("track_overlay", FALSE)
collapse_txs <- get_bool("collapse_txs", TRUE)
col_arg <- get_param("col", "gray70")
group_auto_scale <- get_bool("group_auto_scale", FALSE)
y_max_arg <- get_param("y_max", "")
y_min_arg <- get_param("y_min", "")
txname_arg <- get_param("txname", "")
genename_arg <- get_param("genename", "")
show_axis <- get_bool("show_axis", FALSE)
track_names_arg <- get_param("track_names", "")
track_names_pos <- get_num("track_names_pos", 0)
track_names_to_left <- get_bool("track_names_to_left", FALSE)
gene_fsize <- get_num("gene_fsize", 1)
bw_ord_arg <- get_param("bw_ord", "")
layout_ord_arg <- get_param("layout_ord", "p,b,h,g,c")
regions_bed <- get_param("regions_bed", "")
bw_track_height <- get_num("bw_track_height", 3)
peaks_track_height <- get_num("peaks_track_height", 2)
gene_track_height <- get_num("gene_track_height", 2)
scale_track_height <- get_num("scale_track_height", 2)
chromhmm_track_height <- get_num("chromhmm_track_height", 1)
cytoband_track_height <- get_num("cytoband_track_height", 2)
left_mar <- get_num("left_mar", NA_real_)
boxcol <- get_param("boxcol", "#ffc41a")
boxcolalpha <- get_num("boxcolalpha", 0.4)
peaks_arg <- get_param("peaks", "")
peaks_track_names_arg <- get_param("peaks_track_names", "")
chromhmm_arg <- get_param("chromhmm", "")
chromhmm_names_arg <- get_param("chromhmm_names", "")
chromhmm_cols_arg <- get_param("chromhmm_cols", "")

meta <- fread(file.path(work_dir, "meta.tsv"))
if (!all(c("key","value") %in% colnames(meta)) && ncol(meta) >= 2) {
  colnames(meta)[1:2] <- c("key","value")
}
loci <- meta[key %in% "loci", value][1]

tracks <- fread(file.path(work_dir, "tracks.tsv"))
coldata <- fread(file.path(work_dir, "coldata.tsv"))
setattr(coldata, "is_bw", TRUE)
setattr(coldata, "refbuild", ref_build)

track_summary <- lapply(coldata$bw_sample_names, function(nm) {
  tracks[sample %in% nm, .(chromosome, start, end, size, max)]
})
names(track_summary) <- coldata$bw_sample_names

cyto_path <- file.path(work_dir, "cytoband.tsv")
if (file.exists(cyto_path)) {
  cyto <- fread(cyto_path)
  colnames(cyto)[1:6] <- c("chr","start","end","band","stain","color")
  setkey(cyto, chr, start, end)
} else {
  cyto <- NA
}

etbl_path <- file.path(work_dir, "gene_models.tsv")
if (file.exists(etbl_path)) {
  etbl <- .read_gene_models_tsv(etbl_path)
} else {
  etbl <- NULL
}

attr(track_summary, "meta") <- list(etbl = etbl, cyto = cyto, loci = loci)
summary_list <- list(data = track_summary, colData = coldata)

parse_csv_char <- function(x){
  if (is.null(x) || is.na(x) || !nzchar(x)) return(NULL)
  strsplit(x, ",", fixed = TRUE)[[1]]
}

cols <- parse_csv_char(col_arg)
y_max <- get_csv_num("y_max")
y_min <- get_csv_num("y_min")
txname <- parse_csv_char(txname_arg)
genename <- parse_csv_char(genename_arg)
track_names <- parse_csv_char(track_names_arg)
bw_ord <- parse_csv_char(bw_ord_arg)
layout_ord <- parse_csv_char(layout_ord_arg)
peaks <- parse_csv_char(peaks_arg)
peaks_track_names <- parse_csv_char(peaks_track_names_arg)
chromhmm <- parse_csv_char(chromhmm_arg)
chromhmm_names <- parse_csv_char(chromhmm_names_arg)
chromhmm_cols <- get_named_vec("chromhmm_cols")

regions <- NULL
if (!is.null(regions_bed) && !is.na(regions_bed) && nzchar(regions_bed) && file.exists(regions_bed)) {
  regions <- fread(regions_bed, header = FALSE)
  if (ncol(regions) >= 3) {
    regions <- as.data.frame(regions[,1:3])
  } else {
    regions <- NULL
  }
}

plot_args <- list(
  summary_list = summary_list,
  show_ideogram = show_ideogram,
  draw_gene_track = draw_gene_track,
  track_overlay = track_overlay,
  collapse_txs = collapse_txs,
  groupAutoScale = group_auto_scale,
  show_axis = show_axis,
  track_names_pos = track_names_pos,
  track_names_to_left = track_names_to_left,
  gene_fsize = gene_fsize,
  bw_track_height = bw_track_height,
  peaks_track_height = peaks_track_height,
  gene_track_height = gene_track_height,
  scale_track_height = scale_track_height,
  chromHMM_track_height = chromhmm_track_height,
  cytoband_track_height = cytoband_track_height,
  boxcol = boxcol,
  boxcolalpha = boxcolalpha,
  layout_ord = layout_ord
)
if (!is.null(cols)) plot_args$col <- cols
if (!is.null(y_max)) plot_args$y_max <- y_max
if (!is.null(y_min)) plot_args$y_min <- y_min
if (!is.null(txname)) plot_args$txname <- txname
if (!is.null(genename)) plot_args$genename <- genename
if (!is.null(track_names)) plot_args$track_names <- track_names
if (!is.null(bw_ord)) plot_args$bw_ord <- bw_ord
if (!is.null(regions)) plot_args$regions <- regions
if (!is.na(left_mar)) plot_args$left_mar <- left_mar
if (!is.null(peaks)) plot_args$peaks <- peaks
if (!is.null(peaks_track_names)) plot_args$peaks_track_names <- peaks_track_names
if (!is.null(chromhmm)) plot_args$chromHMM <- chromhmm
if (!is.null(chromhmm_names)) plot_args$chromHMM_names <- chromhmm_names
if (!is.null(chromhmm_cols)) plot_args$chromHMM_cols <- chromhmm_cols

pdf(out_pdf, width = pdf_w, height = pdf_h)
tryCatch({
  do.call(track_plot, plot_args)
}, finally = {
  dev.off()
})
	"##,
    )
    .with_context(|| format!("write R script: {:?}", r_script_path))?;

    // 5) Run Rscript
    let mut cmd = ProcCommand::new(&args.rscript);
    cmd.arg(&r_script_path)
        .arg(trackplot_r)
        .arg(&work_dir)
        .arg(&args.out_pdf)
        .arg(args.pdf_width.to_string())
        .arg(args.pdf_height.to_string());

    let out = cmd.output().context("running Rscript")?;
    if !out.status.success() {
        return Err(anyhow!(
            "Rscript failed (status {}):\n{}",
            out.status,
            String::from_utf8_lossy(&out.stderr)
        ));
    }

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
