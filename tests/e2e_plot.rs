//! End-to-end tests for `tracs plot`: real bigWig data -> `track-extract` ->
//! native rendering -> PDF.
//!
//! Rendering is done entirely in Rust, so these tests need neither R nor a
//! display server. They are skipped (with an explanatory message) when the
//! downloaded fixtures are unavailable, so `cargo test` stays usable on a
//! machine without the GEO data:
//!
//! - bigWigs from `GSE199964` under `localdata/data/GSE199964_RAW/`
//!   (override with `TRACKTOOLS_TESTDATA_DIR`),
//! - a full hg19 `ensGene` GTF under `localdata/data/hg19.ensGene.gtf`
//!   (override with `TRACKTOOLS_TEST_GTF`).

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

use anyhow::{anyhow, Context, Result};

struct TempDir {
    path: PathBuf,
}

impl TempDir {
    fn new(prefix: &str) -> Result<Self> {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap_or_default()
            .as_nanos();
        let path = std::env::temp_dir().join(format!("{prefix}_{}_{}", std::process::id(), nanos));
        fs::create_dir_all(&path).with_context(|| format!("create temp dir: {path:?}"))?;
        Ok(Self { path })
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        // CI sets this so the rendered PDFs and intermediate TSVs survive for
        // the artifact upload step; otherwise the directory is cleaned up.
        if std::env::var("TRACKTOOLS_TEST_KEEP_WORKDIR").is_ok() {
            eprintln!("keeping test work dir: {:?}", self.path);
            return;
        }
        let _ = fs::remove_dir_all(&self.path);
    }
}

fn tracs_exe() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_tracs"))
}

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn test_data_dir() -> PathBuf {
    if let Ok(from_env) = std::env::var("TRACKTOOLS_TESTDATA_DIR") {
        let trimmed = from_env.trim();
        if !trimmed.is_empty() {
            return PathBuf::from(trimmed);
        }
    }
    repo_root()
        .join("localdata")
        .join("data")
        .join("GSE199964_RAW")
}

/// Full hg19 `ensGene` GTF used by the gene-mode end-to-end test.
///
/// A full annotation (rather than the trimmed `testdata/` fixture) is what real
/// usage looks like: `tracs` has to stream a multi-hundred-megabyte GTF and
/// still resolve a single gene correctly.
fn full_gtf() -> PathBuf {
    if let Ok(from_env) = std::env::var("TRACKTOOLS_TEST_GTF") {
        let trimmed = from_env.trim();
        if !trimmed.is_empty() {
            return PathBuf::from(trimmed);
        }
    }
    repo_root()
        .join("localdata")
        .join("data")
        .join("hg19.ensGene.gtf")
}

fn sorted_bigwigs() -> Vec<PathBuf> {
    let dir = test_data_dir();
    let mut found: Vec<PathBuf> = Vec::new();
    let Ok(entries) = fs::read_dir(&dir) else {
        return found;
    };
    for entry in entries.flatten() {
        let path = entry.path();
        let extension = path.extension().and_then(|s| s.to_str()).unwrap_or("");
        if extension.eq_ignore_ascii_case("bigwig") || extension.eq_ignore_ascii_case("bw") {
            found.push(path);
        }
    }
    found.sort();
    found
}



fn read_tsv(path: &Path) -> Result<Vec<Vec<String>>> {
    let text = fs::read_to_string(path).with_context(|| format!("read: {path:?}"))?;
    Ok(text
        .lines()
        .filter(|line| !line.trim().is_empty())
        .map(|line| line.split('\t').map(|f| f.to_string()).collect())
        .collect())
}

fn read_kv_tsv(path: &Path) -> Result<Vec<(String, String)>> {
    let rows = read_tsv(path)?;
    let mut out = Vec::new();
    for row in rows.iter().skip(1) {
        if row.len() >= 2 {
            out.push((row[0].clone(), row[1].clone()));
        }
    }
    Ok(out)
}

fn lookup_kv(rows: &[(String, String)], key: &str) -> Option<String> {
    rows.iter()
        .find(|(k, _)| k == key)
        .map(|(_, v)| v.clone())
}

/// A PDF always starts with the `%PDF-` magic bytes; anything else means the
/// renderer produced a truncated or empty file.
fn assert_valid_pdf(path: &Path) -> Result<()> {
    let bytes = fs::read(path).with_context(|| format!("read pdf: {path:?}"))?;
    if bytes.len() < 5 || &bytes[..5] != b"%PDF-" {
        return Err(anyhow!(
            "not a valid PDF (len={}, head={:?}): {path:?}",
            bytes.len(),
            &bytes[..bytes.len().min(5)]
        ));
    }
    Ok(())
}

/// Shared precondition check: the GEO bigWigs must be present.
///
/// Returns the bigWigs when the test can run, or `None` when it should skip.
fn plot_prerequisites() -> Option<Vec<PathBuf>> {
    let bigwigs = sorted_bigwigs();
    if bigwigs.is_empty() {
        eprintln!(
            "no bigWig found under {:?}; skipping plot end-to-end test.",
            test_data_dir()
        );
        return None;
    }
    Some(bigwigs)
}

/// Plot a fixed hg19 locus for several samples, mirroring the documented
/// `--coldata` workflow, and verify the rendered PDF plus extracted tracks.
#[test]
fn plot_loci_multi_sample_renders_pdf() -> Result<()> {
    let Some(bigwigs) = plot_prerequisites() else {
        return Ok(());
    };
    let selected: Vec<&PathBuf> = bigwigs.iter().take(2).collect();

    let temp = TempDir::new("tracs_e2e_loci")?;
    let coldata_path = temp.path.join("coldata.tsv");
    let work_dir = temp.path.join("work");
    let pdf_path = temp.path.join("tracks.pdf");

    let mut coldata = String::from("bw_files\tbw_sample_names\n");
    let mut sample_names: Vec<String> = Vec::new();
    for (index, bigwig) in selected.iter().enumerate() {
        let sample_name = format!("sample_{}", index + 1);
        coldata.push_str(&format!("{}\t{}\n", bigwig.display(), sample_name));
        sample_names.push(sample_name);
    }
    fs::write(&coldata_path, coldata)?;

    let region = "chr1:158145820-158156686";
    let binsize = "200";
    let output = Command::new(tracs_exe())
        .args([
            "plot",
            "--out",
            pdf_path.to_str().unwrap(),
            "--loci",
            region,
            "--binsize",
            binsize,
            "--coldata",
            coldata_path.to_str().unwrap(),
            "--col",
            "auto",
            "--show-ideogram",
            "false",
            "--draw-gene-track",
            "false",
            "--group-auto-scale",
            "true",
            "--show-axis",
            "true",
            "--work-dir",
            work_dir.to_str().unwrap(),
        ])
        .output()
        .context("run tracs plot")?;
    if !output.status.success() {
        return Err(anyhow!(
            "tracs plot failed (status {}):\nstdout:\n{}\nstderr:\n{}",
            output.status,
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        ));
    }

    assert_valid_pdf(&pdf_path)?;

    let tracks_path = work_dir.join("tracks.tsv");
    let tracks = read_tsv(&tracks_path)?;
    if tracks.len() < 2 {
        return Err(anyhow!("tracks.tsv has no data rows: {tracks_path:?}"));
    }
    // header + 55 bins over the 10,866 bp region at 200 bp bins, per sample.
    let expected_bins = 10866usize.div_ceil(200);
    if tracks.len() - 1 != expected_bins * selected.len() {
        return Err(anyhow!(
            "tracks.tsv row count mismatch: got {} expected {}",
            tracks.len() - 1,
            expected_bins * selected.len()
        ));
    }
    for sample_name in &sample_names {
        if !tracks[1..].iter().any(|row| &row[0] == sample_name) {
            return Err(anyhow!(
                "sample {sample_name} missing from tracks.tsv (first column)"
            ));
        }
    }

    let meta = read_kv_tsv(&work_dir.join("meta.tsv"))?;
    if lookup_kv(&meta, "loci").as_deref() != Some(region) {
        return Err(anyhow!(
            "meta.tsv loci mismatch: {:?} != {region}",
            lookup_kv(&meta, "loci")
        ));
    }

    // Guard against a silently all-zero extraction: GSE199964 has real signal
    // at this locus, so at least one bin must be positive.
    let peak = tracks[1..]
        .iter()
        .filter_map(|row| row.get(5))
        .filter_map(|value| value.parse::<f64>().ok())
        .fold(0.0f64, f64::max);
    if peak <= 0.0 {
        return Err(anyhow!(
            "no positive signal in tracks.tsv for {region}; extraction may be broken"
        ));
    }

    Ok(())
}

/// Resolve a gene through a full hg19 `ensGene` GTF and render it.
///
/// Uses Ensembl gene ids so the lookup is fully offline and deterministic: the
/// GTF branch normalises the query to an ENSG id, and an ENSG input needs
/// neither `org.Hs.eg.db` nor an online symbol lookup.
#[test]
fn plot_gene_with_full_gtf_renders_pdf() -> Result<()> {
    let Some(bigwigs) = plot_prerequisites() else {
        return Ok(());
    };

    let gtf = full_gtf();
    if !gtf.exists() {
        eprintln!(
            "full GTF not found at {gtf:?}; set TRACKTOOLS_TEST_GTF to enable the gene end-to-end test. Skipping."
        );
        return Ok(());
    }
    let bigwig = &bigwigs[0];

    // (query, expected chrom, expected start, expected end, expected strand)
    // Coordinates are for hg19 ensGene; CD1D is on the plus strand and SLC19A1
    // on the minus strand, so both orientations are covered.
    let genes = [
        ("ENSG00000158473", "chr1", 158_149_737u32, 158_154_686u32, "+"),
        ("ENSG00000173638", "chr21", 46_913_486u32, 46_964_325u32, "-"),
    ];

    for (query, expected_chr, expected_start, expected_end, expected_strand) in genes {
        let temp = TempDir::new("tracs_e2e_gene")?;
        let work_dir = temp.path.join("work");
        let pdf_path = temp.path.join(format!("{query}.pdf"));

        let output = Command::new(tracs_exe())
            .args([
                "plot",
                "--out",
                pdf_path.to_str().unwrap(),
                "--gene",
                query,
                "--gtf",
                gtf.to_str().unwrap(),
                "--build",
                "hg19",
                "--binsize",
                "200",
                "--bigwig",
                bigwig.to_str().unwrap(),
                "--sample",
                "sample_1",
                "--show-ideogram",
                "false",
                "--group-auto-scale",
                "true",
                "--show-axis",
                "true",
                "--work-dir",
                work_dir.to_str().unwrap(),
            ])
            .output()
            .with_context(|| format!("run tracs plot --gene {query}"))?;
        if !output.status.success() {
            return Err(anyhow!(
                "tracs plot --gene {query} failed (status {}):\nstdout:\n{}\nstderr:\n{}",
                output.status,
                String::from_utf8_lossy(&output.stdout),
                String::from_utf8_lossy(&output.stderr)
            ));
        }

        assert_valid_pdf(&pdf_path)?;

        let meta = read_kv_tsv(&work_dir.join("meta.tsv"))?;
        let actual_chr = lookup_kv(&meta, "chr")
            .ok_or_else(|| anyhow!("meta.tsv missing chr for gene {query}"))?;
        if actual_chr != expected_chr {
            return Err(anyhow!(
                "gene {query}: chromosome mismatch, got {actual_chr} expected {expected_chr}"
            ));
        }
        let actual_start: u32 = lookup_kv(&meta, "start")
            .ok_or_else(|| anyhow!("meta.tsv missing start for gene {query}"))?
            .parse()?;
        let actual_end: u32 = lookup_kv(&meta, "end")
            .ok_or_else(|| anyhow!("meta.tsv missing end for gene {query}"))?
            .parse()?;
        if actual_start != expected_start || actual_end != expected_end {
            return Err(anyhow!(
                "gene {query}: region mismatch, got {actual_chr}:{actual_start}-{actual_end} expected {expected_chr}:{expected_start}-{expected_end}"
            ));
        }

        // The full GTF must also yield exon models for the gene track.
        let gene_models_path = work_dir.join("gene_models.tsv");
        let gene_models = read_tsv(&gene_models_path)?;
        if gene_models.len() < 2 {
            return Err(anyhow!(
                "gene_models.tsv has no exon rows for {query}: {gene_models_path:?}"
            ));
        }
        let header = &gene_models[0];
        let strand_column = header
            .iter()
            .position(|column| column == "strand")
            .ok_or_else(|| anyhow!("gene_models.tsv missing strand column"))?;
        let exon_start_column = header
            .iter()
            .position(|column| column == "exon_start")
            .ok_or_else(|| anyhow!("gene_models.tsv missing exon_start column"))?;
        let exon_end_column = header
            .iter()
            .position(|column| column == "exon_end")
            .ok_or_else(|| anyhow!("gene_models.tsv missing exon_end column"))?;

        for row in &gene_models[1..] {
            if row[strand_column] != expected_strand {
                return Err(anyhow!(
                    "gene {query}: unexpected strand {:?}, expected {expected_strand}",
                    row[strand_column]
                ));
            }
            // Exons must sit inside the reported gene span.
            let exon_start: u32 = row[exon_start_column].parse()?;
            let exon_end: u32 = row[exon_end_column].parse()?;
            if exon_start < actual_start || exon_end > actual_end || exon_end < exon_start {
                return Err(anyhow!(
                    "gene {query}: exon {exon_start}-{exon_end} outside gene span {actual_start}-{actual_end}"
                ));
            }
        }
    }

    Ok(())
}

/// Resolve a gene by symbol against the full GTF, which requires a
/// symbol -> Ensembl mapping (local `org.Hs.eg.db` sqlite or the online
/// mygene.info lookup). Skipped when neither is configured.
#[test]
fn plot_gene_symbol_with_full_gtf_renders_pdf() -> Result<()> {
    let Some(bigwigs) = plot_prerequisites() else {
        return Ok(());
    };

    let gtf = full_gtf();
    if !gtf.exists() {
        eprintln!("full GTF not found at {gtf:?}; skipping symbol lookup test.");
        return Ok(());
    }

    let local_orgdb = std::env::var("TRACKTOOLS_ORGDB_SQLITE")
        .ok()
        .map(|value| value.trim().to_string())
        .filter(|value| !value.is_empty() && PathBuf::from(value).exists());
    let online_lookup = std::env::var("TRACKTOOLS_GENE_LOOKUP")
        .map(|value| value.eq_ignore_ascii_case("online"))
        .unwrap_or(false);
    if local_orgdb.is_none() && !online_lookup {
        eprintln!(
            "no symbol->Ensembl mapping available (set TRACKTOOLS_ORGDB_SQLITE or TRACKTOOLS_GENE_LOOKUP=online); skipping symbol lookup test."
        );
        return Ok(());
    }

    let bigwig = &bigwigs[0];
    let temp = TempDir::new("tracs_e2e_gene_symbol")?;
    let work_dir = temp.path.join("work");
    let pdf_path = temp.path.join("cd1d_symbol.pdf");

    let mut command = Command::new(tracs_exe());
    command.args([
        "plot",
        "--out",
        pdf_path.to_str().unwrap(),
        "--gene",
        "CD1D",
        "--gtf",
        gtf.to_str().unwrap(),
        "--build",
        "hg19",
        "--binsize",
        "200",
        "--bigwig",
        bigwig.to_str().unwrap(),
        "--sample",
        "sample_1",
        "--show-ideogram",
        "false",
        "--group-auto-scale",
        "true",
        "--show-axis",
        "true",
        "--work-dir",
        work_dir.to_str().unwrap(),
    ]);
    if let Some(orgdb) = local_orgdb.as_deref() {
        command.env("TRACKTOOLS_ORGDB_SQLITE", orgdb);
    }

    let output = command.output().context("run tracs plot --gene CD1D")?;
    if !output.status.success() {
        return Err(anyhow!(
            "tracs plot --gene CD1D failed (status {}):\nstdout:\n{}\nstderr:\n{}",
            output.status,
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        ));
    }

    assert_valid_pdf(&pdf_path)?;

    let meta = read_kv_tsv(&work_dir.join("meta.tsv"))?;
    if lookup_kv(&meta, "chr").as_deref() != Some("chr1") {
        return Err(anyhow!(
            "CD1D resolved to unexpected chromosome: {:?}",
            lookup_kv(&meta, "chr")
        ));
    }

    Ok(())
}
