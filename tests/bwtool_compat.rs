use std::fs;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

use anyhow::{anyhow, Context, Result};
use bigtools::BigWigRead;

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
        let _ = fs::remove_dir_all(&self.path);
    }
}

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn test_data_dir() -> PathBuf {
    if let Ok(p) = std::env::var("TRACKTOOLS_TESTDATA_DIR") {
        let p = p.trim();
        if !p.is_empty() {
            return PathBuf::from(p);
        }
    }
    repo_root()
        .join("localdata")
        .join("data")
        .join("GSE199964_RAW")
}

fn test_bigwig() -> Result<PathBuf> {
    let d = test_data_dir();
    let preferred = d.join("H3K27ac_1.bigWig");
    if preferred.exists() {
        return Ok(preferred);
    }
    let mut cands: Vec<PathBuf> = Vec::new();
    if let Ok(rd) = fs::read_dir(&d) {
        for e in rd.flatten() {
            let p = e.path();
            let ext = p.extension().and_then(|s| s.to_str()).unwrap_or("");
            if ext.eq_ignore_ascii_case("bigwig")
                || ext.eq_ignore_ascii_case("bigWig")
                || ext.eq_ignore_ascii_case("bw")
            {
                cands.push(p);
            }
        }
    }
    cands.sort();
    cands
        .into_iter()
        .next()
        .ok_or_else(|| anyhow!("no bigWig found under test data dir: {:?}", d))
}

fn open_bw(path: &Path) -> Result<bigtools::BigWigRead<bigtools::utils::reopen::ReopenableFile>> {
    BigWigRead::open_file(path).with_context(|| format!("open bigwig: {path:?}"))
}

#[derive(Clone, Debug)]
struct Seg {
    start: u32,
    end: u32,
    value: f32,
}

fn collect_segments(
    bw: &mut bigtools::BigWigRead<impl bigtools::BBIFileRead>,
    chr: &str,
    start: u32,
    end: u32,
) -> Result<Vec<Seg>> {
    if end <= start {
        return Ok(Vec::new());
    }
    let it = bw
        .get_interval(chr, start, end)
        .with_context(|| format!("get_interval({chr},{start},{end})"))?;
    let v = it
        .collect::<std::result::Result<Vec<_>, _>>()
        .map_err(|e| anyhow!("bigwig interval error: {e}"))?;
    Ok(v.into_iter()
        .map(|x| Seg {
            start: x.start,
            end: x.end,
            value: x.value,
        })
        .collect())
}

fn pick_signal_region(_bw_path: &Path) -> Result<(String, u32, u32)> {
    // Keep tests stable: use a loci already used in repository examples.
    Ok(("chr1".to_string(), 158_145_820, 158_156_686))
}

fn approx_eq(a: f64, b: f64, rel: f64) -> bool {
    if a.is_nan() && b.is_nan() {
        return true;
    }
    if a.is_nan() || b.is_nan() {
        return false;
    }
    let scale = a.abs().max(b.abs()).max(1.0);
    (a - b).abs() <= rel * scale
}

fn integrate_sum_and_max(
    bw: &mut bigtools::BigWigRead<impl bigtools::BBIFileRead>,
    chr: &str,
    start: u32,
    end: u32,
    chr_len: u32,
) -> Result<(f64, f64)> {
    let qs = start.min(chr_len);
    let qe = end.min(chr_len);
    if qe <= qs {
        return Ok((0.0, 0.0));
    }
    let segs = collect_segments(bw, chr, qs, qe)?;
    if segs.is_empty() {
        return Ok((0.0, 0.0));
    }
    let mut sum = 0.0f64;
    let mut max = f64::NEG_INFINITY;
    for s in segs {
        let os = s.start.max(qs);
        let oe = s.end.min(qe);
        if oe <= os {
            continue;
        }
        let ol = (oe - os) as f64;
        sum += ol * (s.value as f64);
        max = max.max(s.value as f64);
    }
    if max.is_infinite() {
        max = 0.0;
    }
    Ok((sum, max))
}

fn integrate_avg_in_bin(
    bw: &mut bigtools::BigWigRead<impl bigtools::BBIFileRead>,
    chr: &str,
    bin_start: i64,
    bin_end: i64,
    chr_len: u32,
    binsize: u32,
) -> Result<f64> {
    if bin_end <= bin_start || binsize == 0 {
        return Ok(0.0);
    }
    if bin_end <= 0 || bin_start as u64 >= chr_len as u64 {
        return Ok(0.0);
    }
    let qs = bin_start.max(0) as u32;
    let qe = (bin_end.min(chr_len as i64)) as u32;
    if qe <= qs {
        return Ok(0.0);
    }
    let segs = collect_segments(bw, chr, qs, qe)?;
    let mut sum = 0.0f64;
    for s in segs {
        let os = (s.start as i64).max(bin_start).max(0) as u32;
        let oe = (s.end as i64).min(bin_end).min(chr_len as i64) as u32;
        if oe <= os {
            continue;
        }
        let ol = (oe - os) as f64;
        sum += ol * (s.value as f64);
    }
    Ok(sum / (binsize as f64))
}

fn have_cmd(cmd: &str) -> bool {
    Command::new(cmd)
        .arg("--help")
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .status()
        .is_ok()
}

#[derive(Clone, Debug)]
enum BwtoolInvoker {
    Path,
    CondaEnv(String),
    MicroMambaEnv { env: String, root_prefix: Option<String> },
}

fn bwtool_test_env() -> Option<String> {
    // Preferred generic name, with backwards-compatible fallback.
    std::env::var("TRACKTOOLS_TEST_BWTOOL_ENV")
        .ok()
        .and_then(|s| if s.trim().is_empty() { None } else { Some(s) })
        .or_else(|| {
            std::env::var("TRACKTOOLS_TEST_BWTOOL_CONDA_ENV")
                .ok()
                .and_then(|s| if s.trim().is_empty() { None } else { Some(s) })
        })
}

fn detect_bwtool_invoker() -> Option<BwtoolInvoker> {
    if have_cmd("bwtool") {
        return Some(BwtoolInvoker::Path);
    }
    let env = bwtool_test_env()?;
    // This system typically uses micromamba; prefer it when available.
    if have_cmd("micromamba") {
        let root_prefix = std::env::var("TRACKTOOLS_TEST_MICROMAMBA_ROOT")
            .ok()
            .and_then(|s| if s.trim().is_empty() { None } else { Some(s) })
            .or_else(|| {
                std::env::var("MAMBA_ROOT_PREFIX")
                    .ok()
                    .and_then(|s| if s.trim().is_empty() { None } else { Some(s) })
            });
        return Some(BwtoolInvoker::MicroMambaEnv { env, root_prefix });
    }
    if have_cmd("conda") {
        return Some(BwtoolInvoker::CondaEnv(env));
    }
    None
}

fn run_bwtool(invoker: &BwtoolInvoker, args: &[&str]) -> Result<()> {
    let mut cmd = match invoker {
        BwtoolInvoker::Path => Command::new("bwtool"),
        BwtoolInvoker::CondaEnv(env) => {
            let mut c = Command::new("conda");
            c.args(["run", "-n", env.as_str(), "bwtool"]);
            c
        }
        BwtoolInvoker::MicroMambaEnv { env, root_prefix } => {
            let mut c = Command::new("micromamba");
            let xdg_cache = std::env::var("TRACKTOOLS_TEST_XDG_CACHE_HOME")
                .ok()
                .filter(|s| !s.trim().is_empty())
                .unwrap_or_else(|| "/tmp/tracktools-xdg-cache".to_string());
            c.env("XDG_CACHE_HOME", xdg_cache);
            if let Some(root) = root_prefix.as_deref() {
                c.args(["run", "-r", root, "-n", env.as_str(), "bwtool"]);
            } else {
                c.args(["run", "-n", env.as_str(), "bwtool"]);
            }
            c
        }
    };
    let status = cmd.args(args).status().context("run bwtool")?;
    if !status.success() {
        return Err(anyhow!("bwtool failed with status {status}"));
    }
    Ok(())
}

fn read_tsv(path: &Path) -> Result<Vec<Vec<String>>> {
    let s = fs::read_to_string(path).with_context(|| format!("read: {path:?}"))?;
    Ok(s.lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| l.split('\t').map(|x| x.to_string()).collect())
        .collect())
}

fn header_col_idx(header: &[String], name: &str) -> Option<usize> {
    header.iter().position(|c| {
        let c = c.trim().trim_start_matches('#');
        c.eq_ignore_ascii_case(name)
    })
}

fn tracktools_exe() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_tracktools"))
}

fn mini_gtf() -> PathBuf {
    repo_root()
        .join("testdata")
        .join("hg19.ensGene.CD1D_SLC19A1.mini.gtf")
}

#[test]
fn summary_matches_reference_integrator() -> Result<()> {
    let bw_path = match test_bigwig() {
        Ok(p) => p,
        Err(e) => {
            eprintln!("missing test bigWig; skipping: {e}");
            return Ok(());
        }
    };
    if !bw_path.exists() {
        eprintln!("missing test bigWig; skipping: {:?}", bw_path);
        return Ok(());
    }
    let (chr, region_start, region_end) = pick_signal_region(&bw_path)?;

    let bw = open_bw(&bw_path)?;
    let mut bw = bw.cached();
    let chr_len = bw
        .chroms()
        .iter()
        .find(|c| c.name == chr)
        .map(|c| c.length)
        .ok_or_else(|| anyhow!("missing chrom length for {chr}"))?;

    let td = TempDir::new("tracktools_test_summary")?;
    let bed_path = td.path.join("q.bed");
    let out_path = td.path.join("out.tsv");

    // 0-based half-open.
    let bed = format!(
        "{chr}\t{region_start}\t{}\n{chr}\t{}\t{}\n",
        (region_start + 50).min(region_end),
        region_start.saturating_add(25),
        (region_start + 125).min(region_end)
    );
    fs::write(&bed_path, bed)?;

    let status = Command::new(tracktools_exe())
        .args([
            "summary",
            "-with-sum",
            "-keep-bed",
            "-header",
            bed_path.to_str().unwrap(),
            bw_path.to_str().unwrap(),
            out_path.to_str().unwrap(),
        ])
        .status()
        .context("run tracktools summary")?;
    if !status.success() {
        return Err(anyhow!("tracktools summary failed"));
    }

    let rows = read_tsv(&out_path)?;
    if rows.len() < 2 {
        return Err(anyhow!("expected header + >=1 data row"));
    }
    // header: chromosome start end size sum min max mean
    for (i, r) in rows.iter().enumerate().skip(1) {
        if r.len() < 7 {
            return Err(anyhow!("row {i} too few columns: {:?}", r));
        }
        let bed_rows = read_tsv(&bed_path)?;
        let bed_r = &bed_rows[i - 1];
        let start: u32 = bed_r[1].parse()?;
        let end: u32 = bed_r[2].parse()?;
        let (sum_ref, max_ref) = integrate_sum_and_max(&mut bw, &chr, start, end, chr_len)?;
        let sum_got: f64 = r[4].parse()?;
        let max_got: f64 = r[6].parse()?;
        if !approx_eq(sum_ref, sum_got, 1e-6) {
            return Err(anyhow!(
                "sum mismatch row {i}: ref={sum_ref} got={sum_got}"
            ));
        }
        if !approx_eq(max_ref, max_got, 1e-6) {
            return Err(anyhow!(
                "max mismatch row {i}: ref={max_ref} got={max_got}"
            ));
        }
    }
    Ok(())
}

#[test]
fn matrix_matches_reference_integrator() -> Result<()> {
    let bw_path = match test_bigwig() {
        Ok(p) => p,
        Err(e) => {
            eprintln!("missing test bigWig; skipping: {e}");
            return Ok(());
        }
    };
    if !bw_path.exists() {
        eprintln!("missing test bigWig; skipping: {:?}", bw_path);
        return Ok(());
    }
    let (chr, region_start, region_end) = pick_signal_region(&bw_path)?;

    let bw = open_bw(&bw_path)?;
    let mut bw = bw.cached();
    let chr_len = bw
        .chroms()
        .iter()
        .find(|c| c.name == chr)
        .map(|c| c.length)
        .ok_or_else(|| anyhow!("missing chrom length for {chr}"))?;

    let td = TempDir::new("tracktools_test_matrix")?;
    let bed_path = td.path.join("q.bed");
    let out_path = td.path.join("out.tsv");

    // Two anchors.
    let bed = format!(
        "{chr}\t{region_start}\t{}\n{chr}\t{}\t{}\n",
        (region_start + 10).min(region_end),
        (region_start + 100).min(region_end.saturating_sub(1)),
        (region_start + 120).min(region_end)
    );
    fs::write(&bed_path, bed)?;

    let up = 200u32;
    let down = 200u32;
    let binsize = 50u32;
    let nbins = ((up + down) / binsize) as usize;

    let status = Command::new(tracktools_exe())
        .args([
            "matrix",
            "-starts",
            format!("-tiled-averages={binsize}").as_str(),
            format!("{up}:{down}").as_str(),
            bed_path.to_str().unwrap(),
            bw_path.to_str().unwrap(),
            out_path.to_str().unwrap(),
        ])
        .status()
        .context("run tracktools matrix")?;
    if !status.success() {
        return Err(anyhow!("tracktools matrix failed"));
    }

    let rows = read_tsv(&out_path)?;
    let bed_rows = read_tsv(&bed_path)?;
    if rows.len() != bed_rows.len() {
        return Err(anyhow!(
            "row count mismatch: matrix {} bed {}",
            rows.len(),
            bed_rows.len()
        ));
    }
    for (row_idx, row) in rows.iter().enumerate() {
        if row.len() != nbins {
            return Err(anyhow!(
                "matrix cols mismatch row {row_idx}: got {} expected {nbins}",
                row.len()
            ));
        }
        let start: u32 = bed_rows[row_idx][1].parse()?;
        let anchor = start as i64;
        let region_start = anchor - up as i64;
        for b in 0..nbins {
            let bin_start = region_start + (b as i64) * binsize as i64;
            let bin_end = bin_start + binsize as i64;
            let ref_avg =
                integrate_avg_in_bin(&mut bw, &chr, bin_start, bin_end, chr_len, binsize)?;
            let got: f64 = row[b].parse()?;
            if !approx_eq(ref_avg, got, 1e-6) {
                return Err(anyhow!(
                    "matrix mismatch row {row_idx} bin {b}: ref={ref_avg} got={got}"
                ));
            }
        }
    }
    Ok(())
}

#[test]
fn track_extract_max_matches_reference() -> Result<()> {
    let bw_path = match test_bigwig() {
        Ok(p) => p,
        Err(e) => {
            eprintln!("missing test bigWig; skipping: {e}");
            return Ok(());
        }
    };
    if !bw_path.exists() {
        eprintln!("missing test bigWig; skipping: {:?}", bw_path);
        return Ok(());
    }
    let (chr, region_start, region_end) = pick_signal_region(&bw_path)?;
    let loci = format!("{chr}:{region_start}-{region_end}");

    let td = TempDir::new("tracktools_test_extract")?;
    let out_dir = td.path.join("out");
    fs::create_dir_all(&out_dir)?;

    let status = Command::new(tracktools_exe())
        .args([
            "track-extract",
            "--out-dir",
            out_dir.to_str().unwrap(),
            "--binsize",
            "50",
            "--no-cytoband",
            "--no-gene-models",
            "--loci",
            &loci,
            "--bigwig",
            bw_path.to_str().unwrap(),
            "--sample",
            "S1",
        ])
        .status()
        .context("run tracktools track-extract")?;
    if !status.success() {
        return Err(anyhow!("tracktools track-extract failed"));
    }

    let tracks_path = out_dir.join("tracks.tsv");
    let rows = read_tsv(&tracks_path)?;
    if rows.len() < 2 {
        return Err(anyhow!("tracks.tsv missing data"));
    }
    let bw = open_bw(&bw_path)?;
    let mut bw = bw.cached();
    let chr_len = bw
        .chroms()
        .iter()
        .find(|c| c.name == chr)
        .map(|c| c.length)
        .ok_or_else(|| anyhow!("missing chrom length for {chr}"))?;

    for (i, r) in rows.iter().enumerate().skip(1) {
        if r.len() < 6 {
            return Err(anyhow!("tracks row too few cols at {i}: {:?}", r));
        }
        let start: u32 = r[2].parse()?;
        let end: u32 = r[3].parse()?;
        let segs = collect_segments(&mut bw, &chr, start.min(chr_len), end.min(chr_len))?;
        let mut max_ref = 0.0f64;
        for s in segs {
            if s.end <= s.start {
                continue;
            }
            max_ref = max_ref.max(s.value as f64);
        }
        let got: f64 = r[5].parse()?;
        if !approx_eq(max_ref, got, 1e-6) {
            return Err(anyhow!(
                "track-extract max mismatch row {i}: ref={max_ref} got={got}"
            ));
        }
    }
    Ok(())
}

#[test]
fn bwtool_summary_and_matrix_match_when_available() -> Result<()> {
    let Some(invoker) = detect_bwtool_invoker() else {
        eprintln!(
            "bwtool not found; set TRACKTOOLS_TEST_BWTOOL_ENV=<env> to run bwtool via micromamba/conda. Skipping bwtool compatibility check."
        );
        return Ok(());
    };

    let bw_path = match test_bigwig() {
        Ok(p) => p,
        Err(e) => {
            eprintln!("missing test bigWig; skipping: {e}");
            return Ok(());
        }
    };
    if !bw_path.exists() {
        eprintln!("missing test bigWig; skipping: {:?}", bw_path);
        return Ok(());
    }
    let (chr, region_start, region_end) = pick_signal_region(&bw_path)?;

    let td = TempDir::new("tracktools_test_bwtool")?;
    let bed_path = td.path.join("q.bed");
    let out_tracktools = td.path.join("tracktools.summary.tsv");
    let out_bwtool = td.path.join("bwtool.summary.tsv");
    let mat_tracktools = td.path.join("tracktools.matrix.tsv");
    let mat_bwtool = td.path.join("bwtool.matrix.tsv");

    let bed = format!(
        "{chr}\t{region_start}\t{}\n{chr}\t{}\t{}\n",
        (region_start + 50).min(region_end),
        region_start.saturating_add(25),
        (region_start + 125).min(region_end)
    );
    fs::write(&bed_path, bed)?;

    // summary
    let status = Command::new(tracktools_exe())
        .args([
            "summary",
            "-with-sum",
            "-keep-bed",
            "-header",
            bed_path.to_str().unwrap(),
            bw_path.to_str().unwrap(),
            out_tracktools.to_str().unwrap(),
        ])
        .status()
        .context("run tracktools summary")?;
    if !status.success() {
        return Err(anyhow!("tracktools summary failed"));
    }

    run_bwtool(
        &invoker,
        &[
            "summary",
            "-decimals=9",
            "-with-sum",
            "-keep-bed",
            "-header",
            bed_path.to_str().unwrap(),
            bw_path.to_str().unwrap(),
            out_bwtool.to_str().unwrap(),
        ],
    )?;

    // Compare sum/max columns row-wise.
    let tt = read_tsv(&out_tracktools)?;
    let bw = read_tsv(&out_bwtool)?;
    if tt.len() != bw.len() {
        return Err(anyhow!(
            "summary row count mismatch: tracktools {} bwtool {}",
            tt.len(),
            bw.len()
        ));
    }
    let tt_sum_idx = header_col_idx(&tt[0], "sum").ok_or_else(|| anyhow!("tracktools summary: missing sum column"))?;
    let tt_max_idx = header_col_idx(&tt[0], "max").ok_or_else(|| anyhow!("tracktools summary: missing max column"))?;
    let bw_sum_idx = header_col_idx(&bw[0], "sum").ok_or_else(|| anyhow!("bwtool summary: missing sum column"))?;
    let bw_max_idx = header_col_idx(&bw[0], "max").ok_or_else(|| anyhow!("bwtool summary: missing max column"))?;
    for i in 1..tt.len() {
        let sum_tt: f64 = tt[i][tt_sum_idx].parse()?;
        let max_tt: f64 = tt[i][tt_max_idx].parse()?;
        let sum_bw: f64 = bw[i][bw_sum_idx].parse()?;
        let max_bw: f64 = bw[i][bw_max_idx].parse()?;
        if !approx_eq(sum_tt, sum_bw, 1e-6) || !approx_eq(max_tt, max_bw, 1e-6) {
            return Err(anyhow!(
                "bwtool summary mismatch row {i}: tt(sum={sum_tt},max={max_tt}) bw(sum={sum_bw},max={max_bw})"
            ));
        }
    }

    // matrix
    let up = 200u32;
    let down = 200u32;
    let binsize = 50u32;
    let status = Command::new(tracktools_exe())
        .args([
            "matrix",
            "-starts",
            format!("-tiled-averages={binsize}").as_str(),
            format!("{up}:{down}").as_str(),
            bed_path.to_str().unwrap(),
            bw_path.to_str().unwrap(),
            mat_tracktools.to_str().unwrap(),
        ])
        .status()
        .context("run tracktools matrix")?;
    if !status.success() {
        return Err(anyhow!("tracktools matrix failed"));
    }

    let tiled = format!("-tiled-averages={binsize}");
    let size = format!("{up}:{down}");
    run_bwtool(
        &invoker,
        &[
            "matrix",
            "-decimals=9",
            "-starts",
            tiled.as_str(),
            size.as_str(),
            bed_path.to_str().unwrap(),
            bw_path.to_str().unwrap(),
            mat_bwtool.to_str().unwrap(),
        ],
    )?;

    let tt_m = read_tsv(&mat_tracktools)?;
    let bw_m = read_tsv(&mat_bwtool)?;
    if tt_m.len() != bw_m.len() {
        return Err(anyhow!(
            "matrix row count mismatch: tracktools {} bwtool {}",
            tt_m.len(),
            bw_m.len()
        ));
    }
    for (ri, (a, b)) in tt_m.iter().zip(bw_m.iter()).enumerate() {
        if a.len() != b.len() {
            return Err(anyhow!(
                "matrix col count mismatch row {ri}: tracktools {} bwtool {}",
                a.len(),
                b.len()
            ));
        }
        for (ci, (av, bv)) in a.iter().zip(b.iter()).enumerate() {
            let av: f64 = av.parse()?;
            let bv: f64 = bv.parse()?;
            if !approx_eq(av, bv, 1e-6) {
                return Err(anyhow!(
                    "bwtool matrix mismatch row {ri} col {ci}: tt={av} bw={bv}"
                ));
            }
        }
    }

    Ok(())
}

#[test]
fn track_extract_gene_max_matches_bwtool_summary_when_available() -> Result<()> {
    let enable = std::env::var("CI").is_ok()
        || std::env::var("TRACKTOOLS_TEST_ENABLE_GENE")
            .ok()
            .filter(|v| !v.trim().is_empty())
            .is_some();
    if !enable {
        eprintln!("gene-mode test disabled (set CI=true or TRACKTOOLS_TEST_ENABLE_GENE=1); skipping.");
        return Ok(());
    }

    let Some(invoker) = detect_bwtool_invoker() else {
        eprintln!("bwtool not available; skipping gene-mode compatibility check.");
        return Ok(());
    };

    let bw_path = test_bigwig()?;
    if !bw_path.exists() {
        eprintln!("missing test bigWig; skipping: {:?}", bw_path);
        return Ok(());
    }

    let gtf = mini_gtf();
    if !gtf.exists() {
        return Err(anyhow!("missing mini gtf for gene tests: {:?}", gtf));
    }

    let genes = [("CD1D", "CD1D"), ("SLC19A1", "SLC19A1")];
    for (label, gene) in genes {
        let td = TempDir::new("tracktools_test_gene_extract")?;
        let out_dir = td.path.join("out");
        fs::create_dir_all(&out_dir)?;

        let status = Command::new(tracktools_exe())
            .args([
                "track-extract",
                "--out-dir",
                out_dir.to_str().unwrap(),
                "--binsize",
                "200",
                "--no-cytoband",
                "--gtf",
                gtf.to_str().unwrap(),
                "--gene",
                gene,
                "--bigwig",
                bw_path.to_str().unwrap(),
                "--sample",
                "S1",
            ])
            .status()
            .with_context(|| format!("run tracktools track-extract --gene {gene}"))?;
        if !status.success() {
            return Err(anyhow!("tracktools track-extract failed for gene {gene}"));
        }

        let tracks_path = out_dir.join("tracks.tsv");
        let track_rows = read_tsv(&tracks_path)?;
        if track_rows.len() < 2 {
            return Err(anyhow!("tracks.tsv missing data for gene {gene}"));
        }

        // Build a bed from extracted bins.
        let bed_path = out_dir.join(format!("{label}.bins.bed"));
        {
            let mut w = BufWriter::new(File::create(&bed_path)?);
            for r in track_rows.iter().skip(1) {
                if r.len() < 6 {
                    continue;
                }
                writeln!(&mut w, "{}\t{}\t{}", r[1], r[2], r[3])?;
            }
        }

        let out_bwtool = out_dir.join(format!("{label}.bwtool.summary.tsv"));
        run_bwtool(
            &invoker,
            &[
                "summary",
                "-decimals=9",
                "-keep-bed",
                "-header",
                bed_path.to_str().unwrap(),
                bw_path.to_str().unwrap(),
                out_bwtool.to_str().unwrap(),
            ],
        )?;

        let bw = read_tsv(&out_bwtool)?;
        if bw.len() < 2 {
            return Err(anyhow!("bwtool summary missing data for gene {gene}"));
        }
        let bw_max_idx =
            header_col_idx(&bw[0], "max").ok_or_else(|| anyhow!("bwtool summary: missing max"))?;

        // Compare max per bin row-wise.
        let mut bin_i = 0usize;
        for (i, r) in track_rows.iter().enumerate().skip(1) {
            if r.len() < 6 {
                continue;
            }
            bin_i += 1;
            if bin_i >= bw.len() {
                return Err(anyhow!(
                    "bwtool summary row count mismatch for gene {gene}: tracktools bins {bin_i} bwtool {}",
                    bw.len() - 1
                ));
            }
            let max_tt: f64 = r[5].parse()?;
            let max_bw: f64 = bw[bin_i][bw_max_idx].parse()?;
            if !approx_eq(max_tt, max_bw, 1e-6) {
                return Err(anyhow!(
                    "gene {gene} bin {bin_i} max mismatch: tracktools={max_tt} bwtool={max_bw}"
                ));
            }
        }
    }

    Ok(())
}
