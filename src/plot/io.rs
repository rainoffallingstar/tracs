//! Reading the extraction outputs that feed the renderer.
//!
//! `tracs track-extract` writes `tracks.tsv`, `meta.tsv`, `gene_models.tsv` and
//! `cytoband.tsv`; this module parses them back so the renderer works from the
//! same intermediate files that `--work-dir` exposes for debugging.

use std::collections::BTreeMap;
use std::fs;
use std::path::Path;

use anyhow::{anyhow, Context, Result};

/// One binned signal row from `tracks.tsv`.
#[derive(Clone, Debug, PartialEq)]
pub struct SignalBin {
    pub chromosome: String,
    pub start: u64,
    pub end: u64,
    /// Bin width in bp (the `size` column).
    pub size: u64,
    /// Max signal over the bin.
    pub max: f64,
}

/// All bins for one sample, in file order.
#[derive(Clone, Debug, PartialEq)]
pub struct SampleTrack {
    pub sample: String,
    pub bins: Vec<SignalBin>,
}

impl SampleTrack {
    /// Largest `max` across bins, ignoring non-finite values.
    ///
    /// Mirrors `track_plot()`'s `max(x$max, na.rm = TRUE)`, which drives the
    /// per-track y-axis when `groupAutoScale` is off.
    pub fn max_signal(&self) -> f64 {
        self.bins
            .iter()
            .map(|bin| bin.max)
            .filter(|value| value.is_finite())
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// Smallest `max` across bins, ignoring non-finite values.
    ///
    /// `track_plot()` sets each panel's limits to `c(min(x$max), max(x$max))`, so
    /// the lower bound comes from the same column as the upper one.
    pub fn min_signal(&self) -> f64 {
        self.bins
            .iter()
            .map(|bin| bin.max)
            .filter(|value| value.is_finite())
            .fold(f64::INFINITY, f64::min)
    }
}

/// The plotted region, from `meta.tsv`.
#[derive(Clone, Debug, PartialEq)]
pub struct Region {
    pub chromosome: String,
    pub start: u64,
    pub end: u64,
    pub binsize: u64,
    /// The `loci` string as written by extraction (`chr:start-end`).
    pub loci: String,
}

/// A transcript collapsed from `gene_models.tsv`, ready to draw.
#[derive(Clone, Debug, PartialEq)]
pub struct Transcript {
    pub chromosome: String,
    pub strand: String,
    pub transcript: String,
    pub gene: String,
    /// Transcript span (min exon start .. max exon end).
    pub start: u64,
    pub end: u64,
    /// Exon intervals, sorted by start.
    pub exons: Vec<(u64, u64)>,
}

/// One cytoband row.
#[derive(Clone, Debug, PartialEq)]
pub struct Cytoband {
    pub start: u64,
    pub end: u64,
    pub stain: String,
    pub color: String,
}

/// One chromHMM segment: a genomic interval with a state name.
#[derive(Clone, Debug, PartialEq)]
pub struct ChromHmmSegment {
    pub start: u64,
    pub end: u64,
    /// State label as written by UCSC, e.g. `1_Active_Promoter`.
    pub name: String,
}

/// One chromHMM track (a single segmentation), ready to draw as a row.
#[derive(Clone, Debug, PartialEq)]
pub struct ChromHmmTrack {
    /// Display name; R strips the `wgEncodeBroadHmm`/`HMM` affixes at draw time.
    pub name: String,
    pub segments: Vec<ChromHmmSegment>,
}

/// Reads a TSV into rows of fields, skipping blank lines.
fn read_tsv_rows(path: &Path) -> Result<Vec<Vec<String>>> {
    let text = fs::read_to_string(path).with_context(|| format!("read {path:?}"))?;
    Ok(text
        .lines()
        .filter(|line| !line.trim().is_empty())
        .map(|line| line.split('\t').map(|field| field.to_string()).collect())
        .collect())
}

/// Parses `tracks.tsv` into per-sample tracks.
///
/// Samples are grouped in first-seen order so the panel order matches the input
/// order given to `track-extract`.
pub fn read_tracks(path: &Path) -> Result<Vec<SampleTrack>> {
    let rows = read_tsv_rows(path)?;
    if rows.is_empty() {
        return Err(anyhow!("tracks file is empty: {path:?}"));
    }

    let header = &rows[0];
    let column = |name: &str| -> Result<usize> {
        header
            .iter()
            .position(|field| field == name)
            .ok_or_else(|| anyhow!("tracks file {path:?} is missing column {name:?}"))
    };
    let sample_column = column("sample")?;
    let chromosome_column = column("chromosome")?;
    let start_column = column("start")?;
    let end_column = column("end")?;
    let size_column = column("size")?;
    let max_column = column("max")?;

    let mut tracks: Vec<SampleTrack> = Vec::new();
    for (line_number, row) in rows.iter().enumerate().skip(1) {
        if row.len() <= max_column {
            return Err(anyhow!(
                "tracks row {} has {} columns, expected at least {}",
                line_number + 1,
                row.len(),
                max_column + 1
            ));
        }
        let context = || format!("tracks row {}", line_number + 1);
        let bin = SignalBin {
            chromosome: row[chromosome_column].clone(),
            start: row[start_column].parse().with_context(context)?,
            end: row[end_column].parse().with_context(context)?,
            size: row[size_column].parse().with_context(context)?,
            max: row[max_column].parse().with_context(context)?,
        };

        let sample = row[sample_column].clone();
        match tracks.last_mut() {
            Some(track) if track.sample == sample => track.bins.push(bin),
            _ => {
                // A sample should be contiguous in the file; a repeat would mean
                // the extraction wrote rows out of order.
                if tracks.iter().any(|track| track.sample == sample) {
                    return Err(anyhow!(
                        "sample {sample:?} appears in non-contiguous rows (tracks row {})",
                        line_number + 1
                    ));
                }
                tracks.push(SampleTrack {
                    sample,
                    bins: vec![bin],
                });
            }
        }
    }

    if tracks.is_empty() {
        return Err(anyhow!("tracks file has no data rows: {path:?}"));
    }
    Ok(tracks)
}

/// Parses `meta.tsv` (a two-column key/value file).
pub fn read_region(path: &Path) -> Result<Region> {
    let rows = read_tsv_rows(path)?;
    let mut map: BTreeMap<String, String> = BTreeMap::new();
    for row in rows.iter().skip(1) {
        if row.len() >= 2 {
            map.insert(row[0].clone(), row[1].clone());
        }
    }
    let get = |key: &str| -> Result<String> {
        map.get(key)
            .cloned()
            .ok_or_else(|| anyhow!("meta file {path:?} is missing key {key:?}"))
    };
    Ok(Region {
        chromosome: get("chr")?,
        start: get("start")?.parse().context("parse meta start")?,
        end: get("end")?.parse().context("parse meta end")?,
        binsize: get("binsize")?.parse().context("parse meta binsize")?,
        loci: get("loci")?,
    })
}

/// Parses `gene_models.tsv` and collapses exon rows into transcripts.
///
/// `track_plot()` collapses by transcript when `collapse_txs` is set, which is
/// the default, so transcripts are assembled here and the caller decides whether
/// to draw them all.
pub fn read_transcripts(path: &Path) -> Result<Vec<Transcript>> {
    let rows = read_tsv_rows(path)?;
    if rows.is_empty() {
        return Ok(Vec::new());
    }

    let header = &rows[0];
    let column = |name: &str| -> Result<usize> {
        header
            .iter()
            .position(|field| field == name)
            .ok_or_else(|| anyhow!("gene models file {path:?} is missing column {name:?}"))
    };
    let chromosome_column = column("chr")?;
    let strand_column = column("strand")?;
    let transcript_column = column("tx")?;
    let gene_column = column("gene")?;
    let exon_start_column = column("exon_start")?;
    let exon_end_column = column("exon_end")?;

    // Preserve first-seen order so track ordering is stable.
    let mut order: Vec<String> = Vec::new();
    let mut by_transcript: BTreeMap<String, Transcript> = BTreeMap::new();

    for (line_number, row) in rows.iter().enumerate().skip(1) {
        if row.len() <= exon_end_column {
            return Err(anyhow!(
                "gene models row {} has {} columns, expected at least {}",
                line_number + 1,
                row.len(),
                exon_end_column + 1
            ));
        }
        let context = || format!("gene models row {}", line_number + 1);
        let exon_start: u64 = row[exon_start_column].parse().with_context(context)?;
        let exon_end: u64 = row[exon_end_column].parse().with_context(context)?;
        let (exon_start, exon_end) = if exon_start <= exon_end {
            (exon_start, exon_end)
        } else {
            (exon_end, exon_start)
        };

        let key = row[transcript_column].clone();
        let entry = by_transcript.entry(key.clone()).or_insert_with(|| {
            order.push(key.clone());
            Transcript {
                chromosome: row[chromosome_column].clone(),
                strand: row[strand_column].clone(),
                transcript: row[transcript_column].clone(),
                gene: row[gene_column].clone(),
                start: exon_start,
                end: exon_end,
                exons: Vec::new(),
            }
        });
        entry.start = entry.start.min(exon_start);
        entry.end = entry.end.max(exon_end);
        entry.exons.push((exon_start, exon_end));
    }

    let mut transcripts: Vec<Transcript> = order
        .into_iter()
        .filter_map(|key| by_transcript.remove(&key))
        .collect();
    for transcript in &mut transcripts {
        transcript.exons.sort_unstable();
        transcript.exons.dedup();
    }
    Ok(transcripts)
}

/// Parses `cytoband.tsv` for the chromosome being plotted.
pub fn read_cytobands(path: &Path, chromosome: &str) -> Result<Vec<Cytoband>> {
    let rows = read_tsv_rows(path)?;
    if rows.is_empty() {
        return Ok(Vec::new());
    }
    let header = &rows[0];
    let chromosome_column = header
        .iter()
        .position(|field| field == "chr")
        .ok_or_else(|| anyhow!("cytoband file {path:?} is missing column \"chr\""))?;
    let start_column = header
        .iter()
        .position(|field| field == "start")
        .ok_or_else(|| anyhow!("cytoband file {path:?} is missing column \"start\""))?;
    let end_column = header
        .iter()
        .position(|field| field == "end")
        .ok_or_else(|| anyhow!("cytoband file {path:?} is missing column \"end\""))?;
    let stain_column = header
        .iter()
        .position(|field| field == "stain")
        .ok_or_else(|| anyhow!("cytoband file {path:?} is missing column \"stain\""))?;
    let color_column = header
        .iter()
        .position(|field| field == "color")
        .ok_or_else(|| anyhow!("cytoband file {path:?} is missing column \"color\""))?;

    let mut bands = Vec::new();
    for row in rows.iter().skip(1) {
        if row.len() <= color_column || row[chromosome_column] != chromosome {
            continue;
        }
        bands.push(Cytoband {
            start: row[start_column].parse().unwrap_or(0),
            end: row[end_column].parse().unwrap_or(0),
            stain: row[stain_column].clone(),
            color: row[color_column].clone(),
        });
    }
    Ok(bands)
}

/// Reads a chromHMM segmentation file (4 columns: chr, start, end, name).
///
/// Both the UCSC-derived files written by `tracs plot --ucsc-chromhmm` and
/// user-supplied BED-like files use this shape. The name column may be
/// `1_Active_Promoter` (UCSC state encoding) or any free-form label.
pub fn read_chromhmm(path: &Path, track_name: &str, chromosome: &str) -> Result<ChromHmmTrack> {
    let text = fs::read_to_string(path).with_context(|| format!("read {path:?}"))?;
    let mut segments = Vec::new();
    let mut saw_header = false;

    for (line_number, line) in text.lines().enumerate() {
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 4 {
            return Err(anyhow!(
                "chromHMM file {path:?} line {} has {} columns, expected 4 (chr, start, end, name)",
                line_number + 1,
                fields.len()
            ));
        }
        // The extraction output has no header, but a user-supplied file might.
        if !saw_header && fields[1].eq_ignore_ascii_case("start") {
            saw_header = true;
            continue;
        }
        saw_header = true;

        if fields[0] != chromosome {
            continue;
        }
        let start: u64 = fields[1].trim().parse().with_context(|| {
            format!("chromHMM file {path:?} line {}: bad start {:?}", line_number + 1, fields[1])
        })?;
        let end: u64 = fields[2].trim().parse().with_context(|| {
            format!("chromHMM file {path:?} line {}: bad end {:?}", line_number + 1, fields[2])
        })?;
        segments.push(ChromHmmSegment {
            start,
            end,
            name: fields[3].trim().to_string(),
        });
    }

    segments.sort_by_key(|segment| segment.start);
    Ok(ChromHmmTrack {
        name: track_name.to_string(),
        segments,
    })
}

/// Reads a simple BED-like file into `(start, end)` intervals.
pub fn read_regions(path: &Path) -> Result<Vec<(u64, u64)>> {
    let text = fs::read_to_string(path).with_context(|| format!("read {path:?}"))?;
    let mut regions = Vec::new();
    for line in text.lines() {
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            continue;
        }
        if let (Ok(start), Ok(end)) = (fields[1].parse::<u64>(), fields[2].parse::<u64>()) {
            regions.push((start, end));
        }
    }
    Ok(regions)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_temp(name: &str, contents: &str) -> std::path::PathBuf {
        let mut path = std::env::temp_dir();
        path.push(format!(
            "tracs_io_test_{}_{}_{name}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap_or_default()
                .as_nanos()
        ));
        let mut file = fs::File::create(&path).expect("create temp");
        file.write_all(contents.as_bytes()).expect("write temp");
        path
    }

    #[test]
    fn parses_tracks_grouped_by_sample() -> Result<()> {
        let path = write_temp(
            "tracks.tsv",
            "sample\tchromosome\tstart\tend\tsize\tmax\n\
             s1\tchr1\t100\t200\t100\t1.5\n\
             s1\tchr1\t200\t300\t100\t2.5\n\
             s2\tchr1\t100\t200\t100\t0.5\n",
        );
        let tracks = read_tracks(&path)?;
        assert_eq!(tracks.len(), 2);
        assert_eq!(tracks[0].sample, "s1");
        assert_eq!(tracks[0].bins.len(), 2);
        assert_eq!(tracks[1].sample, "s2");
        assert_eq!(tracks[1].bins.len(), 1);
        assert_eq!(tracks[0].max_signal(), 2.5);
        assert_eq!(tracks[1].max_signal(), 0.5);
        let _ = fs::remove_file(&path);
        Ok(())
    }

    #[test]
    fn rejects_non_contiguous_sample_rows() {
        let path = write_temp(
            "tracks.tsv",
            "sample\tchromosome\tstart\tend\tsize\tmax\n\
             s1\tchr1\t100\t200\t100\t1\n\
             s2\tchr1\t100\t200\t100\t1\n\
             s1\tchr1\t200\t300\t100\t1\n",
        );
        assert!(
            read_tracks(&path).is_err(),
            "interleaved samples should be rejected"
        );
        let _ = fs::remove_file(&path);
    }

    #[test]
    fn parses_region_from_meta() -> Result<()> {
        let path = write_temp(
            "meta.tsv",
            "key\tvalue\nchr\tchr1\nstart\t1000\nend\t2000\nbinsize\t50\nloci\tchr1:1000-2000\n",
        );
        let region = read_region(&path)?;
        assert_eq!(region.chromosome, "chr1");
        assert_eq!(region.start, 1000);
        assert_eq!(region.end, 2000);
        assert_eq!(region.binsize, 50);
        assert_eq!(region.loci, "chr1:1000-2000");
        let _ = fs::remove_file(&path);
        Ok(())
    }

    #[test]
    fn collapses_exons_into_transcripts() -> Result<()> {
        let path = write_temp(
            "gene_models.tsv",
            "chr\tstart\tend\tstrand\ttx\tgene\texon_start\texon_end\n\
             chr1\t100\t300\t+\tT1\tG1\t100\t150\n\
             chr1\t100\t300\t+\tT1\tG1\t250\t300\n\
             chr1\t100\t400\t-\tT2\tG2\t350\t400\n",
        );
        let transcripts = read_transcripts(&path)?;
        assert_eq!(transcripts.len(), 2);

        let first = &transcripts[0];
        assert_eq!(first.transcript, "T1");
        assert_eq!(first.start, 100);
        assert_eq!(first.end, 300);
        assert_eq!(first.exons, vec![(100, 150), (250, 300)]);

        let second = &transcripts[1];
        assert_eq!(second.strand, "-");
        assert_eq!(second.exons, vec![(350, 400)]);
        let _ = fs::remove_file(&path);
        Ok(())
    }

    #[test]
    fn filters_cytobands_by_chromosome() -> Result<()> {
        let path = write_temp(
            "cytoband.tsv",
            "chr\tstart\tend\tband\tstain\tcolor\n\
             chr1\t0\t100\tp1\tgneg\t#FFFFFF\n\
             chr2\t0\t100\tp1\tgpos25\t#C0C0C0\n\
             chr1\t100\t200\tp2\tgpos50\t#808080\n",
        );
        let bands = read_cytobands(&path, "chr1")?;
        assert_eq!(bands.len(), 2);
        assert_eq!(bands[0].color, "#FFFFFF");
        assert_eq!(bands[1].stain, "gpos50");
        let _ = fs::remove_file(&path);
        Ok(())
    }
}
