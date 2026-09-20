//! Principal component analysis of per-region signal matrices.
//!
//! Ported from `trackplot.R`'s `pca_plot()`, which calls R's `prcomp()` on the
//! transpose of an `extract_summary()` table: rows become samples and columns
//! become the regions whose signal is the feature space. Samples that behave
//! alike therefore land near each other, which is the whole point of the panel.
//!
//! R's `prcomp()` centres each column, then takes an SVD of the centred matrix
//! and reports `sdev = d / sqrt(n - 1)` with scores `x = U * d`. That is
//! reproduced here directly rather than through a generic eigen-solver, because
//! the two details that matter downstream are both consequences of the SVD form:
//!
//! - Variance explained is `d^2 / sum(d^2)`, not `d^2 / sum(lambda)` over a
//!   covariance matrix, so a rank-deficient matrix yields zeros in the tail
//!   instead of dividing by a small number.
//! - R's `prcomp()` returns `min(n_samples, n_features)` components (it passes
//!   `nu = 0`, so the decomposition is economy-sized), which means a wide
//!   summary table produces exactly one component per sample.
//!
//! One thing the port deliberately does *not* try to match is the **sign** of
//! each component. R's own documentation says the signs "are arbitrary ... and
//! may differ between different programs for PCA, and even between different
//! builds of R", so a sign is chosen here for stability (the component's
//! largest-magnitude entry is made positive) and the oracle test compares
//! absolute values. Plotting is unaffected, and `--flip-sign` lets a caller
//! match a specific figure by hand.

pub mod panel;

/// Number of components that `prcomp()` would return for a sample-by-feature
/// matrix of this shape.
pub fn component_count(n_samples: usize, n_features: usize) -> usize {
    n_samples.min(n_features)
}

/// One principal component: the direction it points in feature space, plus the
/// spread along it.
#[derive(Clone, Debug, PartialEq)]
pub struct Component {
    /// Singular values scaled the way R scales them: `d / sqrt(n - 1)`.
    pub sdev: f64,
    /// Share of the total variance carried by this component.
    pub variance_explained: f64,
    /// The component's coordinates, one per feature (a column of
    /// `prcomp()$rotation`).
    ///
    /// Unit length for a component that carries variance. A component whose
    /// eigenvalue falls below the rank tolerance has no defined direction, so
    /// every entry is zero rather than an arbitrary unit vector.
    pub loadings: Vec<f64>,
    /// Sample scores along this component (a column of `prcomp()$x`).
    pub scores: Vec<f64>,
}

/// The result of a PCA over samples.
#[derive(Clone, Debug, PartialEq)]
pub struct Pca {
    pub components: Vec<Component>,
    /// The number of sampled features (regions) that went into the fit.
    pub n_features: usize,
    /// The number of samples that went into the fit.
    pub n_samples: usize,
}

impl Pca {
    /// Variance explained by component `index`, or 0.0 when out of range.
    pub fn variance_explained(&self, index: usize) -> f64 {
        self.components
            .get(index)
            .map(|component| component.variance_explained)
            .unwrap_or(0.0)
    }

    /// `sdev` of component `index`, or 0.0 when out of range.
    pub fn sdev(&self, index: usize) -> f64 {
        self.components
            .get(index)
            .map(|component| component.sdev)
            .unwrap_or(0.0)
    }
}

/// A square matrix in row-major order, used for the symmetric eigenvalue solve.
struct SquareMatrix {
    size: usize,
    values: Vec<f64>,
}

impl SquareMatrix {
    fn zeros(size: usize) -> Self {
        Self {
            size,
            values: vec![0.0; size * size],
        }
    }

    fn identity(size: usize) -> Self {
        let mut matrix = Self::zeros(size);
        for index in 0..size {
            matrix.set(index, index, 1.0);
        }
        matrix
    }

    fn get(&self, row: usize, column: usize) -> f64 {
        self.values[row * self.size + column]
    }

    fn set(&mut self, row: usize, column: usize, value: f64) {
        self.values[row * self.size + column] = value;
    }
}

/// Centers each feature (column) across samples, as `prcomp()` does.
///
/// Missing cells stay `NaN` rather than becoming an imputed zero, so a region
/// that is absent for one sample cannot drag the other samples towards it.
fn center_columns(rows: &[Vec<f64>], n_features: usize) -> Vec<Vec<f64>> {
    let mut feature_means = vec![0.0f64; n_features];
    for (column, mean) in feature_means.iter_mut().enumerate() {
        let mut sum = 0.0;
        let mut count = 0usize;
        for row in rows {
            if let Some(value) = row.get(column).filter(|value| value.is_finite()) {
                sum += value;
                count += 1;
            }
        }
        *mean = if count == 0 {
            0.0
        } else {
            sum / count as f64
        };
    }

    rows.iter()
        .map(|row| {
            (0..n_features)
                .map(|column| match row.get(column) {
                    Some(value) if value.is_finite() => value - feature_means[column],
                    _ => f64::NAN,
                })
                .collect()
        })
        .collect()
}

/// Computes the cross-product matrix `X * X^T` of already-centered data.
///
/// Working in sample space (an `n_samples x n_samples` matrix) rather than
/// feature space is what keeps the cost independent of how many regions the
/// summary table holds, which in practice is tens of thousands.
fn centered_gram(centered: &[Vec<f64>], n_features: usize) -> SquareMatrix {
    let n_samples = centered.len();
    let mut gram = SquareMatrix::zeros(n_samples);

    for left in 0..n_samples {
        for right in left..n_samples {
            // Zipping the two centered rows keeps the feature index implicit,
            // which is what the product actually needs: the same column of both.
            let mut sum = 0.0;
            for (left_value, right_value) in centered[left]
                .iter()
                .zip(centered[right].iter())
                .take(n_features)
            {
                // A feature counts only when both samples have it, so a missing
                // cell contributes nothing instead of a fabricated zero.
                if left_value.is_finite() && right_value.is_finite() {
                    sum += left_value * right_value;
                }
            }
            gram.set(left, right, sum);
            gram.set(right, left, sum);
        }
    }

    gram
}

/// Jacobi eigenvalue iteration for a symmetric matrix.
///
/// Returns eigenvalues in descending order with their eigenvectors as columns.
/// Cyclic Jacobi is used because it converges reliably on the small, dense,
/// possibly rank-deficient matrices that a sample-space PCA produces; the
/// sample count is tiny next to the region count, so an O(n^3) sweep is cheap.
fn symmetric_eigen(matrix: &SquareMatrix) -> (Vec<f64>, Vec<Vec<f64>>) {
    let size = matrix.size;
    let mut a = SquareMatrix {
        size,
        values: matrix.values.clone(),
    };
    let mut eigenvectors = SquareMatrix::identity(size);

    // The largest diagonal entry sets the scale for the convergence threshold:
    // an absolute cutoff would never converge on large-magnitude matrices and
    // would stop too early on small ones.
    let scale = (0..size)
        .map(|index| a.get(index, index).abs())
        .fold(0.0f64, f64::max)
        .max(f64::MIN_POSITIVE);
    let tolerance = scale * 1e-15;

    for _ in 0..100 {
        let mut largest_off_diagonal = 0.0f64;
        for row in 0..size {
            for column in (row + 1)..size {
                largest_off_diagonal = largest_off_diagonal.max(a.get(row, column).abs());
            }
        }
        if largest_off_diagonal <= tolerance {
            break;
        }

        for pivot in 0..size {
            for second in (pivot + 1)..size {
                let off_diagonal = a.get(pivot, second);
                if off_diagonal.abs() <= tolerance {
                    continue;
                }
                let difference = a.get(second, second) - a.get(pivot, pivot);
                // The rotation that zeroes this entry; the sign choice keeps the
                // computed angle in the well-conditioned branch.
                let theta = difference / (2.0 * off_diagonal);
                let t = if theta >= 0.0 {
                    1.0 / (theta + (theta * theta + 1.0).sqrt())
                } else {
                    -1.0 / (-theta + (theta * theta + 1.0).sqrt())
                };
                let cosine = 1.0 / (t * t + 1.0).sqrt();
                let sine = t * cosine;

                for index in 0..size {
                    let pivot_value = a.get(index, pivot);
                    let second_value = a.get(index, second);
                    a.set(index, pivot, cosine * pivot_value - sine * second_value);
                    a.set(index, second, sine * pivot_value + cosine * second_value);
                }
                for index in 0..size {
                    let pivot_value = a.get(pivot, index);
                    let second_value = a.get(second, index);
                    a.set(pivot, index, cosine * pivot_value - sine * second_value);
                    a.set(second, index, sine * pivot_value + cosine * second_value);
                }
                for index in 0..size {
                    let pivot_value = eigenvectors.get(index, pivot);
                    let second_value = eigenvectors.get(index, second);
                    eigenvectors.set(
                        index,
                        pivot,
                        cosine * pivot_value - sine * second_value,
                    );
                    eigenvectors.set(
                        index,
                        second,
                        sine * pivot_value + cosine * second_value,
                    );
                }
            }
        }
    }

    let mut order: Vec<usize> = (0..size).collect();
    order.sort_by(|left, right| {
        a.get(*right, *right)
            .partial_cmp(&a.get(*left, *left))
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let eigenvalues = order.iter().map(|index| a.get(*index, *index)).collect();
    let eigenvectors_by_column: Vec<Vec<f64>> = (0..size)
        .map(|column| {
            (0..size)
                .map(|row| eigenvectors.get(row, order[column]))
                .collect()
        })
        .collect();

    (eigenvalues, eigenvectors_by_column)
}

/// Runs the PCA over samples, reproducing `prcomp(t(summary_table))`.
///
/// `rows` holds one entry per sample, each a vector of per-region values in a
/// shared feature order; that is the transpose of the `extract_summary()` table
/// R feeds to `prcomp()`, so it is already the orientation `prcomp()` sees.
///
/// Returns `None` when there is nothing to decompose (fewer than two samples or
/// no features).
pub fn fit(rows: &[Vec<f64>]) -> Option<Pca> {
    let n_samples = rows.len();
    let n_features = rows.iter().map(Vec::len).max().unwrap_or(0);
    if n_samples < 2 || n_features == 0 {
        return None;
    }

    let centered = center_columns(rows, n_features);
    let gram = centered_gram(&centered, n_features);
    let (eigenvalues, eigenvectors) = symmetric_eigen(&gram);

    // `prcomp()` keeps min(n_samples, n_features) components and normalises by
    // the total variance over exactly those components, so the shares sum to 1.
    let kept = component_count(n_samples, n_features);

    // A Jacobi sweep converges to roughly `eps * ||gram||` in absolute terms, so
    // the null space of a rank-deficient matrix comes out as noise around 1e-12
    // rather than as exact zeros. LAPACK (and therefore R) reports exact zeros
    // there, and callers rely on that: `pca_plot()` labels an axis "PC3 [0]" and
    // the scree plot draws a zero bar. Clamping below the same rank tolerance
    // LAPACK uses restores that, and keeps `variance_explained` from reporting a
    // meaningless share for a direction that carries no variance.
    let largest_eigenvalue = eigenvalues
        .iter()
        .take(kept)
        .copied()
        .fold(0.0f64, f64::max);
    let rank_tolerance = largest_eigenvalue * n_samples as f64 * f64::EPSILON;

    let total_variance: f64 = eigenvalues
        .iter()
        .take(kept)
        .copied()
        .filter(|value| *value > rank_tolerance)
        .sum();

    let mut components = Vec::with_capacity(kept);
    for index in 0..kept {
        let raw_eigenvalue = eigenvalues.get(index).copied().unwrap_or(0.0);
        let eigenvalue = if raw_eigenvalue > rank_tolerance {
            raw_eigenvalue
        } else {
            0.0
        };
        let singular = eigenvalue.sqrt();
        // R: `s$d <- s$d / sqrt(max(1, n - 1))`.
        let sdev = singular / (n_samples.saturating_sub(1).max(1) as f64).sqrt();
        let variance_explained = if total_variance > 0.0 {
            eigenvalue / total_variance
        } else {
            0.0
        };

        // Sample scores are `U * d`, and the eigenvectors of the Gram matrix are
        // exactly `U` (up to sign).
        let mut sample_scores: Vec<f64> = eigenvectors
            .get(index)
            .map(|column| column.iter().map(|value| value * singular).collect())
            .unwrap_or_else(|| vec![0.0; n_samples]);

        // Loadings come from projecting the centered data back onto the sample
        // scores, which yields the feature-space direction without ever forming
        // the feature-space covariance matrix.
        //
        // Projecting the *centered* data is what makes a constant region load
        // exactly zero, and it is what keeps loadings consistent with the
        // centered scores they exist to explain.
        let mut loadings = vec![0.0f64; n_features];
        if eigenvalue > 0.0 {
            for (column, loading) in loadings.iter_mut().enumerate() {
                let mut sum = 0.0;
                for (sample, row) in centered.iter().enumerate() {
                    let value = row[column];
                    if value.is_finite() && sample_scores[sample].is_finite() {
                        sum += sample_scores[sample] * value;
                    }
                }
                *loading = sum / eigenvalue;
            }
        }

        // Fix an arbitrary sign convention. R leaves it to LAPACK, and its docs
        // note the signs can differ between builds, so tests compare magnitudes.
        // Anchoring on the largest loading keeps the panel stable when the same
        // data is plotted twice.
        let anchor = loadings
            .iter()
            .enumerate()
            .filter(|(_, value)| value.is_finite())
            .max_by(|left, right| {
                left.1
                    .abs()
                    .partial_cmp(&right.1.abs())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(_, value)| *value)
            .unwrap_or(0.0);
        if anchor < 0.0 {
            for value in loadings.iter_mut() {
                *value = -*value;
            }
            for value in sample_scores.iter_mut() {
                *value = -*value;
            }
        }

        components.push(Component {
            sdev,
            variance_explained,
            loadings,
            scores: sample_scores,
        });
    }

    Some(Pca {
        components,
        n_features,
        n_samples,
    })
}

/// A summary table in the orientation `extract_summary()` produces: one row per
/// region, one column per sample.
///
/// This is what `pca_plot()` receives, and it is the transpose of the
/// sample-major matrix `fit()` takes, so [`samples_from_columns`] bridges them.
#[derive(Clone, Debug, PartialEq)]
pub struct SummaryTable {
    /// Per-region values, one sample wide.
    pub regions: Vec<Vec<f64>>,
    pub sample_names: Vec<String>,
}

impl SummaryTable {
    /// Number of samples (columns) in the table.
    pub fn n_samples(&self) -> usize {
        self.sample_names.len()
    }

    /// Number of regions (rows) in the table.
    pub fn n_regions(&self) -> usize {
        self.regions.len()
    }

    /// Assembles a table from one value column per sample, which is how
    /// `extract_summary()` combines the per-bigWig `bwtool summary` outputs: it
    /// keeps the `sum` column from each and `cbind`s them.
    ///
    /// Columns may differ in length (a region missing from one summary file);
    /// shorter columns are padded with `NaN` so the region count is preserved
    /// and the missing cell stays missing rather than becoming a zero.
    pub fn from_sample_columns(columns: &[(String, Vec<f64>)]) -> Self {
        let n_regions = columns.iter().map(|(_, values)| values.len()).max().unwrap_or(0);
        let sample_names = columns
            .iter()
            .map(|(name, _)| name.clone())
            .collect();
        let regions = (0..n_regions)
            .map(|region| {
                columns
                    .iter()
                    .map(|(_, values)| values.get(region).copied().unwrap_or(f64::NAN))
                    .collect()
            })
            .collect();
        Self {
            regions,
            sample_names,
        }
    }
}

/// R's `sd()`, including its `n - 1` denominator and `na.rm` handling.
///
/// Matching the denominator matters here: `pca_plot()` ranks regions by this
/// value to decide which ones enter the PCA, so a population (n-denominator)
/// standard deviation would pick a different — if usually similar — subset.
pub fn region_sd(values: &[f64]) -> f64 {
    let finite: Vec<f64> = values
        .iter()
        .copied()
        .filter(|value| value.is_finite())
        .collect();
    if finite.len() < 2 {
        return f64::NAN;
    }
    let mean = finite.iter().sum::<f64>() / finite.len() as f64;
    let sum_squares = finite
        .iter()
        .map(|value| (value - mean) * (value - mean))
        .sum::<f64>();
    (sum_squares / (finite.len() - 1) as f64).sqrt()
}

/// Applies `pca_plot()`'s optional log transform in place: `log2(x + 0.1)`.
///
/// The offset keeps zeros finite, and it is what R's default `log2 = FALSE`
/// leaves to the caller. Non-finite entries stay non-finite so they continue to
/// count as missing.
pub fn log2_transform(table: &mut SummaryTable, offset: f64) {
    for region in table.regions.iter_mut() {
        for value in region.iter_mut() {
            if value.is_finite() {
                *value = (*value + offset).log2();
            }
        }
    }
}

/// Orders region indices by descending standard deviation, exactly as
/// `pca_plot()`'s `order(apply(sum_tbl, 1, sd), decreasing = TRUE, na.last = TRUE)`.
///
/// Regions whose standard deviation is `NaN` (fewer than two samples, or all
/// missing) sort last, matching `na.last = TRUE`. Ties keep their original
/// order, since R's `order()` is stable.
pub fn order_regions_by_sd(table: &SummaryTable) -> Vec<usize> {
    let mut indices: Vec<usize> = (0..table.regions.len()).collect();
    indices.sort_by(|left, right| {
        let left_sd = region_sd(&table.regions[*left]);
        let right_sd = region_sd(&table.regions[*right]);
        match (left_sd.is_nan(), right_sd.is_nan()) {
            (true, true) => left.cmp(right),
            (true, false) => std::cmp::Ordering::Greater,
            (false, true) => std::cmp::Ordering::Less,
            (false, false) => right_sd
                .partial_cmp(&left_sd)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then(left.cmp(right)),
        }
    });
    indices
}

/// Keeps the `top` most variable regions, in the order they were ranked.
///
/// R keeps *all* rows when there are fewer than `top` of them, so this is a
/// no-op in that case; otherwise it truncates after ranking. Returning the
/// ranking as well lets a caller report which regions were used.
pub fn select_top_regions(table: &SummaryTable, top: usize) -> (SummaryTable, Vec<usize>) {
    let order = order_regions_by_sd(table);
    if table.regions.len() <= top {
        return (table.clone(), order);
    }
    let kept: Vec<usize> = order.iter().copied().take(top).collect();
    let regions = kept
        .iter()
        .map(|index| table.regions[*index].clone())
        .collect();
    (
        SummaryTable {
            regions,
            sample_names: table.sample_names.clone(),
        },
        order,
    )
}

/// Runs `fit()` on a summary table, applying the same preprocessing
/// `pca_plot()` does: optional log transform, rank regions by variance, keep the
/// top ones, then decompose the transpose.
///
/// Returns the fitted PCA together with the region ranking, so a caller can
/// report which regions were dropped.
pub fn fit_summary_table(
    table: &SummaryTable,
    top: usize,
    log2: bool,
    log2_offset: f64,
) -> Option<(Pca, Vec<usize>)> {
    let mut table = table.clone();
    if log2 {
        log2_transform(&mut table, log2_offset);
    }
    let (selected, ranking) = select_top_regions(&table, top);
    let rows = samples_from_columns(&selected.regions, selected.n_samples());
    fit(&rows).map(|pca| (pca, ranking))
}

/// Turns a regions-by-samples table into the samples-by-regions rows `fit()`
/// wants.
///
/// This is the transpose that `pca_plot()` performs when it hands
/// `t(sum_tbl)` to `prcomp()`: the summary table has one row per region and one
/// column per sample, so reading it "by column" is what produces the
/// sample-major matrix.
pub fn samples_from_columns(columns: &[Vec<f64>], n_samples: usize) -> Vec<Vec<f64>> {
    (0..n_samples)
        .map(|sample| {
            columns
                .iter()
                .map(|column| column.get(sample).copied().unwrap_or(f64::NAN))
                .collect()
        })
        .collect()
}

#[cfg(test)]
// The expected values below are copied verbatim from R's own output at 17
// significant digits. Trimming them to the shortest round-tripping literal
// would obscure that they came from R, so clippy's precision lints are allowed
// here rather than applying their suggestions.
#[allow(clippy::excessive_precision, clippy::approx_constant)]
mod tests {
    use super::*;

    /// Case A of the oracle in sample-major order, which is what `fit()` takes.
    ///
    /// The oracle's matrix A has the regions `[1,2,3,4,5]`, `[5,4,3,2,1]` and
    /// `[2,2,4,4,6]` as its *columns*, because `pca_plot()` feeds
    /// `prcomp(t(sum_tbl))`. Transposing it gives one row per sample.
    fn case_a() -> Vec<Vec<f64>> {
        vec![
            vec![1.0, 2.0, 3.0, 4.0, 5.0],
            vec![5.0, 4.0, 3.0, 2.0, 1.0],
            vec![2.0, 2.0, 4.0, 4.0, 6.0],
        ]
    }

    /// Rebuilds the sample-major matrix from a region-major oracle case.
    fn cases(matrix: &[Vec<f64>]) -> Vec<Vec<f64>> {
        let n_samples = matrix[0].len();
        samples_from_columns(matrix, n_samples)
    }

    #[test]
    fn component_count_matches_prcomp_economy_size() {
        // prcomp() passes nu = 0, so it returns min(n_samples, n_features).
        assert_eq!(component_count(3, 5), 3);
        assert_eq!(component_count(5, 3), 3);
        assert_eq!(component_count(4, 4), 4);
    }

    #[test]
    fn fit_requires_two_samples_and_one_feature() {
        assert!(fit(&[]).is_none());
        assert!(fit(&[vec![1.0, 2.0]]).is_none());
        assert!(fit(&[vec![], vec![]]).is_none());
    }

    #[test]
    fn sdev_and_variance_explained_match_r_case_a() {
        let pca = fit(&case_a()).expect("fit");

        assert_eq!(pca.n_samples, 3);
        assert_eq!(pca.n_features, 5);
        assert_eq!(pca.components.len(), 3);

        // R: sdev = 3.6875367289278058, 0.85755828148397784.
        assert!((pca.sdev(0) - 3.6875367289278058).abs() < 1e-12);
        assert!((pca.sdev(1) - 0.85755828148397784).abs() < 1e-12);
        // R: var explained = 0.94869259026917996, 0.051307409730819993.
        assert!((pca.variance_explained(0) - 0.94869259026917996).abs() < 1e-12);
        assert!((pca.variance_explained(1) - 0.051307409730819993).abs() < 1e-12);
    }

    #[test]
    fn variance_explained_sums_to_one() {
        let pca = fit(&case_a()).expect("fit");
        let total: f64 = pca
            .components
            .iter()
            .map(|component| component.variance_explained)
            .sum();
        assert!(
            (total - 1.0).abs() < 1e-12,
            "shares should sum to 1, got {total}"
        );
    }

    #[test]
    fn scores_match_r_up_to_component_sign() {
        let pca = fit(&case_a()).expect("fit");

        // R reports PC1 as (-2.0035819573047626, 4.2555844005622134, -2.2520024432574508).
        // The overall sign is arbitrary in R, so it is matched here before comparing.
        let expected = [-2.0035819573047626, 4.2555844005622134, -2.2520024432574508];
        let sign = if pca.components[0].scores[0] * expected[0] < 0.0 {
            -1.0
        } else {
            1.0
        };
        for (actual, want) in pca.components[0].scores.iter().zip(expected.iter()) {
            assert!(
                (actual * sign - want).abs() < 1e-12,
                "PC1 score {actual} vs R {want}"
            );
        }

        // R reports PC2 as (-0.8737488873472723, 0.033354472003946034, 0.84039441534332571).
        let expected = [-0.8737488873472723, 0.033354472003946034, 0.84039441534332571];
        let sign = if pca.components[1].scores[0] * expected[0] < 0.0 {
            -1.0
        } else {
            1.0
        };
        for (actual, want) in pca.components[1].scores.iter().zip(expected.iter()) {
            assert!(
                (actual * sign - want).abs() < 1e-12,
                "PC2 score {actual} vs R {want}"
            );
        }
    }

    #[test]
    fn scores_are_centered() {
        let pca = fit(&case_a()).expect("fit");
        for component in &pca.components {
            let mean = component.scores.iter().sum::<f64>() / component.scores.len() as f64;
            assert!(mean.abs() < 1e-12, "scores should be centred, mean {mean}");
        }
    }

    #[test]
    fn loadings_are_unit_length_where_defined() {
        // `prcomp()$rotation` is orthonormal, so each loading vector has norm 1.
        // A component with no variance has no defined direction, so its entries
        // are all zero instead of an arbitrary unit vector.
        let pca = fit(&case_a()).expect("fit");
        for (index, component) in pca.components.iter().enumerate() {
            let norm = component
                .loadings
                .iter()
                .map(|value| value * value)
                .sum::<f64>()
                .sqrt();
            if component.sdev > 1e-12 {
                assert!(
                    (norm - 1.0).abs() < 1e-10,
                    "component {index} loadings have norm {norm}"
                );
            } else {
                assert_eq!(norm, 0.0, "component {index} should be all zeros");
            }
        }
    }

    #[test]
    fn a_constant_region_gets_zero_loading() {
        // Oracle case F: regions 1 and 4 are constant across samples, so they
        // carry no variance and cannot load on any component.
        let matrix = vec![
            vec![4.0, 4.0, 4.0],
            vec![1.0, 2.0, 6.0],
            vec![9.0, 6.0, 2.0],
            vec![5.0, 5.0, 5.0],
        ];
        let pca = fit(&cases(&matrix)).expect("fit");

        assert!((pca.sdev(0) - 4.3650221717900983).abs() < 1e-12);
        assert!((pca.sdev(1) - 0.5290697242464214).abs() < 1e-12);
        for component in &pca.components {
            for constant_index in [0usize, 3] {
                assert!(
                    component.loadings[constant_index].abs() < 1e-12,
                    "constant region {constant_index} loaded {} on a component",
                    component.loadings[constant_index]
                );
            }
        }
    }

    #[test]
    fn collinear_regions_leave_only_one_component_with_variance() {
        // Oracle case C: every region is a fixed offset of the others, so the
        // samples lie on a line and only PC1 carries variance.
        let matrix = vec![
            vec![1.0, 4.0, 7.0],
            vec![2.0, 5.0, 8.0],
            vec![3.0, 6.0, 9.0],
        ];
        let pca = fit(&cases(&matrix)).expect("fit");

        assert!((pca.sdev(0) - 5.196152422706632).abs() < 1e-12);
        assert!(pca.sdev(1) < 1e-12, "PC2 should be degenerate");
        assert!((pca.variance_explained(0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn more_samples_than_regions_keeps_one_component_per_region() {
        // Oracle case E: three regions and five samples, so prcomp keeps three
        // components and the tail is empty rather than zero-padded.
        let matrix = vec![
            vec![1.0, 2.0, 3.0, 4.0, 5.0],
            vec![6.0, 7.0, 8.0, 9.0, 10.0],
            vec![11.0, 12.0, 13.0, 14.0, 15.0],
        ];
        let pca = fit(&cases(&matrix)).expect("fit");

        assert_eq!(pca.n_samples, 5);
        assert_eq!(pca.n_features, 3);
        assert_eq!(pca.components.len(), 3);
        assert!((pca.sdev(0) - 2.7386127875258306).abs() < 1e-12);
        assert!(pca.sdev(1) < 1e-14, "PC2 should be degenerate");
    }

    #[test]
    fn identical_samples_are_degenerate_but_finite() {
        // Oracle case G: two of three samples are identical, so every component
        // after PC1 collapses. The port must report zeros, not NaN or infinity.
        let matrix = vec![
            vec![1.0, 1.0, 9.0],
            vec![3.0, 3.0, 8.0],
            vec![5.0, 5.0, 7.0],
            vec![7.0, 7.0, 6.0],
            vec![2.0, 2.0, 5.0],
        ];
        let pca = fit(&cases(&matrix)).expect("fit");

        assert!((pca.sdev(0) - 5.8594652770823146).abs() < 1e-12);
        assert!(pca.sdev(1) < 1e-12, "PC2 should be degenerate");
        for component in &pca.components {
            for value in component.loadings.iter().chain(component.scores.iter()) {
                assert!(value.is_finite(), "degenerate fit produced {value}");
            }
        }
    }

    #[test]
    fn a_flat_matrix_reports_no_variance_rather_than_dividing_by_zero() {
        let matrix = vec![vec![7.0, 7.0, 7.0], vec![7.0, 7.0, 7.0]];
        let pca = fit(&cases(&matrix)).expect("fit");
        for component in &pca.components {
            assert_eq!(component.variance_explained, 0.0);
        }
    }

    #[test]
    fn two_samples_produce_exactly_two_components() {
        // Oracle case D.
        let matrix = vec![
            vec![10.0, 11.0],
            vec![20.0, 21.0],
            vec![30.0, 31.0],
            vec![40.0, 41.0],
        ];
        let pca = fit(&cases(&matrix)).expect("fit");

        assert_eq!(pca.components.len(), 2);
        assert!((pca.sdev(0) - 1.4142135623730951).abs() < 1e-12);
        assert!((pca.variance_explained(0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn samples_from_columns_transposes_and_pads_missing_values() {
        let columns = vec![vec![1.0, 2.0, 3.0], vec![10.0, 20.0]];
        let samples = samples_from_columns(&columns, 3);

        assert_eq!(samples.len(), 3);
        assert_eq!(samples[0], vec![1.0, 10.0]);
        assert_eq!(samples[1], vec![2.0, 20.0]);
        assert!(samples[2][1].is_nan(), "a short column should pad with NaN");
    }

    #[test]
    fn missing_values_do_not_produce_nan_scores() {
        let columns = vec![
            vec![1.0, 2.0, f64::NAN],
            vec![3.0, 4.0, 5.0],
            vec![6.0, 7.0, 8.0],
        ];
        let pca = fit(&samples_from_columns(&columns, 3)).expect("fit");
        for component in &pca.components {
            for value in &component.scores {
                assert!(value.is_finite(), "missing data produced {value}");
            }
        }
    }

    // ---------------------------------------------------------------------
    // R oracle
    //
    // `testdata/pca_r_oracle.tsv` was produced by R 4.6.0 and holds, per case:
    //
    //   fm <- function(v) paste(sprintf("%.17g", v), collapse = ",")
    //   cases <- list(
    //     A = rbind(c(1, 5, 2), c(2, 4, 2), c(3, 3, 4), c(4, 2, 4), c(5, 1, 6)),
    //     B = rbind(c(0, 7, 1, 2), c(1, 6, 3, 2), c(2, 5, 5, 2), c(3, 4, 7, 2),
    //               c(4, 3, 9, 2), c(5, 2, 11, 2), c(6, 1, 13, 2), c(7, 0, 15, 2)),
    //     C = rbind(c(1, 4, 7), c(2, 5, 8), c(3, 6, 9)),
    //     D = rbind(c(10, 11), c(20, 21), c(30, 31), c(40, 41)),
    //     E = rbind(c(1, 2, 3, 4, 5), c(6, 7, 8, 9, 10), c(11, 12, 13, 14, 15)),
    //     F = rbind(c(4, 4, 4), c(1, 2, 6), c(9, 6, 2), c(5, 5, 5)),
    //     G = rbind(c(1, 1, 9), c(3, 3, 8), c(5, 5, 7), c(7, 7, 6), c(2, 2, 5))
    //   )
    //   for (name in names(cases)) {
    //     m <- cases[[name]]; p <- prcomp(t(m))
    //     cat(sprintf("%s\t%d\t%d\t%s\t%s\t%s\t%s\n", name, nrow(m), ncol(m),
    //       fm(as.vector(t(m))), fm(p$sdev),
    //       fm(p$sdev^2 / sum(p$sdev^2)), fm(as.vector(t(p$x)))))
    //   }
    //
    // Regenerate with the snippet above when the oracle needs extending.
    // ---------------------------------------------------------------------

    /// Splits a comma-separated numeric field.
    fn parse_numbers(field: &str) -> Vec<f64> {
        if field.is_empty() {
            return Vec::new();
        }
        field
            .split(',')
            .map(|token| token.trim().parse::<f64>().unwrap_or(f64::NAN))
            .collect()
    }

    /// Relative tolerance for values R prints with 17 significant digits.
    fn oracle_close(actual: f64, expected: f64) -> bool {
        if actual == expected {
            return true;
        }
        let scale = actual.abs().max(expected.abs()).max(1.0);
        (actual - expected).abs() <= 1e-9 * scale
    }

    #[test]
    fn fit_matches_the_r_prcomp_oracle() {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("testdata")
            .join("pca_r_oracle.tsv");
        let text = std::fs::read_to_string(&path)
            .unwrap_or_else(|error| panic!("read oracle {path:?}: {error}"));

        let mut checked = 0usize;
        for line in text.lines() {
            let line = line.trim_end();
            if line.is_empty() || line.starts_with('#') || line.starts_with("name\t") {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            assert_eq!(fields.len(), 7, "malformed oracle row: {line:?}");

            let name = fields[0];
            let n_regions: usize = fields[1].parse().expect("n_regions");
            let n_samples: usize = fields[2].parse().expect("n_samples");
            let values = parse_numbers(fields[3]);
            let expected_sdev = parse_numbers(fields[4]);
            let expected_variance = parse_numbers(fields[5]);
            let expected_scores = parse_numbers(fields[6]);

            assert_eq!(
                values.len(),
                n_regions * n_samples,
                "{name}: matrix size"
            );

            // R flattens `t(m)` column-major, so values arrive grouped by region.
            let regions: Vec<Vec<f64>> = values.chunks(n_samples).map(<[f64]>::to_vec).collect();
            let fitted = fit(&samples_from_columns(&regions, n_samples))
                .unwrap_or_else(|| panic!("{name}: fit() returned None"));

            assert_eq!(fitted.n_features, n_regions, "{name}: n_features");
            assert_eq!(fitted.n_samples, n_samples, "{name}: n_samples");
            assert_eq!(
                fitted.components.len(),
                expected_sdev.len(),
                "{name}: component count"
            );

            for (index, expected) in expected_sdev.iter().enumerate() {
                assert!(
                    oracle_close(fitted.sdev(index), *expected),
                    "{name}: sdev[{index}] got {}, R has {expected}",
                    fitted.sdev(index)
                );
            }
            for (index, expected) in expected_variance.iter().enumerate() {
                assert!(
                    oracle_close(fitted.variance_explained(index), *expected),
                    "{name}: variance_explained[{index}] got {}, R has {expected}",
                    fitted.variance_explained(index)
                );
            }

            // `t(p$x)` is components-by-samples, flattened sample-major.
            let n_components = expected_sdev.len();
            assert_eq!(
                expected_scores.len(),
                n_samples * n_components,
                "{name}: score count"
            );
            for component in 0..n_components {
                let actual = &fitted.components[component].scores;
                let expected: Vec<f64> = (0..n_samples)
                    .map(|sample| expected_scores[sample * n_components + component])
                    .collect();
                // The component's sign is arbitrary in R, so match it first; the
                // magnitudes are what must agree.
                let dot: f64 = actual.iter().zip(expected.iter()).map(|(a, b)| a * b).sum();
                let sign = if dot < 0.0 { -1.0 } else { 1.0 };
                for (sample, (actual, expected)) in actual.iter().zip(expected.iter()).enumerate() {
                    assert!(
                        oracle_close(actual * sign, *expected),
                        "{name}: score[PC{}][sample {sample}] got {} (sign {sign}), R has {expected}",
                        component + 1,
                        actual * sign
                    );
                }
            }

            checked += 1;
        }

        // Guard against a silently truncated oracle.
        assert!(checked >= 7, "oracle only had {checked} cases");
        eprintln!("validated {checked} PCA cases against the R oracle");
    }
}
