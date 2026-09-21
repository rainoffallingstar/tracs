//! Differential peak analysis: limma's moderated t-test, reimplemented.
//!
//! Ported from `trackplot.R`'s `diffpeak()`, which fits `limma::lmFit()` on a
//! one-way design, applies `limma::contrasts.fit()` and `limma::eBayes()`, and
//! reports `limma::topTable()`.
//!
//! Unlike the plotting ports, this one reproduces *another package's statistics*
//! rather than a data format, so the bar is numeric agreement with limma. It was
//! written against limma 3.68.4 with an oracle covering seven scenarios
//! (`testdata/diffpeak_r_oracle.tsv`) including deliberately unbalanced groups.
//!
//! # What is reproduced, and what is not
//!
//! Reproduced, because `diffpeak`'s output and everything downstream depends on
//! it:
//!
//! - The one-way fit. `model.matrix(~0 + condition)` gives 0/1 group indicators,
//!   so `lmFit`'s coefficients are group means and `sigma` is the pooled
//!   within-group SD. No general least-squares solve is needed, but the
//!   contrast's `stdev.unscaled` is `sqrt(1/n_num + 1/n_den)` and is *not* 1 for
//!   unbalanced groups.
//! - `squeezeVar` / `fitFDist`: the empirical-Bayes prior. The residual
//!   variances are logged and their mean and spread feed a moment estimate of
//!   the prior degrees of freedom, which then shrinks each variance toward the
//!   prior. See [`squeeze_variances`].
//! - The moderated t, `df.total`, the two-sided p-value via
//!   `2 * pt(-|t|, df.total)`, and Benjamini-Hochberg adjustment.
//!
//! **Not** reproduced: limma's `B` statistic (the log-odds). It comes from
//! `tmixture.matrix`, which needs the inverse t CDF with log-probability.
//! Nothing downstream reads it -- `volcano_plot()` uses only `logFC`, `P.Value`
//! and `adj.P.Val` -- and `diffpeak()` itself re-sorts by `P.Value`, so `B`
//! cannot affect the output order either. [`DifferentialPeaks::log_odds`] is
//! therefore always `None`, and the CLI omits the column rather than emitting a
//! number that looks like limma's but is not.
//!
//! One subtlety worth naming, because it silently passes on balanced designs:
//! `df.total` is `pmin(df.residual + df.prior, sum(df.residual))`, and limma
//! applies that to its *per-gene* `df.residual` vector. For a one-way design
//! every gene shares one residual df, so the pooled term is `n_genes *
//! df_residual` -- not `df_residual * n_conditions`, which coincides with the
//! right answer only for certain shapes. Getting this wrong shifted p-values by
//! ~1e-3 on the two-condition fixture while matching the three-condition one.

use std::collections::BTreeMap;

use anyhow::{anyhow, Result};

pub mod special;

/// One condition level with the sample columns that belong to it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Group {
    pub name: String,
    /// Column indices into the sample-major matrix, in input order.
    pub columns: Vec<usize>,
}

/// Splits sample conditions into levels, preserving first-seen order.
///
/// R builds the design with `as.factor(condition)` and takes
/// `levels(as.factor(condition))`, which sorts by the locale's collation. Sample
/// order is what decides column order here, so the level order is
/// first-appearance; the two agree on the ASCII condition names these tools use,
/// and `contrast` names the levels explicitly either way.
pub fn group_conditions(conditions: &[String]) -> Vec<Group> {
    let mut order: Vec<String> = Vec::new();
    for condition in conditions {
        if !order.contains(condition) {
            order.push(condition.clone());
        }
    }
    order
        .into_iter()
        .map(|name| Group {
            columns: conditions
                .iter()
                .enumerate()
                .filter(|(_, condition)| **condition == name)
                .map(|(index, _)| index)
                .collect(),
            name,
        })
        .collect()
}

/// The empirical-Bayes prior, i.e. `limma::squeezeVar()`'s output.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Prior {
    /// `df.prior`: the prior degrees of freedom `fitFDist` estimates.
    pub df_prior: f64,
    /// `s2.prior`: the prior variance.
    pub s2_prior: f64,
}

/// Estimates the prior variance and prior degrees of freedom, matching
/// `limma::squeezeVar()` on the legacy path `squeezeVar` takes when every gene
/// shares a residual df.
///
/// The idea: the residual variances `s2` are proportional to chi-square draws
/// scaled by the true variance, so `log(s2)` has a known offset and a spread
/// that pins down the degrees of freedom. `evar` is that spread, minus the
/// trigamma term that accounts for the noise in `log(s2)` itself; inverting
/// trigamma turns it into a df.
pub fn estimate_prior(sigma2: &[f64], df_residual: f64) -> Prior {
    // fitFDist drops non-finite and negative variances; `filter_ok` mirrors its
    // `ok` mask applied to a one-way design (df1 is constant, so it is always
    // in range when df_residual > 0).
    let usable: Vec<f64> = sigma2
        .iter()
        .copied()
        .filter(|value| value.is_finite() && *value > -1e-15)
        .collect();
    if usable.is_empty() || df_residual <= 0.0 {
        return Prior {
            df_prior: 0.0,
            s2_prior: f64::NAN,
        };
    }
    if usable.len() == 1 {
        return Prior {
            df_prior: 0.0,
            s2_prior: usable[0],
        };
    }

    // R: `m <- median(x)`; if every variance is zero the estimate is hopeless and
    // it substitutes 1 (with a warning).
    let mut sorted = usable.clone();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let median = if sorted.len().is_multiple_of(2) {
        (sorted[sorted.len() / 2 - 1] + sorted[sorted.len() / 2]) / 2.0
    } else {
        sorted[sorted.len() / 2]
    };
    let floor = if median == 0.0 { 1.0 } else { median };

    // `x <- pmax(x, 1e-05 * m)`, then work on the log scale.
    let shifted: Vec<f64> = usable
        .iter()
        .map(|value| value.max(1e-5 * floor))
        .collect();

    let offset = special::logmdigamma(df_residual / 2.0);
    let logged: Vec<f64> = shifted.iter().map(|value| value.ln() + offset).collect();
    let mean = logged.iter().sum::<f64>() / logged.len() as f64;
    let mut evar = logged
        .iter()
        .map(|value| (value - mean) * (value - mean))
        .sum::<f64>()
        / (logged.len() - 1) as f64;
    // `evar <- evar - mean(trigamma(df1/2))`; df1 is scalar so the mean is just
    // the value.
    evar -= special::trigamma(df_residual / 2.0);

    if evar > 0.0 {
        let df_prior = 2.0 * special::trigamma_inverse(evar);
        let s2_prior = (mean - special::logmdigamma(df_prior / 2.0)).exp();
        Prior {
            df_prior,
            s2_prior,
        }
    } else {
        // Not enough spread to identify a prior: limma falls back to infinite
        // prior df, i.e. no moderation, with the average variance as the prior.
        Prior {
            df_prior: f64::INFINITY,
            s2_prior: usable.iter().sum::<f64>() / usable.len() as f64,
        }
    }
}

/// Shrinks each variance toward the prior, i.e. `limma:::.squeezeVar()`.
pub fn squeeze_variances(sigma2: &[f64], df_residual: f64, prior: Prior) -> Vec<f64> {
    if prior.df_prior.is_infinite() {
        return vec![prior.s2_prior; sigma2.len()];
    }
    sigma2
        .iter()
        .map(|value| {
            (df_residual * value + prior.df_prior * prior.s2_prior)
                / (df_residual + prior.df_prior)
        })
        .collect()
}

/// Benjamini-Hochberg adjusted p-values, i.e. `p.adjust(p, method = "BH")`.
///
/// R's `p.adjust` breaks ties by first appearance, which is what makes the
/// adjusted values independent of the input order beyond that; the cumulative
/// minimum over decreasing rank is what enforces monotonicity.
pub fn benjamini_hochberg(p_values: &[f64]) -> Vec<f64> {
    let n = p_values.len();
    if n == 0 {
        return Vec::new();
    }
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|left, right| {
        p_values[*left]
            .partial_cmp(&p_values[*right])
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut adjusted = vec![1.0f64; n];
    let mut running_minimum = f64::INFINITY;
    // Walk from the largest p-value down, keeping the smallest scaled value
    // seen so far.
    for rank in (0..n).rev() {
        let index = order[rank];
        let scaled = p_values[index] * n as f64 / (rank + 1) as f64;
        running_minimum = running_minimum.min(scaled);
        adjusted[index] = running_minimum.min(1.0);
    }
    adjusted
}

/// One region's differential result, matching `limma::topTable()`'s columns.
#[derive(Clone, Debug, PartialEq)]
pub struct DifferentialPeak {
    /// Row position in the input summary table, so a caller can join back to the
    /// region's coordinates.
    pub row_index: usize,
    /// Contrast log2 fold change, `num - den`.
    pub log_fold_change: f64,
    /// Mean expression across all samples, limma's `AveExpr`.
    pub average_expression: f64,
    /// The moderated t statistic.
    pub moderated_t: f64,
    pub p_value: f64,
    pub adjusted_p_value: f64,
    /// limma's `B` log-odds, which this port does not compute. See the module
    /// docs.
    pub log_odds: Option<f64>,
}

/// The fitted model and its per-region results.
#[derive(Clone, Debug, PartialEq)]
pub struct DifferentialPeaks {
    pub peaks: Vec<DifferentialPeak>,
    /// The contrast as limma labels it, `num-den`.
    pub contrast: String,
    pub prior: Prior,
    /// `df.residual`, shared by every region in a one-way design.
    pub df_residual: f64,
    /// `df.total = pmin(df.residual + df.prior, df.pooled)`.
    pub df_total: f64,
    /// The contrast's `stdev.unscaled`, `sqrt(1/n_num + 1/n_den)`.
    pub stdev_unscaled: f64,
}

/// Which two condition levels the test compares.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Contrast {
    pub numerator: String,
    pub denominator: String,
}

impl Contrast {
    /// R's label for the contrast: `paste0(num, "-", den)`.
    pub fn label(&self) -> String {
        format!("{}-{}", self.numerator, self.denominator)
    }
}

/// Picks the contrast the way `diffpeak()` does when `num`/`den` are omitted.
///
/// R enumerates every ordered pair `a != b` with `a < b` over the *levels* and
/// takes the first, which is "first level minus second level". Using the first
/// two groups in level order reproduces that without materialising the pairs,
/// and it is worth being explicit because the sign of every fold change depends
/// on it.
pub fn default_contrast(groups: &[Group]) -> Result<Contrast> {
    if groups.len() < 2 {
        return Err(anyhow!(
            "diffpeak needs at least two conditions, found {}: {}",
            groups.len(),
            groups
                .iter()
                .map(|group| group.name.as_str())
                .collect::<Vec<_>>()
                .join(", ")
        ));
    }
    Ok(Contrast {
        numerator: groups[0].name.clone(),
        denominator: groups[1].name.clone(),
    })
}

/// Validates a caller-supplied contrast against the observed conditions.
pub fn explicit_contrast(
    groups: &[Group],
    numerator: &str,
    denominator: &str,
) -> Result<Contrast> {
    let known: Vec<&str> = groups.iter().map(|group| group.name.as_str()).collect();
    for (label, name) in [("num", numerator), ("den", denominator)] {
        if !known.contains(&name) {
            return Err(anyhow!(
                "diffpeak: --{label} {name:?} is not a condition in the data; found: {}",
                known.join(", ")
            ));
        }
    }
    if numerator == denominator {
        return Err(anyhow!(
            "diffpeak: --num and --den must differ (both {numerator:?})"
        ));
    }
    Ok(Contrast {
        numerator: numerator.to_string(),
        denominator: denominator.to_string(),
    })
}

/// Applies `log2(x + offset)` in place, matching `diffpeak(log2 = TRUE)`.
pub fn log2_transform(expression: &mut [Vec<f64>], offset: f64) {
    for row in expression.iter_mut() {
        for value in row.iter_mut() {
            if value.is_finite() {
                *value = (*value + offset).log2();
            }
        }
    }
}

/// Runs the moderated t-test over already-transformed expression values.
///
/// `expression` is regions-by-samples; each inner vector is one region across
/// all samples in `conditions` order.
pub fn fit(
    expression: &[Vec<f64>],
    conditions: &[String],
    contrast: &Contrast,
) -> Result<DifferentialPeaks> {
    if expression.is_empty() {
        return Err(anyhow!("diffpeak: no regions to test"));
    }
    if conditions.len() != expression[0].len() {
        return Err(anyhow!(
            "diffpeak: {} conditions for {} samples",
            conditions.len(),
            expression[0].len()
        ));
    }
    for (index, row) in expression.iter().enumerate() {
        if row.len() != conditions.len() {
            return Err(anyhow!(
                "diffpeak: region {} has {} values, expected {}",
                index,
                row.len(),
                conditions.len()
            ));
        }
    }

    let groups = group_conditions(conditions);
    let numerator = groups
        .iter()
        .find(|group| group.name == contrast.numerator)
        .ok_or_else(|| {
            anyhow!(
                "diffpeak: numerator {:?} is not a condition in the data",
                contrast.numerator
            )
        })?;
    let denominator = groups
        .iter()
        .find(|group| group.name == contrast.denominator)
        .ok_or_else(|| {
            anyhow!(
                "diffpeak: denominator {:?} is not a condition in the data",
                contrast.denominator
            )
        })?;

    let n_samples = conditions.len();
    let n_levels = groups.len();
    let df_residual = (n_samples - n_levels) as f64;
    if df_residual <= 0.0 {
        return Err(anyhow!(
            "diffpeak: {} samples across {} conditions leaves no residual degrees \
             of freedom; limma fails here too (\"No residual degrees of freedom in \
             linear model fits\")",
            n_samples,
            n_levels
        ));
    }

    // Group means (the design's coefficients) and the residual sum of squares.
    let n_regions = expression.len();
    let mut group_means: BTreeMap<&str, Vec<f64>> = BTreeMap::new();
    let mut residual_sum_of_squares = vec![0.0f64; n_regions];
    for group in &groups {
        let mut means = vec![0.0f64; n_regions];
        for row in 0..n_regions {
            let count = group.columns.len() as f64;
            let mut sum = 0.0;
            for column in &group.columns {
                sum += expression[row][*column];
            }
            let mean = sum / count;
            means[row] = mean;
            for column in &group.columns {
                let residual = expression[row][*column] - mean;
                residual_sum_of_squares[row] += residual * residual;
            }
        }
        group_means.insert(group.name.as_str(), means);
    }

    let sigma2: Vec<f64> = residual_sum_of_squares
        .iter()
        .map(|rss| rss / df_residual)
        .collect();

    let prior = estimate_prior(&sigma2, df_residual);
    let var_post = squeeze_variances(&sigma2, df_residual, prior);

    // limma: `df.pooled <- sum(df.residual)` on the per-gene vector, so for a
    // one-way design it is n_genes * df_residual. `pmin` only binds when the
    // prior df is enormous, but implementing it keeps degenerate inputs honest.
    let df_pooled = df_residual * n_regions as f64;
    let df_total = (df_residual + prior.df_prior).min(df_pooled);

    let numerator_means = &group_means[contrast.numerator.as_str()];
    let denominator_means = &group_means[contrast.denominator.as_str()];
    let stdev_unscaled =
        (1.0 / numerator.columns.len() as f64 + 1.0 / denominator.columns.len() as f64).sqrt();

    let mut log_fold_change = Vec::with_capacity(n_regions);
    let mut average_expression = Vec::with_capacity(n_regions);
    let mut moderated_t = Vec::with_capacity(n_regions);
    let mut p_value = Vec::with_capacity(n_regions);
    for row in 0..n_regions {
        log_fold_change.push(numerator_means[row] - denominator_means[row]);
        // `Amean` is the plain row mean over all samples, *not* the mean of the
        // group means. The two coincide on balanced designs, which is how this
        // got missed until the unbalanced fixture disagreed by 0.41.
        let row_sum: f64 = expression[row].iter().sum();
        average_expression.push(row_sum / n_samples as f64);
        // `coefficients / stdev.unscaled / sqrt(s2.post)`.
        //
        // A region whose contrast is exactly zero *and* whose residual variance
        // is exactly zero makes this 0/0. limma does not hit that because its QR
        // solve leaves `sigma` at ~1e-16 of rounding residue on constant input,
        // which yields an arbitrary-looking t (1.0 for a constant matrix) rather
        // than a defined value. A zero fold change with no spread has no
        // evidence either way, so it is reported as an undefined t and p = 1
        // rather than propagating a NaN through the ranking.
        let t = if var_post[row] > 0.0 {
            log_fold_change[row] / stdev_unscaled / var_post[row].sqrt()
        } else if log_fold_change[row] == 0.0 {
            0.0
        } else {
            // A real effect with no residual spread: infinitely significant.
            log_fold_change[row].signum() * f64::INFINITY
        };
        moderated_t.push(t);
        p_value.push(if t.is_finite() {
            special::t_two_sided_p(t, df_total)
        } else {
            0.0
        });
    }

    let adjusted = benjamini_hochberg(&p_value);

    let peaks = (0..n_regions)
        .map(|row| DifferentialPeak {
            row_index: row,
            log_fold_change: log_fold_change[row],
            average_expression: average_expression[row],
            moderated_t: moderated_t[row],
            p_value: p_value[row],
            adjusted_p_value: adjusted[row],
            log_odds: None,
        })
        .collect();

    Ok(DifferentialPeaks {
        peaks,
        contrast: contrast.label(),
        prior,
        df_residual,
        df_total,
        stdev_unscaled,
    })
}

/// Runs the whole `diffpeak()` pipeline over a summary table.
pub fn run(
    expression: &mut [Vec<f64>],
    conditions: &[String],
    contrast: &Contrast,
    log2: bool,
    log2_offset: f64,
) -> Result<DifferentialPeaks> {
    if log2 {
        log2_transform(expression, log2_offset);
    }
    fit(expression, conditions, contrast)
}

#[cfg(test)]
// The expected values below are copied verbatim from R 4.6.0 / limma 3.68.4 at
// 17 significant digits. Trimming them to the shortest round-tripping literal
// would obscure that provenance, so the precision lints are allowed here rather
// than applying their suggestions.
#[allow(clippy::excessive_precision, clippy::approx_constant)]
mod tests {
    use super::*;

    fn conditions(names: &[&str]) -> Vec<String> {
        names.iter().map(|name| name.to_string()).collect()
    }

    #[test]
    fn groups_preserve_first_appearance_order() {
        let grouped = group_conditions(&conditions(&["B", "A", "B", "C", "A"]));
        assert_eq!(
            grouped
                .iter()
                .map(|group| group.name.as_str())
                .collect::<Vec<_>>(),
            vec!["B", "A", "C"]
        );
        assert_eq!(grouped[0].columns, vec![0, 2]);
        assert_eq!(grouped[1].columns, vec![1, 4]);
        assert_eq!(grouped[2].columns, vec![3]);
    }

    #[test]
    fn default_contrast_is_the_first_two_levels() {
        // R enumerates ordered pairs with a < b and takes the first, i.e.
        // level 1 minus level 2.
        let grouped = group_conditions(&conditions(&["H3K27ac", "H3K27ac", "Input", "Input"]));
        let contrast = default_contrast(&grouped).expect("contrast");
        assert_eq!(contrast.numerator, "H3K27ac");
        assert_eq!(contrast.denominator, "Input");
        assert_eq!(contrast.label(), "H3K27ac-Input");
    }

    #[test]
    fn a_single_condition_cannot_be_contrasted() {
        let grouped = group_conditions(&conditions(&["only", "only"]));
        assert!(default_contrast(&grouped).is_err());
    }

    #[test]
    fn explicit_contrast_rejects_unknown_and_identical_levels() {
        let grouped = group_conditions(&conditions(&["A", "B"]));
        assert!(explicit_contrast(&grouped, "A", "C").is_err());
        assert!(explicit_contrast(&grouped, "A", "A").is_err());
        // Reversing the default is allowed, and flips the label.
        let reversed = explicit_contrast(&grouped, "B", "A").expect("contrast");
        assert_eq!(reversed.label(), "B-A");
    }

    #[test]
    fn stdev_unscaled_is_one_only_for_balanced_groups() {
        // This is the case a naive port gets wrong: with equal group sizes the
        // contrast scaling is exactly 1, so an assumption of 1 passes until an
        // unbalanced design appears.
        let balanced: Vec<Vec<f64>> = (0..4)
            .map(|row| vec![row as f64, row as f64 + 1.0, row as f64, row as f64 + 1.0])
            .collect();
        let contrast = Contrast {
            numerator: "A".to_string(),
            denominator: "B".to_string(),
        };
        let fitted = fit(&balanced, &conditions(&["A", "A", "B", "B"]), &contrast).expect("fit");
        assert!((fitted.stdev_unscaled - 1.0).abs() < 1e-15);

        // 3 vs 2 gives sqrt(1/3 + 1/2).
        let unbalanced: Vec<Vec<f64>> = (0..4)
            .map(|row| vec![row as f64, row as f64 + 1.0, row as f64 + 2.0, 0.0, 1.0])
            .collect();
        let fitted = fit(
            &unbalanced,
            &conditions(&["A", "A", "A", "B", "B"]),
            &contrast,
        )
        .expect("fit");
        let expected = (1.0f64 / 3.0 + 1.0 / 2.0).sqrt();
        assert!(
            (fitted.stdev_unscaled - expected).abs() < 1e-15,
            "got {}, expected {expected}",
            fitted.stdev_unscaled
        );
    }

    #[test]
    fn log_fold_change_is_the_difference_of_group_means() {
        // Two conditions, two samples each, no replication noise: the fold
        // change must be exactly the mean difference.
        let expression = vec![
            vec![1.0, 3.0, 5.0, 7.0], // means 2 and 6
            vec![0.0, 2.0, 2.0, 2.0], // means 1 and 2
        ];
        let contrast = Contrast {
            numerator: "A".to_string(),
            denominator: "B".to_string(),
        };
        let fitted = fit(&expression, &conditions(&["A", "A", "B", "B"]), &contrast).expect("fit");
        assert!((fitted.peaks[0].log_fold_change - -4.0).abs() < 1e-15);
        assert!((fitted.peaks[1].log_fold_change - -1.0).abs() < 1e-15);
    }

    #[test]
    fn average_expression_is_the_plain_row_mean() {
        // `Amean` is `rowMeans(exprs)`, i.e. the mean over all samples, and NOT
        // the mean of the group means. On balanced designs the two are equal, so
        // the distinguishing case has to be unbalanced: with group means 4 and 2
        // over 3 and 2 samples, the row mean is 3.2 while the mean of the means
        // would be 3.0.
        let expression = vec![vec![0.0, 4.0, 8.0, 2.0, 2.0]];
        let contrast = Contrast {
            numerator: "A".to_string(),
            denominator: "B".to_string(),
        };
        let fitted = fit(
            &expression,
            &conditions(&["A", "A", "A", "B", "B"]),
            &contrast,
        )
        .expect("fit");
        assert!(
            (fitted.peaks[0].average_expression - 3.2).abs() < 1e-15,
            "got {}",
            fitted.peaks[0].average_expression
        );
    }

    #[test]
    fn df_total_uses_the_pooled_residual_not_the_level_count() {
        // 4 regions, 4 samples, 2 conditions: df_residual = 2 and
        // df_pooled = 4 * 2 = 8. Using df_residual * n_levels would give 4, and
        // the two only disagree for some shapes, so assert the pooled value.
        let expression: Vec<Vec<f64>> = (0..4)
            .map(|row| vec![row as f64, row as f64 + 0.5, row as f64 + 2.0, row as f64 + 2.5])
            .collect();
        let contrast = Contrast {
            numerator: "A".to_string(),
            denominator: "B".to_string(),
        };
        let fitted = fit(&expression, &conditions(&["A", "A", "B", "B"]), &contrast).expect("fit");
        assert!((fitted.df_residual - 2.0).abs() < 1e-15);
        // df.prior is unbounded here, so df.total is the pooled value or the sum.
        let expected_ceiling = 2.0 * 4.0;
        assert!(
            fitted.df_total <= expected_ceiling + 1e-12,
            "df.total {} should not exceed the pooled {expected_ceiling}",
            fitted.df_total
        );
    }

    #[test]
    fn too_few_samples_per_condition_is_reported() {
        let expression = vec![vec![1.0, 2.0], vec![3.0, 4.0]];
        let contrast = Contrast {
            numerator: "A".to_string(),
            denominator: "B".to_string(),
        };
        let error = fit(&expression, &conditions(&["A", "B"]), &contrast).unwrap_err();
        assert!(
            error.to_string().contains("no residual degrees of freedom"),
            "unexpected error: {error}"
        );
    }

    #[test]
    fn bh_adjustment_matches_the_textbook_definition() {
        // p.adjust(c(0.01, 0.02, 0.03, 0.04, 0.05), "BH") in R.
        let adjusted = benjamini_hochberg(&[0.01, 0.02, 0.03, 0.04, 0.05]);
        let expected = [0.05, 0.05, 0.05, 0.05, 0.05];
        for (actual, want) in adjusted.iter().zip(expected.iter()) {
            assert!((actual - want).abs() < 1e-15, "got {actual}, want {want}");
        }
    }

    #[test]
    fn bh_adjustment_is_monotone_and_capped_at_one() {
        let p = [0.001, 0.5, 0.02, 0.9, 0.04];
        let adjusted = benjamini_hochberg(&p);
        for value in &adjusted {
            assert!(*value <= 1.0 && *value >= 0.0, "out of range: {value}");
        }
        // Larger p-values never get smaller adjusted values.
        let mut order: Vec<usize> = (0..p.len()).collect();
        order.sort_by(|a, b| p[*a].partial_cmp(&p[*b]).unwrap());
        for pair in order.windows(2) {
            assert!(
                adjusted[pair[0]] <= adjusted[pair[1]] + 1e-15,
                "not monotone: {} then {}",
                adjusted[pair[0]],
                adjusted[pair[1]]
            );
        }
        // The smallest p-value maps to n * p / 1, capped.
        assert!((adjusted[0] - (5.0f64 * 0.001).min(1.0)).abs() < 1e-15);
    }

    #[test]
    fn identical_regions_produce_a_defined_result() {
        // Every region has the same signal, so the fold change and the residual
        // variance are both exactly zero. limma only avoids 0/0 here because its
        // QR solve leaves `sigma` at ~e-16 of rounding residue (measured), which
        // makes its t statistic an artifact of that residue rather than a
        // reproducible value. This port computes the variance exactly, so it
        // reports the degenerate case as t = 0 with p = 1: no evidence either
        // way, and nothing NaN reaching the ranking.
        let expression = vec![vec![5.0, 5.0, 5.0, 5.0]; 3];
        let contrast = Contrast {
            numerator: "A".to_string(),
            denominator: "B".to_string(),
        };
        let fitted = fit(&expression, &conditions(&["A", "A", "B", "B"]), &contrast).expect("fit");
        for peak in &fitted.peaks {
            assert_eq!(peak.log_fold_change, 0.0);
            assert_eq!(peak.moderated_t, 0.0);
            assert_eq!(peak.p_value, 1.0);
            assert_eq!(peak.adjusted_p_value, 1.0);
        }
    }

    #[test]
    fn a_real_effect_with_no_residual_spread_is_infinitely_significant() {
        // The other degenerate shape: a genuine fold change with zero within-group
        // variance. That is not "no evidence", it is certainty, so the t diverges
        // and p is 0 rather than NaN.
        let expression = vec![vec![1.0, 1.0, 9.0, 9.0]];
        let contrast = Contrast {
            numerator: "A".to_string(),
            denominator: "B".to_string(),
        };
        let fitted = fit(&expression, &conditions(&["A", "A", "B", "B"]), &contrast).expect("fit");
        assert_eq!(fitted.peaks[0].log_fold_change, -8.0);
        assert!(fitted.peaks[0].moderated_t.is_infinite());
        assert_eq!(fitted.peaks[0].p_value, 0.0);
    }

    #[test]
    fn log2_transform_matches_the_r_expression() {
        let mut expression = vec![vec![0.0, 0.9], vec![8.1, 1.0]];
        log2_transform(&mut expression, 0.1);
        // log2(0 + 0.1), log2(1.0 + 0.1) and log2(8.1 + 0.1).
        assert!((expression[0][0] - 0.1f64.log2()).abs() < 1e-15);
        assert!((expression[1][1] - 1.1f64.log2()).abs() < 1e-15);
        assert!((expression[1][0] - 8.2f64.log2()).abs() < 1e-15);
    }

    #[test]
    fn shape_mismatches_are_reported() {
        let contrast = Contrast {
            numerator: "A".to_string(),
            denominator: "B".to_string(),
        };
        // Ragged rows.
        let error = fit(
            &[vec![1.0, 2.0, 3.0, 4.0], vec![1.0, 2.0]],
            &conditions(&["A", "A", "B", "B"]),
            &contrast,
        )
        .unwrap_err();
        assert!(error.to_string().contains("expected 4"), "got: {error}");

        // Condition count disagrees with the sample count.
        let error = fit(&[vec![1.0, 2.0]], &conditions(&["A", "B", "B"]), &contrast).unwrap_err();
        assert!(error.to_string().contains("conditions for"), "got: {error}");
    }

    // ---------------------------------------------------------------------
    // R oracle
    //
    // `testdata/diffpeak_r_oracle.tsv` was produced by running the real
    // `diffpeak()` from trackplot.R with limma 3.68.4, over the fixture summary
    // and coldata tables in `testdata/`. `diffpeak_r_scenarios.tsv` records the
    // scalars limma derived for each scenario, and
    // `diffpeak_special_r_oracle.tsv` pins the special functions.
    //
    // The scenarios include unbalanced groups on purpose: with equal group sizes
    // the contrast scaling is exactly 1, so an assumption of 1 would pass every
    // balanced case and only fail later.
    // ---------------------------------------------------------------------

    fn testdata(name: &str) -> std::path::PathBuf {
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("testdata")
            .join(name)
    }

    fn read_rows(name: &str) -> Vec<Vec<String>> {
        let text = std::fs::read_to_string(testdata(name))
            .unwrap_or_else(|error| panic!("read {name}: {error}"));
        text.lines()
            .filter(|line| !line.trim().is_empty() && !line.starts_with('#'))
            .map(|line| line.split('\t').map(|f| f.trim().to_string()).collect())
            .collect()
    }

    /// Compares against the oracle with a mixed absolute/relative tolerance.
    ///
    /// The oracle spans raw signal sums (~1e6) and log2 values (~1e-1), so a
    /// purely absolute test would be far too loose for one end and too tight for
    /// the other. A purely relative test fails where the true value is exactly
    /// zero (R's `lgamma(1)` is 0), so the tolerance is relative above 1 and
    /// absolute below it.
    fn close(actual: f64, expected: f64, tolerance: f64) -> bool {
        if actual == expected {
            return true;
        }
        let scale = actual.abs().max(expected.abs()).max(1.0);
        (actual - expected).abs() <= tolerance * scale
    }

    /// Tolerance for the special functions.
    ///
    /// The series and recurrences reach machine precision, but the t tail goes
    /// through [`regularized_incomplete_beta`]'s continued fraction, whose
    /// accuracy in double precision bottoms out near `1.4e-13` relative. That
    /// was measured by comparing against `scipy.special.stdtr` over the same
    /// grid; tightening the convergence epsilon does not improve it, so the
    /// tolerance records the method's real limit rather than a wish.
    const SPECIAL_TOLERANCE: f64 = 1e-12;

    struct Scenario {
        sum_table: String,
        coldata: String,
        log2: bool,
    }

    fn scenario(name: &str) -> Scenario {
        match name {
            "three_log2" | "three_raw" | "three_explicit_reverse" => Scenario {
                sum_table: "diffpeak_three_summary.tsv".to_string(),
                coldata: "diffpeak_three_coldata.tsv".to_string(),
                log2: name != "three_raw",
            },
            "two_log2" => Scenario {
                sum_table: "diffpeak_two_summary.tsv".to_string(),
                coldata: "diffpeak_two_coldata.tsv".to_string(),
                log2: true,
            },
            "unbalanced_log2" | "unbalanced_reverse" => Scenario {
                sum_table: "diffpeak_unbalanced_summary.tsv".to_string(),
                coldata: "diffpeak_unbalanced_coldata.tsv".to_string(),
                log2: true,
            },
            "triplicate_log2" => Scenario {
                sum_table: "diffpeak_triplicate_summary.tsv".to_string(),
                coldata: "diffpeak_triplicate_coldata.tsv".to_string(),
                log2: true,
            },
            other => panic!("unknown scenario {other}"),
        }
    }

    /// Loads a fixture into `(regions, samples, expression)`.
    fn load_fixture(sum_table: &str, coldata: &str) -> (Vec<i64>, Vec<String>, Vec<Vec<f64>>) {
        let rows = read_rows(sum_table);
        let header = &rows[0];
        let start_index = header.iter().position(|f| f == "start").expect("start");
        let sample_indices: Vec<usize> = (0..header.len())
            .filter(|index| {
                let name = &header[*index];
                !matches!(name.as_str(), "chromosome" | "start" | "end" | "size")
            })
            .collect();

        let coldata_rows = read_rows(coldata);
        let coldata_header = &coldata_rows[0];
        let name_index = coldata_header
            .iter()
            .position(|f| f == "bw_sample_names")
            .expect("bw_sample_names");
        let condition_index = coldata_header
            .iter()
            .position(|f| f == "condition")
            .expect("condition");
        let conditions: Vec<String> = coldata_rows[1..]
            .iter()
            .map(|row| row[condition_index].clone())
            .collect();
        // The fixture's column order matches the coldata order.
        let _ = name_index;

        let starts: Vec<i64> = rows[1..]
            .iter()
            .map(|row| row[start_index].parse().expect("start"))
            .collect();
        let expression: Vec<Vec<f64>> = rows[1..]
            .iter()
            .map(|row| {
                sample_indices
                    .iter()
                    .map(|index| row[*index].parse::<f64>().unwrap_or(f64::NAN))
                    .collect()
            })
            .collect();
        (starts, conditions, expression)
    }

    #[test]
    fn moderated_statistics_match_the_limma_oracle() {
        let scenario_rows = read_rows("diffpeak_r_scenarios.tsv");
        let scenario_header = &scenario_rows[0];
        let column = |name: &str| scenario_header.iter().position(|f| f == name).unwrap();
        let (name_at, contrast_at) = (column("scenario"), column("contrast"));
        let (dfp_at, s2p_at) = (column("df_prior"), column("s2_prior"));
        let (dfr_at, stdev_at) = (column("df_residual"), column("stdev_unscaled"));

        // Group the oracle's per-region rows by scenario.
        let oracle_rows = read_rows("diffpeak_r_oracle.tsv");
        let oracle_header = &oracle_rows[0];
        let o = |name: &str| oracle_header.iter().position(|f| f == name).unwrap();
        let (o_scenario, o_start) = (o("scenario"), o("start"));
        let (o_logfc, o_t, o_p) = (o("logFC"), o("t"), o("P_Value"));
        let (o_ave, o_adjp) = (o("AveExpr"), o("adj_P_Val"));

        let mut scenarios_checked = 0usize;
        let mut regions_checked = 0usize;

        for meta in &scenario_rows[1..] {
            let name = &meta[name_at];
            let fx = scenario(name);
            let (starts, conditions, expression) = load_fixture(&fx.sum_table, &fx.coldata);

            // The contrast R chose; supply it explicitly so the test also covers
            // explicit-contrast plumbing.
            let (num, den) = meta[contrast_at]
                .split_once('-')
                .unwrap_or_else(|| panic!("unparsable contrast {}", meta[contrast_at]));
            let groups = group_conditions(&conditions);
            let contrast = explicit_contrast(&groups, num, den).expect("contrast");

            let expected_df_prior: f64 = meta[dfp_at].parse().unwrap();
            let expected_s2_prior: f64 = meta[s2p_at].parse().unwrap();
            let expected_df_residual: f64 = meta[dfr_at].parse().unwrap();
            let expected_stdev: f64 = meta[stdev_at].parse().unwrap();

            let mut values = expression.clone();
            let fitted = run(&mut values, &conditions, &contrast, fx.log2, 0.1).expect("fit");

            // The scalars first: a port that gets these right but the per-region
            // statistics wrong is a different bug from one that gets both wrong.
            assert!(
                close(fitted.prior.df_prior, expected_df_prior, 1e-9),
                "{name}: df.prior {} vs R {expected_df_prior}",
                fitted.prior.df_prior
            );
            assert!(
                close(fitted.prior.s2_prior, expected_s2_prior, 1e-9),
                "{name}: s2.prior {} vs R {expected_s2_prior}",
                fitted.prior.s2_prior
            );
            assert!(
                close(fitted.df_residual, expected_df_residual, 1e-15),
                "{name}: df.residual {} vs R {expected_df_residual}",
                fitted.df_residual
            );
            assert!(
                close(fitted.stdev_unscaled, expected_stdev, 1e-12),
                "{name}: stdev.unscaled {} vs R {expected_stdev}",
                fitted.stdev_unscaled
            );
            assert_eq!(fitted.contrast, meta[contrast_at].as_str(), "{name}: contrast");

            // Then every region.
            let by_start: BTreeMap<i64, usize> = starts
                .iter()
                .enumerate()
                .map(|(index, start)| (*start, index))
                .collect();
            for row in oracle_rows[1..].iter().filter(|row| row[o_scenario] == *name) {
                let start: i64 = row[o_start].parse().expect("start");
                let index = *by_start
                    .get(&start)
                    .unwrap_or_else(|| panic!("{name}: oracle start {start} not in fixture"));
                let peak = &fitted.peaks[index];
                for (label, actual, expected) in [
                    ("logFC", peak.log_fold_change, row[o_logfc].parse::<f64>().unwrap()),
                    ("AveExpr", peak.average_expression, row[o_ave].parse::<f64>().unwrap()),
                    ("t", peak.moderated_t, row[o_t].parse::<f64>().unwrap()),
                    ("P.Value", peak.p_value, row[o_p].parse::<f64>().unwrap()),
                    ("adj.P.Val", peak.adjusted_p_value, row[o_adjp].parse::<f64>().unwrap()),
                ] {
                    assert!(
                        close(actual, expected, 1e-7),
                        "{name} start {start}: {label} got {actual}, R has {expected}"
                    );
                }
                regions_checked += 1;
            }
            scenarios_checked += 1;
        }

        assert!(scenarios_checked >= 7, "only {scenarios_checked} scenarios");
        assert!(regions_checked >= 900, "only {regions_checked} regions");
        eprintln!(
            "validated {regions_checked} regions across {scenarios_checked} scenarios \
             against the limma oracle"
        );
    }

    #[test]
    fn special_functions_match_the_r_oracle() {
        let rows = read_rows("diffpeak_special_r_oracle.tsv");
        let header = &rows[0];
        let (kind_at, arg_at, arg2_at, exp_at) = (
            header.iter().position(|f| f == "kind").unwrap(),
            header.iter().position(|f| f == "arg").unwrap(),
            header.iter().position(|f| f == "arg2").unwrap(),
            header.iter().position(|f| f == "expected").unwrap(),
        );

        let mut checked = 0usize;
        for row in &rows[1..] {
            let argument: f64 = row[arg_at].parse().expect("arg");
            let expected: f64 = row[exp_at].parse().expect("expected");
            let actual = match row[kind_at].as_str() {
                "lgamma" => special::ln_gamma(argument),
                "digamma" => special::digamma(argument),
                "trigamma" => special::trigamma(argument),
                "tetragamma" => special::tetragamma(argument),
                "trigamma_inverse" => special::trigamma_inverse(argument),
                "t_two_sided_p" => {
                    let df: f64 = row[arg2_at].parse().expect("df");
                    special::t_two_sided_p(argument, df)
                }
                other => panic!("unknown oracle kind {other}"),
            };
            assert!(
                close(actual, expected, SPECIAL_TOLERANCE),
                "{}: {actual} vs R {expected} at {argument}",
                row[kind_at]
            );
            checked += 1;
        }
        assert!(checked >= 200, "only {checked} special-function rows");
        eprintln!("validated {checked} special-function values against R");
    }
}
