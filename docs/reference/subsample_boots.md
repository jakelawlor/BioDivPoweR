# Subsample bootstrapped communities to assess power for detecting change across community coverage

Downsample bootstrapped communities to lower sample sizes representing a
gradient of coverage values, and assess the proportion of
correct-direction detections within each coverage interval.

## Usage

``` r
subsample_boots(
  boots,
  pilot,
  method = "single",
  power = c(80),
  target_eff_size = NULL,
  cost_per_sample = NULL,
  seed = NULL,
  analysis_type = "sign",
  effect_minimum = NULL
)
```

## Arguments

- boots:

  tibble of bootstrapped communities across multiple effect sizes.
  Output from \`bootstrap_pilot()\`.

- pilot:

  original pilot study. Identical to the pilot input in
  \`bootstrap_pilot()\`

- method:

  "single" or "two" for single-treatment power analysis (one community
  over time), or two-treatment power analysis (comparing differences
  between two treatment communities)

- power:

  statistical power with which the user wants to detect richness change
  (defaults to 80).

- target_eff_size:

  if supplied, the level of richness difference the user wants to detect
  (in log2 ratio units). Adding a target value will change outputs of
  the function, specifying sampling needs to reach the richness
  difference target.

- cost_per_sample:

  cost per unit sample (linear only), which will add an axis to output
  plots for total cost per unit power.

- seed:

  random seed. Defaults to 1 so repeat runs will be identical, but since
  simulations rely on random draws, changing the seed will result in
  different answers.

- analysis_type:

  "sign" or "minimum". Describes the type of power analysis to conduct.
  "sign" identifies the number of samples needed to detect an effect in
  the same direction as the true effect, "minimum" identifies the number
  of samples needed to detect an effect both in the same direction and
  above user-supplied biologically-relevant effect.

- effect_minimum:

  if \`analysis_type\` is "minimum", a user-defined minimum effect size
  to form a minimum-effect null distribtuion.

## Value

List of outputs: (1) Plot of convergence of detections across community
coverage in multiple effect size bins. (2) Plot of mean detection across
community coverage within multiple effect size bins. (3) Plot of
proportion of correct detections across community coverage within
multiple effect size bins. (4) Plot of functional relationship between
sample size and power to detect richness differences. (5) Tibble of
input and output values.

## Examples

``` r
if(FALSE){subsample_boots(boots, pilot)}
```
