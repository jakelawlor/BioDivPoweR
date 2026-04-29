# Find coverage of ecological community surveys at 1-k samples

Find coverage of ecological community surveys at 1-k samples

## Usage

``` r
find_coverage(occupancy, k = 1000)
```

## Arguments

- occupancy:

  vector of average occupancy of all species in a community (\>0)

- k:

  max number of samples to test. Coverage should asymptote at a large
  number of samples.

## Value

vector of coverage achieved in k samples of an ecological community
