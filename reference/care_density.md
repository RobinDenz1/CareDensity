# Calculate the Care Density for all Patients

This function calculates the classic Care Density Index as defined by
Pollack et al. (2013) for each patient in the supplied dataset. Works
well with large patient-sharing networks.

## Usage

``` r
care_density(data, pat_col=1, data_frame=TRUE)
```

## Arguments

- data:

  A `data.frame` like object containing exactly two columns. One should
  include only patient IDs and the other one only provider IDs. Each row
  should denote one patient-provider contact. Multiple contacts (same
  rows) are allowed but have no effect on the outcome. Both patient and
  provider IDs should be unique, which means that one ID may not be in
  both rows.

- pat_col:

  Specifies which column of `data` includes the patient IDs. If the
  first column contains the patient IDs this should be kept at 1, if the
  second column contains the patient IDs it should be set to 2.

- data_frame:

  Set this argument to `TRUE` to return a `data.frame` instead of the
  `data.table` format that is used under the hood.

## Details

The Care Density (\\C_p\\) is "a patient-level measure that quantifies
the amount of patient-sharing among his or her providers" (DuGoff et al.
2018). Higher care densities have been posited to reflect greater
connections among a patients "care team". Formally, it is defined as:

\$\$C_p = \frac{\sum\_{i = 1}^m w\_{p, j}}{n_p (n_p - 1)/2}\$\$

with \\n_p\\ being the number of providers a patient has visited, \\m\\
defined as the number of all possible combinations of length two and
\\w\_{p, j}\\ being the number of patients that a pair of provider is
sharing. An example is given below and explained more thoroughly in the
vignette of this package.

Under the hood, this function uses the `igraph` package to construct a
patient-sharing network from the provided `data` to calculate the
weights. It then uses the `data.table` package to efficiently calculate
the care densities from a resulting edge list with weights.

## Value

Returns a single `data.frame` (or `data.table`) containing the sum of
all weights (`"sum_weights"`), the number of providers seen by each
patient (`"n"`) and the calculated Care Density (`"care_density"`).

## References

Pollack, Craig Evan, Gary E. Weissman, Klaus W. Lemke, Peter S. Hussey,
and Jonathan P. Weiner. (2013). "Patient Sharing Among Physicians and
Costs of Care: A Network Analytic Approach to Care Coordination Using
Claims Data". Journal of General Internal Medicine 28 (3), pp. 459-465.

DuGoff, Eva H., Sara Fernandes-Taylor, Gary E. Weissman, Joseph H.
Huntley, and Craig Evan Pollack. (2018). "A Scoping Review of
Patient-Sharing Network Studies Using Administrative Data".
Translational Behavioral Medicine 8 (4), pp. 598-625.

## Author

Robin Denz

## See also

[`fragmented_care_density`](https://robindenz1.github.io/CareDensity/reference/fragmented_care_density.md)

## Examples

``` r
library(CareDensity)
library(data.table)
#> 
#> Attaching package: ‘data.table’
#> The following object is masked from ‘package:base’:
#> 
#>     %notin%
library(igraph)
#> 
#> Attaching package: ‘igraph’
#> The following objects are masked from ‘package:stats’:
#> 
#>     decompose, spectrum
#> The following object is masked from ‘package:base’:
#> 
#>     union

# some arbitrary patient-provider contact data
data <- data.frame(PatID=c("1", "1", "1", "2", "2", "3", "3", "4", "5"),
                   ArztID=c("A", "C", "D", "A", "D", "A", "D", "D", "C"))
                   
# calculate the care densities
care_density(data)
#>   PatID sum_weights n care_density
#> 1     1           5 3     1.666667
#> 2     2           3 2     3.000000
#> 3     3           3 2     3.000000
#> 4     4          NA 1           NA
#> 5     5          NA 1           NA
```
