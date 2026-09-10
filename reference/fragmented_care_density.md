# Calculate the Fragmented Care Density for all Patients

This function calculates the Fragmented Care Density Index as defined by
Engels et al. (2024) for each patient in the supplied dataset. Works
well with large patient-sharing networks.

## Usage

``` r
fragmented_care_density(data, pat_col=1, weights, type,
                        by_connection=FALSE, data_frame=TRUE)
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

- weights:

  A `data.frame` containing three columns called `"from"` (character),
  `"to"` (character) and `"weight"` (numeric). The first two columns
  should contain types of providers, thus defining different provider
  connections. All possible non-redundant connections need to be
  specified this way. The `"weight"` column should include the weight
  associated with that connection. When using `by_connection=TRUE` this
  argument can be set to `NULL`, because it won't be needed then. See
  examples for more information.

- type:

  A `data.frame` containing two columns called `"ID"` (containing all
  provider IDs) and `"Type"` (containing the type of the provider). Both
  columns should be character vectors.

- by_connection:

  Either `TRUE` or `FALSE` (default). If `TRUE` this function returns
  the person and connection-specific sums of weights and simple care
  densities instead of returning the fragmented care density directly.
  This may be useful to estimate weights for the `weights` argument.

- data_frame:

  Set this argument to `TRUE` to return a `data.frame` instead of the
  `data.table` format that is used under the hood.

## Details

The Fragmented Care Density is an extension of the classic Care Density
(see
[`care_density`](https://robindenz1.github.io/CareDensity/reference/care_density.md))
and was proposed by Engels et al. (2024). It is also a measure of care
coordination, but it allows a lot more flexibility by using different
weights for different provider-type connections. For example, it may
make sense to weight the amount of patients shared by two general
providers differently than the amount of patients shared by a general
provider and a specialist. Formally, the fragmented care density is
defined as:

\$\$FC_p = \sum\_{j = 1}^{k} w_j \frac{s_j}{n_p(n_p - 1) /2},\$\$

where \\n_p\\ is the number of different providers patient \\p\\ visited
and \\w_j\\ are some connection specific weights. \\k\\ is defined as:

\$\$k = {l \choose 2} + l,\$\$

with \\l\\ being the number of different provider types. Finally,
\\s_j\\ is the sum of the number of patients shared by all doctors of a
specific connection type. See Engels et al. (2024) for more information.

Under the hood, this function uses the `igraph` package to construct a
patient-sharing network from the provided `data` to calculate the
weights. It then uses the `data.table` package to efficiently calculate
the care densities from a resulting edge list with weights.

## Value

Returns a single `data.frame` (or `data.table`) containing output
depending on the specification of the `by_connection` argument.

When `by_connection=FALSE` was used the output only includes the patient
id (`"PatID"`) and the calculated fragmented care densities
(`"fragmented_care_density"`).

When `by_connection=TRUE` was used instead, the output includes the
patient id (`"PatID"`), the connection-type (`"connection"`) the sum of
all weights (`"sum_weights"`), the number of providers seen by each
patient (`"n"`) and the calculated simple care density
(`"care_density"`).

## References

Pollack, Craig Evan, Gary E. Weissman, Klaus W. Lemke, Peter S. Hussey,
and Jonathan P. Weiner. (2013). "Patient Sharing Among Physicians and
Costs of Care: A Network Analytic Approach to Care Coordination Using
Claims Data". Journal of General Internal Medicine 28 (3), pp. 459-465.

Engels, Alexander, Claudia Konnopka, Espen Henken, Martin Härter, and
Hans-Helmut König. (2024). "A Flexible Approach to Measure Care
Coordination Based on Patient-Sharing Networks". BMC Medical Research
Methodology 24 (1), pp. 1-12.

## Author

Robin Denz

## See also

[`care_density`](https://robindenz1.github.io/CareDensity/reference/care_density.md)

## Examples

``` r
library(CareDensity)
library(data.table)
library(igraph)

# some arbitrary patient-provider contact data
data <- data.frame(PatID=c("1", "1", "1", "2", "2", "3", "3", "4", "5"),
                   ArztID=c("A", "C", "D", "A", "D", "A", "D", "D", "C"))

# defining the provider types
d_type <- data.frame(ID=c("A", "C", "D"),
                     Type=c("GP", "GP", "Psychiatrist"))
                     
# defining the connection-specific weights
d_weights <- data.frame(from=c("GP", "GP", "Psychiatrist"),
                        to=c("GP", "Psychiatrist", "Psychiatrist"),
                        weight=c(1.1, 0.8, 1.3))

# calculate the fragmented care densities
fragmented_care_density(data, type=d_type, weights=d_weights)
#>   PatID fragmented_care_density
#> 1     1                1.433333
#> 2     2                2.400000
#> 3     3                2.400000
#> 4     4                      NA
#> 5     5                      NA

# calculate only the connection-specific sums and care-densities per patient
# NOTE: "weights" can be set to NULL here because they won't be used
fragmented_care_density(data, type=d_type, weights=NULL, by_connection=TRUE)
#>   PatID        connection sum_weights n care_density
#> 1     1           GP - GP           1 3    0.3333333
#> 2     1 Psychiatrist - GP           4 3    1.3333333
#> 3     2 Psychiatrist - GP           3 2    3.0000000
#> 4     3 Psychiatrist - GP           3 2    3.0000000
#> 5     4              <NA>          NA 1           NA
#> 6     5              <NA>          NA 1           NA
```
