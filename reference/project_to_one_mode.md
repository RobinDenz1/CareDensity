# Project a Bipartite-Network to a Single Mode

This function takes a bipartite network created using the `igraph`
package and returns an adjacency matrix of one of its underlying modes.
By directly allowing sparse matrices it is faster and more RAM efficient
than some other available versions.

## Usage

``` r
project_to_one_mode(g, mode, sparse=TRUE)
```

## Arguments

- g:

  An `igraph` object with two different node types, defining a bipartite
  graph.

- mode:

  Either `"rows"` or `"cols"`, defining which mode should be calculated.

- sparse:

  Whether to use sparse matrix representations or not. Must be either
  `TRUE` or `FALSE`.

## Details

A bipartite graph only has connections between two types of nodes. For
example, one type of node may be the providers and the other type may be
patients. In health-care databases we would usually see which patient
visited which providers, but we would not see any direct links between
providers and patients. This type of graph can be projected into two
subgraphs. One which only includes the providers and one which only
includes the patients. This function efficiently creates adjacency
matrices of these projections.

The resulting adjacency matrix is a symmetric square matrix that
contains the number of shared patients (in provider level graphs) or the
number of shared providers (in the patient level graph). See the
examples below.

## Value

Returns a sparse adjacency matrix of the specified mode.

## References

Landon, Bruce E., Nancy L. Keating, Michael L. Barnett, Jukka-Pekka
Onnela, Sudeshna Paul, A. James O’Malley, Thomas Keegan, and Nicholas A.
Christakis. (2012). "Variation in Patient-Sharing Networks of Physicians
Across the United States". JAMA 308 (3): 265–73.

## Author

Robin Denz

## See also

`project_to_one_mode`

## Examples

``` r
library(CareDensity)
library(igraph)

# some arbitrary patient-provider contact data
data <- data.frame(PatID=c("1", "1", "1", "2", "2", "3", "3", "4", "5"),
                   ArztID=c("A", "C", "D", "A", "D", "A", "D", "D", "C"))

# create graph
g <- graph_from_data_frame(data, directed=FALSE)

# add type
V(g)$type <- bipartite_mapping(g)$type

## NOTE: we use sparse=FALSE here to show the resulting matrix directly,
#        but it would be more efficient to keep it at sparse=TRUE
# project to patient-level
project_to_one_mode(g, mode="rows", sparse=FALSE)
#>   1 2 3 4 5
#> 1 0 2 2 1 1
#> 2 2 0 2 1 0
#> 3 2 2 0 1 0
#> 4 1 1 1 0 0
#> 5 1 0 0 0 0

# project to provider-level
project_to_one_mode(g, mode="cols", sparse=FALSE)
#>   A C D
#> A 0 1 3
#> C 1 0 1
#> D 3 1 0
```
