\name{project_to_one_mode}
\alias{project_to_one_mode}

\title{
Project a Bipartite-Network to a Single Mode
}
\description{
This function takes a bipartite network created using the \code{igraph} package and returns an adjacency matrix of one of its underlying modes. By directly allowing sparse matrices it is faster and more RAM efficient than some other available versions.
}
\usage{
project_to_one_mode(g, mode, sparse=TRUE)
}
\arguments{
  \item{g}{
An \code{igraph} object with two different node types, defining a bipartite graph.
  }
  \item{mode}{
Either \code{"rows"} or \code{"cols"}, defining which mode should be calculated.
  }
  \item{sparse}{
Whether to use sparse matrix representations or not. Must be either \code{TRUE} or \code{FALSE}.  
  }
}
\details{
A bipartite graph only has connections between two types of nodes. For example, one type of node may be the providers and the other type may be patients. In health-care databases we would usually see which patient visited which providers, but we would not see any direct links between providers and patients. This type of graph can be projected into two subgraphs. One which only includes the providers and one which only includes the patients. This function efficiently creates adjacency matrices of these projections.

The resulting adjacency matrix is a symmetric square matrix that contains the number of shared patients (in provider level graphs) or the number of shared providers (in the patient level graph). See the examples below.
}
\value{
Returns a sparse adjacency matrix of the specified mode.
}
\references{
Landon, Bruce E., Nancy L. Keating, Michael L. Barnett, Jukka-Pekka Onnela, Sudeshna Paul, A. James O’Malley, Thomas Keegan, and Nicholas A. Christakis. (2012). "Variation in Patient-Sharing Networks of Physicians Across the United States". JAMA 308 (3): 265–73.
}
\author{
Robin Denz
}
\seealso{
\code{\link{project_to_one_mode}}
}
\examples{
library(CareDensity)
library(igraph)

# some arbitrary patient-provider contact data
data <- data.frame(PatID=c("1", "1", "1", "2", "2", "3", "3", "4", "5"),
                   ArztID=c("A", "C", "D", "A", "D", "A", "D", "D", "C"))

# create graph
g <- graph.data.frame(data, directed=FALSE)

# add type
V(g)$type <- bipartite.mapping(g)$type

## NOTE: we use sparse=FALSE here to show the resulting matrix directly,
#        but it would be more efficient to keep it at sparse=TRUE
# project to patient-level
project_to_one_mode(g, mode="rows", sparse=FALSE)

# project to provider-level
project_to_one_mode(g, mode="cols", sparse=FALSE)
}
