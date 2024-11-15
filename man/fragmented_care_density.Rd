\name{fragmented_care_density}
\alias{fragmented_care_density}

\title{
Calculate the Fragmented Care Density for all Patients
}
\description{
This function calculates the Fragmented Care Density Index as defined by Engels et al. (2024) for each patient in the supplied dataset. Works well with large patient-sharing networks.
}
\usage{
fragmented_care_density(data, pat_col=1, weights, type,
                        by_connection=FALSE, data_frame=TRUE)
}
\arguments{
  \item{data}{
A \code{data.frame} like object containing exactly two columns. One should include only patient IDs and the other one only provider IDs. Each row should denote one patient-provider contact. Multiple contacts (same rows) are allowed but have no effect on the outcome. Both patient and provider IDs should be unique, which means that one ID may not be in both rows.
  }
  \item{pat_col}{
Specifies which column of \code{data} includes the patient IDs. If the first column contains the patient IDs this should be kept at 1, if the second column contains the patient IDs it should be set to 2.
  }
  \item{weights}{
A \code{data.frame} containing three columns called \code{"from"} (character), \code{"to"} (character) and \code{"weight"} (numeric). The first two columns should contain types of providers, thus defining different provider connections. All possible non-redundant connections need to be specified this way. The \code{"weight"} column should include the weight associated with that connection. When using \code{by_connection=TRUE} this argument can be set to \code{NULL}, because it won't be needed then. See examples for more information.
  }
  \item{type}{
A \code{data.frame} containing two columns called \code{"ID"} (containing all provider IDs) and \code{"Type"} (containing the type of the provider). Both columns should be character vectors.
  }
  \item{by_connection}{
Either \code{TRUE} or \code{FALSE} (default). If \code{TRUE} this function returns the person and connection-specific sums of weights and simple care densities instead of returning the fragmented care density directly. This may be useful to estimate weights for the \code{weights} argument. 
  }
  \item{data_frame}{
Set this argument to \code{TRUE} to return a \code{data.frame} instead of the \code{data.table} format that is used under the hood.
  }
}
\details{
The Fragmented Care Density is an extension of the classic Care Density (see \code{\link{care_density}}) and was proposed by Engels et al. (2024). It is also a measure of care coordination, but it allows a lot more flexibility by using different weights for different provider-type connections. For example, it may make sense to weight the amount of patients shared by two general providers differently than the amount of patients shared by a general provider and a specialist. Formally, the fragmented care density is defined as:

\deqn{FC_p = \sum_{j = 1}^{k} w_j \frac{s_j}{n_p(n_p - 1) /2},}

where \eqn{n_p} is the number of different providers patient \eqn{p} visited and \eqn{w_j} are some connection specific weights. \eqn{k} is defined as:

\deqn{k = {l \choose 2} + l,}

with \eqn{l} being the number of different provider types. Finally, \eqn{s_j} is the sum of the number of patients shared by all doctors of a specific connection type. See Engels et al. (2024) for more information.

Under the hood, this function uses the \code{igraph} package to construct a patient-sharing network from the provided \code{data} to calculate the weights. It then uses the \code{data.table} package to efficiently calculate the care densities from a resulting edge list with weights.
}
\value{
Returns a single \code{data.frame} (or \code{data.table}) containing output depending on the specification of the \code{by_connection} argument.

When \code{by_connection=FALSE} was used the output only includes the patient id (\code{"PatID"}) and the calculated fragmented care densities (\code{"fragmented_care_density"}).

When \code{by_connection=TRUE} was used instead, the output includes the patient id (\code{"PatID"}), the connection-type (\code{"connection"}) the sum of all weights (\code{"sum_weights"}), the number of providers seen by each patient (\code{"n"}) and the calculated simple care density (\code{"care_density"}).
}
\references{
Pollack, Craig Evan, Gary E. Weissman, Klaus W. Lemke, Peter S. Hussey, and Jonathan P. Weiner. (2013). "Patient Sharing Among Physicians and Costs of Care: A Network Analytic Approach to Care Coordination Using Claims Data". Journal of General Internal Medicine 28 (3), pp. 459-465.

Engels, Alexander, Claudia Konnopka, Espen Henken, Martin Härter, and Hans-Helmut König. (2024). "A Flexible Approach to Measure Care Coordination Based on Patient-Sharing Networks". BMC Medical Research Methodology 24 (1), pp. 1-12.
}
\author{
Robin Denz
}
\seealso{
\code{\link{care_density}}
}
\examples{
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

# calculate only the connection-specific sums and care-densities per patient
# NOTE: "weights" can be set to NULL here because they won't be used
fragmented_care_density(data, type=d_type, weights=NULL, by_connection=TRUE)
}
