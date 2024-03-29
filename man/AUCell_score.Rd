% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scLNC_0_genetic.R
\name{AUCell_score}
\alias{AUCell_score}
\title{Calculates lncRNA units activity score.}
\usage{
AUCell_score(object, lnclist = NULL, min.target = 0)
}
\arguments{
\item{object}{A scLNC object.}

\item{lnclist}{lncRNA list. Only keep a subset of lncRNA units, defaults to all lncRNA units.}

\item{min.target}{Selected units with higher number of mRNA. Default is 0.}
}
\value{
Return a scLNC object with AUC slot.
}
\description{
Calculates lncRNA units activity score.
}
\examples{
\dontrun{
data(LNCobject2)
LNCobject <- AUCell_score(object = LNCobject2,lnclist = NULL,min.target = 0)
}

}
