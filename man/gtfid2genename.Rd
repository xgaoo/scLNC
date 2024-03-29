% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scLNC_0_genetic.R
\name{gtfid2genename}
\alias{gtfid2genename}
\title{Gene ensemble ids converts to gene names.}
\usage{
gtfid2genename(id, gtf.info = scLNCgencode)
}
\arguments{
\item{id}{A gene ensemble ids list.}

\item{gtf.info}{Gene annotation from Gencode.}
}
\value{
Return a gene names list.
}
\description{
Gene ensemble ids converts to gene names.
}
\examples{
\dontrun{
data(scLNCgencode)
gtfid2genename(id=c('ENSG00000251191','ENSG00000204054','ENSG00000235304',
'ENSG00000222041'),gtf.info = scLNCgencode)
}

}
