\name{esca-type-data}
\docType{data}
\alias{tfESCA}
\alias{mirESCA}

\title{ESCA type data}
\description{
  The \code{CombinedWeights} object class merges the weight matrices for
  all data sets in a plasma object.
}
\usage{
data(tfESCA)
data(mirESCA)
}
\format{
  Both \code{tfData} and \code{mirESCA} are data frames containng two
  columns. The first column is and \code{ID} column containing the TCGA
  sample barcode for an esophagela cancer sample. The second column,
  called \code{Type} identifies the sample as either "squamous" (for
  likely squamous cell carcinomas that cluster near head and neck
  cancers) or "adeno" (for likley adenocarcinomas that cluster near
  stomach cancers).
}
\source{
  All data supplied here are based upon esophageal cancer data generated
  by the TCGA Research Network (https://www.cancer.gov/tcga).

  The transcription factor classifications of 196 esophageal cancer
  into squamous cell carcinoma or adenocarcinoma are taken from work
  published by Abrams and colleaagues in BMC Genomics.

  The microRNA classifications of 195 esophageal cancer samples
  into squamous cell carcinoma or adenocarcinoma are taken from work
  published by Asiaee and colleaagues in J Comput Biol.

}
\author{Kevin R. Coombes \email{krc@silicovore.com},
  Kyoko Yamaguchi \email{kyoko.yamaguchi@osumc.edu}}
\references{
  Abrams ZB, Zucker M, Wang M, Asiaee Taheri A, Abruzzo LV, Coombes KR.\cr
  \emph{Thirty biologically interpretable clusters of transcription
    factors distinguish cancer type.}\cr
  BMC Genomics. 2018 Oct 11;19(1):738. doi: 10.1186/s12864-018-5093-z.

  Asiaee A, Abrams ZB, Nakayiza S, Sampath D, Coombes KR.\cr
  \emph{Explaining Gene Expression Using Twenty-One MicroRNAs.}\cr
  J Comput Biol. 2020 Jul;27(7):1157-1170. doi: 10.1089/cmb.2019.0321.

}
\keyword{classes}
