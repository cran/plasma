\name{TCGA-ESCA}
\alias{TCGA-ESCA}
\alias{loadESCAdata}
\alias{TCGA-LUSC}
\alias{loadLUSCdata}
\alias{plasmaEnv}
\alias{assemble}
\alias{Outcome}
\alias{m450info}
\docType{data}
\title{
Esophageal carcinoma (ESCA) data or lung squamous cell carcinoma (LUSC)
data from The Cancer Genome Atlas (TCGA).
}
\description{
The \code{TCGA-ESCA} dataset contains the objects \code{assemble},
\code{Outcome}, and \code{m450info} for building the \code{MultiOmics}
object. Because its size exceeds the CRAN limits, the data is stored on
a remote server and must be loaded using the function
\code{loadESCAdata}.

The \code{TCGA-LUSC1}dataset is a parallel object for lung
squamous cell carcinoma (LUSC) data, whihc must be loaded using the
\code{loadLUSCdata} function.
}
\usage{
loadESCAdata(env = plasmaEnv)
loadLUSCdata(env = plasmaEnv)
}
\arguments{
  \item{env}{an environment in which to load the data. The default
    value is a private environment in the package, accessible as
    \code{plasmaEnv}. To make access easier, you can use
    \code{globalenv()} or \code{.GlobalEnv}.}
}
\format{
  The \dQuote{TCGA-ESCA} dataset contains the following:
  \describe{
    \item{\code{assemble}}{A list of 7 different omics dataframes with
      varying numbers of features as rows (D) and varying number of
      patients as columns (N). Note that some of these omics dataframes
      had been manipulated to contain NAs, where these may be complete
      on the GDC Dat Portal from which these data originally came. This
      was done to illustrate the capability of the \code{plasma} package
      on working with missing data.}
    \describe{
    \item{\code{ClinicalBin}}{a dataframe (53x185) of clinical binary values.}
    \item{\code{ClinicalCont}}{a dataframe (6x185) of clinical continuous values.}
    \item{\code{MAF}}{a dataframe (566x184) of minor allele frequencies
      (MAF) that have been converted to binary based on whether they had
      a MAF greater than 0.03 (1) or not (0).}
    \item{\code{Meth450}}{a dataframe (1454x185) of continuous beta
      values from the Illumina Infinium HumanMethylation450 arrays. The
      features in this dataframe have been filtered on mean greater than
      0.15 and a standard deviation greater than 0.3.}
    \item{\code{miRSeq}}{a dataframe (926x166) of continuous counts
      values from microRNA (miRNA) sequencing. The features in this
      dataframe have been filtered on a standard deviation of 0.05.}
    \item{\code{mRNASeq}}{a dataframe (2520x157) of continuous counts
      values from mRNA sequencing data. The features in this dataframe
      have been filtered on a mean greater than 4 and a standard
      deviation greater than 0.7.}
    \item{\code{RPPA}}{a dataframe (192x126) of continuous protein
      expression values from reverse phase protein array (RPPA) assays.}
    }
    \item{\code{Outcome}}{a dataframe (185x5) containing the survival
      outcomes for the patients in \code{assemble}.}
    \item{\code{m450info}}{a dataframe (1454x3) containing gene symbol,
      chromosome number, and genomic coordinate IDs corresponding to the
      features (or \dQuote{probes}) in \code{Meth450}.}
}
}

\source{
  \url{https://portal.gdc.cancer.gov/projects/TCGA-ESCA}
}
\author{Kevin R. Coombes \email{krc@silicovore.com},
  Kyoko Yamaguchi \email{kyoko.yamaguchi@osumc.edu}}
\examples{
fls <- try(loadESCAdata())
if (inherits(fls, "try-error")) {
  stop("Unable to load data from remote server.")
}
}
\keyword{datasets}
