\name{basal}
\alias{basal}
\alias{x}
\alias{y}
\alias{genes}
\alias{group}

\docType{data}

\title{Microarray data set on breast cancer (basal tumors)} 

\usage{data(basal)}

\description{%
  This gene expression data set is freely available, coming from
  Hess \emph{et al}'s paper. Patients were treated with chemotherapy
  prior to surgery. Patient response to the treatment can be classified
  as either a pathologic complete response (pCR) or residual disease
  (not-pCR). The original data set concerns one hundred thirty-three
  patients with stage I--III breast cancer. 
  
  Following Jeanmougin et al. (2011), we restrict to the twenty-nine
  basal tumors from the available samples, divided into 15 pCR and 14
  not-pCR. This particular subtype of breast cancer, harboring
  homogeneous clinical and pathologic features, shows highly variable
  response to chemotherapy. For illustration purpose, we only keep the
  first ten most differentiated probe sets, which exactly match ten
  genes. We also integrate all the probe sets related to these ten
  genes, regardless of their p-values. The \code{basal} data ends up
  with 381 variables (the probe sets) clustered into 172 groups (the
  genes) for a total of 29 samples (the patients).%
}

\format{%
  The four objects are created when loading \code{data(basal)}:
  \describe{
    \item{\code{x}}{matrix with 381 columns and 29
      rows. The \eqn{k}{k}th row gives the expression levels of the 381
      selected probe sets related to the \eqn{k}{k}th patient.
    }
    \item{\code{y}}{binary vector of size 29 indicating the patient status
      (\code{"pcr"} and \code{"not"}).
    }
    \item{\code{group}}{vector of size 381 indicating the group belonging
     (that is, the gene each probe set is related to).
    }
    \item{\code{genes}}{vector of size 172 with the gene names.
    }  
  }
}

\references{%
  Marine Jeanmougin, Mickael Guedj, Christophe Ambroise (preprint, 2011)
  \emph{Defining a robust biological prior from Pathway Analysis to
  drive Network Inference}, \url{http://arxiv.org/abs/1101.3493}.

  K.R. Hess, K. Anderson, W.F. Symmans, V. Valero, N. Ibrahim,
  J.A. Mejia, D. Booser, R.L. Theriault, U.  Buzdar, P.J. Dempsey,
  R. Rouzier, N. Sneige, J.S. Ross, T. Vidaurre, H.L. Gomez,
  G.N. Hortobagyi, and L. Pustzai (2006). Pharmacogenomic predictor of
  sensitivity to preoperative chemotherapy with Paclitaxel and
  Fluorouracil, Doxorubicin, and Cyclophosphamide in breast cancer,
  \emph{Journal of Clinical Oncology}, vol. 24(26), pp. 4236--4244.
}

\examples{
data(basal)
cat("\nmight take about half a minute...")
grplasso <- group.lasso(x, y, group, family="binomial", lambda.min=1e-4)
coolasso <- coop.lasso(x, y, group, family="binomial", lambda.min=1e-4)
par(mfrow=c(1,2))
plot(grplasso, yvar="group")
plot(coolasso, yvar="group")
}
\keyword{datasets}
