\name{est.c}
\alias{est.c}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
calculate the biologically relevant value c in the null hypothesis H0: TS=c, in assumption-free nonparametric test
}
\description{
If there is no INPUT experiment (No control), treat the genes with read counts fewer than 5 as the "null genes". Test statistics were calculated based on the those "null genes" and take the average to obtain the value c, which is used in the null hypothesis H0: TS=c
}
\usage{
est.c(data1, data4, max1 = 5, max4 = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data1}{
Data Matrix (after VST) for condition A
}
  \item{data4}{
Data Matrix (after VST) for condition B
}
  \item{max1}{
Threshold used to decide null genes for condition A. Default as 5
}
  \item{max4}{
Threshold used to decide null genes for condition B. Default as 5
}
}
\details{
Data matrix, the default format is N row by M column. Each row represents for one gene, and each column represents for one bin
}

\value{
\item{tao}{value c in the null hypothesis H0: TS=c}
}


\references{
Qian Wu, Kyoung-Jae Won and Hongzhe Li. (2015) Nonparametric Methods for Identifying Differential Enrichment Regions with ChIP-seq Data. \eqn{Cancer Informatics},14 (Suppl 1), 11-22}

\examples{
data(data1)
data(data4)
Data1=NormTransformation(data1)
Data4=NormTransformation(data4)
tao=est.c(Data1, Data4, max1=5, max4=5)
}
