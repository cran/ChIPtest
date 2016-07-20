\name{TS_kernel}
\alias{TS_kernel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the Test Statistics for kernel-based nonparametric test.
}
\description{
Get the difference between two conditions. Apply Kernel smoothing to fit a smooth curve. Estimate variance for each gene and improve the estimation of variance based on all the genes. Derive test statistics and get the rank list of all the genes. 
}
\usage{
TS_kernel(data, band, quantile)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
difference matrix between two conditions
}
  \item{band}{
bandwidth used in kernel smoothing
}
  \item{quantile}{
threshold used in variance estimation
}
}
\details{
Note 1: Need to chose a bandwidth. Do not recommend to use cross validation (not gene-specific bandwidth) but chose a fixed biological meaningful bandwidth. A fixed bandwidth which can capture the signal profile and smooth out noise would be recommend. The bandwidth used in reference is 20/280.

Note 2: quantile value is based on the distribution of variance estimation of each gene. Recommend to use histogram to double check the distribution. Default 0.9 = 90 \%
}

\value{

\item{TS }{ Kernel based test statistics after WH transformation. Please refer the details in the reference}
\item{TS_sign }{"+" represent for condition B enriched more than condition A; "-" vice versa}
\item{Tmean }{Original test statistics, which is calculated as integral of square of kernel estimator}
}


\references{
Qian Wu, Kyoung-Jae Won and Hongzhe Li. (2015) Nonparametric Methods for Identifying Differential Enrichment Regions with ChIP-seq Data. \eqn{Cancer Informatics},14 (Suppl 1), 11-22}

\examples{
data(data1)
data(data4)
Data1=NormTransformation(data1)
Data4=NormTransformation(data4)
data=Data4-Data1
band=54
TS=TS_kernel(data, band, quantile=0.9)

}
