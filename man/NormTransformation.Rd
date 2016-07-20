\name{NormTransformation}
\alias{NormTransformation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Variance-stabilizing transformation (VST) procedure
}
\description{
Assume observed data approximately follow Poisson Distribution. Apply VST procedure to transform the data as approximate normal with constant variance of 1
}
\usage{
NormTransformation(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Input ChIP-seq data (counts), which approximately follow Poisson Distribution.
}
}
\details{
Please note the input data can be a single value, a vector or a matrix. If it is a matrix, the default format is N row by M column. Each row represents for one gene, and each column represents for one bin.
}
\value{
After VST transformation, the return data matrix would follow Normal Distribution with a constant variance 1
}

\references{
Qian Wu, Kyoung-Jae Won and Hongzhe Li. (2015) Nonparametric Methods for Identifying Differential Enrichment Regions with ChIP-seq Data. \eqn{Cancer Informatics},14 (Suppl 1), 11-22}

\examples{
data(data1)
Data1=NormTransformation(data1)
}
