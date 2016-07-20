\name{ChIPtest_1.0-package}
\alias{ChIPtest_1.0-package}
\alias{ChIPtest_1.0}
\docType{package}
\title{
Nonparametric Methods for Identifying Differential Enrichment Regions with ChIP-seq Data
}
\description{
Nonparametric Tests to identify the differential enrichment region for two conditions or time-course ChIP-seq data. It includes: data preprocessing function, estimation of a small constant used in hypothesis testing, a kernel-based two sample nonparametric test, two assumption-free two sample nonparametric test.
}
\details{
\tabular{ll}{
Package: \tab ChIPtest_1.0\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2016-07-07\cr
License: \tab GPL ( >=2 \cr
LazyLoad: \tab yes \cr
}
}
\author{
Vicky Qian Wu ; Kyoung-Jae Won ; Hongzhe Li <hongzhe@upenn.edu>
Maintainer: Vicky Qian Wu <wuqian7@gmail.com>
}
\references{
Qian Wu, Kyoung-Jae Won and Hongzhe Li. (2015) Nonparametric Methods for Identifying Differential Enrichment Regions with ChIP-seq Data. \eqn{Cancer Informatics},14 (Suppl 1), 11-22
}

\keyword{ ChIPtest, nonparametric test, ChIP-seq, differential enrichment }

\examples{
data(data1)
data(data4)
Data1=NormTransformation(data1)
Data4=NormTransformation(data4)
tao=est.c(Data1, Data4, max1=5, max4=5)
band=54
TS=TS_twosample(Data1, Data4, tao, band, quant=c(0.9,0.9,0.9))
}
