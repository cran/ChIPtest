\name{TS_twosample}
\alias{TS_twosample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Three Nonparametric Test Statistics for two sample ChIP-seq data
}
\description{
It includes three nonparametric test statistics for two sample differential analysis: kernel based nonparametric test, assumption-free nonparametric test with equal variance estimation and unequal variance estimation.  
}
\usage{
TS_twosample(data1, data4, tao, band, quant)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data1}{
data matrix (after VST) for condition A
}
  \item{data4}{
data matrix (after VST) for condition B
}
  \item{tao}{
the biologically relevant value c in the null hypothesis H0: TS=c, in assumption-free nonparametric test
}
  \item{band}{
bandwidth used in kernel smoothing
}
  \item{quant}{
threshold used in variance estimation
}
}

\details{
kernel-based test statistics is the same as "TS_kernel"
}
\value{

\item{TS_kn }{kernel based test statistics}
\item{Deql  }{assumption-free nonparametric test statistics with equal variance}
\item{Dnun  }{assumption-free nonparametric test statistics with unequal variance}
\item{sigma1 }{variance estimation for conditon A under equal variance assumption}
\item{sigma4 }{variance estimation for condition B under unequal variance assumption}
\item{Ts_yvec }{Original statistics, which is calculated as integral of square of kernel estimator}
\item{Dsum }{Original statistics, which is calculated for nonparametric test without smoothing}
\item{Sev  }{ variance estimation under equal variance assumption }
\item{Suv  }{ variance estimation under unequal variance assumption }
\item{Xg   }{ estimation of standard deviation for kernel-based test statistics}

}

\references{
Qian Wu, Kyoung-Jae Won and Hongzhe Li. (2015) Nonparametric Methods for Identifying Differential Enrichment Regions with ChIP-seq Data. \eqn{Cancer Informatics},14 (Suppl 1), 11-22}

\examples{
data(data1)
data(data4)
Data1=NormTransformation(data1)
Data4=NormTransformation(data4)
tao=est.c(Data1, Data4, max1=5, max4=5)
band=54
TS=TS_twosample(Data1, Data4, tao, band, quant=c(0.9,0.9,0.9))
}
