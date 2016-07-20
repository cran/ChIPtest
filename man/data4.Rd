\name{data4}
\alias{data4}
\docType{data}
\title{
Input data (matrix) for condition B
}
\description{
In order to identify the genes that show differential histone modification levels between the two conditions, condition A and condition B, ChIPtest consider the upstream 5000 bp region and downstream 2000 bp region around the transcription start site (TSS) for each gene and divide the regions into 280 bins of 25 bps. Since the two ChIP-seq samples are usually sequenced at different depths (total number of reads), the counts were rescaled according to the sequencing depth ratio. In this example, suppose that there are 5 genes and for each gene, there are 280 observations. The input data matrix has 5 rows and 280 columns. Each row represents for one gene, and each column represents for number of short reads covered at one bin after rescaling. 
}
\usage{data(data4)}
\format{
  The format is:
 num [1:5, 1:280] 0 0 0 0 0 0 0 0 0 0 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:5] "1" "2" "3" "4" ...
  ..$ : chr [1:280] "V3" "V4" "V5" "V6" ...
}
\source{
T.S. Mikkelsen, et al. Comparative Epigenomic Analysis of Murine and Human Adipogenesis. Cell, 143 (156-169): 1156-1166 (2010)
}
\references{
Qian Wu, Kyoung-Jae Won and Hongzhe Li. (2015) Nonparametric Methods for Identifying Differential Enrichment Regions with ChIP-seq Data. \eqn{Cancer Informatics},14 (Suppl 1), 11-22}

\examples{
data(data4)
}
\keyword{datasets}
