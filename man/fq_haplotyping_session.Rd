\name{fq_haplotyping_session}
\alias{fq_haplotyping_session}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Function to assign haploblock alleles in pedigrees using FlexQTL input files
}
\description{
Function fq_haplotyping_session is a user-friendly wrapper for the separate steps
in the imputation and allele assignment process; only the first 2 parameters
need to be specified.
}
\usage{
fq_haplotyping_session(sessionID, mapfile, 
fqparfile = "flexqtl.par", 
pedigreefile = "flexqtl.sort", 
phasedgenofile = "mhaplotypes.csv", 
oldhballelesfile = "", 
messagefile = "messages.txt", 
mrkpolymorphismfile = "mrkpolymorphism.dat", 
hballelesfile = "hballeles.dat", 
HSorighballelesfile = "orig_hballeles_byHS.dat", 
HSfinalhballelesfile = "final_hballeles_byHS.dat", 
hbstatisticsfile = "hbstatistics.dat", 
mrkstatisticsfile = "mrkstatistics.dat", 
origpedimaphbfile = "orighb.ped", 
finalpedimaphbfile = "finalhb.ped", 
origpedimapmrkfile = "origmrk.ped", 
finalpedimapmrkfile = "finalmrk.ped", 
origflexqtlhbfiles = "orighb_flexqtl", 
finalflexqtlhbfiles = "finalhb_flexqtl", 
origflexqtlmrkfiles = "origmrk_flexqtl", 
finalflexqtlmrkfiles = "finalmrk_flexqtl", 
min.allele.freq = 3, 
mv.count = FALSE, 
FQout = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sessionID}{
  (Short) text string, will be prefixed to all output filenames
}
  \item{mapfile}{
  Name of the input file containing the linkage map, must be specified
}
  \item{fqparfile}{
  Name of the FlexQTL parameter file, by default "flexqtl.par"
}
  \item{pedigreefile}{
  Name of the input file containing the pedigree and optionally phenotypic data,
  by default "flexqtl.sort". If this file is not available the datafile 
  specified in the FlexQTL parameter file is used.
}
  \item{phasedgenofile}{
  Name of the input file with phased marker genotypes, by default "mhaplotypes.csv"
}
  \item{oldhballelesfile}{
  Name of an optional input file with haploblock alleles,
  the hballeles file from an earlier PediHaplotyper run
}
  \item{messagefile}{
  Name of the output file that will contain messages
}
  \item{mrkpolymorphismfile}{
  Name of the output file with polymorphism information for each marker
}
  \item{hballelesfile}{
  Name of the output file listing the alleles of each haploblock
}
  \item{HSorighballelesfile}{
  Name of the output file listing all haploblock alleles per Half-Sib
  at the start of the process
}
  \item{HSfinalhballelesfile}{
  Name of the output file listing all haploblock alleles per Half-Sib
  at the end of the process
}
  \item{hbstatisticsfile}{
  Name of the output file that lists the number of changes made to haploblock
  alleles during the process
}
  \item{mrkstatisticsfile}{
  Name of the output file that lists the number of changes made to marker
  alleles during the process
}
  \item{origpedimaphbfile}{
  Name of the output file with Pedimap data for visualizing the haploblock
  allele configurations at the start of the process
}
  \item{finalpedimaphbfile}{
  Name of the output file with Pedimap data for visualizing the haploblock
  allele configurations at the end of the process
}
  \item{origpedimapmrkfile}{
  Name of the output file with Pedimap data for visualizing the marker
  allele configurations at the start of the process
}
  \item{finalpedimapmrkfile}{
  Name of the output file with Pedimap data for visualizing the marker
  allele configurations at the end of the process
}
  \item{origflexqtlhbfiles}{
  Name of the output file with the pedigree and phased haploblock alleles at
  the start of the process, in FlexQTL-readable format
}
  \item{finalflexqtlhbfiles}{
  Name of the output file with the pedigree and phased haploblock alleles at
  the end of the process, in FlexQTL-readable format
}
  \item{origflexqtlmrkfiles}{
  Name of the output file with the pedigree and phased marker alleles at
  the start of the process, in FlexQTL-readable format
}
  \item{finalflexqtlmrkfiles}{
  Name of the output file with the pedigree and phased marker alleles at
  the end of the process, in FlexQTL-readable format
}
  \item{min.allele.freq}{
  A marker is considered sufficiently polymorphic if at least two of its alleles
  each occur at least this number of times (default: 3 times) in the pedigree.
  If not, the marker is not used in the haploblock allele assignment.
}
  \item{mv.count}{
  If FALSE (default), the haploblock allele IDs consist only of a sequential
  number. If TRUE the IDs have the form n(x), where n is the sequential number
  and x is the number of missing marker alleles in this haploblock allele.
}
  \item{FQout}{
  If TRUE the output will be in FlexQTL format, else in generic format.
}
}
\details{
  This function takes phased marker data for individuals in a pedigree, and for
  each haploblock (group of tightly linked markers) assigns haploblock alleles
  to all individuals. In the process conflicts are resolved by deleting or
  correcting marker alleles, and missing marker data are imputed, assuming
  that no within-haploblock recombination occurs over the entire pedigree.
  See the publication and the manual (under References) for further details.
}
\value{
  The function does not return a value but creates a series of output files
  in the current working directory
}
\references{
  Voorrips RE, Bink MCAM, Kruisselbrink JW, Koehorst-van Putten HJJ, Van de Weg WE:
  PediHaplotyper: Software for consistent assignment of marker haplotypes in
  pedigrees. Submitted to Molecular Breeding.
  
  Website with manual and example data: https://www.wageningenur.nl/en/show/PediHaplotyper.htm 
}
\author{
  Roeland E. Voorrips, Wageningen University & Research Centre, email roeland.voorrips@wur.nl
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Download the manual and example files from
## https://www.wageningenur.nl/en/show/PediHaplotyper.htm 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
