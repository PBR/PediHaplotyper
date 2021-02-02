\name{haplotyping_session}
\alias{haplotyping_session}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Function to assign haploblock alleles in pedigrees using generic input files
}
\description{
Function haplotyping_session is a user-friendly wrapper for the separate steps
in the imputation and allele assignment process; only the first 4 parameters
need to be specified.
}
\usage{
haplotyping_session(sessionID, mapfile, pedigreefile, 
	phasedgenofile, oldhballelesfile = "", 
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
	origdatahbfile = "orighb_geno", 
	finaldatahbfile = "finalhb_geno", 
	origdatamrkfile = "origmrk_geno", 
	finaldatamrkfile = "finalmrk_geno", 
	min.allele.freq = 3, 
	mv.count = FALSE, 
	fqparfile = "", 
	FQout = fqparfile != "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sessionID}{
  (Short) text string, will be prefixed to all output filenames
}
  \item{mapfile}{
  Name of the input file containing the linkage map, must be specified
}
  \item{pedigreefile}{
  Name of the input file containing the pedigree and optionally phenotypic data,
  must be specified
}
  \item{phasedgenofile}{
  Name of the input file with phased marker genotypes, must be specified
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
  \item{origdatahbfile}{
  Name of the output file with phased haploblock allele data at the start
  of the process
}
  \item{finaldatahbfile}{
  Name of the output file with phased haploblock allele data at the end
  of the process
}
  \item{origdatamrkfile}{
  Name of the output file with phased marker allele data at the start
  of the process
}
  \item{finaldatamrkfile}{
  Name of the output file with phased marker allele data at the end
  of the process
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
  \item{fqparfile}{
  Name of a FlexQTL parameter file, or "" (default). If a file is specified here
  it is assumed that all input will be in FlexQTL format. This parameter is for
  internal use only; FlexQTL input is handled more conveniently using the function
  fq_haplotyping_session.
}
  \item{FQout}{
  If TRUE and if also a fqparfile is specified the output will be in FlexQTL format,
  else in generic format. By default FQout is TRUE is a fqparfile is specified,
  else FALSE.
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
% ##---- Should be DIRECTLY executable !! ----
% ##-- ==>  Define data, use random,
% ##--	or do  help(data=index)  for the standard data sets.
% 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
