\name{twocols2tworows}
\alias{twocols2tworows}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Function to reshape a data frame or file
}
\description{
  This is a function to help in formatting marker data produced by PediHaplotyper
  for use by other software.
  It takes a data frame or file with two consecutive columns per marker and one 
  row per individual and converts this to two rows per individual and one 
  column per marker. It is more user-friendly than reshape() for this particular 
  case.
}
\usage{
twocols2tworows(source, target = "", sep = "", skip = 0, 
  na.strings = "NA", sep.out = "\t", na.out = "", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{source}{
  Either a data frame or a file name: the data to be transposed. If a file name,
  the file should be readable by read.table() with header=TRUE; all other
  parameters for read.table have the default values of read.table but can be
  specified.
}
  \item{target}{
  Ignored if source is a data frame, else the name of the transposed file
}
  \item{sep}{
  Ignored if source is a data frame, else the separator (passed to read.table)
}
  \item{skip}{
  Ignored if source is a data frame, else the number of lines to skip
  (passed to read.table)
}
  \item{na.strings}{
  Ignored if source is a data frame, else the strings used to represent NA values
  (passed to read.table)
}
  \item{sep.out}{
  Ignored if source is a data frame, else the separator to be used in the
  target file
}
  \item{na.out}{
  Ignored if source is a data frame, else the representation of NA values to be
  used in the target file
}
  \item{\dots}{
  Ignored if source is a data frame, else further parameters to be passed to
  read.table
}
}
\details{
  The source data consist of one column with names of individuals, followed by
  two consecutive columns per marker that contain the marker alleles inherited
  from parent 1 and parent 2. The name of the first column for each marker is 
  the marker name, the name of the second column is ignored, as are the row names.
  This is the format used by PediHaplotyper for its input and output marker data.
  This function reformats the data such that there is only one column per
  marker and two consecutive rows per individual.
  
  If source is a file name the file is read using read.table; it must contain
  a header line, but in contrast to read.table the headers are not modified
  if they are not valid identifiers or if duplicate headers occur.
  The column names of the original and rehaped data frame may be invalid as 
  identifiers; that does not affect any of the PediHaplotyper functions.
}
\value{
The return value is always the reshaped data frame. If source and target are
file names also the target file is created, overwriting any previous file.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(mrkdat) # a data frame with data for two markers and two individuals
s <- transpose(mrkdat)
s # a data frame arranged as the source of twocols2tworows
twocols2tworows(s)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
