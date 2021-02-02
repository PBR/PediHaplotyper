\name{transpose}
\alias{transpose}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Function to transpose a data frame or a file
}
\description{
  This is a function to help in formatting the input and/or output data of 
  PediHaplotyper. It transposes a data frame or file, with some differences
  from t().
}
\usage{
transpose(source, target = "", sep = "", skip = 0, 
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
  This function transposes a data frame such that the original column names go  
  to the first column (which is called "name") and the original contents of
  the first column become the column names. Row names are ignored.
  If source is a file name the file is read using read.table; it must contain
  a header line, but in contrast to read.table the headers are not modified
  if they are not valid identifiers or if duplicate headers occur.
  The column names of the transposed data frame may be invalid as identifiers 
  and/or duplicate column names may occur, depending on the contents of the 
  first column of the source data frame. Invalid column names do not affect  
  any of the PediHaplotyper functions.
}
\value{
The return value is always the transposed data frame. If source and target are
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
mrkdat
transpose(mrkdat)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
