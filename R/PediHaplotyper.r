# PediHaplotyper
# Author: Roeland E. Voorrips,
#         Wageningen University and Research Centre - Plant Breeding
#         P.O. Box 386, 6700AJ Wageningen, the Netherlands
#         email: roeland.voorrips@wur.nl

# Copyright 2015 Roeland E. Voorrips
# PediHaplotyper is distributed under the GNU General Public License (GPL)
# version 2 or later (http://www.gnu.org).
#
# PediHaplotyper is a package for assigning haploblock alleles based on
# phased marker genotypes.
# Development started in November 2013
#
# Citation:
# Voorrips RE, Bink MCAM, Kruisselbrink JW, Koehorst - van Putten HJJ,
#   van de Weg WE (2016) PediHaplotyper: Software for consistent assignment of
#   marker haplotypes in pedigrees. Mol Breeding (2016) 36:119.
#   DOI 10.1007/s11032-016-0539-y

# NOTE:   When the software was originally written, a focalpoint (fp) was a set
#         of tightly linked markers, and a fp haplotype was one of the alleles
#         at a focalpoint. Later on a set of tightly linked markers was called
#         a haploblock (hb). So in this script focalpoint (fp) means the same as
#         haploblock (hb), and a focalpoint haplotype (fphaplotype) the same as
#         a haploblock allele (hballele). The change of names has only partially
#         been implemented, mostly in names of functions and parameters
#         directly used by users



# Function haplotyping_session is a user-friendly wrapper for the separate
# steps that make up the haplotyping process.
# The first 4 parameters specify a sessionID that will be prefixed to
# the output file names, and the 3 mandatory input data files; these 4 parameters
# must be supplied by the user.
# Almost all other parameters are names for input and output files and have
# default values. To suppress an output file (except messagefile), set its name
# to "".
haplotyping_session <- function(
   sessionID, #text string, will be prefixed to all output filenames
   mapfile, #may not be omitted, but in a fq_haplotyping_session it may be ""
   pedigreefile, # may also include phenotypic data
   phasedgenofile,
   oldhballelesfile="", #optional hballeles file from an earlier PediHaplotyper run
   #all output files with default names, will be preceded by sessionID:
   messagefile="messages.txt",
   mrkpolymorphismfile="mrkpolymorphism.dat",
   hballelesfile="hballeles.dat",
   HSorighballelesfile="orig_hballeles_byHS.dat",
   HSfinalhballelesfile="final_hballeles_byHS.dat",
   hbstatisticsfile="hbstatistics.dat",
   mrkstatisticsfile="mrkstatistics.dat",
   origpedimaphbfile="orighb.ped",
   finalpedimaphbfile="finalhb.ped",
   origpedimapmrkfile="origmrk.ped",
   finalpedimapmrkfile="finalmrk.ped",
   origdatahbfile="orighb_geno", #.dat and for FlexQTL .par and .map will be appended
   finaldatahbfile="finalhb_geno", #same
   origdatamrkfile="origmrk_geno", #same
   finaldatamrkfile="finalmrk_geno", #same
   #some non-file parameters, all with default values:
   min.allele.freq=3, #see function read_mhaplotypes - specifies minimum required polymorphism of markers
   mv.count=FALSE, #if TRUE, the number of missing marker scores is shown after the hballele number
   fqparfile="", # if an existing file is specified, all input is assumed to be in flexqtl format;
   #             in that case use function fq_haplotyping_session
   FQout=fqparfile!="" #if TRUE, output is in FlexQTL format including a
   #                    FlexQTL parameter file
  ) {
  # this function performs a complete workflow with default settings
  if (FQout && fqparfile=="") FQout <- FALSE #parfile needed to write FQ output
  if (mv.count) mv.sep <- "()" else mv.sep <- NA
  if (messagefile == "") messagefile <- "messages.txt"
  msgfile <- paste(sessionID, messagefile, sep="_")
  write(paste("haplotyping session started on", date()), file=msgfile)

  cat("reading datafiles ...\n")
  #read map:
  maptype <- checkmapfiletype(mapfile)
  if (fqparfile != "" && maptype == 3) {
    stop("mapfile must be specified")
  }
  map <- switch(maptype,
                read_map_from_table(mapfile),
                read_map_from_extended_FQmap(mapfile),
                read_map_from_mhaplotypes(phasedgenofile))
  if (is.data.frame(map)) {
    #map is a map data.frame, we need to add usemarker:
    usemarker <- rep(TRUE, nrow(map))
  } else {
    #map is a list as returned by read_map_from_table, split:
    usemarker <- map$usemarker
    map <- map$map
  }
  #we delete the not-usemarkers from the map: this will also exclude them
  #from the marker data set later on:
  map <- map[usemarker,]
  rm(usemarker)
  #note that usemarker here is a logical vector indicating for each marker in
  #the map file if it will be used in analysis or not (a user selection, only
  #available from read_map_from_table).
  #This is different from the usemarker returned by read_phasedgenofile: that
  #indicates which markers are sufficiently polymorphic to be used.

  e <- checkMapOrder(map, allowColocalizedMarkers=TRUE)
  if (length(e) > 0)
    stop("map should be organized per chromosome in ascending order")

  chb <- checkHaploblocks(map)
  if (length(chb) > 0) {
    write(paste("\nThe following", length(chb),
                "haploblocks do not form one marker block:"), file=msgfile,
          append=TRUE)
    cat(paste("The following", length(chb),
                "haploblocks do not form one marker block:\n"))
    for (f in chb) {
      write(f, file=msgfile, append=TRUE)
      cat(f); cat("\n")
    }
    stop("not all haploblocks are composed of one marker block")
  }

  #read pedigree with nuisance and trait variables:
  ped <- read_pedigree(pedfilename=pedigreefile, parfilename=fqparfile)
  #Note: this will also correctly read non-flexqtl input
  if (is.character(ped)) stop(paste("error in pedigree:",ped))

  #read oldhballeles if specified:
  if (oldhballelesfile == "") oldhballeles <- NULL else
    oldhballeles <- read_oldhballeles(oldhballelesfile)

  #read phased marker data:
  mhaplo <- read_phasedgenofile(
    map=map, ped=ped, min.allele.freq=min.allele.freq,
    phasedgenofile=phasedgenofile,
    diagnostic_filename=paste(sessionID, mrkpolymorphismfile, sep="_"),
    fq=(fqparfile!=""))
  write(c("", mhaplo$messages), file=msgfile, append=TRUE)

  #calculate the initial haploblock alleles (haplotypes):
  ih <- get_initial_haplotypes(mrkallele.arr=mhaplo$mrkallele.arr,
                               map=map, usemarker=mhaplo$usemarker,
                               oldhballeles=oldhballeles,
                               alleledata=mhaplo$alleledata)
  write(c("", ih$messages), file=msgfile, append=TRUE)
  fp_haplotypes <- ih$fp_haplotypes
  hapnrarr <- ih$hapnrarr
  mhaplo$alleledata <- ih$alleledata
  rm(ih)
  hbmap <- getHaploblockmap(fp_haplotypes)
  msg <- checkMapOrder(hbmap, allowColocalizedMarkers=FALSE)
  if (length(msg) > 0) {
    msg <- c(" ",
             "The following set(s) of haploblocks have identical positions:",
             msg)
    write(msg, file=msgfile, append=TRUE)
    cat(paste(msg, "\n", sep=""))
  }

  #get a list with all HS families:
  HSfam <- get.all.HSfamilies(ped)

  #write all haploblock alleles found initially in each HS family:
  if (HSorighballelesfile != "") {
    cat("writing initial HS family haploblock alleles file ...\n")
    writeHSfamHaplotypes(HSfam=HSfam, fp_haplotypes=fp_haplotypes,
                         hapnrarr=hapnrarr,
                         filename=paste(sessionID, HSorighballelesfile, sep="_"),
                         mv.sep=mv.sep)
  }

  #the actual calculation of haploblock alleles:
  cat("calculation haploblock alleles ...\n")
  procAFP <- processAllFocalpoints(fp_haplotypes, hapnrarr, HSfam)
  write(c("", procAFP$messages), file=msgfile, append=TRUE)

  cat("writing all requested output files ...\n")
  #write all alleles per haploblock:
  if (hballelesfile != "") {
    write.all.haploblockalleles(filename=paste(sessionID, hballelesfile, sep="_"),
                         procAFP$fp_haplotypes, procAFP$hapnrarr,
                         alleledata=mhaplo$alleledata,
                         map=map, na.char="-",
                         mv.sep=mv.sep)
  }

  # For debugging: write the final haploblock alleles grouped:
  # groupNwrite.all.haploblockalleles(procAFP$fp_haplotypes, procAFP$hapnrarr,
  #   filename=paste(sessionID, "hballeles_grouped.dat", sep="_"))

  #write all haploblock alleles found finally in each HS family:
  if (HSfinalhballelesfile != "") {
    writeHSfamHaplotypes(HSfam, procAFP$fp_haplotypes, procAFP$hapnrarr,
                         filename=paste(sessionID, HSfinalhballelesfile, sep="_"),
                         mv.sep=mv.sep)
  }

  #calculate the original and final marker allele data arrays:
  origmarkerarray <- getMarkerarray(fp_haplotypes, hapnrarr,
                                    mhaplo$alleledata,
                                    mhaplo$mrkallele.arr)
  finalmarkerarray <- getMarkerarray(procAFP$fp_haplotypes, procAFP$hapnrarr,
                                     mhaplo$alleledata,
                                     mhaplo$mrkallele.arr)
  mrkconvergence <- procAFP$convergence[match(map$hb, names(fp_haplotypes))]

  #write the haploblock files:
  hbmap <- getHaploblockmap(fp_haplotypes)

  fphapcol <- allelecolors(old=hapnrarr,
                           new=procAFP$hapnrarr,
                           missing=1,
                           convergence=procAFP$convergence)
  allelecolors.statistics(colors=fphapcol, map=hbmap,
                          filename=paste(sessionID, hbstatisticsfile, sep="_"))

  #create a version of missing value symbols for haploblocks for use in writing
  #the output file
  if (is.na(mv.sep)) missing <- 1 else {
    missing <- character(length(fp_haplotypes))
    for (hb in seq_along(fp_haplotypes))
      missing[hb] <- get_hballeleID(hb, 1, fp_haplotypes, mv.sep)
  }

  if (origpedimaphbfile != "" || origdatahbfile != "") {
    #create a version of hapnrarr with the desired type of hballele ID's:
    if (is.na(mv.sep)) {
      hna <- hapnrarr
    } else {
      hna <- conv.hballelenames(x=hapnrarr, fp_haplotypes, mv.sep)
    }
  }
  if (origpedimaphbfile != "") {
    writePedimapFile(filename=paste(sessionID, origpedimaphbfile, sep="_"),
                     map=hbmap,
                     ped=ped,
                     allelearr=hna,
                     missing=missing)
  }
  if (origdatahbfile != "") {
    writeDatafiles(map=hbmap,
                   ped=ped,
                   allelearr=hna,
                   missing=missing,
                   outfiles=paste(sessionID, "_", origdatahbfile, sep=""),
                   origFQparfile=ifelse(FQout, fqparfile, ""),
                   usemarker=TRUE)
  }

  if (finalpedimaphbfile != "" || finaldatahbfile != "") {
    #create a version of procAFP$hapnrarr with the desired type of hballele ID's:
    if (is.na(mv.sep)) {
      hna <- procAFP$hapnrarr
    } else hna <- conv.hballelenames(x=procAFP$hapnrarr, fp_haplotypes, mv.sep)
  }
  if (finalpedimaphbfile != "") {
    writePedimapFile(filename=paste(sessionID, finalpedimaphbfile, sep="_"),
                     map=hbmap,
                     ped=ped,
                     allelearr=hna,
                     missing=missing,
                     allelecolors=fphapcol)
  }
  if (finaldatahbfile != "") {
    writeDatafiles(map=hbmap,
                 ped=ped,
                 allelearr=hna,
                 missing=missing,
                 outfiles=paste(sessionID, "_", finaldatahbfile, sep=""),
                 origFQparfile=ifelse(FQout, fqparfile, ""),
                 usemarker=TRUE)
  }
  if ((origdatahbfile != "" || finaldatahbfile != "") && !FQout) {
    #write the hbmapfile once
    write.table(hbmap, file=paste(sessionID, "_hbmap.dat", sep=""), sep="\t",
                quote=FALSE, na="", col.names=TRUE, row.names=FALSE)
  }

  #write the marker files:
  mrkhapcol <- allelecolors(old=origmarkerarray[[1]],
                            new=finalmarkerarray[[1]],
                            missing=NA,
                            convergence=mrkconvergence,
                            usemarker=finalmarkerarray$usemarker)
  allelecolors.statistics(colors=mrkhapcol, map=map,
                          filename=paste(sessionID, mrkstatisticsfile, sep="_"))
  if (origpedimapmrkfile !="") {
    writePedimapFile(filename=paste(sessionID, origpedimapmrkfile, sep="_"),
                     map=map,
                     ped=ped,
                     allelearr=origmarkerarray[[1]],
                     missing=NA)
  }
  if (finalpedimapmrkfile != "") {
    writePedimapFile(filename=paste(sessionID, finalpedimapmrkfile, sep="_"),
                     map=map,
                     ped=ped,
                     allelearr=finalmarkerarray[[1]],
                     missing=NA,
                     allelecolors=mrkhapcol)
  }
  if (origdatamrkfile != "") {
    writeDatafiles(map=map,
                 ped=ped,
                 allelearr=origmarkerarray[[1]],
                 missing=NA,
                 outfiles=paste(sessionID, "_", origdatamrkfile, sep=""),
                 origFQparfile=ifelse(FQout, fqparfile, ""),
                 usemarker=TRUE)
  }
  if (finaldatamrkfile != "") {
    writeDatafiles(map=map,
                 ped=ped,
                 allelearr=finalmarkerarray[[1]],
                 missing=NA,
                 outfiles=paste(sessionID, "_", finaldatamrkfile, sep=""),
                 origFQparfile=ifelse(FQout, fqparfile, ""),
                 usemarker=TRUE)
  }
  if ((origdatamrkfile != "" || finaldatamrkfile != "") && !FQout) {
    #write the mrkmapfile once
    write.table(map, file=paste(sessionID, "_mrkmap.dat", sep=""), sep="\t",
                quote=FALSE, na="", col.names=TRUE, row.names=FALSE)
  }
  if ((origdatahbfile != "" || finaldatahbfile != "" ||
         origdatamrkfile != "" || finaldatamrkfile != "") && !FQout) {
    #write the pedigreefile once
    write.table(ped, file=paste(sessionID, "_pedigree.dat", sep=""), sep="\t",
                quote=FALSE, na="", col.names=TRUE, row.names=FALSE)
  }

  cat(paste("haplotyping session", sessionID, "finished.\n"))
} #haplotyping_session

# Function fq_haplotyping_session is a user-friendly wrapper for the separate
# steps that make up the haplotyping process, to be used if the input files are
# in FlexQTL format.
# The first 2 parameters specify a sessionID that will be prefixed to
# the output file names, and the map file; these parameters
# must be supplied by the user.
# Almost all other parameters are names for input and output files and have
# default values. To suppress an output file (except messagefile), set its name
# to "".
fq_haplotyping_session <- function(
  sessionID,
  mapfile,
  #three flexqtl input files with default names:
  fqparfile="flexqtl.par",
  pedigreefile="flexqtl.sort", #if not available, the datafile from flexqtl.par is used
  phasedgenofile="mhaplotypes.csv",
  oldhballelesfile="",
  #all output files with default names, will be preceded by sessionID
  messagefile="messages.txt",
  mrkpolymorphismfile="mrkpolymorphism.dat",
  hballelesfile="hballeles.dat",
  HSorighballelesfile="orig_hballeles_byHS.dat",
  HSfinalhballelesfile="final_hballeles_byHS.dat",
  hbstatisticsfile="hbstatistics.dat",
  mrkstatisticsfile="mrkstatistics.dat",
  origpedimaphbfile="orighb.ped",
  finalpedimaphbfile="finalhb.ped",
  origpedimapmrkfile="origmrk.ped",
  finalpedimapmrkfile="finalmrk.ped",
  origflexqtlhbfiles="orighb_flexqtl", #.par, .map and .dat will be appended
  finalflexqtlhbfiles="finalhb_flexqtl",
  origflexqtlmrkfiles="origmrk_flexqtl", #same
  finalflexqtlmrkfiles="finalmrk_flexqtl", #same
  #some non-file parameters, all with default values:
  min.allele.freq=3, #see function read_mhaplotypes
  mv.count=FALSE, #if TRUE, the number of missing marker scores is shown after the hballele number
  FQout=TRUE
) {
  # This function takes flexqtl input with (mostly) default file names
  # and calls the generic haplotyping_session with appropriate parameters
  haplotyping_session(
    sessionID=sessionID, mapfile=mapfile,
    oldhballelesfile=oldhballelesfile,
    pedigreefile=pedigreefile,
    phasedgenofile=phasedgenofile,
    #all output files with default names, will be preceded by sessionID
    messagefile=messagefile,
    mrkpolymorphismfile=mrkpolymorphismfile,
    hballelesfile=hballelesfile,
    HSorighballelesfile=HSorighballelesfile,
    HSfinalhballelesfile=HSfinalhballelesfile,
    hbstatisticsfile=hbstatisticsfile,
    mrkstatisticsfile=mrkstatisticsfile,
    origpedimaphbfile=origpedimaphbfile,
    finalpedimaphbfile=finalpedimaphbfile,
    origpedimapmrkfile=origpedimapmrkfile,
    finalpedimapmrkfile=finalpedimapmrkfile,
    origdatahbfile=origflexqtlhbfiles,
    finaldatahbfile=finalflexqtlhbfiles,
    origdatamrkfile=origflexqtlmrkfiles,
    finaldatamrkfile=finalflexqtlmrkfiles,
    min.allele.freq=min.allele.freq,
    mv.count=mv.count,
    fqparfile=fqparfile,
    FQout=FQout
  )
} #fq_haplotyping_session

read_FQparfile <- function(parfilename="flexqtl.par") {
  #return value: list with items
  # $datafile will only be used (to read pedigree with nuisance and trait names)
  #           if the flexQTL.sort file is not available
  # $mapfile (also not used, we need a mapfile with extra haploblock column)
  # $ncolN number of nuisance trait columns, may be 0
  datafile <- ""
  misvalT <- "-"
  mapfile <- ""
  ncolN <- NA
  ncolT <- NA
  namesN <- character(0)
  namesT <- character(0)
  i <- 1
  line <- readLines(parfilename, warn=F)
  while (i<=length(line) && (is.na(ncolN) || is.na(ncolT) ||
                               datafile == "" || mapfile == "" ||
                               length(namesN) == 0 ||
                               length(namesT) == 0) ) {
    words <- unlist(strsplit(line[i], split="[[:blank:]]"))
    words <- words[nchar(words) > 0]
    if (length(words) > 0 && words[1][1] != ";") {
      #not a comment line
      words[1] <- toupper(words[1])
      if (words[1]=="NCOLN" && !is.na(as.integer(words[2])) )
        ncolN <- as.integer(words[2])
      if (words[1]=="NCOLT" && !is.na(as.integer(words[2])) )
        ncolT <- as.integer(words[2])
      if (words[1]=="DATAFILE" && length(words)>1)
        datafile <- words[2]
      if (words[1]=="MAPFILE" && length(words)>1)
        mapfile <- words[2]
      if (words[1]=="NAMESN" && length(words)>1)
        namesN <- words[-1]
      if (words[1]=="NAMEST" && length(words)>1)
        namesT <- words[-1]
      if (words[1]=="MISVALT" && length(words)>1)
        misvalT <- words[-1]
    }
    i <- i+1
  }
  if (datafile == "") stop(paste("No datafile specified in ",parfilename))
  if (mapfile == "") stop(paste("No mapfile specified in ",parfilename))
  if (is.na(ncolN)) stop(paste("ncolN not found in ",parfilename))
  if (is.na(ncolT)) stop(paste("ncolT not found in ",parfilename))
  #check and correct if nuisance or trait names not defined or incorrect length:
  lN <- length(namesN)
  if (lN < ncolN) {
    namesN[(lN+1):ncolN] <- paste("nuisancevar", (lN+1):ncolN, sep="_")
  }
  if (ncolN > 0) namesN <- namesN[1:ncolN]
  lT <- length(namesT)
  if (lT < ncolT) {
    namesT[(lT+1):ncolT] <- paste("traitvar",(lT+1):ncolT, sep="_")
  }
  if (ncolT > 0) namesT <- namesT[1:ncolT]
  list (datafile=datafile, mapfile=mapfile, ncolN=ncolN, ncolT=ncolT,
        namesN=namesN, namesT=namesT, misvalT=misvalT)
} #read_FQparfile

read_pedigree <- function(pedfilename="flexqtl.sort",
                          parfilename="flexqtl.par") {
  #20150902: if flexqtl parfile is specified and present, then only the
  #pedigree and the nuisance and trait columns are read from flexqtl.sort
  #or the original flexqtl datafile.
  #if no flexqtl parfile, a generic pedfile is read with all columns as trait columns
  #(even if parfilename is specified but no such file present)
  #The generic format has a header line with "name", "parent1", "parent2" and
  #trait names for any further columns, has "" as NA-string and is tab-separated
  #In both the flexqtl and generic format it is possible to exclude individuals
  #by placing a ";" in front of the row
  if (!is.na(parfilename) && !(parfilename=="") && file.exists(parfilename)) {
    #flexQTL parameter and data file available
    FQpar <- read_FQparfile(parfilename)
    if (!file.exists(pedfilename)) pedfilename <- FQpar$datafile
    if (!file.exists(pedfilename)) stop("read_pedigree: pedigreefile not found")
    na.trait <- FQpar$misvalT
    ped <- read.table(pedfilename, header=FALSE, as.is=TRUE, na.strings=na.trait,
                      comment.char=";")
    ped <- ped[,2:(4+FQpar$ncolN+FQpar$ncolT)]
    names(ped) <- c("name", "parent1", "parent2", FQpar$namesN, FQpar$namesT)
    ped$parent1[ped$parent1=="0"] <- NA
    ped$parent2[ped$parent2=="0"] <- NA
  } else {
    #no FQ parfile, we read a file containing only the pedigree and
    #possibly trait values:
    if (!file.exists(pedfilename)) stop("read_pedigree: pedigreefile not found")
    con <- file(pedfilename, "r", blocking = FALSE)
    headers <- strsplit(readLines(con, 1), split="\t")[[1]]
    close(con)
    ped <- read.table(pedfilename, skip=1, header=FALSE, as.is=TRUE, na.strings="",
                      comment.char=";", sep="\t")
    headers[1:3] <- c("name", "parent1", "parent2")
    names(ped) <- headers
  }
  ped <- sortPedigree(ped, colInd=1, colPar1=2, colPar2=3,
                           parentsFirst=TRUE, semiFounders=TRUE)
  ped
} #read_pedigree

read_map_from_mhaplotypes <- function(filename="mhaplotypes.csv") {
  # read the map from line 19-20:
  tmp <- read.table(filename, skip=18, header=F, na.strings="-",sep=",",nrows=2)
  #read the marker names from line 21:
  mhaplo <- read.table(filename, skip=20, header=T, na.strings="-",sep=",",
                       as.is=TRUE, nrows=1)
  map <- data.frame(marker=names(mhaplo)[3:length(mhaplo)],
                     chrom=as.integer(tmp[1,3:length(tmp)]),
                     pos=as.numeric(tmp[2,3:length(tmp)]))
  #add the haploblock here by rounding of position:
  chrchar <- max(nchar(map$chrom))
  chr <- substr(rep("00000", nrow(map)), 1, chrchar - nchar(map$chrom) + 1)
  chr <- paste(substr(chr, 2, chrchar), map$chrom, sep="")
  flp <- as.integer(floor(map$pos))
  pschar <- max(nchar(flp))
  ps <- substr(rep("00000000", nrow(map)), 1, pschar - nchar(flp) + 1)
  ps <- paste(substr(ps, 2, pschar), flp, sep="")
  map$hb <- paste("hb", chr, "_", ps, sep="")
  map
} #read_map_from_mhaplotypes

checkmapfiletype <- function(mapfilename) {
  #if mapfilename=="", mapfiletype=3 (read map from phasedgenotypesfile,
  #valid only for flexqtl input). Else if the first non-empty line that
  #does not start with a ';' has only two words before the first ';' (if that
  #occurs) and the first word is 'group' or 'chrom' the mapfiletype is 2
  #extended FQ mapfile), else mapfiletype=1 (tab-separated table)
  if (mapfilename=="") return(3)
  line <- readLines(mapfilename, warn=FALSE)
  i <- 1 #line number
  words <- character(0)
  while (i<=length(line) && length(words)==0) {
    words <- unlist(strsplit(line[i], split="[[:blank:]]"))
    words <- words[nchar(words)>0]
    if (length(words)>0 && substr(words[1], 1, 1) == ";") words <- character(0)
  }
  if (length(words)<2) stop (paste(mapfilename,": unrecognized map file format"))
  W <- toupper(words[1])
  if (W!="GROUP" && W!="CHROM") return(1)
  if (length(words) == 2) return(2) else {
    if (substr(words[3], 1, 1) ==';') return(2) else return(1)
  }
} #checkmapfiletype

read_map_from_extended_FQmap <- function(mapfilename) {
  #each linkage group starts with a header line with the keyword "GROUP" or
  #"CHROM" (upper/lower case not relevant) and the name of the group
  #following the header lines are one or more lines, each with a locus name,
  #position AND haploblock name (this last is not present in FQ mapfile)
  #The loci must be listed in ascending order of position (although that is
  #not checked here). Words and numbers must be separated by one or more spaces
  #and/or tab characters
  #Empty lines and lines starting with ";" are ignored
  line <- line <- readLines(mapfilename, warn=FALSE)
  i <- 1 #line number
  chromname <- ""
  loc <- 0
  locnames <- character(0)
  chrom <- character(0)
  position <- double(0)
  hb <- character(0)

  while (i<=length(line)) {
    words <- unlist(strsplit(line[i], split="[[:blank:]]"))
    words <- words[nchar(words)>0]
    words
    if (length(words)>0 && words[1][1] != ";") {
      W <- toupper(words[1])
      if (W=="GROUP" || W=="CHROM") {
        if (length(words)<2)
          stop (paste(mapfilename,": GROUP or CHROM without group name on line",i))
        chromname <- words[2]
      }
      else {
        if (chromname=="")
          stop(paste(mapfilename, ": locus",words[1],
                     "appears before GROUP header on line", i))
        if (length(words)<3 ||
              is.na(suppressWarnings(as.double(words[2]))) )
          stop(paste(mapfilename, ": locus name", words[1],
                      "without valid position or haploblock"))
        loc <- loc + 1
        locnames[loc] <- words[1]
        chrom[loc] <- chromname
        position[loc] <- as.double(words[2])
        hb[loc] <- words[3]
        if (loc>1 && chrom[loc] == chrom[loc - 1] &&
            position[loc] < position[loc - 1]) {
          stop(paste(mapfilename, ": locus name", words[1],
                     "position lower than previous locus"))
        }
      }
    } # line not empty or comment
    i <- i + 1
  }
  data.frame(marker=locnames, chrom=chrom, pos=position, hb=hb)
} #read_map_from_extended_FQmap

read_map_from_table <- function(mapfilename,
                                allowDuplicatedMarkers=FALSE,
                                allowColocalizedMarkers=TRUE,
                                roundPosition=3,
                                extraColumns=FALSE) {
  #mapfile is a tab-separated file with a table with at least
  #columns marker, chrom, pos, fp (or alternative names, see code; columns
  #may be in any order); map must be in ascending map order
  #Also a column use or omit, containing 0/1 values specifying marker(s) to
  #use or to omit may be present
  #allowDuplicatedMarkers: if FALSE no two marker names may be identical
  #allowColocalizedMarkers: if FALSE no two markers may be at the same
  #                         position on the same chromosome
  #roundPosition: an integer or NA. If not NA, positions are rounded to this
  #               number of decimals
  #extraColumns: if TRUE: any extra columns besides the required marker,
  #              chromosome, position, haploblock and optional use/omit
  #              are also read, else they are ignored
  #Return value: a list with items
  # $map: a data frame containing the 4 required columns, and any other
  #       columns if extraColumns is TRUE; the map contains only the markers
  #       indicated by the use or omit column is present, and the use or omit
  #       column itself is not present any more.
  # $usemarker: the original use column, or (1-omit) translated to TRUE/FALSE,
  #       or all TRUE's if no use or omit was present; this may be used to
  #       select markers from any other data files if these are in map order.
  #NOTE that markers can also be excluded by commenting out the line
  #     i.e. by putting the comment character ';' in front of the line
  #     In that case the map will be shorter than the original
  #     map, and may not match any more with other data files.
  #     The preferred way to omit markers is therefore the use or omit column.

  #This function is forgiving for some format errors:
  #it will clean leading and trailing blanks in column captions and data
  #and it accepts lines having unequal numbers of fields.
  #All lines are padded with NA fields such that all have equal length,
  #any columns after the one with the last non-blank caption are ignored.
  #This will avoid problems where users have added (invisible) spaces or tabs
  #Comments can be added at the end of any line if preceded by a tab and/or the
  #comment character ';'
  #Lines consisting of only blanks (spaces and tabs) are ignored, and also
  #lines where the first non-blank character is the commentchar

  #first find header line and parse it to count columns (readLines)
  #assign all columns to their correct caption,
  #and delete all columns for which there was no caption
  #convert all column names to their preferred version
  #check for correctly sorted, no NA in required columns
  #check for duplicate marker names and names containing spaces
  #delete all rows where use==0 or omit==1, and then column use or omit
  #round position to 3 decimals
  #optionally: check after rounding for co-localized markers
  sep <- "\t"; commentchar <- ";"
  colnames <- data.frame(
    prefname = factor(c(rep("marker", 5),
                         rep("chrom", 4),
                         rep("pos", 3),
                         rep("use", 4),
                         rep("hb", 4))),
    allowedname = c("mrk", "mark", "marker", "loc", "locus",
                     "chr", "chrom", "chromosome", "group",
                     "pos", "position", "cm",
                     "use", "include", "omit", "exclude", #the last two are treated different from the first two!
                     "hb", "fp", "haploblock", "focalpoint")
  )
  #find and parse header line:
  con <- file(mapfilename, "r")
  line <- 0
  while(TRUE) {
    line <- line + 1
    s <- readLines(con, n=1, warn=FALSE)
    if (length(s) != 1) stop("no header line found")
    s <- gsub("^\\s+|\\s+$", "", s) #trim leading and trailing whitespace
    if (nchar(s) > 0 && substring(s,1,1) != commentchar) break()
  }
  close(con)
  captions <- unlist(strsplit(s, split=sep, fixed=TRUE))
  captions <- gsub("^\\s+|\\s+$", "", captions)
  startcomment <- which(substring(captions, 1, 1) == commentchar)
  if (length(startcomment) > 0) captions <- captions[1:(min(startcomment) - 1)]
  #at least one caption remains
  capchar <- nchar(captions)
  lastcol <- max(which(capchar > 0))
  if (lastcol == 0) stop("header line contains no non-blank column headers")
  captions <- tolower(captions[1:lastcol])
  #remove possible comment included in last caption:
  commentat <- regexpr(paste(" ",commentchar,sep=""), captions[lastcol])[1]
  if (commentat > 0)
    captions[lastcol] <- gsub("^\\s+|\\s+$", "", substring(captions[lastcol], 1, commentat-1))
  # this always has at least one non-blank character left!
  if (sum( (capchar==0)[1:lastcol]) > 0) stop("blank column headers in header line")
  spaceat <- regexpr(" ", captions) #vector with pos of 1st space (or -1) for each caption
  if (sum(spaceat > 0) > 0) stop("some column headers contain spaces")
  # now we have a number of consecutive non-blank column captions without spaces
  #next we check the maximum number of fields in any data line:
  fields <- count.fields(mapfilename, sep=sep, skip=line,
                         blank.lines.skip=TRUE, comment.char=commentchar)
  tmpnames <- paste("V", 1:max(fields), sep="")
  # read the rest of the file without risking the first column to become rownames:
  map <- read.table(mapfilename, skip=line, header=FALSE, col.names=tmpnames,
                    sep=sep, comment.char=commentchar, na.strings="",
                    blank.lines.skip = TRUE, strip.white=TRUE,
                    fill=TRUE)

  if (length(map) < lastcol) stop("more column headers than data columns")
  #read and check the four required columns:
  newmap <- NULL
  mapnames <- c("marker", "chrom", "pos", "hb")
  for (cname in mapnames) {
    fc <- which(captions %in% colnames$allowedname[colnames$prefname==cname])
    if (length(fc) != 1) stop(paste("one column", cname, "required"))
    if (sum(is.na(map[, fc])) > 0)
      stop(paste("no missing values allowed in column", captions[fc]))
    if (is.null(newmap)) newmap <- data.frame(map[, fc]) else
      newmap <- cbind(newmap, map[, fc])
  }
  names(newmap) <- mapnames
  r1 <- 1:(nrow(newmap) - 1); r2 <- r1+1
  # check chromosomes sorted:
  e <- which(newmap$chrom[r2] != newmap$chrom[r1])
  if (length(e) != length(unique(newmap$chrom)) - 1)
    stop("map should be organized per chromosome in ascending order")
  #check markers ascending:
  e <- checkMapOrder(newmap, allowColocalizedMarkers)
  if (length(e) > 0)
    stop("map should be organized per chromosome in ascending order")
  # check for haploblocks forming blocks of successive markers is done elsewhere
  if (!is.na(roundPosition)) newmap$pos <- round(newmap$pos, roundPosition)
  #include any optional columns if requested (this will not copy the
  #use/omit column):
  if (extraColumns) {
    for (cname in setdiff(unique(captions), colnames$allowedname)) {
      fc <- which(captions == cname)
      if (length(fc) != 1) stop(paste("multiple columns", cname))
      if (fc <= length(map)) {
        newmap <- cbind(newmap, map)
        names(newmap)[length(newmap)] <- cname
      }
    }
  }
  #process a possible column use or omit:
  fc <- which (captions %in% colnames$allowedname[colnames$prefname=="use"])
  if (length(fc) > 1) stop("only one column 'use' or 'omit' allowed")
  usemarker <- rep(TRUE, nrow(map))
  if (length(fc) == 1) {
    if (sum(is.na(map[, fc])) > 0 ||
        sum(map[, fc] %in% c(0, 1)) != nrow(map))
      stop(paste("invalid entries in column", captions[fc]))
    if (captions[fc] %in% c("omit", "exclude")) {
      usemarker <- map[, fc] == 0
    } else {
      usemarker <- map[, fc] == 1
    }
    #DO NOT delete the omitted markers:
    #newmap <- newmap[usemarker,]
  }
  list(map=newmap, usemarker=usemarker)
} #read_map_from_table

checkHaploblocks <- function(map) {
  #returns a character vector of haploblocks that do not form one consecutive
  #block of markers, or that group markers over chromosome boundaries
  result <- character(0)
  fp <- as.character(map$hb)
  ufp <- unique(fp)
  for (f in ufp) {
    wf <- which(fp == f)
    if (max(wf) - min(wf) + 1 != length(wf) ||          #not one consec. block
        length(unique(map$chrom[min(wf):max(wf)])) > 1) #not on one chrom
        result <- c(result, f)
  }
  result
} #checkHaploblocks

checkMapOrder <- function(map, allowColocalizedMarkers=TRUE) {
  #map: data frame with at least columns marker, chrom, pos
  #allowColocalizedMarkers: TRUE if successive markers may be at exactly the
  #         same position within a chromosome, FALSE if not
  #(decreasing order of position within a chromosome is never allowed)
  #return value:
  #character vector with in each element the names of two successive
  #markers violating the requirements, separated by a tab.
  #(0 elements if order was ok)
  if (nrow(map) < 1) stop("empty map")
  msg <- character(0)
  if (nrow(map) == 1) return(msg)
  if (allowColocalizedMarkers) {
    x <- which(map$chrom[2:nrow(map)] == map$chrom[1:(nrow(map)-1)] &
                 map$pos[2:nrow(map)] < map$pos[1:(nrow(map)-1)])
  } else {
    x <- which(map$chrom[2:nrow(map)] == map$chrom[1:(nrow(map)-1)] &
                 map$pos[2:nrow(map)] <= map$pos[1:(nrow(map)-1)])
  }
  if (length(x) > 0) {
    markerpairs <- cbind(map$marker[x], map$marker[x+1])
    msg <- paste(map$marker[x], "\t", map$marker[x+1])
  }
  msg
} #checkMapOrder

read_phasedgenofile <- function(map, ped, min.allele.freq=3,
                             phasedgenofile,
                             diagnostic_filename="marker_polymorphism.dat",
                             fq) {
  #map: a data.frame with columns marker, chrom, pos, fp; as returned by
  #     any of the read_map_xxx functions
  #ped: pedigree, a data frame as produced by read_pedigree
  #min.allele.freq: if there are not at least two alleles occurring at least
  #                 min.allele.freq times, usemarker (see below) is set to FALSE
  #phasedgenofile: input file, either in the flexqtl mhaplotypes.csv format
  #                (if fq TRUE) or in generic format (fq FALSE). The generic
  #                format is a tab-separated text file with empty cells for
  #                missing data. It has one row per individual and 2 columns per
  #                marker, with a header row where the first column of each marker
  #                is supposed to have the marker name (the name of the second
  #                marker column is ignored) and the individual names in the
  #                first column.
  #fq: TRUE (FlexQTL input) or FALSE (generic input)
  #return value: a list with items
  # $mrkallele.arr: a 3-dim array; dim 1 = markers, dim 2 = individuals,
  #          dim 3 = 1:2 (parental origin); contains the allele number:
  #          the allele numbers refer to allele names per marker
  #          as specified in the second item:
  # $alleledata: a data frame with columns markername, markernr,
  #              chrom, pos, haploblock, allelename, allelenr,
  #              one line per allele over all markers
  # $usemarker: logical vector in order of map: TRUE for markers to be used,
  #             FALSE for markers to be omitted from analysis (TRUE if at least
  #             2 alleles occur at least min.allele.freq times)
  #side effect: a file is written to diagnostic_filename with info on the
  #             number of observations of each marker allele and with usemarker
  messages <- character(0)
  if (fq) {
    #read flexqtl input:
    mhaplo <- read.table.origheaders(phasedgenofile, skip=20, header=TRUE,
                                     na.strings="-", sep=",", comment.char="",
                                     as.is=TRUE)
    #  no conversion to factors as each column has a mixture of data types
    indcount <- nrow(mhaplo)/5 #as read from mhaplotypes, may be different from pedigree
    # get the haplotypes from the first 2 lines for each individual:
    mhaplo <- mhaplo[sort(c(seq(1, by=5, length.out=indcount),
                            seq(2, by=5, length.out=indcount))),
                     -2] #omit 2nd column VAR (hap1/2, gpo1/2, prob)

  } else {
    #fq FALSE, read generic input:
    mhaplo <- read.table.origheaders(phasedgenofile, header=TRUE, sep="\t",
                                     comment.char="", na.strings=c("", "NA"),
                                     as.is=TRUE)
    #extra checks, since unknown source:
    if (length(mhaplo) < 3 || length(mhaplo) %% 2 == 0)
      return("read_phasedgenofile: invalid column count")
    indcount <- nrow(mhaplo)
    #any column that has only T and F will be interpreted as logical, correct:
    logi <- integer(0)
    for (i in 2:length(mhaplo))
      if (class(mhaplo[,i])=="logical") logi <- c(logi, i)
    if (length(logi) > 0) {
      mh <- read.table(phasedgenofile, skip=1, header=FALSE, sep="\t",
                       comment.char="", colClasses="character")
      for (i in logi) mhaplo[,i] <- mh[,i]
    }
    #now we convert the mhaplo to have 2 columns (instead of 1) for each markers
    #and 2 rows (instead of 1) for each individual:
    mhaplo <- twocols2tworows.dfr(mhaplo)
  } #fq FALSE
  #all markers on map must be in file, not v.v:
  mmrk <- which(names(mhaplo)[-1] %in% map$marker)
  if (length(mmrk) != nrow(map))
      stop("read_phasedgenofile: not all markers on map in file")
  mmrk <- match(map$marker, names(mhaplo)[-1])
  mhaplo <- mhaplo[,c(1, 1+mmrk)] #marker columns in map order
  mrkcount <- length(mmrk)
  #for each column (each marker): if all alleles are actually numeric,
  #convert them to numeric
  for (i in 2:length(mhaplo)) {
    suppressWarnings({
      num <- as.numeric(mhaplo[,i])
      if (sum(is.na(num)) == sum(is.na(mhaplo[,i]))) mhaplo[,i] <- num
    })
  }
  indnames <- mhaplo[seq(1, by=2, length.out=indcount), 1]
  # now each column has all allele names of that marker
  # For efficiency we convert all allele names to integers and keep a
  # translation table, and arrange the data as a 3-dim array:
  haparr <- array(integer(indcount*mrkcount*2), dim=c(mrkcount,indcount,2),
                  dimnames=list(names(mhaplo)[-1], #markers
                                indnames, #individuals
                                c("hap1","hap2"))) # parental origin
  names(dimnames(haparr)) <- c("marker", "individual", "parentalhaplo")
  alleledata <- data.frame()
  for (m in 1:mrkcount) {
    ch <- mhaplo[,m+1]
    chf <- as.factor(ch)
    chi <- as.integer(chf) #levels are now 1,2,3... in numeric or alphabetical order of alleles
    haparr[m,,1] <- chi[seq(1, by=2, length.out=indcount)]
    haparr[m,,2] <- chi[seq(2, by=2, length.out=indcount)]
    mname <- names(mhaplo)[m+1]
    df <- data.frame(markername=rep(mname,length(levels(chf))),
                     markernr=rep(m,length(levels(chf))),
                     chrom=rep(map$chrom[m],length(levels(chf))),
                     pos=rep(map$pos[m],length(levels(chf))),
                     haploblock=rep(map$hb[m],length(levels(chf))),
                     allelename=levels(chf),
                     allelenr=seq_along(levels(chf)))
    alleledata <- rbind(alleledata,df)
  }
  alleledata$allelename <- as.character(alleledata$allelename)
  #all individuals in pedigree must be in file, not v.v:
  pind <- which(dimnames(haparr)[[2]] %in% ped$name)
  if (length(pind) != nrow(ped))
      stop("read_phasedgenofile: not all individuals from pedigree in file")
  mhaplo_not_ped <- setdiff(dimnames(haparr)[[2]], as.character(ped$name))
  if (length(mhaplo_not_ped) > 0) {
    messages <- c(messages,
                  paste("read_phasedgenofile: the following",
                        length(mhaplo_not_ped),
                        "individuals occur in", phasedgenofile,
                        "but not in pedigree;"))
    messages <- c(messages,
                 "they will be ignored and won't occur in the output files:")
    for (i in seq_along(mhaplo_not_ped))
      messages <- c(messages, mhaplo_not_ped[i])
  }
  for (m in seq_along(messages)) cat(paste(messages[m], "\n", sep=""))
  #individuals in pedigree order:
  haparr <- haparr[, match(ped$name,dimnames(haparr)[[2]]),]
  #get the tables of allele counts and the and write the file
  allfrq <- list()
  usemarker <- rep(TRUE, nrow(map))
  max_allele_count <- 0
  for (mrk in seq_along(map$marker)) {
    allfrq[[mrk]] <- table(haparr[mrk,,], useNA="no")
    usemarker[mrk] <- sum(allfrq[[mrk]] >= min.allele.freq) > 1
    #   at least 2 alleles should occur at least min.allele.freq times, that
    #   is more severe than just requiring that the marker is not monomorphic
    #   if min.allele.freq > 1
    if (length(allfrq[[mrk]]) > max_allele_count)
      max_allele_count <- length(allfrq[[mrk]])
  }
  con <- file(diagnostic_filename, "w")
  s <- "marker\tchrom\tpos\thaploblock\tusemarker\tcount_NA"
  for (i in 1:max_allele_count) s <- paste(s,"\tcount_",i,sep="")
  for (i in 1:max_allele_count) s <- paste(s,"\tallele_",i,sep="")
  writeLines(s, con)
  for (mrk in seq_along(map$marker)) {
    s <- paste(map$marker[mrk], map$chrom[mrk], map$pos[mrk], map$hb[mrk],
               as.integer(usemarker[mrk]),
               2*indcount - sum(allfrq[[mrk]]), sep="\t")
    with0s <- c(allfrq[[mrk]], rep(0, max_allele_count - length(allfrq[[mrk]])))
    s <- paste(s, "\t", paste(with0s, sep="", collapse="\t"), sep="")
    with0s <- c(subset(alleledata, markernr==mrk)$allelename,
                rep("", max_allele_count - length(allfrq[[mrk]])))
    s <- paste(s, "\t", paste(with0s, sep="", collapse="\t"), sep="")
    writeLines(s, con)
  }
  close(con)
  list(mrkallele.arr=haparr, alleledata=alleledata, usemarker=usemarker,
       messages=messages)
} #read_phasedgenofile

mhaplotypes2generic <- function(mhaplotypes, generic) {
  #convert the phased genotypes data from flexqtl mhaplotypes.csv format
  #to generic phased genotyp file format
  #mhaplotypes: name of flexqtl file (input)
  #generic: name of new output file (will overwrite existing)
  con <- file(mhaplotypes, "r", blocking = FALSE)
  headers <- readLines(con, 21)
  headers <- strsplit(headers[21], split=",")[[1]]
  close(con)
  #to avoid the unwanted conversion of not allowed characters in marker names
  mhaplo <- read.table(mhaplotypes, skip=20, header=TRUE, na.strings="-",
                       sep=",", comment.char="", as.is=TRUE)
  indcount <- nrow(mhaplo)/5 #as read from mhaplotypes, may be different from pedigree
  # get the haplotypes from the first 2 lines for each individual:
  mhaplo <- mhaplo[sort(c(seq(1, by=5, length.out=indcount),
                          seq(2, by=5, length.out=indcount))),
                   -2] #omit 2nd column VAR (hap1/2, gpo1/2, prob)
  indnames <- mhaplo[,1]
  mhaplo <- data.frame(marker=headers[3:length(headers)],t(mhaplo[,-1]))
  names(mhaplo) <- c("marker", indnames)
  write.table(mhaplo, generic, na="", col.names=TRUE, row.names=FALSE,
              quote=FALSE, sep="\t")
} #mhaplotypes2generic

# function match.haplo returns TRUE if both haplotypes matching (i.e. no
# conflicts at any marker), FALSE if at least one conflict, and a specified
# value (T, F or NA) if at all markers at least one haplotype has NA
# (also if for all markers at least one of the hap is NA)
# hap1, hap2: 2 vectors with marker alleles of equal length > 0 (may be vectors
#             of allele numbers or names)
# match.NA: logical (TRUE or FALsE, not NA); if FALSE, NA in one hap matches NA
#           or any integer in the other; if TRUE, NA matches only NA
#           (i.e. only 2 identical haplotypes return TRUE)
# no.info: logical (TRUE, FALSE or NA): if match.NA is FALSE, the no.info value
#          is returned if at all markers at least one haplotype has NA
match.haplo <- function(hap1, hap2, match.NA, no.info) {
  if (length(hap1)==0 || length(hap2)==0 || length(hap1) != length(hap2)) {
    stop("match.haplo: both haplotypes must have equal length > 0")
  } else {
    neq <- hap1 != hap2
    one.na  <- xor(is.na(hap1), is.na(hap2))
    if (match.NA) {
      return(sum(one.na)==0 &&
        (is.na(sum(neq, na.rm=TRUE)) || sum(neq, na.rm=TRUE)==0))
    } else {
      if (sum(!is.na(neq))==0) {
        return(no.info) #at all markers at least one of the haplotypes has NA
      } else {
        return( sum(neq, na.rm=TRUE)==0)
      }
    }
  }
} #match.haplo

read_oldhballeles <- function(fname) {
  #reads a file with the alleles of all haploblocks from an
  #earlier run, with the aim of assigning the same ID's to the same alleles
  #in the current run.

  #Return value: a list with one item for each haploblock.
  #The items have the names of the haploblock and are character or numeric
  #arrays with markers in columns and haplotypes in rows; the row names are
  #the haplotype ID's (integers!) and the column names are either
  #the marker names (if available from the file) or NULL. The cell contents are
  #the marker allele names or NA.
  oldhballeles <- list()
  hapdata <- read.table(fname, header=FALSE, sep="\t", as.is=TRUE, na.strings="-")
  if (sum(hapdata[1,]=="hballelenr", na.rm=TRUE) != 1 || #this is the ID nr, excl a mv representation
        sum(hapdata[1,]=="markercount", na.rm=TRUE) != 1)
    stop("read_oldhballeles: file invalid")
  startlines <- which(hapdata[,1] == "hbnr")
  #check if there is at least one allele per haploblock (should be true unless
  #file is edited by user) and other errors:
  if (length(startlines) == 0 || startlines[1] != 1 || sum(is.na(startlines)) > 0 ||
        startlines[length(startlines)] == nrow(hapdata))
    stop("read_oldhballeles: file invalid")
  sl <- c(startlines, nrow(hapdata)+1)
  if (min(diff(sl)) <= 0)
    stop("read_oldhballeles: file invalid") #some haploblocks with no alleles
  fpnames <- hapdata[startlines+1, 2]
  markercounts <- as.integer(hapdata[startlines+1, which(hapdata[1,]=="markercount")])

  #hapdata <- hapdata[, c(4, 7:length(hapdata))] #hballeleID and all marker allele columns
  startlines <- c(startlines, nrow(hapdata)+1)
  for (fp in 1:(length(startlines)-1)) {
    hbcol <- which(hapdata[startlines[fp],] != "") #omit the final columns where this haploblock has no markers
    hd <-
      hapdata[(startlines[fp]+1):(startlines[fp+1]-1), hbcol]
    names(hd) <- hapdata[startlines[fp], hbcol]
    hapnrs <- gsub("^\\s+|\\s+$", "", hd[, "hballelenr"]) #whitespace trimmed
    if (sum(is.na(hapnrs)) > 0 || sum(hapnrs == "") > 0)
      stop("read_oldhballeles: missing or empty hballelenr")
    suppressWarnings(
      if (sum(is.na(as.numeric(hapnrs))) > 0)
        stop("read_oldhballeles: all hballelenr must be integers >= 1")
    )
    hapnrs <- as.numeric(hapnrs)
    if (sum(hapnrs < 1) > 0)
      stop("read_oldhballeles: all hballelenr must be integers >= 1")
    if (sum(hapnrs != round(hapnrs)) > 0)
      stop("read_oldhballeles: all hballelenr must be integers >= 1")
    if (length(unique(hapnrs)) != length(hapnrs))
      stop(paste("read_oldhballeles: duplicate hballelenr in haploblock",
                 fpnames[fp]))
    oldhballeles[[fp]] <- as.matrix(hd[, (length(hd)-markercounts[fp]+1):length(hd)])
    rownames(oldhballeles[[fp]]) <- hapnrs
    #colnames(oldhballeles[[fp]]) <- hd[1, -1]
  }
  names(oldhballeles) <- fpnames
  oldhballeles
} #read_oldhballeles

include_oldhballeles <- function(hapmat, fpname, oldhballeles, alleledata,
                                 emptyhapID) {
  #hapmat: a matrix with markers in columns and haplotype alleles in rows;
  #        values are integer (index to alleledata to obtain marker allele names)
  #oldhballeles: a list of character or numeric matrices as returned
  #            by read_oldhballeles, containing the old alleles for all old
  #            haploblocks from an earlier run
  #alleledata: a data frame with at least columns marker, allelenr, allelename;
  #            allelename should be numeric of character, not factor
  #emptyhapID: the ID (an integer) if the allele with only missing marker data
  #Return value: a list with 3 elements:
  #$hapmat: the original hapmat extended with new alleles from oldhballeles
  #         and the names of already present alleles overwritten by those from
  #         oldhballeles if those alleles are in there
  #$alleledata: the original alleledata extended with any extra marker alleles
  #             occurring in oldhballeles (needed to later write the haplotypes
  #             to file with those marker alleles that do not occur in the
  #             present population)
  #$messages: a character vector, empty if no messages
  messages <- character(0)
  if (is.null(oldhballeles) || length(oldhballeles) == 0 || is.null(alleledata))
    return(list(hapmat=hapmat, alleledata=alleledata, messages=messages))
  fp <- which(names(oldhballeles) == fpname)
  if (length(fp) > 1) {
    messages <- c(messages,
                  paste("include_oldhballeles: haploblock", fpname,
                        "occurs more than once; all ignored"))
    return(list(hapmat=hapmat, alleledata=alleledata, messages=messages))
  }
  if (length(fp) == 0)
    #this focalpoint did not occur in the earlier run
    return(list(hapmat=hapmat, alleledata=alleledata, messages=messages))
  oldhballeles <- oldhballeles[[fp]]
  if (nrow(oldhballeles) == 0)
    return(list(hapmat=hapmat, alleledata=alleledata, messages=messages))
  #are the markers in alleledata the same as in hapmat?
  fpmrk <- unique(as.character(alleledata$markername[alleledata$haploblock==fpname]))
  if (!setequal(fpmrk, colnames(hapmat)))
    stop(paste("include_oldhballeles: markers in alleledata don't match hapmat in haploblock",
               fpname))
  #are the markers in oldhballeles the same as in hapmat?
  if (ncol(oldhballeles) != ncol(hapmat))
    stop(paste("include_oldhballeles: markers in oldhballeles don't match in haploblock",
               fpname))
  if (!setequal(colnames(oldhballeles), colnames(hapmat)))
    stop(paste("include_oldhballeles: markers in oldhballeles don't match in haploblock",
               fpname))
  #in case the markers within the focalpoint are rearranged, align them:
  if (sum(colnames(oldhballeles) != colnames(hapmat)) > 0) {
    messages <- c(messages, paste("in oldhballeles markers in different order for haploblock",
                                  fpname, "; rearranged"))
    oldhballeles <- oldhballeles[,match(colnames(oldhballeles), colnames(hapmat))]
  }
  #if we reach this, the markers in all three data sources match.

  #if hballele 1 is present in oldhballeles it should consist of only
  #missing marker data
  empty <- which(rownames(oldhballeles) == emptyhapID)
  if (length(empty) == 1 && sum(!is.na(oldhballeles[empty,])) > 0)
    stop(paste("include_oldhballeles: markers in oldhballele allele", emptyhapID,
               "are not all missing in haploblock", fpname))
  #if any haplotype with only missing data occurs in oldhaballeles it must have
  #the emptyhapID:
  empty <- which(rowSums(!is.na(oldhballeles)) == 0)
  if (length(empty) > 0 && sum(rownames(oldhballeles[empty]) != emptyhapID) > 0)
    stop(paste("include_oldhballeles: the allele in oldhballele with only missing marker values",
               "do not have the hballelenr", emptyhapID,
               "in haploblock", fpname))
  #first add any marker alleles from the previous run to alleledata:
  mcount <- ncol(hapmat)
  for (m in 1:mcount) {
    mrkall <- alleledata[alleledata$markername == colnames(hapmat)[m],]
    oldall <- setdiff(unique(oldhballeles[,m]), mrkall$allelename)
    oldall <- oldall[!is.na(oldall)]
    if (length(oldall) > 0) {
      oldmaxnr <- max(mrkall$allelenr)
      extra <- mrkall[rep(1, length(oldall)),]
      extra$allelename <- oldall
      extra$allelenr <- (oldmaxnr+1):(oldmaxnr+length(oldall))
      alleledata <- rbind(alleledata, extra)
    }
  }

  #next create an array corresponding to oldhballeles with the marker allelenrs
  #instead of the marker allele names:
  hapallelenrs <- matrix(as.integer(NA),
                         nrow=nrow(oldhballeles), ncol=ncol(oldhballeles))
  rownames(hapallelenrs) <- rownames(oldhballeles)
  colnames(hapallelenrs) <- colnames(oldhballeles)
  for (m in 1:mcount) {
    mrkall <- alleledata[alleledata$markername == colnames(hapmat)[m],]
    hapallelenrs[,m] <- mrkall$allelenr[match(oldhballeles[,m], mrkall$allelename)]
  }

  #add the alleles from hapallelenrs sequentially to hapmat, checking if
  #they already existed (this will also remove duplicates within hapallelenrs):
  for (a in seq_len(nrow(hapallelenrs))) {
    b <- nrow(hapmat)
    while (b > 0 && !identical(hapallelenrs[a,], hapmat[b,])) b <- b - 1
    if (b > 0) {
      #allele already in hapmat, assign name from hapallelenrs:
      rownames(hapmat)[b] <- rownames(hapallelenrs)[a]
    } else {
      #this allele did not yet exist and must be added:
      hapmat <- rbind(hapmat, hapallelenrs[a,, drop=FALSE])
    }
  }
  list(hapmat=hapmat, alleledata=alleledata, messages=messages)
} #include_oldhballeles

getnewhapnr <- function(hapmat, index=nrow(hapmat)) {
  #hapmat: a matrix with markers in columns and haploblock alleles in rows;
  #        values are integer (index to alleledata to obtain marker allele names)
  #index: the row in hapmat for which a name must be obtained
  #Currently a new allele gets the next available integer and the numbers are equal to the row numbers
  currnames <- as.integer(rownames(hapmat)[-index])
  if (sum(is.na(currnames)) > 0) stop("getnewhapnr: non-integer numbers occur")
  max(currnames) + 1
} #getnewhapnr

get_initial_haplotypes <- function(mrkallele.arr, map, usemarker,
                                   alleledata,
                                   oldhballeles=NULL,
                                   emptyhapID=1) {
  #obtain the haplotypes for each focal point (with NA as one extra allele)
  #mrkallele.arr: 3-dim array of marker alleles per marker/individual/parent
  #               as returned by read_mhaplotypes
  #map: a data.frame with columns marker, chrom, pos, fp; as returned by
  #     read_map_from_mhaplotypes; used for grouping of markers to focal points
  #usemarker: a logical vector without NAs of length nrow(map), indicating for
  #           each marker if it must be used or not
  #alleledata: as returned my read_mhaplotypes; data frame with at least
  #            columns markername, allelename, allelenr
  #oldhballeles: a data frame with one row per haplotype (as assigned
  #            in a previous run) with at least columns focalpoint (the name),
  #            haploname (name), and N columns with captions marker_1 .. marker_N
  #            containing the alleles (names), "-" for missing alleles or
  #            "%" for unused markers. oldhballeles may not contain duplicate
  #            alleles for the same focalpoint; the number of non-"%" markers
  #            must match the number of markers in the present focalpoint.
  #            The markers in oldhballeles must be the same markers in the same
  #            order as in the current map, but apart from the number of
  #            markers that cannot be checked
  #emptyhapID: the default name of the first haplotype of each focalpoint,
  #            which contains missing data for all markers. This name may be
  #            overwritten if an allele with only missing data is present in
  #            oldhballeles
  #return value: a list with items
  # $fp_haplotypes: a list with one item per focal point; for each focal point
  #                 the following items:
  #    $chrom: character or numeric value, name of the chromosome of this
  #            focalpoint
  #    $pos: numeric value, average position of all markers of this focalpoint
  #    $markers: integer vector of marker numbers in the focal point (numbering
  #              in map order; order of markers within focal point not relevant)
  #    $na.count: integer vector: for each haplotype in this focal point
  #               the number of missing marker alleles;
  #               a haplotype is a sequence of marker alleles, where NA is
  #               considered an allele (e.g. A-A-B, A-B-B and A-NA-B are three
  #               different haplotypes)
  #    $hapmat: integer matrix with markers in columns and haplotypes in rows,
  #             with the marker allele numbers (corresponding to the alleledata
  #             item in the return value of read_mhaplotypes).
  #             colnames are marker names, rownames are haplotype names (by
  #             default consecutive numbers)
  # second item of return value:
  # $hapnrarr: a 3-dim array with per focalpoint / individual / parent the
  #            haplotype number, which in an index into the
  #            fp_haplotypes[[fp.ix]] sub-items
  # third item of return value:
  # $alleledata: the input data frame with additional marker alleles from
  #              oldhballeles
  # $messages: a character vecor

  #first some input checks:
  if (dim(mrkallele.arr)[1] != nrow(map) ||
      length(which(dimnames(mrkallele.arr)[[1]] != map$marker)) > 0) {
    stop("get_initial_haplotypes: markers in mrkallele.arr don't match map")
  }
  if (length(usemarker) == 1) {
    #allow to specify usemarker=TRUE to include all markers
    usemarker <- rep(usemarker, nrow(map))
  }
  if (length(usemarker) != nrow(map) || sum(is.na(usemarker)) > 0) {
    stop("get_initial_haplotypes: usemarker doesn't match map or contains NAs")
  }
  if (sum(usemarker) == 0) {
    stop("get_initial_haplotypes: usemarker excludes all markers")
  }

  #first get the list of chromosome names in supplied map order
  #(not alphabetical / numerical order)
  chromnames <- c(as.character(map$chrom),"")
  chromlist2 <- c("",as.character(map$chrom))
  chromnames <- chromnames[chromnames!=chromlist2]
  length(chromnames) <- length(chromnames)-1
  #next get the list of focalpoints in map order:
  fpnames <- unique(as.character(map$hb))
  empty_fp <- logical(length(fpnames))

  # create fp.df: for each fp its name, chrom, average position)
  fp.df <- data.frame()
  for (fp.ix in seq_along(fpnames)) {
    fpname <- fpnames[fp.ix]
    fpchrom <- unique(map$chrom[map$hb == fpname])
    if (length(fpchrom) != 1)
      stop("get_initial_haplotypes: each focalpoint should be on one chromosome")
    fppos <- round(mean(map$pos[usemarker & map$hb == fpname]), digits=3)
    fp.df <- rbind(fp.df, data.frame(fp=fpname, chrom=fpchrom, pos=fppos))
  }
  fp.df <- fp.df[!is.na(fp.df$pos),] #excludes fps without used markers
  fp.df <- fp.df[order(match(fp.df$chrom, chromnames), fp.df$pos),] #fp.df in map order
  #   NOTE: if this results in a reordering of the fp's, in a later stage the
  #         marker order will be different from the original order, causing
  #         problems in haplo2mrkallnrarr. Therefore we now (20150609) require in
  #         fq_haplotyping_session that the markers are already in map order
  fp.df$chrom <- as.character(fp.df$chrom)
  suppressWarnings(
    if (sum(is.na(as.numeric(fp.df$chrom)))==0) {
      fp.df$chrom <- as.numeric(fp.df$chrom)
    })
  fp.df$fp <- as.character(fp.df$fp)
  suppressWarnings(
    if (sum(is.na(as.numeric(fp.df$fp)))==0) {
    fp.df$fp <- as.numeric(fp.df$fp)
  })
  #next we create and fill fp_haplotypes and hapnrarr
  fpcount <- nrow(fp.df)
  indcount <- dim(mrkallele.arr)[2]
  hapnrarr <- array(integer(fpcount*indcount*2), dim=c(fpcount,indcount,2),
                    dimnames=list(fp.df$fp,                      #focalpoints
                                  dimnames(mrkallele.arr)[[2]],  #individuals
                                  dimnames(mrkallele.arr)[[3]])) #hap1/hap2
  names(dimnames(hapnrarr)) <- c("focalpoint", "individual", "parentalhaplo")

  fp_haplotypes <- list()
  messages <- character(0)
  for (fp.ix in seq_along(fp.df$fp)) {
    fp_haplotypes[[fp.ix]] <- list()
    fp_haplotypes[[fp.ix]]$chrom <- fp.df$chrom[fp.ix]
    fp_haplotypes[[fp.ix]]$pos <- fp.df$pos[fp.ix]
    fp_haplotypes[[fp.ix]]$markers <-
      which(usemarker & map$hb == fp.df$fp[fp.ix]) #index to marker in map
    fp_haplotypes[[fp.ix]]$na.count <- integer(0) #put this before hapmat
    # make an array of haplotypes including NAs for this fp:
    # first row: the haplotype with all marker alleles missing (even if that
    #            doesn't occur in the data)
    fp_haplotypes[[fp.ix]]$hapmat <-
      matrix(as.integer(rep(NA,length(fp_haplotypes[[fp.ix]]$markers))), nrow=1)
    rownames(fp_haplotypes[[fp.ix]]$hapmat) <- emptyhapID
    colnames(fp_haplotypes[[fp.ix]]$hapmat) <-
      map$marker[fp_haplotypes[[fp.ix]]$markers]
    #if earlier (named) haplotypes are given, include them:
    if (!is.null(oldhballeles)) {
      tmp <-
        include_oldhballeles(hapmat=fp_haplotypes[[fp.ix]]$hapmat,
                           fpname=fp.df$fp[fp.ix],
                           oldhballeles=oldhballeles, alleledata=alleledata,
                           emptyhapID=emptyhapID)
      fp_haplotypes[[fp.ix]]$hapmat <- tmp$hapmat
      alleledata <- tmp$alleledata
      messages <- c(messages, tmp$messages)
    }
    #now add all haplotypes from all individuals:
    for (ind in 1:indcount) for (p in 1:2) {
      ihap <- mrkallele.arr[fp_haplotypes[[fp.ix]]$markers,ind,p]
      h <- 1
      while (h<=nrow(fp_haplotypes[[fp.ix]]$hapmat) &&
             !match.haplo(fp_haplotypes[[fp.ix]]$hapmat[h,],
                          ihap, match.NA=TRUE,no.info=TRUE)) {
        h <- h+1
      }
      if (h>nrow(fp_haplotypes[[fp.ix]]$hapmat)) {
        # new haplotype found
        fp_haplotypes[[fp.ix]]$hapmat <-
          rbind(fp_haplotypes[[fp.ix]]$hapmat,ihap)
        rownames(fp_haplotypes[[fp.ix]]$hapmat)[h] <-
          getnewhapnr(fp_haplotypes[[fp.ix]]$hapmat, h)
        #fp_haplotypes[[fp.ix]]$na.count[h] <- sum(is.na(ihap))
      }
      fp_haplotypes[[fp.ix]]$na.count <-
        rowSums(is.na(fp_haplotypes[[fp.ix]]$hapmat))
      hapnrarr[fp.ix,ind,p] <- h
    }
    #rownames(fp_haplotypes[[fp.ix]]$hapmat) <- NULL
    # now we have the complete set of haplotypes_including_NAs in
    # fp_haplotypes[[fp.ix]]$hapmat, and an array hapnrarr has the haplotype
    # numbers corresponding to the rows of that matrix for each focalpoint,
    # individual and hap1/hap2
    msg <- paste("hb ", fp.ix, " = ", fp.df$fp[fp.ix],": ",
                 nrow(fp_haplotypes[[fp.ix]]$hapmat)," initial alleles",
                 sep="")
    messages <- c(messages, msg)
    cat(paste(msg, "\n", sep=""))
  }
  names(fp_haplotypes) <- fp.df$fp #fp_haplotypes is sorted in map order

  list(fp_haplotypes=fp_haplotypes, hapnrarr=hapnrarr,
       alleledata=alleledata, messages=messages)
} # get_initial_haplotypes

get_hballeleID <- function(haploblock, hballele, fp_haplotypes, mv.sep) {
  #haploblock: the index to a haploblock in fp_haplotypes
  #hballele: a vector of *indices* of an allele of the haploblock (i.e. their
  #          row numbers in the $hapmat matrix)
  #fp_haplotypes: a fp_haplotypes object as returned by get_initial_haplotypes
  #mv.sep: if NA the hballeleID is just its number (the row NAME of the $hapmat
  #        matrix);
  #        if a string of length 1 that will be inserted between its number and
  #        the number of missing marker data in the allele;
  #        if a string of length 2, eg "()", these will enclose the number of
  #        missing marker data
  if (haploblock < 1 || haploblock > length(fp_haplotypes) ||
        sum(is.na(hballele)) > 0 ||
        min(hballele) < 1 ||
        max(hballele) > nrow(fp_haplotypes[[haploblock]]$hapmat)) {
    stop("get_hballeleID: invalid parameter values")
  }
  if (is.na(mv.sep)) {
    rownames(fp_haplotypes[[haploblock]]$hapmat)[hballele]
  } else {
    paste(rownames(fp_haplotypes[[haploblock]]$hapmat)[hballele],
          substr(mv.sep, 1, 1),
          #formatC(fp_haplotypes[[haploblock]]$na.count[hballele],
          #width = nchar(length(fp_haplotypes[[haploblock]]$markers)),
          #format = "d", flag = "0"), #padded with zeroes
          fp_haplotypes[[haploblock]]$na.count[hballele], #unpadded
          substr(mv.sep, 2, 2), #"" if nchar(mv.sep)==1
          sep="")
  }
} #get_hballeleID

conv.hballelenames <- function(x, fp_haplotypes, mv.sep) {
  #x: a 3D array as hapnrarr with haploblocks in dim 1, individuals in dim 2
  #   and parental (1/2) in dim 3, containing hballele numbers
  #fp_haplotypes: a list as returned by get.initial.haplotypes
  #mv.sep: a character (pair) used to separate the hballele number from the
  #             number of missing marker scores (in get_hballeleID)
  result <- x
  #we cannot write back to x because after the first hb x is converted
  #to character instead of integer
  for (hb in seq_len(dim(x)[1])) {
      result[hb,,] <- matrix(get_hballeleID(hb, x[hb,,], fp_haplotypes, mv.sep),
                        ncol=2)
  }
  result
} #conv.hballelenames

write.all.haploblockalleles <- function(filename, fp_haplotypes, hapnrarr,
                                     alleledata=NULL, map=NULL, na.char="-",
                                     mv.sep) {
  #fp_haplotypes: a fp_haplotypes object as returned by get_initial_haplotypes
  #hapnrarr: 3-dim array of haplotype numbers (alleles) with dimensions
  #          focalpoint, individual, hap1/hap2, as returned by
  #          get_initial_haplotypes
  #alleledata: if NULL, allele numbers are written, else alleledata must be
  #            an alleledata object as returned by read_mhaplotypes, and then
  #            the marker allele names are written
  #map: data frame
  #na.char: currently works only for allele names, not if alleledata NULL
  maxfpmrk <- 0
  for (fp.ix in seq_along(fp_haplotypes))
    if (length(fp_haplotypes[[fp.ix]]$markers) > maxfpmrk)
      maxfpmrk <- length(fp_haplotypes[[fp.ix]]$markers)
  con <- file(filename, "w")
  for (fp.ix in seq_along(fp_haplotypes)) {
    if (is.null(map)) markerlist <- fp_haplotypes[[fp.ix]]$markers else
      markerlist <- map$marker[fp_haplotypes[[fp.ix]]$markers]
    line <- paste("hbnr\thaploblock\tmarkercount\thballeleID\thballelenr\tnacount\tfreq",
                  paste(markerlist, collapse="\t"), sep="\t")
    filler <- paste(rep("\t", maxfpmrk-length(markerlist)), collapse="")
    line <- paste(line, filler, sep="")
    writeLines(line, con)
#     writeLines(paste("fpnr\thaploblock\tmarkercount\thaplonr\thaploID\tnacount\tfreq",
#                      paste(markerlist, collapse="\t"),
#                      paste(rep("", maxfpmrk-length(markerlist)),
#                            collapse="\t"),
#                      sep="\t"),
#                con)
    #get freq of all haplotypes:
    hapfrq <- table(hapnrarr[fp.ix,,])
    hapfrq <- hapfrq[match(seq_along(fp_haplotypes[[fp.ix]]$na.count),
                           names(hapfrq))] #add the non-occurring haplotypes
    hapfrq[is.na(hapfrq)] <- 0
    startline <- paste(fp.ix,
                       names(fp_haplotypes)[fp.ix],
                       length(fp_haplotypes[[fp.ix]]$markers),
                       sep="\t")
    for (h in seq_along(fp_haplotypes[[fp.ix]]$na.count)) {
      midline <- paste(get_hballeleID(fp.ix, h, fp_haplotypes, mv.sep),
                       rownames(fp_haplotypes[[fp.ix]]$hapmat)[h],
                       fp_haplotypes[[fp.ix]]$na.count[h],
                       hapfrq[h],
                       sep="\t")
      if (is.null(alleledata)) {
        endline <- paste(fp_haplotypes[[fp.ix]]$hapmat[h,],
                         collapse="\t")
      } else {
        endline <- ""
        for (m in seq_along(fp_haplotypes[[fp.ix]]$markers)) {
          hapname <- subset(alleledata,
                            markernr==fp_haplotypes[[fp.ix]]$markers[m] &
                              allelenr==fp_haplotypes[[fp.ix]]$hapmat[h,m]
                           )$allelename
          if (length(hapname)==0) hapname <- na.char
          endline <- paste(endline, hapname, sep="\t")
        }
        endline <- paste(substring(endline,2,nchar(endline)),
                         filler,
                         sep="")
      }
      writeLines(paste(startline,midline,endline,sep="\t"),con)
    }
  }
  close(con)
} #write.all.haploblockalleles

#enumerate all HS families and sort by decreasing size:
get.all.HSfamilies <- function(ped) {
  #ped: pedigree, a data frame as produced by read_pedigree
  #return value: a list with one item per HS family,
  #sorted in order of decreasing family size.
  # items of each family:
  # $parentname
  # $parentnr: the individual number (in the ped) of the parent
  #            to be used for indexing hapnrarr etc
  # $HSf1: the individuals with parent as first parent
  # $HSf2: the individuals with parent as second parent
  # $size: the total size of the HS (HSf1 and HSf2)
  HSparents <- unique(c(as.character(ped$parent1),as.character(ped$parent2)))
  HSparents <- HSparents[!is.na(HSparents)]
  HSfam <- list()
  HSfamSize <- integer(0)
  for (hsp in seq_along(HSparents)) {
    HSf1 <- which(ped$parent1 == HSparents[hsp])
    HSf2 <- which(ped$parent2 == HSparents[hsp])
    # progeny from a selfing is among HSf1 and HSf2
    # that is correct, as both its parental alleles are derived from this parent
    if (length(c(HSf1,HSf2)) == 0) {
      stop("error in get.all.HSfamilies: family size = 0")
    }
    HSnr <- length(HSfam)+1
    HSfam[[HSnr]] <- list()
    HSfam[[HSnr]]$parentname <- HSparents[hsp]
    HSfam[[HSnr]]$parentnr <- which(ped$name == HSparents[hsp])
    HSfam[[HSnr]]$HSf1 <- HSf1 #the individuals with parent as first parent
    HSfam[[HSnr]]$HSf2 <- HSf2 #the individuals with parent as second parent
    HSfam[[HSnr]]$size <- length(c(HSf1,HSf2))
    HSfamSize[HSnr] <- HSfam[[HSnr]]$size
  }
  tmp <- order(HSfamSize, decreasing=TRUE)
  HSfam <- HSfam[tmp]
  HSfam
} #get.all.HSfamilies

getHSParentalHaplotypeFreqs <- function(HSnr, fphapnrarr, HSfam) {
  #HSnr: index to HSfam
  #HSfam: a list as produced by get.all.HSfamilies
  #fphapnrarr: 2-dim array of haplotype numbers (alleles) with dimensions
  #          individual, hap1/hap2: hapnrarr[fp,,] where hapnrarr is a 3-D array
  #          as returned by get_initial_haplotypes
  #return value: a frequency table of all haplotypes inherited from the parent
  #HSf1, HSf2: numbers of :
  HSf1 <- HSfam[[HSnr]]$HSf1 #HS individuals with parent as mother
  HSf2 <- HSfam[[HSnr]]$HSf2 #HS individuals with parent as father
  table(c(fphapnrarr[HSf1, 1], fphapnrarr[HSf2, 2]))
} #getHSParentalHaplotypeFreqs

# get.HSfam.haplotypes collects haplotype info for one haploblock and
# the individuals with parent as mother and as father combined
get.HSfam.haplotypes <- function(HSnr, fphapnrarr, HSfam) {
  #fp: focalpoint number (not ID)
  #HSnr: number of HSfamily (indexes HSfam, the result of get.all.HSfamilies)
  #fphapnrarr: 2-dim array of haplotype numbers (alleles) with dimensions
  #            individual, hap1/hap2: hapnrarr[fp,,] where hapnrarr is a 3-dim
  #            array as returned by get_initial_haplotypes
  #HSfam: a list as produced by get.all.HSfamilies
  #return value: a list with items:
  # $parhaplo: a vector with the 2 haplotypes of the parent
  # $HShaplofrq: a freq.table as returned by getHSParentalHaplotypeFreqs:
  #              a table with the frequencies of the haplotypes in the HS family
  #NOTE that the hballeles in fphapnrarr are the row numbers in fp_haplotypes$hapmat,
  #NOT the hballele names (the rownames of that hapmat) !!!
  result <- list()
  result$parhaplo <- fphapnrarr[HSfam[[HSnr]]$parentnr,]
  result$HShaplofrq <- getHSParentalHaplotypeFreqs(HSnr, fphapnrarr, HSfam)
  result
} #get.HSfam.haplotypes

writeHSfamHaplotypes.faminfo <- function(HSnr, HSfam) {
  #used a.o. by writeHSfamHaplotype, see there for parameters
  paste(HSnr,
        HSfam[[HSnr]]$size,
        HSfam[[HSnr]]$parentname,
        sep="\t")
} #writeHSfamHaplotypes.faminfo

writeHSfamHaplotypes.hapinfo <- function(fp, hap, hapinfo, fp_haplotypes, mv.sep) {
  #used a.o. by writeHSfamHaplotype, see there for further info
  #hap: haplotype number in current focalpoint
  #hapinfo: a list as produced by get.HSfam.haplotypes
  inparent <- (hap %in% hapinfo$parhaplo)
  paste(names(fp_haplotypes)[fp],
        get_hballeleID(fp, hapinfo$parhaplo[1], fp_haplotypes, mv.sep),
        get_hballeleID(fp, hapinfo$parhaplo[2], fp_haplotypes, mv.sep),
        get_hballeleID(fp, hap, fp_haplotypes, mv.sep),
        (0+inparent),
        hapinfo$HShaplofrq[which(names(hapinfo$HShaplofrq) == hap)],
        sep="\t")
} #writeHSfamHaplotypes.hapinfo

writeHSfamHaplotypes <- function(HSfam, fp_haplotypes, hapnrarr, filename,
                                 na.char="-", mv.sep) {
  #HSfam: a list as produced by get.all.HSfamilies
  #hapnrarr: 3-dim array of haplotype numbers (alleles) with dimensions
  #          focalpoint, individual, hap1/hap2, as returned by get_initial_haplotypes
  con <- file(filename, "w")
  writeLines("HSfam\tfamsize\tparent\thaploblock\tparhap1\tparhap2\thballeleID\tinparent\tfreq\tmarkeralleles",con)
  for (HSnr in seq_along(HSfam)) {
    startline <- writeHSfamHaplotypes.faminfo(HSnr, HSfam)
    for (fp in 1:dim(hapnrarr)[1]) {
      hapinfo <- get.HSfam.haplotypes(HSnr, hapnrarr[fp,,], HSfam)
      for (h in seq_along(hapinfo$HShaplofrq)) {
        hap <- as.integer(names(hapinfo$HShaplofrq)[h]) #the current haplotype
        midline <- writeHSfamHaplotypes.hapinfo(fp, hap, hapinfo, fp_haplotypes, mv.sep)
        endline <- fp_haplotypes[[fp]]$hapmat[hap,]
        endline[is.na(endline)] <- na.char
        endline <- paste(endline, sep="\t", collapse="\t")
        writeLines(paste(startline,midline,endline,sep="\t"),con)
      }
    } #for fp
  }
  close(con)
} #writeHSfamHaplotypes

combineHaplotypes <- function(haplotypes,
                              focalpointnr=NULL,
                              fp_haplotypes=NULL) {
  #haplotypes: integer vector of haplotype numbers to combine (if
  #             focalpointnr and fp_haplotypes are not NULL)
  #             or matrix with haplotype sequences in rows
  #focalpointnr: the focalpoint sequential number (not the ID), or NA
  #return value: the consensus haplotype sequence
  #              or NA if not all haplotypes match
  if (length(haplotypes)==0) return(NA)
  if (!is.matrix(haplotypes)) {
    if (is.null(focalpointnr)) return(NA) else {
      if (length(focalpointnr) != 1)
        stop("combineHaplotypes: length(focalpointnr) != 1")
      haplovec <- haplotypes
      haplotypes <- fp_haplotypes[[focalpointnr]]$hapmat[haplotypes[haplovec],, drop=FALSE]
    }
  }
  if (!is.matrix(haplotypes)) stop("combineHaplotypes: not matrix")
  if (nrow(haplotypes)==1) return(haplotypes[1,]) #as vector
  consensus <- haplotypes[1,, drop=FALSE]
  for (i in 2:nrow(haplotypes)) {
    hap <- haplotypes[i,]
    if (!match.haplo(hap,consensus, match.NA=FALSE, no.info=TRUE)) {
      return(NA)
    } else {
      consensus[is.na(consensus)] <- hap[is.na(consensus)]
    }
  }
  return(consensus)
} #combineHaplotypes

calc.matchmatrix <- function(haploseq, names=NULL) {
  # haploseq: matrix with one row per haplotype, with the marker data
  # names: names (numbers) of the haplotypes in haploseq, used as row and column
  #        names of  matchmatrix if given
  #return value: a square matrix with in rows and columns the haplotypes
  #              and in each cell TRUE or FALSE indicating that these 2
  #              haplotypes match each other or not, or NA if at all markers
  #              at least one of both has a NA score
  matchmatrix <- matrix(TRUE, nrow=nrow(haploseq), ncol=nrow(haploseq))
  for (h1 in 1:nrow(haploseq)) for (h2 in h1:nrow(haploseq)) {
    matchmatrix[h1,h2] <- match.haplo(haploseq[h1,], haploseq[h2,],
                          match.NA=FALSE, no.info=NA)
    matchmatrix[h2,h1] <- matchmatrix[h1,h2]
  }
  if (!is.null(names) && length(names) == nrow(haploseq)) {
    colnames(matchmatrix) <- rownames(matchmatrix) <- names
  }
  matchmatrix
} #calc.matchmatrix

groupHaplotypes <- function(haplonr, haploseq, haplofreq, matchmatrix=NULL) {
  # haplonr: integer vector of all haplotype IDs to be grouped, eg. all
  #          haplotypes occurring in one family (ordered by descending frequency)
  # haploseq: matrix with one row per haplonr, with the marker alleles
  # haplofreq: table with frequencies of haplonr in current family
  # matchmatrix: see return value. If NULL it is generated, else used and
  #              returned without checking
  # return value: a list with items
  # $matchmatrix: a square matrix with in rows and columns the haplotypes,
  #               and in each cell T or F indicating if these 2 haplotypes
  #               match each other or NA if at all markers at least one of both
  #               has a NA
  # $groups: a list of vectors each containing a set of haplotype nrs having true
  #          matches (not only absence of conflicts) within the group and
  #          conflicting with the consensus of each other group.
  #          the group list is sorted in order of decreasing group frequency
  # $ungrouped: a vector of haplotype nrs not in one of the groups (because they
  #             match with more than one group, or have only NA matches with all
  #             groups; NOT the $nodata haplotype)
  #             sorted in order of decreasing frequency
  # $nodata: integer(0) or one integer (1: the haplotype with NA at all markers)
  # $consensus: matrix with for each group a row with the consensus marker scores
  # $groupfreq: vector with the total frequencies of each groups, in same order
  #             as $groups
  # $maxgroupfreq: vector with the total frequencies of each group, assuming that
  #                all haplotypes in $ungrouped and $nodata that match (a.o.)
  #                with this group actually belong to this group
  # $ungroupedfreq: vector of frequencies of all haplotypes in $ungrouped
  # $log: (currently, for debugging): a character vector

  #1.After building the groups and removing one item that matches 2 or more groups
  #the remaining items might all belong to the same group (since the removed
  #item might have caused the groups to split). True? Yes, even though it
  #matched some items in the second group it might not have matched all items
  #in the second group; so removing it may change the whole group structure.
  #Therefore each time a haplotype is removed from one of the groups,
  #start over with building the groups and removing one.
  #When no more haplotypes removed: all remaining groups have only haplotypes
  #that match (no conflict, not just no info) all other haplotypes within the groups
  #and that conflict AT LEAST ONE haplotype in every other group.
  #Then we can try to add one by one the removed haplotypes. They can be added if
  #(a) they match all haplotypes in one group and (b) conflict with the consensus
  #in all other groups. We don't create any new groups in that final stage: new
  #groups might cause some haplotypes in other groups now to match the new
  #groups as well
  #
  #Note that within this function we mostly use INDICES to haplonr etc
  #to identify the haplotypes, not the haplotype IDs.
  #this is often shown in names (.hapix, hi = haplotype index)
  log <- "GroupHaplotypes"
  log <- c(log,paste("haplotypes:",paste(haplonr,collapse=" ")))
  if (is.vector(haploseq)) haploseq <- matrix(haploseq, nrow=1)
  result <- list()
  if (is.null(matchmatrix)) {
    result$matchmatrix <- calc.matchmatrix(haploseq, names=haplonr)
  } else {
    result$matchmatrix <- matchmatrix
  }
  nodata.hapix <- which(rowSums(!is.na(haploseq))==0) # the haplotype with all
  #                                                     marker scores unknown
  log <- c(log,paste("nodata.hapix:", paste(nodata.hapix, collapse=" ")))
  hap.nacount <- rowSums(is.na(haploseq))

  #build groups taking the haplotypes in order of increasing number of NA marker scores
  ungrouped.hapix <- integer(0) # haplotypes set (temporarily) aside because they
  #                               match more than 1 group

  #define a nested function assign.hapix that takes one haplotype index hi (into
  #haplonr) and returns the group nr where it belongs:
  #if hi conflicts with all existing groups, the return value is one number above
  #the highest group
  #if there are two or more groups with which hi does not conflict, NA is returned
  #else (just one group that does not conflict with hi): if hi has a true match
  #(other than through NA marker scores) with this group, the group number is
  #returned, else NA
  #Note that this function does not actually create a group or move the haplotype
  #to a group, it just says where the haplotype would fit.
  assign.hapix <- function(hi) {
    if (length(groups.hapix)==0) {
      return(1)
    }
    conflictgrp <- logical(length(groups.hapix)) #for each group, does hi conflict with it?
    matchgrp <- conflictgrp #for each group, does hi have a true match (other than NA)?
    for (gi in seq_along(groups.hapix)) {
      h <- 1 #index to haplotypes within group gi
      while(h<=length(groups.hapix[[gi]]) &&
            !conflictgrp[gi]) {
        if (!is.na(result$matchmatrix[hi,groups.hapix[[gi]][h]])) {
          if (result$matchmatrix[hi,groups.hapix[[gi]][h]]) {
            matchgrp[gi] <- TRUE
          } else {
            conflictgrp[gi] <- TRUE
          }
        }
        h <- h + 1
      } # while h
    }
    if (sum(!conflictgrp)==0) {
      #hi conflicts with all groups, should start new group:
      return(length(groups.hapix) + 1)
    }
    if (sum(!conflictgrp)==1) {
      #hi conflicts with all but one group:
      if (matchgrp[!conflictgrp]) {
        #actual match in only group not conflicting with hi:
        return(which(!conflictgrp))
      }
      #else no actual match:
      return(NA)
    }
    #else hi has no conflicts with more than one group:
    return(NA)
  } #assign.hapix within groupHaplotypes

  removed <- TRUE
  cycle <- 0
  while (removed) {
    removed <- FALSE
    cycle <- cycle + 1
    log <- c(log, paste("cycle:", cycle))
    log <- c(log, paste("ungrouped.hapix:", paste(ungrouped.hapix, collapse=" ")))
    groups.hapix <- list() # same as result$groups, but with index numbers to
    #                        haplonr instead of haplonr ID's)
    again.hapix <- integer(0) #haplotypes temporarily unplaced, to try again later
    for (hi in order(hap.nacount)) {
      if (!(hi %in% c(nodata.hapix,ungrouped.hapix))) {
        #note that since haplonr, haplofreq and haploseq are already sorted by
        #descending haplofreq, we here take the haplotypes with the same number of
        #missing marker scores in order of descending frequency
        hap <- haplonr[hi]
        #compare hi with all existing non-singlet groups and add to existing or new group:
        foundgroup <- assign.hapix(hi)
        if (is.na(foundgroup)) {
          again.hapix <- c(again.hapix, hi)
          log <- c(log, paste("placed in again.hapix:", hi))
        } else {
          if (length(groups.hapix)<foundgroup) {
            groups.hapix[[foundgroup]] <- integer(0)
          }
          groups.hapix[[foundgroup]] <- c(groups.hapix[[foundgroup]],hi)
          log <- c(log, paste("group", foundgroup, ": added", hi))
        }
      }
    } #for h
    # next we try to add the haplotypes set aside because they had no match within
    # their selected group:
    changed <- TRUE
    while (changed && length(again.hapix)>0) {
      log <- c(log, paste("add from again.hapix, contents are",
                          paste(again.hapix, collapse=" ")))
      for (hi in again.hapix) {
        changed <- FALSE
        foundgroup <- assign.hapix(hi)
        if (is.na(foundgroup)) {
          #no change
          log <- c(log, paste("remains in again.hapix:", hi))
        } else {
          if (length(groups.hapix)<foundgroup) {
            groups.hapix[[foundgroup]] <- integer(0)
          }
          groups.hapix[[foundgroup]] <- c(groups.hapix[[foundgroup]],hi)
          again.hapix <- setdiff(again.hapix, hi)
          log <- c(log, paste("group", foundgroup,
                              ": added from again.hapix ", hi))
          changed <- TRUE
        }
      } #for hi
    } #while changed
    #now we know that within a group all haplotypes match, but we must remove the
    #ones that match between groups
    # we do that one by one, in order of decreasing number of missing marker scores
    consensus <- matrix(integer(length(groups.hapix)*ncol(haploseq)),
                        ncol=ncol(haploseq))
    for (gi in seq_along(groups.hapix)) {
      consensus[gi,] <- combineHaplotypes(haploseq[groups.hapix[[gi]],,drop=FALSE])
    }
    remove.order <- order(-hap.nacount, haplofreq) #decreasing nacount, then increasing freq
    log <- c(log, paste("remove.order:", paste(remove.order, collapse=" ")))
    removed <- FALSE
    i <- 1
    while (!removed && i<=length(remove.order)) {
      hi <- remove.order[i] # hi is index to haplonr etc
      logline <- paste("i=",i," hi=",hi," haplonr=",haplonr[hi],sep="")
      if (hi %in% ungrouped.hapix) {
        #ungrouped.hapix contains one haplotype more in each cycle
        log <- c(log,paste(logline," already in ungrouped",sep=""))
      } else if (hi %in% again.hapix) {
        #again.hapix may be different in each cycle
        log <- c(log,paste(logline," in again.hapix",sep=""))
      } else if (hi %in% nodata.hapix) {
        #nodata.hapix contains one or no haplotypes: the one with all markerdata missing
        log <- c(log,paste(logline," in nodata.hapix",sep=""))
      } else {
        #where is this haplotype?
        gi <- length(groups.hapix)
        while (gi>0 && !(hi %in% groups.hapix[[gi]])) gi <- gi-1
        if (gi<=0) stop ("GroupHaplotypes: gi<=0")
        #hi in group gi
        #with which groups does hi not conflict?
        matchgrp <- logical(length(groups.hapix))
        for (gj in seq_along(groups.hapix)) {
          matchgrp[gj] <- match.haplo(haploseq[hi,], consensus[gj,],
                                      match.NA=FALSE, no.info=TRUE)
        }
        if (!(gi %in% which(matchgrp))) stop ("gi not in matchgrp")
        if (sum(matchgrp) > 1) {
          #hi has no conflict with multiple groups and must be removed
          ungrouped.hapix <- c(ungrouped.hapix, hi)
          groups.hapix[[gi]] <- setdiff(groups.hapix[[gi]],hi)
          removed <- TRUE
          log <- c(log,paste(logline," from group ",gi, " also matches group(s)",
                             paste(setdiff(which(matchgrp),gi)),"; removed", sep=""))
        } else {
          log <- c(log,paste(logline," has only matches in its group ",gi, sep=""))
        }

      }
      i <- i+1
    } #while i
  } #while removed

  #There should now not be any empty groups: if in the last cycle a group with one
  #haplotype was created, that same haplotype does not match any other group
  #and therefore will not be removed;
  #but just in case we check for and remove all empty groups:
  emptygrp <- sapply(groups.hapix,"length")==0
  groups.hapix <- groups.hapix[!emptygrp]
  if (sum(emptygrp)>0) log <- c(log,"emptygrp not empty, this should not happen")

  #Finally, add back all from ungrouped.hapix that now match only one remaining group:
  #If we add a new haplotype to an existing group we do not increase its number
  #of NA markers, so it is not possible that haplotypes in other groups will now
  #match this group (although they may match this one new haplotype)
  #If we would create a new group because a haplotype does not match any existing
  #group, it is possible that haplotypes within existing groups do match the new
  #group (with only one haplotype).
  #However it seems impossible that haplotypes not matching any existing group
  #have ended up in ungrouped.hapix.
  #In any case, we now add back haplotypes from ungrouped.hapix ONLY if they
  #do not conflict with exactly one group and have real matches with that group,
  #and we will not create any new groups.
  ungrouped.hapix <- c(ungrouped.hapix,again.hapix)
  rm(again.hapix)
  ungrouped.hapix <- ungrouped.hapix[order(haplofreq[ungrouped.hapix], decreasing=TRUE)]
  changed <- TRUE
  while (changed && length(ungrouped.hapix)>0) {
    log <- c(log, paste("add from ungrouped.hapix, contents are",
                        paste(ungrouped.hapix, collapse=" ")))
    for (hi in ungrouped.hapix) {
      changed <- FALSE
      foundgroup <- assign.hapix(hi)
      if (is.na(foundgroup) || foundgroup>length(groups.hapix)) {
        #no change
        log <- c(log, paste("remains in ungrouped.hapix:", hi))
      } else {
        groups.hapix[[foundgroup]] <- c(groups.hapix[[foundgroup]],hi)
        ungrouped.hapix <- setdiff(again.hapix, hi)
        log <- c(log, paste("group", foundgroup, ": added from ungrouped.hapix ", hi))
        changed <- TRUE
      }
    } #for hi
  } #while changed

  consensus <- matrix(integer(length(groups.hapix)*ncol(haploseq)),
                      ncol=ncol(haploseq))
  for (gi in seq_along(groups.hapix)) {
    consensus[gi,] <- combineHaplotypes(haploseq[groups.hapix[[gi]],,drop=FALSE])
  }

  # calculate the maximum possible freq of each group:
  # we could extend this to check if the different haplotypes from ungrouped hapix
  # that match a given group conflict with each other (in that case they could not
  # all be added to the group). But for now that seems a waste of time.
  try.hapix <- c(ungrouped.hapix, nodata.hapix)
  maxgrpfrq <- integer(length(groups.hapix))
  for (gi in seq_along(groups.hapix)) {
    maxgrp.hapix <- groups.hapix[[gi]]
    for (tryhi in try.hapix) {
      if (sum(!result$matchmatrix[tryhi,maxgrp.hapix],na.rm=TRUE)==0) {
        maxgrp.hapix <- c(maxgrp.hapix,tryhi)
      }
    }
    maxgrpfrq[gi] <- sum(haplofreq[maxgrp.hapix])
  }


  #next we sort the groups in order of decreasing frequency:
  #(note that there are no empty groups, see above)
  grpfrq <- integer(length(groups.hapix))
  for (gi in seq_along(groups.hapix)) {
    grpfrq[gi] <-sum(haplofreq[groups.hapix[[gi]]])
  }
  grporder <- order(grpfrq, decreasing=TRUE)
  result$groups <- list()
  g <- 1
  while (g<=length(grporder)) {
    #now we must convert the haplotype indices to haplotype ID's:
    grp <- grporder[g]
    nwg <- length(result$groups) + 1
    result$groups[[nwg]] <- haplonr[groups.hapix[[grp]]]
    g <- g+1
  }
  result$consensus <- consensus[grporder,,drop=FALSE]
  #ungrouped is already ordered by decreasing frequency
  result$ungrouped <- haplonr[ungrouped.hapix]
  result$nodata <- haplonr[nodata.hapix]
  result$groupfreq <- (grpfrq[grporder])[seq_along(result$groups)] #some groups at end may be empty
  result$maxgroupfreq <- maxgrpfrq[grporder][seq_along(result$groups)]
  result$ungroupedfreq <- sum(haplofreq[ungrouped.hapix])
  result$log <- log
  result
} #groupHaplotypes

#write the components of the return value of groupHaplotypes (for one fp)
write.grouped.haplotypes <- function(con, fp, grouping, fp_haplotypes,
                                     fphaplofreq, haplolen) {
  #grouping: a list as returned by groupHaplotypes
  #fphaplofreq: integer vector with frequencies of all haplotypes for this fp
  startline <- paste(fp,
                     names(fp_haplotypes)[fp],
                     length(fp_haplotypes[[fp]]$markers),
                     sep="\t")
  for (g in seq_along(grouping$groups)) {
    for (hap in grouping$groups[[g]]) {
      midline <- paste(hap,
                       fp_haplotypes[[fp]]$na.count[hap],
                       fphaplofreq[hap], g,
                       #fp_haplotypes[[fp]]$freq[hap], g,
                       grouping$groupfreq[g], grouping$maxgroupfreq[g],
                       sep="\t")
      endline <- paste(fp_haplotypes[[fp]]$hapmat[hap,],
                       sep="\t", collapse="\t")
      if (ncol(fp_haplotypes[[fp]]$hapmat)<haplolen) {
        for (i in 1:(haplolen-ncol(fp_haplotypes[[fp]]$hapmat)))
          endline <- paste(endline,"\t",sep="")
      }
      endline2 <- paste(grouping$consensus[g,],
                        sep="\t", collapse="\t")
      writeLines(paste(startline,midline,endline,"-",endline2,sep="\t"),con)
    }
  }
  for (hap in grouping$ungrouped) {
    midline <- paste(hap,
                     fp_haplotypes[[fp]]$na.count[hap],
                     fphaplofreq[hap], "ungrouped",
                     grouping$ungroupedfreq,"-",
                     sep="\t")
    endline <- paste(fp_haplotypes[[fp]]$hapmat[hap,],
                     sep="\t", collapse="\t")
    writeLines(paste(startline,midline,endline,sep="\t"),con)
  }
  for (hap in grouping$nodata) {
    midline <- paste(hap,
                     fp_haplotypes[[fp]]$na.count[hap],
                     fphaplofreq[hap], "nodata",
                     "-","-",
                     sep="\t")
    endline <- paste(fp_haplotypes[[fp]]$hapmat[hap,],
                     sep="\t", collapse="\t")
    writeLines(paste(startline,midline,endline,sep="\t"),con)
  }
} #write.grouped.haplotypes

groupNwrite.all.haploblockalleles <- function(fp_haplotypes, hapnrarr, filename) {
  #group all haplotypes in data per focalpoint and write to file
  #fp_haplotypes: the fp_haplotypes object returned by get_initial_haplotypes
  #hapnrarr: 3-dim array of haplotype numbers (alleles) with dimensions
  #          focalpoint, individual, hap1/hap2, as returned by
  #          get_initial_haplotypes
  #return value: a list with for each focalpoint a grouping by groupHaplotypes
  #side effect: file filename is written: a tab-separated file with the grouped
  #haplotypes for all focalpoints.
  #Note that this is the overall grouping over all pedigree individuals;
  #the grouping per family will often be different: absence of some haplotypes
  #will affect the grouping of the remaining haplotypes
  res <- list()
  for (fp in seq_along(fp_haplotypes)) {
    haplonr <- 1:length(fp_haplotypes[[fp]]$na.count)
    haplofreq <- table(hapnrarr[fp,,])
    haplofreq <- haplofreq[match(seq_along(fp_haplotypes[[fp]]$na.count),
                                 names(haplofreq))]
    haplofreq[is.na(haplofreq)] <- 0
    haploseq <- fp_haplotypes[[fp]]$hapmat
    or <- order(haplofreq, decreasing=T)
    haplonr <- haplonr[or]
    haplofreq <- haplofreq[or]
    haploseq <- haploseq[or,, drop=FALSE]
    res[[fp]] <- groupHaplotypes(haplonr, haploseq,haplofreq)
    #write(res[[fp]]$log, paste("fp",fp,"_logfile",ver,".txt", sep=""))
  }
  #write all the haplotypes to file:
  con <- file(filename, "w")
  haplolen <- 0
  for (fp in seq_along(fp_haplotypes)) {
    if (haplolen < ncol(fp_haplotypes[[fp]]$hapmat))
      haplolen <- ncol(fp_haplotypes[[fp]]$hapmat)
  }
  captions <- "fpnr\tfocalpoint\tmarkercount\thaplonr\tnacount\tfreq\tgroup\tgrpfrq\tmaxgrpfrq\tmarker_alleles"
  for (i in 1:haplolen) captions <- paste(captions,"\t",sep="")
  captions <- paste(captions,"-\tgrpconsensus",sep="")
  writeLines(captions,con)
  for (fp in 1:length(fp_haplotypes)) {
    #get freq of all haplotypes:
    fphaplofreq <- table(hapnrarr[fp,,])
    fphaplofreq <- fphaplofreq[match(seq_along(fp_haplotypes[[fp]]$na.count),
                                     names(fphaplofreq))]
    fphaplofreq[is.na(fphaplofreq)] <- 0
    write.grouped.haplotypes(con, fp, res[[fp]], fp_haplotypes,
                             fphaplofreq, haplolen)
  }
  close(con)
  invisible(res) #list of the groupings of all focalpoints,
  #               won't print if not assigned
} #groupNwrite.all.haploblockalleles


#   ___________________________________________________________________________
#   Approach of processHSfam
#   VERSION 28-12-2013
#   DONE: omit uninformative / interfering marker scores:
#         where 2 parentals (should) match same group but do not match each other:
#         join haplotypes if possible by setting one marker score to NA in both
#         where 2 parentals identical and match 2 groups, with one group small:
#         join groups if possible by setting one marker score to missing
#   DONE: within each cycle: set clearly incorrect haplotypes or marker scores
#         to NA, and set NA with clear solution to consensus. We dont set a wrong
#         haplotype directly to consensus, to allow filling the NA in the next
#         round by other evidence
#   TODO: idea: if a hi-freq ungrouped matches 2 conflicting groups of which
#         at least one is lo-freq, see if omitting one marker solves problem
#         (esp if this improves fit with parental). We then need a
#         function to find haplotype that fits all haplotypes in vector by
#         setting one marker to NA in all (and still has each haplotype having at
#         least one marker score overlapping with consensus)
#
#   1 group :
#   HSsize>=15, groupfreq incl matching unmatched (so: excl NAs) >
#     HSsize-qbinom(0.02,HSsize,1/3):
#   really one haplotype, no segregation, fill in all HS with consensus
#   (matching parentals also used for consensus and replaced, non-matching
#    parentals replaced by NA; OR: if both parentals match group but not each other,
#    see if setting one marker score to NA in parentals and group solves problem)
#
#   1 group:
#   HSsize>=15, number of NA >= qbinom(0.02,HSsize,1/3), or HSsize<15:
#   there is room for a second haplotype;
#   if both parentals match group and each other: we really have one group, fill
#   in all HS and parentals with consensus incl matching unmatched
#   if both parentals match group but not each other:
#     if HSsize>=15 and group freq>=HSsize/2: non-matching parentals replaced
#       by NA; OR: see if setting one marker score to NA in parentals and group
#       solves problem
#     else: set all HS to consensus and parentals to missing
#   if one parental matches and the other not or NA: second haplotype possible,
#   fill in group 1 (excl matching unmatched) and matching parental with consensus
#   leave remaining HS to NA, leave remaining parental as is (OR: fill in HS NAs
#   with other parental? only one sample supports!)
#   if both parentals not matching: if group large enough at least one parental
#   must be wrong; else group may be wrong
#   ->  group freq excl matching unmatched >2: keep group, leave other HS NA,
#   set parentals to NA; group freq ==2: set parentals and all HS to NA;
#   group freq ==1: set all HS but not parentals to NA
#
#   >=2 groups:
#   HSsize>=15, two largest groups larger than qbinom(0.02,HSsize,1/3): these are
#   the two real haplotypes.
#   if both parentals each match one group: make 2 consensus groups. From unmatched
#     add the haplotypes that match only one of the two groups. Create the overall
#     consensus and set parentals and groups; set all other HS to NA
#   if both parentals missing or one missing and one matching: similar (set missing
#     parental to 2nd group consensus)
#   if one or both parentals not matching: similar, set non-matching parentals
#     to NA
#   if both parentals match the same group: at least one is wrong, set both to NA,
#     and set the two groups incl matching unmatched to consensus per group
#   if parentals identical and match both groups: both groups are real, parentals
#     have at least one missing marker score compared to both groups or one parental
#     is wrong: as previous
#
#   >=2 groups:
#   HSsize>=15, only one group larger than qbinom(0.02,HSsize,1/3): there may be one
#   or two real haplotypes.
#   one parental matches largest group, other parental another group:
#     make 2 consensus groups. From unmatched haplotypes that match only one of
#     the two groups. Create the overall
#     consensus and set parentals and groups; set all other HS to NA
#   one parental matches largest group, other parental missing:
#     make 1 consensus group. add from unmatched all matching.
#     If this total group lets less than qbinom(0.02,HSsize,1/3) outside group: there
#     is only one haplotypes, set missing parental and all other HS also to
#     consensus.
#     Else there may still be two haplotypes: Set one parental plus group to
#     consensus, all other HS to NA
#   both parentals identical and match largest group and one of smaller groups:
#     there is only one haplotype but the 2 groups differ by at least one
#     unreliable marker and the parentals have a missing marker score. See if
#     setting one marker to missing in both groups solves problem. If yes: make
#     consensus score over both groups+parentals+matching unmatched, set all other
#     HS to NA. If no (too many discrepancies between groups): make 1 consensus
#     group, add from unmatched all matching; all other HS and both parentals to NA
#   both parentals missing:
#     make 1 consensus group. add from unmatched all matching. set all other HS
#     to NA.
#   one parental matches a smaller group, other parental missing:
#     assume missing parental matches the largest group, and then as above:
#     make 2 consensus groups. From unmatched haplotypes that match only one of
#     the two groups. Create the overall consensus and set both parentals and
#     groups; set all other HS to NA
#   one parental matches a smaller group, other parental doesnt match a group
#     assume non-matching parental is wrong and actually matches largest group,
#     and then as above:
#     make 2 consensus groups. From unmatched haplotypes that match only one of
#     the two groups. Create the overall consensus, set one parental to group
#     consensus and the other to NA; set all other HS to NA
#   both parentals match a different smaller group or no group:
#     make consensus of largest group. add from unmatched all matching. Set both
#     parentals to NA and also all HS outside largest group
#
#   >= 2 groups
#   HSsize>=15, no groups larger than qbinom(0.02,HSsize,1/3):
#   both parentals identical and match 2 groups:
#     there is only one haplotype but the 2 groups differ by at least one
#     unreliable marker and the parentals have a missing marker score. See if
#     setting one marker to missing in both groups solves problem. If yes: make
#     consensus score over both groups+parentals+matching unmatched, set all other
#     HS to NA. If no (too many discrepancies between groups): set all HS to NA,
#     leave parentals
#   both parentals match only the same group:
#     if parentals match each other: make consensus of parentals and group; add all
#     matching unmatched to group, set all other HS to NA
#     else (parentals dont match each other): make consensus of group, add all
#     matching unmatched: set all other HS and both parentals to NA
#   both parentals each match only one, different group:
#     make two consensus groups, add all unmatched that match only one of these
#     groups, set all other HS to NA
#   one parental matches only one group, other parental missing:
#     make one consensus group, add no unmatched, make all other HS NA
#   two parentals missing, or at least one parental not matching any group:
#     set all HS and both parentals to NA
#
#   >=2 groups
#   HSsize<15
#   both parentals identical and match 2 groups:
#     there is only one haplotype but the 2 groups differ by at least one
#     unreliable marker and the parentals have a missing marker score. See if
#     setting one marker to missing in both groups solves problem. If yes: make
#     consensus score over both groups+parentals+matching unmatched, set all other
#     HS to NA. If no (too many discrepancies between groups): set all HS to NA,
#     leave parentals
#   both parentals match only the same group:
#     if parentals match each other: make consensus of parentals and group; add all
#     matching unmatched to group, set all other HS to NA
#     else (parentals dont match each other): make consensus of group, add all
#     matching unmatched: set all other HS and both parentals to NA
#   both parentals each match only one, different group:
#     make two consensus groups, add all unmatched that match only one of these
#     groups, set all other HS to NA
#   one parental matches only one group, other parental missing:
#     if there are exactly 2 groups and the unmatched group freq >2:
#     make two consensus groups, add no unmatched, make all other HS NA, set one
#     parental to its group consensus and leave other parental NA;
#     else (>2 groups): make one consensus group, add no unmatched, make all other
#     HS NA
#   two parentals missing: if exactly two groups: make two consensus group, add no
#     unmatched, make all other HS NA;
#     else (>2 groups): make all HS missing
#   one parental matches a group, the other not:
#     if there are exactly 2 groups and the unmatched group has freq > 2:
#     make 2 consensus groups, add no unmatched, set non-matching parent to NA
#     else (more than 2 groups or second group has freq <=2): set non-matching
#     parental to NA, make only one consensus, add no unmatched, all other HS to NA
#   both parentals dont match groups: if exactly two groups both with freq>2:
#     make two consensus groups, add no unmatched, make all other HS and both
#     parentals NA;
#     else (>2 groups or freqs<=2): make all HS and parentals missing

#some helper functions for processHSfam:

calcMinHaplofreq <- function(famsize) {
  #famsize: number of individuals in the family
  #return value:
  #the minimum expected frequency of a real haplotype
  #(if this is one of the two haplotypes of parent par and its expected
  # fraction in the progeny is at least max.skewed, then we should observe
  # at least the returned number of progeny with this haplotype at the given
  # significance. The lower the significance (more significant), the lower
  # the minimum expected number of progeny)
  max.skewed = 1/3 # the most extreme skewedness, where unskewed is 0.5
  signif <- 0.02
  qbinom(signif, famsize, max.skewed)
} #calcMinHaplofreq

is.nodata <- function(haploseq) {
  #haploseq: a vector with marker alleles
  #return value: TRUE if all marker alleles are NA, else FALSE
  sum(!is.na(haploseq)) == 0
} #is.nodata

gethaplotypenr <- function(haploseq, fphaplo) {
  #haploseq: a vector with marker alleles
  #fphaplo: fp_haplotypes[[fp]] where fp is the index to the current focalpoint
  #return value: list with
  # $haplotypenr: integer
  # $fphaplo: NULL (meaning not present) if haploseq already in fphaplo,
  #           else fphaplo with haploseq added
  h <- 1
  while (h <= nrow(fphaplo$hapmat) &&
         !match.haplo(haploseq, fphaplo$hapmat[h,],
                     match.NA=TRUE, no.info=FALSE)) {
    h <- h + 1
  }
  result <- list()
  result$haplotypenr <- h
  result$fphaplo <- NULL
  if (h>nrow(fphaplo$hapmat)) {
    # add consensus to fphaplo:
    fphaplo$hapmat <- rbind(fphaplo$hapmat, haploseq)
    rownames(fphaplo$hapmat)[h] <- getnewhapnr(fphaplo$hapmat, h)
    fphaplo$na.count[h] <- sum(is.na(haploseq))
    #fphaplo$freq[h] <- NA
    result$fphaplo <- fphaplo
  }
  result
} #gethaplotypenr

addNArow <- function(m) {
  #m is vector or matrix
  #returns a matrix with a row of NA's rbind-ed to m
  if (is.vector(m)) {
    nas <- rep(NA,length(m))
  } else {
    nas <- rep(NA,ncol(m))
  }
  rbind(m,nas)
} #addNArow

checkNadd1 <- function(consensus, groups, ungrouped, fphaplo, check=TRUE) {
  #checks all haplotypes in ungrouped if they match with the (one) group.
  #matching ungrouped haplotypes are moved from ungrouped to group
  #and the consensus is updated
  #consensus: vector, consensus seq of group before calling
  #groups: integer vector, all haplotype numbers already in group
  #ungrouped: integer vector, all ungrouped haplotypes
  #fphaplo: fp_haplotypes[[fp]] where fp is the index to the current focalpoint
  #check: if TRUE (default) the contents of ungrouped are checked (again)
  #       against groups, else ungrouped and groups remain unchanged and this
  #       function just generates a result list of the unchanged input data
  #return value: list with items as if there were 2 groups,
  #compatible with checkNadd2:
  # $consensus: the updated consensus seq rbinded with a row of NA's
  # $groups: the updated group as groups[[1]], integer(0) as groups[[2]]
  # $ungrouped: the updated ungrouped
  # $consensushaplonr: c(NA, NA) (to be filled in by caller)
  # $addnodata: FALSE, to be filled in by caller (should the nodata
  #             progeny be added to groups 1 or not?)
  #note that fphaplo is not updated here, the initial and final consensus
  #may both not occur in fphaplo but all group and ungrouped haplotypes
  #do occur in fphaplo
  groups <- list(groups, integer(0)) #second group empty, for compatibility
  #                                   with checkNadd2
  if (check) for (h in ungrouped) {
    h.seq <- fphaplo$hapmat[h,]
    if (match.haplo(consensus, h.seq, match.NA=FALSE, no.info=FALSE)) {
      groups[[1]] <- c(groups[[1]], h)
      ungrouped <- setdiff(ungrouped, h)
      consensus <- combineHaplotypes(rbind(consensus, h.seq))
    }
  }
  list(consensus=addNArow(consensus), groups=groups, ungrouped=ungrouped,
       consensushaplonr = c(NA, NA), addnodata=FALSE)
} #checkNadd1

checkNadd2 <- function(consensus, groups, ungrouped, fphaplo, check=TRUE) {
  #checks all haplotypes in ungrouped if they match with only one of the two groups.
  #matching ungrouped haplotypes are moved from ungrouped to their group
  #and the group consensus is updated
  #consensus: 2-row matrix, consensus seq of groups before calling
  #group: list of 2 integer vectors, all haplotype numbers already in groups
  #ungrouped: integer vector, all ungrouped haplotypes
  #fphaplo: fp_haplotypes[[fp]] where fp is the index to the current focalpoint
  #check: if TRUE (default) the contents of ungrouped are checked (again)
  #       against groups, else ungrouped and groups remain unchanged and this
  #       function just generates a result list of the unchanged input data
  #return value: list with items:
  # $consensus: matrix with the two updated group consensus
  # $group: list: the two updated group
  # $ungrouped: the updated ungrouped
  # $consensushaplonr: c(NA, NA) (to be filled in by caller)
  # $addnodata: FALSE, to be filled in by caller (should the nodata
  #             progeny be added to groups 1 or not?)
  #note that fphaplo is not updated here, the initial and final consensus
  #may both not occur in fphaplo but all group and ungrouped haplotypes
  #do occur in fphaplo
  if (check) for (h in ungrouped) {
    h.seq <- fphaplo$hapmat[h,]
    match1 <- match.haplo(consensus[1,], h.seq, match.NA=FALSE, no.info=FALSE)
    match2 <- match.haplo(consensus[2,], h.seq, match.NA=FALSE, no.info=FALSE)
    if (xor(match1, match2)) {
      if (match1) gr <- 1 else gr <- 2
      groups[[gr]] <- c(groups[[gr]], h)
      ungrouped <- setdiff(ungrouped, h)
      consensus[gr,] <- combineHaplotypes(rbind(consensus[gr,], h.seq))
    }
  }
  if (!is.matrix(consensus)) stop("checkNadd2: not matrix")
  list(consensus=consensus, groups=groups, ungrouped=ungrouped,
       consensushaplonr = c(NA, NA), addnodata=FALSE)
} #checkNadd2

setHShaplonrs <- function (dat, HSnr, HSfam, fphapnrarr, fphaplo) {
  #dat: a list as produced by checkNadd1 or 2, specifying which haplonrs are
  #     in groups 1 and 2 and ungrouped, what their consensus seq are,
  #     which haplonrs correspond to the consensus seq (or NA) and whether
  #     nodata should be added to group 1
  #fphapnrarr: a 2-dim array of individuals * parental (1:2)
  #            with the haplotype numbers (but no NA: nodata haplotype is 1)
  #HSnr: number of HS family (index to HSfam)
  #par: 1 or 2, the parental haplotype
  #fphaplo: a list with info on the haplotypes of the current focalpoint
  #return value: a list with
  # $fphaplo = updated version or NULL
  # $fphapnrarr = updated version
  indf1 <- list() #2 elements: individuals in HSf1 belonging to group 1 and 2
  indf2 <- list() #..same for HSf2
  fphaplo.changed <- FALSE
  for (gr in 1:2) {
    if (is.na(dat$consensushaplonr[gr])) {
      #is the group consensus an existing haplotype?
      haplonr <- gethaplotypenr(dat$consensus[gr,], fphaplo)
      h <- haplonr$haplotypenr
      if (!is.null(haplonr$fphaplo)) {
        fphaplo <- haplonr$fphaplo
        fphaplo.changed <- TRUE
      }
      dat$consensushaplonr[gr] <- h
    }
    indf1[[gr]] <- which(fphapnrarr[, 1] %in% dat$groups[[gr]])
    if (dat$addnodata && gr==1) {
      #add the HS indiv with nodata haplotype 1 to group 1
      indnodata <- which(fphapnrarr[, 1] == 1) # 1 is nodata haplotype
      indf1[[gr]] <- union(indf1[[gr]], indnodata)
    }
    indf1[[gr]] <- intersect(indf1[[gr]], HSfam[[HSnr]]$HSf1)
    fphapnrarr[indf1[[gr]], 1] <- rep(dat$consensushaplonr[gr],
                                      length(indf1[[gr]]))
    indf2[[gr]] <- which(fphapnrarr[, 2] %in% dat$groups[[gr]])
    if (dat$addnodata && gr==1) {
      #add the HS indiv with nodata haplotype 1 to group 1
      indnodata <- which(fphapnrarr[, 2] == 1) # 1 is nodata haplotype
      indf2[[gr]] <- union(indf2[[gr]], indnodata)
    }
    indf2[[gr]] <- intersect(indf2[[gr]], HSfam[[HSnr]]$HSf2)
    fphapnrarr[indf2[[gr]], 2] <- rep(dat$consensushaplonr[gr],
                                      length(indf2[[gr]]))
  }
  # ... and set the HS individuals outside the groups to nodata = 1:
  ind <- setdiff(HSfam[[HSnr]]$HSf1, c(indf1[[1]], indf1[[2]]))
  fphapnrarr[ind, 1] <- 1
  ind <- setdiff(HSfam[[HSnr]]$HSf2, c(indf2[[1]], indf2[[2]]))
  fphapnrarr[ind, 2] <- 1
  if (!fphaplo.changed) fphaplo <- NULL
  list(fphaplo=fphaplo, fphapnrarr=fphapnrarr)
} #setHShaplonrs

#Here we implement the above strategy
#processHSfam looks at one HS family HSnr (with parent), for one focalpoint fp,
#and at only the haplotypes inherited from the common parent.
#It tries to determine if one or two haplotypes are inherited from the parent,
#and to resolve conflicts such as multiple haplotypes in the HS progeny,
#discrepancies due to missing marker alleles, conflicts between parent and
#offspring. It replaces conflicting haplotypes by nodata (= haplotype 1),
#and if possible imputes new haplotypes for nodata individuals.
#This function should be called multiple times until no further changes occur,
#or until an earlier configuration is reached again.
processHSfam <- function(fp, HSnr,
                         fphaplo, fphapnrarr, HSfam) {
  # fp: focalpoint (index number, not ID)
  # HSnr: number of HSfamily (indexes HSfam, the result of get.all.HSfamilies)
  # fphaplo: fp_haplotypes[[fp]], part of the object returned by
  #          get_initial_haplotypes
  # fphapnrarr: hapnrarr[fp,,], a 2-dim slice from the array of haplotype
  #             numbers (alleles) returned by get_initial_haplotypes
  # HSfam: a list as produced by get.all.HSfamilies
  #return value: a list with items
  # $fphaplo: possibly changed version of fp_haplotypes[[fp]] (new haplotypes added)
  # $fphaplo.changed: TRUE or FALSE; need to replace fp_haplotypes[[fp]]
  # $fphapnrarr: probably changed version of hapnrarr[fp,,]
  # TODO possible future extension: take into account the number of nodata
  #      haplotypes among HS individuals
  fphaplo.changed <- FALSE
  hapinfo <- get.HSfam.haplotypes(HSnr, fphapnrarr, HSfam)
  HSsize <- sum(hapinfo$HShaplofrq) #nr of HS individuals in HS family
  #        (with parent either as mother or father)
  parhaplo <- hapinfo$parhaplo # the two haplotypes of the parent
  HShaplofreq <- sort(hapinfo$HShaplofrq, decreasing=T) #frequencies of all
  #                haplotypes from parent HSfamily
  HShaplonr <- as.integer(names(HShaplofreq)) #the IDnumbers of these HSfam
  #                                          haplotypes in same order
  parhaploseq <- fphaplo$hapmat[parhaplo,, drop=FALSE]
  parnodata <- parhaplo == 1 # 1 is nodata haplotype
  if (sum(is.na(parnodata)) > 0) stop("processHSfam: parnodata")
  HShaploseq <- fphaplo$hapmat[HShaplonr,,drop=FALSE]
  minfrq <- calcMinHaplofreq(HSsize)
  matchmatrix <- calc.matchmatrix(HShaploseq, names=HShaplonr)
  grouping <- groupHaplotypes(HShaplonr, HShaploseq, HShaplofreq, matchmatrix)
  #####################
  ## from 20140117 - fq_haplotyping_exper9.r we consider a non-overlap between
  ## a parental and a group as a conflict,
  ## else there are so many situations that the approach gets unmanageable
  #####################
  # 20140218: check for non-informative haplotypes if 1 group:
  # If only one group is found, this group may include some
  # haplotypes that are hardly informative (many missing
  # marker scores) but that are grouped only because there
  # are no other haplotypes with which they might also match.
  # If there is at least one parental haplotype that does not
  # match with the group consensus but with which some of the
  # haplotypes in the groups do match, then those haplotypes
  # are not informative.
  # We replace them with nodata=1 in fphapnrarr and adjust the
  # grouping$groups[[1]] and grouping$groupfreq accordingly,
  # before proceeding with the original analysis
  if (length(grouping$groups) == 1) {
    matchpar <- logical(2)
    for (p in 1:2) {
      matchpar[p] <- match.haplo(grouping$consensus[1,], parhaploseq[p,],
                                 match.NA=FALSE, no.info=TRUE)
    }
    for (h in grouping$groups[[1]]) {
      matchfound <- FALSE
      for (p in 1:2) {
        if ((!matchfound) && (!matchpar[p]) &&
              match.haplo(fphaplo$hapmat[h,], parhaploseq[p,],
                          match.NA=FALSE, no.info=TRUE)) {
          #haplotype h from group 1 matches with a parental haplotype
          #that did not match the group consensus; set all HS individuals
          #with this haplotype to nodata=1:
          indgrp <- which(fphapnrarr[, 1] == h)
          ind <- intersect(HSfam[[HSnr]]$HSf1, indgrp)
          fphapnrarr[ind, 1] <- 1
          grouping$groupfreq[1] <- grouping$groupfreq[1] - length(ind)
          indgrp <- which(fphapnrarr[, 2] == h)
          ind <- intersect(HSfam[[HSnr]]$HSf2, indgrp)
          fphapnrarr[ind, 2] <- 1
          grouping$groupfreq[1] <- grouping$groupfreq[1] - length(ind)
          grouping$groups[[1]] <- setdiff(grouping$groups[[1]], h)
          grouping$nodata <- union(grouping$nodata, 1)
          matchfound <- TRUE
        }
      }
    } #for h
    if (xor(length(grouping$groups[[1]]) == 0, grouping$groupfreq[1] == 0) ||
        grouping$groupfreq[1] < 0) {
      stop("processHSfam: error checking one-group situation")
    }
    if (grouping$groupfreq[1] == 0) {
      #the one group is now empty; delete it:
      grouping$groups <- list()
      grouping$groupfreq <- integer(0)
      grouping$consensus <- grouping$consensus[-1,]
    }
  }
  #now we continue with the original analysis:
  if (length(grouping$groups) == 0) {
    #only nodata HS progeny
    #if both parentals exactly the same (match with match.NA TRUE) then all
    #progeny set to this haplotype,
    #else no change
    if (parhaplo[1] == parhaplo[2]) {
      fphapnrarr[HSfam[[HSnr]]$HSf1, 1] <- rep(parhaplo[1],
                                                   length(HSfam[[HSnr]]$HSf1))
      fphapnrarr[HSfam[[HSnr]]$HSf2, 2] <- rep(parhaplo[1],
                                                   length(HSfam[[HSnr]]$HSf2))
    }
  } else if (length(grouping$groups) == 1) {
    #1 group
    if (length(grouping$ungrouped)>0) {
      #in the following we assume that with one group, ungrouped must be empty;
      #but we print a message if this is not true:
      #WRONG: eg if there are two haplotypes that both have marker data but
      #       have no overlapping data, the most frequent will be grouped and
      #       the other left ungrouped
      #So we do not print the warning
      #cat(paste("Warning in processHSfam: 1 group, but ungrouped is not empty",
      #          "\n", sep=""))
    }
    if (HSsize>=15 && HSsize-grouping$groupfreq[1] < minfrq) {
      # there is really only this group
      consensus <- grouping$consensus[1,]
      if (match.haplo(consensus, parhaploseq[1,],
                      match.NA=FALSE, no.info=FALSE) &&
            match.haplo(consensus, parhaploseq[2,],
                        match.NA=FALSE, no.info=FALSE)) {
        #both parentals match group.
        if (match.haplo(parhaploseq[1,], parhaploseq[2,],
                        match.NA=FALSE, no.info=FALSE)) {
          #both parentals also match each other; get consensus with parentals
          consensus <- combineHaplotypes(rbind(consensus, parhaploseq))
        } else {
          #parentals both match group consensus but not each other.
          #due to one or more marker scores that are NA in group consensus
          #but different in parentals.
          #use group consensus also for both parentals
        }
      } else {
        #not both parentals match this group, but if one matches we make a
        #consensus with this one parental:
        for (p in 1:2) {
          if (match.haplo(consensus, parhaploseq[p,],
                          match.NA=FALSE, no.info=FALSE)) {
            consensus <- combineHaplotypes(rbind(consensus, parhaploseq[p,]))
          }
        }
      } #not both parentals match group
      #is the consensus an existing haplotype?
      haplonr <- gethaplotypenr(consensus, fphaplo)
      h <- haplonr$haplotypenr
      if (!is.null(haplonr$fphaplo)) {
        fphaplo <- haplonr$fphaplo
        fphaplo.changed <- TRUE
      }
      #replace the haplotype of all HS individuals with h:
      #(all HS are already in group[[1]] or in nodata)
      fphapnrarr[HSfam[[HSnr]]$HSf1, 1] <- rep(h, length(HSfam[[HSnr]]$HSf1))
      fphapnrarr[HSfam[[HSnr]]$HSf2, 2] <- rep(h, length(HSfam[[HSnr]]$HSf2))
      #replace the parental(s) that do not conflict with consensus with h,
      #and the parental(s) that conflict with nodata = 1:
      for (p in 1:2) {
        if (match.haplo(consensus, parhaploseq[p,],
                        match.NA=FALSE, no.info=TRUE)) {
          fphapnrarr[HSfam[[HSnr]]$parentnr, p] <- h
        } else {
          fphapnrarr[HSfam[[HSnr]]$parentnr, p] <- 1
        }
      }
    } else {
      #1 group, HSsize<15 or nodata>=minfrq, there might be a second group
      #if both parentals match group there is again one group:
      consensus <- grouping$consensus[1,]
      if (match.haplo(consensus, parhaploseq[1,],
                      match.NA=FALSE, no.info=FALSE) &&
            match.haplo(consensus, parhaploseq[2,],
                        match.NA=FALSE, no.info=FALSE)) {
        #both parentals match group.
        consensus <- grouping$consensus[1,]
        if (match.haplo(parhaploseq[1,], parhaploseq[2,],
                        match.NA=FALSE, no.info=FALSE)) {
          #both parentals also match each other; get consensus with parentals
          for (p in 1:2) {
            consensus <- combineHaplotypes(rbind(consensus, parhaploseq[p,]))
          }
        } else {
          #parentals both match group consensus but not each other.
          #due to one or more marker scores that are NA in group consensus
          #but different in parentals.
          #use group consensus also for both parentals
        }
        #is the consensus an existing haplotype?
        haplonr <- gethaplotypenr(consensus, fphaplo)
        h <- haplonr$haplotypenr
        if (!is.null(haplonr$fphaplo)) {
          fphaplo <- haplonr$fphaplo
          fphaplo.changed <- TRUE
        }
        #replace the haplotype of all HS individuals with h:
        #(all HS are already in group[[1]] or in nodata)
        fphapnrarr[HSfam[[HSnr]]$HSf1, 1] <- rep(h, length(HSfam[[HSnr]]$HSf1))
        fphapnrarr[HSfam[[HSnr]]$HSf2, 2] <- rep(h, length(HSfam[[HSnr]]$HSf2))
        #replace the two parentals (both match with group) also with consensus:
        for (p in 1:2) {
          fphapnrarr[HSfam[[HSnr]]$parentnr, p] <- h
        }
      } else {
        #not both parentals match group. If one parental matches and the other
        #is  nodata or not matching we set the group + one parent to overall
        #consensus and leave the other parental and other HS as they were:
        if (match.haplo(consensus, parhaploseq[1,],
                        match.NA=FALSE, no.info=FALSE) ||
            match.haplo(consensus, parhaploseq[2,],
                        match.NA=FALSE, no.info=FALSE)) {
          #one parent matches and the other nodata or not matching
          if (match.haplo(consensus, parhaploseq[1,],
                          match.NA=FALSE, no.info=FALSE)) {
            p <- 1
          } else {
            p <- 2
          }
          consensus <- combineHaplotypes(rbind(grouping$consensus[1,],
                                               parhaploseq[p,]))
          #is the consensus an existing haplotype?
          haplonr <- gethaplotypenr(consensus, fphaplo)
          h <- haplonr$haplotypenr
          if (!is.null(haplonr$fphaplo)) {
            fphaplo <- haplonr$fphaplo
            fphaplo.changed <- TRUE
          }
          #replace the haplotype of the group individuals with h:
          indgrp <- which(fphapnrarr[, 1] %in% grouping$groups[[1]])
          ind <- intersect(HSfam[[HSnr]]$HSf1, indgrp)
          fphapnrarr[ind, 1] <- h
          indgrp <- which(fphapnrarr[, 2] %in% grouping$groups[[1]])
          ind <- intersect(HSfam[[HSnr]]$HSf2, indgrp)
          fphapnrarr[ind, 2] <- h
          #replace the parental that matches with consensus with h,
          #and leave the other parental unchanged:
          fphapnrarr[HSfam[[HSnr]]$parentnr, p] <- h
        } else {
          #no parentals match group: each is nodata or non-matching
          #first get the group consensus haplotype:
          haplonr <- gethaplotypenr(grouping$consensus[1,], fphaplo)
          h <- haplonr$haplotypenr
          if (!is.null(haplonr$fphaplo)) {
            fphaplo <- haplonr$fphaplo
            fphaplo.changed <- TRUE
          }
          #if group large enough at least one parental should match it;
          #else group may be wrong
          if (grouping$groupfreq[1] > 2) {
            #keep group, leave other HS nodata = 1;
            #if one parental is nodata, set this to consensus and leave the other;
            #else set both parentals to nodata = 1
            #is the consensus an existing haplotype?
            haplonr <- gethaplotypenr(grouping$consensus[1,], fphaplo)
            h <- haplonr$haplotypenr
            if (!is.null(haplonr$fphaplo)) {
              fphaplo <- haplonr$fphaplo
              fphaplo.changed <- TRUE
            }
            #replace the haplotype of the group individuals with h:
            indgrp <- which(fphapnrarr[, 1] %in% grouping$groups[[1]])
            ind <- intersect(HSfam[[HSnr]]$HSf1, indgrp)
            fphapnrarr[ind, 1] <- h
            indgrp <- which(fphapnrarr[, 2] %in% grouping$groups[[1]])
            ind <- intersect(HSfam[[HSnr]]$HSf2, indgrp)
            fphapnrarr[ind, 2] <- h
            #set parentals:
            if (xor(parnodata[1], parnodata[2])) {
              #one parental nodata and the other is conflicting
              if (parnodata[1]) {
                fphapnrarr[HSfam[[HSnr]]$parentnr, 1] <- h
              } else {
                fphapnrarr[HSfam[[HSnr]]$parentnr, 2] <- h
              }
            } else {
              #both parentals nodata or both conflicting; set both to nodata = 1:
              for (p in 1:2) {
                fphapnrarr[HSfam[[HSnr]]$parentnr, p] <- 1
              }
            }
          } else {
            #one group with only <=2 individuals;
            #both parentals don't match: each is nodata or non-matching
            #if at least one parent nodata we set the group to consensus,
            #leave rest progeny nodata, and leave parentals as they are,
            #else we set all progeny and both parentals to nodata = 1
            if (parnodata[1] || parnodata[2]) {
              #set the group to consensus, leave rest progeny nodata
              #is the consensus an existing haplotype?
              haplonr <- gethaplotypenr(grouping$consensus[1,], fphaplo)
              h <- haplonr$haplotypenr
              if (!is.null(haplonr$fphaplo)) {
                fphaplo <- haplonr$fphaplo
                fphaplo.changed <- TRUE
              }
              indgrp <- which(fphapnrarr[, 1] %in% grouping$groups[[1]])
              ind <- intersect(HSfam[[HSnr]]$HSf1, indgrp)
              fphapnrarr[ind, 1] <- h
              indgrp <- which(fphapnrarr[, 2] %in% grouping$groups[[1]])
              ind <- intersect(HSfam[[HSnr]]$HSf2, indgrp)
              fphapnrarr[ind, 2] <- h
            } else {
              #set parentals and all HS to nodata = 1:
              fphapnrarr[HSfam[[HSnr]]$HSf1, 1] <- 1
              fphapnrarr[HSfam[[HSnr]]$HSf2, 2] <- 1
              fphapnrarr[HSfam[[HSnr]]$parentnr, 1] <- 1
              fphapnrarr[HSfam[[HSnr]]$parentnr, 2] <- 1
            }
          } #one group with only <=2 individuals
        } #both parentals nodata, not matching or conflicting
      } #not both parentals match group
    } #1 group, HSsize<15 or nodata>=minfrq
  } else {
    # >=2 groups
    if (HSsize >= 15) {
      # >= 2 groups, HSsize >= 15
      dat <- list() #will contain updated info on groups 1 and 2 and ungrouped
      if (grouping$groupfreq[2] >= minfrq) {
        # HSsize >= 15, >= 2 groups, largest 2 >= minfrq
        # probably 2 real haplotypes segregating
        # first add all ungrouped that match only one of the two groups to that group:
        dat <- checkNadd2(consensus=grouping$consensus[1:2,, drop=FALSE],
                          groups=grouping$groups[1:2],
                          ungrouped=grouping$ungrouped,
                          fphaplo=fphaplo)
        #20140314:
        #aim: place all haplotypes that are not in dat$groups into ungrouped
        #but: dat$groups are now changed from grouping$groups (some former
        #ungrouped haplotypes may now be in groups)
        #So: take all haplotypes in grouping$groups (do.call) and move all
        #of those that are not in dat$groups (groupedhap) to dat$ungrouped:
        groupedhap <- union(dat$groups[[1]], dat$groups[[2]])
        dat$ungrouped <- c(dat$ungrouped,
                           setdiff(do.call(c, grouping$groups), groupedhap))
        #next check the parentals:
        matchgrp <- matrix(logical(4), ncol=2)
        for (gr in 1:2) for (p in 1:2) {
          matchgrp[gr,p] <- match.haplo(dat$consensus[gr,], parhaploseq[p,],
                                        match.NA=FALSE, no.info=FALSE)
        }
        #parentals that match more than one group are set to nodata = 1
        #(they have missing data for some important markers)
        #and matchgrp is updated accordingly:
        morematches <- which(colSums(matchgrp) > 1)
        for (p in morematches) {
          fphapnrarr[HSfam[[HSnr]]$parentnr, p] <- 1
          parnodata[p] <- TRUE
          matchgrp[,p] <- rep(FALSE, nrow(matchgrp))
        }
        #in the following situations the parentals are set or kept to nodata:
        if (sum(matchgrp) == 0 || #no matches between any parent and any group
              sum(rowSums(matchgrp)==2) != 0) { #one or both groups match with both parentals
          # no reliable data on parentals, set to nodata = 1:
          fphapnrarr[HSfam[[HSnr]]$parentnr, 1] <- 1
          fphapnrarr[HSfam[[HSnr]]$parentnr, 2] <- 1
        } else if (xor(parnodata[1], parnodata[2])) {
          # one parental known, and it must match exactly one group
          # (if matches 0 groups it would be captured in the first if,
          # and if it matched more groups it is already set to nodata)
          if (parnodata[1]) p <- 2 else p <- 1
          g <- which(matchgrp[,p])
          #calculate overall consensus of parental p and its group g:
          dat$consensus[g,] <- combineHaplotypes(rbind(dat$consensus[g,],
                                                       parhaploseq[p,]))
          #try again to add ungrouped haplotypes:
          dat <- checkNadd2(consensus=dat$consensus,
                            groups=dat$groups,
                            ungrouped=dat$ungrouped,
                            fphaplo=fphaplo)
          #are the group consensus existing haplotypes?
          for (gr in 1:2) {
            haplonr <- gethaplotypenr(dat$consensus[gr,], fphaplo)
            if (!is.null(haplonr$fphaplo)) {
              fphaplo <- haplonr$fphaplo
              fphaplo.changed <- TRUE
            }
            dat$consensushaplonr[gr] <- haplonr$haplotypenr
          }
          #set parental p to group consensus:
          fphapnrarr[HSfam[[HSnr]]$parentnr, p] <- dat$consensushaplonr[g]
          #set unknown parental to other group consensus:
          fphapnrarr[HSfam[[HSnr]]$parentnr, 3-p] <- dat$consensushaplonr[3-g]
        } else {
          #both parentals known
          #cases:
          #- each parental matches a different group: ok
          #- both parentals match the same group: already caught
          #- one parental matches a group and the other not: set the non-
          #  matching to nodata=1 and calculate consensus of matching parental
          #  + group
          #- both parentals don't match group: already caught
          #- one or both parentals match more than 1 group: these were already
          #  set to nodata
          if (sum(colSums(matchgrp)==1) == 2) {
            #both parentals match a group, and these groups must be different
            if (matchgrp[1,1]) p <- 1 else p <- 2
            #now parental 1 matches group g and parental 2 matches group 3-g
            dat$consensus[1,] <- combineHaplotypes(rbind(dat$consensus[1,],
                                                         parhaploseq[p,]))
            dat$consensus[2,] <- combineHaplotypes(rbind(dat$consensus[2,],
                                                         parhaploseq[3-p,]))
            #try again to add ungrouped haplotypes:
            dat <- checkNadd2(consensus=dat$consensus,
                              groups=dat$groups,
                              ungrouped=dat$ungrouped,
                              fphaplo=fphaplo)
            #are the group consensus existing haplotypes?
            for (gr in 1:2) {
              haplonr <- gethaplotypenr(dat$consensus[gr,], fphaplo)
              if (!is.null(haplonr$fphaplo)) {
                fphaplo <- haplonr$fphaplo
                fphaplo.changed <- TRUE
              }
              dat$consensushaplonr[gr] <- haplonr$haplotypenr
            }
            #set parentals to group consensus:
            fphapnrarr[HSfam[[HSnr]]$parentnr, 1] <- dat$consensushaplonr[p]
            fphapnrarr[HSfam[[HSnr]]$parentnr, 2] <- dat$consensushaplonr[3-p]
          } else {
            #only one parent matches a group (and not two groups)
            #this differs from the situation where only one parental known:
            #here we set the known but non-matching parental to nodata = 1,
            #there we set the unknown parental to the other group consensus
            p <- which(colSums(matchgrp) == 1)
            g <- which(matchgrp[,p])
            #calculate overall consensus of this parental p and its group g:
            dat$consensus[g,] <- combineHaplotypes(rbind(dat$consensus[g,],
                                                         parhaploseq[p,]))
            #try again to add ungrouped haplotypes:
            dat <- checkNadd2(consensus=dat$consensus,
                              groups=dat$groups,
                              ungrouped=dat$ungrouped,
                              fphaplo=fphaplo)
            #are the group consensus existing haplotypes?
            for (gr in 1:2) {
              haplonr <- gethaplotypenr(dat$consensus[gr,], fphaplo)
              if (!is.null(haplonr$fphaplo)) {
                fphaplo <- haplonr$fphaplo
                fphaplo.changed <- TRUE
              }
              dat$consensushaplonr[gr] <- haplonr$haplotypenr
            }
            #set parental p to group consensus:
            fphapnrarr[HSfam[[HSnr]]$parentnr, p] <- dat$consensushaplonr[g]
            #set other, conflicting parental to nodata = 1:
            fphapnrarr[HSfam[[HSnr]]$parentnr, 3-p] <- 1
          } #both parentals known, only one matches a group
        } #both parentals known
        #end of different parental situations for HSsize>=15, >= 2 groups, 2 groups >= minfreq
        #finally set the haplonrs of the HS individuals:
        update <- setHShaplonrs(dat, HSnr, HSfam, fphapnrarr, fphaplo)
        fphapnrarr <- update$fphapnrarr
        if (!is.null(update$fphaplo)) {
          fphaplo.changed <- TRUE
          fphaplo <- update$fphaplo
        }
        #end of HSsize >= 15, >= 2 groups, largest 2 >= minfrq
      } else if (grouping$groupfreq[1] >= minfrq) {
        # HSsize >= 15, >= 2 groups, only one >= minfrq
        matchgrp <- matrix(logical(2*length(grouping$groups)), ncol=2)
        for (gr in seq_along(grouping$groups)) for (p in 1:2) {
          matchgrp[gr,p] <- match.haplo(grouping$consensus[gr,], parhaploseq[p,],
                                        match.NA=FALSE, no.info=FALSE)
        }
        #parentals that match more than one group are set to nodata = 1
        #(they have missing data for some important markers)
        #and matchgrp is updated accordingly:
        morematches <- which(colSums(matchgrp) > 1)
        for (p in morematches) {
          fphapnrarr[HSfam[[HSnr]]$parentnr, p] <- 1
          parnodata[p] <- TRUE
          matchgrp[,p] <- rep(FALSE, nrow(matchgrp))
        }
        if (sum(matchgrp) == 0) { #no matches between any parental and any group
          # no reliable data on parentals, set to nodata=1
          # and keep only the large group
          fphapnrarr[HSfam[[HSnr]]$parentnr, 1] <- 1
          fphapnrarr[HSfam[[HSnr]]$parentnr, 2] <- 1
          dat <- checkNadd1(consensus=grouping$consensus[1,],
                            groups=grouping$groups[[1]],
                            ungrouped=grouping$ungrouped,
                            fphaplo,
                            check=FALSE)
          #for (gr in 2:length(grouping$groups)) {
          #  dat$ungrouped <- c(dat$ungrouped, grouping$group[[gr]])
          #}
          #20140314:
          #aim: place all haplotypes that are not in dat$groups into ungrouped
          #but: dat$groups are now changed from grouping$groups (some former
          #ungrouped haplotypes may now be in groups)
          #So: take all haplotypes in grouping$groups (do.call) and move all
          #of those that are not in dat$groups (groupedhap) to dat$ungrouped:
          groupedhap <- union(dat$groups[[1]], dat$groups[[2]])
          dat$ungrouped <- c(dat$ungrouped,
                             setdiff(do.call(c, grouping$groups), groupedhap))

        } else if (xor(parnodata[1], parnodata[2])) {
          #one known parental and this matches with exactly one group.
          #if it matches the largest group, get the consensus, re-check the
          #ungrouped and set the other parental and all F1 indivs outside the
          #largest group to nodata.
          #if it matches a minor group, get the consensus with this minor group,
          #re-check the consensus
          #and assign the consensus of the large group to the other parental
          p <- 3 - which(parnodata)
          g <- which(matchgrp[,p])
          if (g==1) {
            #parental matches large group
            dat <- checkNadd1(consensus=combineHaplotypes(rbind(dat$consensus[1,],
                                                                parhaploseq[p,])),
                              groups=grouping$groups[[1]],
                              ungrouped=grouping$ungrouped,
                              fphaplo)
            #for (gr in 2:length(grouping$groups)) {
            #  dat$ungrouped <- c(dat$ungrouped, grouping$group[[gr]])
            #}
            #20140314:
            #aim: place all haplotypes that are not in dat$groups into ungrouped
            #but: dat$groups are now changed from grouping$groups (some former
            #ungrouped haplotypes may now be in groups)
            #So: take all haplotypes in grouping$groups (do.call) and move all
            #of those that are not in dat$groups (groupedhap) to dat$ungrouped:
            groupedhap <- union(dat$groups[[1]], dat$groups[[2]])
            dat$ungrouped <- c(dat$ungrouped,
                               setdiff(do.call(c, grouping$groups), groupedhap))

            #is new consensus a new haplotype?
            haplonr <- gethaplotypenr(dat$consensus[1,], fphaplo)
            if (!is.null(haplonr$fphaplo)) {
              fphaplo <- haplonr$fphaplo
              fphaplo.changed <- TRUE
            }
            dat$consensushaplonr[1] <- haplonr$haplotypenr
            #set parental p to group consensus:
            fphapnrarr[HSfam[[HSnr]]$parentnr, p] <- dat$consensushaplonr[1]
          } else {
            #parental matches a smaller group
            consensus2 <- combineHaplotypes(rbind(grouping$consensus[g,],
                                                  parhaploseq[p,]))
            dat <- checkNadd2(consensus=rbind(grouping$consensus[1,],
                                              consensus2),
                              groups=grouping$groups[c(1, g)],
                              ungrouped=grouping$ungrouped,
                              fphaplo=fphaplo)
            #for (gr in 2:length(grouping$groups)) {
            #  if (gr != g) dat$ungrouped <- c(dat$ungrouped, grouping$group[[gr]])
            #}
            #20140314:
            #aim: place all haplotypes that are not in dat$groups into ungrouped
            #but: dat$groups are now changed from grouping$groups (some former
            #ungrouped haplotypes may now be in groups)
            #So: take all haplotypes in grouping$groups (do.call) and move all
            #of those that are not in dat$groups (groupedhap) to dat$ungrouped:
            groupedhap <- union(dat$groups[[1]], dat$groups[[2]])
            dat$ungrouped <- c(dat$ungrouped,
                               setdiff(do.call(c, grouping$groups), groupedhap))
            #is new consensus a new haplotype?
            for (gr in 1:2) {
              haplonr <- gethaplotypenr(dat$consensus[gr,], fphaplo)
              if (!is.null(haplonr$fphaplo)) {
                fphaplo <- haplonr$fphaplo
                fphaplo.changed <- TRUE
              }
              dat$consensushaplonr[gr] <- haplonr$haplotypenr
            }
            #set parental p to group consensus:
            fphapnrarr[HSfam[[HSnr]]$parentnr, p] <- dat$consensushaplonr[2]
            #set unknown parental to large group consensus:
            fphapnrarr[HSfam[[HSnr]]$parentnr, 3-p] <- dat$consensushaplonr[1]
          }
        } else {
          #2 parentals known, each matches one or no groups,
          #we have the following possibilities:
          # - both parentals match the large group: only one haplotype segregating,
          #   if both parentals match each other we take the overall consensus,
          #   else we use the group consensus. We set both parentals, HS in group
          #   and all nodata HS to consensus, all other HS to nodata
          # - both parentals match the same smaller group: we keep the smaller and the
          #   large groups (consensus with both parentals if these match each
          #   other) but set both parentals to nodata
          # - one parental matches the large group, one a smaller group:
          #   we keep both groups, consensus with parentals, re-check ungrouped
          # - both parentals match a different smaller group: keep only the
          #   largest group and set both parentals and all other F1 indivs to nodata
          # - one parental matches the large group, the other conflicts with all
          #   groups: keep consensus of parental and large group, re-check
          #   ungrouped, set other parental and all other HS indiv to nodata
          # - one parental matches small group, other conflicts with all groups:
          #   get consensus of small group and parental, re-check ungrouped,
          #   set other parental and all HS outside large and small group
          #   to nodata = 1
          # - both parentals conflict with all groups: already caught
          # - one or both parentals match multiple groups: already caught
          if (sum(colSums(matchgrp)) == 1) {
            #one parental matches a group, the other does not;
            #similar to if one parental known, but now the non-matching parental
            #is set to nodata=1 instead to set to the large group consensus:
            p <- which(colSums(matchgrp) > 0)
            g <- which(matchgrp[,p])
            #set non-matching parental to nodata = 1:
            fphapnrarr[HSfam[[HSnr]]$parentnr, 3-p] <- 1
            if (g==1) {
              #parental matches large group
              dat <- checkNadd1(consensus=combineHaplotypes(rbind(dat$consensus[1,],
                                                                  parhaploseq[p,])),
                                groups=grouping$groups[[1]],
                                ungrouped=grouping$ungrouped,
                                fphaplo)
              dat$addnodata <- TRUE
              #for (gr in 2:length(grouping$groups)) {
              #  dat$ungrouped <- c(dat$ungrouped,grouping$group[[gr]])
              #}
              #20140314:
              #aim: place all haplotypes that are not in dat$groups into ungrouped
              #but: dat$groups are now changed from grouping$groups (some former
              #ungrouped haplotypes may now be in groups)
              #So: take all haplotypes in grouping$groups (do.call) and move all
              #of those that are not in dat$groups (groupedhap) to dat$ungrouped:
              groupedhap <- union(dat$groups[[1]], dat$groups[[2]])
              dat$ungrouped <- c(dat$ungrouped,
                                 setdiff(do.call(c, grouping$groups), groupedhap))
              #is new consensus a new haplotype?
              haplonr <- gethaplotypenr(dat$consensus[1,], fphaplo)
              if (!is.null(haplonr$fphaplo)) {
                fphaplo <- haplonr$fphaplo
                fphaplo.changed <- TRUE
              }
              dat$consensushaplonr[1] <- haplonr$haplotypenr
              #set parental p to group consensus:
              fphapnrarr[HSfam[[HSnr]]$parentnr, p] <- dat$consensushaplonr[1]
            } else {
              #parental matches a smaller group g, the other matches no group
              consensus2 <- combineHaplotypes(rbind(grouping$consensus[g,],
                                                    parhaploseq[p,]))
              dat <- checkNadd2(consensus=rbind(grouping$consensus[1,],
                                                consensus2),
                                groups=grouping$groups[c(1, g)],
                                ungrouped=grouping$ungrouped,
                                fphaplo=fphaplo)
              #for (gr in 2:length(grouping$groups)) {
              #  if (gr != g) dat$ungrouped <- c(dat$ungrouped,
              #                                  grouping$group[[gr]])
              #}
              #20140314:
              #aim: place all haplotypes that are not in dat$groups into ungrouped
              #but: dat$groups are now changed from grouping$groups (some former
              #ungrouped haplotypes may now be in groups)
              #So: take all haplotypes in grouping$groups (do.call) and move all
              #of those that are not in dat$groups (groupedhap) to dat$ungrouped:
              groupedhap <- union(dat$groups[[1]], dat$groups[[2]])
              dat$ungrouped <- c(dat$ungrouped,
                                 setdiff(do.call(c, grouping$groups), groupedhap))
              #is new consensus a new haplotype?
              for (gr in 1:2) {
                haplonr <- gethaplotypenr(dat$consensus[gr,], fphaplo)
                if (!is.null(haplonr$fphaplo)) {
                  fphaplo <- haplonr$fphaplo
                  fphaplo.changed <- TRUE
                }
                dat$consensushaplonr[gr] <- haplonr$haplotypenr
              }
              #set parental p to group consensus:
              fphapnrarr[HSfam[[HSnr]]$parentnr, p] <- dat$consensushaplonr[2]
            }
          } else {
            #both parentals match a group
            g <- integer(2) #group matched by parentals 1 and 2
            g[1] <- which(matchgrp[,1]); g[2] <- which(matchgrp[,2])
            if (g[1] == g[2]) {
              #both parentals match the same group
              #if both parentals match each other we use the consensus of all three,
              #else only the group consensus
              if (match.haplo(parhaploseq[1,], parhaploseq[2,],
                              match.NA=FALSE, no.info=FALSE)) {
                consensus <- combineHaplotypes(rbind(grouping$consensus[g[1],],
                                                     parhaploseq[1,],
                                                     parhaploseq[2,]))
              } else {
                consensus <- grouping$consensus[g[1],]
              }
              if (g[1] == 1) {
                #both parentals match the large group; there is only one haplotype
                dat <- checkNadd1(consensus=consensus,
                                  groups=grouping$groups[[1]],
                                  ungrouped=grouping$ungrouped,
                                  fphaplo)
                #for (gr in 2:length(grouping$groups)) {
                #  dat$ungrouped <- c(dat$ungrouped,grouping$group[[gr]])
                #}
                #20140314:
                #aim: place all haplotypes that are not in dat$groups into ungrouped
                #but: dat$groups are now changed from grouping$groups (some former
                #ungrouped haplotypes may now be in groups)
                #So: take all haplotypes in grouping$groups (do.call) and move all
                #of those that are not in dat$groups (groupedhap) to dat$ungrouped:
                groupedhap <- union(dat$groups[[1]], dat$groups[[2]])
                dat$ungrouped <- c(dat$ungrouped,
                                   setdiff(do.call(c, grouping$groups), groupedhap))
                #add nodata HS to group 1:
                dat$groups[[1]] <- c(dat$groups[[1]], grouping$nodata)
                #is new consensus a new haplotype?
                haplonr <- gethaplotypenr(dat$consensus[1,], fphaplo)
                if (!is.null(haplonr$fphaplo)) {
                  fphaplo <- haplonr$fphaplo
                  fphaplo.changed <- TRUE
                }
                dat$consensushaplonr[1] <- haplonr$haplotypenr
                #set parentals to group consensus:
                fphapnrarr[HSfam[[HSnr]]$parentnr, 1] <- dat$consensushaplonr[1]
                fphapnrarr[HSfam[[HSnr]]$parentnr, 2] <- dat$consensushaplonr[1]
              } else {
                #both parentals match a smaller group; we keep that group and the
                #large group but set both parentals to nodata = 1
                dat <- checkNadd2(consensus=rbind(grouping$consensus[1,],
                                                  consensus),
                                  groups=grouping$groups[c(1, g[1])],
                                  ungrouped=grouping$ungrouped,
                                  fphaplo=fphaplo)
                #for (gr in 2:length(grouping$groups)) {
                #  if (gr != g[1]) dat$ungrouped <- c(dat$ungrouped,
                #                                     grouping$group[[gr]])
                #}
                #20140314:
                #aim: place all haplotypes that are not in dat$groups into ungrouped
                #but: dat$groups are now changed from grouping$groups (some former
                #ungrouped haplotypes may now be in groups)
                #So: take all haplotypes in grouping$groups (do.call) and move all
                #of those that are not in dat$groups (groupedhap) to dat$ungrouped:
                groupedhap <- union(dat$groups[[1]], dat$groups[[2]])
                dat$ungrouped <- c(dat$ungrouped,
                                   setdiff(do.call(c, grouping$groups), groupedhap))
                #set parentals:
                fphapnrarr[HSfam[[HSnr]]$parentnr, 1] <- 1
                fphapnrarr[HSfam[[HSnr]]$parentnr, 2] <- 1
              }
            } else {
              #both parentals match a different group
              if (sum(g == 1) == 0) {
                #both parentals match a different smaller group:
                #keep only the largest group and set both parentals and
                #all other F1 indivs to nodata = 1
                dat <- checkNadd1(consensus=grouping$consensus[1,],
                                  groups=grouping$groups[[1]],
                                  ungrouped=grouping$ungrouped,
                                  fphaplo=fphaplo,
                                  check=FALSE)
                #for (gr in 2:length(grouping$groups)) {
                #  dat$ungrouped <- c(dat$ungrouped,grouping$group[[gr]])
                #}
                #20140314:
                #aim: place all haplotypes that are not in dat$groups into ungrouped
                #but: dat$groups are now changed from grouping$groups (some former
                #ungrouped haplotypes may now be in groups)
                #So: take all haplotypes in grouping$groups (do.call) and move all
                #of those that are not in dat$groups (groupedhap) to dat$ungrouped:
                groupedhap <- union(dat$groups[[1]], dat$groups[[2]])
                dat$ungrouped <- c(dat$ungrouped,
                                   setdiff(do.call(c, grouping$groups), groupedhap))
                #set parental haplotypes:
                fphapnrarr[HSfam[[HSnr]]$parentnr, 1] <- 1
                fphapnrarr[HSfam[[HSnr]]$parentnr, 2] <- 1
              } else {
                #one parental matches the large group, the other parental a small group
                #set both groups with its parental to their consensus
                p <- which(matchgrp[1,]) #the parental matching the large group
                g <- which(matchgrp[,3-p]) #the group matching the other parental
                consensus <- rbind(
                  combineHaplotypes(rbind(grouping$consensus[1,],
                                          parhaploseq[p,])),
                  combineHaplotypes(rbind(grouping$consensus[g,],
                                          parhaploseq[3-p,])))
                dat <- checkNadd2(consensus=consensus,
                                  groups=grouping$groups[c(1, g)],
                                  ungrouped=grouping$ungrouped,
                                  fphaplo=fphaplo)
                #for (gr in 2:length(grouping$groups)) {
                #  if (gr != g) dat$ungrouped <- c(dat$ungrouped,
                #                                  grouping$group[[gr]])
                #}
                #20140314:
                #aim: place all haplotypes that are not in dat$groups into ungrouped
                #but: dat$groups are now changed from grouping$groups (some former
                #ungrouped haplotypes may now be in groups)
                #So: take all haplotypes in grouping$groups (do.call) and move all
                #of those that are not in dat$groups (groupedhap) to dat$ungrouped:
                groupedhap <- union(dat$groups[[1]], dat$groups[[2]])
                dat$ungrouped <- c(dat$ungrouped,
                                   setdiff(do.call(c, grouping$groups), groupedhap))
                #is new consensus a new haplotype?
                for (gr in 1:2) {
                  haplonr <- gethaplotypenr(dat$consensus[gr,], fphaplo)
                  if (!is.null(haplonr$fphaplo)) {
                    fphaplo <- haplonr$fphaplo
                    fphaplo.changed <- TRUE
                  }
                  dat$consensushaplonr[gr] <- haplonr$haplotypenr
                }
                #set parentals to group consensus:
                fphapnrarr[HSfam[[HSnr]]$parentnr, p] <-
                  dat$consensushaplonr[1]
                fphapnrarr[HSfam[[HSnr]]$parentnr, 3-p] <-
                  dat$consensushaplonr[2]
              } #one parental matches the large, the other a small group
            } #both parentals match a different group
          } #both parentals match a group
        } #2 parentals known, each matches one or no groups
        #end of different parental situations for HSsize>=15, >= 2 groups, 1 group >= minfreq
        #finally set the haplonrs of the HS individuals:
        update <- setHShaplonrs(dat, HSnr, HSfam, fphapnrarr, fphaplo)
        fphapnrarr <- update$fphapnrarr
        if (!is.null(update$fphaplo)) {
          fphaplo.changed <- TRUE
          fphaplo <- update$fphaplo
        }
        #end of HSsize >= 15, >= 2 groups, 1 group >= minfreq
      } else {
        # HSsize >= 15, >= 2 groups, no groups >= minfrq
        #TODO: where we have matching parentals and HS groups we now update
        #      both with the overall consensus; but does the evidence from the
        #      small group warrant the filling in of many markers in the parent?
        #      The same question applies in earlier cases where we update
        #      groups with their consensus with a matching parent: one parent
        #      supplies missing marker data for the group.
        #      For the moment we leave it as it is now.
        #      (for the one-group case we have already an alternative check)
        matchgrp <- matrix(logical(2*length(grouping$groups)), ncol=2)
        for (gr in seq_along(grouping$groups)) for (p in 1:2) {
          matchgrp[gr,p] <- match.haplo(grouping$consensus[gr,], parhaploseq[p,],
                                        match.NA=FALSE, no.info=FALSE)
        }
        #parentals that match more than one group are set to nodata = 1
        #(they have missing data for some important markers)
        #and matchgrp is updated accordingly:
        morematches <- which(colSums(matchgrp) > 1)
        for (p in morematches) {
          fphapnrarr[HSfam[[HSnr]]$parentnr, p] <- 1
          parnodata[p] <- TRUE
          matchgrp[,p] <- rep(FALSE,nrow(matchgrp))
        }
        if (sum(matchgrp) == 0) { #no matches between any parental and any group
          # no confirmed data on parentals nor on any group,
          # set both parentals and all F1 progeny to nodata = 1:
          fphapnrarr[HSfam[[HSnr]]$parentnr, 1] <- 1
          fphapnrarr[HSfam[[HSnr]]$parentnr, 2] <- 1
          dat <- list()
          dat$consensus <- rep(NA, ncol(parhaploseq))
          dat$consensus <- rbind(dat$consensus, dat$consensus)
          dat$groups <- list()
          dat$groups[[1]] <- integer(0); dat$groups[[2]] <- integer(0)
          dat$ungrouped <- grouping$ungrouped
          for (gr in seq_along(grouping$groups)) {
            dat$ungrouped <- c(dat$ungrouped, grouping$groups[[gr]])
          }
          dat$consensushaplonr = c(1, 1) # both nodata
          dat$addnodata <- FALSE
        } else if (sum(colSums(matchgrp) == 1) == 1) {
          #one parental matches one group, other parental missing or
          #matches no group:
          #make one consensus group, add no unmatched, make all other HS
          #and other parental nodata = 1
          p <- which(colSums(matchgrp) == 1)
          g <- which(matchgrp[,p])
          consensus <- combineHaplotypes(rbind(grouping$consensus[g,],
                                               parhaploseq[p,]))
          dat <- checkNadd1(consensus=rbind(consensus),
                            groups=grouping$groups[[g]], #g has 1 value, this is a vector
                            ungrouped=grouping$ungrouped,
                            fphaplo=fphaplo,
                            check=FALSE)
          #for (gr in seq_along(grouping$groups)) {
          #  if (gr != g) dat$ungrouped <- c(dat$ungrouped,grouping$group[[gr]])
          #}
          #20140314:
          #aim: place all haplotypes that are not in dat$groups into ungrouped
          #but: dat$groups are now changed from grouping$groups (some former
          #ungrouped haplotypes may now be in groups)
          #So: take all haplotypes in grouping$groups (do.call) and move all
          #of those that are not in dat$groups (groupedhap) to dat$ungrouped:
          groupedhap <- union(dat$groups[[1]], dat$groups[[2]])
          dat$ungrouped <- c(dat$ungrouped,
                             setdiff(do.call(c, grouping$groups), groupedhap))
          #is new consensus a new haplotype?
          haplonr <- gethaplotypenr(dat$consensus[1,], fphaplo)
          if (!is.null(haplonr$fphaplo)) {
            fphaplo <- haplonr$fphaplo
            fphaplo.changed <- TRUE
          }
          dat$consensushaplonr[1] <- haplonr$haplotypenr
          #set parental p to group consensus:
          fphapnrarr[HSfam[[HSnr]]$parentnr, p] <- dat$consensushaplonr[1]
          #set other parental to nodata = 1:
          fphapnrarr[HSfam[[HSnr]]$parentnr, 3-p] <- 1
        } else {
          #2 parentals known, each matches one group,
          #we have the following possibilities:
          # - both parentals match same group:
          #   if parentals match each other we use total consensus and recheck
          #   ungrouped, else use group consensus. Set parentals and group (but
          #   not nodata HS) to consensus, all other HS to nodata
          # - both parentals match a different group:
          #   use both groups + parental consensus, recheck ungrouped
          # - one parental matches a group, the other not: already caught
          # - both parentals match no groups: already caught
          g <- integer(2) #group matched by parentals 1 and 2
          g[1] <- which(matchgrp[,1]); g[2] <- which(matchgrp[,2])
          if (g[1] == g[2]) {
            #both parentals match the same group
            #if both parentals match each other we use the consensus of all three,
            #else only the group consensus
            if (match.haplo(parhaploseq[1,], parhaploseq[2,],
                            match.NA=FALSE, no.info=FALSE)) {
              consensus <- combineHaplotypes(rbind(grouping$consensus[g[1],],
                                                   parhaploseq[1,],
                                                   parhaploseq[2,]))
            } else {
              consensus <- grouping$consensus[g[1],]
            }
            dat <- checkNadd1(consensus=grouping$consensus[g[1],],
                              groups=grouping$groups[[g[1]]],
                              ungrouped=grouping$ungrouped,
                              fphaplo=fphaplo)
            #for (gr in seq_along(grouping$groups)) {
            #  if (gr != g[1]) dat$ungrouped <- c(dat$ungrouped,
            #                                     grouping$group[[gr]])
            #}
            #20140314:
            #aim: place all haplotypes that are not in dat$groups into ungrouped
            #but: dat$groups are now changed from grouping$groups (some former
            #ungrouped haplotypes may now be in groups)
            #So: take all haplotypes in grouping$groups (do.call) and move all
            #of those that are not in dat$groups (groupedhap) to dat$ungrouped:
            groupedhap <- union(dat$groups[[1]], dat$groups[[2]])
            dat$ungrouped <- c(dat$ungrouped,
                               setdiff(do.call(c, grouping$groups), groupedhap))
            #is new consensus a new haplotype?
            haplonr <- gethaplotypenr(dat$consensus[1,], fphaplo)
            if (!is.null(haplonr$fphaplo)) {
              fphaplo <- haplonr$fphaplo
              fphaplo.changed <- TRUE
            }
            dat$consensushaplonr[1] <- haplonr$haplotypenr
            fphapnrarr[HSfam[[HSnr]]$parentnr, 1] <- dat$consensushaplonr[1]
            fphapnrarr[HSfam[[HSnr]]$parentnr, 2] <- dat$consensushaplonr[1]
          } else {
            #both parentals match a different group
            #use the consensuses of both groups
            consensus <- rbind(combineHaplotypes(rbind(grouping$consensus[g[1],],
                                                       parhaploseq[1,])),
                               combineHaplotypes(rbind(grouping$consensus[g[2],],
                                                       parhaploseq[2,])))
            dat <- checkNadd2(consensus=consensus,
                              groups=grouping$groups[g], #g has 2 items, this is a list
                              ungrouped=grouping$ungrouped,
                              fphaplo=fphaplo)
            #old:
            #for (gr in 2:length(grouping$groups)) {
            #  if (!gr %in% g) dat$ungrouped <- c(dat$ungrouped,
            #                                     grouping$group[[gr]])
            #}
            #20140314:
            #aim: place all haplotypes that are not in dat$groups into ungrouped
            #but: dat$groups are now changed from grouping$groups (some former
            #ungrouped haplotypes may now be in groups)
            #So: take all haplotypes in grouping$groups (do.call) and move all
            #of those that are not in dat$groups (groupedhap) to dat$ungrouped:
            groupedhap <- union(dat$groups[[1]], dat$groups[[2]])
            dat$ungrouped <- c(dat$ungrouped,
                               setdiff(do.call(c, grouping$groups), groupedhap))
            #is new consensus a new haplotype?
            for (gr in 1:2) {
              haplonr <- gethaplotypenr(dat$consensus[gr,], fphaplo) #corrected, gr was 1
              if (!is.null(haplonr$fphaplo)) {
                fphaplo <- haplonr$fphaplo
                fphaplo.changed <- TRUE
              }
              dat$consensushaplonr[gr] <- haplonr$haplotypenr
            }
            #set parentals to group consensus:
            fphapnrarr[HSfam[[HSnr]]$parentnr, 1] <- dat$consensushaplonr[1] # corrected, was g[1]
            fphapnrarr[HSfam[[HSnr]]$parentnr, 2] <- dat$consensushaplonr[2] # corrected, was g[2]
          } #both parentals match a different group
        } #2 parentals known, each matches one or no groups
        #end of different parental situations for HSsize>=15, >= 2 groups, no group >= minfreq
        #finally set the haplonrs of the HS individuals:
        update <- setHShaplonrs(dat, HSnr, HSfam, fphapnrarr, fphaplo)
        fphapnrarr <- update$fphapnrarr
        if (!is.null(update$fphaplo)) {
          fphaplo.changed <- TRUE
          fphaplo <- update$fphaplo
        }
      } #end of HSsize >= 15, >= 2 groups, no groups >= minfrq
    } else {
      # HSsize < 15, >= 2 groups
      matchgrp <- matrix(logical(2*length(grouping$groups)), ncol=2)
      for (gr in seq_along(grouping$groups)) for (p in 1:2) {
        matchgrp[gr,p] <- match.haplo(grouping$consensus[gr,], parhaploseq[p,],
                                      match.NA=FALSE, no.info=FALSE)
      }
      #parentals that match more than one group are set to nodata = 1
      #(they have missing data for some important markers)
      #and matchgrp is updated accordingly:
      morematches <- which(colSums(matchgrp) > 1)
      for (p in morematches) {
        fphapnrarr[HSfam[[HSnr]]$parentnr, p] <- 1
        parnodata[p] <- TRUE
        matchgrp[,p] <- rep(FALSE, nrow(matchgrp))
      }
      if (sum(parnodata) == 2) {
        #no parental data; if exactly 2 groups keep
        # these, else set also all HS to nodata
        if (length(grouping$groups) == 2) {
          #exactly 2 groups, keep both
          dat <- checkNadd2(consensus=grouping$consensus,
                            groups=grouping$groups,
                            ungrouped=grouping$ungrouped,
                            fphaplo=fphaplo,
                            check=FALSE)
        } else {
          #more than 2 groups, set all progeny to nodata = 1
          dat <- checkNadd1(consensus=rbind(rep(NA, ncol(parhaploseq))),
                            groups=integer(0),
                            ungrouped=grouping$ungrouped,
                            fphaplo=fphaplo,
                            check=FALSE)
          #for (gr in seq_along(grouping$groups)) {
          #  dat$ungrouped <- c(dat$ungrouped, grouping$groups[[gr]])
          #}
          dat$consensushaplonr <- c(1, 1)
          #20140314:
          #aim: place all haplotypes that are not in dat$groups into ungrouped
          #but: dat$groups are now changed from grouping$groups (some former
          #ungrouped haplotypes may now be in groups)
          #So: take all haplotypes in grouping$groups (do.call) and move all
          #of those that are not in dat$groups (groupedhap) to dat$ungrouped:
          groupedhap <- union(dat$groups[[1]], dat$groups[[2]])
          dat$ungrouped <- c(dat$ungrouped,
                             setdiff(do.call(c, grouping$groups), groupedhap))
        }
      } else if (sum(parnodata) == 1) {
        #one parental known
        p <- 3 - which(parnodata)
        if (sum(matchgrp[, p]) == 1) {
          #the one known parental matches a group:
          #use consensus for group and parental, recheck ungrouped
          g <- which(matchgrp[, p])
          consensus <- combineHaplotypes(rbind(grouping$consensus[g,],
                                               parhaploseq[p,]))
          #if exactly 2 groups, keep also other group,
          #else set all other groups to nodata = 1
          if (length(grouping$groups) == 2) {
            #exactly 2 groups
            dat <- checkNadd2(consensus=rbind(consensus,
                                              grouping$consensus[3-g,]),
                              groups=grouping$groups[c(g, 3-g)],
                              ungrouped=grouping$ungrouped,
                              fphaplo=fphaplo)
            #is new consensus a new haplotype?
            haplonr <- gethaplotypenr(dat$consensus[1,], fphaplo)
            if (!is.null(haplonr$fphaplo)) {
              fphaplo <- haplonr$fphaplo
              fphaplo.changed <- TRUE
            }
            dat$consensushaplonr[1] <- haplonr$haplotypenr
            #set parental p to group consensus:
            fphapnrarr[HSfam[[HSnr]]$parentnr, p] <- dat$consensushaplonr[1]
          } else {
            #there are more than 2 groups, we set all except the one matching
            #the parent to nodata = 1
            dat <- checkNadd1(consensus=consensus,
                              groups=grouping$groups[[g]], #g has 1 value, this is a vector
                              ungrouped=grouping$ungrouped,
                              fphaplo=fphaplo)
            #is new consensus a new haplotype?
            haplonr <- gethaplotypenr(dat$consensus[1,], fphaplo)
            if (!is.null(haplonr$fphaplo)) {
              fphaplo <- haplonr$fphaplo
              fphaplo.changed <- TRUE
            }
            dat$consensushaplonr[1] <- haplonr$haplotypenr
            #set parental p to group consensus:
            fphapnrarr[HSfam[[HSnr]]$parentnr, p] <- dat$consensushaplonr[1]
            #set all other groups to nodata = 1:
            for (gr in seq_along(grouping$groups)) {
              if (gr != g) dat$ungrouped <- c(dat$ungrouped,
                                              grouping$groups[[gr]])
            }
          }
        } else {
          #one known parental doesn't match any group,
          #we set everything to nodata = 1
          dat <- list()
          dat$groups <- list(integer(0),integer(0))
          dat$consensushaplonr <- c(1, 1)
          dat$addnodata <- FALSE
          fphapnrarr[HSfam[[HSnr]]$parentnr, p] <- 1
        }
      } else {
        #both parentals known, each matches 0 or 1 group
        if (sum(matchgrp) == 0) {
          #both parentals match no group, we set everything to nodata = 1
          dat <- checkNadd1(consensus=rbind(rep(NA, ncol(parhaploseq))),
                            groups=integer(0),
                            ungrouped=grouping$ungrouped,
                            fphaplo=fphaplo,
                            check=FALSE)
          #for (gr in seq_along(grouping$groups)) {
          #  dat$ungrouped <- c(dat$ungrouped, grouping$groups[[gr]])
          #}
          #20140314:
          #aim: place all haplotypes that are not in dat$groups into ungrouped
          #but: dat$groups are now changed from grouping$groups (some former
          #ungrouped haplotypes may now be in groups)
          #So: take all haplotypes in grouping$groups (do.call) and move all
          #of those that are not in dat$groups (groupedhap) to dat$ungrouped:
          groupedhap <- union(dat$groups[[1]], dat$groups[[2]])
          dat$ungrouped <- c(dat$ungrouped,
                             setdiff(do.call(c, grouping$groups), groupedhap))
          dat$consensushaplonr <- c(1, 1)
          fphapnrarr[HSfam[[HSnr]]$parentnr, 1] <- 1
          fphapnrarr[HSfam[[HSnr]]$parentnr, 2] <- 1
        } else if (sum(matchgrp) == 1) {
          #one parental matches a group, the other doesn't;
          #we use group&parental consensus,
          #set other parental and groups to nodata = 1
          p <- which(colSums(matchgrp) == 1)
          g <- which(matchgrp[,p])
          dat <- checkNadd1(consensus=combineHaplotypes(rbind(grouping$consensus[g,],
                                                              parhaploseq[p,])),
                            groups=grouping$groups[[g]], #g has 1 value, this is a vector
                            ungrouped=grouping$ungrouped,
                            fphaplo=fphaplo)
          #for (gr in seq_along(grouping$groups)) {
          #  if (gr != g) dat$ungrouped <- c(dat$ungrouped,grouping$group[[gr]])
          #}
          #20140314:
          #aim: place all haplotypes that are not in dat$groups into ungrouped
          #but: dat$groups are now changed from grouping$groups (some former
          #ungrouped haplotypes may now be in groups)
          #So: take all haplotypes in grouping$groups (do.call) and move all
          #of those that are not in dat$groups (groupedhap) to dat$ungrouped:
          groupedhap <- union(dat$groups[[1]], dat$groups[[2]])
          dat$ungrouped <- c(dat$ungrouped,
                             setdiff(do.call(c, grouping$groups), groupedhap))
          #is new consensus a new haplotype?
          haplonr <- gethaplotypenr(dat$consensus[1,], fphaplo)
          if (!is.null(haplonr$fphaplo)) {
            fphaplo <- haplonr$fphaplo
            fphaplo.changed <- TRUE
          }
          dat$consensushaplonr[1] <- haplonr$haplotypenr
          #set parental p to group consensus:
          fphapnrarr[HSfam[[HSnr]]$parentnr, p] <- dat$consensushaplonr[1]
          #set other parental to nodata = 1:
          fphapnrarr[HSfam[[HSnr]]$parentnr, 3-p] <- 1
        } else {
          #2 parentals known, each matches one group,
          #we have the following possibilities:
          # - both parentals match same group:
          #   if parentals match each other we use total consensus and recheck
          #   ungrouped, else use group consensus. Set parentals and group (but
          #   not nodata HS) to consensus, all other HS to nodata = 1
          # - both parentals match a different group:
          #   use both groups + parental consensus, recheck ungrouped
          # - one parental matches a group, the other not: already caught
          # - both parentals match no groups: already caught
          g <- integer(2) #group matched by parentals 1 and 2
          g[1] <- which(matchgrp[,1]); g[2] <- which(matchgrp[,2])
          if (g[1] == g[2]) {
            #both parentals match the same group
            #if both parentals match each other we use the consensus of all three,
            #else only the group consensus
            if (match.haplo(parhaploseq[1,], parhaploseq[2,],
                            match.NA=FALSE, no.info=FALSE)) {
              consensus <- combineHaplotypes(rbind(grouping$consensus[g[1],],
                                                   parhaploseq[1,],
                                                   parhaploseq[2,]))
            } else {
              consensus <- grouping$consensus[g[1],]
            }
            dat <- checkNadd1(consensus=grouping$consensus[g[1],],
                              groups=grouping$groups[[g[1]]],
                              ungrouped=grouping$ungrouped,
                              fphaplo=fphaplo)
            #for (gr in seq_along(grouping$groups)) {
            #  if (gr != g[1]) dat$ungrouped <- c(dat$ungrouped,
            #                                     grouping$group[[gr]])
            #}
            #20140314:
            #aim: place all haplotypes that are not in dat$groups into ungrouped
            #but: dat$groups are now changed from grouping$groups (some former
            #ungrouped haplotypes may now be in groups)
            #So: take all haplotypes in grouping$groups (do.call) and move all
            #of those that are not in dat$groups (groupedhap) to dat$ungrouped:
            groupedhap <- union(dat$groups[[1]], dat$groups[[2]])
            dat$ungrouped <- c(dat$ungrouped,
                               setdiff(do.call(c, grouping$groups), groupedhap))
            #is new consensus a new haplotype?
            haplonr <- gethaplotypenr(dat$consensus[1,], fphaplo)
            if (!is.null(haplonr$fphaplo)) {
              fphaplo <- haplonr$fphaplo
              fphaplo.changed <- TRUE
            }
            dat$consensushaplonr[1] <- haplonr$haplotypenr
            fphapnrarr[HSfam[[HSnr]]$parentnr, 1] <- dat$consensushaplonr[1]
            fphapnrarr[HSfam[[HSnr]]$parentnr, 2] <- dat$consensushaplonr[1]
          } else {
            #both parentals match a different group
            #use the consensuses of both groups
            consensus <- rbind(combineHaplotypes(rbind(grouping$consensus[g[1],],
                                                       parhaploseq[1,])),
                               combineHaplotypes(rbind(grouping$consensus[g[2],],
                                                       parhaploseq[2,])))
            dat <- checkNadd2(consensus=consensus,
                              groups=grouping$groups[g], #g has 2 values, this is a list
                              ungrouped=grouping$ungrouped,
                              fphaplo=fphaplo)
            #due to the order of g, dat$groups 1 and 2 now match parentals 1 and 2
            #for (gr in seq_along(grouping$groups)) {
            #  if (!(gr %in% g)) dat$ungrouped <- c(dat$ungrouped,
            #                                       grouping$group[[gr]])
            #}
            #20140314:
            #aim: place all haplotypes that are not in dat$groups into ungrouped
            #but: dat$groups are now changed from grouping$groups (some former
            #ungrouped haplotypes may now be in groups)
            #So: take all haplotypes in grouping$groups (do.call) and move all
            #of those that are not in dat$groups (groupedhap) to dat$ungrouped:
            groupedhap <- union(dat$groups[[1]], dat$groups[[2]])
            dat$ungrouped <- c(dat$ungrouped,
                               setdiff(do.call(c, grouping$groups), groupedhap))
            #is new consensus a new haplotype?
            for (gr in 1:2) {
              haplonr <- gethaplotypenr(dat$consensus[gr,], fphaplo)
              if (!is.null(haplonr$fphaplo)) {
                fphaplo <- haplonr$fphaplo
                fphaplo.changed <- TRUE
              }
              dat$consensushaplonr[gr] <- haplonr$haplotypenr
            }
            #set parentals to group consensus:
            fphapnrarr[HSfam[[HSnr]]$parentnr, 1] <- dat$consensushaplonr[1]
            fphapnrarr[HSfam[[HSnr]]$parentnr, 2] <- dat$consensushaplonr[2]
          } #both parentals match a different group
        } #both parentals match a group
      } #both parentals known, each matches 0 or 1 group
      #end of different parental situations for HSsize<15, >= 2 groups
      #finally set the haplonrs of the HS individuals:
      update <- setHShaplonrs(dat, HSnr, HSfam, fphapnrarr, fphaplo)
      fphapnrarr <- update$fphapnrarr
      if (!is.null(update$fphaplo)) {
        fphaplo.changed <- TRUE
        fphaplo <- update$fphaplo
      }
    } #end of >= 2 groups, HSsize < 15
  } # >= 2 groups
  list(fphaplo.changed=fphaplo.changed,
       fphaplo=fphaplo,
       fphapnrarr=fphapnrarr,
       grouping=grouping) #for debugging
  #the caller should make sure that the global data structures are updated:
  #hapnrarr[fp,,] <- result$fphapnrarr
  #if (result$fphaplo.changed) fp_haplotypes[[fp]] <- fphaplo
} #processHSfam

calc.oldVSnew <- function(old, new, missing, no.comp=FALSE) {
  #old and new are 2-dim arrays as fphapnrarr or mhaplo$mrkallele.arr[fp,,]
  #            (but could be any pair of vectors or arrays of identical dim)
  #missing is a single value (1 for haplotypes, NA for markers)
  #        that indicates how missing values are represented in old and new
  #no.comp: if TRUE the returned vector has the same names but 5 NA's
  #return value: a vector with the frequencies of 5 different types of changes
  oldVSnew <- rep(NA, 5)
  names(oldVSnew) <- c("scored_unchanged",
                       "missing_missing",
                       "missing_scored",
                       "scored_changed",
                       "scored_missing")
  if (no.comp) return(oldVSnew)
  if (sum(dim(old) != dim(new)) > 0) {
    stop("calc.oldVSnew: old and new do not match")
  }
  missing <- missing[1]
  if (is.na(missing)) {
    # replace NA's with a (hopefully not normally occurring) value:
    missing <- -(.Machine$integer.max) + 4218063
    old[is.na(old)] <- missing
    new[is.na(new)] <- missing
  }
  oldVSnew[1] <- sum(old != missing & new != missing & old == new)
  oldVSnew[2] <- sum(old == missing & new == missing)
  oldVSnew[3] <- sum(old == missing & new != missing)
  oldVSnew[4] <- sum(old != missing & new != missing & old != new)
  oldVSnew[5] <- sum(old != missing & new == missing)
  oldVSnew
} #calc.oldVSnew

processFocalpointHS <- function(fp, fphaplo, fphapnrarr, HSfam, maxcycle=50) {
  # this function calls processHSfam repeatedly such that all HSfamilies in
  # the pedigree are processed for focalpoint fp.
  # In one cycle all families are processed in the order of HSfam (which is
  # ordered by decreasing family size). At the end of each cycle the
  # resulting hapnrarr is compared to all previous versions to detect if the
  # process has converged (final version equal to previous) or is coming
  # back to an earlier version.

  # fp: focalpoint (index number, not ID)
  # fphaplo: fp_haplotypes[[fp]], part of the object returned by
  #          get_initial_haplotypes
  # fphapnrarr: hapnrarr[fp,,], a 2-dim slice from the array of haplotype
  #             numbers (alleles) returned by get_initial_haplotypes
  # HSfam: a list as produced by get.all.HSfamilies
  # maxcycle: the maximum number of cycles
  #return value: a list with items
  # $fphaplo: possibly changed version of fp_haplotypes[[fp]] (new haplotypes added)
  # $fphaplo.changed: TRUE or FALSE; need to replace fp_haplotypes[[fp]]
  # $fphapnrarr: probably changed version of hapnrarr[fp,,]
  # $convergence: "yes", "no" or "cyclical" (no if fphapnearr scores are still
  #               changing when maxcycle is reached; cyclical if the scores keep
  #               repeating themselves in a cyclical way, yes if the scores don't
  #               change in successive cycles)
  # $cycles: the number of cycles used
  # $startVSend: integer(5): frequencies of haplotype changes

  fphaplo.changed <- FALSE
  hapnrarr_versions <- list()
  hapnrarr_versions[[1]] <- fphapnrarr
  testhapnrarr <- hapnrarr_versions
  lastcycle_thisversion <- NA
  cycle <- 0
  testn <- 1
  testversion <- matrix(c(0,0,0,testn),nrow=1)
  colnames(testversion) <- c("fp","cycle","HSnr","testn")
  repeat {
    cycle <- cycle + 1
    changes <- FALSE
    for (HSnr in seq_along(HSfam)) {
      testn <- testn + 1
      testversion <- rbind(testversion, c(fp, cycle, HSnr, testn))
      #cat(paste(fp, cycle, HSnr, testn, "\n"))
      #if (HSnr==12) browser()
      procHS <- processHSfam(fp, HSnr, fphaplo, fphapnrarr, HSfam)
      #if (sum(is.na(procHS$fphapnrarr)) > 0) browser()
      testhapnrarr[[testn]] <- procHS$fphapnrarr
      fphapnrarr <- procHS$fphapnrarr
      if (procHS$fphaplo.changed) {
        changes <- TRUE #new haplotype, non-circular changes
        fphaplo <- procHS$fphaplo
        fphaplo.changed <- TRUE
      }
    }
    if (sum(is.na(fphapnrarr)) > 0) {
      stop(paste("processFocalpoint: NA's generated in fphapnrarr in cycle",
                        cycle))
    }
    hapnrarr_versions[[cycle + 1]] <- fphapnrarr
    found <- FALSE
    if (!changes) {
      cyc <- cycle
      while (cyc >= 1 && !found) {
        diff <- sum(fphapnrarr != hapnrarr_versions[[cyc]])
        found <- diff == 0
        cyc <- cyc - 1
      }
      if (found) {
        lastcycle_thisversion <- cyc
      }
    }
    if (found || cycle >= maxcycle) break
  } #repeat

  if (is.na(lastcycle_thisversion)) {
    convergence <- "no"
  } else if ((cycle - lastcycle_thisversion) == 1) {
    convergence <- "yes"
  } else {
    convergence <- "cyclical"
  }

  if (convergence == "yes") {
    oldVSnew <- calc.oldVSnew(old=hapnrarr_versions[[1]],
                              new=fphapnrarr,
                              missing=1)
  } else {
    fphapnrarr <- hapnrarr_versions[[1]] #back to original version
    oldVSnew <- calc.oldVSnew(no.comp=TRUE)
  }

  list(fphaplo = fphaplo,
       fphaplo.changed = fphaplo.changed,
       fphapnrarr = fphapnrarr,
       convergence = convergence,
       cycles = cycle,
       oldVSnew = oldVSnew,
       hapnrarr_versions = hapnrarr_versions, #for debugging
       testversion = testversion,
       testhapnrarr = testhapnrarr)
} #processFocalpointHS

processAllFocalpoints <- function(fp_haplotypes, hapnrarr, HSfam,
                                  focalpoints=seq_along(fp_haplotypes)) {
  #processAllFocalpoints calls processFocalpointHS for each focalpoint
  #successively, but allows also to select a subset of focalpoints
  #return value: a list with elements
  # $fp_haplotypes: including any new haplotypes created
  # $hapnrarr: updated
  # $convergence: a character vector with one element per focalpoint
  #and the convergence and oldVSnew statistics.
  #By default all focalpoints are analyzed and output, but by supplying
  #an integer vector (with focalpoint numbers) or character vector (with
  #names) also a subset can be selected.
  if (is.character(focalpoints)) {
    fpmap <- getHaploblockmap(fp_haplotypes)
    focalpoints <- match(focalpoints, fpmap$marker)
  }
  focalpoints <- focalpoints[focalpoints %in% seq_along(fp_haplotypes)]
  convergence <- rep("no", length(fp_haplotypes))
  oldVSnew <- matrix(integer(0), nrow=0, ncol=5)
  messages <- character(0)
  for (fp in focalpoints) {
    #cat(paste("fp=", fp, "\n"))
    procFP <- processFocalpointHS(fp, fp_haplotypes[[fp]],
                                  hapnrarr[fp,,], HSfam)
    s <- paste("hb ", fp, " = ", names(fp_haplotypes)[fp], ": convergence = ",
               procFP$convergence, " in ", procFP$cycles, " cycles", sep="")
    convergence[fp] <- procFP$convergence
    #oldVSnew <- rbind(oldVSnew, procFP$oldVSnew)
    if (procFP$fphaplo.changed) {
      s <- paste(s, ", new haplotype(s) created", sep="")
      fp_haplotypes[[fp]] <- procFP$fphaplo
    }
    messages <- c(messages, s)
    cat(paste(s, "\n", sep=""))
    if (procFP$convergence == "yes") hapnrarr[fp,,] <- procFP$fphapnrarr
  }
  list(fp_haplotypes=fp_haplotypes,
       hapnrarr=hapnrarr,
       convergence=convergence,
       messages=messages)
} #processAllFocalpoints


allelecolors <- function(old, new, missing, convergence, usemarker=TRUE) {
  # allelecolors is the name of this function as it originally was only used
  # to determine the color code of each allele in the Pedimap files.
  # Actually these color codes are codes for how the end result of the
  # haplotyping compares with the start: Was convergence reached (or not: codes
  # 8 and 9); was the marker used (or not: code 7); were alleles for focalpoints
  # or markers originally present or not, and filled in, changed or discarded
  # (codes 0-4).
  # Now the same codes are used both for the Pedimap files and for the
  # statistics files.
  #old and new are 3-dim arrays as hapnrarr or mhaplo$mrkallele.arr
  #missing is a single value (1 for haplotypes, NA for markers)
  #        that indicates how missing values are represented in old and new
  #convergence is a character vector with for each focalpoint the
  #            convergence result (yes, no, cyclical)
  #usemarker NULL or a logical vector with length == dim(old)[1]
  #result is an 3-dim array of integers (Pedimap color codes) with the same
  #       dimensions as old and new
  nonconverged <- 9
  cyclical <- 8
  unused <- 7
  unchanged_ok <- 0 #old and new same score, not missing
  unchanged_missing <- 1 #old and new both missing score
  added <- 2    #old was missing, new is a haplotype
  changed <- 3  #old and new different haplotypes, both not missing
  removed <- 4  #old was a haplotype, new is missing
  if (sum(dim(old) != dim(new)) > 0) {
    stop("allelecolors: old and new do not match")
  }
  if (length(convergence) != dim(old)[1]) {
    stop("allelecolors: convergence doesn't match old and new")
  }
  if (length(usemarker) == 1) {
    usemarker <- seq(usemarker, length.out=dim(old)[1])
  } else if (length(usemarker) != dim(old)[1]) {
    stop("allelecolors: usemarker doesn't match old and new")
  }
  colors <- array(integer(prod(dim(old))), dim=dim(old), dimnames=dimnames(old))
  missing <- missing[1]
  if (is.na(missing)) {
    # replace NA's with a (hopefully not normally occurring) value:
    missing <- -(.Machine$integer.max) + 4218063
    old[is.na(old)] <- missing
    new[is.na(new)] <- missing
  }
  colors[new == old & new != missing] <- unchanged_ok
  colors[new == old & new == missing] <- unchanged_missing
  colors[new != old & old == missing] <- added
  colors[new != old & old != missing & new != missing] <- changed
  colors[new != old & new == missing] <- removed
  colors[convergence == "no",,] <- nonconverged
  colors[convergence == "cyclical",,] <- cyclical
  colors[!usemarker,,] <- unused
  colors
} #allelecolors

allelecolors.statistics.mrkrows <- function(colors, map, filename) {
  #The old version of this function, with 1 line per marker and 2 columns
  #per individual
  #The colorvalues calculated by function allelecolors indicate the result
  #of the haplotyping for each marker and each individual allele; here
  #we produce the complete overview of these data.
  #This function takes a 3-dim allelecolors array (colors) with dimensions
  #markers, individuals, hap (1:2) and converts it to a 2-dim array (matrix)
  #with markers in rows and 2 columns per individual (hap1 and hap2),
  #plus top and left margins with totals per color value.
  #in addition the chromosome, positions and in applicable the fp from the map
  #are in columns 2-3 or 2-4.
  #this is written to file and returned (silent)
  if (length(map$marker) != dim(colors)[1] ||
      sum(map$marker != dimnames(colors)[[1]]) > 0) {
    stop("allelecolors.markers: map doesn't match colors array")
  }
  colorvalues <- c(0:4, 7:9) #other color values are not assigned
  indcount <- dim(colors)[2]
  col2d <- matrix(nrow=dim(colors)[1], ncol=2*indcount)
  col2d[,seq(1, by=2, length.out=indcount)] <- colors[,,1]
  col2d[,seq(2, by=2, length.out=indcount)] <- colors[,,2]
  indnames <- dimnames(colors)[[2]]
  colnames(col2d)[seq(2, by=2, length.out=indcount)] <-
    paste(indnames, 2, sep="_")
  colnames(col2d)[seq(1, by=2, length.out=indcount)] <-
    paste(indnames, 1, sep="_")
  rownames(col2d) <- dimnames(colors)[[1]]
  #generate the totals per marker and per individual_hap:
  col1 <- 1+col2d #because 0 is not tabulated
  mrktable <- apply(col1, 1, tabulate, 10) #markers in columns, bins in rows
  haptable <- apply(col1, 2, tabulate, 10) #ind_hap in columns, bins in rows
  #omit entries for color non-existing color values:
  mrktable <- mrktable[colorvalues+1,, drop=FALSE]
  haptable <- haptable[colorvalues+1,, drop=FALSE ]
  #get the diagonal matrix at to left:
  topleft <- matrix(ncol=length(colorvalues), nrow=length(colorvalues))
  for (i in seq_along(colorvalues)) topleft[i, i] <- sum(mrktable[i,])
  #assemble the total matrix:
  result <- matrix(nrow=nrow(topleft) + nrow(col2d),
                   ncol=ncol(topleft) + ncol(col2d))
  result[1:nrow(topleft), 1:ncol(topleft)] <- topleft
  result[(nrow(topleft)+1):nrow(result), (ncol(topleft)+1):ncol(result)] <- col2d
  result[1:nrow(topleft), (ncol(topleft)+1):ncol(result)] <- haptable
  result[(nrow(topleft)+1):nrow(result), 1:ncol(topleft)] <- t(mrktable)
  #header <- c(names(map), colorvalues, colnames(col2d))
  #namecol <- c(colorvalues, rownames(col2d))
  nmap <- data.frame(marker=c(paste("c", colorvalues, sep=""), as.character(map$marker)),
                     chrom=c(rep(NA, length(colorvalues)), map$chrom),
                     pos=c(rep(NA, length(colorvalues)), map$pos))
  if ("hb" %in% names(map)) {
    nmap$hb <- c(rep(NA, length(colorvalues)), as.character(map$hb))
  } else names(nmap)[1] <- "haploblock"
  result <- data.frame(nmap, result)
  names(result) <- c(names(nmap), paste("c", colorvalues, sep=""), colnames(col2d))
  write.table(result, filename, col.names=TRUE, row.names=FALSE, na="",
              quote=FALSE, sep="\t")
  invisible(result)
} #allelecolors.statistics.mrkrows

allelecolors.statistics <- function(colors, map, filename) {
  #New version of the function, with one row per individual and
  #2 columns per marker in the result
  #The colorvalues calculated by function allelecolors indicate the result
  #of the haplotyping for each marker and each individual allele; here
  #we produce the complete overview of these data.
  #This function takes a 3-dim allelecolors array (colors) with dimensions
  #markers, individuals, hap (1:2) and converts it to a 2-dim array (matrix)
  #with markers in rows and 2 columns per individual (hap1 and hap2),
  #plus top and left margins with totals per color value.
  #in addition the chromosome, positions and in applicable the fp from the map
  #are in columns 2-3 or 2-4.
  #this is written to file and returned (silent)
  if (length(map$marker) != dim(colors)[1] ||
        sum(map$marker != dimnames(colors)[[1]]) > 0) {
    stop("allelecolors.markers: map doesn't match colors array")
  }
  colorvalues <- c(0:4, 7:9) #other color values are not assigned
  colornames <- paste("c", colorvalues, sep="")
  mrkcount <- dim(colors)[1]
  mrknames <- dimnames(colors)[[1]]
  indcount <- dim(colors)[2]
  indnames <- dimnames(colors)[[2]]
  col2d <- matrix(nrow=indcount, ncol=2*mrkcount)
  col2d[,seq(1, by=2, length.out=mrkcount)] <- t(colors[,,1])
  col2d[,seq(2, by=2, length.out=mrkcount)] <- t(colors[,,2])
  colnames(col2d)[seq(2, by=2, length.out=mrkcount)] <-
    paste(mrknames, "2nd", sep="_")
  colnames(col2d)[seq(1, by=2, length.out=mrkcount)] <- mrknames
  rownames(col2d) <- indnames
  #generate the totals per individual and per marker_hap:
  col1 <- 1+col2d #because 0 is not tabulated
  indtable <- apply(col1, 1, tabulate, 10) #markers in columns, bins in rows
  haptable <- apply(col1, 2, tabulate, 10) #ind_hap in columns, bins in rows
  #omit entries for color non-existing color values:
  indtable <- indtable[colorvalues+1,, drop=FALSE]
  haptable <- haptable[colorvalues+1,, drop=FALSE ]
  #get the diagonal matrix at top left:
  topleft <- matrix(ncol=length(colorvalues), nrow=length(colorvalues))
  for (i in seq_along(colorvalues)) topleft[i, i] <- sum(indtable[i,])
  #assemble the total matrix:
  result <- matrix(nrow=nrow(topleft) + nrow(col2d),
                   ncol=ncol(topleft) + ncol(col2d))
  result[1:nrow(topleft), 1:ncol(topleft)] <- topleft
  result[(nrow(topleft)+1):nrow(result), (ncol(topleft)+1):ncol(result)] <- col2d
  result[1:nrow(topleft), (ncol(topleft)+1):ncol(result)] <- haptable
  result[(nrow(topleft)+1):nrow(result), 1:ncol(topleft)] <- t(indtable)
  colnames(result) <- c(colornames, colnames(col2d))
  rownames(result) <- c(colornames, indnames)
  #include some lines for chromosome, position and possible haploblock
  #of each marker:
  hormap <- matrix("", nrow=ncol(map)-1, ncol=2*nrow(map))
  hormap[,seq(1, by=2, length.out=mrkcount)] <- t(map[,-1])
  hormap <- cbind(matrix("", ncol=ncol(topleft), nrow=nrow(hormap)), hormap)
  colnames(hormap) <- colnames(result)
  result <- rbind(hormap, result) #same column names, upper is character, lower is integer
  result <- cbind(c(names(map)[-1], colornames, indnames), result)
  colnames(result)[1] <- "name"

  write.table(result, filename, col.names=TRUE, row.names=FALSE, na="",
              quote=FALSE, sep="\t")
  invisible(result)
} #allelecolors.statistics

haplo2mrkallnrarr <- function(fp_haplotypes, hapnrarr,
                              mrkallele.arr=NULL) {
  #haplo2markerarr converts a 3-dim hapnrarr to a 3-dim array of marker
  #allele NUMBERS.
  #mrkallele.arr: if not NULL, a 3-dim array as returned by read_mhaplotypes.
  #               All markers in the haplotypes of fp_haplotypes should occur
  #               in the same order in mrkallele.arr.
  #               If this is the case, after generating a marker allele array
  #               from fp_haplotypes and hapnrarr, any extra markers from
  #               mrkallele.arr will be added.
  #return value: a list with 2 elements:
  # $mrkallnrarr: a 3-dim array of marker allele numbers; if mrkallele.arr
  #                is NULL the first dimension has only the markers in the
  #                fp_haplotypes focalpoints, else the first dimension is
  #                the same as that of mrkallele.arr
  # $usemarker: a logical vector of length dim(mrkallnrarr)[1]; if
  #             mrkallele.arr is NULL all elements are TRUE, else they are
  #             FALSE for the markers that are in mrkallele.arr but not in
  #             fp_haplotypes.
  if (length(fp_haplotypes) != dim(hapnrarr)[1]) {
    stop("haplo2markerarr: fp_haplotypes and hapnrarr don't match")
  }
  #set up the markerarr:
  mrknames <- character(0)
  for (fp in seq_along(fp_haplotypes)) {
    mrknames <- c(mrknames, colnames(fp_haplotypes[[fp]]$hapmat))
  }
  mrkallnrarr <- array(NA, dim=c(length(mrknames),dim(hapnrarr)[2:3]))
  dimnames(mrkallnrarr)[[1]] <- mrknames
  dimnames(mrkallnrarr)[[2]] <- dimnames(hapnrarr)[[2]]
  dimnames(mrkallnrarr)[[3]] <- dimnames(hapnrarr)[[3]]
  names(dimnames(mrkallnrarr))[1] <- "marker"
  names(dimnames(mrkallnrarr))[2:3] <- names(dimnames(hapnrarr))[2:3]

  #fill the mrkallnrarr:
  startmarker <- 1
  for (fp in seq_along(fp_haplotypes)) {
    endmarker <- startmarker + length(fp_haplotypes[[fp]]$markers) - 1
    for (p in 1:2) {
      hapall <- hapnrarr[fp,,p, drop=TRUE] #vector: haplo.allele for each indiv
      mrkalleles <- fp_haplotypes[[fp]]$hapmat[hapall,] #matrix: indiv in rows,
      #                                                  markers in columns
      mrkallnrarr[startmarker:endmarker,,p] <- t(mrkalleles)
    }
    startmarker <- endmarker + 1
  } # for fp
  usemarker <- rep(TRUE, dim(mrkallnrarr)[1])
  if (!is.null(mrkallele.arr)) {
    #check if all new names are in old array in same order:
    ord <- match(dimnames(mrkallnrarr)[[1]], dimnames(mrkallele.arr)[[1]])
    if (sum(is.na(ord)) > 0) {
      stop("haplo2mrkallnrarr: not all markers from fp_haplotypes are in mrkallele.arr")
    }
    if (sum(diff(ord) <= 0) > 0) {
      stop("haplo2mrkallnrarr: markers in fp_haplotypes and mrkallele.arr in different order")
    }
    usemarker <- dimnames(mrkallele.arr)[[1]] %in% dimnames(mrkallnrarr)[[1]]
    mrkallele.arr[usemarker,,] <- mrkallnrarr
    mrkallnrarr <- mrkallele.arr
  }
  list(mrkallnrarr=mrkallnrarr, usemarker=usemarker)
} #haplo2mrkallnrall

mrkallnrarr2mrkallnamearr <- function(mrkallnrarr, alleledata) {
  #mrkallnrarr2mrkallnamearr converts a 3-dim array with marker allele
  #numbers to an array with the same dimensions with marker allele names
  #mrkallnrarr: a 3-dim array with dimensions markers, individuals and parental
  mrkallnamearr <- mrkallnrarr #same dimensions and dimnames
  for (mrknr in seq_along(dimnames(mrkallnrarr)[[1]])) {
    mrkname <- dimnames(mrkallnrarr)[[1]][mrknr]
    allelenames <- alleledata$allelename[alleledata$markername == mrkname]
    mrkallnamearr[mrknr,,] <- allelenames[mrkallnrarr[mrknr,,]]
  }
  mrkallnamearr
} #mrkallnrarr2mrkallnamearr

getMarkerarray <- function(fp_haplotypes, hapnrarr, alleledata,
                           mrkallele.arr=NULL) {
  #produce a 3-dim array of marker allele NAMES from a 3-dim array of
  #focalpoint allele numbers (hapnrarr). If mrkallele.arr is not NULL
  #it is the original 3-dim marker allele (numbers) array, and may include
  #markers that were discarded before analysis (and that therefore are not
  #in fp_haplotypes)
  tmp <- haplo2mrkallnrarr(fp_haplotypes, hapnrarr, mrkallele.arr)
  list(markerarr = mrkallnrarr2mrkallnamearr(tmp$mrkallnrarr, alleledata),
       usemarker = tmp$usemarker)
} #getMarkerarray

getHaploblockmap <- function(fp_haplotypes) {
  map <- data.frame(marker=names(fp_haplotypes),
                    chrom=rep(NA,length(fp_haplotypes)),
                    pos=numeric(length(fp_haplotypes)))
  for (hb in seq_along(fp_haplotypes)) {
    map$chrom[hb] <- fp_haplotypes[[hb]]$chrom
    map$pos[hb] <- fp_haplotypes[[hb]]$pos
  }
  map
} #getHaploblockmap

getUsedMarkermap <- function(map, markernames) {
  subset(map, map$marker %in% markernames)
} ##getUsedMarkermap

get.chromnames <- function(map) {
  #map is a data frame with (at least) columns marker, chrom and pos,
  #sorted in map order
  #result is a list of unique chromosome names in order of the map
  currchrom <- "%&-1$"
  chromnames <- character(0)
  map$chrom <- as.character(map$chrom)
  for (m in 1:nrow(map)) {
    if (map$chrom[m] != currchrom) {
      chromnames <- c(chromnames, map$chrom[m])
      currchrom <- map$chrom[m]
    }
  }
  chromnames
} #get.chromnames

writePedimapFile <- function(filename, map, ped, allelearr, missing,
                             allelecolors=0) {
  #filename: a new file (overwriting existing ones) in Pedimap format.
  #map: data frame with columns (at least) marker, chrom, pos; its column
  #     markers must be exactly identical to dimnames(allelearr)[[1] (so this is
  #     a different map than the original one, if allelearr contains
  #     focalpoint haplotypes and not marker alleles)
  #ped: pedigree, possibly with traits, as before
  #allelearr: a 3-dim array with as first dimension the "markers" in the map
  #           (which are the focalpoints if allelearr is hapnrarr)
  #missing: the value used for missing allele scores in allelearr
  #         (1 for focalpoints, NA for markers)
  #allelecolors: an array of same dimensions as allelearr with the color codes
  #             for all haplotype scores, or integer of length 1
  if (sum(is.na(allelecolors) > 0)) {
    stop("writePedimapFile: allelecolors may not have missing values")
  }
  if (is.vector(allelecolors)) {
    allelecolors <- 0
  } else if (sum(dim(allelearr) != dim(allelecolors)) > 0) {
    stop("writePedimapFile: allelecolors doesn't match allelearr")
  }
  map$marker <- as.character(map$marker)
  if (nrow(map) != dim(allelearr)[1] ||
        sum(map$marker != dimnames(allelearr)[[1]]) > 0) {
    stop("writePedimapFile: map doesn't match allelearr")
  }
  ped$name <- as.character(ped$name)
  if (nrow(ped) != dim(allelearr)[2] ||
      sum(ped$name != dimnames(allelearr)[[2]]) > 0) {
    stop("writePedimapFile: ped doesn't match allelearr")
  }
  if (!is.na(missing[1])) {
    allelearr[allelearr %in% missing] <- NA
  }
  if (length(missing) == 1) missing <- rep(missing, dim(allelearr)[1])

  # write header lines:
  write("POPULATION = mhaplotypes", file=filename)
  write("UNKNOWN = - *", file=filename, append=TRUE)
  write("PLOIDY = 2", file=filename, append=TRUE)
  write("NULLHOMOZ = $", file=filename, append=TRUE)
  write("CONFIRMEDNULL = $$", file=filename, append=TRUE)

  # write pedigree and possible phenotype columns
  write("", file=filename, append=TRUE)
  write("PEDIGREE", file=filename, append=TRUE)
  write("",file=filename, append=TRUE)
  suppressWarnings( #suppress warning about appending column names
    write.table(ped, file=filename, sep="\t", quote=F,
                col.names=TRUE, row.names=FALSE, na="-", append=TRUE)
  )

  #write linkage groups
  currlocus <- 0
  chromnames <- get.chromnames(map)
  for (chr in chromnames) {
    #write chrom to file:
    chrmap <- subset(map, chrom==chr)
    if (nrow(chrmap) > 0) {
      write("", file=filename, append=TRUE)
      write(paste("LINKAGEGROUP ", chr), file=filename, append=TRUE)
      write("", file=filename, append=TRUE)
      write("MAP", file=filename, append=TRUE)
      write.table(chrmap[,c(1,3)], file=filename, sep="\t", quote=FALSE,
                  col.names=FALSE, row.names=FALSE, append=TRUE)
      for (loc in min(which(map$chrom == chr)) : max(which(map$chrom == chr))) {
        currlocus <- currlocus + 1
        write("", file=filename, append=TRUE)
        write(paste("LOCUS ", map$marker[loc]), file=filename, append=TRUE)
        allelenames <- sort(unique(as.vector(allelearr[loc,,])))
        if (!is.na(missing[currlocus]))
          allelenames <- setdiff(allelenames, missing[currlocus])
        write(c("ALLELENAMES", allelenames), file=filename, sep="\t",
              ncolumns=length(allelenames)+1, append=TRUE)
      }
    }
  } # for chr

  # write allele scores
  for (mrk in 1:nrow(map)) {
    if (length(allelecolors) == 1) {
      sc <- cbind(allelearr[mrk,,1],
                  allelearr[mrk,,2])
    } else {
      sc <- cbind(allelearr[mrk,,1],
                  allelecolors[mrk,,1],
                  allelearr[mrk,,2],
                  allelecolors[mrk,,2])
    }
    write("", file=filename, append=TRUE)
    write(paste("ALLELES ", map$marker[mrk],
                " ; chrom ", map$chrom[mrk],
                ", pos ", map$pos[mrk], sep=""),
          file=filename, append=TRUE)
    write.table(sc, file=filename, sep="\t", quote=FALSE,
                col.names=FALSE, row.names=TRUE, na="-", append=TRUE)
  } #for mrk
} #writePedimapFile

writeDatafiles <- function(map, ped, allelearr, missing, outfiles,
                         origFQparfile, usemarker=TRUE) {
  #map: data frame with columns (at least) marker, chrom, pos; its column
  #     marker must be exactly identical to dimnames(allelearr)[[1]] (so this is
  #     a different map than the original one, especially if allelearr contains
  #     focalpoint haplotypes and not marker alleles)
  #ped: as before
  #allelearr: a 3-dim array with as first dimension the "markers" in the map
  #           (which are the focalpoints if allelearr is hapnrarr)
  #           This must be rearranged to a 2-dim array
  #           with rows for individuals and 2 columns per marker.
  #missing: the value used for missing allele scores in allelearr
  #         (1 for focalpoints, NA for markers) or a vector with one symbol for
  #         each haploblock (e.g. 1^5, 1^8, 1^3, ...)
  #outfiles: the basic output filename. if origFQparfile=="" generic output
  #          is written consisting of only the
  #          phased genotypes file outfiles.dat.
  #          Else FlexQTL output is written consisting of the FQ parameter file
  #          outfiles.par, the map file outfiles.map and the data file outfiles.dat
  #          NOTE that for read.table it may be necessary to read the header lines
  #          separately, if trait and/or marker names contain invalid characters
  #          (and anyway the two columns for each marker have the same name so the
  #           second will be changed)
  #origFQparfile: the flexqtl.par used to generate the input data for this
  #               haplotyping script
  #usemarker: logical vector of length 1 or dim(allelearr)[1]: markers with
  #           usemarker=FALSE will not be written (e.g. non-analyzed markers,
  #           non-converged focalpoints)
  map$marker <- as.character(map$marker)
  if (nrow(map) != dim(allelearr)[1] ||
        sum(map$marker != dimnames(allelearr)[[1]])) {
    stop("writeDatafiles: map doesn't match allelearr")
  }
  ped$name <- as.character(ped$name)
  if (nrow(ped) != dim(allelearr)[2] ||
      sum(ped$name != dimnames(allelearr)[[2]]) > 0) {
    stop("writeDatafiles: ped doesn't match allelearr")
  }
  if (length(usemarker) == 1) usemarker <- rep(usemarker, dim(allelearr)[1])
  if (sum(is.na(usemarker)) > 0 || sum(usemarker) == 0) {
    stop("writeDatafiles: usemarker invalid")
  }
  if (length(usemarker) != dim(allelearr)[1]) {
    stop("writeDatafiles: usemarker doesn't match allelearr")
  }
  #select markers according to usemarker and replace missing-symbol by NA:
  map <- map[usemarker,, drop=FALSE]
  allelearr <- allelearr[usemarker,,, drop=FALSE]
  if (!is.na(missing[1])) {
    allelearr[allelearr %in% missing] <- NA
  }
  #make 2-dim array from allelearr:
  FQarr <- matrix(NA, nrow = dim(allelearr)[2], ncol = 2 * dim(allelearr)[1])
  for (mrk in 1:dim(allelearr)[1]) {
    FQarr[,(2*mrk-1):(2*mrk)] <- allelearr[mrk,,]
  }
  rownames(FQarr) <- dimnames(allelearr)[[2]]
  #set the column names to marker names:
  cols <- rep(dimnames(allelearr)[[1]], each=2) #duplicated names
  secondcols <- seq(2, by=2, length.out=dim(allelearr)[1])
  cols[secondcols] <- paste(cols[secondcols], "_2nd", sep="")
  colnames(FQarr) <- as.vector(cols) #? always a vector unless no markers
  #make a data frame in layout of FQ data file:
  FQ.df <- cbind(population=rep(1, nrow(ped)), ped, FQarr)
  names(FQ.df) <- c("population", names(ped), cols) #on creating data frame,
  #                             invalid characters in column names were changed,
  #                             we change back to original names
  if (origFQparfile == "") {
    #generic output:
    write.table(FQ.df[,c(2,(length(ped)+2):length(FQ.df))],
                file=paste(outfiles, ".dat", sep=""), sep="\t",
                quote=FALSE, na="", col.names=TRUE, row.names=FALSE)
  } else {
    #FlexQTL output:
    FQparfile <- paste(outfiles, ".par", sep="")
    FQmapfile <- paste(outfiles, ".map", sep="")
    FQdatafile <- paste(outfiles, ".dat", sep="")
    names(FQ.df)[1] <- paste(";", names(FQ.df)[1], sep="") #prepend a ;
    #    (without space after it) to make header a comment line for FQ
    write.table(FQ.df, file=FQdatafile, sep="\t",
                quote=FALSE, na="-", row.names=FALSE, col.names=TRUE)
    #write FQ map file:
    chromnames <- get.chromnames(map)
    mrkPerChr <- integer(length(chromnames))
    con <- file(FQmapfile, "w")
    for (chr in seq_along(chromnames)) {
      writeLines(paste("group",chromnames[chr], sep="\t"), con)
      chrmap <- map[map$chrom==chromnames[chr], c("marker","pos")]
      mrkPerChr[chr] <- nrow(chrmap)
      write.table(chrmap, con, append=TRUE, quote=FALSE,
                  sep="\t", row.names=FALSE, col.names=FALSE)
    }
    close(con)
    #next create a new flexqtl par-file for this data and map file, based on the original:

    trim.leading <- function (x)  sub("^\\s+", "", x) #removes leading whitespace

    find.second.word <- function(x) {
      i <- 1
      while (i <= nchar(x) && substr(x,i,i) %in% c(" ", "\t")) i <- i + 1
      while (i <= nchar(x) && !(substr(x,i,i) %in% c(" ", "\t"))) i <- i + 1
      while (i <= nchar(x) && substr(x,i,i) %in% c(" ", "\t")) i <- i + 1
      if (i > nchar(x)) stop("find.second.word: not found")
      i
    }

    chromnames <- get.chromnames(map)
    FQpar <- readLines(origFQparfile)
    fqp <- substr(trim.leading(FQpar), 1, 6) #6 chars from first non-blank
    i <- which(fqp == "datafi"); j <- find.second.word(FQpar[i])
    FQpar[i] <- paste(substr(FQpar[i], 1, j-1), FQdatafile, sep="")
    i <- which(fqp == "mapfil"); j <- find.second.word(FQpar[i])
    FQpar[i] <- paste(substr(FQpar[i], 1, j-1), FQmapfile, sep="")
    i <- which(fqp == "nchrom"); j <- find.second.word(FQpar[i])
    FQpar[i] <- paste(substr(FQpar[i], 1, j-1), length(chromnames), sep="")
    i <- which(fqp == "nmrkrC"); j <- find.second.word(FQpar[i])
    FQpar[i] <- paste(substr(FQpar[i], 1, j-1), paste(mrkPerChr, collapse=" "),
                      sep="")
    write(FQpar, file=FQparfile)
  }
} #writeDatafiles

calc.recombs <- function(hap1,hap2) {
  #function currently not used
  #hap1, hap2: 2 integer vectors of equal length > 0
  #result: matrix with all possible single recombinants in rows,
  #        excluding recombinants identical to one of the parental haplotypes,
  #        with NA treated as one extra allele
  if (length(hap1)==0 || length(hap2)==0 || length(hap1) != length(hap2)) {
    stop("calc.recombs: both haplotypes must have equal length > 0")
  } else {
    result <- matrix(integer(0),ncol=length(hap1))
    for (pos in 2:length(hap1)) {
      # recomb between marker (pos-1) and marker pos
      if (match.haplo(hap1[1:(pos-1)],hap2[1:(pos-1)],TRUE) ||
          match.haplo(hap1[pos:length(hap1)],hap2[pos,length(hap2)],TRUE)) {
        #on one or both sides identical subhaplotypes, no recombinants
        #that are different from one of the parental haplotypes
      } else {
        result <- rbind(result,
                        c(hap1[1:(pos-1)],hap2[pos:length(hap2)]),
                        c(hap2[1:(pos-1)],hap1[pos:length(hap1)]))
      }
    }
  }
} #calc.recombs

# Functions to re-order a pedigree data frame

sortPedigree <- function(dat, colInd, colPar1, colPar2,
                         parentsFirst=TRUE, semiFounders=FALSE) {
  #based on sortPedigree in Pedimap
  #function sorts a data frame containing pedigree info such
  #that parents always occur before (default) or always after their offspring
  #with minimal shuffling
  #dat: the data frame to be ordered
  #colInd, colPar1, colPar2: the column numbers with the names of the
  #        individuals, first and second parents. Individuals may not have
  #        missing values, but for (semi-)founders one or both parents are NA.
  #        All individuals must be different.
  #parentsFirst: if TRUE (default), parents will be sorted before any progeny,
  #              if FALSE, after all progeny
  #semiFounders: if FALSE (default) no semiFounders are allowed: parent1 and
  #              parent2 must either both be NA or both be non-NA
  #Return value: character (error message) if the input is incorrect,
  #              else the dat data frame with the rows reordered if necessary,
  #              and with the colInd, colPar1 and colPar2 as factors with
  #              the same set of levels.

  #first checks of the input (except circular pedigrees):
  if (nrow(dat) == 0) return(dat)
  dat[, colInd] <- as.character(dat[, colInd])
  dat[, colPar1] <- as.character(dat[, colPar1])
  dat[, colPar2] <- as.character(dat[, colPar2])
  if (sum(is.na(dat[, colInd])) > 0)
    return ("missing individual names not allowed")
  if (length(unique(dat[, colInd])) < nrow(dat))
    return("duplicate individual names not allowed")
  if (length(setdiff(dat[, colPar1], c(NA_character_, dat[, colInd]))) > 0 ||
        length(setdiff(dat[, colPar2], c(NA, dat[, colInd]))) > 0  )
    return("all parents should also be listed as individuals")
  if (!semiFounders &&
        sum(xor(is.na(dat[, colPar1]), is.na(dat[, colPar2]))) > 0)
    return("semi-founders not allowed")
  if (sum(dat[, colPar1] == dat[, colInd], na.rm=TRUE) > 0 ||
        sum(dat[, colPar2] == dat[, colInd], na.rm=TRUE) > 0  )
    return("some individuals are listed as their own parents")
  if (!parentsFirst) dat <- dat[nrow(dat):1,] #reverse dat:
  #we sort in parents first order, in this case starting with the reversed
  #pedigree so that if it is already sorted as parentsLast no changes will
  #be made (at the end the pedigree is again reversed for parentsLast)
  #double size of dat for sorting:
  pedsize <- nrow(dat)
  dat <- rbind(dat, dat) #double number of rows, for sorting
  for (i in 1:length(dat)) dat[(pedsize+1):nrow(dat), i] <-
    rep(NA, pedsize) #lower half empty

  #first step: place first founder at line 1
  founders <- which(is.na(dat[1:pedsize, colPar1]) &
                      is.na(dat[1:pedsize, colPar2]))
  if (length(founders) == 0)
    return ("pedigree contains no founders")
  if (pedsize == 1) {
    dat <- dat[1,]
    dat[, colInd] <- factor(dat[, colInd]) #both parents are NA
    return(dat)
  }
  if (founders[1] > 1) {
    #move all individuals from 1 to founders[1]-1 to end
    #and then move up all indivs to 1
    dat[(pedsize+1):(pedsize + founders[1] - 1),] <- dat[1:(founders[1] - 1),]
    dat[1:pedsize,] <- dat[founders[1]:(pedsize + founders[1] - 1),]
  }
  #second step: loop, placing in each pass the next individual whose parents
  #(if there are any) are already placed
  pass <- 2
  while(pass < pedsize) {
    #for sorting, pass<pedsize-1 would be sufficient, but in that case
    #the last Indiv might have a non-existing parent which would not be noticed
    #find first indiv i that can be placed at position pass,
    #i.e. whose parents (if present) are already placed
    #For all ind from position pass to end:
    #find position of parents (0 if parent NA, NA if parent not yet placed):
    posP1 <- match(dat[pass:pedsize, colPar1], dat[1:(pass-1), colInd])
    posP1[which(is.na(dat[pass:pedsize, colPar1]))] <- 0
    posP2 <- match(dat[pass:pedsize, colPar2], dat[1:(pass-1), colInd])
    posP2[which(is.na(dat[pass:pedsize, colPar2]))] <- 0
    #find the next individuals (from position pass) that could now be placed
    nxt <- which((posP1 < pass) & (posP2 < pass))
    if (length(nxt) == 0)
      return("circular references in pedigree")
    if (nxt[1] == 1) {
      #individual dat[pass,] is already in place and possibly a number of the
      #next individuals as well, find out which:
      last.ok <- max(which(nxt == 1:length(nxt)))
      #we set pass to the next individual after that:
      pass <- pass + last.ok
    } else {
      #nxt[1] > 1: the next individual that can now be placed is
      #dat[pass + nxt[1] - 1,]
      #move all individuals from pass to pass + nxt[1] - 2 to end
      #and then move up all indivs to pass
      #dat[(pedsize+1):(2*pedsize - pass - nxt[1] + 1),] <-
      dat[(pedsize + 1):(pedsize + nxt[1] - 1),] <-  #nxt[1]-1 indiv
        dat[pass:(pass + nxt[1] - 2),]               #nxt[1]-1 indiv
      dat[pass:pedsize,] <-                          #pedsize-pass+1 indiv
        dat[(pass + nxt[1] - 1):
              (pedsize + nxt[1] - 1),]                 #pedsize-pass+1 indiv
      #the individual at pass is now ok, move to next:
      pass <- pass + 1
    }
  } #while pass
  dat <- dat[1:pedsize,]
  if (!parentsFirst) dat <- dat[nrow(dat):1,] #reverse dat back again
  dat[, colInd] <- factor(dat[, colInd])
  dat[, colPar1] <- factor(dat[, colPar1], levels=levels(dat[, colInd]))
  dat[, colPar2] <- factor(dat[, colPar2], levels=levels(dat[, colInd]))
  dat
} #sortPedigree

isPedigreeSorted <- function(dat, colInd, colPar1, colPar2,
                             parentsFirst=TRUE, semiFounders=FALSE) {
  #parameters as sortPedigree
  #return value: error message if input incorrect (see sortPedigree),
  #              else TRUE or FALSE
  dat[[colInd]] <- as.character(dat[[colInd]])
  dat[[colPar1]] <- as.character(dat[[colPar1]])
  dat[[colPar2]] <- as.character(dat[[colPar2]])
  if (sum(is.na(dat[[colInd]])) > 0)
    return ("missing individual names not allowed")
  if (length(unique(dat[[colInd]])) < nrow(dat))
    return("duplicate individual names not allowed")
  if (length(setdiff(dat[[colPar1]], dat[[colInd]])) > 0 ||
        length(setdiff(dat[[colPar2]], dat[[colInd]])) > 0  )
    return("all parents should be listed as individuals")
  if (!semiFounders &&
        sum(xor(is.na(dat[[colPar1]]), is.na(dat[[colPar2]]))) > 0)
    return("semi-founders not allowed")
  if (parentsFirst) {
    parline <- match(dat[[colPar1]], dat[[colInd]])
    if (sum(parline > (1:nrow(dat)), na.rm=TRUE) > 0) return(FALSE)
    parline <- match(dat[[colPar2]], dat[[colInd]])
    if (sum(parline > (1:nrow(dat)), na.rm=TRUE) > 0) return(FALSE)
    return(TRUE)
  } else {
    parline <- match(dat[[colPar1]], dat[[colInd]])
    if (sum(parline < (1:nrow(dat)), na.rm=TRUE) > 0) return(FALSE)
    parline <- match(dat[[colPar2]], dat[[colInd]])
    if (sum(parline < (1:nrow(dat)), na.rm=TRUE) > 0) return(FALSE)
    return(TRUE)
  }
} #isPedigreeSorted

#############################################################
#some file manipulation tools

read.table.origheaders <- function(file, header=FALSE, sep="", skip=0, ...) {
  #read a data frame from a file and restore the original headers if present,
  #even if some or all are invalid identifiers
  if (!header) {
    df <- read.table(file=file, header=FALSE, sep=sep, skip=skip, ...)
  } else {
    #first read the header line:
    con <- file(file, "r", blocking = FALSE)
    headers <- readLines(con, n=skip+1)[skip+1]
    close(con)
    if (sep == "") headers <- strsplit(headers, split="\\s+")[[1]] else
      headers <- strsplit(headers, split=sep)[[1]]
    #then read the table without the header:
    df <- read.table(file, header=FALSE, sep=sep, skip=skip+1, ...)
    if (length(headers) != length(df))
      stop("read.table.origheaders: count of column headers does not match count of data columns")
    #add the original headers:
    names(df) <- headers
  }
  df
} #read.table.origheaders

testSourceTarget <- function(source, target, funcname) {
  message <- ""
  if (class(source) != "character" || length(source) != 1) {
    message <- paste(funcname, ": source must be a character of length 1")
  } else if (!file.exists(source)) {
    message <- paste(funcname, ": source '", source, "' not found", sep="")
  } else if (class(target) != "character" || length(target) != 1) {
    message <- paste(funcname, ": target must be a character of length 1")
  } else {
    suppressWarnings(if (!file.create(target))
    message <- paste(funcname, ": target '", target, "' cannot be created", sep=""))
  }
  message
} #testSourceTarget

transpose.dfr <- function(dfr) {
  #transposes a data frame such that the column names go to the first column
  #and the first column becomes the column names,
  #even if some of those column names would not be valid identifiers.
  #the name of the first column becomes "name"
  nwcolnames <- as.character(dfr[,1])
  dfr <- data.frame(name=names(dfr)[-1], t(dfr[,-1]))
  names(dfr) <- c("name", nwcolnames)
  rownames(dfr) <- NULL
  dfr
} #transpose.dfr

transpose.file <- function(source, target, sep="", skip=0,
                           na.strings="NA",
                           sep.out="\t", na.out="", ...) {
  #transposes the file source using function transpose.dfr and writes it to
  #file target.
  #sep, skip and na.strings are as for read.table; any other parameters for
  #read.table can be added optionally.
  #sep.out and na.out are the separator and the missing value symbol used in
  #the output file; the default values are those expected by haplotyping_session
  msg <- testSourceTarget(source, target, "transpose.file")
  if (msg != "") stop (msg)
  dfr <- read.table.origheaders(source, header=TRUE, sep=sep, skip=skip,
                                na.strings=na.strings, ...)
  dfr <- transpose(dfr)
  write.table(dfr, file=target, quote=FALSE, sep=sep.out, na=na.out,
              col.names=TRUE, row.names=FALSE)
  dfr
} #transpose.file

transpose <- function(source, target="", sep="", skip=0,
                           na.strings="NA",
                           sep.out="\t", na.out="", ...) {
  #if source is a data. frame, transpose.dfr is called and all other parameters
  #are ignored, else if source is a character transpose.file is called
  if (is.data.frame(source)) {
    transpose.dfr(source)
  } else if (is.character(source)) {
    transpose.file(source, target, sep, skip,
                   na.strings,
                   sep.out, na.out, ...)
  } else stop("transpose: source must be a data frame or a file name")
} #transpose


twocols2tworows.dfr <- function(dfr) {
  #dfr is a data frame that has one column of marker (or individual) names
  #and then a pair of columns for each individual (or marker). These two
  #columns store the two alleles for that individual and that marker. There is
  #one row for each marker (or individual)
  #The result is a similar data frame that now has two rows for each marker
  #(or individual) and one column for each individual (or marker).
  #Note that this is not a transposition of the data frame: if the markers
  #were in the rows and the individuals in the columns (or the other way round)
  #that does not change, only the two alleles in one marker and one individual
  #end up one above the other instead of side by side.
  #The marker (or individual) names in column 1 are duplicated (twice the same
  #name) and the first column name of each pair of columns is kept.
  if (length(dfr) %% 2 == 0) stop("twocols2tworows: number of columns should be odd")
  ncolitems <- (length(dfr) - 1) / 2
  firstrowitems <- seq(2, by=2, length.out=ncolitems)
  secondrowitems <- firstrowitems + 1;
  nwcolnames <- names(dfr)[c(1, firstrowitems)]
  #names(dfr) <- NULL
  for (i in 2:length(dfr))
    if (class(dfr[,i])=="factor") dfr[,i] <- as.character(dfr[,i])
  result1 <- dfr[, c(1, firstrowitems)]
  result2 <- dfr[, c(1, secondrowitems)]
  names(result2) <- names(result1)
  result <- rbind (result1, result2)
  rownrs <- rep(1:nrow(dfr), each=2)
  rownrs <- rownrs + c(0, nrow(dfr))
  result <- result[rownrs,]
  rownames(result) <- NULL
  result
} #twocols2tworows.dfr

twocols2tworows.file <- function(source, target, sep="", skip=0,
                           na.strings="NA",
                           sep.out="\t", na.out="", ...) {
  #changes the file source using function twocols2tworows.dfr and writes it to
  #file target.
  #sep, skip and na.strings are as for read.table; any other parameters for
  #read.table can be added optionally.
  #sep.out and na.out are the separator and the missing value symbol used in
  #the output file; the default values are those expected by haplotyping_session
  msg <- testSourceTarget(source, target, "twocols2tworows.file")
  if (msg != "") stop (msg)
  dfr <- read.table.origheaders(source, header=TRUE, sep=sep, skip=skip,
                                na.strings=na.strings, ...)
  dfr <- twocols2tworows.dfr(dfr)
  write.table(dfr, file=target, quote=FALSE, sep=sep.out, na=na.out,
              col.names=TRUE, row.names=FALSE)
} #twocols2tworows.file

twocols2tworows <- function(source, target="", sep="", skip=0,
                      na.strings="NA",
                      sep.out="\t", na.out="", ...) {
  #if source is a data.frame, twocols2tworows.dfr is called and all other parameters
  #are ignored, else if source is a character twocols2tworows.file is called
  if (is.data.frame(source)) {
    twocols2tworows.dfr(source)
  } else if (is.character(source)) {
    twocols2tworows.file(source, target, sep, skip,
                   na.strings,
                   sep.out, na.out, ...)
  } else stop("twocols2tworows: source must be a data frame or a file name")
} #twocols2tworows

tworows2twocols.dfr <- function(dfr) {
  #dfr is a data frame that has one column of marker (or individual) names
  #and then one column for each individual (or marker). There is a pair of rows
  #for each marker (or individual). These two rows store the two alleles for
  #that individual and that marker.
  #The result is a similar data frame that now has one row for each marker
  #(or individual) and two columns for each individual (or marker).
  #Note that this is not a transposition of the data frame: if the markers
  #were in the rows and the individuals in the columns (or the other way round)
  #that does not change, only the two alleles in one marker and one individual
  #end up side by side instead of one above the other.
  #The first of each pair of marker (or individual) names in column 1 is used,
  #and the column names are duplicated with "_2nd" appended to the second copy
  #(no checking if that results in duplicate column names).
  if (nrow(dfr) %% 2 == 1) stop("tworows2twocols: number of rows should be even")
  nrowitems <- nrow(dfr) / 2
  firstcolitems <- seq(1, by=2, length.out=nrowitems)
  secondcolitems <- firstcolitems + 1;
  result <- cbind(dfr[firstcolitems,], dfr[secondcolitems, -1])
  names(result) <- c(names(dfr), paste(names(dfr)[-1], "_2nd", sep=""))
  rownames(result) <- NULL
  colnrs <- rep(2:length(dfr), each=2)
  colnrs <- c(1, colnrs + c(0, length(dfr)-1))
  result[, colnrs]
} #tworows2twocols.dfr

tworows2twocols.file <- function(source, target, sep="", skip=0,
                                 na.strings="NA",
                                 sep.out="\t", na.out="", ...) {
  #changes the file source using function twocols2tworows.dfr and writes it to
  #file target.
  #sep, skip and na.strings are as for read.table; any other parameters for
  #read.table can be added optionally.
  #sep.out and na.out are the separator and the missing value symbol used in
  #the output file; the default values are those expected by haplotyping_session
  msg <- testSourceTarget(source, target, "tworows2twocols.file")
  if (msg != "") stop (msg)
  dfr <- read.table.origheaders(source, header=TRUE, sep=sep, skip=skip,
                                na.strings=na.strings, ...)
  dfr <- tworows2twocols.dfr(dfr)
  write.table(dfr, file=target, quote=FALSE, sep=sep.out, na=na.out,
              col.names=TRUE, row.names=FALSE)
} #tworows2twocols.file

tworows2twocols <- function(source, target="", sep="", skip=0,
                            na.strings="NA",
                            sep.out="\t", na.out="", ...) {
  #if source is a data.frame, tworows2twocols.dfr is called and all other parameters
  #are ignored, else if source is a character tworows2twocols.file is called
  if (is.data.frame(source)) {
    tworows2twocols.dfr(source)
  } else if (is.character(source)) {
    tworows2twocols.file(source, target, sep, skip,
                         na.strings,
                         sep.out, na.out, ...)
  } else stop("tworows2twocols: source must be a data frame or a file name")
} #tworows2twocols
