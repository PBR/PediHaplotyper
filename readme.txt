PediHaplotyper

PediHaplotyper is software for assigning haploblock alleles
to individuals in a pedigree, based on observed marker
alleles that have already been phased by other software.
A haploblock is defined as a group of very closely linked
markers among which there are (almost) no recombinations
over the pedigree. For more information see the
PediHaplotyper manual.pdf

This is an R package. To install it, use the command

devtools::install_github("PBR/PediHaplotyper")

(for this, the devtools package must be installed first;
this is available from CRAN).

This is the version published (see citation below), 
with only a few minor changes in the comments and the 
DESCRIPTION file.

Citation:
Voorrips RE, Bink MCAM, Kruisselbrink JW, 
    Koehorst - van Putten HJJ, van de Weg WE (2016) 
....PediHaplotyper: Software for consistent assignment of
    marker haplotypes in pedigrees. 
....Mol Breeding (2016) 36:119.
    DOI 10.1007/s11032-016-0539
