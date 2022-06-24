## -----------------------------------------------------------------------------
##============================================================================##
## A    Load the library
##============================================================================##
library(Rgff)

## -----------------------------------------------------------------------------
##============================================================================##
## B    Load first example data
##============================================================================##

dir <- system.file("extdata", package="Rgff")
gffFile1 <- file.path(dir,"AthSmall.gff3")


## ----message=FALSE, warning=FALSE---------------------------------------------
##============================================================================##
## C    Check the consistency and order of the GFF file 
##============================================================================##

check_gff(gffFile1)


## -----------------------------------------------------------------------------
##============================================================================##
## D    Load and check second example data
##============================================================================##

gffFile2 <- file.path(dir,"eden.gff3")

check_gff(gffFile2)


## -----------------------------------------------------------------------------

# read the first lines of "eden.gff3" file
head(read.table(gffFile2,sep="\t",header=FALSE), n=7L)

## -----------------------------------------------------------------------------
##============================================================================##
## E    Obtain the stats of the GFF file
##============================================================================##
gff_stats(gffFile1)

## -----------------------------------------------------------------------------
##============================================================================##
## F    Obtain the stats of the GFF file, disaggregated by chromosome
##============================================================================##
print(gff_stats_by_chr(gffFile1), n=50)

## -----------------------------------------------------------------------------
##============================================================================##
## G    Extract the feature organization of the GFF file as a tree
##============================================================================##
get_features(gffFile1)


## ----out.height="100%", out.width="100%", message=FALSE, warning=FALSE--------
##============================================================================##
## H    Plot the dependency tree of the GFF file
##============================================================================##

#install DiagrammeR if you do not have it installed (you need to do this only once)
# install.packages("DiagrammeR")

#load DiagrammeR
library("DiagrammeR")

#plot the features tree
plot_features(gffFile1)

## ----out.height="100%", out.width="100%"--------------------------------------
##=================================================================================##
## I  Plot the dependency tree of the GFF file in PNG format (default format)
##    and include the number of items of each feature
##=================================================================================##

plot_features(gffFile1, includeCounts = TRUE)

## ---- results=FALSE, message=FALSE, warning=FALSE-----------------------------
##=================================================================================##
## J  Plot the dependency tree of the GFF file in PDF format 
##=================================================================================##

# get the plot in a PDF file
outPlot1 <- file.path(dir,"treeplot_from_gff3_file.pdf")
plot_features(gffFile1, outPlot1, exportFormat = "pdf", includeCounts = FALSE)


## ----results=FALSE, message=FALSE, warning=FALSE, eval=FALSE------------------
#  
#  # installing and loading the required packages for svg format
#  # install.packages("DiagrammeRsvg")
#  # install.packages("rsvg")
#  
#  library("DiagrammeRsvg")
#  library("rsvg")
#  
#  # get the plot in a svg file
#  outPlot2 <- file.path(dir,"outplot_from_gff3.svg")
#  plot_features(gffFile1, outPlot2, exportFormat = "svg", includeCounts = TRUE)
#  

## ----message=FALSE, warning=FALSE---------------------------------------------
##============================================================================##
## K    Extract the feature organization of the GFF file in data.frame format
##============================================================================##
get_features(gffFile1, outFormat = 'data.frame', includeCounts = TRUE)


## ----message=FALSE, warning=FALSE---------------------------------------------
##=================================================================================##
## L    Extract the feature organization of the GFF as JSON
##=================================================================================##
gffFile1_json_features <- get_features(gffFile1, outFormat = 'JSON')
strsplit(gffFile1_json_features,"\\n");


## ----message=FALSE, warning=FALSE---------------------------------------------
##=================================================================================##
## M  Sort an unsorted GFF file
##=================================================================================##

#sorts the unsorted file gffFile2 (eden.gff3)
gffFile2_sorted <- sort_gff(gffFile2)

# check if the sorted file is well-formatted
check_gff(gffFile2_sorted)

# let's take a look to the sorted file
head(read.table(gffFile2_sorted,sep="\t"), n=10L)


## ----message=FALSE, warning=FALSE---------------------------------------------

##============================================================================##
## N    Convert a GFF file to SAF format, only the "gene" feature
##============================================================================##

safFileConverted <- saf_from_gff(gffFile1, features = c("gene"))

read.table(safFileConverted,sep="\t",header=TRUE)


## ----message=FALSE, warning=FALSE---------------------------------------------

##================================================================================##
## O    Convert a GFF file to SAF format, both "gene" and "ncRNA_gene" features
##================================================================================##

safFileConverted2 <- saf_from_gff(gffFile1, features = c("gene","ncRNA_gene"))

read.table(safFileConverted2,sep="\t",header=TRUE)


## ----message=FALSE, warning=FALSE---------------------------------------------
##============================================================================##
## P    Convert a GFF file to SAF format, compiling "exons by gene"
##============================================================================##

safFileConverted3 <- saf_from_gff(gffFile1, features = c("gene > exon"))

read.table(safFileConverted3,sep="\t",header=TRUE)


## ----message=FALSE, warning=FALSE---------------------------------------------

safFileConverted4 <- saf_from_gff(gffFile1)


## ----message=FALSE, warning=FALSE---------------------------------------------

safFileConverted5 <- saf_from_gff(gffFile1, features = c("gene : exon"), sep = ':')


## ----message=FALSE, warning=FALSE---------------------------------------------

##==============================================================================##
## Q    Convert a GFF file to SAF format, compiling "exons by gene"  
##      and "exons by non-coding RNA genes"
##==============================================================================##

safFileConverted6 <- saf_from_gff(gffFile1, features = c("gene > exon","ncRNA_gene > exon"))

read.table(safFileConverted6,sep="\t",header=TRUE)


## ----message=FALSE, warning=FALSE---------------------------------------------
##============================================================================##
## R    Convert from GTF to GFF3
##============================================================================##

# load and show our example GTF file
gtfFile1 <- file.path(dir,"AthSmall.gtf")
head(read.table(gtfFile1,sep="\t"))


## ----echo=TRUE, results='hide', message=FALSE, warning=FALSE------------------

# convert the GTF format to GFF3 format
gffFileConverted <- gtf_to_gff3(gtfFile1, forceOverwrite = TRUE)



## ----message=FALSE, warning=FALSE---------------------------------------------
# show the results of the conversion
head(read.table(gffFileConverted,sep="\t"))


