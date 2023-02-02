# aries

Functions for loading the
[http://www.ariesepigenomics.org.uk/](ARIES)
DNA methylation profiles.

## Accessible Resource for Integrated Epigenomic Studies (ARIES)

The Accessible Resource for Integrated Epigenomic Studies (ARIES) is a
BBSRC-funded resource of epigenomic information on a range of human
tissues, including DNA methylation data on peripheral blood at
multiple time points across the lifecourse.

> Caroline L Relton, Tom Gaunt, Wendy McArdle, Karen Ho, Aparna
> Duggirala, Hashem Shihab, Geoff Woodward, Oliver Lyttleton, David M
> Evans, Wolf Reik, Yu-Lee Paul, Gabriella Ficz, Susan E Ozanne, Anil
> Wipat, Keith Flanagan, Allyson Lister, Bastiaan T Heijmans, Susan M
> Ring, and George Davey Smith. Data Resource Profile: Accessible
> Resource for Integrated Epigenomic Studies (ARIES)
> Int. J. Epidemiol. (2015) 44 (4): 1181-1190

http://www.ariesepigenomics.org.uk/

Access to the ARIES dataset can be obtained by submitting a project proposal
to the Avon Longitudinal Study of Parents and Children (ALSPAC)
executive committee (http://www.bristol.ac.uk/alspac/researchers/access/).

Once you have a copy of the dataset, this R package can be used
to load the ARIES dataset in R.
 
## Installation

Install the `devtools` package if it is not already installed:
         install.packages("devtools")

Load `devtools` and then install `aries`.

         library(devtools)
         install_github("MRCIEU/aries")

## Examples

The following example shows how to load the methylation profiles
for blood samples collected around age 15.

         library(aries)
         aries.dir <- "path/to/aries"
         aries <- aries.select(aries.dir, time.point="15up")
         aries$meth <- aries.methylation(aries)

A more complete example including a basic EWAS
is also [available](https://mrcieu.github.io/aries/tutorial/tutorial.html).

## ARIES by number

|          | Illumina 450k| Illumina EPIC | Description |
|:---------|----:|----:|:----|
|cord      |  914|    0| Study child cord blood |
|c43m      |   85|    0| Study child age 3.5y |
|c61m      |   34|    0| Study child age 5y |
|F7        |  978|    0| Study child age 7y |
|F9        |  374|    0| Study child age 9y |
|15up      |  981| 1920| Study child age 15-17y |
|F24       |    0|  828| Study child age 24y |
|antenatal |  987|    0| Mothers during the study child pregnancy |
|FOM       |  992|    0| Mothers approximately 18 years after study pregnancy |
|FOF       |  588|    0| Fathers approximately 18 years after study pregnancy |

These numbers can be obtained from the data using following code:

         library(aries)
         aries.dir <- "path/to/aries"
         info <- aries.info(path="")
         with(info$common$samples, table(time_point,chip))


