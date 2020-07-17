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
is also [available](http://htmlpreview.github.io/?https://github.com/MRCIEU/aries/docs/tutorial/tutorial.html).

