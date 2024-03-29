### Training Day 3 ###
## By Cole Nawrocki ## 

## Topic: Packages and package managers

## CRAN
# The packages that you install via `install.packages()` come from CRAN by default. CRAN is the de facto 
# repository for R packages. Basically, a packages gets reviewed by developers before it can be uploaded 
# to CRAN. The advantage here is that if a package is on CRAN, then it is "stable." This means that it will 
# download on any device running R without a problem. The disadvantage is that it takes forever for CRAN 
# packages to be updated. So, if the people who made the package decide to change something, then it won't be 
# in the CRAN download of their package for a while. Sometimes, a package gets updated to depend on the newer, 
# non-CRAN version of a package. So, to keep using the functions from the updated package, you have two 
# options: 1) Uninstall packages and re-install the specific versions that are compatible or 2) Install the 
# non-CRAN version of the package. In my experience, option 2 is best. There are many ways to do this, but two 
# simple ways that work for many packages is to download from the "r-project" or "r-universe" repositories or 
# to download from GitHub. Note: anyone can upload a package to GitHub without it being reviewed/vetted, so be 
# careful--who knows what mistakes exist in people's code? To do these things: 

# r-project
install.packages("BiocManager", repos = "https://cloud.r-project.org")

# GitHub
install.packages("remotes")
remotes::install_github("mojaveazure/seurat-disk") # format is: "user-name/repository-name"

# Sometimes a package is only available from one of these other installations.

## Source vs. Binary and C Compilers
# Furthermore, it is my understanding that `install.packages()` defaults to download the binary version of
# a package. The advantage here is that you do not need a C compiler for this. The disadvantage here is 
# that this is code that has been translated from its original form. So, sometimes, this can cause problems, 
# for reasons that I do not understand. Often, to solve these problems, you simply need to get a C compiler 
# and then install from source. R and Python are based on C. That is, they themselves and many of the 
# packages that they run were originally written in C. So a C compiler allows you to download the source
# code, which is C, and then "compile" it into R. There are two popular C compilers for linux-based systems
# as far as I know. They are gcc and clang. Once you have one or both of them, you can install packages 
# from source. 

## Package Managers
# But, how do you install these? These are not R packages; they are modules for your system. So, you use a 
# package manager. Popular package managers are homebrew, apt-get, conda, and pip. Homebrew is for MacOS. 
# apt-get is for Ubuntu. Conda and pip are for python packages, mostly. Let's get homebrew by following these
# instructions: https://brew.sh/. Then you use it to install modules from your terminal like this: 
# $ brew install wget

## Getting the gcc C Compiler and Installling from Source
# $ brew install gcc
# To test, try the following: 
install.packages("irlba", type="source") 
# If it doesn't work, then we will have to get technical. Likely, the "Makeconf" file in your R directory is
# pointing to the wrong compiler or to the wrong place on your machine. If this is the case, open your 
# "Makeconf" located here: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/etc/Makeconf, which is 
# just a txt file. Find the FLIBS variable and change it to the following: 
# FLIBS =  -L/usr/local/Cellar/gcc/13.2.0/lib/gcc/13/gcc/x86_64-apple-darwin23/13 -L/usr/local/Cellar/gcc/13.2.0/lib/gcc/13 -lgfortran -lquadmath -lm
# This should be right if you installed gcc from homebrew. Note: you may have to correct the versions in the 
# line above. To do this, navigate to /usr/local/Cellar/gcc and follow the directories to find your version. 
# Then, correct the line above. Now, R should be able to locate the compiler. Try again to install a package 
# from source:  
install.packages("irlba", type="source") 

## Bioconductor 
# Bioconductor is essentially CRAN but for Bioinformatics packages. It comes with its own package installer 
# that you can use to install the official, reviewed packages from the Bioconductor repository. You installed
# Bioconductor from the r-project above. Now, we can use it to install the scDblFinder package that we were 
# having trouble installing in the last training script. 
BiocManager::install("scDblFinder")

# In fact, when I first tried to use the package, I ran into the issue that the irlba package was not 
# compatible anymore with the scDblFinder package. To fix it, I had to re-install the irlba package from 
# source, which we just did above. Now, the scDblFinder package ought to work correctly. Go ahead and 
# test it out.

## A note on this process
# This exact process above won't work for every single scenario. You may have to install a package from
# a different package manager. Or, you may need the clang compiler. The important lesson here is that you do 
# not have to get freaked out when the package you want is failing to install (usually it says "package had
# non-zero exit status" or something like that). Just find the line where it says, "can't find the following
# directory" or "this package requires a different version of the following dependency" and then use your 
# package manager knowledge and your googling skills to find the right command line calls to copy and 
# paste into your terminal. Pro tip: Stack exchange and GitHub issue forums are SUPER HELPFUL for this stuff. 
# Also, developers tend to get antsy when the packages they made do not install correctlt and they will often 
# help you or other people super quickly when they get messaged. So... message them.  

