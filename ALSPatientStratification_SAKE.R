#####################################  SAKE GUI Install for Patient Stratification  ########################################################

#Notes written by: Jarrett Eshima
#Date: October, 2020
#For: Use by the Dr. Barbara Smith Lab at Arizona State University

#############################################  Installation INSTRUCTIONS   #################################################################
setRepositories() #Run this command and copy/paste the following into the console prompt below: 1 2 3 4 5

#Steps 1-6 only need to be completed once (per computing environment) before running step 7.
#Note: Be extremely careful editing the System Environment (step 2), any pasted paths are added to the end of the environment path, without the ability to remove it.
#Make sure to copy and paste the old PATH to a .txt file to save any potentially R breaking mistakes.


#Step 1 - Download sake.master.zip from https://genome.cshlp.org/content/28/9/1353.short (supplemental files)
    #Extract zip files to known folder
    #Open the DESCRIPTION file in the main Sake folder using .txt editor (notepad or similar program)
    #Edit the second line under the "Remotes:" heading
        #Replace the line with the following: bioc::release/biomaRt,


#Step 2 - Rtools40 install
#https://cran.r-project.org/bin/windows/Rtools/
writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron") #DO NOT DELETE .Renviron file from documents folder
Sys.which("make")
#Console Output: "C:\\PROGRA~1\\rtools40\\usr\\bin\\make.exe"
install.packages("jsonlite", type = "source") #This should work, if not delete the Rtools folder and reinstall


#Step 3 - Download dependent packages
setRepositories() #make sure 1-5 are selected or else the next line will fail
install.packages(c("annotate", "AnnotationHub", "biomaRt", "DESeq2", "gage", "gageData", "GO.db", "pathview", "plotly", "DT"))
install.packages("devtools")
library(devtools)
#Optional
install.packages("githubinstall")


#Step 4 - Pkgmaker Installation - both lines below are equivalent
githubinstall::gh_install_packages('renozao/pkgmaker', ref = 'develop')
devtools::install_github('renozao/pkgmaker')
library(pkgmaker)


#Step 5 - SAKE Installation - both lines are equivalent
library(githubinstall)
githubinstall::gh_install_packages("naikai/sake")
devtools::install_github('naikai/sake')


#Step 6 - Roll back package versions for: plotly, shiny, shinydashboard, htmltools, htmlwidgets, promises, crosstalk, and DT
remove.packages("shinythemes", lib="~/R/win-library/4.0") #Do this for all packages listed above
#Package Versions for Original Working Software:  shiny 1.5.0, shinydashboard 0.6.1, htmltools 0.5.0, htmlwidgets 1.5.2, promises 1.1.1, crosstalk 1.1.0.1, plotly 4.9.2, and DT 0.16
#Resource: https://www.rdocumentation.org/packages/DT/versions/0.16
sessionInfo() #can be used to check package versions
packageVersion("shinyBS") #this works too
getOption("repos") #List of options for repository links

library(devtools)
install_version("shinythemes", version = "1.1.2", repos = "https://cran.rstudio.com/")


##########################################  After Successful Installation  #####################################################


#Step 7 - RUN SAKE

library(sake)
shiny::runApp(system.file("sake", package="sake"))



################################################################################################################################
#Save the SAKE NMF output (zip) from each run for robust subtype assignment
#Save the full SAKE NMF feature table for feature selection