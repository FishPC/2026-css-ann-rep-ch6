# call/install packages ---------------------------------------------------
## package list
packFontHandler <- function(){
  packages <- c(
    'ggplot2',
    'flextable',
    'officer',
    'extrafont',
    'dplyr',
    'metafor',
    'esc',
    'dmetar',
    'tidyverse',
    'boot',
    'rms',
    'plotly',
    'scatterplot3d',
    'qpcR',
    'MuMIn',
    'parallel',
    'openair',
    'EnvStats',
    'equatags',
    'js',
    'ggfun',
    'processx',
    'reticulate',
    'ggpubr',
    'ggmagnify',
    'rjags',
    'tidyr',
    'splitstackshape',
    'janitor',
    'lubridate',
    'gginnards'
  )
  
  ## install or load packages
  if (!require(install.load)) {
    install.packages('install.load')
  }
  
  install.load::install_load(packages)
  
  devtools::install_github('ropensci/plotly', upgrade = "never")
  
  ## dmetar must be installed from github
  if (!require('remotes')) {
    install.packages('remotes')
  }
  remotes::install_github('MathiasHarrer/dmetar')
  remotes::install_github("hughjonesd/ggmagnify")
  
  eval(metafor:::.MuMIn,envir = globalenv())
  
  ## installs for exporting static Plotly figures
  ### uncomment to install miniconda if not already on machine
  # reticulate::install_miniconda(force = TRUE)
  # reticulate::conda_install('r-reticulate', 'python-kaleido==0.1.*')
  # reticulate::conda_install('r-reticulate', 'plotly', channel = 'plotly')
  # reticulate::use_miniconda('r-reticulate')
  
# reconcile fonts ---------------------------------------------------------
  remotes::install_version('Rttf2pt1', version = '1.3.8')
  font_import(prompt = FALSE, pattern = 'calibri')
  fonts()
  loadfonts(device = 'win')
  windowsFonts()
}