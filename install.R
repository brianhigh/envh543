# Install all of the packages used in this repository.

mylib <- Sys.getenv("R_LIBS_USER")
if (!dir.exists(mylib)) dir.create(mylib, recursive = T, showWarnings = F)

if (!suppressWarnings(require("pacman"))) install.packages("pacman")

pacman::p_load(tinytex)
if (!dir.exists(tinytex_root(error = F))) tinytex::install_tinytex()

pacman::p_load(reshape, doBy, knitr, ggplot2, Hmisc, mc2d, miscTools)
pacman::p_load(deSolve, boot, cubature, lattice, foreign, openxlsx)
pacman::p_load(stats4, car, gsl)
