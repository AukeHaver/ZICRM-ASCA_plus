### ZICRM-ASCA+ supplementary Functions
### Author: Auke Haver
### BDA GROUP SILS Amsterdam
### Project: Zero-inflated GLMM ASCA

library("tidyverse") # General IO
library("ggplot2")
library("ggpubr") # Extension on ggplot2
library("glmmTMB") # zero-inflated GLMM fitting
library("parallel") # Multi-core processing
library("writexl")

# Virids color objects 
dark_purple <- rgb(72/255,40/255,120/255)
teal <- rgb(38/255,130/255,142/255)