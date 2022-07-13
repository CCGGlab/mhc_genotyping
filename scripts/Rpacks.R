# Rpacks
#########

# Install additional R packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

install.packages("devtools")
devtools::install_github("thomasp85/patchwork")
