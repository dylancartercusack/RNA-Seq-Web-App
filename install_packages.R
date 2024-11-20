# Package names
packages <- c(
   "shiny",
   "ggplot2",
   "tidyverse",
   "corrr",
   "shinythemes",
   "plotly",
   "forcats",
   "ggrepel",
   "dplyr",
   "ggdark",
   "rsconnect",
   "ggtext"
)
# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
   install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))