suppressPackageStartupMessages(library(tools, quietly = TRUE)) # text formatting

# format data-y field names into printing english
MakeTitleCase <- function(st) {
  st <- gsub("_"," ", st)
  st <- gsub("\\."," ", st)
  st <- toTitleCase(st)
  return(st)
}