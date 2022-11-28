#!/usr/bin/env Rscript

library(knitr)

args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
html <- args[2]

md <- sub(".*\\/(.*)Rmd", "\\1md", input)

knit2html(input = input, output = html)

# remove the md file that is put in project root directory
unlink(md)