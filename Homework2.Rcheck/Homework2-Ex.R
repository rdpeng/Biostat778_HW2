pkgname <- "Homework2"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('Homework2')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("hw2_data")
### * hw2_data

flush(stderr()); flush(stdout())

### Name: hw2_data
### Title: Data description of homework2
### Aliases: hw2_data

### ** Examples

data(hw2_data)



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
