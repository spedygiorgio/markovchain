library("Rcpp")

# this will do compileAttributes() automatically
Rcpp.package.skeleton("markovchain", example_code = FALSE,
			cpp_files = c("1_functions4Fitting.cpp"))


