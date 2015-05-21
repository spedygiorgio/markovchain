# TODO: Add Display Message
# 
# Author: Giorgio
###############################################################################

# loading the markovchain package

.onLoad <- function(libname, pkgname) {
  packageStartupMessage('This is markovchain package')
}

# for unloading dynamic libraries

.onUnload <- function (libpath) {
  library.dynam.unload("markovchain", libpath)
}