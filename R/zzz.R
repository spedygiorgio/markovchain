# TODO: Add Display Message
# 
# Author: Giorgio
###############################################################################

# loading the markovchain package

.onAttach <- function(libname, pkgname) {
  desc  <- packageDescription(pkgname, libname)
  packageStartupMessage('Package:  ', desc$Package, '\n',
                        'Version:  ', desc$Version, '\n', 
                        'Date:     ', desc$Date, '\n',
                        'BugReport: ', desc$BugReports, '\n')
}

# for unloading dynamic libraries

.onUnload <- function (libpath) {
  library.dynam.unload("markovchain", libpath)
}