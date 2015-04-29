# TODO: Add Display Message
# 
# Author: Giorgio
###############################################################################



# for unloading dynamic libraries

.onUnload <- function (libpath) {
  library.dynam.unload("markovchain", libpath)
}