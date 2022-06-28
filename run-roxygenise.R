#!/usr/bin/env Rscript
setwd("eigenmap")
system("rm man/*")
system("R CMD SHLIB -o src/eigenmap.so src/*.c")
roxygen2::roxygenise()
system("rm -f src/*.so")
system("rm -f src/*.o")
system("find -type f \\( -not -name \"MD5\" \\) -exec md5sum '{}' \\; > MD5")
