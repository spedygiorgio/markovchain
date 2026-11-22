# 1. Pulizia degli oggetti compilati e degli output delle vignette
devtools::clean_dll()          # rimuove src/*.o, *.so, *.dll
devtools::clean_vignettes()    # rimuove gli output delle vignette
unlink(c("inst/doc", "vignettes/.build.timestamp"), recursive = TRUE)

# 2. Rigenerazione della documentazione e ricompilazione
devtools::document()           # roxygen2 (controlla la versione, vedi sotto)
Rcpp::compileAttributes()      # se il cambio di versione ha toccato i moduli C++
devtools::load_all()           # ricompila il codice C++ da zero

# 3. Test
devtools::test()

# 4. Vignette, isolate dal resto
devtools::build_vignettes()