# 1. Scarica il pacchetto (se caricato) e cancella i compilati (.o, .dll)
devtools::unload("markovchain")
devtools::clean_dll()        # ri-esportata da pkgbuild: cancella i file compilati

# 2. (Consigliato) riavvia la sessione R: Ctrl+Shift+F10 in RStudio
#    oppure .rs.restartR()

# 3. Rigenera RcppExports, NAMESPACE e man/
devtools::document()         # chiama anche compileAttributes() per Rcpp

# 4. Ricompila e carica
devtools::load_all()         # ricompila in automatico se i sorgenti sono cambiati
# per forzare la ricompilazione completa:
pkgbuild::compile_dll(force = TRUE)

# 5. Test
devtools::test()
devtools::test(filter = "verifyMarkovProperty")   # un solo file

# 6. Check in ambiente pulito
devtools::check()                         # costruisce il tar.gz ed esegue R CMD check
devtools::check(args = "--as-cran")       # più vicino a CRAN