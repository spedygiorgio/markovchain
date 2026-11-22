# Cambi strutturali nei download di {markovchain} - v2
# - serie storica lunga (dal 2014) per stimare la stagionalita' annuale (STL)
# - giorni con 0 download trattati come log mancanti di cranlogs (non come zeri veri)
# - spike isolati smorzati (clipping dei residui a 3 MAD) prima di cercare i break
# - break multipli (Bai-Perron) sulla serie destagionalizzata:
#     A) storia completa  B) finestra dal 2025-07-01
# - confronto con i download dei reverse dependencies
# Uso: Rscript analisi-download-markovchain-v2.R

pkgs <- c("cranlogs", "strucchange", "ggplot2", "jsonlite")
miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) install.packages(miss)
suppressPackageStartupMessages({
  library(cranlogs); library(strucchange); library(ggplot2)
})

pkg        <- "markovchain"
hist_from  <- as.Date("2014-01-01")   # serie usata per stimare la stagionalita'
win_from   <- as.Date("2025-07-01")   # finestra di interesse
to         <- Sys.Date() - 1
DO_REVDEPS <- TRUE

## 1. Download giornalieri ----------------------------------------------------
d <- cran_downloads(packages = pkg, from = hist_from, to = to)
d$date <- as.Date(d$date)
stopifnot(nrow(d) > 0)
write.csv(d, "markovchain_downloads_daily_full.csv", row.names = FALSE)

# cranlogs restituisce 0 (non NA) quando mancano i log: li trattiamo come mancanti
out <- d[d$count == 0, ]
cat(sprintf("Giorni con 0 download (probabili buchi nei log): %d\n", nrow(out)))
if (nrow(out)) print(head(out[, c("date", "count")], 30))

## 2. Serie settimanale (lun-dom, solo settimane di calendario complete) --------
d$week <- as.Date(cut(d$date, "week", start.on.monday = TRUE))
d$ok   <- d$count > 0
wk <- data.frame(week = sort(unique(d$week)))
wk$n     <- as.vector(tapply(d$count, d$week, length))
wk$n_ok  <- as.vector(tapply(d$ok,    d$week, sum))
wk$count <- as.vector(tapply(d$count, d$week, sum))
wk <- wk[wk$n == 7, ]                   # scarta prima/ultima settimana parziale
wk$ok <- wk$n_ok == 7                   # settimana senza buchi nei log
wk$y  <- ifelse(wk$ok, log1p(wk$count), NA_real_)
cat(sprintf("Settimane: %d (%s -> %s); con buchi nei log: %d\n",
            nrow(wk), min(wk$week), max(wk$week), sum(!wk$ok)))

## 3. Destagionalizzazione STL (periodo ~52 settimane) ---------------------------
idx   <- seq_len(nrow(wk))
y_fil <- approx(idx[wk$ok], wk$y[wk$ok], xout = idx, rule = 2)$y   # solo per STL
st    <- stl(ts(y_fil, frequency = 52), s.window = "periodic", robust = TRUE)
wk$seas  <- as.numeric(st$time.series[, "seasonal"])
wk$trend <- as.numeric(st$time.series[, "trend"])
wk$rem   <- as.numeric(st$time.series[, "remainder"])
cl <- 3 * mad(wk$rem[wk$ok])
wk$adj <- wk$trend + pmax(pmin(wk$rem, cl), -cl)   # destagionalizzata, spike smorzati
wk$adj[!wk$ok] <- wk$trend[!wk$ok]                 # settimane con buchi: solo trend
cat(sprintf("Ampiezza stagionale (max/min): x%.2f ; spike smorzati: %d settimane\n",
            exp(diff(range(wk$seas))), sum(abs(wk$rem) > cl & wk$ok)))
cat("Attenzione: la stagionalita' e' stimata come costante su tutta la storia.\n")

## 4. Break multipli ---------------------------------------------------------------
find_breaks <- function(x, dates, h, maxbr) {
  n  <- length(x)
  mb <- min(maxbr, floor(n / h) - 1)
  bp <- breakpoints(x ~ 1, h = h, breaks = mb)
  print(summary(bp))
  b  <- breakpoints(bp)$breakpoints
  b  <- b[!is.na(b)]
  seg <- cut(seq_len(n), c(0, b, n), labels = FALSE)
  fit <- expm1(ave(x, seg))
  reg <- do.call(rbind, lapply(sort(unique(seg)), function(s) {
    i <- which(seg == s)
    data.frame(regime = s, dal = dates[min(i)], al = dates[max(i)],
               settimane = length(i), livello_sett = round(expm1(mean(x[i]))))
  }))
  ci <- NULL
  if (length(b)) {
    ci <- tryCatch(confint(bp, breaks = length(b))$confint, error = function(e) NULL)
    if (!is.null(ci)) {
      ci <- data.frame(nuovo_regime_dal = dates[b + 1],
                       ic95_da = dates[pmin(pmax(ci[, 1], 1), n - 1) + 1],
                       ic95_a  = dates[pmin(pmax(ci[, 3], 1), n - 1) + 1])
    } else ci <- data.frame(nuovo_regime_dal = dates[b + 1])
  }
  list(b = b, fit = fit, regimi = reg, break_ci = ci)
}

cat("\n==== A) Storia completa (serie destagionalizzata, segmenti >= 26 settimane) ====\n")
A <- find_breaks(wk$adj, wk$week, h = 26, maxbr = 10)
print(A$break_ci); print(A$regimi)
wk$fitA <- A$fit

cat("\n==== B) Finestra dal", format(win_from), "(segmenti >= 8 settimane) ====\n")
w <- wk[wk$week >= win_from, ]
B <- find_breaks(w$adj, w$week, h = 8, maxbr = 5)
print(B$break_ci); print(B$regimi)
w$fitB <- B$fit

# robustezza (opzionale)
if (requireNamespace("changepoint", quietly = TRUE)) {
  cp <- changepoint::cpt.mean(w$adj, method = "PELT", penalty = "MBIC")
  cat("PELT (MBIC) sulla finestra, nuovi regimi dal:\n")
  print(w$week[changepoint::cpts(cp) + 1])
}

## 5. Release CRAN ---------------------------------------------------------------
rel <- tryCatch({
  tl <- jsonlite::fromJSON(sprintf("https://crandb.r-pkg.org/%s/all", pkg))$timeline
  r  <- data.frame(version = names(tl), date = as.Date(unlist(tl)))
  r[r$date >= hist_from, ]
}, error = function(e) NULL)
if (!is.null(rel)) { cat("\nRelease:\n"); print(tail(rel, 12)) }

## 6. Reverse dependencies (opzionale) -------------------------------------------
rev_w <- NULL
if (DO_REVDEPS) {
  rd <- tryCatch({
    db <- available.packages(repos = "https://cloud.r-project.org")
    tools::package_dependencies(pkg, db = db, reverse = TRUE,
                                which = c("Depends", "Imports", "LinkingTo"))[[pkg]]
  }, error = function(e) NULL)
  if (length(rd)) {
    cat(sprintf("\nReverse dependencies (Depends/Imports/LinkingTo): %d\n", length(rd)))
    rdd <- cran_downloads(packages = rd, from = win_from, to = to)
    rdd$date <- as.Date(rdd$date)
    tot <- aggregate(count ~ package, rdd, sum)
    print(head(tot[order(-tot$count), ], 10), row.names = FALSE)
    rdd$week <- as.Date(cut(rdd$date, "week", start.on.monday = TRUE))
    rev_w <- aggregate(count ~ week, rdd, sum)
    names(rev_w)[2] <- "revdeps"
    m <- merge(w[, c("week", "count")], rev_w, by = "week")
    cat(sprintf("Correlazione settimanale (log) markovchain vs somma revdeps: %.2f\n",
                cor(log1p(m$count), log1p(m$revdeps))))
  }
}
# Per capire CHI scarica (CI/bot): i log grezzi di Posit (cran-logs.rstudio.com)
# contengono ip_id, r_version, os per ogni download; non e' incluso qui.

## 7. Grafici -----------------------------------------------------------------------
vl_rel <- function(p, xr) {
  if (!is.null(rel) && nrow(rel[rel$date >= xr[1] & rel$date <= xr[2], ]))
    p + geom_vline(data = rel, aes(xintercept = as.numeric(date)),
                   colour = "darkgreen", alpha = 0.5) +
      geom_text(data = rel[rel$date >= xr[1] & rel$date <= xr[2], ],
                aes(date, Inf, label = version), angle = 90, vjust = -0.3,
                hjust = 1.1, size = 3, colour = "darkgreen")
  else p
}

pA <- ggplot(wk, aes(week + 3.5)) +
  geom_line(aes(y = count), colour = "grey70") +
  geom_line(aes(y = expm1(adj)), colour = "steelblue") +
  geom_step(aes(y = fitA), colour = "firebrick", linewidth = 1, direction = "mid") +
  geom_vline(xintercept = as.numeric(wk$week[A$b + 1]), linetype = 2) +
  labs(title = sprintf("{%s}: download settimanali, storia completa", pkg),
       subtitle = "Grigio: osservato / Blu: destagionalizzato / Rosso: livello per regime",
       x = NULL, y = "download a settimana") + theme_minimal()
pA <- vl_rel(pA, range(wk$week))
ggsave("markovchain_A_storia_completa.png", pA, width = 11, height = 5, dpi = 150)

pB <- ggplot(w, aes(week + 3.5)) +
  geom_line(aes(y = count), colour = "grey60") +
  geom_line(aes(y = expm1(adj)), colour = "steelblue") +
  geom_step(aes(y = fitB), colour = "firebrick", linewidth = 1, direction = "mid") +
  geom_vline(xintercept = as.numeric(w$week[B$b + 1]), linetype = 2) +
  labs(title = sprintf("{%s}: dal %s, serie destagionalizzata", pkg, format(win_from)),
       subtitle = "Grigio: osservato / Blu: destagionalizzato / Rosso: livello per regime",
       x = NULL, y = "download a settimana") + theme_minimal()
pB <- vl_rel(pB, range(w$week))
ggsave("markovchain_B_finestra.png", pB, width = 10, height = 5, dpi = 150)

if (!is.null(rev_w)) {
  m2 <- merge(w[, c("week", "count")], rev_w, by = "week")
  pC <- ggplot(m2, aes(week + 3.5)) +
    geom_line(aes(y = count, colour = "markovchain")) +
    geom_line(aes(y = revdeps, colour = "somma reverse deps")) +
    scale_y_log10() + labs(x = NULL, y = "download a settimana (log)", colour = NULL,
                           title = "markovchain vs reverse dependencies") + theme_minimal()
  ggsave("markovchain_C_revdeps.png", pC, width = 10, height = 4.5, dpi = 150)
}
print(pA); print(pB)
