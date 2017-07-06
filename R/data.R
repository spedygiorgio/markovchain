#' @name sales
#' 
#' @title Sales Demand Sequences
#' 
#' @description Sales demand sequences of five products (A, B, C, D, E).
#'              Each row corresponds to a sequence. First row corresponds to Sequence A, 
#'              Second row to Sequence B and so on.
#' 
#' @usage data("sales")
#' 
#' @details The example can be used to fit High order multivariate
#'          markov chain.
#' 
#' @examples 
#' data("sales")
#' # fitHighOrderMultivarMC(seqMat = sales, order = 2, Norm = 2)
#' 
"sales"
#' @name blanden
#' 
#' @title Mobility between income quartiles
#' 
#' @description This table show mobility between income quartiles for father and sons for the 1970 cohort born
#' 
#' @usage data(blanden)
#' 
#' @details The rows represent fathers' income quartile when the son is aged 16, whilst the columns represent sons' income quartiles when he is aged 30 (in 2000).
#' 
#' @source Personal reworking
#' 
#' @references Jo Blanden, Paul Gregg and Stephen Machin, Intergenerational Mobility in Europe and North America, Center for Economic Performances (2005)
#' 
#' @examples
#' data(blanden)
#' mobilityMc<-as(blanden, "markovchain")
"blanden"
#' @name craigsendi
#' 
#' @title CD4 cells counts on HIV Infects between zero and six month 
#' 
#' @description This is the table shown in Craig and Sendi paper showing zero and six month CD4 cells count in six brakets
#' 
#' @usage data(craigsendi)
#' 
#' @format 
#' The format is:
#' table [1:3, 1:3] 682 154 19 33 64 19 25 47 43
#' - attr(*, "dimnames")=List of 2
#'   ..$ : chr [1:3] "0-49" "50-74" "75-UP"
#'   ..$ : chr [1:3] "0-49" "50-74" "75-UP"
#'   
#' @details Rows represent counts at the beginning, cols represent counts after six months.
#' 
#' @source Estimation of the transition matrix of a discrete time Markov chain, Bruce A. Craig and Peter P. Sendi, Health Economics 11, 2002.
#' 
#' @references see source
#' 
#' @examples 
#' data(craigsendi)
#' csMc<-as(craigsendi, "markovchain")
#' steadyStates(csMc)
"craigsendi"
#' @name holson
#' 
#' @title Holson data set
#' 
#' @description A data set containing 1000 life histories trajectories and a categorical status (1,2,3) observed on eleven evenly spaced steps.
#' 
#' @usage data(holson)
#' 
#' @format 
#' A data frame with 1000 observations on the following 12 variables.
#' \describe{
#' \item{\code{id}}{unique id}
#' \item{\code{time1}}{observed status at i-th time}
#' \item{\code{time2}}{observed status at i-th time}
#' \item{\code{time3}}{observed status at i-th time}
#' \item{\code{time4}}{observed status at i-th time}
#' \item{\code{time5}}{observed status at i-th time}
#' \item{\code{time6}}{observed status at i-th time}
#' \item{\code{time7}}{observed status at i-th time}
#' \item{\code{time8}}{observed status at i-th time}
#' \item{\code{time9}}{observed status at i-th time}
#' \item{\code{time10}}{observed status at i-th time}
#' \item{\code{time11}}{observed status at i-th time}
#' }
#' 
#' @details The example can be used to fit a \code{markovchain} or a \code{markovchainList} object.
#' 
#' @source Private communications
#' 
#' @references Private communications
#' 
#' @examples 
#' data(holson)
#' head(holson)
"holson"
#' @name kullback
#' 
#' @title Example from Kullback and Kupperman Tests for Contingency Tables
#' 
#' @format A list containing two 6x6 non - negative integer matrices
#' 
#' @usage data(kullback)
#' 
#' @description A list of two matrices representing raw transitions between two states
"kullback"
#' @name rain
#' 
#' @title Alofi island daily rainfall
#' 
#' @description Rainfall measured in Alofi Island
#' 
#' @usage data(rain)
#' 
#' @format 
#' A data frame with 1096 observations on the following 2 variables.
#' \describe{
#'   \item{\code{V1}}{a numeric vector, showing original coding}
#'   \item{\code{rain}}{a character vector, showing daily rainfall millilitres brackets}
#' }
#' 
#' @source Avery Henderson
#' 
#' @references Avery Henderson, Fitting markov chain models on discrete time series such as DNA sequences
#' 
#' @examples 
#' data(rain)
#' rainMc<-markovchainFit(data=rain$rain)
"rain"
#' @name preproglucacon
#' 
#' @title Preprogluccacon DNA protein bases sequences
#' 
#' @description Sequence of bases for preproglucacon DNA protein
#' 
#' @usage data(preproglucacon)
#' 
#' @format 
#' A data frame with 1572 observations on the following 2 variables.
#' \describe{
#'   \item{\code{V1}}{a numeric vector, showing original coding}
#'   \item{\code{preproglucacon}}{a character vector, showing initial of DNA bases (Adenine, Cytosine, Guanine, Thymine)}
#' }
#' 
#' @source Avery Henderson
#' 
#' @references Averuy Henderson, Fitting markov chain models on discrete time series such as DNA sequences
#' 
#' @examples 
#' data(preproglucacon)
#' preproglucaconMc<-markovchainFit(data=preproglucacon$preproglucacon)
"preproglucacon"
#' @name tm_abs
#' 
#' @title Single Year Corporate Credit Rating Transititions
#' 
#' @description Matrix of Standard and Poor's Global Corporate Rating Transition Frequencies 2000 (NR Removed)
#' 
#' @usage data(tm_abs)
#' 
#' @format 
#' The format is:
#' num [1:8, 1:8] 17 2 0 0 0 0 0 0 1 455 ...
#' - attr(*, "dimnames")=List of 2
#' ..$ : chr [1:8] "AAA" "AA" "A" "BBB" ...
#' ..$ : chr [1:8] "AAA" "AA" "A" "BBB" ...
#' 
#' @references 
#' European Securities and Markets Authority, 2016 
#' https://cerep.esma.europa.eu/cerep-web/statistics/transitionMatrice.xhtml
#' 
#' @examples 
#' data(tm_abs)
"tm_abs"
