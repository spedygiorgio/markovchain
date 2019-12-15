#helper functions


#helper function for checkMP

.findNijPjk<-function(Nijk = Nijk, Nij = Nij, trans, row = 1){
  i <- Nijk[row,1]
  j <- Nijk[row,2]
  k <- Nijk[row,3]
  
  fromCh <- as.character(j)
  toCh <- as.character(k)
  Pjk <- trans[fromCh,toCh]
  
  m1 <- which(Nij[, 1] == i)
  m2 <- which(Nij[, 2] == j)
  m <- c(m1, m2)
  return(Nij[m[anyDuplicated(m)], 3] * Pjk)
}




#' @name verifyMarkovProperty
#' 
#' @rdname statisticalTests
#' @family statisticalTests
#' 
#' @title Various functions to perform statistical inference of DTMC
#' @description These functions verify the Markov property, assess 
#'              the order and stationarity of the Markov chain.
#' 
#' @param sequence An empirical sequence.
#' @param verbose Should test results be printed out?
#' @param nblocks Number of blocks.
#' 
#' @return Verification result
#' 
#' @references Anderson and Goodman.
#' 
#' @author Tae Seung Kang, Giorgio Alfredo Spedicato
#' 
#' @seealso \code{markovchain}
#' 
#' @examples 
#' sequence <- c("a", "b", "a", "a", "a", "a", "b", "a", "b",
#'               "a", "b", "a", "a", "b", "b", "b", "a")
#' mcFit <- markovchainFit(data = sequence, byrow = FALSE)
#' verifyMarkovProperty(sequence)
#' assessOrder(sequence)
#' assessStationarity(sequence, 1)
#' 
#' 
#' @export

# check if the sequence holds the Markov property
verifyMarkovProperty <- function(sequence, verbose = TRUE) {
  warning("The accuracy of the statistical inference functions has been questioned. It will be thoroughly investigated in future versions of the package.")
  #fitting the markovchain
  
  transMatrix <- markovchainFit(data = sequence)$estimate@transitionMatrix
  
  #make the (n-2)x3 matrix for observations
  subSample<-sequence[1:(length(sequence) - (length(sequence)%%3))]
  
  seqSet1<-matrix(c(subSample[1:(length(subSample) - 2)],
                    subSample[2:(length(subSample) - 1)],
                    subSample[3:(length(subSample))]
  ),ncol = 3) #fill the matrix in reverse order so position 11 is the first obersvation,12 second and 13 third
  #compute row frequencies
  temp<-as.data.frame(seqSet1)
  Nijk<-aggregate(temp, by = temp, length)[1:(ncol(temp) + 1)]
  
  seqSet2 <- seqSet1[, -3] #make matrix of couples
  temp2 <- as.data.frame(seqSet2)
  Nij <- aggregate(temp2, by = temp2, length)[1:(ncol(temp2) + 1)] #rowfrequencies included
  
  
  test<-c(length = dim(Nijk)[1])
  #compute the test statistic
  invisible(lapply(seq_len(dim(Nijk)[1]),function(i)
    {
    foundNijPjk <- .findNijPjk(Nijk = Nijk, Nij = Nij, trans = transMatrix, row = i)
    test[i] <<- ((Nijk[i,4]-foundNijPjk)^2)/foundNijPjk
  })
  )
  statistic <- sum(test)
  #return value of the test statistic and test at confience level 95% and 99%
  
  #dof
  
  dof = length(unique(sequence))^3
  
  
  pvalue <- 1-pchisq(q = statistic,df = dof)
  
  out <- list(statistic = statistic,dof = dof,p.value = pvalue)
  
  # previous version
  # n <- length(sequence)
  # u <- unique(sequence)
  # stateNames <- u
  # nelements <- length(stateNames)
  # mat <- matlab::zeros(nrow = nelements, ncol = 3)
  # 
  # # SSO: state sequence occurrences
  # # TSO: two state occurences
  # dimnames(mat) <- list(stateNames, c("SSO", "TSO", "TSO-SSO"))
  # 
  # # numeric vector initialized with zero for all states
  # SSO <- numeric()
  # for(i in 1:nelements) {
  #   sname <- stateNames[i]
  #   SSO[sname] <- 0
  # }
  # 
  # # numeric vector initialized with zero for all states
  # TSO <- SSO
  # 
  # # store the output to be returned
  # out <- list()
  # 
  # for(present in stateNames) {
  #   for(future in stateNames) {
  #     
  #     for(i in 1:nelements) {
  #       TSO[i] <- SSO[i] <- 0  
  #     }
  #     
  #     # populate TSO and SSO vector
  #     for(i in 1:(n-1)) {
  #       # let the ith state as past state
  #       past <- sequence[i]
  #       
  #       # if next state in the sequence is present state
  #       if(sequence[i+1] == present) {
  #         TSO[past] <- TSO[past] + 1
  #         
  #         # if next to next state in the sequence is future state
  #         if((i < n - 1) && (sequence[i+2] == future)) {
  #           SSO[past] <- SSO[past] + 1
  #         }
  #       }
  #     }
  #     
  #     # populate the matrix
  #     # first column corresponds to SSO, second to TSO and
  #     # third to their difference
  #     
  #     for(i in 1:(length(SSO))) {
  #       mat[i, 1] <- SSO[i]
  #       mat[i, 2] <- TSO[i]
  #       mat[i, 3] <- TSO[i] - SSO[i]
  #     }
  #     
  #   }
  # }
  #     
  #     # chi-squared test
  #     
  #     # between SSO and TSO-SSO
  #     table <- as.data.frame(mat[, c(1, 3)])
  #     
  #     # an object of class htest
  #     res <- chisq.test(table)
  #     
  #     # extract all information from htest object
  #     # and stored the result in the form of list
  #     res <- c(res)
  #     
  #     # SSO and TSO
  #     table <- as.data.frame(mat[ , c(1, 2)])
  #     
  #     # stored the table in the list
  #     res[["table"]] <- table
  #     
  #     # store the result corresponding to present state and future state
  #     out[[paste0(present, future)]] <- res
  
  if (verbose == TRUE) {
    cat("Testing markovianity property on given data sequence\n")
    cat("Chi - square statistic is:", statistic, "\n")
    cat("Degrees of freedom are:", dof, "\n")
    cat("And corresponding p-value is:", pvalue, "\n")  
  }
  
  invisible(out)
}



#' @rdname statisticalTests
#' @export

# check if sequence is of first order or of second order
assessOrder <- function(sequence, verbose = TRUE) {
  warning("The accuracy of the statistical inference functions has been questioned. It will be thoroughly investigated in future versions of the package.")
  # length of sequence
  n <- length(sequence)
  
  # unique states
  states <- unique(sequence)
  
  # number of unique states
  nelements <- length(states)
  
  TStat <- 0
  for(present in states) {
    # going to be a transition matrix 
    mat <- matlab::zeros(nelements)
    dimnames(mat) <- list(states, states)
    
    # populate transition matrix
    for(i in 1:(n - 2)) {
      if(present == sequence[i + 1]) {
        past <- sequence[i]
        future <- sequence[i + 2]
        mat[past, future] <- mat[past, future] + 1 
      }
    }
    
    # chi-squared test
    res <- chisq.test(mat)
    TStat <- TStat + res$statistic
  }
  
  k <- nelements
  df <- k * (k - 1)^2
  pvalue <- 1-pchisq(q = TStat, df)
  out <- list(statistic = TStat[[1]], p.value = pvalue[[1]])
  
  # returning the output
  if (verbose == TRUE) {
    cat("The assessOrder test statistic is: ", TStat, "\n")
    cat("The Chi-Square d.f. are: ", df, "\n")
    cat("The p-value is: ", pvalue, "\n")
  }
  
  invisible(out)
}

#' @rdname statisticalTests
#' @export

# check if sequence is stationary
assessStationarity <- function(sequence, nblocks, verbose = TRUE) {
  warning("The accuracy of the statistical inference functions has been questioned. It will be thoroughly investigated in future versions of the package.")
  # length of sequence
  n <- length(sequence)
  
  # size of each block
  blocksize <- n / nblocks
  
  # vector of unique states
  states <- unique(sequence)
  
  # number of states
  nstates <- length(states) 
  
  # sum of the statistics
  TStat <- 0 
  
  # chi-squared test for each state
  for(i in states) {
    
    # init matrix
    mat <- matlab :: zeros(nblocks, nstates)
    dimnames(mat) <- list(1:nblocks, states)
    
    # compute the transition matrix from sequence
    for(j in 1:(n - 1)) {
      if(sequence[j] == i) {
        # row index
        b <- ceiling(j / blocksize) 
        
        # next state
        future <- sequence[j+1] 
        
        # update transition matrix
        mat[b, future] <- mat[b, future] + 1 
      }
    }
    
    # vector to store row sum of matrix
    rowsums <- rowSums(mat)
    
    # store the indices with zero row sum
    indices <- which(rowsums == 0) 
    
    # update rows with zero sum
    for(k in indices) mat[k, ] <- 1/nstates 
    
    # update row sum after checking zero sum row
    rowsums <- rowSums(mat)
    
    # row-wise normalize. 
    mat <- mat/rowsums 
    
    # Some columns may still be all zeros. This causes NaN for chi-squared test.
    # chi-squared test
    res <- chisq.test(mat)
    TStat <- TStat + res$statistic
  }
  
  k <- nstates
  
  # degree of freedom
  df <- k * (nblocks - 1) * (k - 1)
  pvalue <- 1 - pchisq(TStat, df)
  
  # returning the output
  
  if (verbose==TRUE) {
    cat("The assessStationarity test statistic is: ", TStat, "\n")
    cat("The Chi-Square d.f. are: ", df, "\n")
    cat("The p-value is: ", pvalue, "\n")
  }
  
  out <- list(statistic = TStat[[1]], p.value = pvalue[[1]])
  
  invisible(out)
}

# sequence to transition frequencey matrix
.seq2mat <- function(sequence) {
  
  # basic requirement to create transition matrix
  n <- length(sequence)
  states <- unique(sequence)
  nstates <- length(states)
  
  # create transition matrix
  mat <- matlab::zeros(nstates)
  dimnames(mat) <- list(states, states)
  
  # populate transition matrix
  for(i in 1:(n - 1)) {
    from <- sequence[i]
    to <- sequence[i+1]
    mat[from, to] <- mat[from, to] + 1
  }
  
  return (mat)
}




#' @title test whether an empirical transition matrix is compatible to a theoretical one
#' 
#' @description This function tests whether an empirical transition matrix is statistically compatible
#' with a theoretical one. It is a chi-square based test
#' 
#' @rdname statisticalTests
#' @family statisticalTests
#'
#' @param data matrix, character or list to be converted in a raw transition matrix
#' @param object a markovchain object
#'
#' @return a list with following slots: statistic (the chi - square statistic), dof (degrees of freedom), and corresponding p-value
#' @export
#'
#' @examples
#' 
#' #Example taken from Kullback Kupperman Tests for Contingency Tables and Markov Chains
#' 
#' sequence<-c(0,1,2,2,1,0,0,0,0,0,0,1,2,2,2,1,0,0,1,0,0,0,0,0,0,1,1,
#' 2,0,0,2,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,2,1,0,
#' 0,2,1,0,0,0,0,0,0,1,1,1,2,2,0,0,2,1,1,1,1,2,1,1,1,1,1,1,1,1,1,0,2,
#' 0,1,1,0,0,0,1,2,2,0,0,0,0,0,0,2,2,2,1,1,1,1,0,1,1,1,1,0,0,2,1,1,
#' 0,0,0,0,0,2,2,1,1,1,1,1,2,1,2,0,0,0,1,2,2,2,0,0,0,1,1)
#' 
#' mc=matrix(c(5/8,1/4,1/8,1/4,1/2,1/4,1/4,3/8,3/8),byrow=TRUE, nrow=3)
#' rownames(mc)<-colnames(mc)<-0:2; theoreticalMc<-as(mc, "markovchain")
#' 
#' verifyEmpiricalToTheoretical(data=sequence,object=theoreticalMc)
#' 

verifyEmpiricalToTheoretical <- function(data, object, verbose = TRUE) {
  warning("The accuracy of the statistical inference functions has been questioned. It will be thoroughly investigated in future versions of the package.")  
  if (!class(object) == 'markovchain') stop("Error! Object should belong to the markovchain class")
  if (missing(data) | missing(object)) stop("Error! Required inputs missing")
  if (!(class(data) %in% c("matrix","character","numeric"))) stop("Error! Data should be either a raw transition matrix or 
                                                                  either a character or a numeric element")
  
  if (class(data) %in% c("character","numeric")) data<-createSequenceMatrix(stringchar = data)
  
  if (length(setdiff(names(data),names(object))) > 0) stop("Error! Empirical and theoretical tm have different support")
  
  # (possibly rearrange columns and rownames)
  
  data <- data[match(rownames(data),names(object)),] #matching rows
  data <- data[,match(colnames(data),names(object))] #matching cols
  
  f_i_dot <-colSums(data)
  
  statistic <- 0
  
  for (i in 1:dim(object)) {
    for (j in 1:dim(object)) {
      if (data[i, j]>0&object[i, j]>0) statistic <- statistic + data[i, j]*log(data[i, j]/(f_i_dot[i]*object[i, j]))
    }
  }
  
  statistic <- statistic * 2
  
  null_elements <- sum(object@transitionMatrix == 0)
  
  dof <- dim(object) * (dim(object) - 1) - null_elements #r(r-1) - c, c null element ob objects
  
  p.value <- 1 - pchisq(q = statistic,df = dof)
  
  if (verbose == TRUE) {
    cat("Testing whether the\n");print(data);cat("transition matrix is compatible with\n");print(object@transitionMatrix);print("theoretical transition matrix")
    cat("ChiSq statistic is",statistic,"d.o.f are",dof,"corresponding p-value is",p.value,"\n")  
  }

  #return output
  out <- list(statistic = statistic, dof = dof,pvalue = p.value)
  
  return(out)
}




.checkMatrix4Homogeneity<-function(matr) {
  out<-TRUE
  if (length(colnames(matr)) == 0) {message("Error! No colnames in input matrix"); out = FALSE}
  if (length(rownames(matr)) == 0) {message("Error! No rownames in input matrix"); out = FALSE}
  if (!all.equal(rownames(matr),colnames(matr))) {message("Error! Colnames <> Rownames")}
  if (any(matr<0)) {message("Error! Negative elements"); out = FALSE}
  return(out)
}

.addNamedColumns <- function(matr, fullnames) {
  if ( length( setdiff(names(matr),fullnames) )>0)  stop("Error! Names in matr not in fullnames")
  fullnames<-sort(fullnames)
  newMatr<-matrix(0,nrow = length(fullnames),ncol = length(fullnames),dimnames = list(fullnames,fullnames))
  
  current_support = colnames(matr)
  current_dim = dim(matr)
  
  for (i in 1:current_dim[1]) { #cycle on row
    for (j in 1:current_dim[2]) { #cycle on cols
      item<-matr[i,j] #take the element
      which_row_trans<-current_support[i] #define current row and cols
      which_col_trans<-current_support[j]
      # lookup element in the pooled table
      row_to_write <-match(x=which_row_trans,table = fullnames)
      col_to_write <-match(x=which_col_trans,table = fullnames)
      # write element into the pooled table
      newMatr[row_to_write,col_to_write] <- newMatr[row_to_write,col_to_write] + item
    }
  }

  return(newMatr)
}


#' @title Verify Homogeneity across transition matrices
#' 
#' @description Verifies that the s elements in the input list belongs to the same DTMC
#' 
#' @rdname statisticalTests
#' @family statisticalTests
#'
#' @param inputList A list of items that can coerced to transition matrices
#'
#' @return a list of transition matrices?
#' @export
#'
#' @examples
#' 
#' data(kullback)
#' verifyHomogeneity(inputList=kullback,verbose=TRUE)
#' 
verifyHomogeneity<-function(inputList, verbose = TRUE) {
  warning("The accuracy of the statistical inference functions has been questioned. It will be thoroughly investigated in future versions of the package.")  
  if (class(inputList) != "list") stop("Error! inputList should be a string")
  if (length(inputList) < 2) stop("Error! inputList length lower than 2")
  
  #checks whether all inputs can be put as transition matrices
  
  for (i in 1:length(inputList)) {
    if (is.matrix(inputList[[i]]) == TRUE) {
      checks<-.checkMatrix4Homogeneity(inputList[[i]])
      if (!checks) stop("Error! Element ", i, " to be checked")
    } else {
      inputList[[i]]<-createSequenceMatrix(stringchar = inputList[[i]]) #convert all elements into transition matrices
    }
  }
  
  # create the pooled raw transition matrix and the matrix of rowsums
  all.names<-character()
  for (i in 1:length(inputList)) {
    all.names<-c(all.names, rownames(inputList[[i]]))
  }
  all.names<-sort(unique(all.names))
  ##initialize
  PooledRawTransitionMatrix <- matrix(0,nrow = length(all.names),ncol = length(all.names),dimnames = list(all.names, all.names))
  RowSumsMatrix <- matrix(0, nrow = length(inputList),ncol=length(all.names),dimnames = list(1:length(inputList),all.names))
  ##sum for each element in the list
  
  for (k in 1:length(inputList)) {
    current_support = rownames(inputList[[k]])
    current_dim = dim(inputList[[k]])
    
    for (i in 1:current_dim[1]) { #cycle on row
      for (j in 1:current_dim[2]) { #cycle on cols
        num_trans<-inputList[[k]][i, j] #take the element
        which_row_trans <- current_support[i] #define current row and cols
        which_col_trans <- current_support[j]
        # lookup element in the pooled table
        row_to_write <-match(x = which_row_trans,table = all.names)
        col_to_write <-match(x = which_col_trans,table = all.names)
        # write element into the pooled table
        PooledRawTransitionMatrix[row_to_write,col_to_write]=PooledRawTransitionMatrix[row_to_write,col_to_write]+num_trans
      }
    }
  }
  
  #create the matrix of rowsums fij.
  for (k in 1:length(inputList)) {
    my_row_sums <- rowSums(inputList[[k]])
    current_support = names(my_row_sums)
    for (i in 1:length(current_support)) {
      my_element<-my_row_sums[i]
      col_to_write<-match(x=current_support[i],table = all.names)
      RowSumsMatrix[k, col_to_write]<-RowSumsMatrix[k, col_to_write] + my_element
    }
  }
  
  # compute the chi - square statistic
  
  statistic <- 0
 # degreesOfFreedomLess <- 0
  newInputList <- lapply(inputList, .addNamedColumns,fullnames = all.names)
  number_of_transitions <- sapply(newInputList,sum)
  total_transitions <- sum(number_of_transitions)
  
  for (s in 1:length(inputList)) { #cycle across inputs
    for (j in 1:length(all.names)) { #cycle across rows
      for (k in 1:length(all.names)) { #cycle across cols
       if (any(newInputList[[s]][j,k] == 0, number_of_transitions[s] == 0, PooledRawTransitionMatrix[j,k] == 0)) {
         statistic <- statistic + 0 # zero element in log expr does not contribute to statistics 
#         degreesOfFreedomLess <- degreesOfFreedomLess +1
       } else {
         statistic <- statistic + newInputList[[s]][j, k] * log((total_transitions*newInputList[[s]][j, k])/(number_of_transitions[s]*PooledRawTransitionMatrix[j,k]))
       }
    }
    }
  }
  
  statistic <- statistic * 2
  #dof (s-1)*(r^2-1)-#zeros
  degrees_of_freedom <- (length(inputList) - 1)*(length(all.names)^2 - 1)#-degreesOfFreedomLess
  
  p.value <- 1 - pchisq(q = statistic,df = degrees_of_freedom)
  
  if (verbose == TRUE) {
    cat("Testing homogeneity of DTMC underlying input list \n")
    cat("ChiSq statistic is",statistic,"d.o.f are",degrees_of_freedom,"corresponding p-value is",p.value,"\n")  
  }
  
  #return output
  out <- list(statistic = statistic, dof = degrees_of_freedom,pvalue = p.value)
  return(out)
}
