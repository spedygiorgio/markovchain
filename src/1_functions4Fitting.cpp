// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include <ctime>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace RcppArmadillo;
using namespace RcppParallel;
using namespace std;

#include "mathHelperFunctions.h"
#include "mapFitFunctionsSAI.h"
#include "sampler.h"   
#include <math.h>


// [[Rcpp::export(.markovchainSequenceRcpp)]]
CharacterVector markovchainSequenceRcpp(int n, S4 markovchain, CharacterVector t0,
                                        bool include_t0 = false) {
  
  // character vector to store the result
  CharacterVector chain(n);
  
  // transition mastrix
  NumericMatrix transitionMatrix = markovchain.slot("transitionMatrix");
  
  // possible states
  CharacterVector states = markovchain.slot("states");
  
  // current state
  CharacterVector state = t0;
  
  NumericVector rowProbs(states.size());
  CharacterVector outstate;
  
  
  
  for(int i = 0;i < n;i++) {
    
    // extracting row probabilties for the given state from transition matrix
    int row_no = 0;
    for(int j = 0;j < states.size();j++) {
      
      /* last element of state character vector
         because of markovchainListRcpp, a seq of length greater than 1
         is also passed whose end state is the beginning state here */
      
      if(states[j] == state[state.size()-1]) {
        row_no = j;
        break;
      }
    }
    
    for(int j = 0; j < states.size(); j++) {
      rowProbs[j] = transitionMatrix(row_no, j);
    }
    
    // calculate next state
    outstate = sample(states, 1, false, rowProbs);
    chain[i] = outstate[0];  
    state = outstate;
    
  }
  
  if (include_t0)
    chain.push_front(t0[0]);
  
  return chain;
}

bool checkSequenceRcpp(List object) {
  bool out = true;
  int nob = object.size();
  
  // if there is only one markovchain object return true
  if (nob == 1)
    return(true);
  
  S4 ob0, ob1;
  CharacterVector statesN1, statesN, intersection;
  
  for(int i = 1; i < nob;i++) {
    
    CharacterVector statesNm1;
    
    ob0 = S4(object[i-1]);
    ob1 = S4(object[i]);
    
    statesN1 = ob0.slot("states"); 
    statesN = ob1.slot("states");
    
    NumericMatrix matr = ob0.slot("transitionMatrix");
    double csum = 0;
    
    for(int j = 0;j < matr.ncol();j++) {
      csum = 0;
      for(int k = 0;k < matr.nrow();k++) {
        csum += matr(k, j);
      }

      if(csum != 0) statesNm1.push_back(statesN1[j]);
    }
    
    intersection = intersect(statesNm1, statesN);
    if(not setequal(intersection, statesNm1)) {
      out = false;
      break;
    }
  }
  return(out);
}

// [[Rcpp::export(.markovchainListRcpp)]]
List markovchainListRcpp(int n, List object, bool include_t0 = false, CharacterVector t0
                           = CharacterVector()) {
  
  bool verify = checkSequenceRcpp(object);
  
  if (not verify) {
    warning("Warning: some states in the markovchain sequences are not contained in the following states!");
  }
  
  // size of result vector
  int sz = n*object.size();
  if(include_t0) sz += n;
  
  int vin = 0; // useful in filling below vectors
  NumericVector iteration(sz);
  CharacterVector values(sz); 
  
  
  S4 ob(object[0]);
  
  CharacterVector sampledValues, newVals;
  
  // Initial State selection if not passed
  //-----------------------------------------------------------------------------
      CharacterVector ustates = ob.slot("states");
      NumericVector rowProbs;
  
      for(int i = 0;i < ustates.size();i++) {
        rowProbs.push_back(1.0 / ustates.size());
      }
  
      bool rselect = (t0.size() == 0);
  
      if(rselect) {
        t0 = sample(ustates, 1, false, rowProbs);
      }
  //------------------------------------------------------------------------------
  
  // check whether t0 is in unique states or not 
  for(int i = 0;i < ustates.size();i++) {
    if(ustates[i] == t0[0]) break;
    
    else if(i == ustates.size()-1) stop("Error! Initial state not defined");
  }
  
  // generate n sequence
  for(int i = 0;i < n;i++) {
    
    // random selection of initial state if not passed to the function
    if(rselect) {
      t0 = sample(ustates, 1, false, rowProbs);
    }
    
    sampledValues = markovchainSequenceRcpp(1, object[0], t0, include_t0);
    
    if(object.size() > 1) {
      for(int j = 1;j < object.size();j++) {
        newVals = markovchainSequenceRcpp(1, object[j], sampledValues);
        sampledValues.push_back(newVals[0]);
      }
    }
    
    for(int k = 0;k < sampledValues.size();k++) {
      iteration[vin] = i + 1;
      values[vin++] = sampledValues[k];
    }
  }
  
  return(List::create(iteration, values));
}

struct MCList : public Worker
{   
  // 3-D matrix where each slice is a transition matrix
  const arma::cube mat;
  
  // number of transition matrices
  const int num_mat;
  
  // matrix where ith row vector store the list of states 
  // names present in ith transition matrix 
  const vector<vector<string> > names;
  
  // vector ehose ith element stor the dimension of ith
  // transition matrix
  const vector<int> size_emat;
  
  // whether to include first state
  const bool include_t0;
  
  // info about initial state
  const bool init; // whether initial state is passed to the method
  const string init_state; // if yes what's the name
  
  // each element of list is a sequence
  list<vector<string> > output;
  
  // constructor for initialization
  MCList(const arma::cube &pmat, const int &pnum_mat, const vector<vector<string> > &pnames, 
         const vector<int> psize_emat, const bool &pinclude_t0, const bool &pinit, const string &pinit_state) : 
    mat(pmat), num_mat(pnum_mat), names(pnames), size_emat(psize_emat), include_t0(pinclude_t0), init(pinit), 
    init_state(pinit_state) {}
  
  
  MCList(const MCList& mclist, Split) : 
    mat(mclist.mat), num_mat(mclist.num_mat), names(mclist.names), size_emat(mclist.size_emat) , 
    include_t0(mclist.include_t0), init(mclist.init), init_state(mclist.init_state) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    
    // to take care of include_t0
    unsigned int ci = 0;
    if(include_t0) ci = 1;
    
    // to store single sequence generated each time
    vector<string> temp(num_mat+ci);
    
    // initial probability and states indices
    arma::vec in_probs(size_emat[0]);
    arma::vec in_states(size_emat[0]);
    
    // assume equal chances of selection of states for the first time
    for(unsigned int i=0;i<in_probs.size();i++) {
      in_probs[i] = 1.0/size_emat[0];
      in_states[i] = i;
    }
    
    // to store the index of the state selected
    arma::vec istate;
    string t0;
    
    // every time generate one sequence
    for(unsigned int p = begin; p < end; p++) {
      
      if(not init) {
        
        // randomly selected state
        istate = rsample(in_states, 1, false, in_probs);
        t0 = names[0][istate[0]];  
      }
      
      else {
        t0 = init_state;
      }
      
      // include the state in the sequence
      if(include_t0) temp[0] = t0;
      
      // to generate one sequence
      for(unsigned int i = 0; i < num_mat; i++) {
        
        // useful for generating rows probabilty vector
        unsigned int j = 0;
        for(j = 0; j < size_emat[i]; j++) {
          if(names[i][j] == t0) break;
        }
        
        // vector to be passed to rsample method
        arma::vec probs(size_emat[i]);
        arma::vec states(size_emat[i]);
        
        for(unsigned int k=0; k < probs.size(); k++) {
          probs[k] = mat(j, k, i);
          states[k] = k;
        }
        
        
        // new state selected
        arma::vec elmt = rsample(states, 1, false, probs);
        t0 = names[i][elmt[0]];
        
        // populate sequence
        temp[i+ci] = t0;
        
      }
      
      // insert one sequence to the output
      output.push_back(temp);  
    }
    
  }
  
  void join(const MCList& rhs) { 
    
    // constant iterator to the first element of rhs.output  
    list<vector<string> >::const_iterator it = rhs.output.begin();
    
    // merge the result of two parallel computation
    for(;it != rhs.output.end();it++) {
      output.push_back(*it);
    }
    
  }
  
};


// Function to generate a list of sequence of states in parallel from non-homogeneous Markov chains.
// 
// Provided any markovchainList object, it returns a list of sequence of states coming 
// from the underlying stationary distribution. 
// 
// @param listObject markovchainList object
// @param n Sample size
// @param include_t0 Specify if the initial state shall be used
// 
// @return A List
// @author Giorgio Spedicato, Deepak Yadav
//   
// @examples
// statesNames <- c("a")
// mcA <- new("markovchain", states = statesNames, transitionMatrix = 
//        matrix(c(1), nrow = 1, byrow = TRUE, 
//        dimnames = list(statesNames, statesNames)))
//  
// statesNames <- c("a","b")
// mcB <- new("markovchain", states = statesNames, transitionMatrix = 
//        matrix(c(0.5, 0.5, 0.3, 0.7), nrow = 2, byrow = TRUE, 
//        dimnames = list(statesNames, statesNames)))
//        
// statesNames <- c("a","b","c")       
// mcC <- new("markovchain", states = statesNames, transitionMatrix = 
//        matrix(c(0.2, 0.5, 0.3, 0, 0.2, 0.8, 0.1, 0.8, 0.1), nrow = 3, 
//        byrow = TRUE, dimnames = list(statesNames, statesNames)))  
//
// mclist <- new("markovchainList", markovchains = list(mcA, mcB, mcC))   
// 
// markovchainSequenceParallelRcpp(mclist, 99999, TRUE)
// 
// 
// [[Rcpp::export(.markovchainSequenceParallelRcpp)]]
List markovchainSequenceParallelRcpp(S4 listObject, int n, bool include_t0 = false,
                                     CharacterVector init_state = CharacterVector()) {
  
  // list of markovchain object
  List object = listObject.slot("markovchains");
  
  bool verify = checkSequenceRcpp(object);
  if (not verify) {
    warning("Warning: some states in the markovchain sequences are not contained in the following states!");
  }
  
  // store number of transition matrices
  int num_matrix = object.size();
  
  // maximum of dimension of all transition matrices
  int max_dim_mat = 0;
  
  // to store the dimension of each transition matrices
  vector<int> size_emat(num_matrix);
  
  // to store list of states in a transition matrix
  CharacterVector states;
  
  // calculate max_dim_mat and populate size_emat
  for(int i = 0; i < num_matrix; i++) {
    
    // extract ith markovchain object
    S4  ob = object[i];
    
    // list of states in ith markovchain object
    states = ob.slot("states");
    
    // keep track of maximun dimension
    if(states.size() > max_dim_mat) max_dim_mat = states.size();
    
    // size of ith transition matrix
    size_emat[i] = states.size();
  }
  
  // Matrix with ith row store the states in ith t-matrix
  vector<vector<string> > names(num_matrix, vector<string>(max_dim_mat));
  
  // to store all t-matrix
  arma::cube mat(max_dim_mat, max_dim_mat, num_matrix);
  mat.fill(0);
  
  for(int i = 0; i < num_matrix;i++) {
    
    // ith markovchain object
    S4  ob = object[i];
    
    // t-matrix and states names
    NumericMatrix tmat = ob.slot("transitionMatrix"); 
    CharacterVector stat_names = ob.slot("states");
    
    // populate 3-D matrix
    for(int j = 0;j < tmat.nrow();j++) {
      for(int k = 0; k < tmat.ncol();k++) {
        
        mat(j, k, i) = tmat(j, k);
        
      }
      
      // populate names of states
      names[i][j] = stat_names[j];
    }
    
  }
  
  // initial state is passed or not
  bool init = false;
  string ini_state;
  
  if(init_state.size() != 0) {
    init = true;
    ini_state = as<string>(init_state[0]);
  }
  
  // create an object of MCList class
  MCList mcList(mat, num_matrix, names, size_emat, include_t0, init, ini_state);
  
  // start parallel computation
  parallelReduce(0, n, mcList);
  
  // list of sequences  
  return wrap(mcList.output);
  
}


// convert a frequency matrix to a transition probability matrix
NumericMatrix _toRowProbs(NumericMatrix x, bool sanitize = false) {
  int nrow = x.nrow(), ncol = x.ncol();
  NumericMatrix out(nrow);

  for (int i = 0; i < nrow; i++) {
    double rowSum = 0;
    for (int j = 0; j < ncol; j++) 
      rowSum += x(i, j);
    for (int j = 0; j < ncol; j++) {
      
      if(sanitize == true) {
        if(rowSum == 0) {
          out(i, j) = 1.0/ncol;  
        } else {
          out(i, j) = x(i, j) / rowSum;
        }
      }
      
      else {
        if(rowSum == 0) {
          out(i, j) = 0;  
        } else {
          out(i, j) = x(i, j) / rowSum;
        }
      }
      
    }
  }
  out.attr("dimnames") = List::create(rownames(x), colnames(x)); 
  return out;
}

// Create a frequency matrix
//' @rdname markovchainFit
//' 
//' @export
// [[Rcpp::export]]
NumericMatrix createSequenceMatrix(SEXP stringChar, bool toRowProbs = false, bool sanitize = false,
                                   CharacterVector possibleStates = CharacterVector()) {
  
  CharacterVector stringchar = as<CharacterVector>(stringChar);
  CharacterVector elements = unique(union_(stringchar, possibleStates)).sort();
  int sizeMatr = elements.size();
  
  NumericMatrix freqMatrix(sizeMatr);
  freqMatrix.attr("dimnames") = List::create(elements, elements); 
  CharacterVector rnames = rownames(freqMatrix);
  
  if(Rf_isMatrix(stringChar)) {
    
    // coerce SEXP to CharacterMatrix
    CharacterMatrix seqMat = as<CharacterMatrix>(stringChar);
    
    // number of columns must be 2
    if(seqMat.ncol() != 2) {
      stop("Number of columns must be 2");
    }
    
    // populate frequency matrix
    int posFrom = 0, posTo = 0;
    for(long long i = 0; i < seqMat.nrow(); i ++) {
      for (int j = 0; j < rnames.size(); j ++) {
        if(seqMat(i, 0) == rnames[j]) posFrom = j;
        if(seqMat(i, 1) == rnames[j]) posTo = j;
      }
      freqMatrix(posFrom, posTo)++;
    }  
  } 
  
  else {
    
    int posFrom = 0, posTo = 0;
    for(long long i = 0; i < stringchar.size() - 1; i ++) {
      for (int j = 0; j < rnames.size(); j ++) {
        if(stringchar[i] == rnames[j]) posFrom = j;
        if(stringchar[i + 1] == rnames[j]) posTo = j;
      }
      freqMatrix(posFrom, posTo)++;
    }  
  }
  
  // sanitizing if any row in the matrix sums to zero by posing the corresponding diagonal equal to 1/dim
  if(sanitize == true)
  {
    for (int i = 0; i < sizeMatr; i++) {
      double rowSum = 0;
      for (int j = 0; j < sizeMatr; j++) 
        rowSum += freqMatrix(i, j);
      if(rowSum == 0)
        for (int j = 0; j < sizeMatr; j++) 
          freqMatrix(i, j) = 1;
    }
  }
  
  if(toRowProbs == true)
    return _toRowProbs(freqMatrix, sanitize);
  
  return (freqMatrix);
  
}

// log-likelihood
double _loglikelihood(CharacterVector seq, NumericMatrix transMatr) {
  
  // to store the result
  double out = 0;
  
  // states names
  CharacterVector rnames = rownames(transMatr);
  
  // caculate out
  int from = 0, to = 0; 
  for(long long i = 0; i < seq.size() - 1; i ++) {
    for(int r = 0; r < rnames.size(); r ++) {
      if(rnames[r] == seq[i]) from = r; 
      if(rnames[r] == seq[i + 1]) to = r; 
    }    
    out += log(transMatr(from, to));
  }
  
  return out;
}

// Fit DTMC using MLE
List _mcFitMle(CharacterVector stringchar, bool byrow, double confidencelevel, bool sanitize = false, 
               CharacterVector possibleStates = CharacterVector()) {
  
  // unique states
  CharacterVector elements = unique(union_(stringchar, possibleStates)).sort();
  // matrix size = nrows = ncols
  int sizeMatr = elements.size();
  
  // initial matrix = transition matrix
  NumericMatrix initialMatr(sizeMatr);
  // frequencey matrix 
  NumericMatrix freqMatr(sizeMatr);
  
  // set names of states as rows name and columns name
  initialMatr.attr("dimnames") = List::create(elements, elements); 

  // populate frequency matrix
  int posFrom = 0, posTo = 0; 
  for(long long i = 0; i < stringchar.size() - 1; i++) {  
    for (int j = 0; j < sizeMatr; j++) {            
      if(stringchar[i] == elements[j]) posFrom = j;
      if(stringchar[i + 1] == elements[j]) posTo = j;
    }
    
    freqMatr(posFrom, posTo)++;
  }

  // take care of rows with all entries 0 
  for (int i = 0; i < sizeMatr; i++) {  
  	double rowSum = 0;
  	for (int j = 0; j < sizeMatr; j++) { 
  	  rowSum += freqMatr(i, j);
  	}
  		
  	// calculate rows probability
  	for (int j = 0; j < sizeMatr; j++) { 
  	  if(rowSum == 0) {
  	    initialMatr(i, j) = (sanitize ? 1.0/sizeMatr : 0);
  	  }
  	  else {
  	    initialMatr(i, j) = freqMatr(i, j)/rowSum;
  	  }
  	}
  }

  // transpose the matrix if byrow is false
  if(byrow == false) {
    initialMatr = _transpose(initialMatr); 
  }

  // z score for given confidence interval
  double zscore = stats::qnorm_0(confidencelevel, 1.0, 0.0);
  
  // store dimension of matrix
  int nrows = initialMatr.nrow();
  int ncols = initialMatr.ncol();
  
  // matrix to store end results
  NumericMatrix lowerEndpointMatr(nrows, ncols);
  NumericMatrix upperEndpointMatr(nrows, ncols);
  NumericMatrix standardError(nrows, ncols);

  // populate above defined matrix 
  double marginOfError, lowerEndpoint, upperEndpoint;
  for(int i = 0; i < nrows; i++) {
    for(int j = 0; j < ncols; j++) {
      if(freqMatr(i, j) == 0) {
        
        // whether entire ith row is zero or not
        bool notrans = true;
        for(int k = 0; k < ncols; k++) {
          
          // if entire ith row is not zero then set notrans to false  
          if(freqMatr(i, k) != 0) {
            standardError(i, j) = lowerEndpointMatr(i, j) = upperEndpointMatr(i, j) = 0;
            notrans = false;
            break;
          }
          
        }
        
        // if entire ith row is zero 
        if(notrans) 
          standardError(i, j) = lowerEndpointMatr(i, j) = upperEndpointMatr(i, j) = 1;
      } 
      else {
        // standard error calculation
        standardError(i, j) = initialMatr(i, j) / sqrt(freqMatr(i, j));
        
        // marginal error calculation
        marginOfError = zscore * standardError(i, j);
        
        // lower and upper end point calculation
        lowerEndpoint = initialMatr(i, j) - marginOfError;
        upperEndpoint = initialMatr(i, j) + marginOfError;
        
        // taking care that upper and lower end point should be between 0(included) and 1(included)
        lowerEndpointMatr(i, j) = (lowerEndpoint > 1.0) ? 1.0 : ((0.0 > lowerEndpoint) ? 0.0 : lowerEndpoint);
        upperEndpointMatr(i, j) = (upperEndpoint > 1.0) ? 1.0 : ((0.0 > upperEndpoint) ? 0.0 : upperEndpoint);
      }
    }
  }
  
  // set the rows and columns name as states names
  standardError.attr("dimnames") = upperEndpointMatr.attr("dimnames") 
          = lowerEndpointMatr.attr("dimnames") = List::create(elements, elements); 

  // create markov chain object
  S4 outMc("markovchain");
  outMc.slot("transitionMatrix") = initialMatr;
  outMc.slot("name") = "MLE Fit";  
  
  // return a list of important results
  return List::create(_["estimate"] = outMc,
                      _["standardError"] = standardError,
		                  _["confidenceInterval"] = List::create(_["confidenceLevel"] = confidencelevel, 
		                                                         _["lowerEndpointMatrix"] = lowerEndpointMatr, 
		                                                         _["upperEndpointMatrix"] = upperEndpointMatr)							
	       );
}

// Fit DTMC using Laplacian smooth
List _mcFitLaplacianSmooth(CharacterVector stringchar, bool byrow, double laplacian = 0.01, bool sanitize = false,
                           CharacterVector possibleStates = CharacterVector()) {
  
  // create frequency matrix
  NumericMatrix origNum = createSequenceMatrix(stringchar, false, sanitize, possibleStates);
  
  // store dimension of frequency matrix
  int nRows = origNum.nrow(), nCols = origNum.ncol();
  
  // convert frequency matrix to transition matrix
  for(int i = 0; i < nRows; i ++) {
	  double rowSum = 0;
	  
	  // add laplacian correction to each entry
	  // also calculate row's sum
	  for(int j = 0; j < nCols; j ++) {
    		origNum(i,j) += laplacian;
    		rowSum += origNum(i,j);
  	}
	  
  	// get a transition matrix and a DTMC
	  for(int j = 0; j < nCols; j ++) { 
	    if(rowSum == 0)
	      origNum(i,j) = sanitize ? origNum(i,j)/rowSum : 0;
	    else 
	      origNum(i,j) = origNum(i,j)/rowSum;
	  }
  }
  
  // transpose transition matrix = columnwise storage 
  if(byrow == false) {
    origNum = _transpose(origNum);  
  }
 
  // create markovchain object
  S4 outMc("markovchain");
  outMc.slot("transitionMatrix") = origNum;
  outMc.slot("name") = "Laplacian Smooth Fit";  

  return List::create(_["estimate"] = outMc);
}

// bootstrap a sequence to produce a list of sample sequences
List _bootstrapCharacterSequences(CharacterVector stringchar, int n, long long size = -1, 
                                  CharacterVector possibleStates = CharacterVector()) {
  
  // store length of sequence
  if(size == -1) {
    size = stringchar.size();  
  }
  
  // frequency matrix
  NumericMatrix contingencyMatrix = createSequenceMatrix(stringchar, true, true, possibleStates);
  
  // many samples from a given a sequence :: bootstrap
  // res list is helper list
  List samples, res;
  
  // state names
  CharacterVector itemset = rownames(contingencyMatrix);
  
  // number of distinct states
  int itemsetsize = itemset.size();

  // access R sample function
  Function sample("sample");
  
  
  for(int i = 0; i < n; i ++) {
    // charseq to store a fresh sequence
    CharacterVector charseq, resvec;	
    
    // select a random number between 0 ans itemsize
    int rnd = (int)(runif(1)(0) * itemsetsize);
    
    // random state
    String ch = itemset[rnd];
    
    // push start state to charseq
    charseq.push_back(ch);
    
    
    for(long long j = 1; j < size; j ++) {
      // store row probability
      NumericVector probsVector;
      
      // populate probsVector
      for(int k = 0; k < itemsetsize; k ++) {
        if((std::string)itemset[k] == (std::string) ch) {
          probsVector = contingencyMatrix(k, _);	
          break;
        }
      }
      
      // select next state of sequence
      res = sample(itemset, 1, true, probsVector);
      resvec = res[0];
      
      // current state
      ch = resvec[0];
      charseq.push_back(ch);
    }
    
    // every add one sequence
    samples.push_back(charseq);
  }

  // return a list of n sequence of same length as of given sequence
  return samples;
}


// estimate from the list of bootstrapped matrices
List _fromBoot2Estimate(List listMatr) {
  
  // number of transition matrices 
  int sampleSize = listMatr.size();
  
  // first transition matrix
  NumericMatrix firstMat = listMatr[0];
  
  // dimension of matrix
  int matrDim = firstMat.nrow();
  
  // matrix to store mean and standard deviation
  NumericMatrix matrMean(matrDim), matrSd(matrDim);

  // populate mean and sd matrix
  for(int i = 0; i < matrDim; i ++) { 
  	for(int j = 0; j < matrDim; j ++) { 
  	  
  	  NumericVector probsEstimated;
  	  for(int k = 0; k < sampleSize; k ++) {
        NumericMatrix mat = listMatr[k];
  	    probsEstimated.push_back(mat(i,j));
  	  }
  	  
  	  matrMean(i, j) = mean(probsEstimated);
  	  matrSd(i, j) = Rcpp::sd(probsEstimated);
  	}
  }
  
  // set rows and columns names = states names
  matrMean.attr("dimnames") = List::create(rownames(firstMat), colnames(firstMat)); 
  matrSd.attr("dimnames") = matrMean.attr("dimnames");
  
  // return list of estimated mean transition matrix and estimated sd matrix 
  return List::create(_["estMu"] = matrMean, _["estSigma"] = matrSd);
}

// Fit DTMC using bootstrap method
List _mcFitBootStrap(CharacterVector data, int nboot, bool byrow, bool parallel, double confidencelevel, bool sanitize = false,
                     CharacterVector possibleStates = CharacterVector()) {
  
  // list of sequence generated using given sequence
  List theList = _bootstrapCharacterSequences(data, nboot, data.size());
  
  // number of new sequence
  int n = theList.size();
  
  // to store frequency matrix for every sequence
  List pmsBootStrapped(n);

  // populate pmsBootStrapped 
  if(parallel)
    for(int i = 0; i < n; i++)
      pmsBootStrapped[i] = createSequenceMatrix(theList[i], true, sanitize, possibleStates);
  
  else 
    for(int i = 0; i < n; i++) 
      pmsBootStrapped[i] = createSequenceMatrix(theList[i], true, sanitize, possibleStates);
  
  
  List estimateList = _fromBoot2Estimate(pmsBootStrapped);
  
  // transition matrix
  NumericMatrix transMatr = _toRowProbs(estimateList["estMu"], sanitize);

  // markovchain object
  S4 estimate("markovchain");
  estimate.slot("transitionMatrix") = transMatr;
  estimate.slot("byrow") = byrow;
  estimate.slot("name") = "BootStrap Estimate";  

  // z score for given confidence interval
  double zscore = stats::qnorm_0(confidencelevel, 1.0, 0.0);
  
  // store dimension of matrix
  int nrows = transMatr.nrow();
  int ncols = transMatr.ncol();
  
  // matrix to store end results
  NumericMatrix lowerEndpointMatr(nrows, ncols), upperEndpointMatr(nrows, ncols);
  NumericMatrix sigma = estimateList["estSigma"], standardError(nrows, ncols);
  
  // populate above defined matrix 
  double marginOfError, lowerEndpoint, upperEndpoint;
  for(int i = 0; i < nrows; i ++) {
    for(int j = 0; j < ncols; j ++) {
      
      // standard error calculation
      standardError(i, j) = sigma(i, j) / sqrt(double(n));
      
      // marginal error calculation
      marginOfError = zscore * standardError(i, j);
      
      // lower and upper end point calculation
      lowerEndpoint = transMatr(i, j) - marginOfError;
      upperEndpoint = transMatr(i, j) + marginOfError;
      
      // taking care that upper and lower end point should be between 0(included) and 1(included)
      lowerEndpointMatr(i, j) = (lowerEndpoint > 1.0) ? 1.0 : ((0.0 > lowerEndpoint) ? 0.0 : lowerEndpoint);
      upperEndpointMatr(i, j) = (upperEndpoint > 1.0) ? 1.0 : ((0.0 > upperEndpoint) ? 0.0 : upperEndpoint);
    }
  }
  
  // set the rows and columns name as states names
  standardError.attr("dimnames") = upperEndpointMatr.attr("dimnames") 
              = lowerEndpointMatr.attr("dimnames") = transMatr.attr("dimnames"); 
  
  // return a list of important results
  List out = List::create(_["estimate"] = estimate,
                          _["standardError"] = standardError,
                          _["confidenceInterval"] = List::create(_["confidenceLevel"] = confidencelevel, 
                                                                 _["lowerEndpointMatrix"] = lowerEndpointMatr,
                                                                 _["upperEndpointMatrix"] = upperEndpointMatr),
                          _["bootStrapSamples"] = pmsBootStrapped
		                     ); 

  return out;
}

// convert matrix data to transition probability matrix
S4 _matr2Mc(CharacterMatrix matrData, double laplacian = 0, bool sanitize = false) {
  
  // dimension of input matrix
  long long nRows = matrData.nrow(), nCols = matrData.ncol();
  
  // set of states
  std::set<std::string> uniqueVals;
  
  // populate uniqueVals set
  for(long long i = 0; i < nRows; i++) 
    for(long long j = 0; j < nCols; j++) 
      uniqueVals.insert((std::string)matrData(i, j));	
  
  // unique states
  int usize = uniqueVals.size();
  
  // matrix of dimension usize
  NumericMatrix contingencyMatrix (usize);
  
  // state names as rows name and columns name
  contingencyMatrix.attr("dimnames") = List::create(uniqueVals, uniqueVals); 
  
  // iterator for set of states
  std::set<std::string>::iterator it;
  
  // populate contingency matrix
  int stateBegin = 0, stateEnd = 0;
  for(long long i = 0; i < nRows; i ++) {
    for(long long j = 1; j < nCols; j ++) {
      
      // row and column number of begin state and end state
      int k = 0;
      for(it = uniqueVals.begin(); it != uniqueVals.end(); ++it, k++) {
        if(*it == (std::string)matrData(i, j-1)) {
          stateBegin = k;
        }
        
        if(*it == (std::string)matrData(i,j)) {
          stateEnd = k;
        }
      }
      
      contingencyMatrix(stateBegin,stateEnd)++;
    }
  }
  
  // add laplacian correction if needed
  for(int i = 0; i < usize; i++) {
    double rowSum = 0;
    for(int j = 0; j < usize; j++) {
      contingencyMatrix(i,j) += laplacian;
      rowSum += contingencyMatrix(i, j);
    }
    
    // get the transition matrix and a DTMC
    for(int j = 0; j < usize; j ++) {
      
      if(sanitize == true) {
        if(rowSum == 0)  {
          contingencyMatrix(i,j) = 1.0/usize;
        } else {
          contingencyMatrix(i,j) /= rowSum;
        }
      } 
      
      else {
        if(rowSum == 0) {
          contingencyMatrix(i,j) = 0;
        } else {
          contingencyMatrix(i,j) /= rowSum;
        }
      }
      
    }
    
  }
  
  // markovchain object
  S4 outMc("markovchain");
  outMc.slot("transitionMatrix") = contingencyMatrix;
  
  return(outMc);
}

//' @name inferHyperparam
//' @title Function to infer the hyperparameters for Bayesian inference from an a priori matrix or a data set
//' @description Since the Bayesian inference approach implemented in the package is based on conjugate priors, 
//'              hyperparameters must be provided to model the prior probability distribution of the chain 
//'              parameters. The hyperparameters are inferred from a given a priori matrix under the assumption 
//'              that the matrix provided corresponds to the mean (expected) values of the chain parameters. A 
//'              scaling factor vector must be provided too. Alternatively, the hyperparameters can be inferred 
//'              from a data set. 
//'              
//' @param transMatr A valid transition matrix, with dimension names.
//' @param scale A vector of scaling factors, each element corresponds to the row names of the provided transition 
//'              matrix transMatr, in the same order. 
//' @param data A data set from which the hyperparameters are inferred.  
//' 
//' @details transMatr and scale need not be provided if data is provided.
//' @return Returns the hyperparameter matrix in a list.
//' 
//' @note The hyperparameter matrix returned is such that the row and column names are sorted alphanumerically, 
//'       and the elements in the matrix are correspondingly permuted. 
//' 
//' @references Yalamanchi SB, Spedicato GA (2015). Bayesian Inference of First Order Markov Chains. R
//'             package version 0.2.5       
//'             
//' @author Sai Bhargav Yalamanchi, Giorgio Spedicato
//' @seealso \code{\link{markovchainFit}}, \code{\link{predictiveDistribution}}
//' 
//' @examples
//' data(rain, package = "markovchain")
//' inferHyperparam(data = rain$rain)
//'  
//' weatherStates <- c("sunny", "cloudy", "rain")
//' weatherMatrix <- matrix(data = c(0.7, 0.2, 0.1, 
//'                                  0.3, 0.4, 0.3, 
//'                                  0.2, 0.4, 0.4), 
//'                         byrow = TRUE, nrow = 3, 
//'                         dimnames = list(weatherStates, weatherStates))
//' inferHyperparam(transMatr = weatherMatrix, scale = c(10, 10, 10))
//'  
//' @export
//'  
// [[Rcpp::export]]
List inferHyperparam(NumericMatrix transMatr = NumericMatrix(), NumericVector scale = NumericVector(), CharacterVector data = CharacterVector()) {
  
  // stop if there is only one element in the matrix and size of data sequence is zero
  if(transMatr.nrow() * transMatr.ncol() == 1 && data.size() == 0)
    stop("Provide the prior transition matrix or the data set in order to infer the hyperparameters");
  
  // to store final result
  List out;
  
  // Number of elements are greater than 1
  if(transMatr.nrow() * transMatr.ncol() != 1) {
    if(scale.size() == 0) {
      stop("Provide a non-zero scaling factor vector to infer integer hyperparameters");
    }
      
    // --------begin validity checks for the transition matrix---------
    if(transMatr.nrow() != transMatr.ncol()) {
      stop("Transition matrix dimensions are inconsistent");
    }
      
    //  number of rows in transition matrix
    int sizeMatr = transMatr.nrow();
    
    // if any element is greater than 1 or less than 0 then raise error
    // sum of each rows must lie between 1 - eps and 1 + eps
    for(int i = 0; i < sizeMatr; i++) {
      double rowSum = 0., eps = 1e-10;
      
      for(int j = 0; j < sizeMatr; j++) {
        if(transMatr(i, j) < 0. || transMatr(i, j) > 1.)
          stop("The entries in the transition matrix must each belong to the interval [0, 1]");
        else
          rowSum += transMatr(i, j);
      }   
        
      if(rowSum <= 1. - eps || rowSum >= 1. + eps)
        stop("Each rows of the transition matrix must sum to 1");
    }
    
    // rows and columns name of transition matrix 
    List dimNames = transMatr.attr("dimnames");
    CharacterVector colNames = dimNames[1];
    CharacterVector rowNames = dimNames[0];
    
    // sorted rows and columns names
    CharacterVector sortedColNames(sizeMatr), sortedRowNames(sizeMatr);
    for(int i = 0; i < sizeMatr; i++) {
      sortedColNames(i) = colNames(i), sortedRowNames(i) = rowNames(i);
    }
    sortedColNames.sort();
    sortedRowNames.sort();
    
    // rows names vector and columns name vector must be same
    // and no names in names vectors should be same
    for(int i = 0; i < sizeMatr; i++) {
      if(i > 0 && (sortedColNames(i) == sortedColNames(i-1) || sortedRowNames(i) == sortedRowNames(i-1)))  
        stop("The states must all be unique");
      else if(sortedColNames(i) != sortedRowNames(i))
        stop("The set of row names must be the same as the set of column names");  
    } 
    // --------end of validity checks for the transition matrix---------  
      
    
    // --------beginning of validity checks for the scale factor vector---------  
    
    // length of scale vector must be equal to number of rows in transition matrix
    if(scale.size() != sizeMatr) 
      stop("The dimensions of the scale vector must match the number of states in the chain");
      
    // if any value in the scale vector is zero  
    for(int i = 0; i < sizeMatr; i++) {
      if(scale(i) == 0)
        stop("The scaling factors must be non-zero!");
    }
    // --------end of validity checks for the scale factor vector---------  
      
    
    // Creation of output matrix i.e. hyper param matrix
    NumericMatrix hpScaled(sizeMatr);
    hpScaled.attr("dimnames") = List::create(rowNames, colNames);
    
    // populate hyper param matrix
    for(int i = 0; i < sizeMatr; i++)
      for(int j = 0; j < sizeMatr; j++)
        hpScaled(i, j) = scale(i) * transMatr(i, j);
    
    /* shift rows and columns so that names of rows
       and columns names will be in sorted order */
    hpScaled = sortByDimNames(hpScaled);
    
    // store list of hyper param scaled matrix
    out = List::create(_["scaledInference"] = hpScaled);
  }
  
  else if(data.size() != 0) {
    
    // to store unique states in sorted order
    CharacterVector elements = data;
    for(int i = 0; i < data.size(); i++)
      elements.push_back(data[i]);
    elements = unique(elements).sort();
    
    // size of hyperparam matrix
    int sizeMatr = elements.size();
    
    // create hyperparam matrix
    NumericMatrix hpData(sizeMatr);
    hpData.attr("dimnames") = List::create(elements, elements); 
    std::fill(hpData.begin(), hpData.end(), 1);
    
    // populate hyper param matrix
    int posFrom = 0, posTo = 0;
    for(long long i = 0; i < data.size() - 1; i ++) {
      for (int j = 0; j < sizeMatr; j ++) {
        if(data[i] == elements[j]) posFrom = j;
        if(data[i + 1] == elements[j]) posTo = j;
      }
      hpData(posFrom,posTo)++;
    }
    
    // ouput data
    out = List::create(_["dataInference"] = hpData);
  }
  
  return out;
}


//' @name markovchainFit
//' @title Function to fit a discrete Markov chain
//' @description Given a sequence of states arising from a stationary state, 
//'  it fits the underlying Markov chain distribution using either MLE (also using a 
//'  Laplacian smoother), bootstrap or by MAP (Bayesian) inference.
//'  
//' @param data A character list.
//' @param method Method used to estimate the Markov chain. Either "mle", "map", "bootstrap" or "laplace"
//' @param byrow it tells whether the output Markov chain should show the transition probabilities by row.
//' @param nboot Number of bootstrap replicates in case "bootstrap" is used.
//' @param laplacian Laplacian smoothing parameter, default zero. It is only used when "laplace" method 
//'                  is chosen.  
//' @param name Optional character for name slot. 
//' @param parallel Use parallel processing when performing Boostrap estimates.
//' @param confidencelevel \deqn{\alpha} level for conficence intervals width. 
//'                        Used only when \code{method} equal to "mle".
//' @param hyperparam Hyperparameter matrix for the a priori distribution. If none is provided, 
//'                   default value of 1 is assigned to each parameter. This must be of size kxk 
//'                   where k is the number of states in the chain and the values should typically 
//'                   be non-negative integers.                        
//' @param stringChar Equivalent to data. Either a nx2 matrix or a character vector.
//' @param toRowProbs converts a sequence matrix into a probability matrix
//' @param sanitize put 1 in all rows having rowSum equal to zero
//' @param possibleStates Possible states which are not present in the given sequence
//' 
//' @return A list containing an estimate, log-likelihood, and, when "bootstrap" method is used, a matrix 
//'         of standards deviations and the bootstrap samples. When the "mle", "bootstrap" or "map" method 
//'         is used, the lower and upper confidence bounds are returned along with the standard error. 
//'         The "map" method also returns the expected value of the parameters with respect to the 
//'         posterior distribution.
//' @references A First Course in Probability (8th Edition), Sheldon Ross, Prentice Hall 2010
//'             
//'             Inferring Markov Chains: Bayesian Estimation, Model Comparison, Entropy Rate, 
//'             and Out-of-Class Modeling, Christopher C. Strelioff, James P. Crutchfield, 
//'             Alfred Hubler, Santa Fe Institute
//' 
//'             Yalamanchi SB, Spedicato GA (2015). Bayesian Inference of First Order Markov Chains. R
//'             package version 0.2.5          
//'             
//' @author Giorgio Spedicato, Tae Seung Kang, Sai Bhargav Yalamanchi
//' @note This function has been rewritten in Rcpp. Bootstrap algorithm has been defined "euristically". 
//'       In addition, parallel facility is not complete, involving only a part of the bootstrap process.
//'       When \code{data} is either a \code{data.frame} or a \code{matrix} object, only MLE fit is 
//'       currently available.
//'       
//' @seealso \code{\link{markovchainSequence}}, \code{\link{markovchainListFit}}
//' @examples
//' sequence <- c("a", "b", "a", "a", "a", "a", "b", "a", "b", "a", "b", "a", "a", 
//'               "b", "b", "b", "a")        
//' sequenceMatr <- createSequenceMatrix(sequence, sanitize = FALSE)
//' mcFitMLE <- markovchainFit(data = sequence)
//' mcFitBSP <- markovchainFit(data = sequence, method = "bootstrap", nboot = 5, name = "Bootstrap Mc")
//'
//' @rdname markovchainFit
//' 
//' @export
//' 
// [[Rcpp::export]]
List markovchainFit(SEXP data, String method = "mle", bool byrow = true, int nboot = 10, double laplacian = 0
            , String name = "", bool parallel = false, double confidencelevel = 0.95, NumericMatrix hyperparam
             = NumericMatrix(), bool sanitize = false, CharacterVector possibleStates = CharacterVector()) {
  
  // list to store the output
  List out;
  
  // if input data is data frame or matrix
  if(Rf_inherits(data, "data.frame") || Rf_isMatrix(data)) { 
    
    // store input data in mat
    CharacterMatrix mat;
    
    // if data is a data.frame force it to matrix
  	if(Rf_inherits(data, "data.frame")) {
  	  DataFrame df(data);
  	  
  	  // matrix : no of rows = no of rows in df : same for number of columns
  	  mat = CharacterMatrix(df.nrows(), df.size());
  	  for(long long i = 0; i < df.size(); i++) {
  	    mat(_,i) = CharacterVector(df[i]);
  	  }
 	  } 
  	else {
		  mat = data;
	  }
  	
  	// byrow assumes distinct observations (trajectiories) are per row
  	// otherwise transpose
  	if(!byrow) {
  	  mat = _transpose(mat); 
  	}
  	
   	S4 outMc =_matr2Mc(mat, laplacian, sanitize);
   	out = List::create(_["estimate"] = outMc);
  } 
  else {
    
    if(method == "mle") {
      out = _mcFitMle(data, byrow, confidencelevel, sanitize, possibleStates); 
    }
    
    if(method == "bootstrap") {
      out = _mcFitBootStrap(data, nboot, byrow, parallel, confidencelevel, sanitize, possibleStates);
    }
  
    if(method == "laplace") {
      out = _mcFitLaplacianSmooth(data, byrow, laplacian, sanitize, possibleStates);
    }
    
    if(method == "map") {
      out = _mcFitMap(data, byrow, confidencelevel, hyperparam, sanitize, possibleStates);
    }
  }
  
  // markovchain object
  S4 estimate = out["estimate"];
  if(name != "") {
    estimate.slot("name") = name; 
  }
  
  // transition matrix
  NumericMatrix transMatr = estimate.slot("transitionMatrix");
  
  // data is neither data frame nor matrix
  if(!Rf_inherits(data, "data.frame") && !Rf_isMatrix(data)) 
    out["logLikelihood"] = _loglikelihood(data, transMatr);
  
  estimate.slot("states") = rownames(transMatr);
  out["estimate"] = estimate;
  
  return out;
}