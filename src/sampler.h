// sampler.h: Rcpp/Armadillo equivalent to R's sample().  
// This is intended for use in C++ functions, and should *not* be called from R.
// It should yield identical results to R in most cases, 
// and stop with errors when results are expected to differ.
//
// Copyright (C)  2012 - 2014  Christian Gunning
// Copyright (C)  2013  Romain Francois
//
// This file is part of RcppArmadillo.
//
// RcppArmadillo is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// RcppArmadillo is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with RcppArmadillo.  If not, see <http://www.gnu.org/licenses/>.

#ifndef SAMPLE_H
#define SAMPLE_H

arma::vec rsample(const arma::vec &x, const int size, const bool replace, arma::vec prob_ );
void RSampleNoReplace(arma::vec &index, int nOrig, int size);
void RSampleReplace(arma::vec &index, int nOrig, int size);
void RProbSampleNoReplace(arma::vec &index, int nOrig, int size, arma::vec &prob);
void RProbSampleReplace(arma::vec &index, int nOrig, int size, arma::vec &prob);
void RWalkerProbSampleReplace(arma::vec &index, int nOrig, int size, arma::vec &prob);
void RFixProb(arma::vec &prob, const int size, const bool replace);

#endif
