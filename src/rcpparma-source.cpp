// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]


//' Partition of a graph given its centers
//' 
//' Each center in order takes up all its free neighbors. The process is repeated 
//' until there are no more free regions.
//' @title Graph partition
//' @param neighbors A matrix defining a the links between the graph nodes.
//' @param centers The desired centers to start the algorithm
//' @return a vector indicating the cluster center for each graph node.
//' @author Marcelo Costa
// [[Rcpp::export]]
Rcpp::NumericVector RcppPartition(const arma::mat & neighbors, const Rcpp::NumericVector & centers) {
  arma::mat viz = neighbors;
  Rcpp::NumericVector center = centers;
  int n_vizrows = viz.n_rows;
  int n_vizcols = viz.n_cols;
  int n_center  = center.size();
  int i=0, j=0, k, sum=0, x, a;
  arma::mat viznew = viz;
  
  int n, max;
  max = viz(i,j);
  for(i=0;i<n_vizrows;i++){
    for(j=0;j<n_vizcols;j++)
      if(viz(i,j)>max) max = viz(i,j);
  }
  n = max;
  
  Rcpp::NumericVector cluster(n);
  Rcpp::NumericVector filled(n);
  
  for (i=0;i<n;i++){
    cluster[i] = 0;
    filled[i] = 0;
  }
  
  while(sum!=n){
    sum = 0;
    for(i = 0;i < n_center;i++){
      x = center[i];
      cluster[x-1] = x;
      filled[x-1] = 1;
      for(j=0;j<n_vizrows;j++){
        if((viz(j,0)==x) && (filled[viz(j,1)-1] == 0)) {
          a = viz(j,1);
          cluster[viz(j,1)-1] = x;
          filled[viz(j,1)-1]  = 1;
          for (k=0;k<n_vizrows;k++){
            if (viz(k,0)==a) {
              viznew(k,0) = x;
            }
          }
        }
        
      }
    }
    viz = viznew;
    for(i=0;i<n;i++){
      sum = sum + filled[i];
    }
  }

  return cluster;
}

//' Generate frequency matrix
//' 
//' Frequency matrix
//' @title Frequency matrix generation
//' @param partitions A vector indicating the cluster center for each graph node.
//' @return a matrix.
//' @author Leandro Mineti
// [[Rcpp::export]]
Rcpp::NumericMatrix RcppFreqMatrix(const Rcpp::NumericVector & partitions) {
  Rcpp::NumericVector part = partitions;
  int n = part.size(), mat_index, mat_test;
  arma::mat mat = arma::zeros(n,n);
  
  for(mat_index=0; mat_index<n; mat_index++){
    for(mat_test=mat_index+1; mat_test<n; mat_test++){
      if(part[mat_index] == part[mat_test]) mat(mat_index,mat_test) = 1;
    }
  }
  return Rcpp::wrap(mat);
}
