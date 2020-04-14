#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat roncalli(mat M, int b){
  int n = M.n_rows;
  int p = M.n_cols;

  mat dummy(p,p), out(p,p);
  dummy.zeros();
  out.zeros();

  for(int s = 1; s < b; s++){
    mat place_h(n,p);
    place_h.zeros();
    place_h.tail_rows(n-s) = M.rows(0,n-s-1);

    dummy =  trans(place_h)*M;
    out += dummy;
  }

  out += trans(M)*M;
  return(out/n);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat usual(mat M, int b){
  int n = M.n_rows;
  int p = M.n_cols;

  mat dummy(p,p), out(p,p);
  dummy.zeros();
  out.zeros();

  for(int s = 1; s < b; s++){
    dummy =  trans(M.rows(0,n-s-1))*M.rows(s,n-1);
    out += dummy;
  }

  out += trans(M)*M;
  return(out/n);

}