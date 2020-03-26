#include <RcppArmadillo.h>
#include <math.h> 

using namespace arma;
using namespace Rcpp;

// returns the lag window value for the corresponding window
double lag(int s, double b, String method)
{
  if(method == "bartlett")
  {
    return(1 - s/b);
  }
  else
  {
    return((1 + cos(PI * s/b))/2 );
  }
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

mat msveC_new_2(const arma::mat& chain, double b, String method = "bartlett"){
  int n = chain.n_rows;
  int p = chain.n_cols;
  
  vec wt(n); // wt = rep(0,n)
  wt.zeros();
  for(int s=0; s<b ; ++s){
    wt[s] = lag(s, b, method);
  }
  
  vec wt_m(2*n); 
  wt_m.zeros();
  
  wt_m.head(n) = wt;
  
  for(int s=n+1; s< 2*n; s++){
    wt_m[s] = wt[2*n-s];
  }
  // eigen values 
  vec eivals = real(fft(wt_m));  // no need for complex here
  
  mat chain_new(2*n,p);
  chain_new.zeros();
  chain_new.head_rows(n) =  chain; //join_vert(chain, m);
  
  cx_mat chain_cx(2*n,p);
  chain_cx.zeros();
  chain_cx = fft(chain_new);
  
  for(int i=0;i<p;i++){
    chain_cx.col(i) = chain_cx.col(i) % eivals;
  }
  
  chain_cx = ifft(chain_cx);
  
  chain_cx = trans(chain) * chain_cx.rows(0,n-1);
  return(real(chain_cx/n));
  
  
}