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

mat msveC_new(const arma::mat& chain, double b, String method = "bartlett"){
  int n = chain.n_rows;
  int p = chain.n_cols;
  
  vec wt(n); // wt = rep(0,n)
  for(int s=0; s<b ; ++s){
     wt[s] = lag(s, b, method);
  }

  
  vec wt_m(2*n); 
 
  for(int s=0; s<n ;s++){
    wt_m[s] = wt[s];
  }
  
  wt_m[n] = 0;
  for(int s=n+1; s< 2*n; s++){
    wt_m[s] = wt[2*n-1-s];
  }
  // eigen values 
  vec evals = real(fft(wt_m));

  mat m(n,p); 
  m.zeros();

  mat chain_new = join_vert(chain,m); // rbind(mcond, matrix(0,n,q))

  cx_mat chain_fft = fft(chain_new);

  mat eval_mat(2*n,2*n);
  eval_mat.zeros();
  for(int i=0; i< (2*n); i++){
    eval_mat(i,i) = evals[i];
  }

  chain_fft = eval_mat*chain_fft;

  mat chain_ifft = real(ifft(chain_fft));

  mat s_chain = chain_ifft.rows(0,n-1);

  mat result = trans(chain)*s_chain;
  return(result/n);

}


