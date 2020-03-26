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

mat msveC_new_1(const arma::cx_mat& chain, double b, String method = "bartlett"){
  int n = chain.n_rows;
  int p = chain.n_cols;
  
  vec wt(n); // wt = rep(0,n)
  wt.zeros();
  for(int s=0; s<b ; ++s){
     wt[s] = lag(s, b, method);
  }

  
  vec wt_m(2*n); 
  wt_m.zeros();
 
  for(int s=0; s<n ;s++){
    wt_m[s] = wt[s];
  }
  
  for(int s=n+1; s< 2*n; s++){
    wt_m[s] = wt[2*n-s];
  }
  // eigen values 
  cx_vec eivals = fft(wt_m);

  cx_mat m(n,p); 
  m.zeros();

  cx_mat chain_new(2*n,p);
  chain_new.zeros();

  
  chain_new = join_vert(chain, m);
  chain_new = fft(chain_new);

  cx_mat eival_mat(2*n,p);
  eival_mat.zeros();

  for(int i=0;i<p;i++){
     eival_mat.col(i) = eivals;
  }

  chain_new = chain_new % eival_mat; 



  chain_new = ifft(chain_new);

  mat result(n,p);
  result.zeros(); 

  result = trans(real(chain)) * real(chain_new.rows(0,n-1));
  return(result/n);

  //return(((real(chain_new.rows(0,n-1).t()))*real(chain_new.rows(0,n-1)))/n);
   
}