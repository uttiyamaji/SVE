#include <RcppArmadillo.h>
#include <math.h> 

using namespace arma;
using namespace Rcpp;


vec weight(int m, int n){
  vec weights(n);
  weights.zeros();

  //weights.head(2*m) = 1;

  for(int i=0; i < 2*m; i++){
    weights[i] = 1;
  }

  return weights;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

mat HAC(mat mcond, int m){
  int n = mcond.n_rows;
  int p = mcond.n_cols;
  
  vec wt = weight(m, n);
  
  vec wt_m(2*n); 
  wt_m.zeros();
  
  wt_m.head(n) = wt;
  
  for(int s=n+1; s< 2*n; s++){
    wt_m[s] = wt[2*n-s];
  }
  // eigen values 
  vec eivals = real(fft(wt_m));  // no need for complex here
  
  mat mcond_new(2*n,p);
  mcond_new.zeros();
  mcond_new.head_rows(n) =  mcond; //join_vert(chain, m);
  
  cx_mat mcond_cx(2*n,p);
  mcond_cx.zeros();
  mcond_cx = fft(mcond_new);
  
  for(int i=0;i<p;i++){
    mcond_cx.col(i) = mcond_cx.col(i) % eivals;
  }
  
  mcond_cx = ifft(mcond_cx);
  
  mcond_cx = trans(mcond) * mcond_cx.rows(0,n-1);
  return(real(mcond_cx/n));
  
  
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List initseq_new(mat M){
  int  m, trun;
  mat mu=mean(M);
  //center the rows
  M.each_row() -= mu;
  int n=M.n_rows, p=M.n_cols;
  //Dtm is the vector of det(Sig)'s
  NumericVector Dtm;
  //gam_0 and gam_1 are for gam_2m and gam_2m+1, resp.
  mat Sig(p,p), Sig1(p,p);

  Sig.zeros();
  Sig1.zeros();

  double dtm;

  int sn = n/2;

  for(m = 0; m < n/2; m++){
    Sig = HAC(M ,m );

    if (eig_sym(Sig)(0)>0)
    {
      sn=m;
      break;
    }
  }
  if (sn>n/2-1) 
  {
    stop("Not enough samples.");
  }


  Dtm=det(Sig);
  
  for(m = sn+1; m < n/2; ++m){
    Sig1 = HAC(M, m);

    dtm = det(Sig1); 

    if (dtm<=Dtm(m-sn-1)) break;
    //update Sig
    Sig=Sig1;
    //record dtm
    Dtm.push_back(dtm);
  }

  trun = Dtm.size()-1+sn;
  List res; res["Sig"]=Sig; res["Dtm"]=Dtm; res["trunc"]= trun; res["sn"]=sn;

  return res;

}