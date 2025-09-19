#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

#define ZERO 2e-16 // practical zero

// ______________________________________________________
// L1 NORM A.K.A LASSO

// [[Rcpp::export]]
arma::vec elt_norm_L1(arma::vec x) {
  return(arma::abs(x));
}

// [[Rcpp::export]]
double pen_norm_L1(arma::vec x, arma::vec w) {
  return(arma::accu(arma::abs(w % x)));
}

// [[Rcpp::export]]
double dual_norm_L1(arma::vec x) {
  return(arma::max(arma::abs(x))) ;
}

// [[Rcpp::export]]
arma::vec proximal_L1(const arma::vec& x, const arma::vec& w, double lambda) {
  return(arma::max(arma::zeros(x.n_elem), arma::abs(x) - lambda*w ) % sign(x));
}

// ______________________________________________________
// LINF NORM A.K.A BOUNDED REGRESSION

// [[Rcpp::export]]
arma::vec elt_norm_LINF(const arma::vec& x) {
  return(arma::abs(x));
}

// [[Rcpp::export]]
double pen_norm_LINF(const arma::vec& x, const arma::vec& w) {
  return(arma::max(arma::abs(w % x)));
}

// [[Rcpp::export]]
double dual_norm_LINF(const arma::vec& x) {
  return(arma::accu(arma::abs(x))) ;
}

// [[Rcpp::export]]
arma::vec proximal_LINF(const arma::vec& x, const arma::vec& w, double lambda) {
  arma::uword p = x.n_elem;
  arma::vec u, proj;
  arma::vec res = arma::zeros<arma::vec>(p);
  
  if ( arma::accu(arma::abs(x) / (lambda * w)) >= 1) {
    
    // Reordering absolute values
    u = arma::sort(arma::abs(x), "descend");
    
    // values of the projected coordinate if non zero (dual problem)
    proj = (arma::cumsum(u) - lambda*w)/arma::linspace<arma::vec>(1,p,p);
    
    // selecting non-null entries (dual)
    arma::uvec maxs = arma::sort(arma::find(u-proj > ZERO), "descend") ;
    double thresh = proj[maxs[0]];
    
    // solving primal problem
    // We keep the smallest values and threshold the common values to +- thresh
    for (arma::uword k=0; k < p; k++) {
      if (fabs(x(k)) > ZERO) {
        if (x(k) > 0) {
          res(k) = fmin(fabs(x(k)),thresh);
        } else {
          res(k) = -fmin(fabs(x(k)),thresh);
        }
      }
    }
  }
  return(res);
}

// ______________________________________________________
// L1/L2 NORM A.K.A GROUP-LASSO

// [[Rcpp::export]]
arma::vec elt_norm_L1L2(const arma::vec& x, const arma::uvec& pk) {
  
  arma::vec res = zeros<arma::vec> (pk.n_elem) ; // output with group norms
  arma::uword ind = 0 ; // index to go through the groups
  
  for (arma::uword k=0; k<pk.n_elem; k++) {
    res(k) = norm(x.subvec(ind, ind + pk(k) - 1), 2);
    ind += pk(k);
  }
  
  return(res);
}

// [[Rcpp::export]]
arma::vec elt_norm_L1L2_rep(const arma::vec& x, const arma::uvec& pk) {
  
  arma::vec res = zeros<arma::vec> (x.n_elem) ; // output with group norms
  arma::uword ind = 0 ; // index to go through the groups
  
  for (arma::uword k=0; k<pk.n_elem; k++) {
    double norm_current = norm(x.subvec(ind, ind + pk(k) - 1), 2);
    res.subvec(ind, ind + pk(k) - 1).fill(norm_current);
    ind += pk(k);
  }
  
  return(res);
}

// [[Rcpp::export]]
double pen_norm_L1L2(const arma::vec& x, const arma::uvec& pk, const arma::vec& w) {
  return(accu(w % elt_norm_L1L2(x, pk)));
}

// [[Rcpp::export]]
double dual_norm_L1L2(const arma::vec& x, const arma::uvec& pk) {
  return(max(elt_norm_L1L2(x, pk))) ;
}

// [[Rcpp::export]]
arma::vec proximal_L1L2(const arma::vec& x, const arma::uvec& pk, const arma::vec& w, double lambda) {
  
  arma::vec res = zeros<arma::vec>(x.n_elem);
  arma::uword ind = 0 ;
  
  arma::vec tmp = max(zeros(pk.n_elem), 1-(lambda*w)/elt_norm_L1L2(x, pk)) ;
  
  for (arma::uword k=0; k<pk.n_elem; k++) {
    res.subvec(ind, ind + pk(k) - 1) = tmp(k) * x.subvec(ind, ind + pk(k) - 1);
    ind += pk(k);
  }
  
  return(res);
}

// ______________________________________________________
// L1/LINF NORM A.K.A GROUP-LASSO type 2

// [[Rcpp::export]]
arma::vec elt_norm_L1LINF(const arma::vec& x, const arma::uvec& pk) {
  
  arma::vec res = zeros<arma::vec> (pk.n_elem) ; // output with group norms
  arma::uword ind = 0 ; // index to go through the groups
  
  for (arma::uword k=0; k<pk.n_elem; k++) {
    res(k) = max(abs(x.subvec(ind, ind + pk(k) - 1))) ;
    ind += pk(k);
  }
  
  return(res);
}

// [[Rcpp::export]]
arma::vec elt_norm_L1LINF_rep(const arma::vec& x, const arma::uvec& pk) {
  
  arma::vec res = zeros<arma::vec> (x.n_elem) ; // output with group norms
  arma::uword ind = 0 ; // index to go through the groups

  for (arma::uword k=0; k<pk.n_elem; k++) {
    double norm_current = max(abs(x.subvec(ind, ind + pk(k) - 1)));
    res.subvec(ind, ind + pk(k) - 1).fill(norm_current);
    ind += pk(k);
  }

  return(res);
}

// [[Rcpp::export]]
double pen_norm_L1LINF(const arma::vec& x, const arma::uvec& pk, const arma::vec& w) {
  return(accu(w % elt_norm_L1LINF(x, pk)));
}

// [[Rcpp::export]]
double dual_norm_L1LINF(const arma::vec& x, const arma::uvec& pk) {
  
  arma::vec  res = zeros<arma::vec> (pk.n_elem) ; // output with group norms
  arma::uword ind = 0 ; // index to go through the groups
  
  for (arma::uword k=0; k<pk.n_elem; k++) {
    res(k) = accu(abs(x.subvec(ind, ind + pk(k) - 1)));
    ind += pk(k);
  }
  
  return(max(res)) ;
}

// [[Rcpp::export]]
arma::vec proximal_L1LINF(const arma::vec& x, const arma::uvec& pk, const arma::vec& w, double lambda) {
  
  uword ind = 0, p ;
  arma::vec u, v, proj;
  arma::vec res = zeros<arma::vec>(sum(pk));
  
  for (arma::uword k=0; k<pk.n_elem; k++) {
    v = x.subvec(ind,ind+pk(k)-1) ;
    p = v.n_elem ;
    
    // proximal l-inf
    if ( accu(abs(v) / (lambda*w)) >= 1) {
      // Reordering absolute values
      u = sort(abs(v), "descend");
      
      // values of the projected coordinate if non zero (dual problem)
      proj = (cumsum(u) - lambda*w)/linspace<arma::vec>(1,p,p);
      
      // selecting non-null entries (dual)
      arma::uvec maxs = sort(find(u-proj>ZERO), "descend") ;
      double thresh = proj[maxs[0]];
      
      // solving primal problem
      // We keep the smallest values and threshold the common values to +- thresh
      for (arma::uword j=0; j<p ;j++) {
        if (fabs(v(j)) > ZERO) {
          v(j) = fmin(fabs(v(j)),thresh);
        } else {
          v(j) = -fmin(fabs(v(j)),thresh);
        }
      }
    }
    
    res.subvec(ind,ind+pk(k)-1) =  v ;
    ind += pk(k);
  }
  return(res);
}

// ______________________________________________________
// COOP(ERATIVE) NORM A.K.A COOPERATIVE-LASSO

// [[Rcpp::export]]
arma::vec elt_norm_COOP(const arma::vec& x, const arma::uvec& pk) {
  
  arma::vec  res = zeros<arma::vec> (pk.n_elem) ; // output with group norms
  arma::uword ind = 0 ; // index to go through the groups
  
  for (arma::uword k=0; k<pk.n_elem; k++) {
    res(k) = 
      norm(max(zeros(pk(k)),x.subvec(ind, ind + pk(k) - 1)), 2) + 
      norm(min(zeros(pk(k)),x.subvec(ind, ind + pk(k) - 1)), 2);
    ind += pk(k);
  }
  
  return(res);
}

// [[Rcpp::export]]
arma::vec elt_norm_COOP_rep(const arma::vec& x, const arma::uvec& pk) {
  
  arma::vec  res = zeros<arma::vec> (x.n_elem) ; // output with group norms
  arma::uword ind = 0 ; // index to go through the groups
  
  for (arma::uword k=0; k<pk.n_elem; k++) {
    double norm_current = 
      norm(max(zeros(pk(k)), x.subvec(ind, ind + pk(k) - 1)), 2) + 
      norm(min(zeros(pk(k)), x.subvec(ind, ind + pk(k) - 1)), 2) ;
    res.subvec(ind, ind + pk(k) - 1).fill(norm_current);
    ind += pk(k);
  }

  return(res);
}

// [[Rcpp::export]]
double pen_norm_COOP(const arma::vec& x, const arma::uvec& pk, const arma::vec& w) {
  return(accu(w % elt_norm_COOP(x, pk)));
}

// [[Rcpp::export]]
double dual_norm_COOP(arma::vec x, arma::uvec pk) {
  return(max(elt_norm_COOP(x, pk)));
}

// [[Rcpp::export]]
arma::vec proximal_COOP(const arma::vec& x, const arma::uvec& pk, const arma::vec& w, double lambda) {
  
  arma::vec res = zeros<arma::vec>(x.n_elem);
  arma::uword ind = 0 ;
  
  arma::vec tmp = max(zeros(pk.n_elem), 1-(lambda*w)/elt_norm_COOP(x, pk)) ;
  
  for (arma::uword k=0; k<pk.n_elem; k++) {
    res.subvec(ind, ind + pk(k) - 1) = tmp(k) * x.subvec(ind, ind + pk(k) - 1);
    ind += pk(k);
  }
  return(res);
}

