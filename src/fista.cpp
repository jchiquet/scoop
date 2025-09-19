#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

#include <cmath>
#include <algorithm>
#include "utils.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
void ISTA_LM(
  const arma::ivec& Ks,
  const arma::ivec& pk,
  const arma::mat& XtX,
  const arma::vec& Xty,
  const arma::vec& lambda,
  double L,
  int max_it,
  double eps,
  arma::vec& x0,
  arma::vec& xk,
  int penalty, // 1 = Lasso, 2 = Group-Lasso, ...
  int& iter
) {
  double delta = 2 * eps;
  
  while ((delta > eps) && (iter < max_it)) {
    
    xk = x0 - (XtX * x0 - Xty) / L;
    
    // proximal operator selon le penalty
    switch (penalty) {
      case 1: proximal_las(lambda, L, xk); break;
      case 2: proximal_grp(lambda, L, xk, pk); break;
      case 3: proximal_coo(lambda, L, xk, pk); break;
      case 4: proximal_tree_grp(Ks, lambda, L, xk, pk); break;
      case 5: proximal_tree_coo(Ks, lambda, L, xk, pk); break;
    }
    
    // assess convergence
    delta = arma::norm(x0 - xk, 2);
    
    // parameter update
    x0 = xk;
    iter++;
    
    R_CheckUserInterrupt();
  }
}

// [[Rcpp::export]]
void FISTA_LM(
  const arma::ivec& Ks,
  const arma::ivec& pk,
  const arma::mat& XtX,
  const arma::vec& Xty,
  const arma::vec& lambda,
  double L,
  int max_it,
  double eps,
  arma::vec& x0,
  arma::vec& xk,
  int penalty,
  int& iter
) {
  double delta = 2 * eps;
  double t0 = 1.0, tk;
  arma::vec s = x0;
  
  while ((delta > eps) && (iter < max_it)) {
    
    xk = s - (XtX * s - Xty) / L;
    
    // apply proximal operator
    switch (penalty) {
      case 1: proximal_las(lambda, L, xk); break;
      case 2: proximal_grp(lambda, L, xk, pk); break;
      case 3: proximal_coo(lambda, L, xk, pk); break;
      case 4: proximal_tree_grp(Ks, lambda, L, xk, pk); break;
      case 5: proximal_tree_coo(Ks, lambda, L, xk, pk); break;
    }
    
    // FISTA auxiliary variables 
    tk = 0.5 * (1 + std::sqrt(1 + 4 * t0 * t0));
    s = xk + (t0 - 1.0) / tk * (xk - x0);
    
    // Convergence
    delta = arma::norm(x0 - xk, 2);
    
    // parameters update
    t0 = tk;
    x0 = xk;
    iter++;
    
    R_CheckUserInterrupt();
  }
}

// ISTA pour la régression logistique avec recherche de pas

// [[Rcpp::export]]
void ISTA_LRM(
  const arma::ivec& Ks,
  const arma::ivec& pk,
  const arma::mat& X,
  const arma::vec& y,
  const arma::vec& lambda,
  double L0,
  int max_it,
  double eps,
  arma::vec& x0,
  arma::vec& xk,
  int penalty,
  int& iter
) {
  double delta = 2 * eps;
  double normDiff2 = 0.0, scalDiffd = 0.0;
  double fx0 = 0.0, fxk = 0.0;
  double L = 2.0, maxL = L0, ratio = 1.0;
  bool found = false;
  
  arma::vec eta0(y.n_elem), etak(y.n_elem), df(x0.n_elem);
  
  while ((delta > eps) && (iter < max_it)) {
    
    eta0 = X * x0;
    
    // df = - X * (y - sigmoid(eta0))
    arma::vec prob = 1.0 / (1.0 + arma::exp(-eta0));
    df = - X.t() * (y - prob);
  
    // fx0 = -sum(y*eta0 - log(1+exp(eta0)))
    fx0 = -arma::accu(y % eta0 - arma::log(1 + arma::exp(eta0)));
    
    // récupération de L antérieur
    L = std::fmin(ratio, L / 2.0);
    
    // recherche linéaire
    found = false;
    while (!found) {
      // tentative : descente
      xk = x0 - df / L;
      
      // proximal
      switch (penalty) {
        case 1: proximal_las(lambda, L, xk); break;
        case 2: proximal_grp(lambda, L, xk, pk); break;
        case 3: proximal_coo(lambda, L, xk, pk); break;
        case 4: proximal_tree_grp(Ks, lambda, L, xk, pk); break;
        case 5: proximal_tree_coo(Ks, lambda, L, xk, pk); break;
      }
      
      // etak = X^T * xk
      etak = X * xk;
      
      // norme et produit scalaire
      arma::vec diff = xk - x0;
      normDiff2 = arma::dot(diff, diff);
      scalDiffd = arma::dot(diff, df);
      
      // fxk = -sum(y*etak - log(1+exp(etak)))
      fxk = -arma::accu(y % etak - arma::log(1 + arma::exp(etak)));
      
      // ratio = 2 * (fxk - fx0 - scalDiffd) / normDiff2
      ratio = 2.0 * (fxk - fx0 - scalDiffd) / normDiff2;
      
      if (L >= ratio || L >= maxL) {
        found = true;
      } else {
        L = std::fmin(std::fmax(2 * L, ratio), maxL);
      }
      
      R_CheckUserInterrupt();
    }
    
    // préparation itération suivante
    delta = std::sqrt(normDiff2);
    x0 = xk;
    found = false;
    iter++;
    
    R_CheckUserInterrupt();
  }
}

// FISTA pour la régression logistique avec recherche de pas

// [[Rcpp::export]]
void FISTA_LRM(
    const arma::ivec& Ks,
    const arma::ivec& pk,
    const arma::mat& X,
    const arma::vec& y,
    const arma::vec& lambda,
    double L0,
    int max_it,
    double eps,
    arma::vec& x0,
    arma::vec& xk,
    int penalty,
    int& iter
) {
  double delta = 2 * eps;
  double normDiff2 = 0.0, scalDiffd = 0.0;
  double fy0 = 0.0, fyk = 0.0;
  double t0 = 1.0, tk;
  double L = 2.0, maxL = L0, ratio = 1.0;
  bool found = false;
  
  arma::vec s = x0;
  arma::vec eta0(y.n_elem), etak(y.n_elem), df(x0.n_elem);
  
  while ((delta > eps) && (iter < max_it)) {
    
    eta0 = X * s;
    
    // df = -X.t() * (y - sigmoid(eta0))
    arma::vec prob = 1.0 / (1.0 + arma::exp(-eta0));
    df = - X.t() * (y - prob);
    
    // fy0
    fy0 = -arma::accu(y % eta0 - arma::log(1 + arma::exp(eta0)));
    
    // récupération L
    L = std::fmin(ratio, L / 2.0);
    
    // line search
    found = false;
    while (!found) {
      // tentative descente
      xk = s - df / L;
      
      // proximal
      switch (penalty) {
        case 1: proximal_las(lambda, L, xk); break;
        case 2: proximal_grp(lambda, L, xk, pk); break;
        case 3: proximal_coo(lambda, L, xk, pk); break;
        case 4: proximal_tree_grp(Ks, lambda, L, xk, pk); break;
        case 5: proximal_tree_coo(Ks, lambda, L, xk, pk); break;
      }
      
      // etak = X^T * xk
      etak = X * xk;
      
      arma::vec diff = xk - s;
      normDiff2 = arma::dot(diff, diff);
      scalDiffd = arma::dot(diff, df);
      
      fyk = -arma::accu(y % etak - arma::log(1 + arma::exp(etak)));
      
      ratio = 2.0 * (fyk - fy0 - scalDiffd) / normDiff2;
      
      if (L >= ratio || L >= maxL) {
        found = true;
      } else {
        L = std::fmin(std::fmax(2 * L, ratio), maxL);
      }
      
      R_CheckUserInterrupt();
    }
    
    // mise à jour t
    tk = 0.5 * (1 + std::sqrt(1 + 4 * t0 * t0));
    
    // mise à jour s (extrapolation Nesterov)
    s = xk + (t0 - 1.0) / tk * (xk - x0);
    
    // préparation
    delta = std::sqrt(normDiff2);
    t0 = tk;
    x0 = xk;
    found = false;
    iter++;
    
    R_CheckUserInterrupt();
  }
}
