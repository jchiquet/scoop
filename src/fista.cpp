#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

#include <cmath>
#include <algorithm>
#include "utils.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
void ISTA_LM(
  int p,
  int H,
  int K,
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
    // u = x0 - (XtX * x0 - Xty)/L
    xk = x0 - (XtX * x0 - Xty) / L;
    
    // proximal operator selon le penalty
    switch (penalty) {
      case 1: proximal_las(lambda, L, xk); break;
      case 2: proximal_grp(p, K, pk, lambda, L, xk); break;
      case 3: proximal_coo(p, K, pk, lambda, L, xk); break;
      case 4: proximal_tree_grp(p, H, K, Ks, pk, lambda, L, xk); break;
      case 5: proximal_tree_coo(p, H, K, Ks, pk, lambda, L, xk); break;
    }
    
    // calcul de la différence
    delta = arma::norm(x0 - xk, 2);
    
    // mise à jour
    x0 = xk;
    iter++;
    
    R_CheckUserInterrupt(); // si tu compiles toujours avec Rcpp
  }
}

// [[Rcpp::export]]
void FISTA_LM(
  int p,
  int H,
  int K,
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
    // u = s - (XtX*s - Xty)/L
    xk = s - (XtX * s - Xty) / L;
    
    // proximal
    switch (penalty) {
      case 1: proximal_las(lambda, L, xk); break;
      case 2: proximal_grp(p, K, pk, lambda, L, xk); break;
      case 3: proximal_coo(p, K, pk, lambda, L, xk); break;
      case 4: proximal_tree_grp(p, H, K, Ks, pk, lambda, L, xk); break;
      case 5: proximal_tree_coo(p, H, K, Ks, pk, lambda, L, xk); break;
    }
    
    // mise à jour t
    tk = 0.5 * (1 + std::sqrt(1 + 4 * t0 * t0));
    
    // mise à jour yk
    s = xk + (t0 - 1.0) / tk * (xk - x0);
    
    // delta
    delta = arma::norm(x0 - xk, 2);
    
    // update
    t0 = tk;
    x0 = xk;
    iter++;
    
    R_CheckUserInterrupt();
  }
}

// ISTA pour la régression logistique avec recherche de pas

// [[Rcpp::export]]
void ISTA_LRM(
  int p,
  int H,
  int K,
  const arma::ivec& Ks,
  const arma::ivec& pk,
  int n,
  const arma::mat& X,   // taille p x n
  const arma::vec& y,   // taille n
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
  
  arma::vec eta0(n), etak(n), df(p);
  
  while ((delta > eps) && (iter < max_it)) {
    // eta0 = X^T * x0  (NB: X est p x n, donc X.t() est n x p)
    eta0 = X.t() * x0;
    
    // df = - X * (y - sigmoid(eta0))
    arma::vec prob = 1.0 / (1.0 + arma::exp(-eta0));
    df = - X * (y - prob);
    
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
        case 2: proximal_grp(p, K, pk, lambda, L, xk); break;
        case 3: proximal_coo(p, K, pk, lambda, L, xk); break;
        case 4: proximal_tree_grp(p, H, K, Ks, pk, lambda, L, xk); break;
        case 5: proximal_tree_coo(p, H, K, Ks, pk, lambda, L, xk); break;
      }
      
      // etak = X^T * xk
      etak = X.t() * xk;
      
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
    int p,
    int H,
    int K,
    const arma::ivec& Ks,
    const arma::ivec& pk,
    int n,
    const arma::mat& X,   // taille p x n
    const arma::vec& y,   // taille n
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
  arma::vec eta0(n), etak(n), df(p);
  
  while ((delta > eps) && (iter < max_it)) {
    // eta0 = X^T * s
    eta0 = X.t() * s;
    
    // df = -X * (y - sigmoid(eta0))
    arma::vec prob = 1.0 / (1.0 + arma::exp(-eta0));
    df = - X * (y - prob);
    
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
        case 2: proximal_grp(p, K, pk, lambda, L, xk); break;
        case 3: proximal_coo(p, K, pk, lambda, L, xk); break;
        case 4: proximal_tree_grp(p, H, K, Ks, pk, lambda, L, xk); break;
        case 5: proximal_tree_coo(p, H, K, Ks, pk, lambda, L, xk); break;
      }
      
      // etak = X^T * xk
      etak = X.t() * xk;
      
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
