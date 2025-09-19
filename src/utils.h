#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

#include <cmath>
#include <algorithm>

using namespace Rcpp;
using namespace arma;

#define ZERO 2e-16 // practical zero

// =====================================================
// Prox operators
// =====================================================

// Lasso proximal operator

// [[Rcpp::export]]
void proximal_las(const arma::vec& lambda, double L, arma::vec& u) {
  for (arma::uword j = 0; j < u.n_elem; j++) {
    double shrink = std::max(0.0, 1.0 - lambda[j] / (L * std::abs(u[j])));
    u[j] *= shrink;
  }
}

// Group Lasso

// [[Rcpp::export]]
void proximal_grp(const arma::vec& lambda, double L,
                  arma::vec& u, const arma::ivec& pk) {
  int K = pk.n_elem;
  int ind = 0;
  for (int k = 0; k < K; k++) {
    arma::vec sub = u.subvec(ind, ind + pk[k] - 1);
    double normk = arma::norm(sub, 2);
    double shrink = (normk > 0) ? std::max(0.0, 1.0 - lambda[k] / (L * normk)) : 0.0;
    u.subvec(ind, ind + pk[k] - 1) = sub * shrink;
    ind += pk[k];
  }
}

// Cooperative Lasso

// [[Rcpp::export]]
void proximal_coo(const arma::vec& lambda, double L,
                  arma::vec& u, const arma::ivec& pk) {
  int K = pk.n_elem;
  int ind = 0;
  for (int k = 0; k < K; k++) {
    arma::vec sub = u.subvec(ind, ind + pk[k] - 1);
    
    arma::vec sub_pos = arma::clamp(sub, 0, arma::datum::inf);
    arma::vec sub_neg = arma::clamp(-sub, 0, arma::datum::inf);
    
    double norm_pos = arma::norm(sub_pos, 2);
    double norm_neg = arma::norm(sub_neg, 2);
    
    double shrink_pos = (norm_pos > 0) ? std::max(0.0, 1.0 - lambda[k] / (L * norm_pos)) : 0.0;
    double shrink_neg = (norm_neg > 0) ? std::max(0.0, 1.0 - lambda[k] / (L * norm_neg)) : 0.0;
    
    for (int j = 0; j < pk[k]; j++) {
      if (sub[j] > 0) {
        u[ind + j] = sub[j] * shrink_pos;
      } else {
        u[ind + j] = sub[j] * shrink_neg;
      }
    }
    ind += pk[k];
  }
}

// Tree Group Lasso

// [[Rcpp::export]]
void proximal_tree_grp(const arma::ivec& Ks, const arma::vec& lambda,
                       double L, arma::vec& u, const arma::ivec& pk) {
  int k0 = 0;
  for (int h = 0; h < Ks.n_elem; h++) {
    int Kc = Ks[h];
    arma::ivec pkc = pk.subvec(k0, k0 + Kc - 1);
    arma::vec lambdac = lambda.subvec(k0, k0 + Kc - 1);
    
    proximal_grp(lambdac, L, u, pkc);
    k0 += Kc;
  }
}

// Tree Cooperative Lasso

// [[Rcpp::export]]
void proximal_tree_coo(const arma::ivec& Ks, const arma::vec& lambda,
                       double L, arma::vec& u, const arma::ivec& pk) {
  int k0 = 0;
  for (int h = 0; h < Ks.n_elem; h++) {
    int Kc = Ks[h];
    arma::ivec pkc = pk.subvec(k0, k0 + Kc - 1);
    arma::vec lambdac = lambda.subvec(k0, k0 + Kc - 1);
    
    proximal_coo(lambdac, L, u, pkc);
    k0 += Kc;
  }
}
