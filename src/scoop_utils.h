// #include <R.h>
// #include <Rmath.h>
// #include <math.h>
// #include <R_ext/Utils.h>
// 
// void groupnorm(int     *p,
// 	       int    *K,
// 	       int    pk[*K],
// 	       double beta[*p],
// 	       double norm[*K]) {
//   
//   int ind=0;
//   
//   for (int k=0; k<*K; k++) {
//     norm[k] = 0.0;
//     for (int i=0; i<pk[k]; i++) {
//       norm[k] += pow(beta[ind],2);
//       ind++;
//     }
//     norm[k] = sqrt(norm[k]);
//   }
// }
// 
// double coopnorm(int     *p,
// 		int    *K,
// 		int    pk[*K],
// 		double beta[*p],
// 		double norm[*K]) {
// 
//   int ind=0;
//   double norm_pos, norm_neg;
// 
//   for (int k=0; k<*K; k++) {
//     norm_pos = 0;
//     norm_neg = 0;
//     for (int i=0; i<pk[k]; i++) {
//       norm_pos += pow(fmax(0, beta[ind]),2);
//       norm_neg += pow(fmax(0,-beta[ind]),2);
//       ind++;
//     }
//     norm[k] = sqrt(norm_pos) + sqrt(norm_neg);
//   }
// }
// 
// void groupnormrep(int   *p,
// 		  int    *K,
// 		  int    pk[*K],
// 		  double beta[*p],
// 		  double norm[*p]) {
//   
//   int ind=0, ind0;
//   double normk;
//   
//   for (int k=0; k<*K; k++) {
//     ind0=ind;
//     normk=0;
//     for (int i=0; i<pk[k]; i++) {
//       normk += pow(beta[ind],2);
//       ind++;
//     }
//     for (int j=ind0; j<(ind0+pk[k]); j++) {
//       norm[j] = sqrt(normk);
//     }
//   }
// }
// 
// void proximal_las(int     p,
// 		  double lambda[p],
// 		  double L,
// 		  double u[p]) {
//   
//   for (int j=0; j<p; j++) {
//     u[j] = fmax(0,1-lambda[j]/(L*fabs(u[j]))) * u[j];
//   }
// }
// 
// void proximal_grp(int     p,
// 		  int     K,
// 		  int     pk[K],
// 		  double lambda[K],
// 		  double L,
// 		  double u[p]) {
//   
//   int ind=0, ind0;
//   double normk, shrink;
//   
//   for (int k=0; k<K; k++) {
//     ind0=ind;
//     normk=0;
//     for (int i=0; i<pk[k]; i++) {
//       normk += pow(u[ind],2);
//       ind++;
//     }
//     shrink = lambda[k] / (L * sqrt(normk));
//     for (int j=ind0; j<(ind0+pk[k]); j++) {
//       u[j] = fmax(0,1-shrink) * u[j];
//     }
//   }
// }
// 
// void proximal_coo(int     p,
// 		  int     K,
// 		  int     pk[K],
// 		  double lambda[K],
// 		  double L,
// 		  double u[p]) {
//   
//   int ind=0, ind0;
//   double normk, normk_pos, normk_neg, shrink_pos, shrink_neg;
//   
//   for (int k=0; k<K; k++) {
//     ind0=ind;
//     normk_pos=0;
//     normk_neg=0;
//     for (int i=0; i<pk[k]; i++) {
//       if(u[ind] > 0) {
// 	normk_pos += pow(u[ind],2);
//       } else {
// 	normk_neg += pow(u[ind],2);
//       }
//       ind++;
//     }
//     shrink_pos = lambda[k] / (L * sqrt(normk_pos));
//     shrink_neg = lambda[k] / (L * sqrt(normk_neg));
//     for (int j=ind0; j<(ind0+pk[k]); j++) {
//       if(u[j] > 0) {
// 	u[j] = fmax(0,1-shrink_pos) * u[j];
//       } else {
// 	u[j] = fmax(0,1-shrink_neg) * u[j];	
//       }
//     }
//   }
// }
// 
// void proximal_tree_grp(
// 		   int    p        ,
//  		   int    H        ,
// 		   int    K        ,
// 		   int    Ks[H]    ,
// 		   int    pk[K]    ,
// 		   double lambda[K],
// 		   double L,
// 		   double u[p]     ) {
//   
//   int h=0, k=0, k0=0, Kc;
//   
//   // For each level of the tree (starting from the bottom)
//   for (h=0; h<H; h++) {
//     Kc = Ks[h];
//     int pkc[Kc];
//     double lambdac[Kc];
//     for (k=0; k<Kc; k++) {
//       pkc[k] = pk[k+k0];
//       lambdac[k] = lambda[k+k0];
//     }
// 
//     proximal_grp(p, Kc, pkc, lambdac, L, u);
//     k0 = k0+Kc;
//   }
// }
// 
// void proximal_tree_coo(
// 		   int    p        ,
//  		   int    H        ,
// 		   int    K        ,
// 		   int    Ks[H]    ,
// 		   int    pk[K]    ,
// 		   double lambda[K],
// 		   double L,
// 		   double u[p]     ) {
//   
//   int h=0, k=0, k0=0, Kc;
//   
//   // For each level of the tree (starting from the bottom)
//   for (h=0; h<H; h++) {
//     Kc = Ks[h];
//     int pkc[Kc];
//     double lambdac[Kc];
//     for (k=0; k<Kc; k++) {
//       pkc[k] = pk[k+k0];
//       lambdac[k] = lambda[k+k0];
//     }
// 
//     proximal_coo(p, Kc, pkc, lambdac, L, u);
//     k0 = k0+Kc;
//   }
// }
// 
// void proximal_tree_grp_standalone(
// 		   int    *p        ,
// 		   int    *H        ,
// 		   int    *K        ,
// 		   int    Ks[*H]    ,
// 		   int    pk[*K]    ,
// 		   double lambda[*K],
// 		   double *L,
// 		   double u[*p]     ) {
//   
//   int h=0, k=0, k0=0, Kc;
//   
//   for (h=0; h<*H; h++) {
//     Kc = Ks[h];
//     int pkc[Kc];
//     double lambdac[Kc];
//     for (k=0; k<Kc; k++) {
//       pkc[k] = pk[k+k0];
//       lambdac[k] = lambda[k+k0];
//     }
//     proximal_grp(*p, Kc, pkc, lambdac, *L, u);
//     k0 = k0+Kc;
//   }
// }
// 
// void proximal_tree_coo_standalone(
// 		   int    *p        ,
// 		   int    *H        ,
// 		   int    *K        ,
// 		   int    Ks[*H]    ,
// 		   int    pk[*K]    ,
// 		   double lambda[*K],
// 		   double *L,
// 		   double u[*p]     ) {
//   
//   int h=0, k=0, k0=0, Kc;
//   
//   for (h=0; h<*H; h++) {
//     Kc = Ks[h];
//     int pkc[Kc];
//     double lambdac[Kc];
//     for (k=0; k<Kc; k++) {
//       pkc[k] = pk[k+k0];
//       lambdac[k] = lambda[k+k0];
//     }
//     proximal_coo(*p, Kc, pkc, lambdac, *L, u);
//     k0 = k0+Kc;
//   }
// }
