#include "scoop_utils.h"

void ISTA_LM(int *p,
	     int *H,
	     int *K,
	     int Ks[*H],
	     int pk[*K],
	     double XtX[*p][*p],
	     double Xty[*p],
	     double lambda[*K],
	     double *L,
	     int *max_it,
	     double *eps,
	     double x0[*p],
	     double xk[*p],
	     int *penalty, // 1 = Lasso, 2 = Group-Lasso, 3 = Coop-Lasso, 4 = Tree-Group, 5 = Tree-Coop
	     int *i) {
  
  double delta = 2 * *eps;
  int j,l;

  while ( (delta > *eps) && (*i < *max_it) ) {

    // arg1 proximal: u = x0 - (XtX %*% x0 - Xty)/L
    for (j=0; j<*p; j++) {
      xk[j] = 0;
      for (l=0;l<*p; l++) {
	xk[j] += XtX[j][l] * x0[l];
      }
      xk[j] = x0[j] - (xk[j] - Xty[j]) / *L;
    }

    switch(*penalty) {
      case 1: proximal_las(*p, lambda, *L, xk); 
	break;
      case 2: proximal_grp(*p, *K, pk, lambda, *L, xk);
	break;
      case 3: proximal_coo(*p, *K, pk, lambda, *L, xk);
	break;
      case 4: proximal_tree_grp(*p, *H, *K, Ks, pk, lambda, *L, xk);
	break;
      case 5: proximal_tree_coo(*p, *H, *K, Ks, pk, lambda, *L, xk);
    }

    // preparing next iterate
    delta = 0;
    for(j=0;j<*p;j++) {
      delta += pow(x0[j]-xk[j],2);
    }
    delta = sqrt(delta);
    
    for(j=0;j<*p;j++) {
      x0[j] = xk[j];
    }
    *i += 1;

    R_CheckUserInterrupt();
  }
}

void FISTA_LM(int *p,
	      int *H,
	      int *K,
	      int Ks[*H],
	      int pk[*K],
	      double XtX[*p][*p],
	      double Xty[*p],
	      double lambda[*K],
	      double *L,
	      int *max_it,
	      double *eps,
	      double x0[*p],
	      double xk[*p],
	      int *penalty, // 1 = Lasso, 2 = Group-Lasso, 3 = Coop-Lasso, 4 = Tree-Group, 5 = Tree-Coop
	      int *i) {
  
  double delta = 2 * *eps ;
  double t0 = 1.0, tk;
  double s[*p];
  int j,l;

  // initializing y0
  for (j=0; j<*p; j++) {
    s[j] = x0[j];
  }
  
  while ( (delta > *eps) && (*i < *max_it) ) {
        
    // arg1 proximal: u = y0 - (XtX %*% y0 - Xty)/L
    for (j=0; j<*p; j++) {
      xk[j] = 0;
      for (l=0;l<*p; l++) {
	xk[j] += XtX[j][l] * s[l];
      }
      xk[j] = s[j] - (xk[j] - Xty[j]) / *L;
    }

    // updating xk
    switch(*penalty) {
      case 1: proximal_las(*p, lambda, *L, xk); 
	break;
      case 2: proximal_grp(*p, *K, pk, lambda, *L, xk);
	break;
      case 3: proximal_coo(*p, *K, pk, lambda, *L, xk);
	break;
      case 4: proximal_tree_grp(*p, *H, *K, Ks, pk, lambda, *L, xk);
	break;
      case 5: proximal_tree_coo(*p, *H, *K, Ks, pk, lambda, *L, xk);
    }

    // updating t
    tk = 0.5 * (1+sqrt(1+4* t0*t0)) ;
    // updating yk
    for (j=0; j<*p; j++) {
      s[j] = xk[j] + (t0-1.0)/tk * (xk[j] - x0[j]);
    }    

    // preparing next iterate
    delta = 0;
    for(j=0;j<*p;j++) {
      delta += pow(x0[j]-xk[j],2);
    }
    delta = sqrt(delta);
    t0 = tk;
    for(j=0;j<*p;j++) {
      x0[j] = xk[j];
    }
    *i += 1;

    R_CheckUserInterrupt();
  }
}

void ISTA_LRM(int *p,
	      int *H,
	      int *K,
	      int Ks[*H],
	      int pk[*K],
	      int *n,
	      double X[*p][*n],
	      double y[*n],
	      double lambda[*K],
	      double *L0,
	      int *max_it,
	      double *eps,
	      double x0[*p],
	      double xk[*p],
	      int *penalty, // 1 = Lasso, 2 = Group-Lasso, 3 = Coop-Lasso, 4 = Tree-Group, 5 = Tree-Coop
	      int *i) {
  
  double normDiff2, scalDiffd, fx0, fxk,  delta = 2 * *eps;
  double eta0[*n], etak[*n], df[*p] ;
  int j,k,l, found=0;
  double L = 2, maxL=*L0, ratio = 1.0;
  
  while ( (delta > *eps) && (*i < *max_it) ) {
    
    // eta0 = X * x0 + mu
    for (l=0; l<*n; l++) {
      eta0[l] = 0.0;
      for (j=0;j<*p;j++) {
	eta0[l] += X[j][l] * x0[j];
      }
    }
    
    // df = -t(X) (y - 1/(1+exp(eta0)))
    for (j=0; j<*p; j++) {
      df[j] = 0;
      for (l=0;l<*n;l++) {
	df[j] -= X[j][l] * (y[l] - 1/(1+exp(-eta0[l])));
      }
    }
    fx0=0;
    for (l=0;l<*n;l++) {
      fx0 -= y[l] * eta0[l] - log(1+exp(eta0[l]));
    }
    
    // recovering the previous feasible estimate of the Lipschitz constant
    L = fmin(ratio,L/2);
    
    // THE LINE SEARCH FOR L
    while (!found) {
      for (j=0; j<*p; j++) {
	xk[j] = x0[j] - df[j] / L;
      }

      // xk = proxi(x0 - nabla.f/L;lambda/L)
      switch(*penalty) {
        case 1: proximal_las(*p, lambda, L, xk); 
	  break;
        case 2: proximal_grp(*p, *K, pk, lambda, L, xk);
	  break;
        case 3: proximal_coo(*p, *K, pk, lambda, L, xk);
	  break;
        case 4: proximal_tree_grp(*p, *H, *K, Ks, pk, lambda, L, xk);
	  break;
        case 5: proximal_tree_coo(*p, *H, *K, Ks, pk, lambda, L, xk);
      }

      // etak = X * xk
      for (l=0; l<*n; l++) {
	etak[l] = 0.0;
	for (j=0;j<*p;j++) {
	  etak[l] += X[j][l] * xk[j];
	}
      }

      // normDiff2 = ||x0-xk||^2, 
      normDiff2 = 0.0;
      scalDiffd = 0.0;
      for(j=0;j<*p;j++) {
	normDiff2 += pow(x0[j] - xk[j],2);
	scalDiffd += (xk[j] - x0[j]) * df[j];
      }
      
      // fx0 = -sum (y * etak - log(1+exp(etak))
      fxk=0;
      for (l=0;l<*n;l++) {
	fxk -= y[l] * etak[l] - log(1+exp(etak[l]));
      }

      // ratio = 2 * (f.bin.xk - f.bin.x0 - t(xk-x0) %*% nabla.f) / normDiff2
      ratio = 2 * (fxk - fx0 - scalDiffd) / normDiff2 ;

      if (L >= ratio || L >= maxL) {
	found = 1;
      } else {
	L = fmin(fmax(2*L, ratio),maxL);
      }

      R_CheckUserInterrupt();

    } // end line search
    
    // preparing next iterate
    delta = sqrt(normDiff2);
    for(j=0;j<*p;j++) {
      x0[j] = xk[j];
    }
    found = 0;
    *i += 1;

    R_CheckUserInterrupt();
  }
}

void FISTA_LRM(int *p,
	       int *H,
	       int *K,
	       int Ks[*H],
	       int pk[*K],
	       int *n,
	       double X[*p][*n],
	       double y[*n],
	       double lambda[*K],
	       double *L0,
	       int *max_it,
	       double *eps,
	       double x0[*p],
	       double xk[*p],
	       int *penalty, // 1 = Lasso, 2 = Group-Lasso, 3 = Coop-Lasso, 4 = Tree-Group, 5 = Tree-Coop
	       int *i) {
  
  double normDiff2, scalDiffd, fy0, fyk, delta = 2 * *eps, t0=1.0, tk;
  double s[*p], eta0[*n], etak[*n], df[*p] ;
  int j,k,l, found=0; 
  double L = 2, maxL = *L0, ratio=1.0;

  // initializing y0
  for (j=0; j<*p; j++) {
    s[j] = x0[j];
  }
  
  while ( (delta > *eps) && (*i < *max_it) ) {
    
    // eta0 = X * y0
    for (l=0; l<*n; l++) {
      eta0[l] = 0;
      for (j=0;j<*p;j++) {
	eta0[l] += X[j][l] * s[j];
      }
    }
    
    // nabla.f = -t(X) (y - 1/(1+exp(eta0)))
    for (j=0; j<*p; j++) {
      df[j] = 0;
      for (l=0;l<*n;l++) {
	df[j] -= X[j][l] * (y[l] - 1/(1+exp(-eta0[l])));
      }
    }
    
    // fx0 = -sum (y * eta0 - log(1+exp(eta0))
    fy0=0;
    for (l=0;l<*n;l++) {
      fy0 -= y[l] * eta0[l] - log(1+exp(eta0[l]));
    }
    
    // recovering the previous feasible estimate of the Lipschitz constant
    L = fmin(ratio,L/2);

    // THE LINE SEARCH FOR L
    while (!found) {
      // u = x0 - nabla.f/L
      for (j=0; j<*p; j++) {
	xk[j] = s[j] - df[j] / L;
      }
      // xk = proxi(x0 - nabla.f/L;lambda/L)
      switch(*penalty) {
        case 1: proximal_las(*p, lambda, L, xk); 
	  break;
        case 2: proximal_grp(*p, *K, pk, lambda, L, xk);
	  break;
        case 3: proximal_coo(*p, *K, pk, lambda, L, xk);
	  break;
        case 4: proximal_tree_grp(*p, *H, *K, Ks, pk, lambda, L, xk);
	  break;
        case 5: proximal_tree_coo(*p, *H, *K, Ks, pk, lambda, L, xk);
      }

      // etak = X * xk
      for (l=0; l<*n; l++) {
	etak[l] = 0;
	for (j=0;j<*p;j++) {
	  etak[l] += X[j][l] * xk[j];
	}
      }

      // normDiff2 = ||x0-xk||^2, 
      normDiff2 = 0.0;
      scalDiffd = 0.0;
      for(j=0;j<*p;j++) {
	normDiff2 += pow(s[j] - xk[j],2);
	scalDiffd += (xk[j] - s[j]) * df[j] ;
      }

      // fx0 = -sum (y * etak - log(1+exp(etak))
      fyk=0;
      for (l=0;l<*n;l++) {
	fyk -= y[l] * etak[l] - log(1+exp(etak[l]));
      }

      // ratio = 2 * (f.bin.xk - f.bin.x0 - t(xk-x0) %*% nabla.f) / normDiff2
      ratio = 2 * (fyk - fy0 - scalDiffd) / normDiff2 ;

      if (L >= ratio || L >= maxL) {
	found = 1;
      } else {
	L = fmin(fmax(2*L, ratio),maxL);
      }

      R_CheckUserInterrupt();

    } // end line search

    // updating t
    tk = 0.5 * (1+sqrt(1+4* t0*t0)) ;
    // updating yk
    for (j=0; j<*p; j++) {
      s[j] = xk[j] + (t0-1.0)/tk * (xk[j] - x0[j]);
    }
    
    // preparing next iterate
    delta = sqrt(normDiff2);
    t0 = tk;
    for(j=0;j<*p;j++) {
      x0[j] = xk[j];
    }
    found = 0;
    *i += 1;

    R_CheckUserInterrupt();
  }
}
