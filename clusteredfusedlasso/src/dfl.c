#include <R.h>
#include <Rmath.h>

void dfl1D(int *np, double *y, int *Kp, double *b, double *lamp, int *g);


// Parameters:
// np - Pointer to the number of data points
// y - Pointer to the data 
// Kp - Pointer to the number of groups
// b - Pointer to the group levels 
// lamp - Pointer to lambda
// g - Pointer to the group assigments (this is empty, and will be 
//     filled by the data). This has the same length as y
void dfl1D(int *np, double *y, int *Kp, double *b, double *lamp, int *g) {
  int n=*np;
  int K=*Kp;
  double lam=*lamp;
  int i;

  // create a 2D array of back-pointers
  int **bp = malloc((n-1)*sizeof(int*));
  for (i=0; i<n-1; i++) {
    bp[i] = malloc(K*sizeof(int));
  }

  // create 2 1D arrays of min function values
  double *f = malloc(K*sizeof(double));
  double *ff = malloc(K*sizeof(double));
  for (i=0; i<K; i++) {
    ff[i] = 0;
  }

  int m,j;
  double a;

  // iterate over g's
  for (m=0; m<n-1; m++) {
    // iterate over values for next g
    for (i=0; i<K; i++) {
      f[i] = 0.5*(y[m]-b[0])*(y[m]-b[0]) + lam*fabs(b[i]-b[0]) + ff[0];
      bp[m][i] = 0;
      // iterate over values for the current g
      for (j=1; j<K; j++) {
	a = 0.5*(y[m]-b[j])*(y[m]-b[j]) + lam*fabs(b[i]-b[j]) + ff[j];
	if (a < f[i]) {
	  f[i] = a;
	  bp[m][i] = j;
	}
      }
    }

    // save the current min function values
    for (i=0; i<K; i++) {
      ff[i] = f[i];
    }
  }
    
  // now solve for the last g
  double amin = 0.5*(y[n-1]-b[0])*(y[n-1]-b[0]) + ff[0];
  g[n-1] = 0;
  for (i=1; i<K; i++) {
    a = 0.5*(y[n-1]-b[i])*(y[n-1]-b[i]) + ff[i];
    if (a < amin) {
      amin = a;
      g[n-1] = i;
    }
  }
  
  // and use our back-pointers to solve for the rest
  for (i=n-2; i>=0; i--) {
    g[i] = bp[i][g[i+1]];
  }

  // free allocated memory
  free(ff);
  free(f);
  for (i=0; i<n-1; i++) free(bp[i]);
  free(bp);
}
