// Curve Fit

// g++ main.cpp -o executable_name

// See problemDescription/curveFit.pdf

// This program uses Gaussian elimination with and without (partial) pivoting to solve
// SUM_(j=0)^(n) M(i,j) A(j) = B(i), where M(i,j) = 1/(i+j+1) and
// B(i) = int_0^1 dx x^i sin(x), to give a least squares n degree polynomial
// approximation sin(x) = SUM_(i=0)^(n) A(i) x^i on the interval x = 0 to 1.

// const and pass by reference?

#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

// can probably combine these two
void gaussNoPivot(vector <vector <double> > M, vector <double>& A, vector <double> B, int dim);
void gaussPivot(vector <vector <double> > M, vector <double>& A, vector <double> B, int dim);
void sub(vector <vector <double> > M, vector <double>& A, vector <double> B, int dim);

int main() {
  int n = 10;       // polynomial degree
  int dim = n + 1;  // dimension of arrays

  // Set B arrays
  vector <double> Bnp;
  vector <double> Bp;

  Bnp.resize(dim);
  Bp.resize(dim);

  Bnp[0] = 1.0 - cos(1.0);
  Bp[0] = Bnp[0];

  Bnp[1] = sin(1.0) - cos(1.0);
  Bp[1] = Bnp[1];
 
  for (int i = 2; i < dim; i++) {
    Bnp[i] = -cos(1.0) + ((double) i)*sin(1.0) - ((double) i)*((double) i - 1.0)*Bnp[i-2];
    Bp[i] = Bnp[i];
  }

  // Set M arrays
  vector <vector <double> > Mnp;
  vector <vector <double> > Mp;

  Mnp.resize(dim);
  Mp.resize(dim);
  
  for (int i = 0; i < dim; i++) {
    Mnp[i].resize(dim);
    Mp[i].resize(dim);
    for (int j = 0; j < dim; j++) {
      Mnp[i][j] = 1.0/((double) i + (double) j + 1.0);
      Mp[i][j] = Mnp[i][j];
    }
  }

  // Determine A arrays
  vector <double> Anp;
  vector <double> Ap;

  Anp.resize(dim);
  Ap.resize(dim);

  gaussNoPivot(Mnp, Anp, Bnp, dim);

  cout << "A without pivoting:" << endl;
  for (int i = 0; i < dim; i++) {
    cout << "Anp[" << i << "] = " << Anp[i] << endl;
  }

  cout << "sin(x) = ";
  for (int i = 0; i <= n; i++) {
    cout << Anp[i] << " * x^" << i << " + ";
  }

  cout << endl;

/*
  gaussPivot(Mp, Ap, Bp, dim);

  cout << "A with pivoting:" << endl;
  for (int i = 0; i < dim; i++) {
    cout << "Ap[" << i << "] = " << Ap[i] << endl;
  }

  cout << "sin(x) = ";
  for (int i = 0; i <= n; i++) {
    cout << Ap[i] << " * x^" << i << " + ";
  }

  cout << "Difference in A with and without pivoting:" << endl;
  for (int i = 0; i < dim; i++) {
    cout << "Anp[" << i << "] - Ap[" << i << "] = " << Anp[i] - Ap[i] << endl;
  }
*/

  return 0;
}


void sub(vector <vector <double> > M, vector <double>& A, vector <double> B, int dim) {
  int i;
  double sum;
  A[dim-1] = B[dim-1]/M[dim-1][dim-1];

  // Loop backwards, starting at bottom row.
  for (int i = dim-2; i >= 0; i--) { 
    sum = B[i];
    for (int j = i+1; j < dim; j++) {
      sum = sum - M[i][j]*A[j];
    }
    A[i] = sum/M[i][i];
  }
}

/* FUNCTION FOR GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING	
! Subroutine for Gaussian elimination with partial pivoting
	! --------------------------------------------------
	SUBROUTINE gausspivot(M, A, B, dim)
		implicit none
		integer dim, maxdim, i, j, k
		integer maxindex
		parameter (maxdim=50)
		double precision M(maxdim,maxdim), A(maxdim), B(maxdim)
		double precision maxvalue, swap

		DO k = 1, (dim-1)					! ELIMINATION
			maxvalue = 0.d0
			DO i = k, dim
				IF (M(i,k).gt.maxvalue) THEN
					maxvalue = M(i,k)
					maxindex = i
				END IF
			END DO
			IF (maxindex.ne.k) THEN
				DO j = 1, dim
					swap = M(maxindex,j)
					M(maxindex,j) = M(k,j)
					M(k,j) = swap
				END DO
				swap = B(maxindex)
				B(maxindex) = B(k)
				B(k) = swap
			END IF
			DO i = (k+1), dim				! Loops through rows.
				M(i,k) = M(i,k)/M(k,k)		! This would be zero, so we'll store ratio here.
				DO j = (k+1), dim			! Loops through columns.
					M(i,j) = M(i,j) - M(i,k)*M(k,j)
				END DO
				B(i) = B(i) - M(i,k)*B(k)
			END DO
		END DO
		CALL sub(M, A, B, dim)				! SUBSTITUTION
	RETURN
	END
*/

void gaussNoPivot(vector <vector <double> > M, vector <double>& A, vector <double> B, int dim) {
  for (int k = 0; k < dim; k++) {       // Elimination
    for (int i = k+1; i < dim; i++) {   // Loops through rows
      M[i][k] = M[i][k]/M[k][k];        // This would be zero, so we'll store ratio here.
      for (int j = k+1; j < dim; j++) { // Loops through columns.
        M[i][j] = M[i][j] - M[i][k]*M[k][j];
      }
      B[i] = B[i] - M[i][k]*B[k];
    }
  }
  sub(M, A, B, dim);
}
