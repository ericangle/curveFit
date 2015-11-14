#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

// can probably combine these two
void gaussNoPivot(vector <vector <double> > M, vector <double>& A, vector <double> B, int dim);
void gaussPivot(vector <vector <double> > M, vector <double>& A, vector <double> B, int dim);
void sub(vector <vector <double> > M, vector <double>& A, vector <double> B, int dim);

int main() {
  int n = 12;       // polynomial degree
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
  cout << endl;

  gaussPivot(Mp, Ap, Bp, dim);

  cout << "A with pivoting:" << endl;
  for (int i = 0; i < dim; i++) {
    cout << "Ap[" << i << "] = " << Ap[i] << endl;
  }

  cout << "sin(x) = ";
  for (int i = 0; i <= n; i++) {
    cout << Ap[i] << " * x^" << i << " + ";
  }

  cout << endl;
  cout << endl;

  cout << "Percent difference in A with and without pivoting:" << endl;
  for (int i = 0; i < dim; i++) {
    cout << "i = " << i << ": " << 100.0*abs(1.0 - Ap[i]/Anp[i]) << endl;
  }

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

void gaussPivot(vector <vector <double> > M, vector <double>& A, vector <double> B, int dim) {
  double maxvalue, swap;
  int maxindex;
  for (int k = 0; k < dim; k++) {  // Elimination
    maxvalue = 0.0;
    for (int i = k; i < dim; i++) {
      if (M[i][k] > maxvalue) {
        maxvalue = M[i][k];
        maxindex = i;
      }
    }
    if (maxindex != k) {
      for(int j = 0; j < dim; j++) {
        swap = M[maxindex][j];
	M[maxindex][j] = M[k][j];
	M[k][j] = swap;
      }
      swap = B[maxindex];
      B[maxindex] = B[k];
      B[k] = swap;
    }
    for (int i = k+1; i < dim; i++) {  // Loops through rows. go to NoPivot now?
      M[i][k] = M[i][k]/M[k][k];        // This would be zero, so we'll store ratio here.
      for (int j = k+1; j < dim; j++) { // Loops through columns.
        M[i][j] = M[i][j] - M[i][k]*M[k][j];
      }
      B[i] = B[i] - M[i][k]*B[k];
    }
  }
  sub(M, A, B, dim);  // substitution
}

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
  sub(M, A, B, dim);  // substitution
}
