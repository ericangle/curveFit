// Curve Fit

// g++ main.cpp -o executable_name

// Note about PDF explaining background

// This program ...

#include <iostream>
using namespace std;

int main() {
  cout << "Hello world." << endl;
  return 0;
}

/*
	! LEAST.F
	!
	! This program uses Gaussian elimination with and without (partial) pivoting to solve
	! SUM_(j=1)^(n+1) M(i,j) A(j) = B(i), where M(i,j) = 1/(i+j-1) is the Hilbert matrix and
	! B(i) = int_0^1 dx x^i sin(x), to give a least squares n degree polynomial
	! approximation sin(x) = SUM_(i=1)^(n+1) A(i) x^(i-1) on the interval x = 0 to 1.
	! ---------------------------------------------------------------------

	PROGRAM least

	implicit none

	integer i, j, k						! Loop variables
	integer n, dim						! dim = n + 1, n = polynomial degree
	integer maxdim					! Maximum potential value of dim
	parameter (maxdim=50)
	double precision Mnp(maxdim,maxdim), Anp(maxdim), Bnp(maxdim)	! No pivot arrays
	double precision Mp(maxdim,maxdim), Ap(maxdim), Bp(maxdim)		! Pivot arrays
	
	PRINT *, 'Polynomial degree?'
	READ *, n
	dim = n + 1

	Bnp(1) = 1.d0 - dcos(1.d0)			! With Bnp(1) & Bnp(2), Bnp(k) for k > 2 is given
	Bnp(2) = dsin(1.d0) - dcos(1.d0)		! recursively in the below DO loop.

	DO i = 3, dim						! Defines elements of Bnp(i>2)
		Bnp(i) = - dcos(1.d0) + dble(i-1)*dsin(1.d0) - dble(i*i-3*i+2)*Bnp(i-2)
	END DO

	DO i = 1, dim						! Defines elements of Mnp
		DO j = 1, dim
			Mnp(i,j) = 1.d0/dble(i+j-1)
			Mp(i,j) = Mnp(i,j)			! Creates duplicate M for use with pivoting
		END DO
		Bp(i) = Bnp(i)					! Creates duplicate B for use with pivoting
	END DO

	CALL gaussnopivot(Mnp, Anp, Bnp, dim)	! Perform Gaussian elimination without pivoting
	PRINT *, ' '
	PRINT *, 'A(i = 1 to',dim,') without pivoting'
	DO i = 1, dim						! Print Anp
		WRITE (*,10) Anp(i)
10		FORMAT(1PG22.14)
	END DO

	CALL gausspivot(Mp, Ap, Bp, dim)		! Perform Gaussian elimination with pivoting
	PRINT *, ' '
	PRINT *, 'A(i = 1 to',dim,') with pivoting'
	DO i = 1, dim						! Print Ap
		WRITE (*,11) Ap(i)
11		FORMAT(1PG22.14)
	END DO

	PRINT *, ' '
	PRINT *, 'Difference in A(i = 1 to',dim,') with and without pivoting'
	DO i = 1, dim						! Print Anp
		WRITE (*,12) Anp(i) - Ap(i)
12		FORMAT(1PG22.14)
	END DO

	PRINT *, ' '

	STOP
	END

	! Subroutine for Gaussian elimination without pivoting
	! --------------------------------------------------
	SUBROUTINE gaussnopivot(M, A, B, dim)
		implicit none
		integer dim, maxdim, i, j, k
		parameter (maxdim=50)
		double precision M(maxdim,maxdim), A(maxdim), B(maxdim)

		DO k = 1, (dim-1)					! ELIMINATION
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

	! Subroutine for substitution
	! --------------------------------------------------
	SUBROUTINE sub(M, A, B, dim)
		implicit none
		integer dim, maxdim, i, j, k
		parameter (maxdim=50)
		double precision M(maxdim,maxdim), A(maxdim), B(maxdim), sum

		A(dim) = B(dim)/M(dim,dim)
		DO k = 1, (dim-1)
			i = dim - k					! Loop backwards, starting at bottom row.
			sum = B(i)
			DO j = (i+1),dim
				sum = sum - M(i,j)*A(j)
			END DO
			A(i) = sum/M(i,i)
		END DO
	RETURN
	END
*/ 
