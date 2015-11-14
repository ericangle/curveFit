.. role:: raw-math(raw)
    :format: latex html

curveFit
=============

This program uses Gaussian elimination with and without (partial) pivoting to solve
SUM_(j=0)^(n) M(i,j) A(j) = B(i), where M(i,j) = 1/(i+j+1) and
B(i) = int_0^1 dx x^i sin(x), to give a least squares n degree polynomial
approximation sin(x) = SUM_(i=0)^(n) A(i) x^i on the interval x = 0 to 1.

Run Tests
------------------

To run the tests, ...
