/* 2.2 & 2.3
2.2 backsubstitution is just reducing the matrix to a upper triangular matrix
2.3 LU decomp reduces a matrix A to lower and upper triangular 
such that L * U = A
We can use that decomposition to solve 
A * x = (L * U) * x = L * (U * x) = b
by first solving L * y = b and then U * x = y
*/ 