#Adaptive truncated SVD of A for the given error tolerance. (adpRandSVD.m)  
Syntax:  
&#8195 s= adpRandSVD(A, tol)  
 s= adpRandSVD(A, tol, P)  
 s= adpRandSVD(A, tol, P, bsize)
 [U, S, V]= adpRandSVD(A, tol)
 [U, S, V]= adpRandSVD(A, tol, P)
 [U, S, V]= adpRandSVD(A, tol, P, bsize)
 [Q, B]= adpRandSVD(A, tol)
 [Q, B]= adpRandSVD(A, tol, P)
 [Q, B]= adpRandSVD(A, tol, P, bsize)
 -tol is the error tolerance for norm(A-U*S*V')/norm(A). We measure the
  Frobenius norm. tol is a value less than 1, like 0.1, 0.85.
 -P is an optional parameter to balance time and accuacy (default value 2).
  With large P, the accuracy increases with runtime overhead.
 -bsize is the number of rows/columns processed as a block (default value
  20).
 -Q and B forms a rank-k approximation of A, where Q is an mxk orthonormal 
  matrix, and B is a kxn matrix (A is of mxn).
Algorithm: a blocked randomized algorithm randQB_EI for auto-rank problem.
  The resulting low-rank approximation guarantees the error is within tol.
  To avoid large computing cost and inaccuracy, it may terminate after a 
  maximum rank is attained. 
                  --by Wenjian Yu, Nov. 2016