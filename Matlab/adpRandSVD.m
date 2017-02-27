function [Q, B, V] = adpRandSVD(A, tol, p, block_size)
%Adaptive truncated SVD of A for the given error tolerance. 
%Suitable for a dense/sparse A, can be faster than the randQB for sparse A.
%Syntax:
% s= adpRandSVD(A, tol)
% s= adpRandSVD(A, tol, P)
% s= adpRandSVD(A, tol, P, bsize)
% [U, S, V]= adpRandSVD(A, tol)
% [U, S, V]= adpRandSVD(A, tol, P)
% [U, S, V]= adpRandSVD(A, tol, P, bsize)
% [Q, B]= adpRandSVD(A, tol)
% [Q, B]= adpRandSVD(A, tol, P)
% [Q, B]= adpRandSVD(A, tol, P, bsize)
% -tol is the error tolerance for norm(A-U*S*V')/norm(A). We measure with
%  Frobenius norm. tol is a value less than 1, like 0.1, 0.5.
% -P is an optional parameter to balance time and accuacy (default value 1).
%  With large P, the accuracy increases with runtime overhead.
% -bsize is the number of rows/columns processed as a block (default value
%  20).
% -Q and B form a rank-k approximation of A, where Q is an mxk orthonormal 
%  matrix, and B is a kxn matrix (A is of mxn).
%Algorithm: a blocked randomized algorithm randQB_EI for the fixed-precision 
%  problem of low-rank matrix approximation.
%  The result guarantees the relative error of approximation is within tol.
%  To avoid large computing cost and inaccuracy, it may terminate after a 
%  maximum rank is attained. 
%                  --by Wenjian Yu, Nov. 2016

if nargin<2,
    disp('Less than 2 input arguments');
    return;
elseif nargin<3
    p= 1;
    block_size=20;
elseif nargin<4
    block_size=20;
end

    [m, n]  = size(A);
    
    Q = zeros(m, 0);
    B = zeros(0, n);
    E= norm(A, 'fro')^2; E0= E;
    threshold= tol^2*E;
    maxiter= ceil(min(m,n)/3/block_size);   % avoid so many iterations ?
    flag= false;

    for i=1:maxiter,
        Omg = randn(n, block_size);           % nxb matrix
        q = A * Omg - (Q * (B * Omg));        % mxb matrix
        [q, ~] = qr(q, 0);
        
        for j = 1:p
            [Omg, ~] = qr(A'*q - B'*(Q'*q), 0);  % yuwj corrected it
            [q, ~] = qr(A*Omg - Q*(B*Omg), 0);
        end
        
        [q, ~] = qr(q - Q * (Q' * q), 0);
        
        %b = q' * A;
        %b = q' * A - q' * Q * B;
        Omg= A'*q- B'*(Q'*q);  % b', a nxb matrix
        
        % variant
        % bv = qv' * A - qv' * Qv * Bv;
        
        Q = [Q, q];
        B = [B; Omg'];
        
        temp = E- norm(Omg, 'fro')^2;
        
        if temp< threshold,
            for j=1:block_size,
                E= E-norm(Omg(:, j))^2;
                if E< threshold,
                    flag= true;
                    break;
                end
            end
        else
            E= temp;
        end
        if flag,
            k= (i-1)*block_size+j;
            break;
        end
    end
    
    if ~flag,
        k= i*block_size;
        disp(sprintf('Error = %f. Fail to attain the accuracy demand within rank %d!', sqrt(E/E0), k));
    end

if nargout==2,
    Q= Q(:, 1:k);
    B= B(1:k, :);
elseif nargout==1,
    Q= svd(B, 'econ');
    Q= Q(1:k);
elseif nargout==3,
    [U1, B, V]= svd(B, 'econ');
    V= V(:, 1:k);
    B= diag(B);
    B= B(1:k);
    Q= Q*U1(:,1:k);
end

end
