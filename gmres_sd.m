function [x, error, its, flag] = gmres_sd( A, x, b, L, U, restrt, max_it, tol)
%GMRES_SD   Left-preconditioned GMRES in single/double precision
%   Solves Ax=b by solving the preconditioned linear system (LU)^{-1}Ax=(LU)^{-1}b
%   using the Generalized Minimal residual ( GMRES ) method.
%   Currently uses (preconditioned) relative residual norm to check for convergence 
%   (same as Matlab GMRES)
%   Single precision used throughout, except in applying (U\L\A) to a vector 
%   which is done in double precision
%
%   input   A        REAL nonsymmetric positive definite matrix
%           x        REAL initial guess vector
%           b        REAL right hand side vector
%           L        REAL L factor of lu(A)
%           U        REAL U factor of lu(A)
%           restrt   INTEGER number of iterations between restarts
%           max_it   INTEGER maximum number of iterations
%           tol      REAL error tolerance
%
%   output  x        REAL solution vector
%           error    REAL error norm
%           iter     INTEGER number of (inner) iterations performed
%           flag     INTEGER: 0 = solution found to tolerance
%                             1 = no convergence given max_it

flag = 0;
its = 0;

%Ensure single working precision
A = single(A);
b = single(b);
x = single(x);

%Cast half precision L and U factors to working precision
L = single(L);
U = single(U);


rtmp = b-A*x;

r = double(L)\double(rtmp);
r = double(U)\double(r);
r = single(r);

bnrm2 = norm(r );
if  ( bnrm2 == 0.0 ), bnrm2 = 1.0; end


error(1) = norm( r ) / bnrm2;
if ( error(1) < tol ) return, end

[n,~] = size(A);                                  % initialize workspace
m = restrt;
V(1:n,1:m+1) = single(zeros(n,m+1));
H(1:m+1,1:m) = single(zeros(m+1,m));
cs(1:m) = single(zeros(m,1));
sn(1:m) = single(zeros(m,1));
e1    = single(zeros(n,1));
e1(1) = single(1.0);

for iter = 1:max_it,                              % begin iteration
    rtmp = single(b-A*x);
    r = double(U)\(double(L)\double(rtmp));
    r = single(r);
    
    V(:,1) = r / norm( r );
    s = norm( r )*e1;
    for i = 1:m,                     % construct orthonormal basis via GS
        its = its+1;
        vcur = V(:,i);      
        
        vcur = double(U)\(double(L)\(double(A)*double(vcur)));
        
        w = single(vcur);
      
        for k = 1:i,
            H(k,i)= w'*V(:,k);
            w = w - H(k,i)*V(:,k);
        end
        H(i+1,i) = norm( w );
        V(:,i+1) = w / H(i+1,i);
        for k = 1:i-1,                              % apply Givens rotation
            temp     =  cs(k)*H(k,i) + sn(k)*H(k+1,i);
            H(k+1,i) = -sn(k)*H(k,i) + cs(k)*H(k+1,i);
            H(k,i)   = temp;
        end
        [cs(i),sn(i)] = rotmat( H(i,i), H(i+1,i) ); % form i-th rotation matrix
        temp   = cs(i)*s(i);                        % approximate residual norm
        s(i+1) = -sn(i)*s(i);
        s(i)   = temp;
        H(i,i) = cs(i)*H(i,i) + sn(i)*H(i+1,i);
        H(i+1,i) = 0.0;
        
        error((iter-1)*m+i+1)  = abs(s(i+1)) / bnrm2;
        if ( error((iter-1)*m+i+1) <= tol ),                        % update approximation
            y = H(1:i,1:i) \ s(1:i);                 % and exit
            addvec = V(:,1:i)*y;
            x = x + addvec;
            break;
        end
    end
    
    if ( error(end) <= tol ), break, end
    y = H(1:m,1:m) \ s(1:m);
    addvec = V(:,1:m)*y;
    x = x + addvec;                            % update approximation
    rtmp = b-A*x;
    r = double(U)\(double(L)\double(rtmp));           % compute residual
    r = single(r);
    s(i+1) = norm(r);
    error = [error,s(i+1) / bnrm2];                        % check convergence
    if ( error(end) <= tol ), break, end;
end

if ( error(end) > tol ) flag = 1; end;                 % converged




function [ c, s ] = rotmat( a, b )

%
% Compute the Givens rotation matrix parameters for a and b.
%
if ( b == 0.0 ),
    c = 1.0;
    s = 0.0;
elseif ( abs(b) > abs(a) ),
    temp = a / b;
    s = 1.0 / sqrt( 1.0 + temp^2 );
    c = temp * s;
else
    temp = b / a;
    c = 1.0 / sqrt( 1.0 + temp^2 );
    s = temp * c;
end