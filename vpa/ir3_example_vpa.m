%IR3_EXAMPLE_VPA    Example script for running iterative refinement with 3 precisions
%   Generates data for plots in Figure 10.1 of MIMS Eprint 2017.# 
%   Uses vpa for extended precision computations

n = 100;
maxit = 10;

rng(1);
A = gallery('randsvd',n,1e9,2);
b = randn(n,1);


%Run standard IR with uf = single, u = double, ur = double
fprintf('Running LU-IR\n'); 
sir3_vpa(A,b,1,2,2,maxit);
drawnow

%Run standard IR with uf = single, u = double, ur = quad
fprintf('Running LU-IR\n'); 
sir3_vpa(A,b,1,2,4,maxit);
drawnow

%Run standard IR with uf = double, u = double, ur = quad
fprintf('Running LU-IR\n'); 
sir3_vpa(A,b,2,2,4,maxit);
drawnow

%Run GMRES-IR with uf = single, u = double, ur = quad
fprintf('Running GMRES-IR\n'); 
gmresir3_vpa(A,b,1,2,4,maxit,1e-6);
drawnow