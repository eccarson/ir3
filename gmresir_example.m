%GMRESIR_EXAMPLE    Example script for comparing LU-IR and GMRES-IR (with 2 precisions)
%   Generates data for plots in Figures 5.1 and 5.2 of MIMS EPrint 2017.12
%   Note: Requires Advanpix multiprecision toolbox

n = 100;
maxit = 10;

condnums_single = [1e7,1e8,1e9,1e10];
condnums_double = [1e15,1e16,1e17,1e18];

%Run single/double tests
for i = 1:numel(condnums_single)
    
    rng(1);
    A = gallery('randsvd',n,condnums_single(i),3);
    b = randn(n,1);
    
    %Run standard IR with uf = u = single, ur = double
    fprintf('Running LU-IR\n');
    sir3(A,b,1,1,2,maxit);
    drawnow
    
    %Run GMRES-IR with uf = u = single, ur = double
    fprintf('Running GMRES-IR\n'); 
    gmresir3(A,b,1,1,2,maxit,1e-4);
    drawnow
    
end

%Run double/quad tests
for i = 1:numel(condnums_double)
    
    rng(1);
    A = gallery('randsvd',n,condnums_double(i),3);
    b = randn(n,1);
    
    %Run standard IR with uf = u = double, ur = quad
    fprintf('Running LU-IR\n');
    sir3(A,b,2,2,4,maxit);
    drawnow
    
    %Run GMRES-IR with uf = u = double, ur = quad
    fprintf('Running GMRES-IR\n'); 
    gmresir3(A,b,2,2,4,maxit,1e-4);
    drawnow
    
end
