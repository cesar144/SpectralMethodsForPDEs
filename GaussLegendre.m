function [x,w] = GaussLegendre(N)
% GAUSS-LEGENDRE weights and nodes
% nodes x (Legendre points) and weights w for quadrature
% Node x_j corresponds to weight w_j
alpha = (1:N-1) ./ sqrt(4*(1:N-1).^2 -1);
A = sparse(1:N-1, 2:N, alpha, N ,N) + sparse(2:N, 1:N-1, alpha, N, N);
[w,x] = eigs(A,N);
x = diag(x);
aux =  1 ./ (w(1,1:N)*sqrt(2));
aux = ones(N,N) .* aux;
w = w .* aux;
w = w.^2;
w = 1./sum(w);
end