function [x,w] = LobattoGaussLegendre(N)
% LOBATTO-GAUSS-LEGENDRE weights and nodes
% nodes x (Legendre points) and weights w for quadrature
% Node x_j corresponds to weight w_j
alpha = sqrt(((1:N-2).*((1:N-2)+2))./ (((1:N-2)+1/2).*((1:N-2)+3/2)));
A = sparse(1:N-2, 2:N-1, alpha, N-1 ,N-1) + sparse(2:N-1, 1:N-2, alpha, N-1, N-1);
[w,x] = eigs(A,N);
x = diag(x)/2;
x = x';
aux =  sqrt(3) ./ (w(1,1:N-1)*2);
aux = ones(N-1,N-1) .* aux;
w = w .* aux;
w = w.^2;
w = ((1-x.^2).^(-1)) ./ sum(w);
w = [2/(N*(N+1)), w, 2/(N*(N+1))];
x = [-1, x, 1];
[x ,s] = sort(x, 'descend');
w = w(s)';
end
