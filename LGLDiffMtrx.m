function D = LGLDiffMtrx(N)
% LOBATTO GAUSS LEGENDRE Diferentation Matrix
x = LobattoGaussLegendre(N)';
X = repmat(x,1,N+1);
dX = X-X';
dX(logical(eye(N+1))) = NaN;
a = prod(dX', 'omitnan')';
a = repmat(a, 1, N+1);
a = a.*((1./a)');
dX = 1 ./dX;
D = dX .* a;
d = sum(dX' , 'omitnan');
D(logical(eye(N+1))) = d';
end