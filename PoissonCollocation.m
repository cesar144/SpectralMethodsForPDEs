%Poisson Problem on [-1,1]x[-1,1] with u=0 on boundary
% Set up grids and tensor product Laplacian and solve for u:
N = 64; x = LobattoGaussLegendre(N)';
D = LGLDiffMtrx(N);
y = x;
[xx,yy] = meshgrid(x(2:N),y(2:N));
xx = xx(:); yy = yy(:);
% stretch 2D grids to 1D vectors
fun = @(x, y) sin(x);
f = 10*sin(8*xx.*(yy-1));
tic
D2 = D^2; D2 = D2(2:N,2:N); I = eye(N-1);
L = kron(I,D2) + kron(D2,I);
% Laplacian
figure(1), clf, spy(L), drawnow
u = L\f; toc
% solve problem and watch the clock
% Reshape long 1D results onto 2D grid:
uu = zeros(N+1,N+1); uu(2:N,2:N) = reshape(u,N-1,N-1);
[xx,yy] = meshgrid(x,y);
v = uu(N/4+1,N/4+1);
% Interpolate to finer grid and plot:
[xxx,yyy] = meshgrid(-1:.04:1,-1:.04:1);
uuu = interp2(xx,yy,uu,xxx,yyy,'spline');
figure(2), clf, mesh(xxx,yyy,uuu), colormap([0 0 0])
text(.4,-.3,-.3,sprintf('u(2^{-1/2},2^{-1/2}) = %14.11f',v))
xlabel x, ylabel y, zlabel u