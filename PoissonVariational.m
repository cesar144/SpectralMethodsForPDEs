%Poisson Problem on [-1,1]x[-1,1] with u=0 on boundary
N = 64; 
%Nodes and Weights
[x,w] = LobattoGaussLegendre(N);
x = x'; w = w';
%Differentiation Matrix
D = LGLDiffMtrx(N);
% stretch 2D grids to 1D vectors
nPts = (N-1)^2;
Aijv = NaN(((N-1)^2)*(2*N-3), 3); %contains the indexes and values of the matrix A with non−zero values
y = x;
[xx,yy] = meshgrid(x(2:N),y(2:N));
xx = xx(:); yy = yy(:);
f = 10*sin(8*xx.*(yy-1));
r = w;
[ww, rr] = meshgrid(w(2:N),r(2:N));
ww = ww(:); rr = rr(:);
tic
% Computute mass matrix
M = ww .*rr;
b = M .* f; 
% Compute stiffness matrix
aux = 0;
for pt = 1 : nPts
    [i,j] = ind2sub([N-1, N-1], pt);
    aux = aux + 1;
    dx = D(:,i+1);
    dy = D(:, j+1);
    Aijv(aux,:) = [pt, pt, w(j+1)*(w*dx.^2) + w(i+1)*(w*dy.^2)];
    for s = 1:N-1
        if s ~= i
            aux  = aux + 1;
            auxd = D(:, s+1);
            auxpt = sub2ind([N-1, N-1], s, j);
            Aijv(aux,:) = [pt, auxpt,w(j+1)*(w*(auxd.*dx))];
        end
    end
    for t = 1:N-1
        if  t ~= j
            aux  = aux + 1;  
            auxd = D(:, t+1);
            auxpt = sub2ind([N-1, N-1], i, t);
            Aijv(aux,:) = [pt, auxpt, w(i+1)*(w*(auxd.*dy))];
        end
    end
end
% Solving
A = sparse(Aijv(1:aux,1),Aijv(1:aux,2),Aijv(1:aux,3));
figure(1), clf, spy(A), drawnow
u = -A\b; toc
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