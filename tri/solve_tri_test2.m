clear all; close all; clc;

%% Specify parameters and problem
% Koornwinder poly params (unit weight)
a = 1/2; b = 1/2; c = 1/2;
% legendre analog quadrature rule on triangle
load 'triquadLeg_16_27.mat';
X = Zk(1:N); Y = Zk(N+1:2*N); W = Zk(2*N+1:3*N);
% specify problem
%%%syms x y
%%%U(x,y) = x.*y.*(1-x-y).*(x+y).^8;
%U(x,y) = x.*y.*(1-x-y).*exp(y.*sin(x+y));
%U(x,y) = x.*y.*(1-x-y).*exp(-(x.^2+y.^2));
%%%dxU(x,y) = diff(U,'x',1);
%%%dyU(x,y) = diff(U,'y',1);
%%%lapU(x,y) = diff(dxU,'x',1)+diff(dyU,'y',1);
U = @(x,y) x.*y.*(1-x-y).*(x+y).^8;
%U(x,y) = x.*y.*(1-x-y).*exp(y.*sin(x+y));
%U(x,y) = x.*y.*(1-x-y).*exp(-(x.^2+y.^2));
dxU = @(x,y) -(y.*(x + y).^7.*(10*x.^2 + (-1 + y).*y + x.*(-9 + 11*y)));
dxxU = @(x,y) -2*y.*(x + y).^6.*(45*x.^2 + 18*x.*(-2 + 3*y) + y.*(-8 + 9*y));
dyU = @(x,y) -(x.*(x + y).^7.*(x.^2 + y.*(-9 + 10*y) + x.*(-1 + 11*y)));
dyyU = @(x,y) -2*x.*(x + y).^6.*(9*x.^2 + 9.*y.*(-4 + 5*y) + x.*(-8 + 54*y));
lapU = @(x,y) dxxU(x,y) + dyyU(x,y);

% convert to matlab funcs
f = @(x,y) lapU(x,y);%matlabFunction(lapU);
u = @(x,y) U(x,y);%matlabFunction(U);

%% Precompute

% highest poly degree is m-1, and there are M total polys
m = 11; n = m+1; M = m*(m+1)/2; 

% normalization under (a,b,c), (a+1,b,c) etc.
H_abc = structure_factors_tri(n+1,a,b,c);
H_a1bc = structure_factors_tri(n+1,a+1,b,c);
H_a1b1c = structure_factors_tri(n+1,a+1,b+1,c);
H_a1b1c1 = structure_factors_tri(n+1,a+1,b+1,c+1);
H_ab1c = structure_factors_tri(n+1,a,b+1,c);
H_abc1 = structure_factors_tri(n+1,a,b,c+1);
H_a1bc1 = structure_factors_tri(n+1,a+1,b,c+1);
H_ab1c1 = structure_factors_tri(n+1,a,b+1,c+1);
H_a2bc2 = structure_factors_tri(n+1,a+2,b,c+2);
H_a2b1c2 = structure_factors_tri(n+1,a+2,b+1,c+2);
H_ab2c2 = structure_factors_tri(n+1,a,b+2,c+2);
H_a1b2c2 = structure_factors_tri(n+1,a+1,b+2,c+2);
H_a2b2c2 = structure_factors_tri(n+1,a+2,b+2,c+2);
H_a2b1c1 = structure_factors_tri(n+1,a+2,b+1,c+1);
H_a2b2c1 = structure_factors_tri(n+1,a+2,b+2,c+1);

% vandermonde under (a,b,c) 
V_abc = jPoly_tri(X,Y,H_abc,n-1,a,b,c);
% weighted vandermonde under (a+1,b+1,c+1)
V_a1b1c1w = jPoly_tri_weighted(X,Y,H_a1b1c1(1:n,1:n),n-2,a+1,b+1,c+1);

% promotion matrices (a+2,b,c+2) -> (a+2,b+1,c+2) etc.
K_a2bc2_a2b1c2 = promotion_mat_tri(a+2,b,c+2,H_a2bc2,H_a2b1c2,1);
K_a2b1c2_a2b2c2 = promotion_mat_tri(a+2,b+1,c+2,H_a2b1c2,H_a2b2c2,1);
K_ab2c2_a1b2c2 = promotion_mat_tri(a,b+2,c+2,H_ab2c2,H_a1b2c2,0);
K_a1b2c2_a2b2c2 = promotion_mat_tri(a+1,b+2,c+2,H_a1b2c2,H_a2b2c2,0);
K_abc_a1bc = promotion_mat_tri(a,b,c,H_abc,H_a1bc,0);
K_a1bc_a1b1c = promotion_mat_tri(a+1,b,c,H_a1bc,H_a1b1c,1);
K_a1b1c_a1b1c1 = promotion_mat_tri(a+1,b+1,c,H_a1b1c,H_a1b1c1,2);
K_a1bc1_a1b1c1 = promotion_mat_tri(a+1,b,c+1,H_a1bc1,H_a1b1c1,1);
K_ab1c1_a1b1c1 = promotion_mat_tri(a,b+1,c+1,H_ab1c1,H_a1b1c1,0);
K_a1b1c1_a2b1c1 = promotion_mat_tri(a+1,b+1,c+1,H_a1b1c1,H_a2b1c1,0);
K_a2b1c1_a2b2c1 = promotion_mat_tri(a+2,b+1,c+1,H_a2b1c1,H_a2b2c1,1);
K_a2b2c1_a2b2c2 = promotion_mat_tri(a+2,b+2,c+1,H_a2b2c1,H_a2b2c2,2);

% unweighted derivative matrices
Dx_abc_a1bc1 = D1_tri(a,b,c,H_abc,H_a1bc1,0);
Dx_a1bc1_a2bc2 = D1_tri(a+1,b,c+1,H_a1bc1,H_a2bc2,0);
Dy_abc_ab1c1 = D1_tri(a,b,c,H_abc,H_ab1c1,1);
Dy_ab1c1_ab2c2 = D1_tri(a,b+1,c+1,H_ab1c1,H_ab2c2,1);

% weighted derivative amatrices
Wx_a1b1c1_ab1c = D1_tri_weighted(a+1,b+1,c+1,H_a1b1c1,H_ab1c,0);
Wy_a1b1c1_a1bc = D1_tri_weighted(a+1,b+1,c+1,H_a1b1c1,H_a1bc,1);

% lowering matrices
Lx_ab1c_abc = lowering_mat_tri(a,b+1,c,H_ab1c,H_abc,1);
Ly_a1bc_abc = lowering_mat_tri(a+1,b,c,H_a1bc,H_abc,0);

% laplacian in weighted basis (3/2,/3/2,3/2)
Dxx = K_a1bc1_a1b1c1*Dx_abc_a1bc1*Lx_ab1c_abc*Wx_a1b1c1_ab1c; 
Dyy = K_ab1c1_a1b1c1*Dy_abc_ab1c1*Ly_a1bc_abc*Wy_a1b1c1_a1bc;
% remove zero rows and columns of weighted laplacian
Lap_a1b1c1w_a1b1c1 = Dxx(1:M,1:M) + Dyy(1:M,1:M);

%% solve with weighted laplacian
% coefficients of RHS for poisson in (1/2,1/2,1/2)
cf_abc = V_abc(:,1:M)'*(f(X,Y).*W);
% promote coefficients to (3/2,3/2,3/2)
cf_a1b1c1 = K_a1b1c_a1b1c1(1:M,1:M)*(K_a1bc_a1b1c(1:M,1:M)*(K_abc_a1bc(1:M,1:M)*cf_abc));
% solution coefficients
cu_a1b1c1w = Lap_a1b1c1w_a1b1c1\cf_a1b1c1;
% 2-norm relative error
fprintf('weighted laplacian error: \t %1.8e\n', ...
  norm(V_a1b1c1w*cu_a1b1c1w-u(X,Y))/norm(u(X,Y)))

%% let's try adding boundary conditions to unweighted laplacian

% total number of equations
M = n*(n+1)/2;
% laplacian in unweighted basis
Lap_abc_a2b2c2 = K_a2b1c2_a2b2c2*K_a2bc2_a2b1c2*Dx_a1bc1_a2bc2*Dx_abc_a1bc1 + ...\
                 K_a1b2c2_a2b2c2*K_ab2c2_a1b2c2*Dy_ab1c1_ab2c2*Dy_abc_ab1c1;  

               %%
%%%%%  Find quad nodes that are on the boundary of the triangle.
Nsamp = 10; NN = 0:(Nsamp-1);
xtmp = (-cos(pi*NN/(Nsamp-1))'+1)/2;              % 2nd-kind Chebyshev points.
% use gjQuad
xtmp = xtmp(2:end-1);

%%%
xx_bottom = [xtmp';0*xtmp'];
xx_left = [0*xtmp'; xtmp'];
xx_hyp = [xtmp'; 1-xtmp'];


% eval polynomials on the boundary
Vx0 = jPoly_tri(xx_left(1,:).',xx_left(2,:).',H_abc,n-1,a,b,c);
Vy0 = jPoly_tri(xx_bottom(1,:).',xx_bottom(2,:).',H_abc,n-1,a,b,c);
V1mx = jPoly_tri(xx_hyp(1,:).',xx_hyp(2,:).',H_abc,n-1,a,b,c);
% append horizontally and find how many equations are missing from laplacian
%%%%V_bnd_h = [Vx0 Vy0 V1mx];

% to fill the matrix I just choose to eliminate one of the equations from
% above.
MM = [Vx0; Vy0; V1mx];
MM = MM(2:end,:);

rr = rank(Lap_abc_a2b2c2);

%%%% this is how many rows are zero in laplacian (rank(L)+nbnd_eq = size(L,1))
%%%nbnd_eq = rank(V_bnd_h); %
% BUT the size is wrong - we need it to be (M x 3*length(X))
% so append vertically and transpose
%%%V_bnd_v = [Vx0; Vy0; V1mx]';
% NOW (rank(L)+nbnd_eq_v > size(L,1))
%%%nbnd_eq_v = rank(V_bnd_v);
% get QR pivots
%%%[~,~,Jv] = qr(V_bnd_v,0); Jvs = Jv(1:nbnd_eq);
% slice out the needed rows 
% (even though nbnd_eq < nbnd_eq_v)
%%%V_bnd_vs = V_bnd_v(:,Jvs)';
% replace deficient rows in laplacian
%%%%Lap_abc_a2b2c2(M-nbnd_eq+1:end,:) = V_bnd_vs;
Lap_abc_a2b2c2(rr+1:end,:) = MM;

%% solve with unweighted laplacian
% coefficients of RHS for poisson in (a,b,c)
cf_abc = V_abc'*(f(X,Y).*W);
% promote coeffs to (a+2,b+2,c+2)
cf_a2b2c2 = K_a2b2c1_a2b2c2*(K_a2b1c1_a2b2c1*(K_a1b1c1_a2b1c1*(K_a1b1c_a1b1c1* ...\
            (K_a1bc_a1b1c*(K_abc_a1bc*cf_abc)))));
% solution coefficients
cu_abc = Lap_abc_a2b2c2 \ ...
        [cf_a2b2c2(1:rr);zeros(M-rr,1)];
% 2-norm relative error
fprintf('unweighted laplacian error: \t %1.8e\n', ...
        (norm(V_abc*cu_abc-u(X,Y))/norm(u(X,Y))))
