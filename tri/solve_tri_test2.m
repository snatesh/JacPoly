clear all; close all; clc;
set(groot, 'defaultLineLineWidth', 2);
set(groot,'defaultLineMarkerSize',12);
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultAxesFontSize',25);
set(groot,'defaultAxesTitleFontSizeMultiplier',1.1);
set(groot,'defaultLegendFontSize',25);
%% Specify parameters and problem
% Koornwinder poly params (unit weight)
a = 1/2; b = 1/2; c = 1/2;
% legendre analog quadrature rule on triangle
load 'triquadLeg_16_27.mat';
R = Zk(1:N); S = Zk(N+1:2*N); W = Zk(2*N+1:3*N);
% specify problem
syms x y
U(x,y) = x.*y.*(1-x-y).*(x+y).^8;
Us{1} = U;
U(x,y) = x.*y.*(1-x-y).*exp(y.*sin(x+y));
Us{2} = U;
U(x,y) = x.*y.*(1-x-y).*exp(-(x.^2+y.^2));
Us{3} = U;
% vertices of triangle
tri_ref = [0 1 0 0 0 1];
tri = [0.5 0.7 0.6 0 0.05 0.5];
%tri = [0 1 0 0 0 0.5];
Atri = [tri(2)-tri(1) tri(3)-tri(1);tri(5)-tri(4) tri(6)-tri(4)];
[X_fun,Y_fun,L_fun] = get_Param();

% get general tri coords
Xrs = X_fun(R,S,tri(1),tri(2),tri(3));
Yrs = Y_fun(R,S,tri(4),tri(5),tri(6));
%%
ms = 3:15;
for iF = 1:length(Us)
disp(iF);
U(x,y) = Us{iF};
dxU(x,y) = diff(U,'x',1);
dyU(x,y) = diff(U,'y',1);
lapU(x,y) = diff(dxU,'x',1)+diff(dyU,'y',1);
dyF(x,y) = diff(lapU,'y',1);
% convert to matlab funcs
f = matlabFunction(lapU);
u = matlabFunction(U);
dxyF(x,y) = diff(dyF,'x',1);
dyf = matlabFunction(dyF);
dxyf = matlabFunction(dxyF);
%% Precompute
for j = 1:length(ms)
% highest poly degree is m-1, and there are M total polys
m = ms(j); n = m+1; M = m*(m+1)/2; 
%  Find quad nodes that are on the boundary of the triangle.
if j == length(ms)
rr = m*(m-1)/2; 
nbnd_eq = 2*m+1;%n*(n+1)/2 - rr; 
Nsamp_lr = floor(nbnd_eq/(2+sqrt(2))); 
Nsamp_h = nbnd_eq-2*(Nsamp_lr);
xtmp_lr = (1+gjQuad(Nsamp_lr,a-1/2,c-1/2))/2;
xtmp_h = (1+gjQuad(Nsamp_h,a-1/2,c-1/2))/2;
% for nbnd_eq = 7, put 2 on bottom and left, 3 on hyp
xx_bottom_ref = [xtmp_lr';0*xtmp_lr'];
xx_left_ref = [0*xtmp_lr'; xtmp_lr'];
% TRYME : put higher order on hyp
xx_hyp_ref = [xtmp_h'; 1-xtmp_h'];

xx_bottom = Atri*xx_bottom_ref + [tri(1);tri(4)];
xx_left = Atri*xx_left_ref + [tri(1);tri(4)];
xx_hyp = Atri*xx_hyp_ref + [tri(1);tri(4)];

figure(3)
plot(xx_bottom_ref(1,:),xx_bottom_ref(2,:),'k-'); hold on;
plot(xx_left_ref(1,:),xx_left_ref(2,:),'k-')
plot(xx_hyp_ref(1,:),xx_hyp_ref(2,:),'k-')



Xrs = X_fun(R,S,tri(1),tri(2),tri(3));
Yrs = Y_fun(R,S,tri(4),tri(5),tri(6));
plot(Xrs,Yrs,'.');
plot(xx_bottom(1,:),xx_bottom(2,:),'r-');
plot(xx_left(1,:),xx_left(2,:),'r-')
plot(xx_hyp(1,:),xx_hyp(2,:),'r-')
end

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
H_a1b1c2 = structure_factors_tri(n+1,a+1,b+1,c+2);
H_a1b2c2 = structure_factors_tri(n+1,a+1,b+2,c+2);
H_a2b2c2 = structure_factors_tri(n+1,a+2,b+2,c+2);
H_a2b1c1 = structure_factors_tri(n+1,a+2,b+1,c+1);
H_a2b2c1 = structure_factors_tri(n+1,a+2,b+2,c+1);

% vandermonde under (a,b,c) 
V_abc = jPoly_tri(R,S,H_abc,n-1,a,b,c);
% weighted vandermonde under (a+1,b+1,c+1)
V_a1b1c1w = jPoly_tri_weighted(R,S,H_a1b1c1(1:n,1:n),n-2,a+1,b+1,c+1);

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
K_a1b1c2_a2b1c2 = promotion_mat_tri(a+1,b+1,c+2,H_a1b1c2,H_a2b1c2,0);


K_a2b1c1_a2b2c1 = promotion_mat_tri(a+2,b+1,c+1,H_a2b1c1,H_a2b2c1,1);
K_a2b2c1_a2b2c2 = promotion_mat_tri(a+2,b+2,c+1,H_a2b2c1,H_a2b2c2,2);

% unweighted derivative matrices
Dx_abc_a1bc1 = D1_tri(a,b,c,H_abc,H_a1bc1,0);
Dx_a1bc1_a2bc2 = D1_tri(a+1,b,c+1,H_a1bc1,H_a2bc2,0);
Dx_ab1c1_a1b1c2 = D1_tri(a,b+1,c+1,H_ab1c1,H_a1b1c2,0);

Dy_abc_ab1c1 = D1_tri(a,b,c,H_abc,H_ab1c1,1);
Dy_ab1c1_ab2c2 = D1_tri(a,b+1,c+1,H_ab1c1,H_ab2c2,1);

% weighted derivative amatrices
Wx_a1b1c1_ab1c = D1_tri_weighted(a+1,b+1,c+1,H_a1b1c1,H_ab1c,0);
Wy_a1b1c1_a1bc = D1_tri_weighted(a+1,b+1,c+1,H_a1b1c1,H_a1bc,1);

% lowering matrices
Lx_ab1c_abc = lowering_mat_tri(a,b+1,c,H_ab1c,H_abc,1);
Ly_a1bc_abc = lowering_mat_tri(a+1,b,c,H_a1bc,H_abc,0);

% laplacian in weighted basis (3/2,/3/2,3/2)
Dxx_w = K_a1bc1_a1b1c1*Dx_abc_a1bc1*Lx_ab1c_abc*Wx_a1b1c1_ab1c; 
Dyy_w = K_ab1c1_a1b1c1*Dy_abc_ab1c1*Ly_a1bc_abc*Wy_a1b1c1_a1bc;
% remove zero rows and columns of weighted laplacian
Lap_a1b1c1w_a1b1c1 = Dxx_w(1:M,1:M) + Dyy_w(1:M,1:M);
conds_w(j) = cond(Lap_a1b1c1w_a1b1c1);

%% solve with weighted laplacian
% coefficients of RHS for poisson in (1/2,1/2,1/2)
cf_abc = V_abc(:,1:M)'*(f(R,S).*W);
% promote coefficients to (3/2,3/2,3/2)
cf_a1b1c1 = K_a1b1c_a1b1c1(1:M,1:M)*(K_a1bc_a1b1c(1:M,1:M)*(K_abc_a1bc(1:M,1:M)*cf_abc));
% solution coefficients
cu_a1b1c1w = Lap_a1b1c1w_a1b1c1\cf_a1b1c1;
% 2-norm relative error
errs_w(j,iF) = norm(V_a1b1c1w*cu_a1b1c1w-u(R,S))/norm(u(R,S));

% laplacian in unweighted basis
Dxx = K_a2b1c2_a2b2c2*K_a2bc2_a2b1c2*Dx_a1bc1_a2bc2*Dx_abc_a1bc1;
Dyy = K_a1b2c2_a2b2c2*K_ab2c2_a1b2c2*Dy_ab1c1_ab2c2*Dy_abc_ab1c1;  
Dxy = K_a2b1c2_a2b2c2*K_a1b1c2_a2b1c2*Dx_ab1c1_a1b1c2*Dy_abc_ab1c1;

% truncate to match order
Dxx = Dxx(1:M,1:M);
Dyy = Dyy(1:M,1:M);
Dxy = Dxy(1:M,1:M);
Lap_abc_a2b2c2 = Dxx + Dyy;
Lap_abc_a2b2c2_new = L_fun(tri(1),tri(2),tri(3),...
                           tri(4),tri(5),tri(6),Dxx,Dyy,Dxy);

%  Find quad nodes that are on the boundary of the triangle.
rr = rank(Lap_abc_a2b2c2);%m*(m-1)/2; %rank(Lap_abc_a2b2c2); 
nbnd_eq = 2*m-1;%M-rr;%2*m+1;%M - rr; 
Nsamp_lr = floor(nbnd_eq/(2+sqrt(2))); 
Nsamp_h = nbnd_eq-2*(Nsamp_lr);

xtmp_lr = (1+gjQuad(Nsamp_lr,a-1/2,c-1/2))/2;
xtmp_h = (1+gjQuad(Nsamp_h,a-1/2,c-1/2))/2;
xx_bottom_ref = [xtmp_lr';0*xtmp_lr'];
xx_left_ref = [0*xtmp_lr'; xtmp_lr'];
xx_hyp_ref = [xtmp_h'; 1-xtmp_h'];
xx_bottom = Atri*xx_bottom_ref + [tri(1);tri(4)];
xx_left = Atri*xx_left_ref + [tri(1);tri(4)];
xx_hyp = Atri*xx_hyp_ref + [tri(1);tri(4)];

% eval polynomials on the boundary
Vx0 = jPoly_tri(xx_left_ref(1,:).',xx_left_ref(2,:).',H_abc,n-2,a,b,c);
Vy0 = jPoly_tri(xx_bottom_ref(1,:).',xx_bottom_ref(2,:).',H_abc,n-2,a,b,c);
V1mx = jPoly_tri(xx_hyp_ref(1,:).',xx_hyp_ref(2,:).',H_abc,n-2,a,b,c);

MM = [Vx0; Vy0; V1mx];
Lap_abc_a2b2c2(rr+1:end,:) = MM;
Lap_abc_a2b2c2_new(rr+1:end,:) = MM;
conds_uw(j) = cond(Lap_abc_a2b2c2_new);
%% solve with unweighted laplacian
% coefficients of RHS for poisson in (a,b,c)
cf_abc = V_abc(:,1:M)'*(f(R,S).*W);
% promote coeffs to (a+2,b+2,c+2)
cf_a2b2c2 = K_a2b2c1_a2b2c2(1:M,1:M)*(K_a2b1c1_a2b2c1(1:M,1:M)*(K_a1b1c1_a2b1c1(1:M,1:M)*(K_a1b1c_a1b1c1(1:M,1:M)* ...\
            (K_a1bc_a1b1c(1:M,1:M)*(K_abc_a1bc(1:M,1:M)*cf_abc)))));
% solution coefficients
cu_abc = Lap_abc_a2b2c2 \ ...
        [cf_a2b2c2(1:rr);zeros(nbnd_eq,1)];
% 2-norm relative error
errs_uw(j,iF) = norm(V_abc(:,1:M)*cu_abc-u(R,S))/norm(u(R,S));
%fprintf('unweighted laplacian error: \t %1.8e\n', errs_uw(j))
Ns(j) = M;

% get coeffs of f(x(r,s),y(r,s))= h(r,s)
ch_abc = V_abc(:,1:M)'*(f(Xrs,Yrs).*W);
% promote coefs
ch_a2b2c2 = K_a2b2c1_a2b2c2(1:M,1:M)*(K_a2b1c1_a2b2c1(1:M,1:M)*(K_a1b1c1_a2b1c1(1:M,1:M)*(K_a1b1c_a1b1c1(1:M,1:M)* ...\
            (K_a1bc_a1b1c(1:M,1:M)*(K_abc_a1bc(1:M,1:M)*ch_abc)))));
% sol coefs
cv_abc = Lap_abc_a2b2c2_new \ ...
             [ch_a2b2c2(1:rr);u(xx_left(1,:)',xx_left(2,:)');...
              u(xx_bottom(1,:)',xx_bottom(2,:)');...
              u(xx_hyp(1,:)',xx_hyp(2,:)')];
         
errs_uw_new(j,iF) = norm(V_abc(:,1:M)*cv_abc-u(Xrs,Yrs))/norm(u(Xrs,Yrs));

%norm(V_abc(:,1:M)*cu_abc_tri-u(Xtri,Ytri))

cu_abc = V_abc(:,1:M)'*(u(Xrs,Yrs).*W);
clap = Lap_abc_a2b2c2_new*cu_abc;
V_a2b2c2 = jPoly_tri(R,S,H_a2b2c2,n-1,a+2,b+2,c+2);
errs_new(j,iF) = norm(V_a2b2c2(:,1:M)*clap-f(Xrs,Yrs));


% if iF == 2
% cdxyf_a2b2c2 = Dxy*cf_abc;
% end


end
end

%%
figure(1); clf
subplot(1,2,1)
spy(Lap_a1b1c1w_a1b1c1);
title('Weighted Laplacian');
subplot(1,2,2)
spy(Lap_abc_a2b2c2);
title('Unweighted Laplacian');
figure(2); clf
subplot(1,2,1)
semilogy(Ns,errs_w(:,1),'o--','displayname','$xy(1-x-y)(x+y)^8$'); hold on;
semilogy(Ns,errs_w(:,2),'o--','displayname','$xy(1-x-y)\exp(y\sin(x+y))$');
semilogy(Ns,errs_w(:,3),'o--','displayname','$xy(1-x-y)\exp(-(x^2+y^2))$');
title('Weighted Laplacian with zero Dirichlet BCs')
xlabs = {'(6,2)','(10,3)','(15,4)','(21,5)','(28,6)','(36,7)','(45,8)','(55,9)','(66,10)','(78,11)','(91,12)','(105,13)','(120,14)'};
xticks(Ns);
xticklabels(xlabs);
ax = gca;
ax.XAxis.FontSize = 11;
xlabel('$(N,n)$');
ax.XLabel.FontSize = 25;
ylabel('Releative error $\frac{||u-\hat{u}||}{||u||}$');
legend show; legend boxoff;
subplot(1,2,2)
semilogy(Ns,errs_uw_new(:,1),'o--'); hold on;
semilogy(Ns,errs_uw_new(:,2),'o--');
semilogy(Ns,errs_uw_new(:,3),'o--');
title('Unweighted Laplacian with zero Dirichlet BCs')
xticks(Ns);
xticklabels(xlabs);
ax = gca;
ax.XAxis.FontSize = 11;
xlabel('$(N,n)$');
ax.XLabel.FontSize = 25;


function [X_fun,Y_fun,L_fun] = get_Param()
% (x,y) in T, (r,s) in T_ref
% [x0,y0,x1,y1,x2,y2] are vertices of T counter clockwise
syms x y x0 x1 x2 y0 y1 y2 A r s;

% map from T to Tref
A = [x1-x0,x2-x0;y1-y0,y2-y0]; 
T2Tref = inv(A) * [x-x0;y-y0];
Tref2T = A*[r;s]+[x0;y0];

R_eq(x,y,x0,x1,x2,y0,y1,y2) = T2Tref(1);
S_eq(x,y,x0,x1,x2,y0,y1,y2) = T2Tref(2);
X_eq(r,s,x0,x1,x2) = Tref2T(1);
Y_eq(r,s,y0,y1,y2) = Tref2T(2);
dxR  = diff(R_eq,x);
dyR  = diff(R_eq,y);
dxS  = diff(S_eq,x);
dyS  = diff(S_eq,y);
% % diff ops on T_ref
syms Drr_eq Dss_eq Drs_eq; 
syms Drr_eq Dss_eq Drs_eq; 
% % get laplace operator on T
Dxx = dxR^2 .* Drr_eq + 2 * dxR .* dxS .* Drs_eq + dxS^2 .* Dss_eq;
Dyy = dyR^2 .* Drr_eq + 2 * dyR .* dyS .* Drs_eq + dyS^2 .* Dss_eq;
L_eq(x0,x1,x2,y0,y1,y2,Drr_eq,Dss_eq,Drs_eq) = Dxx+Dyy;

X_fun = matlabFunction(X_eq);
Y_fun = matlabFunction(Y_eq);
L_fun = matlabFunction(L_eq);
end






