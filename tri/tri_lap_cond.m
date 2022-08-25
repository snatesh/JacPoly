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
% vertices of tr0iangle
tri_ref = [0 1 0 0 0 1];
tri = [0 1 0 0 0 0.001];
[X_fun,Y_fun,L_fun] = get_Param();

%% Precompute
% highest poly degree is m-1, and there are M total polys
m = 14; n = m+1; M = m*(m+1)/2;
% normalization under (a,b,c), (a+1,b,c) etc.
H_abc = structure_factors_tri(n+1,a,b,c);
H_a1bc = structure_factors_tri(n+1,a+1,b,c);
H_a1b1c = structure_factors_tri(n+1,a+1,b+1,c);
H_a1b1c1 = structure_factors_tri(n+1,a+1,b+1,c+1);
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

% promotion matrices (a+2,b,c+2) -> (a+2,b+1,c+2) etc.
K_a2bc2_a2b1c2 = promotion_mat_tri(a+2,b,c+2,H_a2bc2,H_a2b1c2,1);
K_a2b1c2_a2b2c2 = promotion_mat_tri(a+2,b+1,c+2,H_a2b1c2,H_a2b2c2,1);
K_ab2c2_a1b2c2 = promotion_mat_tri(a,b+2,c+2,H_ab2c2,H_a1b2c2,0);
K_a1b2c2_a2b2c2 = promotion_mat_tri(a+1,b+2,c+2,H_a1b2c2,H_a2b2c2,0);
K_abc_a1bc = promotion_mat_tri(a,b,c,H_abc,H_a1bc,0);
K_a1bc_a1b1c = promotion_mat_tri(a+1,b,c,H_a1bc,H_a1b1c,1);
K_a1b1c_a1b1c1 = promotion_mat_tri(a+1,b+1,c,H_a1b1c,H_a1b1c1,2);
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

% laplacian in unweighted basis
Dxx = K_a2b1c2_a2b2c2*K_a2bc2_a2b1c2*Dx_a1bc1_a2bc2*Dx_abc_a1bc1;
Dyy = K_a1b2c2_a2b2c2*K_ab2c2_a1b2c2*Dy_ab1c1_ab2c2*Dy_abc_ab1c1;
Dxy = K_a2b1c2_a2b2c2*K_a1b1c2_a2b1c2*Dx_ab1c1_a1b1c2*Dy_abc_ab1c1;

% truncate to match order
Dxx = Dxx(1:M,1:M);
Dyy = Dyy(1:M,1:M);
Dxy = Dxy(1:M,1:M);
% bottom length
a1 = norm([tri(2),tri(5)]-[tri(1),tri(4)]);
% hyp length
b1 = norm([tri(2),tri(5)]-[tri(3),tri(6)]);
% left length
c1 = norm([tri(1),tri(4)]-[tri(3),tri(6)]);
perim = a1+b1+c1;
s = perim/2;

%  Distribute quad nodes on the boundary of the triangle.
rr = (m-2)*(m-1)/2;
nbnd_eq = 2*m-1;

Nsamp_l = ceil(c1/perim*nbnd_eq);
Nsamp_h = floor(b1/perim*nbnd_eq);
Nsamp_b = nbnd_eq-(Nsamp_l+Nsamp_h);

Lap_abc_a2b2c2 = L_fun(tri(1),tri(2),tri(3),...
    tri(4),tri(5),tri(6),Dxx,Dyy,Dxy);

as = linspace(0.3,0.4,100); bs = as;
mincond = inf;
for ia = 1:length(as)
    for ib = 1:length(bs)
        xtmp_b = (1+gjQuad(Nsamp_b,as(ia),bs(ib)))/2;
        xtmp_l = (1+gjQuad(Nsamp_l,as(ia),bs(ib)))/2;
        xtmp_h = (1+gjQuad(Nsamp_h,as(ia),bs(ib)))/2;
        
        
        xx_bottom_ref = [xtmp_b';0*xtmp_b'];
        xx_left_ref = [0*xtmp_l'; xtmp_l'];
        xx_hyp_ref = [xtmp_h'; 1-xtmp_h'];
        
        % eval polynomials on the boundary
        Vx0 = jPoly_tri(xx_left_ref(1,:).',xx_left_ref(2,:).',H_abc,n-2,a,b,c);
        Vy0 = jPoly_tri(xx_bottom_ref(1,:).',xx_bottom_ref(2,:).',H_abc,n-2,a,b,c);
        V1mx = jPoly_tri(xx_hyp_ref(1,:).',xx_hyp_ref(2,:).',H_abc,n-2,a,b,c);
        
        MM = [Vx0; Vy0; V1mx];
        Lap_abc_a2b2c2(rr+1:end,:) = MM;
        conds_uw = cond(Lap_abc_a2b2c2);
        if conds_uw < mincond
            mincond = conds_uw
            astar=ia; bstar=ib;
        end
    end
end

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






