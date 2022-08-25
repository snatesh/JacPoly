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
syms x y t
U(x,y,t) = x.*y.*(t*(1-x)-y).*(x+y).^8;
lapU(x,y,t) = diff(U,'x',2)+diff(U,'y',2);
Us{1} = matlabFunction(U);
lapUs{1} = matlabFunction(lapU);
U(x,y,t) = x.*y.*(t*(1-x)-y).*exp(y.*sin(x+y));
lapU(x,y,t) = diff(U,'x',2)+diff(U,'y',2);
Us{2} = matlabFunction(U);
lapUs{2} = matlabFunction(lapU);
U(x,y,t) = x.*y.*(t*(1-x)-y).*exp(-(x.^2+y.^2));
lapU(x,y,t) = diff(U,'x',2)+diff(U,'y',2);
Us{3} = matlabFunction(U);
lapUs{3} = matlabFunction(lapU);

% vertices of triangle
tri_ref = [0 1 0 0 0 1];
tops = [linspace(0.001,0.09,80),linspace(0.1,1,80)];

[X_fun,Y_fun,L_fun] = get_Param();

ms = 3:14;
% Precompute
for j = 1:length(ms)
    % highest poly degree is m-1, and there are M total polys
    m = ms(j); n = m+1; M = m*(m+1)/2;    
    % normalization under (a,b,c), (a+1,b,c) etc.
    H_abc = structure_factors_tri(n+1,a,b,c);
    H_a1bc = structure_factors_tri(n+1,a+1,b,c);
    H_a1b1c = structure_factors_tri(n+1,a+1,b+1,c);
    H_a1b1c1 = structure_factors_tri(n+1,a+1,b+1,c+1);
    H_ab1c = structure_factors_tri(n+1,a,b+1,c);
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
    K_abc_a1bc = promotion_mat_tri(a,b,c,H_abc,H_a1bc,0);
    K_a1bc_a1b1c = promotion_mat_tri(a+1,b,c,H_a1bc,H_a1b1c,1);
    K_a1b1c_a1b1c1 = promotion_mat_tri(a+1,b+1,c,H_a1b1c,H_a1b1c1,2);
    K_a1bc1_a1b1c1 = promotion_mat_tri(a+1,b,c+1,H_a1bc1,H_a1b1c1,1);
    K_ab1c1_a1b1c1 = promotion_mat_tri(a,b+1,c+1,H_ab1c1,H_a1b1c1,0);
    % unweighted derivative matrices
    Dx_abc_a1bc1 = D1_tri(a,b,c,H_abc,H_a1bc1,0);
    Dy_abc_ab1c1 = D1_tri(a,b,c,H_abc,H_ab1c1,1);
    % weighted derivative amatrices
    Wx_a1b1c1_ab1c = D1_tri_weighted(a+1,b+1,c+1,H_a1b1c1,H_ab1c,0);
    Wy_a1b1c1_a1bc = D1_tri_weighted(a+1,b+1,c+1,H_a1b1c1,H_a1bc,1);
    % lowering matrices
    Lx_ab1c_abc = lowering_mat_tri(a,b+1,c,H_ab1c,H_abc,1);
    Ly_a1bc_abc = lowering_mat_tri(a+1,b,c,H_a1bc,H_abc,0);
    % laplacian in weighted basis (3/2,/3/2,3/2)
    Dxx_w = K_a1bc1_a1b1c1*Dx_abc_a1bc1*Lx_ab1c_abc*Wx_a1b1c1_ab1c;
    Dyy_w = K_ab1c1_a1b1c1*Dy_abc_ab1c1*Ly_a1bc_abc*Wy_a1b1c1_a1bc;
    Dxy_w = K_a1bc1_a1b1c1*Dx_abc_a1bc1*Ly_a1bc_abc*Wy_a1b1c1_a1bc;
    Dxx_w = Dxx_w(1:M,1:M);
    Dyy_w = Dyy_w(1:M,1:M);
    Dxy_w = Dxy_w(1:M,1:M);
    
    for iT = 1:length(tops)
        tri = [0 1 0 0 0 tops(iT)];
        Atri = [tri(2)-tri(1) tri(3)-tri(1);tri(5)-tri(4) tri(6)-tri(4)];
        a1 = norm([tri(2),tri(5)]-[tri(1),tri(4)]);
        b1 = norm([tri(2),tri(5)]-[tri(3),tri(6)]);
        c1 = norm([tri(1),tri(4)]-[tri(3),tri(6)]);
        abc = a1 * b1 * c1; s = (a1 + b1 + c1) / 2;
        aratios(iT) = abc/(8 * (s - a1) * (s - b1) * (s - c1));
        % remove zero rows and columns of weighted laplacian
        Lap_a1b1c1w_a1b1c1 = L_fun(tri(1),tri(2),tri(3),...
                                   tri(4),tri(5),tri(6),Dxx_w,Dyy_w,Dxy_w);
        conds_w(iT,j) = cond(Lap_a1b1c1w_a1b1c1);
        
        % get general tri coords
        Xrs = X_fun(R,S,tri(1),tri(2),tri(3));
        Yrs = Y_fun(R,S,tri(4),tri(5),tri(6));

        for iF = 1:length(Us)
            f = @(x,y) lapUs{iF}(x,y,tops(iT));
            u = @(x,y) Us{iF}(x,y,tops(iT));
            % solve with weighted laplacian
            % coefficients of RHS for poisson in (1/2,1/2,1/2)
            cf_abc = V_abc(:,1:M)'*(f(Xrs,Yrs).*W);
            % promote coefficients to (3/2,3/2,3/2)
            cf_a1b1c1 = K_a1b1c_a1b1c1(1:M,1:M)*(K_a1bc_a1b1c(1:M,1:M)*(K_abc_a1bc(1:M,1:M)*cf_abc));
            % solution coefficients
            cu_a1b1c1w = Lap_a1b1c1w_a1b1c1\cf_a1b1c1;
            % 2-norm relative error
            errs_w(iT,j,iF) = norm(V_a1b1c1w*cu_a1b1c1w-u(Xrs,Yrs))/norm(u(Xrs,Yrs));
        end
    end
    disp(j);            
    Ns(j) = M;
end
save('weighted.mat');

%%

figure(1)
clf
for j = 1:length(tops)
pl1=semilogy(Ns,errs_w(j,:,1),'k--','displayname','$(x+y)^8$'); hold on;
pl2=semilogy(Ns,errs_w(j,:,2),'r--','displayname','$\exp(y\sin(x+y))$'); hold on;
pl3=semilogy(Ns,errs_w(j,:,3),'b--','displayname','$\exp(-(x^2+y^2))$'); hold on;
if j > 1
  set(get(get(pl1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
  set(get(get(pl2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
  set(get(get(pl3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end
end
title('Weighted Laplacian with zero Dirichlet BCs');
xlim([0,105])
xlabs = {'(6,2)','(10,3)','(15,4)','(21,5)','(28,6)','(36,7)','(45,8)','(55,9)','(66,10)','(78,11)','(91,12)','(105,13)'};
xticks(Ns);
xticklabels(xlabs);
ax = gca;
ax.XAxis.FontSize = 10;
xlabel('$(N,n)$');
ylabel('Releative error $\frac{||u-\hat{u}||}{||u||}$');
ax.XLabel.FontSize = 25;
legend show
legend boxoff
legend('location','southwest')
axes('Position',[.6 .6 .2 .2])
box on
semilogy(aratios,conds_w(:,12),'.-');
xlabel('aspect ratio','fontsize',14); 
title('Condition number of Laplacian for $n=13$','fontsize',14);

figure(2)
for j = 1:length(tops)
semilogy(Ns,conds_w(j,:,1),'k--'); hold on;
end
xlabel('$N$','fontsize',25); 
title('Condition number of weighted Laplacian','fontsize',25);



% figure(1)
% clf
% pl(1)=semilogy(aratios,errs_w(:,12,1),'o--','displayname','$(x+y)^8$'); hold on;
% pl(2)=semilogy(aratios,errs_w(:,12,2),'s--','displayname','$\exp(y\sin(x+y))$'); hold on;
% pl(3)=semilogy(aratios,errs_w(:,12,3),'^--','displayname','$\exp(-(x^2+y^2))$'); hold on;
% for j=1:3
%     pl(j).MarkerFaceColor = pl(j).Color;
% end
% legend show
% legend boxoff
% legend('location','southeast')
% ax = gca;
% ax.XAxis.FontSize = 14.5;
% xlabel('aspect ratio');
% ax.XLabel.FontSize = 25;
% ylabel('Releative error $\frac{||u-\hat{u}||}{||u||}$');
% title('Order n=13');
% 
% axes('Position',[.65 .45 .2 .2])
% box on
% semilogy(aratios,conds_w(:,12),'.-');
% xlabel('aspect ratio','fontsize',10); title('Condition number of Laplacian','fontsize',10);
% 
% 
% 
% figure(2)
% clf
% pl(1)=semilogy(Ns,errs_w(55,:,1),'o--','displayname','$(x+y)^8$'); hold on;
% pl(2)=semilogy(Ns,errs_w(55,:,2),'s--','displayname','$\exp(y\sin(x+y))$'); hold on;
% pl(3)=semilogy(Ns,errs_w(55,:,3),'^--','displayname','$\exp(-(x^2+y^2))$'); hold on;
% for j=1:3
%     pl(j).MarkerFaceColor = pl(j).Color;
% end
% xlabs = {'(6,2)','(10,3)','(15,4)','(21,5)','(28,6)','(36,7)','(45,8)','(55,9)','(66,10)','(78,11)','(91,12)','(105,13)','(120,14)'};
% xticks(Ns);
% xticklabels(xlabs);
% ax = gca;
% ax.XAxis.FontSize = 10;
% xlabel('$(N,n)$');
% ylabel('Releative error $\frac{||u-\hat{u}||}{||u||}$');
% ax.XLabel.FontSize = 25;
% legend show
% legend boxoff
% legend('location','southwest')
% title('aspect ratio = 3.2');
% 
% axes('Position',[.6 .6 .2 .2])
% box on
% semilogy(Ns,conds_w(55,:),'.-');
% xlabel('$N$','fontsize',10); title('Condition number of Laplacian','fontsize',10);
% 
% 

% %%
% figure(1);
% semilogy(Ns,errs_w(1,:,1),'o--','displayname','$xy(1-x-y)(x+y)^8$'); hold on;
% semilogy(Ns,errs_w(1,:,2),'o--','displayname','$xy(1-x-y)\exp(y\sin(x+y))$');
% semilogy(Ns,errs_w(1,:,3),'o--','displayname','$xy(1-x-y)\exp(-(x^2+y^2))$');
% title('Weighted Laplacian with zero Dirichlet BCs')
% xlabs = {'(6,2)','(10,3)','(15,4)','(21,5)','(28,6)','(36,7)','(45,8)','(55,9)','(66,10)','(78,11)','(91,12)','(105,13)'};
% xticks(Ns);
% xticklabels(xlabs);
% ax = gca;
% ax.XAxis.FontSize = 14.5;
% xlabel('$(N,n)$');
% ax.XLabel.FontSize = 25;
% ylabel('Releative error $\frac{||u-\hat{u}||}{||u||}$');
% legend show; legend boxoff;



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






