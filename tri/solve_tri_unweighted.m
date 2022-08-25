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
U(x,y) = (x+y).^8;
lapU(x,y) = diff(U,'x',2)+diff(U,'y',2);
lapUs{1} = matlabFunction(lapU);
Us{1} = matlabFunction(U);
U(x,y) = exp(y.*sin(x+y));
lapU(x,y) = diff(U,'x',2)+diff(U,'y',2);
lapUs{2} = matlabFunction(lapU);
Us{2} = matlabFunction(U);
U(x,y) = exp(-(x.^2+y.^2));
lapU(x,y) = diff(U,'x',2)+diff(U,'y',2);
lapUs{3} = matlabFunction(lapU);
Us{3} = matlabFunction(U);
% vertices of triangle
tri_ref = [0 1 0 0 0 1];
%tops = [0.1, 0.125, 0.15, 0.175, 0.2, 0.3, 0.4, 0.5, 0.7, 1];
%tops = [0.1 0.15 0.2 0.4 0.7 1];
tops = [linspace(0.001,0.09,80),linspace(0.1,1,80)];

[X_fun,Y_fun,L_fun] = get_Param();

%% Precompute
ms = 3:14;
for j = 1:length(ms)
    % highest poly degree is m-1, and there are M total polys
    m = ms(j); n = m+1; M = m*(m+1)/2;    
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
    

    for iF = 1:length(Us)
        u = Us{iF};
        f = lapUs{iF};
        for iT = 1:length(tops)
            tri = [0 1 0 0 0 tops(iT)];
            %tri = tris(iT,:);
            Atri = [tri(2)-tri(1) tri(3)-tri(1);tri(5)-tri(4) tri(6)-tri(4)];
            %aratios(iT) = cond(Atri);
            
            % bottom length
            a1 = norm([tri(2),tri(5)]-[tri(1),tri(4)]);
            % hyp length
            b1 = norm([tri(2),tri(5)]-[tri(3),tri(6)]);
            % left length
            c1 = norm([tri(1),tri(4)]-[tri(3),tri(6)]);
            perim = a1+b1+c1; 
            s = perim/2;
            aratios(iT) = a1*b1*c1/(8*(s-a1)*(s-b1)*(s-c1));
            
            %  Distribute quad nodes on the boundary of the triangle.
            rr = (m-2)*(m-1)/2;
            nbnd_eq = 2*m-1;
            
            Nsamp_l = ceil(c1/perim*nbnd_eq);
            Nsamp_h = floor(b1/perim*nbnd_eq);
            Nsamp_b = nbnd_eq-(Nsamp_l+Nsamp_h);
                
             
            %Nsamp_lb = floor(nbnd_eq/(2+sqrt(2)));
            %Nsamp_h = nbnd_eq-2*(Nsamp_lb);
            
            xtmp_b = (1+gjQuad(Nsamp_b,a-1/2,c-1/2))/2;
            xtmp_l = (1+gjQuad(Nsamp_l,a-1/2,c-1/2))/2;
            xtmp_h = (1+gjQuad(Nsamp_h,a-1/2,c-1/2))/2;
%             astar = 1.941414141414141;
%             bstar = 2.012121212121212;
%             xtmp_b = (1+gjQuad(Nsamp_b,astar,bstar))/2;
%             xtmp_l = (1+gjQuad(Nsamp_l,astar,bstar))/2;
%             xtmp_h = (1+gjQuad(Nsamp_h,astar,bstar))/2;

            %xtmp_bl = (1+gjQuad(Nsamp_lb,a-1/2,c-1/2))/2;
            %xtmp_h = (1+gjQuad(Nsamp_h,a-1/2,c-1/2))/2;
            
            xx_bottom_ref = [xtmp_b';0*xtmp_b'];
            xx_left_ref = [0*xtmp_l'; xtmp_l'];
            xx_hyp_ref = [xtmp_h'; 1-xtmp_h'];
            
            
            % eval polynomials on the boundary
            Vx0 = jPoly_tri(xx_left_ref(1,:).',xx_left_ref(2,:).',H_abc,n-2,a,b,c);
            Vy0 = jPoly_tri(xx_bottom_ref(1,:).',xx_bottom_ref(2,:).',H_abc,n-2,a,b,c);
            V1mx = jPoly_tri(xx_hyp_ref(1,:).',xx_hyp_ref(2,:).',H_abc,n-2,a,b,c);
            
            MM = [Vx0; Vy0; V1mx];
            
            
            
            Lap_abc_a2b2c2 = L_fun(tri(1),tri(2),tri(3),...
                                   tri(4),tri(5),tri(6),Dxx,Dyy,Dxy);
            
            Lap_abc_a2b2c2(rr+1:end,:) = MM;
            conds_uw(iT,j) = cond(Lap_abc_a2b2c2);
            % get general tri coords
            Xrs = X_fun(R,S,tri(1),tri(2),tri(3));
            Yrs = Y_fun(R,S,tri(4),tri(5),tri(6));
            xx_bottom = Atri*xx_bottom_ref + [tri(1);tri(4)];
            xx_left = Atri*xx_left_ref + [tri(1);tri(4)];
            xx_hyp = Atri*xx_hyp_ref + [tri(1);tri(4)];
            % get coeffs of f(x(r,s),y(r,s))= h(r,s)
            ch_abc = V_abc(:,1:M)'*(f(Xrs,Yrs).*W);
            % promote coefs
            ch_a2b2c2 = K_a2b2c1_a2b2c2(1:M,1:M)*(K_a2b1c1_a2b2c1(1:M,1:M)*(K_a1b1c1_a2b1c1(1:M,1:M)*(K_a1b1c_a1b1c1(1:M,1:M)* ...\
                (K_a1bc_a1b1c(1:M,1:M)*(K_abc_a1bc(1:M,1:M)*ch_abc)))));
            % sol coefs
            cv_abc = Lap_abc_a2b2c2 \ ...
                [ch_a2b2c2(1:rr);u(xx_left(1,:)',xx_left(2,:)');...
                u(xx_bottom(1,:)',xx_bottom(2,:)');...
                u(xx_hyp(1,:)',xx_hyp(2,:)')];
            
            errs_uw(iT,j,iF) = norm(V_abc(:,1:M)*cv_abc-u(Xrs,Yrs))/norm(u(Xrs,Yrs));
            Ns(j) = M;
        end
    end
    disp(j)
end
save('unweighted.mat');
%%
load weighted.mat 
load unweighted.mat
%%
figure(1)
clf
subplot(1,2,1)
for j = 1:length(tops)

pl1=semilogy(Ns,errs_uw(j,:,1),'k--','displayname','$(x+y)^8$'); hold on;
pl2=semilogy(Ns,errs_uw(j,:,2),'r--','displayname','$\exp(y\sin(x+y))$'); hold on;
pl3=semilogy(Ns,errs_uw(j,:,3),'b--','displayname','$\exp(-(x^2+y^2))$'); hold on;
if j > 1
  set(get(get(pl1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
  set(get(get(pl2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
  set(get(get(pl3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end
end
title('Unweighted Laplacian with non-zero Dirichlet BCs');
xlim([0,105])
xlabs = {'(6,2)','(10,3)','(15,4)','(21,5)','(28,6)','(36,7)','(45,8)','(55,9)','(66,10)','(78,11)','(91,12)','(105,13)'};
xticks(Ns);
xticklabels(xlabs);
ax = gca;
ax.XAxis.FontSize = 12;
xlabel('$(N,n)$');
ylabel('Releative error $\frac{||u-\hat{u}||}{||u||}$');
ax.XLabel.FontSize = 25;

axes('Position',[.28 .7 .14 .14])
box on
semilogy(aratios,conds_uw(:,12),'.-');
xlabel('aspect ratio','fontsize',12); 
title('cond(L) $n=13$','fontsize',12);

subplot(1,2,2)
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
legend show
legend boxoff
legend('location','southwest')
legend('fontsize',20)
xlim([0,105])
xlabs = {'(6,2)','(10,3)','(15,4)','(21,5)','(28,6)','(36,7)','(45,8)','(55,9)','(66,10)','(78,11)','(91,12)','(105,13)'};
xticks(Ns);
xticklabels(xlabs);
ax = gca;
ax.XAxis.FontSize = 12;
xlabel('$(N,n)$');
ax.XLabel.FontSize = 25;
axes('Position',[.745 .7 .14 .14])
box on
semilogy(aratios,conds_w(:,12),'.-');
xlabel('aspect ratio','fontsize',12); 
title('cond(L) $n=13$','fontsize',12);


%%
figure(2)
clf
for j = 1:length(tops)
pl1=semilogy(Ns,conds_uw(j,:,1),'b--','displayname','unweighted'); hold on;
pl2=semilogy(Ns,conds_w(j,:,1),'r--','displayname','weighted');
if j > 1
  set(get(get(pl1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
  set(get(get(pl2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end
end
xlabel('$N$','fontsize',25); 
title('Condition number of Laplacian','fontsize',25);
legend show;
legend boxoff;
legend('location','northwest');
% figure(1)
% clf
% pl(1)=semilogy(aratios,errs_uw(:,12,1),'o--','displayname','$(x+y)^8$'); hold on;
% pl(2)=semilogy(aratios,errs_uw(:,12,2),'s--','displayname','$\exp(y\sin(x+y))$'); hold on;
% pl(3)=semilogy(aratios,errs_uw(:,12,3),'^--','displayname','$\exp(-(x^2+y^2))$'); hold on;
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
% axes('Position',[.6 .4 .2 .2])
% box on
% semilogy(aratios,conds_uw(:,12),'.-');
% xlabel('aspect ratio','fontsize',10); title('Condition number of Laplacian','fontsize',10);
% 
% 
% 
% figure(2)
% clf
% pl(1)=semilogy(Ns,errs_uw(1,:,1),'o--','displayname','$(x+y)^8$'); hold on;
% pl(2)=semilogy(Ns,errs_uw(1,:,2),'s--','displayname','$\exp(y\sin(x+y))$'); hold on;
% pl(3)=semilogy(Ns,errs_uw(1,:,3),'^--','displayname','$\exp(-(x^2+y^2))$'); hold on;
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
% title('aspect ratio = 50');
% 
% axes('Position',[.6 .6 .2 .2])
% box on
% semilogy(Ns,conds_uw(1,:),'.-');
% xlabel('$N$','fontsize',10); title('Condition number of Laplacian','fontsize',10);

% 
% figure(4)
% clf
% for j = 1:length(tops)
% tri = [0 1 0 0 0 tops(j)];
% plot([tri(1),tri(2)], [tri(4),tri(5)], 'k-','markersize',14); hold on;
% plot([tri(2),tri(3)], [tri(5),tri(6)], 'k-','markersize',14);
% plot([tri(3),tri(1)], [tri(6),tri(4)], 'k-','markersize',14);
% end

%pl(2)=semilogy(aratios,errs_uw(:,12,2),'s--','displayname','$\exp(y\sin(x+y))$'); hold on;
%pl(3)=semilogy(aratios,errs_uw(:,12,3),'^--','displayname','$\exp(-(x^2+y^2))$'); hold on;

% 
% figure(2)
% plot(xx_bottom_ref(1,:),xx_bottom_ref(2,:),'k-'); hold on;
% plot(xx_left_ref(1,:),xx_left_ref(2,:),'k-')
% plot(xx_hyp_ref(1,:),xx_hyp_ref(2,:),'k-')
% plot(Xrs,Yrs,'.');
% plot(xx_bottom(1,:),xx_bottom(2,:),'r-');
% plot(xx_left(1,:),xx_left(2,:),'r-')
% plot(xx_hyp(1,:),xx_hyp(2,:),'r-')

% 
% semilogy(Ns,errs_uw(iT,:,2),'o--','displayname','$\exp(y\sin(x+y))$');
% semilogy(Ns,errs_uw(iT,:,3),'o--','displayname','$\exp(-(x^2+y^2))$');
%%


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






