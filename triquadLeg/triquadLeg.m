addpath('../');
clear all; close all; clc;
set(groot, 'defaultLineLineWidth', 2);
set(groot,'defaultLineMarkerSize',12);
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultAxesFontSize',25);
set(groot,'defaultAxesTitleFontSizeMultiplier',1.1);
set(groot,'defaultLegendFontSize',25);


a = 1/2; b = 1/2; c = 1/2;
load 'triquadLeg_16_27.mat';
Hm = structure_factors_tri(m,a,b,c);
ag = randn(2,1); bg = randn(2,1);
ftest{1} = @(X,Y) cos(2*pi*bg(1) + ag(1)*X + ag(2)*Y);
ftest{2} = @(X,Y) 1./(((1/ag(1))^2+(X-bg(1)).^2).*((1/ag(2)).^2+(Y-bg(2)).^2));
ftest{3} = @(X,Y) exp(-(ag(1)^2*(X-bg(1)).^2+ag(2)^2*(Y-bg(2)).^2));
ftest{4} = @(X,Y) (ag(1)*X+ag(2)*Y).^9;
ftest{5} = @(X,Y) (ag(1)*X+ag(2)*Y).^23;
nms = [2 3; 3 5; 4 6; 5 8; 6 10; 7 12; 8 13; 9 15; 10 17; 11 18; 12 20; 13 22; 14 23; 15 24];

errs = zeros(size(nms,1),size(ftest,2));

for ff = 1:5
f = ftest{ff};
Xk = Zk(1:N); Yk = Zk(N+1:2*N); Wk = Zk(2*N+1:end)';
Vm = jPoly_tri(Xk,Yk,Hm,m-1,a,b,c);
ref = Wk*f(Xk,Yk);


figure(1)
tri = [0 1 0 0 0 1];
for j = 1:(size(nms,1))
n = nms(j,1); m = nms(j,2); 
fname = strcat('triquadLeg_',num2str(n),'_',num2str(m),'.mat');
load(fname);
Xk = Zk(1:N); Yk = Zk(N+1:2*N); Wk = Zk(2*N+1:end)';
Vm = jPoly_tri(Xk,Yk,Hm,m-1,a,b,c);
errs(j,ff) = abs(ref-Wk*f(Xk,Yk))/abs(ref);
Ns(j) = N; Ms(j) = M;
condVm(j) = cond(Vm);
subplot(2,7,j)
title(strcat('$N=',num2str(N),'$'));
plot_tri(tri,'k-'); hold on;
plot(Xk,Yk,'b.');
end
end

errs(find(errs == 0)) = 1e-17;
xlabs = {'(3,2)','(6,4)','(10,5)','(15,7)','(21,9)','(28,11)','(36,12)','(45,14)','(55,16)','(66,17)','(78,19)','(91,21)','(105,22)','(120,23)'};
ms = nms(:,2); ns = nms(:,1);

figure(2)
pl = semilogy(Ns, errs,'o--');
pl(1).DisplayName = '$\cos(2\pi b_1 + \sum_{i=1}^2a_ix_i)$';
pl(2).DisplayName = '$\prod_{i=1}^2(a_i^{-2}+(x_i-b_i)^2)^{-1}$';
pl(3).DisplayName = '$\exp(-\sum_{i=1}^2 a_i(x_i-b_i)^2)$';
pl(4).DisplayName = '$(\sum_{i=1}^2 a_ix_i)^{9}$';
pl(5).DisplayName = '$(\sum_{i=1}^2 a_ix_i)^{23}$';
xticks(Ns);
xticklabels(xlabs);
ax = gca;
ax.XAxis.FontSize = 14.5;
xlabel('$(N,m)$');
ax.XLabel.FontSize = 25;
ylabel('Releative error $\frac{|I-\hat{I}|}{|I|}$');
ylim([1e-17,1]);
legend show;


%%
fprintf('%3g ',ns-1)
fprintf('\n');
fprintf('%3g ',ms-1)
fprintf('\n')
fprintf('%3g ',Ns)
fprintf('\n')
fprintf('%3g ',Ms)
fprintf('\n')
fprintf('%.3g ',condVm)
fprintf('\n')

function plot_tri(tri,col)
plot([tri(1),tri(2)], [tri(4),tri(5)], col); hold on;
plot([tri(2),tri(3)], [tri(5),tri(6)], col);
plot([tri(3),tri(1)], [tri(6),tri(4)], col);
end

