clear all; close all; clc;
set(groot, 'defaultLineLineWidth', 2);
set(groot,'defaultLineMarkerSize',12);
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultAxesFontSize',25);
set(groot,'defaultAxesTitleFontSizeMultiplier',1.1);
set(groot,'defaultLegendFontSize',25);
%%
% testing subroutines

% Jacoby polymomial parameters
a = 3; b = 4;
% create dense grid for visualization and analysis
[xs,ws] = gjQuad(50,a,b);
% rhs and sol
fs{1} = @(x) x; sols{1} = @(x) (1/6)*x.*(x.^2-1);
fs{2} = @(x) x.^2; sols{2} = @(x)  (1/12)*(x.^4 - 1);
fs{3} = @(x) x.^3; sols{3} = @(x) (1/20)*x.*(x.^4-1);
fs{4} = @(x) exp(x).*sin(x); sols{4} = @(x) 0.25 * ((-x + exp(1)^2*(x+1)+1)*cos(1)/exp(1) - 2*exp(x).*cos(x));
fs{5} = @(x) atan(x); sols{5} = @(x) 0.5*(x.*(log(2)-log(x.^2+1))+(x.^2-1).*atan(x));
Ns = 3:40; 
legs = {'$x$','$x^{2}$','$x^{3}$','$e^x\sin(x)$','$\tan^{-1}(x)$'};

for j = 1:length(fs)
[cu1,cu2,cu3,err1,err2,err3,Vs] = sol_verify(a,b,xs,ws,fs{j},sols{j},Ns);
leg = legs{j};
figure(1);
subplot(1,3,1)
semilogy(Ns,err1,'--','displayname',leg); hold on;
subplot(1,3,2)
semilogy(Ns,err2,'--','displayname',leg); hold on;
subplot(1,3,3)
semilogy(Ns,err3,'--','displayname',leg); hold on;
figure(2)
pl=plot(xs,Vs*cu1,'-','displayname',leg); hold on;
end
figure(1)
ax1 = subplot(1,3,1);
title('Vandermonde inversion');
ylabel('$||\hat{u}-u||_{L_2^w}$');
xlabel('$N$');
ylim([1e-18,1e-1]);
ax2 = subplot(1,3,2);
title('Weighted $L_2$ projection CCJ');
xlabel('$N$');
ylim([1e-18,1e-1]);
legend show; legend boxoff;
ax3 = subplot(1,3,3);
title('Weighted $L_2$ projection GJ');
xlabel('$N$');
ylim([1e-18,1e-1]);
legend show; legend boxoff;
linkaxes([ax1,ax2,ax3],'y');
figure(2)
xlabel('$x$');
ylabel('$u(x)$');
legend show; legend boxoff;

%% N point Gauss quad rule for jacoby poly with params a,b

N = 10; a = 0; b = 0;
[X,w] = gjQuad(N,a,b)
f = @(x) exp(x).*sin(x);
If = (cos(1) + sin(1) + exp(1)^2 * (-cos(1) + sin(1)))/(2*exp(1));
w*f(X)-If



%%
N = 20;
as = linspace(-0.99,1,10000); bs = as(randperm(length(as)));
as1 = linspace(1.1,20,10000); bs1 = as1(randperm(length(as1)));
as2 = linspace(-0.99,20,10000); bs2 = as2;

f = @(x) exp(sqrt(1-x.^2)).*sin(x);
l2err_new = zeros(length(as),1); ierr_new = l2err_new; gerr_new = l2err_new;
l2err_new1 = zeros(length(as1),1); ierr_new1 = l2err_new1; gerr_new1 = l2err_new1;
l2err_new2 = zeros(length(as2),1); ierr_new2 = l2err_new2; gerr_new2 = l2err_new2;
for k = 1:length(as)
a = as(k); b = bs(k); 
[x,w] = ccjQuad(2*N-1,a,b);
[V,n] = jPoly(x,N,a,b);
c = (V.' * (f(x) .* w.')) ./ n;
l2err_new(k) = sqrt(abs(w*(V*c-f(x)).^2));
[x1,w1] = ccjQuad(N,a,b);
[V1,n] = jPoly(x1,N,a,b);
c = V1 \ f(x1);
ierr_new(k) = sqrt(abs(w*(V*c-f(x)).^2));
[X,W] = gjQuad(N,a,b);
[V2,n] = jPoly(X,N,a,b);
c = (V2.' * (f(X) .* W.'))./n;
gerr_new(k) = sqrt(abs(w*(V*c-f(x))).^2);
end

for k = 1:length(as1)
a = as1(k); b = bs1(k); 
[x,w] = ccjQuad(2*N-1,a,b);
[V,n] = jPoly(x,N,a,b);
c = (V.' * (f(x) .* w.')) ./ n;
l2err_new1(k) = sqrt(abs(w*(V*c-f(x)).^2));
[x1,w1] = ccjQuad(N,a,b);
[V1,n] = jPoly(x1,N,a,b);
c = V1 \ f(x1);
ierr_new1(k) = sqrt(abs(w*(V*c-f(x)).^2));
[X,W] = gjQuad(N,a,b);
[V2,n] = jPoly(X,N,a,b);
c = (V2.' * (f(X) .* W.'))./n;
gerr_new1(k) = sqrt(abs(w*(V*c-f(x))).^2);

end

for k = 1:length(as2)
a = as2(k); b = bs2(k); 
[x,w] = ccjQuad(2*N-1,a,b);
[V,n] = jPoly(x,N,a,b);
c = (V.' * (f(x) .* w.')) ./ n;
l2err_new2(k) = sqrt(abs(w*(V*c-f(x)).^2));
[x1,w1] = ccjQuad(N,a,b);
[V1,n] = jPoly(x1,N,a,b);
c = V1 \ f(x1);
ierr_new2(k) = sqrt(abs(w*(V*c-f(x)).^2));

[X,W] = gjQuad(N,a,b);
[V2,n] = jPoly(X,N,a,b);
c = (V2.' * (f(X) .* W.'))./n;
gerr_new2(k) = sqrt(abs(w*(V*c-f(x))).^2);
end


ax1 = subplot(1,2,1);
semilogy(as+bs,l2err_new,'r.','displayname','$a+b \in (-2,2)$');hold on;
semilogy(as1+bs1,l2err_new1,'b.','displayname','$a,b>1$');
semilogy(as2+bs2,l2err_new2,'k.','displayname','$a=b$','Markersize',7);
title('CCJ');
xlabel('$a+b$');
ylabel('$||\hat{f}-f||_{L_2^w}$, $N=20$');
legend show; legend boxoff;
ylim([1e-20,1e-2]);
ax2 = subplot(1,2,2);
semilogy(as+bs,gerr_new,'r.','displayname','$a+b \in (-2,2)$');hold on;
semilogy(as1+bs1,gerr_new1,'b.','displayname','$a,b>1$');
semilogy(as2+bs2,gerr_new2,'k.','displayname','$a=b$','Markersize',7);
ylim([1e-20,1e-2]);
title('GJ');
xlabel('$a+b$');
linkaxes([ax1,ax2],'y');

%%
function [cu1,cu2,cu3,err1,err2,err3,Vs] = sol_verify(a,b,xs,ws,f,sol,Ns)
err1 = zeros(length(Ns),1); err2 = err1; err3 = err1;
for j = 1:length(Ns)
N = Ns(j);
[cu1,cu2,cu3] = solve(a,b,N,f);
% interpolation op for sol coeffs onto resolved grid
[Vs,~] = jPoly(xs,N,a,b);
% weighted l2 error
err1(j) = sqrt(ws*(Vs*cu1-sol(xs)).^2);
err2(j) = sqrt(ws*(Vs*cu2-sol(xs)).^2);
err3(j) = sqrt(ws*(Vs*cu3-sol(xs)).^2);
end
end

function [cu1,cu2,cu3] = solve(a,b,N,f)
% solve u'' = f with homogeneous Dirichlet BCs
% a,b - Jacoby poly params
% interpolation op for rhs coeffs
[xf,~] = ccjQuad(N,a,b);
[Vf,~] = jPoly(xf,N,a+2,b+2);
% l2 projection op (doubled grid)
%[xl2,wl2] = ccjQuad(2*N-1,a+2,b+2);
%[Vl2,facl2] = jPoly(xl2,N,a+2,b+2);

[xl2,wl2] = ccjQuad(2*N-1,a,b);
[Vl2,facl2] = jPoly(xl2,N,a,b);

% l2 proj with gauss-jacobi
%[xgj,wgj] = gjQuad(N,a+2,b+2);
%[Vgj,facgj] = jPoly(xgj,N,a+2,b+2);

[xgj,wgj] = gjQuad(N,a,b);
[Vgj,facgj] = jPoly(xgj,N,a,b);
% promotion operator
K = promotion_mat(a,b,N);
K1 = promotion_mat(a+1,b+1,N);
K2 = K1 * K;

% rhs coefs using vandermonde inversion
c1 = Vf \ f(xf);
% rhs coefs using l2 projection clencurt-jacobi
c2 = (Vl2.' * (f(xl2) .* wl2.')) ./ facl2;
c2 = K2 * c2;
% rhs coefs using l2 projection gauss-jacobi
c3 = (Vgj.' * (f(xgj) .* wgj.')) ./ facgj;
% promote rhs coefs 
c3 = K2 * c3;
% second deriv op
D2 = sparseSecD(a,b,N-1);
% add BCs
D2_hat = bcRows(D2,a,b,0);
% solve for sol coeffs
cu1 = D2_hat \ [c1(1:end-2);0;0];
cu2 = D2_hat \ [c2(1:end-2);0;0];
cu3 = D2_hat \ [c3(1:end-2);0;0];
end

function D2 = sparseSecD(a,b,N)
% Modal second derivative matrix for Jacoby poly (a,b,N)
d = 0.25 * (a + b + 1 + (2:N)) .* (a + b + 2 + (2:N));
D2 = sparse(diag(d,2));
end

function D2 = bcRows(D2,a,b,type)
% Boundary condition rows for second order 1D eq
% type = 0 - Dirichlet
N = size(D2,1);
if type == 0
    j = 0:N-1;
    jpb = j + b;
    jpa = j + a;
    D2(N-1,:) = (-1).^j .* gamma(jpb+1)./(gamma(jpb-j+1).*gamma(j+1));
    D2(N,:) = gamma(jpa+1)./(gamma(jpa-j+1).*gamma(j+1));
end
end

function [x, w] = ccjQuad(n, a, b)
%CCJQUADWTS   Clenshaw-Curtis-Jacobi quadrature weights.
%   [W, X] = CCJQUAD(N, A, B) returns the N-point Clenshaw-Curtis-Jacobi
%   quadrature nodes, X = CHEBPTS(N), and weights, W, corresponding to the
%   weight function w(t) = (1-t)^A * (1+t)^B on the interval [-1,1].
if ( a == b && a == 0 ) % Clenshaw-Curtis

    c = 2./[1, 1-(2:2:(n-1)).^2];          % Standard Chebyshev moments
    c = [c, c(floor(n/2):-1:2)];           % Mirror for DCT via FFT 
    w = ifft(c);                           % Interior weights
    w([1, n]) = w(1)/2;                    % Boundary weights

elseif ( a == b )       % Gegenbauer
    
    l = a + .5;                            % Gegenbauer parameter
    g0 = gamma(l+.5)*sqrt(pi)/gamma(l+1);
    k = 1:floor((n-1)/2); 
    c = g0*[1, cumprod((k-l-1)./(k+l))];   % Chebyshev moments for (1-x)^a(1+x)^b
    c = [c, c(floor(n/2):-1:2)];           % Mirror for DCT via FFT 
    w = ifft(c);                           % Interior weights
    w([1, n]) = w(1)/2;                    % Boundary weights
    
else                    % Jacobi
    
    c = [1, (a-b)/(a+b+2), zeros(1, n-2)]; % Initialise moments
    for r = 1:n % Recurrence relation for 3F2([r, -r, b +1 ; .5, a+b+2] ; 1 ):
        c(r+2) = - (2*(b-a)*c(r+1) + (a+b+2-r)*c(r)) / (a+b+2+r);
    end
    c = 2^(a+b+1)*gamma(a+1)*gamma(b+1)/gamma(a+b+2) * c; % Moments (with const)
    v = ifft([c(1:n), c(n-1:-1:2)]);       % Mirror for DCT via FFT 
    w = [v(1), 2*v(2:n-1), v(n)];          % Rescale interior weights

end

if ( nargout > 1 )
    x = -cos(pi*(0:n-1)/(n-1))';              % 2nd-kind Chebyshev points.
end
end


%%
% syms a b a1 b1 c1
% syms A [3,3] matrix
% assume(a > -1); assumeAlso(a,'real');
% assume(b > -1); assumeAlso(b,'real');
% A(1,:) = [(1/2)*(1+b)*(2+b), -(1+b), 1];
% A(2,:) = [(1/2)*(1+a)*(2+a), (1+a), 1];
% A(3,:) = [(a+1)*(a+2)/2-(a+2)*(a+b+3)/2+(a+b+3)*(a+b+4)/8, 1/2*(a-b),1];
% coeffs = simplify(symmatrix2sym(A \ [1+b,1+a,0]'));
% N = 100;
% J = jMat(N,2,2);
% X1 = sort(eigs(J,N));
% J = jMat(N,-1/2,-1/2);
% X2 = sort(eigs(J,N));
% cheb = sort(sin((N - 2*(0:N-1) - 1)*pi/(2*N)))';
% plot(X1,sqrt(1-X1.^2),'r.-'); hold on;
% plot(X2,sqrt(1-X2.^2),'bo'); hold on;
% plot(cheb,sqrt(1-cheb.^2),'mp');