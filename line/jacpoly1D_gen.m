clear all; close all; clc;
set(groot, 'defaultLineLineWidth', 2);
set(groot,'defaultLineMarkerSize',12);
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultAxesFontSize',25);
set(groot,'defaultAxesTitleFontSizeMultiplier',1.1);
set(groot,'defaultLegendFontSize',25);

% define problem p(x)u'' + r(x)u' + q(x)u = f(x)
% and solution u(x) (for manufacturing rhs and bc)
syms x
u(x) = exp(x .* sin(x));
up(x) = diff(u); 
upp(x) = diff(up);
P(x) = sin(x);
R(x) = cos(x);
Q(x) = atan(x);
F(x) = P(x) * upp(x) + R(x) * up(x) + Q(x) * u(x);
% convert to matlab funcs
f = matlabFunction(F);
p = matlabFunction(P);
r = matlabFunction(R);
q = matlabFunction(Q);
sol = matlabFunction(u);
clear x u up upp P R Q F
%%
% Jacoby polymomial parameters
a = 3; b = 4; mp = 10; mr = 10; mq = 10;
% create dense grid for visualization and analysis
[xs,ws] = gjQuad(40,a,b);
Ns = 10:1:40; 

[cu,err,Vs,L_hat] = sol_verify(a,b,xs,ws,f,p,r,q,sol,Ns,mp,mr,mq);

figure(1);
semilogy(Ns,err,'--'); hold on;
xlabel('$N$');
ylabel('$||\hat{u}-u||_{L_2^w}$');
figure(2)
pl=plot(xs,Vs*cu,'-'); hold on;
xlabel('$x$');
ylabel('$u(x)$');
%%
function [cu,err,Vs,L_hat] = sol_verify(a,b,xs,ws,f,p,r,q,sol,Ns,mp,mr,mq)
err = zeros(length(Ns),1); 
for j = 1:length(Ns)
N = Ns(j); 
[cu,L_hat] = solve(a,b,N,mp,mr,mq,f,p,r,q,[sol(-1),sol(1)]);
% interpolation op for sol coeffs onto resolved grid
[Vs,~] = jPoly(xs,N,a,b);
% weighted l2 error
err(j) = sqrt(ws*(Vs*cu-sol(xs)).^2);
figure(3);
spy(L_hat);
end
end

function [cu,varargout] = solve(a,b,n,mp,mr,mq,f,p,r,q,bc)
% solve pu''+ru' + qu= f with BCs

% Inputs: 
% a,b - Jacoby poly params
% n-1 - degree of f
% m[prq]-1 - degere of p,r,q 
% f   - rhs
% p,r,q - as in ODE
% bc  - [bc(1),bc(2)] boundary condition at -1,1

% Outputs:
% cu - coefficients of sol

% f lives on n point grid
[Xf,Wf] = gjQuad(n,a,b);
[Vf,hf] = jPoly(Xf,n,a,b);
% p lives on mp pt grid 
[Xp,Wp] = gjQuad(mp,a,b); 
[Vp,hp] = jPoly(Xp,mp,a,b);
% r lives on mr pt grid 
[Xr,Wr] = gjQuad(mr,a,b);
[Vr,hr] = jPoly(Xr,mr,a,b);
% q lives on mq pt grid 
[Xq,Wq] = gjQuad(mq,a,b);
[Vq,hq] = jPoly(Xq,mq,a,b);

% promotion operator
Kf = promotion_mat(a,b,n);
K1f = promotion_mat(a+1,b+1,n);
Kp = promotion_mat(a,b,mp);
K1p = promotion_mat(a+1,b+1,mp);
Kr = promotion_mat(a,b,mr);
K1r = promotion_mat(a+1,b+1,mr);
Kq = promotion_mat(a,b,mq);
K1q = promotion_mat(a+1,b+1,mq);

% N coeffs of f
cf = (Vf.' * (f(Xf) .* Wf.')) ./ hf;
% promote f coeffs 
cf = K1f * (Kf * cf); 

% TODO: compute N coeffs of pqr
% and filter based on tol

% mp coeffs of p
cp = (Vp.' * (p(Xp) .* Wp.'))./hp;
% promote p coeffs
cp = K1p * (Kp * cp);
% mr coeffs of r
cr = (Vr.' * (r(Xr) .* Wr.'))./hr;
% promote r coeffs
cr = K1r * (Kr * cr);
% mq coeffs of q
cq = (Vq.' * (q(Xq) .* Wq.'))./hq;
% promote r coeffs
cq = K1q * (Kq * cq);

% multiplication by x operator
Jx = multByX(Xf,n,a+2,b+2); 
% multiplication by "any func" operators
Mop_p = multBy(Jx,mp,a+2,b+2); 
Mop_r = multBy(Jx,mr,a+2,b+2);
Mop_q = multBy(Jx,mq,a+2,b+2);

% multiplication by p,r,q operators
Mp = multByFunc(cp,Mop_p); 
Mr = multByFunc(cr,Mop_r); 
Mq = multByFunc(cq,Mop_q); 

% first deriv op
D1 = sparseFirstD(a,b,n-1);
% second deriv op
D2 = sparseSecD(a,b,n-1);
% complete differential operator
L = Mp * D2 + Mr * (K1f * D1) + Mq * (K1f * Kf); 
% add BCs
L_hat = bcRows(L,a,b,0);
% solve for sol coeffs
cu = L_hat \ [cf(1:end-2);bc(1);bc(2)];
if nargout == 2
   varargout{1} = L_hat; 
end
end

function D1 = sparseFirstD(a,b,N)
d = 0.5 * (a + b + 1 + (1:N));
D1 = sparse(diag(d,1));
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

function Jx = multByX(x,N,a,b)
% Generate the multipilcation operator for 
% coeffs of Jacobi poly expansion

% Input:
% N-1 - max poly deg
% a,b - Jacoby poly parameters
% x - grid from -1 to 1

% Output:
% Multiplication operator M s.t. if f is vector 
% of function values on grid X with f = V*c 
% for Vandermonde V and coeffs c, then
% X.*f = V*(M*c)
% i.e. M*c are the Jacoby coeffs of xf(x)

% Jacobi matrix with p_n contrib
[J, anm1] = jMatOG(N,a,b); 
% vandermonde up to Nth deg poly
[V,~] = jPoly(x,N+1,a,b); 
% p_n contrib for mult mat
en = V(:,end)*anm1; 
Jx = J'; Jx(:,end) = Jx(:,end) + V(:,1:end-1) \ en; 
end

function Mop = multBy(Jx,M,a,b)
% Generate the operator which can be used 
% to create the multiplication by f operator
% for any g represented by at most a degree M-1 poly
%
% Inputs:
% Jx - NxN multiplication by X op (for deg(f)<= N-1)
% M-1 - max poly degree for g (for f*g, with deg(g) <= M-1)
% a,b - Jacobi poly params
%
% Output:
% Mop - Multiplication by "any g" with deg(g) <= M-1
%       i.e. if c_g are M coeffs of g, 
%               c_f are N coeffs of f, 
%               and I is NxN identity,
%       then M = kron(c_g,I)' * Mop is such that 
%            M*g = c_{gf}, the N coeffs of g*f
N = size(Jx,1);
L = clenshawOM(Jx,M,a,b); 
In = eye(N); e0 = (1:M==1)'; 
Mop = L \ kron(e0,In);
end

function Mg = multByFunc(gm,Mop)
N = size(Mop,2);
In = eye(N); 
Mg = kron(gm,In)' * Mop;
end 

function L = clenshawOM(Jx,M,a,b)
% Generate the block lower triangular matrix
% for Clenshaw's algorithm for operator valued input

% Inputs:
% Jx - NxN matrix
% M-1 - max poly degree
% a,b - Jacobi params

% Output:
% L - MN x MN block lower triangular matrix
%     s.t. if c are the M coeffs of f, then
%     f(Jx) = kron(c,I)' * (L \ kron(e0,I))
%     with e0 1st basis vec in R^M 
%     and I the NxN identity    

N = size(Jx,1);
[J,~] = jMatOG(M,a,b); I = speye(N);
L = zeros(N*M,N*M);

cvec = diag(J,-1); 
bvec = diag(J,0);
avec = diag(J,1);
avec = [1;avec]; cvec = [cvec;0];

% fill main block diagonal (NxN blocks)
for j = 1:M
    L((j-1)*N+1:j*N,(j-1)*N+1:j*N) = avec(j) * I;
end

% fill sub1 block diagonal 
for j = 1:(M-1)
    L(j*N+1:(j+1)*N,(j-1)*N+1:j*N) = bvec(j)*I - Jx;
end

% fill sub2 block diagonal
for j = 1:(M-2)
   L((j+1)*N+1:(j+2)*N,(j-1)*N+1:j*N) = cvec(j) * I; 
end
L = sparse(L);
end