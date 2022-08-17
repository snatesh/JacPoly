 function [J,anm1] = jMatOG(N,a,b)
% Make the Jacobi coefficient matrix for 
% orthogonal Jacobi polys
% Input:
% N - max poly deg
% a,b - Jacoby poly parameters
% Output:
% Jacobi matrix J s.t. xP(x) = JP(x) + a_{n-1}p_n(x)e_n,
% where P = (p_0,p_1,...,p_{n-1})
NN = (0:N)';
An = (2*NN+a+b+1).*(2*NN+a+b+2)./(2*(NN+1).*(NN+a+b+1)); An(1) = 0.5*(a+b)+1;
Bn = (a^2-b^2).*(2*NN+a+b+1)./(2*(NN+1).*(NN+a+b+1).*(2*NN+a+b)); Bn(1) = 0.5*(a-b);
Cn = (NN+a).*(NN+b).*(2*NN+a+b+2)./((NN+1).*(NN+a+b+1).*(2*NN+a+b));

avec = 1./An(1:N);
bvec = -Bn(1:N).*avec; 
cvec = Cn(2:end)./An(2:end);
J = spdiags([cvec bvec [0;avec(1:N-1)]], -1:1, N,N); 
anm1 = avec(end);
end

% apb = a + b; aa  = a * a; bb  = b * b;
% avec = zeros(N,1); bvec = avec; cvec = avec;
% A0 = 0.5 * apb + 1; B0 = 0.5 * (a-b); C1 = (1+a) * (1+b) * (2+apb+2) / (2*(1+apb+1)*(2+apb));
% A1 = (2+apb+1)*(2+apb+2)/(4*(1+apb+1));
% avec(1) = 1 / A0; 
% bvec(1) = -B0/A0; 
% cvec(1) = C1/A1;
% for k = 1:(N-1)
%     k2 = 2*k; kp1 = k+1; k21 = 2*kp1;
%     k2apb = k2 + apb; k2apb1 = k21+apb;
%     Ak = (k2apb+1)*(k2apb+2)/(2*(k+1)*(k+apb+1));
%     Ak1 = (k2apb1+1)*(k2apb1+2)/(2*(kp1+1)*(kp1+apb+1));    
%     Bk = (k2apb+1)*(aa-bb)/(2*(k+1)*(k+apb+1)*k2apb);
%     Ck = (kp1+a)*(kp1+b)*(k2apb1+2)/((kp1+1)*(kp1+apb+1)*(k2apb1));
%     cvec(k+1) = Ck/Ak1;
%     bvec(k+1) = -Bk/Ak;
%     avec(k+1) = 1/Ak;
% end
% J = spdiags([cvec bvec [0;avec(1:N-1)]], -1:1, N,N); 
% anm1 = avec(end);