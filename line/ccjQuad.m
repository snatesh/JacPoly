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