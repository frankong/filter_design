function x = fdr(N, a, b, v)
% Design a zero-phase order-N filter
% that minimizes ripples, ie |X(exp(jw)) - v|
%
% Inputs:
%       N - Filter order
%       a - Frequency bands start point
%       b - Frequency bands end point
%       v - value at frequency bands
% Output:
%       x - Order-N filter
%
% Example:
%       Low pass filter with cutoff pi/2 and 1.0% passband ripple
%       x = filt_design(20, -pi/2, pi/2, 0.99, 1.01);

assert(length(a) == length(b) & length(b) == length(v));
M = length(a);

% Construct matrices
A = amat(N);
Bs = cell(M,1);
for m = 1:M
    Bs{m} = bmat(N, a(m), b(m));
end
delta = sparse(N+1, 1, 1.0, 2*N+1, 1);

% Solve
cvx_begin SDP
    variable X(N+1, N+1)        hermitian
    variable Fl(N+1, N+1, M)    hermitian
    variable Gl(N, N, M)        hermitian
    variable Fu(N+1, N+1, M)    hermitian
    variable Gu(N, N, M)        hermitian
    variable r
    minimize r
    subject to
        for m = 1:M
            A*vec(X) - (v(m)-r)*delta == A*vec(Fl(:,:,m)) + Bs{m}*vec(Gl(:,:,m))
            (v(m)+r)*delta - A*vec(X) == A*vec(Fu(:,:,m)) + Bs{m}*vec(Gu(:,:,m))
            Fl(:,:,m) >= 0
            Gl(:,:,m) >= 0
            Fu(:,:,m) >= 0
            Gu(:,:,m) >= 0
        end
        X >= 0
cvx_end

x = A*vec(X);

end

function A = amat(N)
si = zeros((N+1)^2,1);
sj = zeros((N+1)^2,1);
ss = zeros((N+1)^2,1);

c = 1;
for i = 0:N
    for j = 0:N
        si(c) = i - j + N + 1;
        sj(c) = (i+1) + j * (N+1);
        ss(c) = 1.0;
        c = c + 1;
    end
end
A = sparse(si, sj, ss, 2*N+1, (N+1)^2);
end


function B = bmat(N, alpha, beta)

a = tan(alpha/2);
b = tan(beta/2);

if (pi == beta && -pi == alpha)
    d0 = 0;
    d1 = 0;
elseif (pi == beta)
    d0 = -a/2;
    d1 = -a/4 - 1j/4;
elseif (-pi == alpha)
    d0 = b/2;
    d1 = b/4 + 1j/4;
else
    d0 = -(a*b+1)/2;
    d1 = (1-a*b)/4 - 1j*(a+b)/4;
end

si = zeros(3*N^2,1);
sj = zeros(3*N^2,1);
ss = zeros(3*N^2,1);

c = 1;
for i = 0:N-1
    for j = 0:N-1
        si(c) = i - j + N + 1;
        sj(c) = (i+1) + j * N;
        ss(c) = d0;
        c = c + 1;
        
        
        si(c) = i - j + N + 1 - 1;
        sj(c) = (i+1) + j * N;
        ss(c) = d1;
        c = c + 1;
        
        
        si(c) = i - j + N + 1 + 1;
        sj(c) = (i+1) + j * N;
        ss(c) = d1';
        c = c + 1;
    end
end
B = sparse(si, sj, ss, 2*N+1, N^2);
end