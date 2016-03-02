function x = fdtr(N, a, b, v)
% Design a order-N polynomial
% that minimizes ripples, ie |X(t) - v|
% where X(t) = sum_{i=0}^N x_t t^i
%
% Inputs:
%       N - Polynomial degree. Supports only even N
%       a - Interval start points
%       b - Interval end points
%       v - value at intervals
% Output:
%       x - Order-N polynomial
%

assert(mod(N,2) == 0)
N = N / 2;

assert(length(a) == length(b) & length(b) == length(v));
M = length(a);

% Construct matrices
A = amat(N);
Bs = cell(M,1);
for m = 1:M
    Bs{m} = bmat(N, a(m), b(m));
end
delta = sparse(1, 1, 1.0, 2*N+1, 1);

% Solve
cvx_begin SDP
    variable x(2*N+1)             complex
    variable Fl(N+1, N+1, M)    hermitian
    variable Gl(N, N, M)        hermitian
    variable Fu(N+1, N+1, M)    hermitian
    variable Gu(N, N, M)        hermitian
    variable r
    minimize r
    subject to
        for m = 1:M
            x - (v(m)-r)*delta == A*vec(Fl(:,:,m)) + Bs{m}*vec(Gl(:,:,m))
            (v(m)+r)*delta - x == A*vec(Fu(:,:,m)) + Bs{m}*vec(Gu(:,:,m))
            Fl(:,:,m) >= 0
            Gl(:,:,m) >= 0
            Fu(:,:,m) >= 0
            Gu(:,:,m) >= 0
        end
cvx_end

end

function A = amat(N)
si = zeros((N+1)^2,1);
sj = zeros((N+1)^2,1);
ss = zeros((N+1)^2,1);

c = 1;
for i = 0:N
    for j = 0:N
        si(c) = i + j + 1;
        sj(c) = (i+1) + j * (N+1);
        ss(c) = 1.0;
        c = c + 1;
    end
end
A = sparse(si, sj, ss, 2*N+1, (N+1)^2);
end


function B = bmat(N, a, b)

si = zeros(3*N^2,1);
sj = zeros(3*N^2,1);
ss = zeros(3*N^2,1);

c = 1;
for i = 0:N-1
    for j = 0:N-1
        si(c) = i + j + 1;
        sj(c) = (i+1) + j * N;
        ss(c) = -a*b;
        c = c + 1;
        
        
        si(c) = i + j + 1 + 1;
        sj(c) = (i+1) + j * N;
        ss(c) = a+b;
        c = c + 1;
        
        
        si(c) = i + j + 1 + 2;
        sj(c) = (i+1) + j * N;
        ss(c) = -1;
        c = c + 1;
    end
end
B = sparse(si, sj, ss, 2*N+1, N^2);
end