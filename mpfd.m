function [x, X] = mpfd(N, a, b, m)
% Design a minimum-phase order-N filter that minimizes maximum absolute deviation from the specified magnitude response
% 
% Inputs:
% 
%       N - Filter order.
%       a - Frequency band start points. From -pi to pi.
%       b - Frequency band end points. From -pi to pi.
%       m - Magnitude response at frequency bands.
% 
%       a, b, and m are length-B vectors where B is the number of bands.
% Output:
% 
%       x - Order-N filter
%       X - Gram matrix of x. Ideally should be x*x'

assert(length(a) == length(b) & length(b) == length(m));
M = length(a);

% Construct matrices
A = amat(N);
Bs = cell(M,1);
for i = 1:M
    Bs{i} = bmat(N, a(i), b(i));
end
dirac = sparse(N+1, 1, 1.0, 2*N+1, 1);

% Solve
cvx_begin SDP
    variable X(N+1, N+1)        semidefinite hermitian     
    variable Fl(N+1, N+1, M)    hermitian
    variable Gl(N, N, M)        hermitian
    variable Fu(N+1, N+1, M)    hermitian
    variable Gu(N, N, M)        hermitian
    variable d
    minimize d - 0.0001*X(1,1)
    subject to
        for i = 1:M
            A*vec(X) - (m(i)^2 - d)*dirac == A*vec(Fl(:,:,i)) + Bs{i}*vec(Gl(:,:,i))
            (m(i)^2 + d)*dirac - A*vec(X) == A*vec(Fu(:,:,i)) + Bs{i}*vec(Gu(:,:,i))
            Fl(:,:,i) >= 0
            Gl(:,:,i) >= 0
            Fu(:,:,i) >= 0
            Gu(:,:,i) >= 0
        end
        d >= 0
cvx_end

[v,s] = eig(X);
[y,i] = max(diag(s));
x = sqrt(s(i,i)) * v(:,i);

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