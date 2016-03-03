function x = pd(N, a, b, m)
% Design a degree-N polynomial that minimizes maximum absolute deviation from the specified magnitude response
% 
% Inputs:
% 
%       N - Polynomial degree.
%       a - Bands start points. From 0 to 1.
%       b - Bands end points. From 0 to 1.
%       v - Magnitude response at intevals.
% 
%       a, b, v are length-B vectors where B is the number of bands.
% Output:
% 
%       x - Order-N filter

assert(mod(N,2) == 0)
N = N / 2;

assert(length(a) == length(b) & length(b) == length(m));
M = length(a);

% Construct matrices
A = amat(N);
Bs = cell(M,1);
for c = 1:M
    Bs{c} = bmat(N, a(c), b(c));
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
        for c = 1:M
            x - (m(c)-r)*delta == A*vec(Fl(:,:,c)) + Bs{c}*vec(Gl(:,:,c))
            (m(c)+r)*delta - x == A*vec(Fu(:,:,c)) + Bs{c}*vec(Gu(:,:,c))
            Fl(:,:,c) >= 0
            Gl(:,:,c) >= 0
            Fu(:,:,c) >= 0
            Gu(:,:,c) >= 0
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