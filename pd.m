function x = pd(N, a, b, m)
% Design a degree-2N polynomial that minimizes maximum absolute deviation from the specified magnitude response
% 
% Inputs:
% 
%       N - Polynomial degree.
%       a - Bands start points. From 0 to 1.
%       b - Bands end points. From 0 to 1.
%       m - Magnitude response at intevals.
% 
%       a, b, m are length-B vectors where B is the number of bands.
% Output:
% 
%       x - Degree N polynomial
%
% Reference:
% Page 11 eq 1.30 and 1.31 in Dumitrescu 2007

assert(length(a) == length(b) & length(b) == length(m));
M = length(a);
n = floor(N/2);
e = mod(N,2) == 0;

% Construct matrices
As = cell(M,1);
Bs = cell(M,1);
for i = 1:M
    if e
        As{i} = aemat(n);
        Bs{i} = bemat(n, a(i), b(i));
    else
        As{i} = aomat(n, a(i));
        Bs{i} = bomat(n, b(i));
    end
end
dirac = sparse(1, 1, 1.0, N+1, 1);

% Solve
cvx_begin SDP
    variable x(N+1)             complex
    variable Fl(n+1, n+1, M)    hermitian
    variable Gl(n-e+1, n-e+1, M)        hermitian
    variable Fu(n+1, n+1, M)    hermitian
    variable Gu(n-e+1, n-e+1, M)        hermitian
    variable r
    minimize r
    subject to
        for i = 1:M
            x - (m(i)-r)*dirac == As{i}*vec(Fl(:,:,i)) + Bs{i}*vec(Gl(:,:,i))
            (m(i)+r)*dirac - x == As{i}*vec(Fu(:,:,i)) + Bs{i}*vec(Gu(:,:,i))
            Fl(:,:,i) >= 0
            Gl(:,:,i) >= 0
            Fu(:,:,i) >= 0
            Gu(:,:,i) >= 0
        end
cvx_end

end

function A = aemat(n)
    si = zeros((n+1)^2,1);
    sj = zeros((n+1)^2,1);
    ss = zeros((n+1)^2,1);
    
    c = 1;
    for i = 0:n
        for j = 0:n
            si(c) = i + j + 1;
            sj(c) = (i+1) + j * (n+1);
            ss(c) = 1.0;
            c = c + 1;
        end
    end
    A = sparse(si, sj, ss, 2*n+1, (n+1)^2);
end

function B = bemat(n, a, b)
si = zeros(3*n^2,1);
sj = zeros(3*n^2,1);
ss = zeros(3*n^2,1);

c = 1;
for i = 0:n-1
    for j = 0:n-1
        si(c) = i + j + 1;
        sj(c) = (i+1) + j * n;
        ss(c) = -a*b;
        c = c + 1;
        
        
        si(c) = i + j + 1 + 1;
        sj(c) = (i+1) + j * n;
        ss(c) = a+b;
        c = c + 1;
        
        
        si(c) = i + j + 1 + 2;
        sj(c) = (i+1) + j * n;
        ss(c) = -1;
        c = c + 1;
    end
end
B = sparse(si, sj, ss, 2*n+1, n^2);
end


    
function A = aomat(n, a)
    si = zeros(2*(n+1)^2,1);
    sj = zeros(2*(n+1)^2,1);
    ss = zeros(2*(n+1)^2,1);
    
    c = 1;
    for i = 0:n
        for j = 0:n
            si(c) = i + j + 1;
            sj(c) = (i+1) + j * (n+1);
            ss(c) = -a;
            c = c + 1;
            
            
            si(c) = i + j + 1 + 1;
            sj(c) = (i+1) + j * (n+1);
            ss(c) = 1;
            c = c + 1;
        end
    end
    A = sparse(si, sj, ss, 2*n+2, (n+1)^2);
end


function B = bomat(n, b)
    si = zeros(2*(n+1)^2,1);
    sj = zeros(2*(n+1)^2,1);
    ss = zeros(2*(n+1)^2,1);
    
    c = 1;
    for i = 0:n
        for j = 0:n
            si(c) = i + j + 1;
            sj(c) = (i+1) + j * (n+1);
            ss(c) = b;
            c = c + 1;
            
            
            si(c) = i + j + 1 + 1;
            sj(c) = (i+1) + j * (n+1);
            ss(c) = -1;
            c = c + 1;
        end
    end
    B = sparse(si, sj, ss, 2*n+2, (n+1)^2);
end