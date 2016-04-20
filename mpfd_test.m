%%
clc
clear
close all

%% Minimum Phase Filter Design
N = 25;

% a = [-0.5,     -1.0,       0.7] * pi;
% b = [0.5,      -0.7,       1.0] * pi;

a = [-1,     -pi,        1.5];
% a = [-pi/2,     -pi,        2*pi/3];
b = [1,      -1.5,    pi];
m = [1.0,       0.0,        0.0];

[x, G] = mpfd(N, a, b, m);

%% Plot
n = 0:N;
figure, plot(n, real(x), n, imag(x));
xlabel('n');
title('Filter');
legend('real','imaginary');

w = linspace(-pi, pi, 1000);
% X = fftshift(fft(fftshift(zpad(x, 1000))));
X = fftshift(fft(x, 1000));
figure, plot(w, abs(X));
axis([-4, 4, 0, 1.2])
xlabel('\omega');
ylabel('Magnitude');
title('Magnitude response');

figure, plot(w, -unwrap(angle(X)));
xlabel('\omega');
ylabel('Magnitude');
title('Phase response');