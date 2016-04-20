%%
clc
clear
close all

%% Polynomial Design
N = 12;

a = [-1,     -pi,        1.5];
% a = [-pi/2,     -pi,        2*pi/3];
b = [1,      -1.5,    pi];
% b = [pi/2,      -2*pi/3,    pi];
m = [1.0,       0.0,        0.0];

x = fd(N, a, b, m);

%% Plot Magnitude Response 
n = -N:N;
figure, plot(n, real(x), n, imag(x));
xlabel('n');
title('Filter');
legend('real','imaginary');

w = linspace(-pi, pi, 1000);
X = fftshift(fft(fftshift(zpad(x, 1000))));
figure, plot(w, abs(X));
axis([-4, 4, 0, 1.2])
xlabel('\omega');
ylabel('Magnitude');
title('Magnitude response');

figure, plot(w, angle(X));
xlabel('\omega');
ylabel('Magnitude');
title('Phase response');