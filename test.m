%%
clc
clear
close all

%%
N = 20;

a = [-pi/2,     -pi,        2*pi/3];
b = [pi/2,      -2*pi/3,    pi];
l = [0.95,      0,          0];
u = [1.05,      0.05,       0.05];
v = [1.0,       0.0,        0.0];

x = fdr(N, a, b, v);



w = linspace(-pi, pi, 1000);
X = fftshift(fft(x, 1000));
figure, plot(w, abs(X));
xlabel('\omega');
ylabel('Magnitude');
title('Magnitude response');
