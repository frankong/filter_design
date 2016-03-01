%%
clc
clear
close all

%%
N = 30;

a = [-pi/2];
b = [pi/2];
l = [0.99];
u = [1.01];

x = filt_design(N, a, b, l, u);

w = linspace(-pi, pi, 1000);
X = fftshift(fft(x, 1000));
figure, plot(w, abs(X));
xlabel('\omega');
ylabel('Magnitude');
title('Magnitude response');
