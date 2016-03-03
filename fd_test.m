%%
clc
clear
close all

%% Polynomial Design
N = 12;

a = [-pi/2,     -pi,        2*pi/3];
b = [pi/2,      -2*pi/3,    pi];
m = [1.0,       0.0,        0.0];

x = fd(N, a, b, m);

%% Plot Magnitude Response 
w = linspace(-pi, pi, 1000);
X = fftshift(fft(x, 1000));
figure, plot(w, abs(X));
xlabel('\omega');
ylabel('Magnitude');
title('Magnitude response');
