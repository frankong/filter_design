%%
clc
clear
close all

%%
N = 12;

a = [-pi/2,     -pi,        2*pi/3];
b = [pi/2,      -2*pi/3,    pi];
l = [0.98,      0,          0];
u = [1.02,      0.02,       0.02];
v = [1.0,       0.0,        0.0];

xr = fdzr(N, a, b, v);
xe = fdze(N, a, b, l, u);

w = linspace(-pi, pi, 1000);
Xr = fftshift(fft(xr, 1000));
Xe = fftshift(fft(xe, 1000));
figure, plot(w, abs(Xr), w, abs(Xe));
xlabel('\omega');
ylabel('Magnitude');
title('Magnitude response');
legend('Minimum ripple', 'Minimum energy');
