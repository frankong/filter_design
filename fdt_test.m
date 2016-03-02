%%
clc
clear
close all

%%
N = 12;

a = [0,        0.8];
b = [0.1,     1.0];
l = [-0.01,      0.59];
u = [0.01,      0.61];
v = [0.0,      0.6];

xr = fdtr(N, a, b, v);
xe = fdte(N, a, b, l, u);

i = 0:N;
t = linspace(0,1,1000)';
T = bsxfun(@power, t, i);
Xr = T * xr;
Xe = T * xe;

figure, plot(t, Xr, t, Xe);
xlabel('t');
ylabel('Reponse');
title('Response');
legend('Minimum ripple', 'Minimum energy');
