%%
clc
clear
close all

%% Polynomial Design
N = 12;

a = [0.0,    0.7];
b = [0.3,    1.0];
v = [0.0,    0.1];

x = pd(N, a, b, v);


%% Plot Magnitude Response 
i = 0:N;
t = linspace(0,1,1000)';
T = bsxfun(@power, t, i);
X = T * x;

figure, plot(t, X);
xlabel('t');
ylabel('Reponse');
title('Response');
