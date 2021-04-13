clear all;
close all;

% Simple recreation of Scatchard plots. Uses equations from MVH paper for
% cooperative and noncooperative systems.

K = 1;  %kinetic association constant
n = 3;  %protein length

w = 1;  %cooperativity constant

x = [0:(1/n)/1000:1/n];    %x-values for scatchard plot

Y_NonCoop = K.*(1-(n.*x)).*((1-(n.*x))/(1-(n-1).*x)).^(n-1);

figure(1);
plot(x,Y_NonCoop,'b');
hold on;
%plot(v,Y_Coop,'r');
xlim([0 max(x)]);
xlabel('v');
%ylim([0 (max(Y_Coop)+0.1*max(Y_Coop))]);
ylabel('v/L');
title('Scatchard Plots');
legend('\omega = 1',['\omega = ', num2str(w)]);