clear all;
close all;

% Simple recreation of Scatchard plots. Uses equations from MVH paper for
% cooperative and noncooperative systems.

K = 1;  %kinetic association constant

n_values = [1,2,3,5,10];  %protein length

for n = n_values
    v = 0:(1/n)/1000:1/n;    %x-values for scatchard plot
    Y_N = K*(1-(n*v)).*(((1-(n*v))./(1-((n-1)*v))).^(n-1));

    figure(1);
    plot(v,Y_N,'LineWidth',2);
    hold on;
end
figure(1);
xlim([0 max(1./n_values)]);
xlabel('v');
ylim([0 K]);
ylabel('v/L');
title('Non-Cooperative Scatchard Plot');
box on;

Legend_n = cell(length(n_values),1);
for a = 1:length(n_values)
    Legend_n{a} = ['n = ', num2str(n_values(a))];
end
figure(1);
legend(Legend_n);

w_values = [0.1,0.5,1,3,5];  %cooperativity constants
n_C = 1;    %protein size for cooperative Scatchard plot

for w = w_values
    x = 0:(1/n_C)/1000:(1/n_C);
    
    if w ~= 1
        R = sqrt(((1-((n_C+1)*x)).^2)+4*w*x.*(1-(n_C*x)));
        Y_C = K*(1-(n_C*x)).*(((((2*w)-1).*(1-(n_C*x))+x-R)./(2*(w-1).*(1-(n_C*x)))).^(n_C-1)).*(((1-((n_C+1)*x)+R)./(2*(1-x))).^2);
    else
        Y_C = K*(1-(n_C*x)).*(((1-(n_C*x))./(1-((n_C-1)*x))).^(n_C-1));
    end

    figure(2);
    plot(x,Y_C,'LineWidth',2);
    hold on;
end
figure(2);
xlim([0 1/n_C]);
xlabel('v');
ylim([0 3]);
ylabel('v/L');
title('Cooperative Scatchard Plots');
box on;

Legend_w = cell(length(w_values),1);
for b = 1:length(w_values)
    Legend_w{b} = ['\omega = ', num2str(w_values(b))];
end
figure(2);
legend(Legend_w);