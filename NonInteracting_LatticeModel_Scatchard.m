clear all;
close all;

% This code will randomly bind and unbind proteins to a DNA lattice
% represented by an array of zeros and ones. A zero represents a free
% location on the lattice while a one represents a location where a protein
% is bound. The length of the lattice and the proteins are variable. The
% probabilities of binding and unbinding are based on k_on and k_off values
% as well as chemical kinetics. Proteins cannot interact on the lattice.
% The code creates a Scatchard plot (v/L vs. v) with each point
% representing an equilibrium binding density for the lattice saturation.
% k_on, k_off, L, and dt values must be selected such that k_on*L*dt < 1 as
% well as k_off*dt < 1.

N = 100;    %length of DNA lattice
n = 3;      %protein length
k_on = 1;   %kinetic rate constant for binding
k_off = 1;  %kinetic rate constant for unbinding
K = k_on/k_off; %calculates equilibrium constant

currentBound = zeros(1,N);  %allocate memory for currentBound array
Iterations = 500;  %how many binding/unbinding loops to run
t(1) = 0;

for i = 1:1000    %for loop to vary L (free protein concentration)
    DNA = zeros(1,N);   %array to model DNA lattice
    
    vLoop(1) = 0;
    xAB(1) = 0;
    
    L(i) = i/40;   %stores values for L
    dt = 0.02;     %small time step required for probability calculations
    BindCounter = 0;    %counts how many binding events occur
    UnbindCounter = 0;   %counts how many unbinding events occur
    for a = 1:Iterations   %loops binding and unbinding runs multiple times
        for j = randperm(N-(n-1))   %for loop checks each location on DNA in random order
            if DNA(j:j+(n-1)) == 0  %checks availiability at location
                if rand <= k_on*L(i)*dt    %checks probability based on kinetics
                    DNA(j:j+(n-1)) = 1; %space matches all requirements so protein binds
                    BindHist(BindCounter+1) = j;  %array to record history of all past binding locations
                    BindCounter = BindCounter+1;    %increases counter of total binding events
                    currentBound(j) = 1;    %shows which positions a protein is bound to currently
                end
            end
        end
        for m = find(currentBound == 1) %picks each location where a protein is currently bound
            if rand <= k_off*dt    %checks probability based on kinetics
                DNA(m:m+(n-1)) = 0; %unbinds protein from location
                UnbindHist(UnbindCounter+1) = m;    %stores location in history of all past unbinding events
                UnbindCounter = UnbindCounter+1;    %increases counter of total unbinding events
                currentBound(m) = 0;    %shows no more protein bound to that location
            end
        end
        vLoop(a+1) = (sum(DNA)/n)/N;  %calculates binding density after each iteration
        t(a+1) = t(a)+dt;
    end
    v(i) = mean(vLoop(Iterations-50:Iterations));   %records the equilibrium (?) binding density
    ScatchY(i) = v(i)./L(i);   %calculates values of v/L for the model
end

x = 0:(1/n)/1000:1/n;
TheoreticalScatchY = K.*(1-(n.*x)).*(((1-(n.*x))./(1-((n-1).*x))).^(n-1));   %calculates theoretical values of v/L using Eq. 10 of McGhee paper


% figure(1);
% plot(t,vLoop);          %plots binding density with each iteration of binding and unbinding to show equilibrium
% xlabel('Iterations');
% xlim([0 Iterations]);
% ylabel('Binding Density');
% ylim([0 max(vLoop+0.01)]);

figure(2);
scatter(v,ScatchY, 5,'r','filled');        %plots Scatchard plot with both real data and theoretical values
hold on;
plot(x,TheoreticalScatchY,'k');
xlabel('v');
xlim([0 1/n]);
ylabel('v/L');
ylim([0 K+0.25*K]);
title('Scatchard Plot');
legend('Real Data','Theoretical');

% figure();
% scatter(L,v,5,'r','filled');          %plots v vs. L just to see possible relationships
% xlabel('L (Free Protein Concentration)');
% ylabel('v (Binding Density)');
% title('Equilibrium Density vs. Protein Concentration');

% cftool;                    %opens cftool to fit a curve to the data points of the Scatchard plot