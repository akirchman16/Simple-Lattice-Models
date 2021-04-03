clear all;
close all;

% This code will bind and unbind proteins to a DNA lattice without any
% cooperativity. No Scatchard plot will be produced but the fractional
% coverage of the DNA over time can be plotted.

N = 1000;    %length of DNA lattice
n = 3;      %protein length
k_on = 1;   %kinetic rate constant for binding
k_off = 1;  %kinetic rate constant for unbinding
L = 0.5;    %concentration of free proteins

dt = 0.001; %time step between loops
K = k_on/k_off; %calculates equilibrium constant

Iterations = 1000;  %how many binding/unbinding loops to run
currentBound = zeros(1,N);  %allocate memory for currentBound array
xAB = 0; %initially an empty lattice
t(1) = 0;
ProteinCount(1) = xAB;
FracCover(1) = xAB/N;

DNA = zeros(1,N);   %array to model DNA lattice

BindCounter = 0;    %counts how many binding events occur
UnbindCounter = 0;   %counts how many unbinding events occur

for a = 1:Iterations   %loops binding and unbinding runs multiple times
    for j = randperm(N-(n-1))   %for loop checks each location on DNA in random order
        if DNA(j:j+(n-1)) == 0  %checks availiability at location
            if rand <= k_on*L*dt    %checks probability based on kinetics
                DNA(j:j+(n-1)) = 1; %space matches all requirements so protein binds
                BindHist(BindCounter+1) = j;  %array to record history of all past binding locations
                BindCounter = BindCounter+1;    %increases counter of total binding events
                currentBound(j) = 1;    %shows which positions a protein is bound to currently
                xAB = xAB+1;
            end
        end
    end
    for m = find(currentBound == 1) %picks each location where a protein is currently bound
        if rand <= k_off*dt    %checks probability based on kinetics
            DNA(m:m+(n-1)) = 0; %unbinds protein from location
            UnbindHist(UnbindCounter+1) = m;    %stores location in history of all past unbinding events
            UnbindCounter = UnbindCounter+1;    %increases counter of total unbinding events
            currentBound(m) = 0;    %shows no more protein bound to that location
            xAB = xAB-1;
        end
    end
    t(a+1) = t(a)+dt;
    ProteinCount(a+1) = sum(DNA)/n;
    FracCover(a+1) = xAB/N;
end

figure();
scatter(t,FracCover,2,'b','filled');    %plots fractional coverage over time
hold on;
xlabel('Time, t (s)');
xlim([0 max(t)]);
ylabel('Fractional Coverage');
ylim([0 1]);
title(['Fractional Coverage (K = ' num2str(K) ', N = ' num2str(N) ', n = ' num2str(n) ')']);
legend([ num2str(Iterations) ' Events']);