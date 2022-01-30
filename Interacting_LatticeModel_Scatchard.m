clear all;
close all;

% This code will randomly bind and unbind proteins to a DNA lattice
% represented by an array of zeros and ones. A zero represents a free
% location on the lattice while a one represents a location where a protein
% is bound. The length of the lattice and the proteins are variable. The
% probabilities of binding and unbinding are based on k_on and k_off values
% as well as chemical kinetics. Proteins can interact on the lattice and
% binding probabilities are affected by these interactions. The code 
% creates a Scatchard plot (v/L vs. v) with each point representing an 
% equilibrium binding density for the lattice saturation.
% k_on, k_off, L, and dt values must be selected such that k_on*L*dt < 1 as
% well as k_off*dt < 1.
% There are 'dummy zeros' on either end of the lattice (positions 1 and
% N+1) in order for the code to check the binding location type.


N = 8660;    %length of DNA lattice
n = 5;      %protein length
k_on = 1;   %kinetic rate constant for binding
k_off = 1;  %kinetic rate constant for unbinding
K = k_on/k_off; %calculates equilibrium constant
w = 0.1;  %co-operativity parameter

Iterations = 200;  %how many binding/unbinding loops to run

dt = 0.01;     %small time step required for probability calculations

for i = 1:1000    %for loop to vary L (free protein concentration)
    DNA = zeros(1,N+2);   %array to model DNA lattice with 'dummy zeros' on each end
    currentBound = zeros(1,N+2);  %allocate memory for currentBound array with 'dummy zeros'
    L(i) = i/400;   %stores values for L
    BindHist = 0;   %resets BindHist for each loop of binding/unbinding events
    BindCounter = 0;    %counts how many binding events occur
    isolBindCounter = 0;    %counts isolated binding events
    scBindCounter = 0;  %counts singly contiguous binding events
    dcBindCounter = 0;  %counts doubly contiguous binding events
    UnbindCounter = 0;   %counts how many unbinding events occur
    for a = 1:Iterations   %loops binding and unbinding runs multiple times
        for j = randperm(N-(n-1))+1   %for loop checks each location on DNA in random order (binding events)
            if DNA(j:j+(n-1)) == 0  %checks availiability at location
                if DNA(j-1) == 0 && DNA(j+n) == 0  %checks if binding site is isolated
                    if rand <= k_on*L(i)*dt    %checks probability based on kinetics
                        DNA(j:j+(n-1)) = 1; %space matches all requirements so protein binds
                        BindHist(1,BindCounter+1) = j;  %array to record history of all past binding locations
                        BindCounter = BindCounter+1;    %increases counter of total binding events
                        isolBindHist(isolBindCounter+1) = j;    %stores location of all isolated binding events
                        isolBindCounter = isolBindCounter+1;    %increases counter of isolated binding events
                        currentBound(j) = 1;    %shows which positions a protein is bound to currently
                    end
                end
                if (DNA(j-1) == 0 && DNA(j+n) == 1) || (DNA(j-1) == 1 && DNA(j+n) == 0)    %checks if binding site is singly contiguous (sc)
                    if rand <= k_on*L(i)*w*dt    %checks probability based on kinetics
                        DNA(j:j+(n-1)) = 1; %space matches all requirements so protein binds
                        BindHist(2,BindCounter+1) = j;  %array to record history of all past binding locations
                        BindCounter = BindCounter+1;    %increases counter of total binding events
                        scBindHist(scBindCounter+1) = j;    %stores location of all sc binding events
                        scBindCounter = scBindCounter+1;    %increases counter of all sc binding events
                        currentBound(j) = 1;    %shows which positions a protein is bound to currently
                    end
                end
                if DNA(j-1) == 1 && DNA(j+n) == 1  %checks if binding site is doubly contiguous (dc)
                    if rand <= k_on*L(i)*(w^2)*dt    %checks probability based on kinetics
                        DNA(j:j+(n-1)) = 1; %space matches all requirements so protein binds
                        BindHist(3,BindCounter+1) = j;  %array to record history of all past binding locations
                        BindCounter = BindCounter+1;    %increases counter of total binding events
                        dcBindHist(dcBindCounter+1) = j;    %stores location of all dc binding events
                        dcBindCounter = dcBindCounter+1;    %increases counter of all dc binding events
                        currentBound(j) = 1;    %shows which positions a protein is bound to currently
                    end 
                end
           
            end
        end
        for m = find(currentBound == 1) %picks each location where a protein is currently bound (unbinding events)
            if rand <= k_off*dt    %checks probability based on kinetics
                DNA(m:m+(n-1)) = 0; %unbinds protein from location
                UnbindHist(UnbindCounter+1) = m;    %stores location in history of all past unbinding events
                UnbindCounter = UnbindCounter+1;    %increases counter of total unbinding events
                currentBound(m) = 0;    %shows no more protein bound to that location
            end
        end
        vLoop(a) = (sum(DNA)/n)/N;  %calculates binding density after each iteration
    end
    v(i) = mean(vLoop(Iterations-50:Iterations));   %records the equilibrium (?) binding density
    R(i) = sqrt(((1-((n+1).*v(i))).^2)+(4*w.*v(i).*(1-(n.*v(i))))); %R value necessary for theoretical values calculation
    ScatchY(i) = v(i)./L(i);   %calculates values of v/L for the model
    if w == 1
        TheoreticalScatch(i) = K*(1-(n.*v(i))).*(((1-(n.*v(i)))./(1-((n-1).*v(i)))).^(n-1));   %calculates theoretical values of v/L using Eq. 10 of McGhee paper
    else
        TheoreticalScatch(i) = K*(1-n.*v(i))*(((((2*w)-1)*(1-(n.*v(i)))+v(i)-R(i))/(2*(w-1)*(1-n.*v(i)))).^(n-1))*(((1-((n+1).*v(i))+R(i))/(2*(1-(n.*v(i))))).^2);   %calculates theoretical values of v/L using Eq. 15 of McGhee paper
    end
end

% x = 1:Iterations;
% figure();              %plots binding density vs. iteration to show equilibrium
% plot(x,vLoop);          %plots binding density with each iteration of binding and unbinding
% xlabel('Iterations');
% xlim([0 Iterations]);
% ylabel('Binding Density');
% ylim([0 max(vLoop+0.01)]);

figure();
scatter(v,ScatchY, 5, 'r', 'filled');        %plots Scatchard plot with both real data and theoretical values
hold on;
plot(v,TheoreticalScatch,'b');
xlabel('v');
xlim([0 1/n]);
ylabel('v/L');
ylim([0 2*K]);
title('Scatchard Plot'); box on;
legend('Simulation','McGhee & von Hippel Model');

% cftool;         %opens cftool to fit a curve to the data points in Scatchard plot
% curve fit equation (copy/paste into 'custom fit' in cftool)???:
% (K-(K*n*x))*(((((2*w)-1)*(1-(n*x))+x-sqrt(((1-(n+1)*x)^(2))+4*w*x*(1-(n*x))))/(2*(w-1)*(1-(n*x))))^(n-1))*(((1-((n+1)*x)+sqrt(((1-(n+1)*x)^(2))+4*w*x*(1-(n*x))))/(2*(1-(n*x))))^2)

% figure();
% scatter(L,v,5,'r','filled');          %plots v vs. L just to see possible relationships
% xlabel('L (Free Protein Concentration)');
% ylabel('v (Binding Density)');
% title('Equilibrium Density vs. Protein Concentration');