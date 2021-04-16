clear all;
close all;

% This code will randomly bind and unbind proteins to a DNA lattice
% represented by an array of zeros and ones. A zero represents a free
% location on the lattice while a one represents a location where a protein
% is bound. The length of the lattice and the proteins are variable. The
% probabilities of binding and unbinding are based on k_on and k_off values
% as well as chemical kinetics. Proteins can interact on the lattice and
% binding probabilities are affected by these interactions. The code
% creates a growth profile over time of the given situation based on the
% given parameters.
% k_on, k_off, L, and dt values must be selected such that k_on*L*dt < 1 as
% well as k_off*dt < 1.
% There are 'dummy zeros' on either end of the lattice (positions 1 and
% N+1) in order for the code to check the binding location type.


N = 8660;    %length of DNA lattice
n = 3;      %protein length
k_on = 0.1;   %kinetic rate constant for binding
k_off = 1;  %kinetic rate constant for unbinding
w = 1;  %co-operativity parameter
L = 8;  %concentration of free proteins

dt = 0.01;     %small time step required for probability calculations

Probabilities = [k_on*L*dt,k_on*L*w*dt,k_on*L*(w^2)*dt,k_off*dt];   %probabilites used for comparisons
if max(Probabilities) > 1
    disp('ERROR: dt TOO BIG');
elseif max(Probabilities) < 0.001
    disp('ERROR: dt TOO SMALL');
end

K = k_on/k_off; %calculates equilibrium constant

Iterations = 500;  %how many binding/unbinding loops to run

DNA = zeros(1,N+2);   %array to model DNA lattice with 'dummy zeros' on each end
currentBound = zeros(1,N+2);  %allocate memory for currentBound array with 'dummy zeros'

BindHist = 0;   %resets BindHist for each loop of binding/unbinding events
BindCounter = 0;    %counts how many binding events occur
isolBindCounter = 0;    %counts isolated binding events
scBindCounter = 0;  %counts singly contiguous binding events
dcBindCounter = 0;  %counts doubly contiguous binding events
UnbindCounter = 0;   %counts how many unbinding events occur

t = zeros(1,Iterations);    %memory allocation of t and FracCover
FracCover = zeros(1,Iterations);

t(1) = 0;       %initial conditions
FracCover(1) = 0;
for a = 2:Iterations   %loops binding and unbinding runs multiple times
    for j = randperm(N-(n-1))+1   %for loop checks each location on DNA in random order (binding events)
        if DNA(j:j+(n-1)) == 0  %checks availiability at location
            if DNA(j-1) == 0 && DNA(j+n) == 0  %checks if binding site is isolated
                if rand <= k_on*L*dt    %checks probability based on kinetics
                    DNA(j:j+(n-1)) = 1; %space matches all requirements so protein binds
                    BindHist(1,BindCounter+1) = j;  %array to record history of all past binding locations
                    BindCounter = BindCounter+1;    %increases counter of total binding events
                    isolBindHist(isolBindCounter+1) = j;    %stores location of all isolated binding events
                    isolBindCounter = isolBindCounter+1;    %increases counter of isolated binding events
                    currentBound(j) = 1;    %shows which positions a protein is bound to currently
                end
            end
            if (DNA(j-1) == 0 && DNA(j+n) == 1) || (DNA(j-1) == 1 && DNA(j+n) == 0)    %checks if binding site is singly contiguous (sc)
                if rand <= k_on*L*w*dt    %checks probability based on kinetics
                    DNA(j:j+(n-1)) = 1; %space matches all requirements so protein binds
                    BindHist(2,BindCounter+1) = j;  %array to record history of all past binding locations
                    BindCounter = BindCounter+1;    %increases counter of total binding events
                    scBindHist(scBindCounter+1) = j;    %stores location of all sc binding events
                    scBindCounter = scBindCounter+1;    %increases counter of all sc binding events
                    currentBound(j) = 1;    %shows which positions a protein is bound to currently
                end
            end
            if DNA(j-1) == 1 && DNA(j+n) == 1  %checks if binding site is doubly contiguous (dc)
                if rand <= k_on*L*(w^2)*dt    %checks probability based on kinetics
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
    
    FracCover(a) = sum(DNA)/N;  %fractional coverage of DNA lattice
    
    t(a) = t(a-1)+dt;
end

figure(1);
scatter(t,FracCover,3,'r','filled');
hold on;
xlabel('Time, t');
ylabel('Saturation');
xlim([0 dt*Iterations]);
ylim([0 1]);
title('Saturation of DNA Lattice');
box on;