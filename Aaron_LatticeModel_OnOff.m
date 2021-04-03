%June 2nd, 2020 - Proteins binding on to lattice with Scatchard plot
%This program models a protein and lattice (both variable lengths). The
%proteins bind to the lattice at random locations and unbind from random
%locations based on probabilities. The program will then create a Scatchard
%plot based on Eq. 10 in the McGhee paper. The initial concentration of the
%free protein and the association constant are also included.

N = 1000;           %variable length of the DNA lattice
n = 5;              %variable length of the protein
L = 1;              %initial concentration of free proteins
K = 1;              %K, the association constant
BindProb = 0.8;     %binding probability (between 0 and 1)
UnbindProb = 0.3;   %unbinding probability (between 0 and 1)

DNA = zeros(1,N);           %matrix to represent the DNA lattice (0=empty, 1=filled)
Iterations = 1000;         %number of iterations of the code

BindHist = zeros(1,Iterations);        %matrix to observe history of random binding locations
UnbindHist = zeros(1,Iterations);      %matrix to observ history of random unbinding locations
RandBHist = zeros(1,Iterations);       %matrix to observe history of random binding numbers
RandUHist = zeros(1,Iterations);       %matrix to observe history of random unbinding numbers
ScatchY = zeros(1,Iterations);         %allocates memory for y-axis of Scatchard plot
v = zeros(1,Iterations);               %allocates memory for binding density

for i = 1:Iterations
    RandB = rand;             %random number to check binding probability
    RandBHist(i) = RandB;     %stores RandB value in RandBHist for referencing
    RandU = rand;             %random number to check unbinding probability
    RandUHist(i) = RandU;     %stores RandU value in RandUHist for referencing
    if RandB <= BindProb               %only runs if RandB is <= to binding probability
        Bind = randi((N-n),1);         %Bind is a random location on the DNA for proteins to bind
        BindHist(i) = Bind;            %stores Bind location in BindHist for referencing
        if DNA(Bind:Bind+(n-1)) == 0   %checks if random locations are empty
            DNA(Bind:Bind+(n-1)) = 1;  %fills corresponding sections with protein
        end
    end
    Go = 1;             %variable for while statement
    while Go == 1       %while statement to guarantee random unbinding location is start of a protein
        pos = randi(length(BindHist(1:i)));       %picks position a protein has already bound to
        Unbind = BindHist(pos);                   %sets Unbind to that position
        if Unbind == 1               %if Unbind == 1
            Go = 0;                  %move forward to check if full protein is present (should be)
        else
            if Unbind == 0          %if the unbind location selected is zero
                Go = 1;             %pick a new location
            else
                A = rem(sum(DNA(1:(Unbind-1))),n);  %otherwise calculate remainder of sum of DNA prior to this location
                if A == 0                 %if remainder equals 0
                    Go = 0;               %we've picked the start of a protein              
                end
            end
        end
    end
    if RandU <= UnbindProb           %only runs if RandU is less than Unbinding probability
        if (Unbind == 1) & (DNA(1:n) == 1)   %if a protein is present at the random location   
            DNA(1:n) = 0;                    %unbind the protein
        else
            if DNA(Unbind:Unbind+(n-1)) == 1    %if a protein exists at random location
                DNA(Unbind:Unbind+(n-1)) = 0;   %unbind the protein
            end
        end
    end
    B = (sum(DNA))/n;         %total number of proteins bound to the lattice
    v(i) = B/N;               %calculates binding density
    ScatchY = K*(1-(n.*v)).*(((1-(n.*v))./(1-((n-1).*v))).^(n-1));  %y-axis of Scatchard plot (Eq. 10 in McGhee paper)
end

figure();
plot(v,ScatchY, 'b');   %generates Scatchard plot
hold on;
xlabel('v');
xlim([0,1/n]);
ylabel('v/L');
ylim([0,K]);
title('Scatchard Plot (K = 1, N = 1000)');

% X = 1:Iterations;
% 
% figure();
% scatter(X,v,5,'r','filled');    %plots binding density over time
% xlabel('"time"');
% ylabel('Binding Density');