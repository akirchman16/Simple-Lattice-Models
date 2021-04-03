%A demonstration of chemical kinetics in a discrete-time system

%In a math or chemistry class, we'd start by taking the concentration of
%"reacted" stuff and giving it a name like "A".  In this case, our reaction
%is "binding", and we don't have a concentration so much as a number of
%bound molecules.  That's fine, we'll still call it A.  How does A change
%over time?  Well, at the beginning, if nothing is bound, we can define its
%rate of change by 

%  dA/dt = k_on * L         (1)

%k_on is the kinetic rate constant for the molecules to bind, and L is the
%concentration of free molecules.  This should make sense; if you have more
%molecules to begin with, they should be finding free binding locations and
%binding on faster.  However, this raises a problem: we have to take into
%account how many remaining locations there are to bind, since we now
%that this number will change over time.  So we get

%  dA/dt = k_on * L * (number of possible binding locations)        (2)

%If all the proteins have an n=1, then the remaining locations would just 
%be the total number of possible binding locations (let's call that N)
%minus the number currently occupied (again, equal to A for n=1).  

%  dA/dt = k_on * L * (N-A)                 (3)

%Normally, chemists talk about these equations in terms of the
%concentrations being in Molar (moles/liter); here, it's possible we could
%take L to be a whole number of molecules (if we're simulating a small
%enough number of them).  But the goal is not to keep track of every single
%protein when it isn't bound, but instead to just imagine that "there's a
%lot out there" in the solution.  For this reason, when a protein enters or
%exits the lattice, L is not appreciably changed.  We will be treating it
%as a "constant" from here on out (though really it's more of a parameter,
%since we'll be changing it a lot for different runs of the code to see how
%it affects things).

%So what we're left with is a rate constant that has units something like
%"protein units per second per Molar of L per binding location"
%That sounds clunky, but it's just to get it to cancel out with the L*(N-A) and
%produce the units of dA/dt (which are protein units per second).  The most
%important part is the "per second" bit, which is the key to understanding
%how our kinetic model will be working.

%That's that for binding.  How about the other half?  
%Once a lot of proteins have bound, there must be some chance for them to
%dissociate.  This chance should be proportional to the number of proteins
%that are bound (each of them has some small likelihood) and...that's it.
%The size of L (in solution) does not affect a bound protein to fall off.
%This means that if there were no binding going on, the description of our
%unbinding would be simply:

% dA/dt = - k_off * A                       (4)

%(the negative sign is because every unbinding event decreases A)
%(The units of k_off are just "per second")
%Say, that's the differential equation describing an exponential decay!
%However, in real life we have a competition between binding and unbinding
%at any given moment, so the total change in bound protein, A, is:

% dA/dt = k_on * L * (N-A) - k_off * A              (5)

%This is already pretty cool.  You could plug that into a differential
%equation solver if you wanted.  A few more things we can say about it:

%The equilibrium constant K describes the limiting state of the reaction.
%In other words, you let this thing go long enough, it balances out to a
%stable population of A, and we can calculate this using the equilibrium
%quotient:

%      K = (A)/((N-A)*L)    (6)

%But what's even cooler is that, in the way we've written it, K is
%dependent on only the rate constants from earlier:

%     K = k_on/k_off                (7)

%You can rearrange equations (6) and (7) and plug them into
%equation (5) and see that they make dA/dt = 0; a steady state.  

%The fact that it works this way is one of the coolest results in physical
%chemistry, in that it hints at deep symmetries and fundamental
%relationships that are not at all obvious from looking at the separate
%pieces.  That's good physics.

%Some of the above equations won't work perfectly when n is not 1 (the
%number of available spots is not just N-A, as McGhee pointed out and you
%well know by now) but the concepts are still correct.


% This has been a brief introduction into kinetics.  But how do we get this
% to happen in our code?  We must make the jump from continuous quantities
% to discrete ones; time happens in chunks, and you can't add a fraction of
% a protein.  We can take our guidance from the above equations.

%First, consider how to make something go onto the lattice.  The rate must
%be proportional to L, some constant k_on, but also to the number of
%available slots.  How do we get that to happen?  Well, in each time step,
%you could simply iterate over the whole lattice, check each spot to see if
%it's available, and if it is available, give it the opportunity to bind a
%protein with a certain chance.  There is one flaw in this, and that is
%that it biases certain spots to receive proteins more often than others.
%t's a better idea to choose a random order to go over the lattice each
%time, which you can do quickly with a for loop that looks like this:
lengthLattice = 500;
for i = randperm(lengthLattice)
    %i is now the index of a spot on the lattice
    %i will be every number from 1 to 500 in a random order
    %check to see if this spot is unoccupied/has a space open enough for a
    %protein to bind
    %if it's available, give a protein a chance to bind
    %if something bound, update whatever lists you use to keep track of
    %what spots are free/have a space open enough for something to bind
end

%randperm is just a fun function for things like this.  

%but how do we "give something a chance to bind"? Well, here's where the
%beauty of our rate constants come in.  Set a time step to be some small
%number and call it dt.  If at each time point the chance of a protein
%coming onto an unoccupied spot is equal to dt*(k_on * L), then the average
%rate of change over a larger period of time is equal to exactly what we
%asked it to be in (2).  How do we get that randomness?  rand<P!
t = 0; %make a time variable
dt = .0001; %very small
k_on = 1;
L = 1;
boundAtSpot = 0; %imagine we're focused on just one binding location
while ~boundAtSpot %going to repeat this until something binds
    if rand<(L*k_on)*dt %the important line
        boundAtSpot = 1;
    end
    t = t+dt;
end

%Note: dt has to be small, such that (L*k_on)*dt << 1.  Why?  This is where
%the discrete nature of the problem comes in.  If dt is too big, our
%approximate breaks down.  I don't remember the proof for that, but it's
%obvious that something's not working right if (L*k_on)*dt is ever bigger
%than 1.

%Anyway, to prove that the above is generating the behavior we want, I'm
%going to run it a bunch of times to see how long, in general, it takes our
%hypothetical molecule to bind to that one spot:


for j = 1:10000
    boundAtSpot = 0; %reset the stuff
    t = 0;
    while ~boundAtSpot %going to repeat this until something binds
        if rand<(L*k_on)*dt
            boundAtSpot = 1;
        end
        t = t+dt;
    end
    tListOn(j) = t;
end

%Okay, let's see what we got.

histogram(tListOn,100); 
xlabel('Time to bind (s)')
meanOnTime = mean(tListOn)
stdOnTime = std(tListOn)

%This is exactly the behavior we expect.  Stochastic events like this
%always have an exponential distribution of time lengths.  It happens in
%radioactivity, kinetics, earthquakes, whatever.  Look up the wikipedia
%page on the exponential distribution for more examples and insights.  The
%important thing is that we see that the mean time to bind on is 1 second.
%This implies an average rate of 1 per second.  And this is exactly what we
%expect because our equation (per binding location) was dA/dt = L * k_on.

%In the above code, change the k_on to 2, or the L to 1/2, and see what
%happens to the time it takes to bind.

%That's the mechanics of making something bind.  But once they've bound, we
%must give each molecule the opportunity to leave.  This is in accordance
%with the equation

%   dA/dt = - k_off * A

%   The number which leave must be proportional to the number which are on
%   the lattice.  How do we do this?  Iterate over everybody again!

boundLocations = [5, 17, 28, 34] %four proteins are bound right now
for i = boundLocations
    if rand<k_off*dt
        %remove the protein from the lattice
        %update variables
    end
end

%(Doing the unbinding in a random order is not necessary because the
%events are independent of each other; the binding events are not
%independent, because if something binds at position 1, binding at position
%2 is no longer possible if n>1)

%This strategy again enusres that we get an off rate proportional to k_off
%and to the total number of A.  It will result in that nice exponential
%distribution, and as long as dt is small it will be a very good
%approximation.

%So a few tips:

%dt has to be small for another reason: we don't want too many proteins to
%be leaving/entering on any given time step, since that could also bias the
%final result.  This just means that, ideally, both A*k_off and
%L*k_on*(number of available locations) never goes too far above 1.  But at
%some point of making dt super small your code is just going to be taking a
%long time to run, and if that happens talk to me and we'll see if there's 
%a good way to speed it up

%But until your code is taking more than a minute to run, don't bother
%making it faster.  The rule is, don't spend more than a minute to improve
%your code's performance by a few seconds.  Unless your code is actually
%taking multiple minutes (or hours) to run, the time spent making it
%"faster" by removing if statements or making things tighter is always
%going to be longer than the time saved during run.

%That said, if your code is taking an hour, there's some paralellization
%and vectorization stuff I can go over to help it be a bit faster.


%These techniques are used, in some variations, throughout the world of
%computational physics.  Understanding the mechanisms behind this stuff is
%a cool skill.

%I've probably missed something.  Try things and do some reading as well.
%If you're stuck, ask me, or if you just want to run things by someone.

