clear all; %clear all variables
close all; %close all open figure windows
%DEFINE NETWORK PARAMETERS
NumPatt = 1 %Number of memory patterns that network will attempt to store and retrieve
N = 100 %Number of neurons in the network (and number of elements in a memory pattern)
NumTimeSteps = 20 %length of run
%MEMORY PATTERNS
Mem_mat = 2*round(rand(N,NumPatt))-1; %random strings  of 1's and -1's
%DEFINE SYNAPTIC WEIGHT MATRIX
W_mat = Mem_mat*Mem_mat';
%INITIALIZE NETWORK ACTIVITY
qStart = 0.7; %average overlap of initial activity pattern with first memory (the one to be retrieved)
for i=1:N
 if rand < qStart %assign x_i equal to Mem(i,1)
 x_vect(i)=Mem_mat(i,1);
 else %assign random value for x_i
 if rand < .5
 x_vect(i) = 1;
 else
 x_vect(i) = -1;
 end
 end
end
%compute initial overlap
q_vect = zeros(1,NumTimeSteps);
q_vect(1) = x_vect*Mem_mat(:,1)/N; 

%RUN MODEL AND COMPUTE OVERLAP AT EACH TIME STEP
for tStep=2:NumTimeSteps
 x_vect = (2*(W_mat*x_vect' >= 0) - 1)'; %update rule
 q_vect(tStep) = x_vect*Mem_mat(:,1)/N; %compute overlap of newly updated x_vect with Mem1
end

figure(1)
subplot(3,1,1)
plot(0:NumTimeSteps-1,q_vect,'-o')
xlabel('Time step')
ylabel('Overlap of x(t) with Mem^{(1)}')
subplot(3,1,2)
image(50*Mem_mat(:,1)')
xlabel('Memory pattern element')
ylabel('blue=-1,red=+1')
subplot(3,1,3)
image(50*x_vect)
xlabel(strcat('x(t=',num2str(NumTimeSteps-1),')'))
ylabel('blue=-1,red=+1')