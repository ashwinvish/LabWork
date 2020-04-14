clc;
clear all; 
TauS = 1000; % in msec
TauU = 100; % in msec
C = 0.5666; % in nF
Q = 0.1885/1000 ; % in seconds
Vthre = -20; % in mV
GammaTh = 2; % in mV

Vj = -60:1:60; % in mV

for i = 1:length(Vj)
        Sigma(i) = 1/(1+exp((Vthre-Vj(i))/GammaTh));
end
%% Plot presynaptic activation
figure;
plot(Vj,Sigma);

xlabel('Presynaptic potential mV');
ylabel('\sigma(V_j) Activation state');

%% Steady State Activation Curves

%Rf = -9:1:60; % Point of inflection
%Theta = 0.5:0.5:35;
r = 0:1:80;
Rf = 14;
Theta = [6];
%Temp =  zeros(1:length(r));
clear i;
figure;

  
    for j= 1:length(Theta)
            Ainf(j) = 1/(1+exp(Rf/Theta(j)));
            Binf(j) = 1/(1-Ainf(j)); 
        for i=1:length(r)
            Temp(i) =  Binf(j)*(1/(1+exp((Rf-r(i))/Theta(j)))-Ainf(j));
        end
        %Sinf(j,:) = Temp(:);
        plot(r,Temp);
        hold all;
        clear Temp;
    end