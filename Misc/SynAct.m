% synaptic activation parameters


function [INa, IK, IKt, Ileak, Iin, V, time, m, h, n, bk] = SynAct(Vstart,Tmax,Iamp)
clc;
% Iamp in nA

TauS = 1000; % in msec
TauU = 100; % in msec
C = 0.5666; % in nF
Q = 0.1885/1000 ; % in seconds
Vthre = -20; % in mV
GammaTh = 2; % in mV

% Leak current
Gleak = 0.1490; % in uS
Eleak = -65; %in mV

%Sodium Current
GNabar = 56.66; %us
VNa = 55 ;%in mV
TauH = 0.1; % in msec

% Potassium Current
GKbar = 39.66; % in uS
VK = -80; % in mV
TauN = 0.1; % in mSec

% Transient Potassium
GKtbar = 4.264 ; % in uS
TauB = 30; % in mSec

% Spiking Model Shriki-Hansld-Sompolinsky model
dt= .01;
T = Tmax;
N = floor(T/dt);
V = zeros(1,N);
INa = zeros(1,N);
IK = zeros (1,N);
IKt = zeros(1,N);
Ileak = zeros(1,N);
Inoise = 1000*[0.23+0.06*randn(1,N)];

slope = Iamp/(N/2); % slope of ramp
Iin = zeros(1,N);
for i =1:N
    if i<=N/2
        Iin(i) = slope*(i-1);
    else
        Iin(i) = (Iamp*2)-(slope*(i-1));
    end
    
end

%Iin = Iamp*[zeros(1,10) ones(1,(N-20)) zeros(1,10)];
%Iin = Iamp*ones(1,N);
m= zeros(1,N);
h = zeros(1,N);
n = zeros(1,N);
bk = zeros(1,N);

V(1) = Vstart ; % in mV

for i = 2:N
    %Activation Variables
    m(i) = AlphaM(V(i-1)) / (AlphaM(V(i-1))+BetaM(V(i-1))) ;
    h(i) = h(i-1) + dt*((1/TauH)*(AlphaH(V(i-1))*(1-h(i-1))-(BetaH(V(i-1))*h(i-1))));
    n(i) = n(i-1) + dt*((1/TauN)*(AlphaN(V(i-1))*(1-n(i-1))-(BetaN(V(i-1))*n(i-1))));
    bk(i) = bk(i-1) + dt*((1/TauB)*(bKt(V(i-1)-bk(i-1))));
    % Ionic Conductances
    GNa = GNabar*m(i)^3*h(i);
    GK = GKbar*n(i)^4;
    GKt = GKtbar*(aKt(V(i-1)))^3*bk(i);
    %Ionic Currents
    INa(i) = GNa*(V(i-1)-VNa);
    IK(i) = GK*(V(i-1)-VK);
    IKt(i) = GKt*(V(i-1)-VK);
    Ileak(i) = Gleak*(V(i-1)-Eleak);
    % Membrane Voltage
    V(i) = V(i-1)+dt*(1/C*(-Ileak(i)-INa(i)-IK(i)-IKt(i)+Iin(i)+Inoise(i)));
end

time = [0:dt:Tmax-dt];
subplot(4,1,1) , plot(time, Iin);
subplot(4,1,1), ylabel ('Input Current nA');
subplot (4,1,2), plot(time,V);
subplot (4,1,2), ylabel('Membrane Voltage mV');
subplot(4,1,3), plot(time, INa, 'r');
hold on
subplot(4,1,3), plot(time, IK, 'g');
subplot(4,1,3), plot(time, IKt, 'b');
subplot(4,1,3), plot(time, Ileak, 'c');
subplot(4,1,3), ylabel ('Ionic Current nA');
xlabel ('Time msec');
subplot(4,1,4), plot(time, Inoise);
subplot(4,1,4), ylabel('Inoise nA');;

%subplot(3,1,1), legend('INa, IK, IKt,ILeak');


    function a = AlphaM(V)
        a = ((V+30)*0.1) / (1-exp(-(V+30)*0.1));
    end

    function b = BetaM(V)
        b = 4*exp(-(V+55)/18);
    end

    function c = AlphaH(V)
        c = 0.07*exp(-(V+44)/20);
    end

    function d = BetaH(V)
        d = 1/(exp(-(V+14)*0.1)+1);
    end



    function e = AlphaN(V)
        e = (0.1*((V+34)*0.1)) / (1-exp(-(V+34)*0.1));
    end

    function f = BetaN(V)
        f = 0.125*exp(-(V+44)/80);
    end


    function g = aKt(V)
        g = 1/(exp(-(V+50)/20)+1);
    end

    function h = bKt(V)
        h = 1/(exp((V+80)/60)+1);
    end


end



%% Synaptic Activation Parameters

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

Xlabel('Presynaptic potential mV');
Ylabel('\sigma(V_j) Activation state');

%% Steady State Activation Curves

%Rf = -9:1:60; % Point of inflection
%Theta = 0.5:0.5:35;
r = 0:1:80;
Rf = 14;
Theta = [6,14,22,30];
clear i;

  
    for j= 1:length(Theta)
            Ainf(j) = 1/(1+exp(Rf/Theta(j)));
            Binf(j) = 1/(1-Ainf(j)); 
        for i=1:length(r)
            Temp(i) =  Binf(j)*(1/(1+exp((Rf-r(i))/Theta(j)))-Ainf(j));
        end
        Sinf(j,:) = Temp(:);
        plot(r,Temp);
        hold all;
    end
  




    