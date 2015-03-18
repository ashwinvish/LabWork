% 100ms acquisition with 20ms steps of varying amplitude
% Script to acquire data for short periods of time



%stimulation=zeros(1000,1);
%%
time_step = 100; % ms , time of acquisition
base = 20; %ms , length of latency before acquisition
stim = 50; %ms, length of stimulus
Amp = 100; % mV or nA , depending on IC or VC
freq = 10000; % Hz, freq of data acquisition
rem = time_step-(base+stim); % ms, remaining time


%%
% for single pulse stimulus

stimulation = Amp*[zeros(base*10,1); ones(stim*10,1); zeros(rem*10,1)];
outputdata = aoi(stimulation,[0],freq]; % sturcture [stimulus, recording channel(0,1,2..), freq of data acq]
figure;
subplot(211), plot(outputdata);
subplot(212), plot(stimulation);

%% for single stimulus with varying amplitude and plotting IV curve


Amp=[-60:10:60]; % 10 mv of nA steps
maxim = zeros(length(Amp),1)
figure;
stimulation = [zeros(base*10,1); ones(stim*10,1);zeros(rem*10,1)];

for i = 1: length(Amp)
temp_input=Amp(i)*stimulation;
outputdata=aoi(temp_input,[0],freq);
maxim(i) = find(max(-1*outputdata(base*10:(base+stim)*10))); % determining the peak of the output

hold on;
subplot(211), plot(temp_input);
subplot(212),plot(outputdata)
hold all;
end;

figure;
plot(Amp,maxim);

%%
%outputdata=aoi(stimulation,[0],10000);

%% sample for minis over 10 sec

stimulation=[zeros(100000,1)];
outputdata=aoi(stimulation,[0],10000);
plot(outputdata);
