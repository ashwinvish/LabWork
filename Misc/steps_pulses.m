% 100ms acquisition with 20ms steps of varying amplitude

%stimulation=zeros(1000,1);
%%
amplitude=[10:10:100];

for i = 1:10
stimulation=amplitude(i)*[zeros(1000,1); ones(20,1); zeros(100,1); ones(20,1);zeros(1000,1)];
outputdata=aoi(stimulation,[0],10000);
hold all;
subplot(211), plot(stimulation);
subplot(212),plot(outputdata)
%hold on;
end;
%%
%outputdata=aoi(stimulation,[0],10000);

%% sample for minis
stimulation=[zeros(500000,1)];
outputdata=aoi(stimulation,[0],10000);
plot(outputdata);
