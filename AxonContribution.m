colors = cbrewer('qual', 'Set1', 10);
startup;

% DBX - red
% ALX - blue
% Barhl - green

if ismac
    addpath(genpath('/Users/ashwin/Documents/LabWork'));
    df = readtable('/Users/ashwin/Documents/SynapseDetector/09202018.csv');
else
    addpath(genpath('/usr/people/ashwinv/seungmount/research/Ashwin/Scripts'));
    df = readtable('/usr/people/ashwinv/seungmount/research/Ashwin/SynapseDetector/09202018.csv');
end

%load('FiringRates.mat');
load('TAU.mat');
load('STA.mat');
load('AllCells.mat');
Firing = A;

DbxCells = [76182,76183,76185,76186,76188,76189,76191,76199,76200];
OriginalCellOrderDBX = [8,9,10,11,21,12,15,2,3];
DbxTimeConstants = TAU(OriginalCellOrderDBX);

t = [-2:0.05:7];


%%

corrAndSizeDBX = [];
for i = 1:length(DbxCells)
    for j = 1:length(DbxCells)
        if i == j
            corrAndSizeDBX = [corrAndSizeDBX];
        else
            [preSynapticInputs1,preSynapticInputs1PSD] = SynapticPartners(DbxCells(i),1,df);
            [preSynapticInputs2,preSynapticInputs2PSD] = SynapticPartners(DbxCells(j),1,df);
            [commonInputs,loc1,loc2] = intersect(preSynapticInputs1, preSynapticInputs2);
            psdSize1 = df.size(preSynapticInputs1PSD(loc1));
            psdSize2 = df.size(preSynapticInputs2PSD(loc2));
            corr = corrcoef(Firing(:,OriginalCellOrderDBX(i)), Firing(:,OriginalCellOrderDBX(j)));
            corrAndSizeDBX = [corrAndSizeDBX; commonInputs, corr(1,2).*ones(size(psdSize1,1),1), psdSize1, psdSize2];
            
            clear {'preSynapticInputs1','preSynapticInputs1PSD', 'psdSize1','loc1'};
            clear {'preSynapticInputs2', 'preSynapticInputs2PSD','psdSize2', 'loc2'};
            clear commonInputs;
        end
    end
end

%%

uniqueCommonInputs = unique(corrAndSizeDBX(:,1));
index  = 1;

for i =1:size(uniqueCommonInputs)
    [A,B] = SynapticPartners(uniqueCommonInputs(i),2,df);
    [lia,lib] = ismember(DbxCells,A);
    if sum(lia)>1
        averagedFiring(index,:) = mean(Firing(:,OriginalCellOrderDBX(find(lia))),2);
        index = index+1;
        axonID(index)  = uniqueCommonInputs(i);
    end
end


averageFiringDBX = mean(Firing(:,OriginalCellOrderDBX),2);

for i = 1:size(axonID,2)
    subplot(12,12,i)
   % plot(t,(averagedFiring(i,:)'- averageFiringDBX));
      plot(t,(averagedFiring(i,:)));
      hold on;
      plot(t,averageFiringDBX);
   % DiffFromAverage(i,:) = (averagedFiring(i,:)'- averageFiringDBX);
    line([-2,7], [0,0],'color','k');
    %plot(t,mean(averagedFiring(:,:)),'k');
   % plot(t,mean(Firing(:,OriginalCellOrderDBX),2),'r');
   % YLim = [0,1];
    XLim = [-2,7];
    box off;
    %axis square;
    axis off;
    title(axonID(i));
end

%% sort the traces

[coeffs, score] = pca(DiffFromAverage);
%plot3(score(:,1),score(:,2),score(:,3),'o')
%box on;
% for i = 1:size(axonID,2)
%     text(score(i,1)+0.01,score(i,2)+0.01,score(i,3)+0.01, num2str(axonID(i)));
% end

figure;
 [idx,c] = kmeans(averagedFiring,2,'Distance','cityblock');
 
 temp = find(idx == 1);
 
%  for i = 1:length(temp)
%      plot(t,DiffFromAverage(temp(i),:),'Color',[0,0,1,0.8]);
%      hold on;
%  end
 hold on;
 plot(t, mean(DiffFromAverage(temp,:)),'b');

 clear temp
 
 temp = find(idx == 2);
%   for i = 1:length(temp)
%      plot(t,DiffFromAverage(temp(i),:),'r');
%      hold on;
%   end
  
 plot(t, mean(DiffFromAverage(temp,:)),'r');

clear temp


%%

[idx1,c1] = kmeans(Firing(:,OriginalCellOrderDBX)',2);

