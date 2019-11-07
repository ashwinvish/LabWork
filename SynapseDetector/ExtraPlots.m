load scores.mat

colorSchemes;

figure;
subplot(4,4,1);
histogram(1-mnx,0:0.1:1,'Normalization','probability','FaceColor',ABDcolor,'FaceAlpha',0.2,'EdgeColor','none');
hold on;
histogram(1-diO,0:0.1:1,'Normalization','probability','FaceColor','r','FaceAlpha',0.2,'EdgeColor','none');
yyaxis ('right');
histogram(1-mnx,0:0.05:1,'DisplayStyle','stairs','Normalization','cdf','LineWidth',2,'EdgeColor',ABDcolor,'LineStyle',':');
histogram(1-diO,0:0.05:1,'DisplayStyle','stairs','Normalization','cdf','LineWidth',2,'EdgeColor','r','LineStyle',':');
%axis square;
box off;
legend({'mnx','diO'});
%offsetAxes(gca);

%% Eye traces

addpath(genpath('/Users/ashwin/Documents/LabWork'));
load viewRightDirSTATraces2-Tonic.mat
startup

figure;
for i = 1:7
subplot(7,1,i);
plot(TcellBB{i}{1}(:,1), EcellBB{i}{1}(:,1),'color','b','LineWidth',2); % eye position traces in degrees
hold on;
plot(TcellBB{i}{1}(:,2), EcellBB{i}{1}(:,2),'color','r','LineWidth',2);
yyaxis('right');
plot(TFBB{i}{1},FBB{i}{1},'color','k','LineWidth',2) % fluorescence in dF/F 
box off;
end


figure;
for i = 1:7
subplot(7,1,i);
[~,a] = min(abs(TcellBB{i}{1}(:,1)-50));
[~,b] = min(abs(TcellBB{i}{1}(:,2)-50));
[~,c] = min(abs(TcellBB{i}{1}(:,1)-100));
[~,d] = min(abs(TcellBB{i}{1}(:,2)-100));
plot(TcellBB{i}{1}(a:c,1), EcellBB{i}{1}(a:c,1),'color','b','LineWidth',2); % eye position traces in degrees
hold on;
plot(TcellBB{i}{1}(b:d,2), EcellBB{i}{1}(b:d,2),'color','r','LineWidth',2);
yyaxis('right')
plot(TFBB{i}{1}(50:100,1),FBB{i}{1}(50:100,1),'color','k','LineWidth',2) % fluorescence in dF/F 
box off;
end

% figure;
% plot(TcellBB{2}{1},EcellBB{2}{1}) % fluorescence in dF/F 
% yyaxis('right');
% plot(TFBB{2}{1},FBB{2}{1}) % fluorescence in dF/F \
% 
% 
% figure;
% plot(TcellBB{2}{1},EcellBB{2}{1}) % fluorescence in dF/F 
% yyaxis('right');
% plot(TFBB{2}{1},FBB{2}{1}) % fluorescence in dF/F \
% 
% 
% figure;
% plot(TcellBB{3}{1},EcellBB{3}{1}) % fluorescence in dF/F 
% yyaxis('right');
% plot(TFBB{3}{1},FBB{3}{1}) % fluorescence in dF/F \
% 
% 
% figure;
% plot(TcellBB{3}{2},EcellBB{3}{2}) % fluorescence in dF/F 
% yyaxis('right');
% plot(TFBB{3}{2},FBB{3}{2}) % fluorescence in dF/F 


%% Velocity traces

load viewRightDirSTATraces-Burst.mat

figure;
for i = 1:4
    
subplot(4,1,i)
[~,a] = min(abs(TcellBB{i}{1}(:,1)-50));
[~,b] = min(abs(TcellBB{i}{1}(:,2)-50));
[~,c] = min(abs(TcellBB{i}{1}(:,1)-200));
[~,d] = min(abs(TcellBB{i}{1}(:,2)-200));
plot(TcellBB{i}{1}(a:c,1), EcellBB{i}{1}(a:c,1),'color','b','LineWidth',2); % eye position traces in degrees
hold on;
plot(TcellBB{i}{1}(b:d,2), EcellBB{i}{1}(b:d,2),'color','r','LineWidth',2);
yyaxis('right')
plot(TFBB{i}{1}(50:200,1),FBB{i}{1}(50:200,1),'color','k','LineWidth',2) % fluorescence in dF/F 
box off;
% 
%     plot(TcellBB{i}{1}(:,1), EcellBB{i}{1}(:,1),'color','b','LineWidth',2); % eye position traces in degrees
% hold on
% plot(TcellBB{i}{1}(:,2), EcellBB{i}{1}(:,2),'color','r','LineWidth',2); % eye position traces in degrees
% 
% %plot(TcellBB{1}{1},EcellBB{1}{1}) % eye position traces in degrees
% yyaxis('right');
% plot(TFBB{1}{1},FBB{1}{1}) % fluorescence in dF/F 
end



%% 

colorSchemes

load ABDPutativeSaccadic.mat
load ABDiPutativeSaccadic.mat


figure;
subplot(4,4,1)
errorbar(ABDPutativeSaccadic.meanABDgradient,ABDPutativeSaccadic.stdABDgradient./sqrt(29),...
    '-o','color',ABDcolor,'LineWidth',2,'MarkerFaceColor','w');
hold on;
errorbar(ABDiPutativeSaccadic.meanABDigradient,ABDiPutativeSaccadic.stdABDigradient./sqrt(29),...
    '-o','color',ABDicolor,'LineWidth',2,'MarkerFaceColor','w');
axis square;
set(gca,'XTickLabels',[0,0.5,1]);
xlabel('Norm. pathlength');
ylabel('Norm. count');
box off;
offsetAxes(gca);

subplot(4,4,2)

temp = isMotor(ABDPutativeSaccadic.cellIDs,df);
temp2 = isMotor(ABDiPutativeSaccadic.cellIDs,df);
plot(sum(temp(:,2:3),2))




[a,b] = getABDgradient('ABDr',ABDPutativeSaccadic.cellIDs',true,false);
[c,d] = getABDgradient('ABDc',ABDPutativeSaccadic.cellIDs',true,false);

[e,f] = getABDgradient(ABDIr,ABDiPutativeSaccadic.cellIDs',true,false);
[g,h] = getABDgradient(ABDIc,ABDiPutativeSaccadic.cellIDs',true,false);
