load('101112 _files1_4.mat')
FLUOR(1).ROI=SPT;
FLUOR=stimSummaryEM(FLUOR);
%%
delt=.05; % global threshold
STA=FLUOR(1).STA;
SPT=FLUOR(1).ROI;
sta=[STA.staf];err=[STA.err];
er=(sta./err);er=mean(er);
rsqs=[SPT.rsqs];
% for i=1:size(sta,2);
%     confplot(sta(:,i),err(:,i));
%     title(num2str(er(i)));
%     pause;
% end
staE=[STA.staE];
stc=diag(corr(sta(1:80,:),staE(1:80,:)))';
cor=max(abs([SPT.cor]));

[Y,I] =  sort(stc);
I = I(I>50);
noOfCells = 143-50; % removing cells from plane1 that was not in the EM volume
a=find(stc>.6);
a=intersect(a,find(rsqs>.25));
a=intersect(a,find(cor>.2)); % a describes the conditions that were used to filter for VPNI cells
figure(1);
histogram(stc(I),20);

figure(2);
for i = 1:noOfCells
    subplot(10,10,i)
    h = plotyy(1:180,sta(1:180,I(i)), 1:180, staE(1:180,I(i)));
    %line([80, 80],[0, max(sta(1:180,I(i)))], 'color', 'k', 'LineStyle', '--');
    plotTitle = sprintf('%1.2f,%1.2f,%1.2f',stc(I(i)),rsqs(I(i)),cor(I(i)));
    %axis(h(1),'tight','off');
    %axis(h(2),'tight','off');
    if ~ismember(I(i),a)
        title(plotTitle,'color','k');
    else
        title(plotTitle,'color','r');
    end
end


