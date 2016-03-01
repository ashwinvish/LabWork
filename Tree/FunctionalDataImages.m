% LM functional plots

load('101112 _files1_4.mat')
FLUOR(1).ROI=SPT;
FLUOR=stimSummaryEM(FLUOR);
STA=FLUOR(1).STA;
SPT=FLUOR(1).ROI;
close all;

for i=2:length(SPT); % starts from 2 since we do not have EM for 1st of length(SPT) planes.
    % plot the imaging plane and ROIs
    ROI=SPT(i);
    sta=STA(i).staf; % spike triggered average fluorescence
    stae=STA(i).staE; % spike triggered average error
    rsqs=SPT(i).rsqs;
    int=SPT(i).intCorrect;
    err=STA(i).err;
    cor=max(abs(SPT(i).cor));
    stc=diag(corr(sta(1:80,:),stae(1:80,:)))';
    cls=find(stc>.6); % coorelation threshold for staf with stae
    cls=intersect(cls,find(rsqs>.25)); % CIRF as reported in Miri 2011
   cls=intersect(cls,find(cor>.2)); % threshold for correlation of fl and eye pos.       
    if i==2;
        cls=[cls 17 9];
    end
    areaI=SPT(i).areaI;
    figure();
    im=ROI.mnImage; % LM image
    ct=ROI.centroid; % centroid of cells of ineterst
    ct=ct(:,cls)';
    im=repmat(ROI.mnImage,[1 1 3]);
    im=im/max(im(:));
    mx=find(im>max(im(:))*.7); 
    mn=find(im<max(im(:))*.0);
    im(mn)=0;im(mx)=.8;
    image(im+.1);
    str=['image',num2str(i)];
    figure();
    ROIcellLabel(ROI,cls,1,[1 0 0]);
    set(gcf,'units','inches','position',[3 3 4 4]);
    set(gcf,'paperposition',[3 3 4 4]);
    
    % plot spike triggered average fluorescence with errors
    sc=.75;
    n=length(cls);
    t=1:size(sta,1);

    figure();
    % eye pos
    subplot(n+1,1,1);
    plot(stae(:,1),'b','linewidth',2);
    xlim([0 t(end)+40]);
    set(gca,'color','none');
    str=['sta',num2str(i)];
    pbaspect([1,1,1]);
    
    for j=1:n;
        a=sta(:,cls(j)); % get sta for the cells in cls
        e=err(:,cls(j)); % get error for cells in cls
        l=max(abs(e));     
        subplot(n+1,1,j+1);
        confplot(t,a,e);hold on; % continious error boundaries
        plot(t,0*t,'k--'); % dashed line at eye pos 0;
        set(gca,'color','none');
        ylim([-.05 l+max(a)]);
        g=text(t(end)+10,a(end)+.05,num2str(j));
        xlim([1 t(end)+40])
        set(g,'color','r','fontsize',14,'fontweight','bold');
        pbaspect([1,1,1]);
    end

    
    n=length(cls);
    figure();
    
    subplot(n+1,1,1);
    th=SPT(i).theta;
    plot(th(:,1),medfilt1(th(:,2),11),'b','linewidth',4);
    hold on;
    xlim([0 th(end,1)+30]);
    ylim([min(medfilt1(th(:,2),11)) max(medfilt1(th(:,2),11))]);
    set(gca,'color','none');
    str=['timeseries',num2str(i)];
    pbaspect([10,1,1]);
    
    for j=1:n;
        t=SPT(i).time(:,cls(j)); % time
        a=int(:,cls(j)); % deltaf/f
        a = a/max(a);
        subplot(n+1,1,j+1);
        plot(t,a,'k','linewidth',4);
        set(gca,'color','none');
        %ylim([min(a/max(a)) 1]);
        g=text(t(end)+10,a(end)+.05,num2str(j));
        xlim([1 t(end)+30])
        set(g,'color','r','fontsize',14,'fontweight','bold');
        pbaspect([10,1,1]);
    end

end