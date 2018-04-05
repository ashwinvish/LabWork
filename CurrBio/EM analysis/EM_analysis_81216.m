clear
% addpath(genpath('C:\Users\Daie\Desktop\EM analysis'))
addpath(genpath('C:\Users\Daie\Documents\EM\EM analysis\EM analysis\'))
load('101112 _files1_4.mat')
FLUOR(1).ROI=SPT;
FLUOR=stimSummaryEM(FLUOR);
%%
delt=.05;
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
load('calibrate1');
cor=max(abs([SPT.cor]));
% a=find(er>3 & cor>.2);
ind=size(SPT(1).intensity,2);
% for i=2:length(SPT);
%     areaI=[areaI SPT(i).areaI+ind];
%     b=size(SPT(i).intensity,2);
%     ind=ind+b;
% end
x=[SPT.centroid];
z=[SPT.z];
int=[SPT.intCorrect];

a=find(stc>.6);
a=intersect(a,find(rsqs>.25));
a=intersect(a,find(cor>.2));
a=sort([a size(SPT(1).intensity,2)+[17 9]]);
% a=intersect(a,areaI);
a=intersect(a,find(z~=0));
sta=sta(22:end,a);
x=x(:,a);x/512*calibrate1(3);
x=[x;z(a)];
n=length(a);
int=int(:,a);
for i=1:n;
    for j=1:n;
        d(i,j)=sqrt(sum((x(:,i)-x(:,j)).^2));
    end
end

load('modes_stafs');
fs=[fs sta];
t=0:delt:delt*(size(fs,1)-1);t=t-2;
r=svd_rates(fs,t,1.9,[3,3],-100);
r=r(:,end-(n-1):end);
firing=r;fluorescence=sta(1:end-1,:);T=t(1:end-1);
r=r(40:end,:);
t=0:.05:.05*(size(r,1)-1);
pwCor=corr(int);
o=find(eye(n,n)==0);
% [Xs,Ys,stds,Cs,Ps]=mean_bin_plot(d(o),pwCor(o),5,1,1);
staE=[STA.staE];eyes=mean(staE,2);
%%
for i=1:n;k=ezfit(t,r(:,i),'a*exp(-x/t)');tau(i)=k.m(2);end
tau(tau>100)=100;tau(tau<1)=1;
[ds,ps]=PairwiseSpatialStructure2(x',tau',1);
[Xs,Ys,stds,Cs,Ps]=mean_bin_plot(ds,ps,5,1,1);
%%
% folder='C:\Users\kpd7\Documents\data\EM\101112\zfigures\';
folder = 'C:\Users\Daie\Documents\EM\figures\';
STA=FLUOR(1).STA;
SPT=FLUOR(1).ROI;
TAU = [];
close all;
for i=2:length(SPT);
    ROI=SPT(i);
    sta=STA(i).staf;stae=STA(i).staE;
    rsqs=SPT(i).rsqs;
    int=SPT(i).intCorrect;
    err=STA(i).err;
    cor=max(abs(SPT(i).cor));
    stc=diag(corr(sta(1:80,:),stae(1:80,:)))';
    cls=find(stc>.6);
    cls=intersect(cls,find(cor>.2));
    cls=intersect(cls,find(rsqs>.25));
    if i==2;cls=[cls 17 9];end
    areaI=SPT(i).areaI;
%     cls=intersect(cls,areaI);
%     cls=1:size(sta,2);
    figure(3*(i-1)+1);
    im=ROI.mnImage;
    ct=ROI.centroid;
    ct=ct(:,cls)';
    
    im=repmat(ROI.mnImage,[1 1 3]);
    im=im/max(im(:));
    mx=find(im>max(im(:))*.7);
    mn=find(im<max(im(:))*.0);
    im(mn)=0;im(mx)=.8;
    image(im+.1);
    str=['image',num2str(i)];
    hold on;
%     plot(ct(:,1),ct(:,2),'r.')
    ROIcellLabel(ROI,cls,1,[1 0 0])
    set(gcf,'units','inches','position',[3 3 4 4]);
    set(gcf,'paperposition',[3 3 4 4]);
        print(gcf,'-depsc',[folder,str,'.eps'])
    
    
    load('modes_stafs');
    t=0:delt:delt*(size(fs,1)-1);t=t-2;
    r=svd_rates([fs sta(22:end,:)],t,1.9,[3,3],-100);
    r=r(:,end-(size(sta,2)-1):end);
    firing=r;fluorescence=sta(1:end-1,:);T=t(1:end-1);
    r=r(40:end,:);
    t=0:.05:.05*(size(r,1)-1);
    
    clear tau
    for tau_i=1:size(sta,2);k=ezfit(t,r(:,tau_i),'a*exp(-x/t)');tau(tau_i)=k.m(2);end
    tau(tau>100)=100;tau(tau<1)=1;
    
    TAU = [TAU;tau(cls)'];
    
    sc=.75;
    figure(3*(i-1)+2);figure_initialize;
    set(gcf,'position',[4 .1 1.2 9]*sc,'paperposition',[4 .1 1.2 9]*sc);
    n=length(cls);t=1:size(sta,1);set(gcf,'color','w');
    for j=1:n;
        a=sta(:,cls(j));
        e=err(:,cls(j));
        l=max(abs(e));
        KDsubplot(n+1,1,[j+1,1],0);        
        confplot(t,a,e);hold on;
        plot(t,0*t,'k--');
        set(gca,'visible','off');
        ylim([-.05 l+max(a)]);
        g=text(t(end)+10,a(end)+.05,num2str(j));
        xlim([1 t(end)+40])
        set(g,'color','r','fontsize',14,'fontweight','bold');
        g = text(t(end) + 10,a(end) - .05,['\tau = ',num2str(round(tau(cls(j))))])
        set(g,'color','r','fontsize',10,'fontweight','bold');
    end
     KDsubplot(n+1,1,[1,1],0);
    plot(stae(:,1),'b','linewidth',2);
    xlim([0 t(end)+40]);
    set(gca,'visible','off');
    str=['sta',num2str(i)];
    print(gcf,'-depsc',[folder,str,'.eps'])
    
    figure(3*(i-1)+3);figure_initialize;
    set(gcf,'position',[4 .1 4 9]*sc,'paperposition',[4 .1 4 9]*sc);
    n=length(cls);set(gcf,'color','w');
    for j=1:n;
        t=SPT(i).time(:,cls(j));
        a=int(:,cls(j));
        KDsubplot(n+1,1,[j+1,1],0);
        plot(t,a,'k','linewidth',2);
        set(gca,'visible','off');
        ylim([min(a) max(a)]);
        g=text(t(end)+10,a(end)+.05,num2str(j));
        xlim([1 t(end)+30])
        set(g,'color','r','fontsize',14,'fontweight','bold');
    end
    KDsubplot(n+1,1,[1,1],0);
    th=SPT(i).theta;
    plot(th(:,1),medfilt1(th(:,2),11),'b','linewidth',2);
    xlim([0 th(end,1)+30]);
    set(gca,'visible','off');
    str=['timeseries',num2str(i)];
    print(gcf,'-depsc',[folder,str,'.eps'])
end
