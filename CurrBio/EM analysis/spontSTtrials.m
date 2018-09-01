function STA=spontSTtrials(j,t,cl,e,eyeThr,fixLength)
global SPT delt
delt = 0.05;
%% Saccade triggered averages of the spontaneous data
E=[e;-e];%loop through twice, but with the eye position flipped the second time.
%This has the effect of finding cells with positive eye correlation first and
%then finding those with negative correlation on the second iteration
% if max(max(cl))>10;
%     cl=cl./repmat(mean(cl),length(cl),1);cl=cl-repmat(mean(cl),length(cl),1);%normalize fluorescence
% end
% x=linspace(-1,1,length(cl))';xmat=[ones(size(x)) x];clear x;cl=cl-xmat*pinv(xmat)*cl;%subtract the baseline
% zer=sort(cl);cl=cl-repmat(mean(zer(1:round(length(zer)*.1),:)),length(cl),1);%find lowest 20% of fluorescence points (baseline)

%%%%%%%%%%%%use this to subtract low freq baseline%%%%%%%%%%%%%%%%%%%%%
% per=.2;
% minutes=60/delt;
% for ii=1:size(cl,2);
%     tt=[];y=[];    
%             for i=1:floor(t(end)/60);
%                 a=(i-1)*minutes+1;
%                 a=a:a+minutes;a(a>length(t))=[];
%                 a=round(a);
%                 T=t(a);c=cl(a,ii);
%                 [a,b]=sort(c);
%                 ind=sort(b(1:floor(length(b)*per)));
%                 y=[y;c(ind)];
%                 tt=[tt T(ind)];
%             end
%             x=tt'-mean(tt);x=x/max(x);
%             %         x=[x'.^2 x' ones(size(x'))];
%             T=t-mean(t);T=T/max(T);
%             G=cl(:,ii);
%             A=[];
% %             FUN=inline('p(1)*t+p(2)+p(3)*(t-p(4)).^2','p','t');
% %             b = LSQCURVEFIT(FUN,[1 1 1 0],x,y,...
% %                 [-100 -100 -100 -.3333],[100 100 100 .3333]);
%             FUN=inline('p(1)*t+p(2)+p(3)*(t).^2+p(4)*t.^3','p','t');
%             b = LSQCURVEFIT(FUN,[1 1 1 1],x,y);
%             g=FUN(b,T);A=[A g'];
%             cl(:,ii)=cl(:,ii)-g';
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%This is an alternative method suggested by reviewer 2%%%
% for ii=1:size(cl,2)
%     [b,a]=butter(2,delt/t(end)*2,'high');
%     y=filter(b,a,cl(:,ii));
%     k=sort(y);mn=mean(k(1:floor(length(k)*.2)));
%     cl(:,ii)=y-mn;
% end
%%%%%%%%%%%%%%%%%%%%%



c=corrcoef([e' cl]);cor=c(2:end,1);%find correlation between all cells and the eye
%     f=figure;subplot(1,5,1:4);stack_plot(t/.05,cl);set(gcf,'position',[360 93 500 829]);
%     ex=input('exclude?');subplot(155);plot(t/.05,e);close(f);
ex=[];
fix=round(fixLength/delt);pre=round(2/delt);len=pre+fix+1;T=linspace(-pre*delt,delt*fix,len);%fixation length,length before fixation,time of fixation
staf=zeros(length(T),size(cl,2));
dosu=zeros(size(staf,1),size(staf,2),20);
for W=1;
    e=E(W,:);%f=figure;plot(e);ex=input('exclude?');close(f);
    [c,pp]=corrcoef([e' cl]);COR=c(2:end,1);pval=pp(2:end,1); %find correlation between all cells and the eye
    COR=ones(size(COR));
    if length(find(COR>0))>0
        CL=cl(:,COR>0);%CL=[zeros(pre,size(CL,2));CL;zeros(fix,size(CL,2))];
        vel=[0 diff(medfilt1(e,1+round(1/delt)))]/delt;%calculate eye velocity
        k=sort(e);R=mean(k(1:round(length(k)*.4)));%find the bottom 40% of eye positions
        sacc=find(abs(vel)>50);aa=find(diff(sacc)<fix+floor(1/delt));sacc(aa)=[];
        sacc(sacc+fix>length(t))=[];%keep saccades that are separated by at least 6s.
        sacc([1 diff(vel(sacc)>0)]==0 & vel(sacc)<0)=[];sacL=sacc(vel(sacc)>0);sacR=sacc(vel(sacc)<0);%only keep contralateral saccades that are preceded by an ipsilateral saccade
        sacL=repmat(sacL,len,1)+repmat((-pre:fix)',1,length(sacL));%create a matrix of fixation indices, columns correspond to individual fixations and rows are the time points
        sacR=repmat(sacR,len,1)+repmat((-pre:fix)',1,length(sacR));
        ll=ismember(sacL,ex);rr=ismember(sacR,ex);sacL(:,any(ll)==1)=[];sacR(:,any(rr)==1)=[];
        LL=reshape(sacL,numel(sacL),1);%these are all of the data points for sta's
        RR=reshape(sacR,numel(sacR),1);%off sta's
        LL(LL<1)=1;RR(RR<1)=1;
        
        %editing
%         eyeThr=0;
        Lpos=max(reshape(e(LL),len,size(sacL,2)));%percent=abs((mean(Lpos)-Lpos)./(mean(Lpos)-R));%saccade magnitudes
        sacL(:,Lpos<eyeThr(1) | Lpos>eyeThr(2))=[];%get rid of leftward saccades whose amplitude is >20% smaller than average
        EE(W)=mean(Lpos(Lpos>eyeThr(1) & Lpos<eyeThr(2)));
        
        %
        
        RR=reshape(sacR,numel(sacR),1);LL=reshape(sacL,numel(sacL),1);RR(RR<1)=1;LL(LL<1)=1;RR(RR<1)=1;
        upL=reshape(CL(LL,:),len,size(sacL,2),size(CL,2));%all on saccades
        tempEye=repmat(e',1,size(CL,2));
        eyL=reshape(tempEye(LL,:),len,size(sacL,2),size(CL,2));
        dnL=reshape(CL(RR,:),len,size(sacR,2),size(CL,2));%all off saccades
        %find the saccades following off saccades, find the avg fluor immediately
        %before the cell turns back on, call this zero. I am setting zero
        %separately for each fixation because occasionally there is some drift in
        %the fluorescence so it seems more appropriate to set local rather than
        %global zeros.
        [a,b]=intersect(sacc,sacR(1,:)+pre);b=b+1;b(b>length(sacc))=length(sacc);
        sacRz=sacc(b);p3=round(3/delt);p1=round(1/delt);p05=round(.5/delt);
        sacRz=repmat(sacRz,pre,1)+repmat((-p3:-p1-1)',1,length(sacRz));LL=reshape(sacRz,numel(sacRz),1);
        dnz=repmat(mean(reshape(CL(LL,:),pre,size(sacRz,2),size(CL,2)),1),len,1);
        %         dnz=dnz*0;
        stadL=squeeze(mean(dnL-dnz,2));stad(:,COR>0)=stadL;
        clear tauL RL conInt;
        %%this is the old cirf fit, doesn't have conf intervals
        for i=1:size(stadL,2);
            k=ezfit(T(pre+1:end),stadL(pre+1:end,i),'a*exp(-x/t)');
            tauL(i)=k.m(2);RL(i)=k.r;
        end;
%         tauL=ones(1,size(stadL,2));RL=tauL;
        %%%
        %%this is the new cirf fit, does have conf intervals
%         for i=1:size(stadL,2);
%             decay=['p(1)*exp(-t/p(2))'];
%             ConExp=inline(decay,'p','t');
%             try
%             [tauf,resid,Jacobian,sigma] = NLINFIT(T(pre+1:end)',stadL(pre+1:end,i)...
%                 ,ConExp,[.1 1]);
%             CI = nlparci(tauf,resid,'covar',sigma,'alpha',.1);
%             conInt(i)=CI(2,2)-tauf(2);
%             fittt(:,1)=ConExp(tauf,T(pre+1:end));
%             RL(i)=1-sum((stadL(pre+1:end,i)-fittt(:,1)).^2)./...
%                 sum((stadL(pre+1:end,i)-mean(stadL(pre+1:end,i))).^2);
%             tauL(i)=tauf(2);
%             catch
%                 RL(i)=nan;tauL(i)=nan;conInt(i)=nan;
%             end
%         end;
%         %%%%%%%%
        TAU(COR>0)=tauL;cirfR(COR>0)=RL;%cirfCI(COR>0)=conInt;
        for qq=1:20;
        if size(dnL,2)<2;
            dosR=NaN(size(dnL,1),size(dnL,3));%dnL=NaN(size(dnL));
        else;
            ind=randperm(size(dnL,2));ind=ind(1:2);
            dosR(:,:,qq)=1/sqrt(2*size(dnL,2))*...
                squeeze(dnL(:,ind(1),:)-dnL(:,ind(2),:));
        end
        end
        stadL=squeeze(mean(dnL,2));stad(:,COR>0)=stadL;dosd(:,COR>0,:)=dosR;        
        %find the fluorescence immediately before the cell turns on and subtract.
        %For on saccades that are not preceded by off saccades, set zero as the
        %average fluorescence before the on saccade immediately following the most
        %recent off saccade.
        [a,b]=intersect(sacL(1,:)+pre,sacRz(1,:)+p3);sacLz=sacc(b);
        %         upz(:,b,:)=repmat(mean(upL(1:1+pre,b,:),1),len,1);b=setdiff(1:size(upL,2),b);b(b==1)=3;
        %         upz(:,b,:)=repmat(mean(upL(1:1+pre,b-1,:),1),len,1);
        upz=upL*0;
        upL=upL-upz;
        %subtract the burst so that all on saccades start with fluorescence=0, this allows us
        %to average saccades that start from different positions (i.e. fluorescence is different,
        %but underlying firing rate is the same)
        TT=repmat(T'-T(pre+p05+1),[1,size(upL,2),size(upL,3)]);TT(TT<0)=inf;
        b=repmat(upL(pre+p05+1,:,:),size(upL,1),1).*...
            exp(-TT./reshape(repmat(tauL,[size(upL,1)*size(upL,2),1]),[size(upL)]));
        b=0*b;
        upL=upL-b;
        for qq=1:20;
        if size(upL,2)<2;
            dosL(:,:,qq)=NaN(size(upL,1),size(upL,3));%upL=NaN(size(upL));
        else;
            ind=randperm(size(upL,2));ind=ind(1:2);
            dosL(:,:,qq)=1/sqrt(2*size(upL,2))*...
                squeeze(upL(:,ind(1),:)-upL(:,ind(2),:));
        end
        end
        staf(:,COR>0)=squeeze(mean(upL,2));
        err(:,COR>0)=squeeze(std(upL,0,2))/sqrt(size(upL,2));
        dosu(:,COR>0,:)=dosL;
        stf(W).all_trials=zeros(size(upL,1),size(upL,2),size(cl,2));
        stf(W).all_trials(:,:,COR>0)=upL;
        staE(:,COR>0)=squeeze(mean(eyL,2));
        Ts=T;clear upz dosL dosR
    end;
end
% sp=fieldnames(SPT);
STA.staf=staf;STA.staE=staE;STA.stad=stad;
STA.err=err;
STA.stf=stf;
STA.dosu=dosu;STA.dosd=dosd;
STA.intcor=cor;STA.pval=pval;STA.T=Ts;
STA.cirfTau=TAU;STA.cirfRsq=cirfR;%STA.cirfCI=cirfCI;
[c,pp]=corrcoef([e' cl]);
STA.cor=(c(2:end,1))';