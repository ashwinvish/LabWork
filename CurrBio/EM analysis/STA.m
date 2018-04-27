% spike triggered averages

time = t ;                  % time
delt = 0.05;                % interpolated time

% interpolated Calcium signals
for i = 1:size(Cinterp,2)
CinterpMat(i,:) = cell2mat(Cinterp(i));
end
CinterpMat = CinterpMat'; 

EyePos = e;                         % eye positons
E=[e;-e];                           %loop through twice, but with the eye position flipped the second time.
eyeThr=[3 140];                     % Eye thresholds
fixLength=7;                        % fixation length in seconds
k=sort(e);
R=mean(k(1:round(length(k)*.4)));   %find the bottom 40% of eye positions

%find correlation between all cells and the eye
[c,pp] =corrcoef([e' CinterpMat]);  % coorelation and Pvalues
cor=c(2:end,1);                     
pval=pp(2:end,1);

ex=[];
fix=round(fixLength/delt); 
pre=round(2/delt);
len=pre+fix+1;
T=linspace(-pre*delt,delt*fix,len); %fixation length,length before fixation,time of fixation
staf=zeros(length(T),size(Cinterp,2));
dosu=zeros(size(staf,1),size(staf,2),20);

% eye velocity 
vel=[0 diff(medfilt1(e,1+round(1/delt)))]/delt;                            %calculate eye velocity

% saccades
sacc=find(abs(vel)>50);aa=find(diff(sacc)<fix+floor(1/delt));sacc(aa)=[];  % find saccades
sacc(sacc+fix>length(t))=[];                                               %keep saccades that are separated by at least 6s.
sacc([1 diff(vel(sacc)>0)]==0 & vel(sacc)<0)=[];
sacL=sacc(vel(sacc)>0);
sacR=sacc(vel(sacc)<0);                                                    %only keep contralateral saccades that are preceded by an ipsilateral saccade
sacL=repmat(sacL,len,1)+repmat((-pre:fix)',1,length(sacL));                %create a matrix of fixation indices, columns correspond to individual fixations and rows are the time points
sacR=repmat(sacR,len,1)+repmat((-pre:fix)',1,length(sacR));
ll=ismember(sacL,ex);rr=ismember(sacR,ex);
sacL(:,any(ll)==1)=[];sacR(:,any(rr)==1)=[];
LL=reshape(sacL,numel(sacL),1);                                            %these are all of the data points for sta's
RR=reshape(sacR,numel(sacR),1);                                            %off sta's
LL(LL<1)=1;RR(RR<1)=1;

Lpos=max(reshape(e(LL),len,size(sacL,2)));                                 %saccade magnitudes
sacL(:,Lpos<eyeThr(1) | Lpos>eyeThr(2))=[];                                %get rid of leftward saccades whose amplitude is >20% smaller than average

upL=reshape(CinterpMat(LL,:),len,size(sacL,2),size(CinterpMat,2));         %all on saccades
tempEye=repmat(e',1,size(CinterpMat,2));
eyL=reshape(tempEye(LL,:),len,size(sacL,2),size(CinterpMat,2));
dnL=reshape(CinterpMat(RR,:),len,size(sacR,2),size(CinterpMat,2));

for W=1;
    e=E(W,:);  
    %[c,pp]=corrcoef([e' cl]);COR=c(2:end,1);pval=pp(2:end,1);%find correlation between all cells and the eye
    COR=ones(size(COR));
    if length(find(COR>0))>0
        CL=cl(:,COR>0);%CL=[zeros(pre,size(CL,2));CL;zeros(fix,size(CL,2))];       

        EE(W)=mean(Lpos(Lpos>eyeThr(1) & Lpos<eyeThr(2)));

        %all off saccades
        %find the saccades following off saccades, find the avg fluor immediately
        %before the cell turns back on, call this zero. I am setting zero
        %separately for each fixation because occasionally there is some drift in
        %the fluorescence so it seems more appropriate to set local rather than
        %global zeros.
        [a,b]=intersect(sacc,sacR(1,:)+pre);
        b=b+1;
        b(b>length(sacc))=length(sacc);
        sacRz=sacc(b);
        p3=round(3/delt);
        p1=round(1/delt);
        p05=round(.5/delt);
        sacRz=repmat(sacRz,pre,1)+repmat((-p3:-p1-1)',1,length(sacRz));LL=reshape(sacRz,numel(sacRz),1);
        dnz=repmat(mean(reshape(CL(LL,:),pre,size(sacRz,2),size(CL,2)),1),len,1);
        %         dnz=dnz*0;
        stadL=squeeze(mean(dnL-dnz,2));
        stad(:,COR>0)=stadL;
        clear tauL RL conInt;
        %%this is the old cirf fit, doesn't have conf intervals
        for i=1:size(stadL,2);
            k=ezfit(T(pre+1:end),stadL(pre+1:end,i),'a*exp(-x/t)');
            tauL(i)=k.m(2);RL(i)=k.r;
        end;

        TAU(COR>0)=tauL;cirfR(COR>0)=RL;
        for qq=1:20;
        if size(dnL,2)<2;
            dosR=NaN(size(dnL,1),size(dnL,3));
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