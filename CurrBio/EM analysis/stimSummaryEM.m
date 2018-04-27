function FLUOR=stimSummary3(FLUOR);
% cd C:\Users\kpd7\Documents\MATLAB\2011stim
if isfield(FLUOR,'STA')==1;
    FLUOR=rmfield(FLUOR,'STA');FLUOR=rmfield(FLUOR,'stimTA')
end
for QQ=1:length(FLUOR)
    for Q=1:length(FLUOR(QQ).ROI);
        global delt
        if size(FLUOR(QQ).ROI(Q).theta,2)==4;
            stimm=1;
        else
            stimm=0;
        end
        ROI=FLUOR(QQ).ROI(Q);
        a=[ROI.time(end,:) ROI.theta(end,1)];
        a(a<2)=[];strt=max([ROI.time(1,:) ROI.theta(1,1)]);
        tend=min(a);delt=.05;t=strt:delt:tend;
        if stimm==0;
            E=interp1(ROI.theta(:,1),ROI.theta(:,2),t,'linear');
           % E=interp1(ROI.theta(:,1),ROI.theta(:,3),t,'linear'); % changed by AV
            k=sort(E);l=length(k);
            bot=k(floor(l*.05));top=k(floor(l*.95));
            E=E-((top-bot)/2+bot);
            e=(E)*30/(top-bot);
        else
            E=ROI.theta(:,4);dll=floor(1/mean(diff(ROI.theta(:,1))));
            E=interp1(ROI.theta(1:end,1),E,t,'linear');
            e=max(E)-E;
        end
        a=find(isnan(E)==1);E(a)=[];t(a)=[];e(a)=[];
        cl=zeros(length(t),size(ROI.intensity,2));rsqs=zeros(1,size(cl,2));
        for i=1:size(ROI.intensity,2);
            if isnan(ROI.intensity(1,i))==0
                if any(ROI.intensity(:,i))==1
                    a=[ROI.intensity(:,i)];
                    a=(a-mean(a))/mean(a);
                    sp=a;
                    cl(:,i)=interp1(ROI.time(:,i),sp,t,'linear');
                    [cl(:,i),ey,bs,rsqs(i)]=baseline_regress(cl(:,i),e,t,e*0);
                end;end;end;
        FLUOR(QQ).ROI(Q).rsqs=rsqs;
        eyeThr=[3 140];fixLength=7;
        if stimm==0
            STA(Q)=...
                spontSTtrials(Q,t,cl,e,eyeThr,fixLength);
        else
            FLUOR(QQ).stimTA(Q)=stimTrials(t,cl,e);
        end
    end
    if stimm==0;
        FLUOR(QQ).STA=STA;
    end
end
ST=interpStim(e);
FLUOR(2).ST=ST;
%%