function [cl,eyfit,baseline,rsq,cor,beta]=baseline_regress(cl,e,t,st);
delt=mean(diff(t));
e=medfilt1(e,51);
e=e';st=st';
v=[0;diff(e)/delt];
vl=v;vr=-v;
vl(vl<0)=0;vr(vr<0)=0;
el=e;er=-e;
el(el<0)=0;er(er<0)=0;
x=linspace(-1,1,length(e))';
K=ones(size(e));
sl=st;sr=-st;
sl(sl<0)=0;sr(sr<0)=0;
xmat=[e el er vl vr sl sr];
xmat=xmat-repmat(xmat(1,:),size(xmat,1),1);
for i=1:size(xmat,2);
    a=conv(xmat(:,i),exp(-t/1.9));
    xmat(:,i)=a(1:length(t));
end
stp=size(xmat,2)+1;
xmat=[xmat K x x.^2 x.^3];
for i=1:size(cl,2);
    beta=pinv(xmat)*cl(:,i);
    fit(:,i)=xmat*beta;
end
baseline=xmat(:,stp:end)*beta(stp:end);
eyfit=xmat(:,1:stp-1)*beta(1:stp-1);
cl=cl-baseline;
beta=beta(stp:end);
rsq=1-sum((cl-eyfit).^2)/sum((cl-mean(cl)).^2);
k=sort(cl);mn=k(floor(length(k)*.2));a = find(cl<mn);
cor=corr(cl,eyfit,'type','spearman');
%k=sort(eyfit);mn=k(floor(length(k)*.2));a=find(eyfit<mn);
%a = find(cl<mn);
mn=mean(cl(a));
cl=cl-mn;