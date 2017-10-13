t=0:.01:5;n=20;taus=linspace(1,20,n);for i=1:n;r(:,i)=exp(-t/taus(i));end;
int=diff(r)/.01*.1+r(2:end,:);beta=pinv(r(2:end,:))*int;
FUN=inline('xdata*x','x','xdata');

xdata=r(2:end,:);
for i=1:n
X(:,i) = lsqcurvefit(@(x,xdata) FUN(x,r(2:end,:))...
    ,.1*ones(n,1),r(2:end,:),int(:,i),zeros(n,1),ones(n,1));end
[u,s,v]=svd(r(2:end,:));win=v'*X;x=win(:,1);y=int(:,1);nr=u*s;
for K=1:n
    nx(:,K)=linspace(x(K,1)-,x(K,1)+x(K,1)*.1,20);
    x=win(:,1);
    for i=1:length(nx);
        x(K)=nx(i,K);
        chi(i,K)=sum((int(:,1)-nr*x).^2);
    end
end
rlchi=sum((int(:,1)-nr*win(:,1)).^2);
ind=1:size(chi,1);rl=win(:,1);

hess=2*xdata'*xdata;jac=2*((int-xdata*X)'*(-xdata));

[u,s,v]=svd(hess+jac);
win=v'*beta;
win(5:end,:)=rand(16,n);
new=v*win;