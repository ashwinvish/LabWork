function r=svd_rates(f,t,tau,num,cut);

delt=mean(diff(t));
[u,s,v]=svd(f);
fs=f*v(:,1:num(1))*v(:,1:num)';
rs=tau*diff(fs)/delt+fs(1:end-1,:);
a=find(t(2:end)>=cut);
rs=rs(a,:);
[u,s,v]=svd(rs);
r=rs*v(:,1:num(2))*v(:,1:num(2))';