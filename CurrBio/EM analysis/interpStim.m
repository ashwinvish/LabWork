function ST=interpStim(e);
delt=.05;
t=-2:delt:7;
T=-2:1:7;
ST=[];
for i=1:100;
    strt=round(rand*20);strt(strt<1)=1;
    a=e(strt:20:end);
    b=find(a>.5);
    clear st
    b(b+7>length(a))=[];b(b-2<1)=[];
    for j=1:length(b);
        st=interp1(T,a(b-2:b+7),t,'linear');
        ST=[ST st'];
    end
end
        
ST=mean(ST,2);
    