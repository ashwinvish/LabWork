% plot fits and eye position
function [h] = eyeFits(fits, eyePos, C, time)
% fits are the regressed or fitted eye positions
% eyePos is the positon of the eye
% C is the interpolated eye df/f
% time is the time period
% rsqs is the rsquare error between C and fits

Cells =  size(fits,2);
Cmax = cellfun(@max,C);

%rsq=1-sum((C-fits).^2)/sum((C-mean(C)).^2);
%cor=corr(C,fits,'type','spearman');
delt = mean(diff(time)); % time steps
v=[0,diff(eyePos)/delt]; % eye velocity

vl=v; % left eye velocity
vr=-v; % right eye velocity

vl(vl<0)=0;
vr(vr<0)=0;

el=eyePos;
er=-eyePos;

el(el<0)=0;
er(er<0)=0;



[m,n] = sort(corr);

cols = colorcet('r3');
subplot(1,2,1)
for i = 1:Cells
    col = datasample(cols,1); 
    plot(time, i+ (C{n(i)}/max(C{n(i)})),'Color', cols(i,:), 'LineWidth',1.5);%[col,0.7]);
    hold on;
    text (max(time)+10,i, num2str(n(i)));
    axis off;
end
plot(time, (5+i)+ (eyePos./max(eyePos)),'color','k','LineWidth',1.5);

[~,locU] = findpeaks(diff(eyePos),'MinPeakHeight',1.5);
[~,locD] = findpeaks(diff(-1*eyePos),'MinPeakHeight',1.5);

xU = [locU locU];
xU = repmat(time(locU),2,1);
yU = [ones(size(locU,1),1) 5+i.*ones(size(locU,1),1)];
yU = repmat(yU,size(xU,2),1);

xD = [locD locD];
xD = repmat(time(locD),2,1);
yD = repmat(yU(1,:),size(xD,2),1);

l1 = line(xU,yU','Color',[0,0,1,0.7]);
l2 = line(xD,yD','Color', [1,0,0,0.7]);


subplot(1,2,2)
for i = 1:Cells
    plot(time, i+(fits{n(i)}/max(fits{n(i)})),'color',cols(i,:), 'LineWidth',1.5);% col);
    hold on;
    text (max(time)+10, i, num2str(n(i)));
    axis off;
end
plot(time, (5+i)+ (eyePos./max(eyePos)),'color','k','LineWidth',1.5);

l1 = line(xU,yU','Color',[0,0,1,0.7]);
l2 = line(xD,yD','Color', [1,0,0,0.7]);

end
