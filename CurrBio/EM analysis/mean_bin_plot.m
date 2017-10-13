function [X,Y,stdEr,C,P]=mean_bin_plot(xx,yy,col,plt,A,color,g)
%function for binning pairwise distance vs. pairwise something data. xx is
%a matrix of distances, each column corresponds to a direction. yy is a
%matrix of corresponding pairwise values (i.e. pairwise difference in
%persistence), plt=1 if you want to output a plot, plt=any number
%otherwise. Col is number of data ponits to reshape the matrix.
if nargin==5;color='bbbbbbb';end
for i=1:size(xx,2);
    x=xx(:,i);y=yy(:,i)/A;
%     rem=.05;
%     [a,b]=sort(abs(y));ind=b([(1:floor(rem*length(b)))...
%         (floor((1-rem)*length(b)):length(b))]);
    a=find(isnan(x)==1);x(a)=[];y(a)=[];
    b=find(x>0);
%     b=1:length(x);
    [c,p]=corr(x(b),y(b),'type','spearman');
    if nargin==7;
        [c,p]=corr(x(b),y(b),'type','pearson');
    end
    [a,b]=sort(x);x=x(b);y=y(b);
    row=floor(length(x)/col);len=row*col;
    x=reshape(x(1:len),row,col);y=reshape(y(1:len),row,col);
    fivePercent=round(.0*size(x,1));
    [a,b]=sort(y);
    for ll=1:size(x,2);
        x(:,ll)=x(b(:,ll),ll);
    end
    y=sort(y);
    x=x(max([fivePercent 1]):end-fivePercent,:);
    y=y(max([fivePercent 1]):end-fivePercent,:);
    X(:,i)=mean(x);Y(:,i)=mean(y);stdEr(:,i)=std(y)/sqrt(row);
    % plot(x,y,'.','color',.6*[1 1 1]);hold on;
    if plt==1
        if size(xx,2)>1;subplot(1,size(xx,2),i);end
        plot(X(:,i),Y(:,i),'k.-','markersize',25,'color','r', 'LineWidth', 2);
        hold on;
        plot(x,y,'k.','markersize',25,'color',color(i), 'LineWidth', 2);
        %errorbar(X(:,i),Y(:,i),stdEr(:,i),'k.-','markersize',25,'color',color(i), 'LineWidth', 2);
        set(gca,'FontName', 'Arial','FontSize',40, 'LineWidth', 2)
        box off;
    end
    C(:,i)=diag(c);P(:,i)=diag(p);
    hold on;
end