function ROIcellLabel(ROI,cls,clLabel,col)
%ROI is a  structure of fluorescence time series created using
%FluorescenceAnalysis.m, cls are the cells within this structure that you
%choose to label, set clLabel to 1 if you want to label each ROI with the
%cell's number
k=double(max(ROI.mask(:,:,cls),[],3));
[FX,FY] = gradient(k);
mask=sqrt(FX.^2+FY.^2);mask(mask~=0)=1;
im=repmat(ROI.mnImage,[1 1 3]);
im=im/max(im(:));
mx=find(im>max(im(:))*.7);
mn=find(im<max(im(:))*.0);
im(mn)=0;im(mx)=.8;
image(im+.1);colormap(gray);
col(col>.9)=.9;
for i=1:3;
    a=im(:,:,i);
    a(mask==1)=col(i);
    im(:,:,i)=a;
end
image(im+.1);
set(gca,'position',[0 0 1 1]);set(gcf,'toolbar','none');

s=size(im);
cellnum=1;
if clLabel==1
    for i=1:length(ROI.centroid);
        if ismember(i,cls)
            xpos=ROI.centroid(1,i);ypos=ROI.centroid(2,i);
            aa=annotation('textbox','units','normalized','position',[xpos(1)/s(2) .95-ypos(1)/s(1)+.05 .05 .05]);
            set(aa,'linestyle','none','color','w','fontweight','bold','string',num2str(cellnum));
            cellnum=cellnum+1;
        end
    end
end
% set(gcf,'position',[190 100 500 500])
set(gca,'xtick',[],'ytick',[])