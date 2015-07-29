function [tempXY,tempXZ,tempYZ] = threeView( volume, cmap , varargin)
%threeView will genetrate the three view of a volume
%   volume is the 3D array
%   cmap is the colormap, e.g. jet
%   varargin are the variable arguments, e.g specific planes to be plotted

if nargin == 2
    %plot XY
    tempXY = max(volume,[],3);
    h1 = imagesc(0,90,imrotate(tempXY,-90));
    %set(h1,'AlphaData', trans);
    set(gca,'YDir','normal');
    colormap(cmap);
    hold on;
    axis image;
    
    %plot XZ
    tempXZ = max(volume,[],2);
    for i = 1:size(tempXZ,3)
        B(:,i) = tempXZ(:,:,i);
    end
    h2 = imagesc(0,0,imrotate(flipdim(B,2),-90));
    %set(h2,'AlphaData',trans);
    colormap(cmap);
    axis image;
    
    %plot YZ
    tempYZ = max(volume,[],1);
    for i = 1:size(tempYZ,3)
        C(i,:) = tempYZ(:,:,i);
    end
    h3 = imagesc(160,90,imrotate(flipdim(C,1),-90));
    %set(h3,'AlphaData',trans);
    set(gca,'YDir','normal');
    colormap(cmap);
    axis image
    
    axis off;
    axis vis3d;
    
else if nargin> 2
        vol = 0;
        tempXY = varargin{1};
        h1 = imagesc(0,90,imrotate(tempXY,-90));
        %set(h1,'AlphaData', trans);
        set(gca,'YDir','normal');
        colormap(cmap);
        hold on;
        axis image;
        
        %plot XZ
        tempXZ = varargin{2};
        for i = 1:size(tempXZ,3)
            B(:,i) = tempXZ(:,:,i);
        end
        h2 = imagesc(0,0,imrotate(flipdim(B,2),-90));
        %set(h2,'AlphaData',trans);
        colormap(cmap);
        axis image;
        
        %plot YZ
        tempYZ = varargin{3};
        for i = 1:size(tempYZ,3)
            C(i,:) = tempYZ(:,:,i);
        end
        h3 = imagesc(160,90,imrotate(flipdim(C,1),-90));
        %set(h3,'AlphaData',trans);
        set(gca,'YDir','normal');
        colormap(cmap);
        axis image
        
        axis off;
        axis vis3d;

    end
end


