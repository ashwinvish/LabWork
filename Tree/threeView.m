function [area] = threeView(volume, Soma, cmap, varargin)
%threeView will genetrate the three view of a volume
%   volume is the 3D array
%   Soma is the cartesian coordinates of the somas
%   cmap is the colormap, e.g. jet
%   varargin are the variable arguments, e.g specific planes to be plotted


trans = 0.5; % if transparency is needed
if nargin == 3
    %plot XY
    tempXY = max(volume,[],3);
    [m,n] = size(tempXY);
    h1 = imagesc(0,70,imrotate(tempXY,-90));
    set(h1,'AlphaData', trans);
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
    set(h2,'AlphaData', trans);
    colormap(cmap);axis image;
    
    %plot YZ
    tempYZ = max(volume,[],1);
    for i = 1:size(tempYZ,3)
        C(i,:) = tempYZ(:,:,i);
    end
    h3 = imagesc(130,70,imrotate(flipdim(C,1),-90));
    set(h3,'AlphaData', trans);
    set(gca,'YDir','normal');
    colormap(cmap);
    axis image;
    axis vis3d;
    
    %plot XZ soma location
    scatter(140-Soma(:,1), 60 -Soma(:,3),150,'Marker','o', 'MarkerFaceColor','w', 'MarkerEdgeColor','k','LineWidth', 2);
    %plot XY soma location
    scatter(140-Soma(:,1), 10 + Soma(:,2),150,'Marker','o', 'MarkerFaceColor','w', 'MarkerEdgeColor','k','LineWidth', 2);
    %plot YZ soma location
    scatter(130+Soma(:,3), 10 + Soma(:,2),150,'Marker','o', 'MarkerFaceColor','w', 'MarkerEdgeColor','K','LineWidth', 2);
    
    axis off ;
    
else if nargin> 3
        vol = 0;
        tempXY = varargin{1};
        imagesc(0,70,imrotate(tempXY,-90));
        h1 = gcf;
        set(h1.CurrentAxes, 'LineWidth', 2);
        colormap(cmap);
        hold on;
        axis image;
        
        %plot XZ
        tempXZ = varargin{2};
        for i = 1:size(tempXZ,3)
            B(:,i) = tempXZ(:,:,i);
        end
        imagesc(0,0,imrotate(flipdim(B,2),-90));
        h2 = gcf;
        set(h2.CurrentAxes, 'LineWidth', 2);
        colormap(cmap);
        axis image;
        
        %plot YZ
        tempYZ = varargin{3};
        for i = 1:size(tempYZ,3)
            C(i,:) = tempYZ(:,:,i);
        end
        imagesc(130,70,imrotate(flipdim(C,1),-90));
        h3 = gcf;
        set(h3.CurrentAxes, 'LineWidth', 2);
        set(gca,'YDir','normal');
        colormap(cmap);
        axis image;
        axis vis3d;
        
         axis off;
         axis vis3d;
        
        % calculate the overlap area in downsampled coordinates
        area = sum(tempXY(:)); % should be same for all planes.
        
        %
        Soma2 =  varargin{4};
        
        %plot XZ soma location
        scatter(140-Soma(:,1), 60 -Soma(:,3),50,'Marker','o',   'MarkerFaceColor','none','MarkerEdgeColor','red','LineWidth', 2);
        %plot XY soma location
        scatter(140-Soma(:,1), 10 + Soma(:,2),50, 'Marker','o',  'MarkerFaceColor','none','MarkerEdgeColor','red','LineWidth', 2);
        %plot YZ soma location
        scatter(130+Soma(:,3), 10 + Soma(:,2),50, 'Marker','o',  'MarkerFaceColor','none','MarkerEdgeColor','red','LineWidth', 2);
        
        if ~isempty(Soma2)
        
        %plot XZ soma location
        scatter(140-Soma2(:,1), 60 -Soma2(:,3),50, 'Marker','o', 'MarkerFaceColor','none','MarkerEdgeColor','green','LineWidth', 2);
        %plot XY soma location
        scatter(140-Soma2(:,1), 10 + Soma2(:,2),50,'Marker','o',  'MarkerFaceColor','none','MarkerEdgeColor','green','LineWidth', 2);
        %plot YZ soma location
        scatter(130+Soma2(:,3), 10 + Soma2(:,2),50,'Marker','o',  'MarkerFaceColor','none','MarkerEdgeColor','green','LineWidth', 2);
        
        end
        
        
    end
    
end


