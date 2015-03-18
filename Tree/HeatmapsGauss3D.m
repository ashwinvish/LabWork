clear vol;
res =1000; % resolution of image
ksize = 18; % size of Kernel
Gwin = gausswin(ksize+1);
K = Gwin*Gwin';
event = 1;

for i = 1:(ksize+1)
    K3D(:,:,i) = K(:,:)*Gwin(i);
end

for kk = 1:length(allPost)
    volPost = zeros((15e4-0)/res,(2.5e5-0)/res,(8e4-0)/res); 
    % add 10000/res pixels to the zaxis to avoid edge effects
    for n = 1:length(allPost{kk});
        volPost(round(allPost{kk}(n,1)/res) + (-ksize/2:ksize/2), round(allPost{kk}(n,2)/res) + (-ksize/2:ksize/2), round(allPost{kk}(n,3)/res + 10000/res)+ (-ksize/2:ksize/2)) = ...
            volPost(round(allPost{kk}(n,1)/res) + (-ksize/2:ksize/2), round(allPost{kk}(n,2)/res) + (-ksize/2:ksize/2),round(allPost{kk}(n,3)/res + 10000/res)+ (-ksize/2:ksize/2))  + K3D;
    end
   
    subplot(3,8,kk);
%     vol3d('CData',volPost,'xdata',[0 250], 'ydata', [0 150] , 'zdata', [0 80]);
%     hold on;
%     plot3(160-CellSoma(kk,2)/res,90 + CellSoma(kk,1)/res,(CellSoma(kk,3))/res,'Marker','o', 'MarkerFaceColor','k');
%     colormap(jet);
%     alphamap([0, linspace(0.1,0,255)]);
%     
%     daspect([1,1,1]);
%     view([-50,50]);
%     box on;
    
    
    
    %plot XZ
    tempXZ = max(volPost,[],2);
    for i = 1:size(tempXZ,3)
        B(:,i) = tempXZ(:,:,i);
    end
    imagesc(0,0,imrotate(B,-90)); colormap(jet);
    title(cellIDs{kk});
    hold on;
    plot(160-CellSoma(kk,2)/res,80-(CellSoma(kk,3))/res,'Marker','o', 'MarkerFaceColor','w');
    clear B;
    
    %plot XY
    tempXY = max(volPost,[],3);
    h2 = imagesc(0,90,imrotate(tempXY,-90));set(gca,'YDir','normal'); colormap(jet);
    title(cellIDs{kk});
    plot(160-CellSoma(kk,2)/res, 90 + CellSoma(kk,1)/res,'Marker','o', 'MarkerFaceColor','w');

    %plot yz
    tempYZ = max(volPost,[],1);
    for i = 1:size(tempYZ,3)
        C(i,:) = tempYZ(:,:,i);
    end
    imagesc(160,90,imrotate(C,-90));set(gca,'YDir','normal');
    colormap(jet);
    title(cellIDs{kk});
    plot(160+(CellSoma(kk,3))/res, 90+ CellSoma(kk,1)/res,'Marker','o', 'MarkerFaceColor','w');
    
    axis off;
    axis image;
    clear C;
    
    [m,I] = max(volPost(:));
    maxPostDesnity(kk) = m;
    [x,y,z] = ind2sub(size(volPost),I);
    XPost(kk) = x; YPost(kk) = y; ZPost(kk) = z-10000/res;
end
%%
%Presynapse Heatmap
figure(14);
clear volPre;
clear('tempXY','tempXZ','tempYZ');
clear ('B','C');
jj=1;
for kk = 1:length(allPreSynapse)
    if cellfun('isempty',allPreSynapse(1,kk)) == 1
        continue;
    else
        volPre = zeros((15e4-0)/res,(2.6e5-0)/res,(8e4-0)/res);
        for n = 1:size(allPreSynapse{kk},1);
            volPre(round(allPreSynapse{kk}(n,1)/res) + (-ksize/2:ksize/2), round(allPreSynapse{kk}(n,2)/res) + (-ksize/2:ksize/2), round(allPreSynapse{kk}(n,3)/res + 10000/res)+ (-ksize/2:ksize/2)) = ...
                volPre(round(allPreSynapse{kk}(n,1)/res) + (-ksize/2:ksize/2), round(allPreSynapse{kk}(n,2)/res) + (-ksize/2:ksize/2),round(allPreSynapse{kk}(n,3)/res + 10000/res)+ (-ksize/2:ksize/2))  + K3D;
        end
        subplot(3,8,jj);
        %plot XY
        tempXY = max(volPre,[],3);
        %figure('Units','pixels','Position', [0,0,m,n]);
        imagesc(0,90,imrotate(tempXY,-90));set(gca,'YDir','normal'); colormap(jet);
        title(cellIDs{kk});
        axis off
        hold on;
        plot(160-CellSoma(kk,2)/res, 90 + CellSoma(kk,1)/res,'Marker','o', 'MarkerFaceColor','w');
        axis image;
        %plot XZ
        tempXZ = max(volPre,[],2);
        for i = 1:size(tempXZ,3)
            B(:,i) = tempXZ(:,:,i);
        end
        imagesc(0,0,imrotate(B,-90));set(gca,'YDir','normal'); colormap(jet);
        title(cellIDs{kk});
        axis off
        hold on;
        plot(160-CellSoma(kk,2)/res,80-(CellSoma(kk,3))/res,'Marker','o', 'MarkerFaceColor','w');
        axis image;
        clear B;
        %plot yz
        tempYZ = max(volPre,[],1);
        for i = 1:size(tempYZ,3)
            C(i,:) = tempYZ(:,:,i);
        end
        imagesc(160,90,imrotate(C,-90));set(gca,'YDir','normal'); colormap(jet);
        title(cellIDs{kk});
        axis off
        hold on;
        plot(160+(CellSoma(kk,3))/res, 90+ CellSoma(kk,1)/res,'Marker','o', 'MarkerFaceColor','w');
        axis image;
        clear C;
        jj = jj+1;
    end
    
    axis off;
    axis image;
    clear C;
    
    
    [m,I] = max(volPre(:));
    maxPreDesnity(kk) = m;
    [x,y,z] = ind2sub(size(volPre),I);
    XPre(kk) = x; YPre(kk) = y; ZPre(kk) = z-10000/res;
      
end
