%
load('./101112 _files1_4.mat')

% view the eyevars data

for i = 1:size(SPT,2)
    figure(i)
    subplot(5,1,1)
    plot(SPT(i).eyevars(:,1));
    title('raw eye pos');
    subplot(5,1,2)
    plot(SPT(i).eyevars(:,2));
    title('ipsi side');
    subplot(5,1,3)
    plot(SPT(i).eyevars(:,3));
    title('contra side');
    subplot(5,1,4)
    plot(SPT(i).eyevars(:,4));
    title('ipsi eye vel?');
    subplot(5,1,5)
    plot(SPT(i).eyevars(:,5));
    title('contra eye vel?');
    suptitle(SPT(i).filename);
end

%view tmap
for i = 1:size(SPT,2)
    figure(i)
    subplot(2,2,1)
    imagesc(SPT(i).tMap(:,:,1));
    colormap gray;
    axis square;
    subplot(2,2,2)
    imagesc(SPT(i).tMap(:,:,2));
    colormap gray;
    axis square;
    subplot(2,2,3)
    imagesc(SPT(i).tMap(:,:,3));
    colormap gray;
    axis square;
    subplot(2,2,4)
    imagesc(SPT(i).tMap(:,:,4));
    colormap gray;
    axis square;
    suptitle(SPT(i).filename);
end

% view mnImage (Median?)
figure
for i = 1:size(SPT,2)  
    subplot(2,2,i)
    imagesc(SPT(i).mnImage);
    hold on
    scatter(SPT(i).centroid(1,:),SPT(i).centroid(2,:) ,'r.');
    colormap gray;
    axis square;
    axis off
    title(SPT(i).filename);
end

% plot correct fluorescence intensities for all cells
for i = 1:size(SPT,2)
    figure(i);
    n = round(SPT(i).numROIs/5);
    subplot(n,5,1)
    plot(SPT(i).eyevars(:,1));
    set(gca,'XLim',[0,300]);
    for j = 1:SPT(i).numROIs
        subplot(n,6,j+1)
        plot(SPT(i).intCorrect(:,j));
        set(gca,'XLim',[0,300]);
    end
end
        
        






