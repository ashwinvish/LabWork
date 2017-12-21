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
    title('ipsi eye pos');
    subplot(5,1,3)
    plot(SPT(i).eyevars(:,3));
    title('contra eye pos');
    subplot(5,1,4)
    plot(SPT(i).eyevars(:,4));
    title('ipsi eye vel');
    subplot(5,1,5)
    plot(SPT(i).eyevars(:,5));
    title('contra eye vel');
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




