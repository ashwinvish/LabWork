%% compare 2 skeletons with each other 

% get file names
%d2 = [pwd '\skeletons\ashleigh\'];
%d1 = [pwd '\skeletons\ashwin\'];

d2 = [pwd '/skeletons/ashleigh/'];
d1 = [pwd '/skeletons/ashwin/'];

[f1,f2,id]=find_samecell_swcfiles(d1,d2);
% boolean variable for files containing axon initiation points
AwinInitBool = false(length(f1),1);
AwinInitBool([[4:7],21],1)=true;

N = length(f1);
for j=5:5
    % get skeletons from swc files
    q = getPoints(f1{j});
    close(gcf)
    q2 = getPoints(f2{j});
    close(gcf)
    
    % calculate the nearest neighbor distance to every point in Ashleigh's
    % file and vice-versa
    [~,D2] = knnsearch(q,q2);
    [~,D1] = knnsearch(q2,q);
        
    
    figure
    highlight_both(q,q2,D1,D2,200,'au')
     legend('ashwin','ashleigh');
    title(num2str(id(j)))
    axis vis3d
 
    if AwinInitBool(j)
        showAxonInit(f1{j},j);
   end
  %    set(gcf,'PaperPosition',[0 0 8.5 11]);
    %    print(gcf,'-dpdf',num2str(id(j)))
end