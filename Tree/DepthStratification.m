%
clc
postSynapses = [];
preSynapses = [];
figure();
subplot(1,2,1);
for i = 1:numel(cellIDs)
    if ismember (cellIDs{i}, cellIDsAlx)==1
        DisplayTree(allTrees{i},[1],false, [eval([cellIDs{i},'_axon'])],[rand,rand, rand],allPreSynapse{i}, allPost{i});
        postSynapses =  [postSynapses; allPost{i}];
        preSynapses = [preSynapses; allPreSynapse{i}];
    else
        continue;
    end
end
axis vis3d;
view(-180,0);

h2 = subplot(1,2,2);
%figure();
postSynapses = sortrows(postSynapses,3);
preSynapses = sortrows(preSynapses,3);
[hist1,edge1] = histcounts(postSynapses(:,3)/1000,'BinWidth',1);%,'Orientation','horizontal', 'FaceColor',[0.9,0,0]);
%hold on;
[hist2,edge2] = histcounts(preSynapses(:,3)/1000,'BinWidth',1);%,'Orientation','horizontal', 'FaceColor',[0,0.8,0] );
plot(edge1,hist1,'-');
plot(edge2,hist2,'-');

axis normal;
box off;
set(h2,'Units','Normalized','Position',[0.5,0.41,0.2,0.2],'YDir','reverse');
hold off;

% figure();
% postSynapses = [];
% preSynapses = [];
% figure();
% subplot(1,2,1);
% for i = 1:numel(cellIDs)
%     if ismember (cellIDs{i}, cellIDsDbx)==1
%         DisplayTree(allTrees{i},[1],false, [eval([cellIDs{i},'_axon'])],[rand,rand, rand],allPreSynapse{i}, allPost{i});
%         postSynapses =  [postSynapses; allPost{i}];
%         preSynapses = [preSynapses; allPreSynapse{i}];
%     else
%         continue;
%     end
% end
% axis vis3d;
% view(-180,0);
% 
% h2 = subplot(1,2,2);
% %figure();
% postSynapses = sortrows(postSynapses,3);
% preSynapses = sortrows(preSynapses,3);
% histogram(h2,postSynapses(:,3)/1000,'BinWidth',5,'DisplayStyle','stairs','Orientation','horizontal', 'FaceColor',[0.9,0,0]);
% hold on;
% histogram(h2,preSynapses(:,3)/1000,'BinWidth',5,'DisplayStyle','stairs','Orientation','horizontal', 'FaceColor',[0,0.8,0] );
% axis normal;
% box off;
% set(h2,'Units','Normalized','Position',[0.5,0.41,0.2,0.2],'YDir','reverse');
% hold off;
% 
% figure();
% postSynapses = [];
% preSynapses = [];
% figure();
% subplot(1,2,1);
% for i = 1:numel(cellIDs)
%     if ismember (cellIDs{i}, cellIDsL)==1
%         DisplayTree(allTrees{i},[1],false, [eval([cellIDs{i},'_axon'])],[rand,rand, rand],allPreSynapse{i}, allPost{i});
%         postSynapses =  [postSynapses; allPost{i}];
%         preSynapses = [preSynapses; allPreSynapse{i}];
%     else
%         continue;
%     end
% end
% axis vis3d;
% view(-180,0);
% 
% h2 = subplot(1,2,2);
% %figure();
% postSynapses = sortrows(postSynapses,3);
% preSynapses = sortrows(preSynapses,3);
% histogram(h2,postSynapses(:,3)/1000,'BinWidth',5,'DisplayStyle','stairs','Orientation','horizontal', 'FaceColor',[0.9,0,0]);
% hold on;
% histogram(h2,preSynapses(:,3)/1000,'BinWidth',5,'DisplayStyle','stairs','Orientation','horizontal', 'FaceColor',[0,0.8,0] );
% axis normal;
% box off;
% set(h2,'Units','Normalized','Position',[0.5,0.41,0.2,0.2],'YDir','reverse', 'YLim',[0 ,60]);
% hold off;
