function [ postSynapses, preSynapses,hist1,hist2 ] = DepthStratification( cellID, allTrees, cellIDs )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

load CellAxons
subplot(1,2,1);

DisplayTree(allTrees{cellID},[1],false, [eval([cellIDs{cellID},'_axon'])],[rand,rand, rand],allPreSynapse{cellID}, allPost{cellID});
postSynapses =  [postSynapses; allPost{i}];
preSynapses = [preSynapses; allPreSynapse{i}];

axis vis3d;
view(-180,0);

h2 = subplot(1,2,2);

postSynapses = sortrows(postSynapses,3);
preSynapses = sortrows(preSynapses,3);
hist1 = histogram(h2,postSynapses(:,3)/1000,'BinWidth',1,'Orientation','horizontal', 'FaceColor',[0.9,0,0]);
hold on;
hist2 = histogram(h2,preSynapses(:,3)/1000,'BinWidth',1,'Orientation','horizontal', 'FaceColor',[0,0.8,0] );
axis normal;
box off;
set(h2,'Units','Normalized','Position',[0.5,0.41,0.2,0.2],'YDir','reverse');
hold off;


end

