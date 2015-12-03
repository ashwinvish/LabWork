% plot all LM traces in same orintation with axon in blue, dendrtic tree is
% not complete, so in some cases the root node is misleading.

clear LMPathLen;
ha = tight_subplot(2,3,[.05 .05],[.05 .1],[.01 .01]);

tree{1} = treeVisAV('fish3003_ch1-exported-000.swc',1);
%LMPathLen{1} = LMPathLengths('fish3003_ch1-exported-000.swc',1,2)
axes(ha(1));

tree{2} = treeVisAV('fish2017_ch2.swc',1);
%LMPathLen{2} = LMPathLengths('fish2017_ch2.swc',1,2)
axes(ha(2));

tree{3} = treeVisAV('fish1059_ch2-axons.swc',1);
%LMPathLen{3} = LMPathLengths('fish1059_ch2-axons.swc',1,2)
axes(ha(3));

tree{4} = treeVisAV('fish3075_118-axons.swc',1);
%LMPathLen{4} = LMPathLengths('fish3075_118-axons.swc',1,2)
axes(ha(4));

tree{5} = treeVisAV('fish1006-axon.swc',1);
%LMPathLen{5} = LMPathLengths('fish1006-axon.swc',1,2)
axes(ha(5));

tree{6} = treeVisAV('fish1013_ch2.swc',1);
%LMPathLen{6} = LMPathLengths('fish1013_ch2.swc',1,2)
axes(ha(6));

  