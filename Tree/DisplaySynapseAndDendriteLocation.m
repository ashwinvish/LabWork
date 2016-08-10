function [] = DisplaySynapseAndDendriteLocation(group, cellIDs, allTrees, allPost, groupColor)


GroupTreePathLenghts = [];
GroupPostSynapsePathLength = [];
GroupTreeNormalizedPathLength = [];
GroupNormalizedPostSynapsePathLength = [];

for i = 1:numel(cellIDs)
    if ismember (cellIDs{i}, group) ==1
        TreeNodePoints =  NodePoints(allTrees{i});
        TreepathLength = findPathLength_new([cellIDs{i} , '_WithTags.swc'],allTrees{i},[5,5,45],TreeNodePoints);
        PostSynapsePathLength = findPathLength_new([cellIDs{i} , '_WithTags.swc'],allTrees{i},[5,5,45],allPost{i});
        GroupTreePathLenghts = [GroupTreePathLenghts;TreepathLength];
        GroupTreeNormalizedPathLength = [GroupTreeNormalizedPathLength; NormalizedLength(TreepathLength)];
        GroupPostSynapsePathLength = [GroupPostSynapsePathLength;   PostSynapsePathLength];
        GroupNormalizedPostSynapsePathLength = [GroupNormalizedPostSynapsePathLength; NormalizedLength(PostSynapsePathLength)];
        clear TreepathLength;
        clear PostSynapsePathLength
        clear den;
    end
end

[h,p] = kstest2(GroupTreeNormalizedPathLength,GroupNormalizedPostSynapsePathLength);
sprintf('P value is %f',p)


subplot(1,2,1)
histogram(GroupTreePathLenghts/1000, 'BinWidth', 10, 'Normalization', 'Probability');
hold on;
histogram(GroupPostSynapsePathLength/1000, 'BinWidth', 10, 'Normalization', 'Probability');
axis square
box off;
ylabel('Probability');
xlabel('Pathlength in \mum');
set(gca,'FontName','Arial', 'FontSize', 40, 'LineWidth',2);
set(gca,'color', [groupColor,0.2]);
legend({'Dendrite','PostSynapse'},'FontName','Arial', 'FontSize', 40, 'box', 'off');


subplot(1,2,2);
histogram(GroupTreeNormalizedPathLength, 'BinWidth', 0.1, 'Normalization', 'Probability');
hold on;
histogram(GroupNormalizedPostSynapsePathLength, 'BinWidth', 0.1, 'Normalization', 'Probability');
axis square
box off;
ylabel('Probability');
xlabel('Norm. pathlength');
set(gca,'FontName','Arial', 'FontSize', 40, 'LineWidth',2);
set(gca,'color', [groupColor,0.2]);


end



