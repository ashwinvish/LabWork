% use termial nodes to calculate axonal pathlength

for i = 1:size(cellIDs,2)
    if eval([cellIDs{i},'_axon'])>0
        AxnNodes = eval([cellIDs{i},'_axon']);
        TerminalNodes = find(Terminals{i}==0);
        CommonNodes = AxnNodes(ismember(AxnNodes,TerminalNodes)==1);
        for jj = 1:numel(CommonNodes)
            treeNodes(jj,1:3) = allTrees{i}{CommonNodes(jj)}{1,3};
        end
        allLengthToPreNode{i} = findPathLengthNew([cellIDs{i} , '_WithTags.swc'],allTrees{i},[5,5,45],treeNodes);
        treeNodes = [];
        %temp(1:size(allLengthToPreNode{i},1)-1,i) = diff(sort(allLengthToPreNode{i}));
        %axLength = [axLength,sum(temp(:,i))];
        
    else
        %axLength = [axLength,0];
        continue;
    end
end