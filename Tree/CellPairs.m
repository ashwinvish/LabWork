A = 22;
B = 22;
clear temp1;
clear temp2;
temp1 = allPreSynapse{A}; temp2 = allPostSynapse{A};
treeVisualizer(allTrees{A}, [1],[],[{temp2} {temp1}],false,{[1,0.5,0]}, 1:numel(allTrees{A}), false); 
% temp3 = allPreSynapse{B}; temp4 = allPostSynapse{B};
% treeVisualizer(allTrees{B}, [1],[],[{temp4} {temp3}],false,{[1,0.5,0]}, 1:numel(allTrees{B}), false); 

%%
%for i=1:numel(allTrees)
    
    A = 1;
    B = 6;
    %subplot(3,8,i);
    clear temp1;
    clear temp2;
    temp1 = allPreSynapse{A}; temp2 = allPostSynapse{A};
    [C1] = TreeConvexHull(allTrees{A}, [1],[],[{temp2} {temp1}],false,{[rand,rand,rand]},'green',1:numel(allTrees{A}));
    temp3 = allPreSynapse{B}; temp4 = allPostSynapse{B};
    [C2] = TreeConvexHull(allTrees{B}, [1],[],[{temp4} {temp3}],false,{[rand,rand,rand]},'red',1:numel(allTrees{B}));
%end
