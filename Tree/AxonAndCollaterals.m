% plot all Colaterals of Int1_5
%subplot(1,2,1);
DisplayTree(allTrees{5},[1],true,eval([cellIDs{5},'_axon']),[1 0.5 0.3]);
title('Complete axon of Int1_5');

%subplot(1,2,2);
% recurrent collaterals
DisplayTree(allTrees{5},[1],true,[117,124,128,131,132,134,136,138,140,141,142,143,144,152,154,155,156,157,158],[1 0.5 0.3]);
DisplayTree(allTrees{5},[1],false,[103,117,123,125,126,127,137,139,150,151,153],[1 0.5 0.3]);
DisplayTree(allTrees{5},[1],false,[92,99,100,103,115],[1 0.5 0.3]);
DisplayTree(allTrees{5},[1],false,[77,92,106,110,114,116,133,135,145,146,147,148,149],[1 0.5 0.3]);
DisplayTree(allTrees{5},[1],false,[159,161,163],[1 0.5 0.3]);
DisplayTree(allTrees{5},[1],false,[163,165,168],[1 0.5 0.3]);
DisplayTree(allTrees{5},[1],false,[163,165,166,167,169,170],[1 0.5 0.3]);
DisplayTree(allTrees{5},[1],false,[77,80,81,82,84,85,86],[1 0.5 0.3]);
% long range axons
% DisplayTree(allTrees{5},[1],false,[56,107,159],[1 0.5 0.3]);
% DisplayTree(allTrees{5},[1],false,[56,77],[1 0.5 0.3]);


% plot all colaterals for Int1_4














