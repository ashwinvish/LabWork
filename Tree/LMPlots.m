% plot all LM traces in same orintation with axon in blue
    ha = tight_subplot(2,3,[.05 .05],[.05 .1],[.01 .01]);

    treeVisAV('fish3003_ch1-exported-000.swc',1);
    axes(ha(1));
    treeVisAV('fish2017_ch2.swc',1);
    axes(ha(2));
    treeVisAV('fish1059_ch2-axons.swc',1);
    axes(ha(3));
    treeVisAV('fish3075_118-axons.swc',1);
    axes(ha(4));
    treeVisAV('fish1006-axon.swc',1);
    axes(ha(5));
    treeVisAV('fish1013_ch2.swc',1);
    axes(ha(6));
    
% set(ha(1:6),'BoxStyle','full');
% %title(sprintf('CellID %s', cellIDs{kk}));
% set(gcf,'color','w');
%     