% Plots to generate the major populations

colors = cbrewer('qual','Set1',10);
subplot(1,4,1)
transform_swc_AV(TVNs,colors(1,:),[238],false);
subplot(1,4,2)
transform_swc_AV(MVNs,colors(2,:),[186],false);
subplot(1,4,3)
transform_swc_AV(CELL_Id,colors(3,:),[186],false);
subplot(1,4,4)
transform_swc_AV(vestibularCellIds,colors(4,:),[],false);

