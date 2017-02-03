AlxPostSites = [];
AlxPreSites = [];
TransPostSites = [];
TransPreSites = [];
DbxPostSites = [];
BarhlPostSites = [];

AlxSoma = [];
TransSoma = [];
BarhlSoma = [];
DbxSoma = [];


for kk = 1:numel(cellIDs)
    if ismember(cellIDs{kk},cellIDsAlx)==1
        AlxPostSites  = [AlxPostSites; allPost{kk}];
        AlxPreSites = [AlxPreSites; allPreSynapse{kk}];
        AlxSoma = [AlxSoma; CellSoma(kk,:)];
    elseif ismember(cellIDs{kk},cellIDsTrans)==1
        TransPostSites  = [TransPostSites; allPost{kk}];
        TransPreSites = [TransPreSites; allPreSynapse{kk}];
        TransSoma = [TransSoma; CellSoma(kk,:)];
    elseif ismember(cellIDs{kk},cellIDsDbx)==1
        DbxPostSites = [DbxPostSites; allPost{kk}];
        DbxSoma = [DbxSoma; CellSoma(kk,:)];
    else
        BarhlPostSites = [BarhlPostSites; allPost{kk}];
        BarhlSoma = [BarhlSoma; CellSoma(kk,:)];
    end
end

res = 1000;                                                                                          % downsampling factor
ksize = 20000;

figure('units','normalized','outerposition',[0 0 1 1]); 
AlxPostVol = HeatMapFish(ksize,res, AlxPostSites,AlxSoma,cellIDsAlx,true);
title('Ipsilateral postsynaptic sites', 'FontName', 'Arial', 'FontSize', 40);
colormap(jet);
c1 = colorbar;
c1.Label.String = 'sites/\mum';
c1.Location = 'westoutside';
c1.FontSize = 40;
set(gcf, 'color', 'none');
%export_fig('/Users/admin/Dropbox/ZFishPresentation/IpsiPostHeatMap.png');

figure('units','normalized','outerposition',[0 0 1 1]);
AlxPreVol = HeatMapFish(ksize,res, AlxPreSites,AlxSoma,cellIDsAlx,true);
title('Ipsilateral presynaptic sites', 'FontName', 'Arial', 'FontSize', 40);
colormap(jet);
c2 = colorbar;
c2.Label.String = 'sites/\mum';
c2.FontSize = 40;
set(gcf,'color', 'none');
%export_fig('/Users/admin/Dropbox/ZFishPresentation/IpsiPreHeatMap.png');


figure('units','normalized','outerposition',[0 0 1 1]);
HeatMapFish(ksize,res, TransPostSites,TransSoma,cellIDsTrans,true);
title('Ipsi-Contra postsynaptic sites','FontName', 'Arial', 'FontSize', 40);
colormap(jet);
c1 = colorbar;
c1.Label.String = 'sites/\mum';
c1.Location = 'westoutside';
c1.FontSize = 40;
set(gcf, 'color', 'none');
%export_fig('/Users/admin/Dropbox/ZFishPresentation/IpsiContraPostHeatMap.png');

figure('units','normalized','outerposition',[0 0 1 1]);
HeatMapFish(ksize,res, TransPreSites,TransSoma,cellIDsTrans,true);
title('Ipsi-Contra presynaptic sites', 'FontName', 'Arial', 'FontSize', 40);
colormap(jet);
c2 = colorbar;
c2.Label.String = 'sites/\mum';
c2.FontSize = 40;
set(gcf, 'color', 'none');
%export_fig('/Users/admin/Dropbox/ZFishPresentation/IpsiContraPreHeatMap.png');


figure('units','normalized','outerposition',[0 0 1 1]);
HeatMapFish(ksize,res, DbxPostSites,DbxSoma,cellIDsDbx,true);
title('Contra postsynaptic sites','FontName', 'Arial', 'FontSize', 40);
colormap(jet);
c1 = colorbar;
c1.Label.String = 'sites/\mum';
c1.Location = 'westoutside';
c1.FontSize = 40;
set(gcf, 'color', 'none');
%export_fig('/Users/admin/Dropbox/ZFishPresentation/ContraPostHeatMap.png');


figure('units','normalized','outerposition',[0 0 1 1]);
HeatMapFish(ksize,res, BarhlPostSites,BarhlSoma,cellIDsL,true);
title('Unknown postsynaptic sites','FontName', 'Arial', 'FontSize', 40);
colormap(jet);
c1 = colorbar;
c1.Label.String = 'sites/\mum';
c1.Location = 'westoutside';
c1.FontSize = 40;
set(gcf, 'color', 'none');
%export_fig('/Users/admin/Dropbox/ZFishPresentation/UnknownPostHeatMap.png');

%% All Sites


HeatMapFish(ksize,res, [AlxPostSites; TransPostSites;DbxPostSites;BarhlPostSites],[AlxSoma;TransSoma; DbxSoma; BarhlSoma],cellIDs,true);

%%
figure();
HeatMapFish(ksize,res, [AlxPreSites; TransPreSites],[AlxSoma;TransSoma; DbxSoma; BarhlSoma],cellIDs,true);

