% To calculate the overlapping area between cells

ksize = 12000;
res = 1000;

% compute overlapping area for all cell pairs

for i = 1:size(cellIDs,2)
    if ~cellfun('isempty',allPreSynapse(1,i)) == 1					% do only for cells that have presynaptic sites
        figure;
        for ii = 1:size(cellIDs,2)
            h = subplot(3,8,ii);
            [area] = dotVol(volPre{i},volPost{ii},CellSoma(i,:),CellSoma(ii,:),res);
            if area == 0
                delete (h);
            end
            IntArea(i,ii) = area;
            str = sprintf('Presynaptic cell (red): %s \nPostSynaptic cell (green): %s',cellIDs{i},cellIDs{ii});
            title(str,'FontSize',5);
        end
    else
        continue;
    end
end

AlxAlxArea = [];
AlxTransArea = [];
AlxDbxArea = [];
AlxBarhlArea = [];


for i = 1:numel(cellIDs)
    for j = 1:numel(cellIDs)
        if ismember(cellIDs{i}, cellIDsAlx) == 1
            if ismember(cellIDs{j}, cellIDsAlx) == 1
                AlxAlxArea= [AlxAlxArea, IntArea(i,j)];
            end
            if ismember(cellIDs{j}, cellIDsDbx) == 1
                AlxDbxArea = [AlxDbxArea, IntArea(i,j)];
            end
            if ismember(cellIDs{j}, cellIDsTrans) ==1
                AlxTransArea = [AlxTransArea, IntArea(i,j)];
            end
            if ismember(cellIDs{j}, cellIDsL) ==1
                AlxBarhlArea = [AlxBarhlArea, IntArea(i,j)];
            end
        end
    end
end

TransAlxArea = [];
TransTransArea = [];
TransDbxArea = [];
TransBarhlArea = [];


for i = 1:numel(cellIDs)
    for j = 1:numel(cellIDs)
        if ismember(cellIDs{i}, cellIDsTrans) == 1
            if ismember(cellIDs{j}, cellIDsAlx) == 1
                TransAlxArea= [TransAlxArea, IntArea(i,j)];
            end
            if ismember(cellIDs{j}, cellIDsDbx) == 1
                TransDbxArea = [TransDbxArea, IntArea(i,j)];
            end
            if ismember(cellIDs{j}, cellIDsTrans) ==1
                TransTransArea = [TransTransArea, IntArea(i,j)];
            end
            if ismember(cellIDs{j}, cellIDsL) ==1
                TransBarhlArea = [TransBarhlArea, IntArea(i,j)];
            end
        end
    end
end

% remove those pairs with zero overlap
IpsiToIpsi = AlxAlxArea;
IpsiToIpsi = IpsiToIpsi(find(IpsiToIpsi~=0));

IpsiContraToIpsiContra = TransTransArea;
IpsiContraToIpsiContra = IpsiContraToIpsiContra(find(IpsiContraToIpsiContra~=0));

IpsiToContra = [AlxDbxArea';TransDbxArea'];
IpsiToContra = IpsiToContra(find(IpsiToContra~=0));
IpsiToContra = IpsiToContra';

IpsiToUnknown =[AlxBarhlArea';TransBarhlArea'];
IpsiToUnknown = IpsiToUnknown(find(IpsiToUnknown~=0));
IpsiToUnknown = IpsiToUnknown';

h = boxplot([IpsiToIpsi./1e6;IpsiContraToIpsiContra./1e6; IpsiToContra./1e6;IpsiToUnknown./1e6],...
    [ones(size(IpsiToIpsi,1),1); 2*ones(size(IpsiContraToIpsiContra,1),1); 3*ones(size(IpsiToContra,1),1); 4*ones(size(IpsiToUnknown,1),1)],...
    'Notch','off', 'Symbol', 'ko','Colors',[calx;ctrans;cdbx;cbarhl],'OutlierSize', 25);
set(h,{'linew'},{4});
h1 = findobj(gca,'tag','Median');
set(h1,'Color','k');
set(gca,'FontName', 'Arial', 'FontSize', 40, 'LineWidth',4);
box off;

hold on;

plot([ones(size(IpsiToIpsi,1),1); 2*ones(size(IpsiContraToIpsiContra,1),1); 3*ones(size(IpsiToContra,1),1); 4*ones(size(IpsiToUnknown,1),1)]...
    ,[IpsiToIpsi./1e6;IpsiContraToIpsiContra./1e6; IpsiToContra./1e6;IpsiToUnknown./1e6], ...
     'Marker', 'o', 'MarkerFaceColor','none' ,'MarkerSize',25,'LineStyle','none', 'MarkerEdgeColor','k', 'LineWidth', 4);
axis square;

ylabel('Overlapping area (\mum)', 'FontName', 'Arial', 'FontSize', 40);




%% for all cells

AlxPreVol = HeatMapFish(ksize,res, AlxPre,AlxSoma,cellIDsAlx,false);
AlxPostVol = HeatMapFish(ksize,res, AlxPost,AlxSoma,cellIDsAlx,false);
TransPreVol = HeatMapFish(ksize,res, TransPre,TransSoma,cellIDsTrans,false);
TransPostVol = HeatMapFish(ksize,res, TransPost,TransSoma,cellIDsTrans,false);
DbxPostVol = HeatMapFish(ksize,res, DbxPost,DbxSoma,cellIDsDbx,false);
BarhlPostVol = HeatMapFish(ksize,res, BarhlPost,BarhlSoma,cellIDsL,false);

%Ipsi axons with Ipsi dendrites
figure();
AlxAlxArea = dotVol(AlxPreVol,AlxPostVol,AlxSoma, AlxSoma, res);
suptitle('Ipis Axons with Ipsi dendrites');


% IpsiContra axons with IpsiContraDendrites
figure();
TransArea = dotVol(TransPreVol,TransPostVol,TransSoma, TransSoma, res);
suptitle('IpisContra Axons with IpsiContra dendrites');

% IpsiAxons with Contra dendrites
figure();
subplot(1,2,1);
AlxDbxArea = dotVol(AlxPreVol,DbxPostVol,AlxSoma, DbxSoma, res);
subplot(1,2,2);
TransDbxArea = dotVol(TransPreVol,DbxPostVol,TransSoma, DbxSoma, res);
suptitle('Ipis Axons with Contra dendrites');


% IpsiAxons with Unknown dendrites
figure();
subplot(1,2,1);
AlxBarhlArea = dotVol(AlxPreVol,BarhlPostVol,AlxSoma, BarhlSoma, res);
subplot(1,2,2);
TransBarhlArea = dotVol(TransPreVol,BarhlPostVol,TransSoma, BarhlSoma, res);
suptitle('Ipis Axons with Unkown dendrites');


figure();
bar([ AlxAlxArea/1e6,TransArea/1e6,(AlxDbxArea+TransDbxArea)/1e6, (AlxBarhlArea+TransBarhlArea)/1e6] ,...
    'color',{calx, ctrans, cdbx, cbarhl});




