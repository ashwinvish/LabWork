AxonInitiation{1} = NaN;
AxonInitiation{2} = [62850,120545,6255];
AxonInitiation{3} = [55725,104395,12330];
AxonInitiation{4} = [89520,174240,15750; 36535,166090,46305];
AxonInitiation{5} = [77820,162670,45315]
AxonInitiation{6} = [65180,195125,18900;59795,194995,31005 ]
AxonInitiation{7} = [60315,171125,26145]
AxonInitiation{8} = [73095,116140,23175]
AxonInitiation{9} = [63115,116325,14670]
AxonInitiation{10} = [62175,117240,13005]
AxonInitiation{11} = [61015,125075,12735]
AxonInitiation{12} = [58120,104615,16875]
AxonInitiation{13} = [44945,189525,27450]
AxonInitiation{14} = NaN
AxonInitiation{15} = [60220,167570,29070]
AxonInitiation{16} = [43015,190305,43515];
AxonInitiation{17} = NaN
AxonInitiation{18} = NaN
AxonInitiation{19} = NaN
AxonInitiation{20} = NaN
AxonInitiation{21} = [54245,157650,17730]
AxonInitiation{22} = [63500,176335,47925]

%%

AlxAxonInitialLength = [];
TransAxonInitialLength = [];
DbxAxonInitialLength = [];
BarhlAxonInitiationLength = [];

for i = 1:numel(cellIDs)
    if ismember(cellIDs{i}, cellIDsAlx) ==1
        temp = findPathLength_new([cellIDs{i} , '_WithTags.swc'],allTrees{i},[5,5,45],AxonInitiation{i});
        AlxAxonInitialLength =  [AlxAxonInitialLength; temp];
        AxonInitiationLength(i) = temp(1);
        clear temp;
    elseif ismember(cellIDs{i}, cellIDsTrans) ==1
        temp = findPathLength_new([cellIDs{i} , '_WithTags.swc'],allTrees{i},[5,5,45],AxonInitiation{i});
        TransAxonInitialLength =  [TransAxonInitialLength; temp];
        AxonInitiationLength(i) = temp(1);
        clear temp;
    elseif ismember(cellIDs{i}, cellIDsDbx) ==1
        temp = findPathLength_new([cellIDs{i} , '_WithTags.swc'],allTrees{i},[5,5,45],AxonInitiation{i});
        DbxAxonInitialLength =  [DbxAxonInitialLength; temp];
        AxonInitiationLength(i) = temp(1);
        clear temp;
    else
        BarhlAxonInitiationLength = 0;
        AxonInitiationLength(i) = 0;
    end 
    
end

plot(ones(1,size(AlxAxonInitialLength,1)), AlxAxonInitialLength/1000, 'o', 'MarkerSize', 35, 'MarkerFaceColor', calx, 'MarkerEdgeColor', 'k');
hold on;
plot(2*ones(1,size(TransAxonInitialLength,1)), TransAxonInitialLength/1000, 'o', 'MarkerSize', 35, 'MarkerFaceColor', ctrans, 'MarkerEdgeColor', 'k');
plot(3*ones(1,size(DbxAxonInitialLength,1)), DbxAxonInitialLength/1000, 'o', 'MarkerSize', 35, 'MarkerFaceColor', cdbx, 'MarkerEdgeColor', 'k');
set(gca, 'XLim',[0.5,3.5],'LineWidth',2, 'FontName', 'Arial', 'FontSize', 40);
ylabel('Axon initiation (\mum)', 'FontName', 'Arial', 'FontSize', 40);
box off;
axis square;










