
AlxPostLength = [];
DbxPostLength = [];
BarhlPostLength = [];
for i = 1:length(cellIDs)
    if ismember(cellIDs{i},cellIDsAlx)==1
        [AlxPostLength] = [AlxPostLength; allLengthToPostNode{i}];
    elseif ismember(cellIDs{i},cellIDsDbx)==1
        [DbxPostLength] = [DbxPostLength; allLengthToPostNode{i}];
    else
        [BarhlPostLength] = [BarhlPostLength; allLengthToPostNode{i}];
    end
end

figure();
histogram(AlxPostLength/1000,'FaceColor',calx,'Normalization','probability','BinWidth', 10);
hold on;
histogram(DbxPostLength/1000,'FaceColor', cdbx,'Normalization','probability','BinWidth', 10);
histogram(BarhlPostLength/1000,'FaceColor', cbarhl,'Normalization','probability','BinWidth', 10);

figure();
subplot(2,3,2);
h1 = histogram(allPostSynapticLength/1000,'FaceColor', [0.9,0,0]); % dimensions in microns
h1.Normalization = 'probability';
h1.BinWidth = 10;
h1.FaceAlpha = 0.2;
hold on;
h2 = histogram(allPreSynapticLength/1000,'FaceColor',[0,0.8,0]); % dimensions in microns
h2.Normalization = 'probability';
h2.BinWidth = 10;
h2.FaceAlpha = 0.2;
axis square
box off;
ylabel('Probability');
xlabel('Pathlength in \mum');
%set(gca,'XLim', [0,250],'FontName','Arial', 'FontSize', 40, 'LineWidth',2);

%figure();
subplot(2,3,4);
h1 = histogram(AlxPostLength/1000,'FaceColor',[0.9,0,0], 'Normalization','probability','BinWidth', 10);
h1.FaceAlpha = 0.2;
hold on;
h2 = histogram(allPreSynapticLength/1000,'FaceColor',[0,0.8,0],'Normalization','probability', 'BinWidth',10); % dimensions in microns
h2.FaceAlpha = 0.2;
axis square
box off;
ylabel('Probability');
xlabel('Pathlength in \mum');
%set(gca,'XLim', [0,250],'FontName','Arial', 'FontSize', 40, 'LineWidth',2);
title('Alx Cells','FontName','Arial');

%figure();
subplot(2,3,5);
h1 = histogram(DbxPostLength/1000,'FaceColor',[0.9,0,0],'Normalization','probability','BinWidth', 10);
h1.FaceAlpha = 0.2;
hold on;
h2 = histogram(allPreSynapticLength/1000,'FaceColor',[0,0.8,0],'Normalization','probability','BinWidth',10); % dimensions in microns
h2.FaceAlpha = 0.2;
axis square
box off;
ylabel('Probability');
xlabel('Pathlength in \mum');
%set(gca,'XLim', [0,250],'FontName','Arial', 'FontSize', 40, 'LineWidth',2);
title('Dbx Cells','FontName','Arial');


%figure();
subplot(2,3,6);
h1 = histogram(BarhlPostLength/1000,'FaceColor',[0.9,0,0],'Normalization','probability','BinWidth', 10);
h1.FaceAlpha = 0.2;
hold on;
h2 = histogram(allPreSynapticLength/1000,'FaceColor',[0,0.8,0],'Normalization','probability', 'BinWidth',10); % dimensions in microns
h2.FaceAlpha = 0.2;
axis square
box off;
ylabel('Probability');
xlabel('Pathlength in \mum');
%set(gca,'XLim', [0,250],'FontName','Arial', 'FontSize', 40, 'LineWidth',2);
title('Barhl Cells','FontName','Arial');




