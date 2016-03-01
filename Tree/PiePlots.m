% Pie plots
DenLengthPie = [sum(AlxDenLenght), sum(TransDenLength), sum(DbxDenLength), sum(BarhlDenLength)];
figure;
%labels = {'group1','group2','group3','group4'};
h = pie(DenLengthPie);
set(h(1),'FaceColor',calx, 'LineWidth',2);
set(h(2),'FontName', 'Arial', 'FontSize', 40);
set(h(3),'FaceColor',ctrans, 'LineWidth',2);
set(h(4),'FontName', 'Arial', 'FontSize', 40);
set(h(5),'FaceColor',cdbx, 'LineWidth',2);
set(h(6),'FontName', 'Arial', 'FontSize', 40);
set(h(7),'FaceColor',cbarhl, 'LineWidth',2);
set(h(8),'FontName', 'Arial', 'FontSize', 40);
set(gcf,'color','w');
title('Total dendritic length','FontName', 'Arial', 'FontSize', 40) ;

%Normalized Plot
DenLengthPie = [sum(AlxDenLenght)/numel(AlxDenLenght), sum(TransDenLength)/numel(TransDenLength), sum(DbxDenLength)/numel(DbxDenLength), sum(BarhlDenLength)/numel(BarhlDenLength)];
figure;
%labels = {'group1','group2','group3','group4'};
h = pie(DenLengthPie);
set(h(1),'FaceColor',calx, 'LineWidth',2);
set(h(2),'FontName', 'Arial', 'FontSize', 40);
set(h(3),'FaceColor',ctrans, 'LineWidth',2);
set(h(4),'FontName', 'Arial', 'FontSize', 40);
set(h(5),'FaceColor',cdbx, 'LineWidth',2);
set(h(6),'FontName', 'Arial', 'FontSize', 40);
set(h(7),'FaceColor',cbarhl, 'LineWidth',2);
set(h(8),'FontName', 'Arial', 'FontSize', 40);
set(gcf,'color','w');
title('Normalized dendritic length','FontName', 'Arial', 'FontSize', 40) ;

figure();

AxLengthPie = [sum(AlxAxnLength), sum(TransAxnLength), sum(DbxAxnLength), sum(BarhlAxnLenght)];
h = pie(AxLengthPie);
set(h(1),'FaceColor',calx, 'LineWidth',2);
set(h(2),'FontName', 'Arial', 'FontSize', 40);
set(h(3),'FaceColor',ctrans);
set(h(4),'FontName', 'Arial', 'FontSize', 40);
set(h(5),'FaceColor',cdbx, 'LineWidth',2);
set(h(6),'FontName', 'Arial', 'FontSize', 40);
set(gcf,'color','w');
title('Total axonal length','FontName', 'Arial', 'FontSize', 40) ;

figure();

AxLengthPie = [sum(AlxAxnLength)/numel(AlxAxnLength), sum(TransAxnLength)/ numel(TransAxnLength), sum(DbxAxnLength)/numel(DbxAxnLength), sum(BarhlAxnLenght)/numel(BarhlAxnLenght)];
h = pie(AxLengthPie);
set(h(1),'FaceColor',calx, 'LineWidth',2);
set(h(2),'FontName', 'Arial', 'FontSize', 40);
set(h(3),'FaceColor',ctrans);
set(h(4),'FontName', 'Arial', 'FontSize', 40);
set(h(5),'FaceColor',cdbx, 'LineWidth',2);
set(h(6),'FontName', 'Arial', 'FontSize', 40);
set(gcf,'color','w');
title('Normalized axonal length','FontName', 'Arial', 'FontSize', 40) ;

figure();

SynapsePre = [sum(AlxPre),sum(TransPre),0,0];
h = pie(SynapsePre);
set(h(1),'FaceColor',calx, 'LineWidth',2);
set(h(2),'FontName', 'Arial', 'FontSize', 40);
set(h(3),'FaceColor',ctrans, 'LineWidth',2);
set(h(4),'FontName', 'Arial', 'FontSize', 40);
set(gcf,'color','w');
title('Total presynaptic sites','FontName', 'Arial', 'FontSize', 40) ;

SynapsePre = [sum(AlxPre)/numel(cellIDsAlx),sum(TransPre)/numel(cellIDsTrans),0,0];
h = pie(SynapsePre);
set(h(1),'FaceColor',calx, 'LineWidth',2);
set(h(2),'FontName', 'Arial', 'FontSize', 40);
set(h(3),'FaceColor',ctrans, 'LineWidth',2);
set(h(4),'FontName', 'Arial', 'FontSize', 40);
set(gcf,'color','w');
title('Normalized presynaptic sites','FontName', 'Arial', 'FontSize', 40) ;

figure();
SynapsePost = [sum(AlxPost), sum(TransPost), sum(DbxPost), sum(BarhlPost)];
h = pie(SynapsePost);
set(h(1),'FaceColor',calx, 'LineWidth',2);
set(h(2),'FontName', 'Arial', 'FontSize', 40);
set(h(3),'FaceColor',ctrans, 'LineWidth',2);
set(h(4),'FontName', 'Arial', 'FontSize', 40);
set(h(5),'FaceColor',cdbx, 'LineWidth',2);
set(h(6),'FontName', 'Arial', 'FontSize', 40);
set(h(7),'FaceColor',cbarhl, 'LineWidth',2);
set(h(8),'FontName', 'Arial', 'FontSize', 40);
set(gcf,'color','w');
title('Total postsynaptic sites','FontName', 'Arial', 'FontSize', 40) ;

figure();
SynapsePost = [sum(AlxPost)/numel(cellIDsAlx), sum(TransPost)/numel(cellIDsTrans), sum(DbxPost)/numel(cellIDsDbx), sum(BarhlPost)/numel(cellIDsL)];
h = pie(SynapsePost);
set(h(1),'FaceColor',calx, 'LineWidth',2);
set(h(2),'FontName', 'Arial', 'FontSize', 40);
set(h(3),'FaceColor',ctrans, 'LineWidth',2);
set(h(4),'FontName', 'Arial', 'FontSize', 40);
set(h(5),'FaceColor',cdbx, 'LineWidth',2);
set(h(6),'FontName', 'Arial', 'FontSize', 40);
set(h(7),'FaceColor',cbarhl, 'LineWidth',2);
set(h(8),'FontName', 'Arial', 'FontSize', 40);
set(gcf,'color','w');
title('Normalized postsynaptic sites','FontName', 'Arial', 'FontSize', 40) ;