load ConnMatrixPre.mat
% load connMatCKT_withNegativeExamples.mat
% load MatOrder_withNegativeExamples.mat
load AllCells.mat

newNeurons  = [77061,76392,76422,77062,76475,76571,76455,76886,76606,76404,77050,...
    77044,77079,76897,76892,76411,76863,76349,76386,76414,76412,76415,76837,76859,...
    76541,76352,76472,76879,76400,76831,76909,76496,76875,77075,76406,77063,76918,...
    76894,77068,76493,76612,76887,76862,77056,76463,76410,76614,76470,77080,76402,...
    76354,76497,77091,77032,76568,76607,76556,76893,77089,77053,76401,76869,76471,...
    76912,77072,76391,76388,76418,76492,76468,77066,76552,76382,77048,76903,76555,...
    76531,76870,77052,76452,76896,76844,77070,76610,76880,76907,76898,76485,77094,...
    76851,76551,77098,76539,76384,77031,77036,76889,77092,76403,77096,77086,76865,...
    76530,76466,77069,77088,76569,76293,76489,76405,76454,77060,76560,76563,77081,...
    76462,76465,76393,76464,76387,77085,76868,76866,76890,76860,76348,76381,76827,...
    76474,77074,76856,77095,76538,77042,76491,77054,77078,77084,76564,76567,76902,...
    76608,76451,77093,76570,76861,77058,76486,76457,76845,76885,76871,76490,76419,...
    76473,76833,76609,76469,76553,76857,77082,76456,77083,76613,76872,77090,76407,...
    77059,76565,77041,76899,77071,77076,76554,77039,76453,76858,76409,77057,76390,...
    76561,76537,76864,76834,76835,76900,76964,76972,76948,76952,76936,76940,76970,...
    76933,77022,76943,76969,76953,77024,76957,76942,76971,76938,76923,77021,77020,...
    76955,76921,76924,76929,76950,76931,76922,77029,77027,76926,76956,77025,76930,76968,76928];


if ismac
    addpath(genpath('/Users/ashwin/Documents/'));
else
    addpath(genpath('/usr/people/ashwinv/seungmount/research/Ashwin/BCT'));
end

[m,n,v] = find(ConnMatrixPre); % m -postSynapse,row, n - column, v - number of synapses,value
%[m,n,v] = find(connMatCKT)
% for i = 1:size(intPartners,1)
%     AllCellNames{i} = num2str(intPartners(i));
% end
AllCellNames = AllCells;
%AllCellNames = MatOrder;

G = digraph(n',m',v',size(AllCellNames,1));                 % digraph(source, target, weight, numNodes)
%G.Nodes.Names = AllCellNames';
G.Nodes.Names = num2str(AllCellNames);  % nodeNames

% for i = 1:size(AllCellNames)
%     temp = find(AllCellNames(i) == mlOrdered);
%     if ~isempty(temp)
%         AllCellManualID(i) = mlOrderedTypes(temp);
%     else
%         AllCellManualID(i) = cellstr(sprintf('%d',AllCellNames(i)));
%     end
% end

%%
inDeg = indegree(G);
%notConnected = find(inDeg < 2);                       % weakly connected nodes
%inDeg(notConnected) = [];                             % drop them from deg
%Gin = rmnode(G, notConnected);
[~, ranking] = sort(inDeg,'descend');                 % get ranking by degree
top3in = G.Nodes.Names(ranking(1:10));                  % get top 3 node names

outDeg = outdegree(G);
%notConnected = find(outDeg < 2);                       % weakly connected nodes
%outDeg(notConnected) = [];
%Gout = rmnode(G, notConnected);
[~, ranking] = sort(outDeg,'descend');                 % get ranking by degree
top3out = G.Nodes.Names(ranking(1:10));                  % get top 3 node names

%%

%Adj = adjacency(G,'Weights',G.Edges.Weight);
nn = numnodes(G);
[s,t] = findedge(G);
A = sparse(s,t,G.Edges.Weight,nn,nn);
%Inc = incidence(G);
%Lap = laplacian(G);

%% try clustering
% Iterations = 10;
% parfor i = 1:1:Iterations
%     [Mmany(:,i),Q(:,i)] = community_louvain(A);
%     disp(i)
% end

Q0 = -1; Q1 = 0; % initialize modularity values
i =1 ; 
while abs(Q1-Q0)>1e-5;           % while modularity increases
    Q0 = Q1;                % perform community detection
    [M, Q1] = community_louvain(A,[],1:size(A,1));
    plot(i,Q1-Q0,'o');
    hold on;
    drawnow;
    i = i+1;
end
%%
[Vsort,idx] = sort(M);
figure; cspy(A(idx,idx),'Colormap',colorcet('R3'),'Levels',255,'MarkerSize',15);
blocks = find(diff(round(Vsort)));
hold on;
for i = 1:length(blocks)
    line([0,size(ConnMatrixPre,1)],[blocks(i),blocks(i)],'color','k');
    line([blocks(i),blocks(i)],[0,size(ConnMatrixPre,1)],'color','k');
end

AllCellNamesOrdered = AllCellNames(idx);
load allApprovedIntegrators.mat



% locater78Integrators = find(ismember(AllCellNamesOrdered,ALX.cellIDs));
% locater456Integrators = find(ismember(AllCellNamesOrdered,[SaccadicProjectingToABDexclusively.cellID'; SaccadicProjectingToABDiexclusively.cellID']));
% locater456Bino = find(ismember(AllCellNamesOrdered,SaccadicProjectingToABD_ABDi.cellID));
locateIntegrators = find(ismember(AllCellNamesOrdered,allApprovedIntegrators));
locateNewNeurons = find(ismember(AllCellNamesOrdered,newNeurons));

% text(locater78Integrators,repmat(-10,1,numel(locater78Integrators)), sprintf('\\downarrow'), 'HorizontalAlignment','center', 'FontWeight','bold','Color','r');
% text(locater456Integrators,repmat(-10,1,numel(locater456Integrators)), sprintf('\\downarrow'), 'HorizontalAlignment','center', 'FontWeight','bold','Color','b');
% text(locater456Bino,repmat(-10,1,numel(locater456Bino)), sprintf('\\downarrow'), 'HorizontalAlignment','center', 'FontWeight','bold','Color','g');
text(locateNewNeurons,repmat(-10,1,numel(locateNewNeurons)), sprintf('\\downarrow'), 'HorizontalAlignment','center', 'FontWeight','bold','Color','k');
text(locateIntegrators,repmat(-10,1,numel(locateIntegrators)), sprintf('\\downarrow'), 'HorizontalAlignment','center', 'FontWeight','bold','Color','r');




% set(gca,'YTick',1:size(A),'YTickLabel',AllCellNames(idx),'XTick',1:size(A),...
%     'XTickLabel',AllCellNames(idx),'XTickLabelRotation',45,'XAxisLocation','top');

%%
% find int cells in idx sapce
[c,ia,ib] = intersect(1:18,idx);

figure; subplot(1,2,1); cspy(IntConnMatrixPre,'Colormap',colorcet('L8'),'Levels',255,'MarkerSize',15);
text(repmat(-50,1,18), 1:18, sprintf('\\rightarrow'), 'HorizontalAlignment','center', 'FontWeight','bold');
text(1:18,repmat(-50,1,18), sprintf('\\downarrow'), 'HorizontalAlignment','center', 'FontWeight','bold');
box off;

subplot(1,2,2); cspy(A(idx,idx),'Colormap',colorcet('L8'),'Levels',255,'MarkerSize',15);
text(repmat(-50,1,size(ib,1)), ib, sprintf('\\rightarrow'), 'HorizontalAlignment','center', 'FontWeight','bold');
text( ib,repmat(-50,1,size(ib,1)), sprintf('\\downarrow'), 'HorizontalAlignment','center', 'FontWeight','bold')
box off;

blocks = find(diff(Vsort));
for i = 1:length(blocks)
    line([0,size(IntConnMatrixPre,1)],[blocks(i),blocks(i)],'color','k');
    line([blocks(i),blocks(i)],[0,size(IntConnMatrixPre,1)],'color','k');
end

