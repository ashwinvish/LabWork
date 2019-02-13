load('IntConnMatrixPre.mat');
load('IntPartners.mat');
if ismac
    addpath(genpath('/Users/admin/Documents/ComDetTBv091'));
    addpath(genpath('/Users/admin/Documents/BCT'));
else
    addpath(genpath('/usr/people/ashwinv/seungmount/research/Ashwin/BCT'));
end

[m,n,v] = find(IntConnMatrixPre); % m -postSynapse,row, n - column, v - number of synapses,value
% for i = 1:size(intPartners,1)
%     AllCellNames{i} = num2str(intPartners(i));
% end
AllCellNames = intPartners;
G = digraph(n',m',v',size(AllCellNames,1));                 % digraph(source, target, weight, numNodes)
%G.Nodes.Names = AllCellNames';
G.Nodes.Names = num2str(intPartners);  % nodeNames
for i = 1:size(AllCellNames)
    temp = find(AllCellNames(i) == mlOrdered);
    if ~isempty(temp)
        AllCellManualID(i) = mlOrderedTypes(temp);
    else
        AllCellManualID(i) = cellstr(sprintf('%d',AllCellNames(i)));
    end
end

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
Iterations = 10;
parfor i = 1:1:Iterations
    [Mmany(:,i),Q(:,i)] = community_louvain(A);
    disp(i)
end

% Q0 = -1; Q1 = 0;            % initialize modularity values
% while Q1-Q0>1e-5;           % while modularity increases
%     Q0 = Q1;                % perform community detection
%     [M, Q1] = community_louvain(A,[],1:size(A,1));
% %     plot(Q1-Q0);
%     %hold on;
% end
%%
[Vsort,idx] = sort((Mmany(:,10)));
figure; cspy(A(idx,idx),'Colormap',colorcet('L8'),'Levels',255,'MarkerSize',15);
blocks = find(diff(round(Vsort)));
hold on;
for i = 1:length(blocks)
    line([0,size(IntConnMatrixPre,1)],[blocks(i),blocks(i)],'color','k');
    line([blocks(i),blocks(i)],[0,size(IntConnMatrixPre,1)],'color','k');
end


set(gca,'YTick',1:size(A),'YTickLabel',AllCellManualID(idx),'XTick',1:size(A),...
    'XTickLabel',AllCellManualID(idx),'XTickLabelRotation',45,'XAxisLocation','top');

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

