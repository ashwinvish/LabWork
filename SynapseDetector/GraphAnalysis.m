load('IntConnMatrix.mat');
load('IntPartners.mat');
addpath(genpath('/Users/admin/Documents/ComDetTBv091'));
addpath(genpath('/Users/admin/Documents/BCT'));
[m,n,v] = find(IntConnMatrixPre); % m -postSynapse, n - presynapse, v - number of synapses
for i = 1:size(intPartners,1)
    AllCellNames{i} = num2str(intPartners(i));
end
G = digraph(n',m',v',size(AllCellNames,2));                 % digraph(source, target, weight, numNodes)
G.Nodes.Names = AllCellNames';                              % nodeNames
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
Iterations = 10000;
for i = 1:1:Iterations
[M(:,i),Q(:,i)] = community_louvain(A);
end
%%
[Vsort,idx] = sort(M);
figure; cspy(A(idx,idx),'Colormap',colorcet('L8'),'Levels',255,'MarkerSize',15);
blocks = find(diff(Vsort));
hold on;
for i = 1:length(blocks)
    line([0,size(IntConnMatrixPre,1)],[blocks(i),blocks(i)],'color','k');
    line([blocks(i),blocks(i)],[0,size(IntConnMatrixPre,1)],'color','k');
end
%%
% find int cells in idx sapce
[c,ia,ib] = intersect(1:18,idx);

figure; subplot(1,2,1); cspy(IntConnMatrixPre,'Colormap',colorcet('L8'),'Levels',255,'MarkerSize',15);
text(repmat(-50,1,18), 1:18, sprintf('\\rightarrow'), 'HorizontalAlignment','center', 'FontWeight','bold');
text(1:18,repmat(-50,1,18), sprintf('\\downarrow'), 'HorizontalAlignment','center', 'FontWeight','bold');
box off;

subplot(1,2,2); cspy(IntConnMatrixPre(idx,idx),'Colormap',colorcet('L8'),'Levels',255,'MarkerSize',15);
text(repmat(-50,1,size(ib,1)), ib, sprintf('\\rightarrow'), 'HorizontalAlignment','center', 'FontWeight','bold');
text( ib,repmat(-50,1,size(ib,1)), sprintf('\\downarrow'), 'HorizontalAlignment','center', 'FontWeight','bold')
box off;

blocks = find(diff(Vsort));
for i = 1:length(blocks)
    line([0,size(IntConnMatrixPre,1)],[blocks(i),blocks(i)],'color','k');
    line([blocks(i),blocks(i)],[0,size(IntConnMatrixPre,1)],'color','k');
end

