function [] = AxonalTree(AxnTree,treeno, cellIDs, Display)
% AXONALTREE is to extract all the nodes of the dendrites of a tree
%   Tree is the tree whos dendrites are being extracted
%   treeno is the ID associated with the tree
%   cellIDs are the IDs of all cells in the dataset
%   col is a string assocaited with the color. eg. 'b'
%   Display true or false

col = [rand rand rand];
load CellAxons.mat
axonTree = [];
inducingNodes = eval([cellIDs{treeno},'_axon']);
subplot(1,2,1)

for kk=1:numel(AxnTree)
    children = AxnTree{kk}{2};
    for mm = 1:numel(children)
        if ismember(kk,inducingNodes)
            Axtempx=[AxnTree{children(mm)}{3}(1); AxnTree{children(mm)}{4}{1}(:,1); AxnTree{kk}{3}(1)];
            Axtempy=[AxnTree{children(mm)}{3}(2); AxnTree{children(mm)}{4}{1}(:,2); AxnTree{kk}{3}(2)];
            Axtempz=-[AxnTree{children(mm)}{3}(3); AxnTree{children(mm)}{4}{1}(:,3); AxnTree{kk}{3}(3)];
            if Display == true
                plot3(Axtempx,Axtempy,Axtempz,'color',col,'lineWidth',2);
                hold on;
            end
            axonTree = [axonTree; Axtempx Axtempy Axtempz];
        else
            continue;
        end
    end
end
% xy view
view(-180,0);
set (gca,'XTick',[], 'YTick',[],'ZTick', [], 'ZLim',[-60000, -0],'XLim',[20000 , 140000]);

% plot histogram
subplot(1,2,2)
h = histogram(-1*axonTree(:,3)/1000, 'BinWidth',1,'Visible','off');%,'Normalization','probability');
x= h.Values;
y = h.BinEdges;
plot([0,x],y,'-','Color',col);
%plot(h.Values,1:2:60-1);
hold all;
set(gca,'YLim',[0 ,60],'YDir','reverse');
box off;


