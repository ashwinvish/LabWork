function BranchOrderVisualizer(tree,highlightedNodes,inducingNodes,specialNodes,newFigure,validNodes)
%UNTITLED14 Summary of this function goes here
%   Detailed explanation goes here

rndclr = jet(length(unique(inducingNodes))+1);
symCell={'o','s','v','x','d','*'};
synapseColor = [[0.9,0,0]; [0,0.8,0]];
MEdgeColor = [[0.5,0,0]; [0,0.5,0]];
hold on;

if nargin < 7
    validNodes = [1:numel(tree)];
    if nargin < 6
        newFigure = true;
        if nargin < 5
            specialNodes = [];
            if nargin < 3
                inducingNodes = 1:numel(tree);
                if nargin < 2
                    highlightedNodes = [];
                end
            end
        end
    end
end


for kk=1:numel(tree)
    children = tree{kk}{2};
    for mm=1:numel(children)
        tempx=[tree{children(mm)}{3}(1); tree{children(mm)}{4}{1}(:,1); tree{kk}{3}(1)];
        tempy=[tree{children(mm)}{3}(2); tree{children(mm)}{4}{1}(:,2); tree{kk}{3}(2)];
        tempz=-[tree{children(mm)}{3}(3); tree{children(mm)}{4}{1}(:,3); tree{kk}{3}(3)];
        h =  plot3(tempx,tempy,tempz,'color',rndclr(inducingNodes(kk)+1,:),'lineWidth',2);
    end
end


for kk=1:numel(highlightedNodes)
    plot3(tree{highlightedNodes(kk)}{3}(1),tree{highlightedNodes(kk)}{3}(2),-tree{highlightedNodes(kk)}{3}(3),'Marker','o','MarkerSize' , 10, 'LineWidth', 3, 'MarkerFaceColor','k','MarkerEdgeColor' , 'k' )
end

if ~isempty(specialNodes) % special nodes to be marked and hilighted
    
    for mm = 1
        for kk=1:size(specialNodes{mm}, 1)
            hSyn1 = plot3(specialNodes{mm}(kk,1),specialNodes{mm}(kk,2),-specialNodes{mm}(kk,3),'Marker',symCell{rem(mm-1,numel(symCell))+1}, 'MarkerSize' , 5, 'LineWidth', 0.1 , 'MarkerFaceColor',synapseColor(mm,:) , 'MarkerEdgeColor' , MEdgeColor(mm,:));
        end
    end
    
    for mm = 2
        for kk=1:size(specialNodes{mm}, 1)
            hSyn2 = plot3(specialNodes{mm}(kk,1),specialNodes{mm}(kk,2),-specialNodes{mm}(kk,3),'Marker',symCell{rem(mm-1,numel(symCell))+1}, 'MarkerSize' , 10, 'LineWidth', 0.1 , 'MarkerFaceColor', synapseColor(mm,:) , 'MarkerEdgeColor' , MEdgeColor(mm,:));
            
        end
    end
    
end
    
box on;
axis([ 20000 140000 60000 250000 -60000 0]);
plot( [20000, 40000], [70000, 70000],'-k' ) % insert 20um sclaebar
daspect([1 1 1]); % make aspect ratio [1 1 1]
set (gca,'XTick',[], 'YTick',[],'ZTick', [], 'Ydir','reverse');
view([-180,90]); % xy view
%colorbar ;
%set(gca,'CLim',[min(inducingNodes) max(inducingNodes)]);
end



