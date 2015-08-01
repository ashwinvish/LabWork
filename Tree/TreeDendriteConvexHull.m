function [ ConHull,ConVol, X,Y,Z,htri1] = TreeDendriteConvexHull(tree,highlightedNodes,inducingNodes,specialNodes,newFigure,colorString,HullColor,validNodes,pixelUnits )
%TreeConvexHull Plots the convex hull of tree, similar to treeVisualizer.m

relativeRes = [5 5 45]; % in nm
rndclr = colorString;
% 
% list = lineage(tree,1,[]);
% list(inducingNodes) = [];
% inducingNodes = list;

symCell={'o','s','v','x','d','*'};
synapseColor = {'g', 'r'}; % red - presynaptic; g- postsynaptic
%synapseColor = {'green', 'red'}; % red - presynaptic; g- postsynaptic

if nargin < 9
    pixelUnits = false;
    if nargin < 8
        validNodes = [1:numel(tree)];
        if nargin <7
            HullColor= 'cyan';
            if nargin < 6
                rndclr = colorString;
                if nargin < 5
                    newFigure = true;
                    if nargin < 3
                        inducingNodes = 1:numel(tree);
                        if nargin < 2
                            highlightedNodes = [];
                        end
                    end
                end
            end
        end
    end
    
    if newFigure
        figure;hold;
    else
        hold on;
    end
    
    X = [];  Y = []; Z = [];
  
    for kk=1:numel(tree)
        children = tree{kk}{2};
        
        for mm=1:numel(children)
            if pixelUnits
                tempx=[tree{children(mm)}{3}(1); tree{children(mm)}{4}{1}(:,1); tree{kk}{3}(1)]*relativeRes(1);
                tempy=[tree{children(mm)}{3}(2); tree{children(mm)}{4}{1}(:,2); tree{kk}{3}(2)]*relativeRes(2);
                tempz=-[tree{children(mm)}{3}(3); tree{children(mm)}{4}{1}(:,3); tree{kk}{3}(3)]*relativeRes(3);
            else
                tempx=[tree{children(mm)}{3}(1); tree{children(mm)}{4}{1}(:,1); tree{kk}{3}(1)];
                tempy=[tree{children(mm)}{3}(2); tree{children(mm)}{4}{1}(:,2); tree{kk}{3}(2)];
                tempz=-[tree{children(mm)}{3}(3); tree{children(mm)}{4}{1}(:,3); tree{kk}{3}(3)];
            end
            
            if ismember(kk,inducingNodes) && ismember(children(mm),inducingNodes)
                plot3(tempx,tempy,tempz,'color', [1,0.3,0]);
                 X = vertcat(X,tempx); Y= vertcat(Y,tempy); Z = vertcat(Z,tempz);
            else 
                if ismember(kk,validNodes)
                    plot3(tempx,tempy,tempz,'color', rndclr{1});
                end
            end
        end 
    end
    % plot convexhull of tree
    [ConHull ConVol] = convhull(X,Y,Z,'simplify', true);
    htri1 = trimesh(ConHull,X,Y,Z,'faceColor',HullColor,'FaceAlpha',0.2,'EdgeColor','black','EdgeAlpha',0.1);
       
    
    for kk=1:numel(highlightedNodes)
        if pixelUnits
            plot3(tree{highlightedNodes(kk)}{3}(1)*relativeRes(1),tree{highlightedNodes(kk)}{3}(2)*relativeRes(2),-tree{highlightedNodes(kk)}{3}(3)* relativeRes(3),'Marker','o','MarkerSize',10 , 'MarkerFaceColor',rndclr{1},'MarkerEdgeColor','k' )
        else
            plot3(tree{highlightedNodes(kk)}{3}(1),tree{highlightedNodes(kk)}{3}(2),-tree{highlightedNodes(kk)}{3}(3),'Marker','o','MarkerSize' , 10,  'MarkerFaceColor',rndclr{1},'MarkerEdgeColor' , 'k' )
        end
    end
    
    for mm = 1: numel(specialNodes)
        for kk=1:size(specialNodes{mm}, 1)
            if pixelUnits
                plot3(specialNodes{mm}(kk,1)*relativeRes(1),specialNodes{mm}(kk,2)*relativeRes(2),-specialNodes{mm}(kk,3)*relativeRes(3),'Marker',symCell{rem(mm-1,numel(symCell))+1}, 'MarkerSize' , 4,  'MarkerFaceColor', synapseColor{mm} , 'MarkerEdgeColor' , 'none');
            else
                plot3(specialNodes{mm}(kk,1),specialNodes{mm}(kk,2),-specialNodes{mm}(kk,3),'Marker',symCell{rem(mm-1,numel(symCell))+1}, 'MarkerSize' , 4,  'MarkerFaceColor', synapseColor{mm} , 'MarkerEdgeColor' , 'none');
            end
        end
    end

    daspect([1 1 1]);
    axis vis3d;
    axis([20000 140000  60000 250000 -60000 0]);

    plot( [20000, 40000], [70000, 70000],'-k' ) % insert 20um sclaebar

    
    box off;
    XColor = [1,1,1]; YColor = [1,1,1];
    set (gca, 'XTick',[], 'YTick',[],'ZTick', [], 'Ydir','reverse');
    view([-180,90]); 
end

