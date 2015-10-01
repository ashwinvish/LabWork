function treeVisualizer(tree,highlightedNodes,inducingNodes,specialNodes,newFigure,colorString,validNodes,pixelUnits)
%relativeRes = [5 5 45]; % in nm
rndclr = colorString;
symCell={'o','s','v','x','d','*'};
%synapseColor = {[0 0.4 0.4], [0.2 0 0.8]}; % red - presynaptic; g- postsynaptic [1 0.2 0.2]
synapseColor = [[0,1,0];[1,0,0]];
if nargin < 8
    pixelUnits = false;
    if nargin < 7
        validNodes = [1:numel(tree)];
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
    figure;
    %figure('units','normalized','outerposition',[0 0 1 1]);
    hold;
else
    hold on;
end

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
            h =  plot3(tempx,tempy,tempz,'color',rndclr{2},'lineWidth',1);
            %h.Color(4) = 0.4; %  transparency 0-1
        else
            if ismember(kk,validNodes)
                h = plot3(tempx,tempy,tempz,'color', rndclr{1},'lineWidth',0.5);
                % h.Color(4) = 0.4;
            end
        end
    end
end

for kk=1:numel(highlightedNodes)
    if pixelUnits
        plot3(tree{highlightedNodes(kk)}{3}(1)*relativeRes(1),tree{highlightedNodes(kk)}{3}(2)*relativeRes(2),-tree{highlightedNodes(kk)}{3}(3)* relativeRes(3),'Marker','o','MarkerSize',10 , 'MarkerFaceColor',rndclr{1},'MarkerEdgeColor','k' )
    else
        plot3(tree{highlightedNodes(kk)}{3}(1),tree{highlightedNodes(kk)}{3}(2),-tree{highlightedNodes(kk)}{3}(3),'Marker','o','MarkerSize' , 10, 'LineWidth', 0.1, 'MarkerFaceColor',rndclr{1},'MarkerEdgeColor' , 'k' )
    end
end



for mm = 1: numel(specialNodes)
    for kk=1:size(specialNodes{mm}, 1)
        if pixelUnits
            hSyn = plot3(specialNodes{mm}(kk,1)*relativeRes(1),specialNodes{mm}(kk,2)*relativeRes(2),-specialNodes{mm}(kk,3)*relativeRes(3),'Marker',symCell{rem(mm-1,numel(symCell))+1}, 'MarkerSize' , 3,  'MarkerFaceColor', synapseColor(mm,:) , 'MarkerEdgeColor' , [0.8 0.8 0.8]);
%             drawnow;
%             hMarkers = hSyn.MarkerHandle;
%             hMarkers.FaceColorData = uint8(255*[synapseColor(mm,:),0.5])';  % Alpha=0.5 => 50% transparent 
        else
            hSyn = plot3(specialNodes{mm}(kk,1),specialNodes{mm}(kk,2),-specialNodes{mm}(kk,3),'Marker',symCell{rem(mm-1,numel(symCell))+1}, 'MarkerSize' , 3, 'LineWidth', 0.1 , 'MarkerFaceColor',synapseColor(mm,:) , 'MarkerEdgeColor' , [0.8 0.8 0.8] );
%             drawnow;
%             hMarkers = hSyn.MarkerHandle;
%             hMarkers.FaceColorData =  uint8(255*[synapseColor(mm,:),0.5])';  % Alpha=0.5 => 50% transparent red
        end
    end
end




%axis vis3d;
box on;
axis([ 20000 140000 60000 250000 -60000 0]);
plot( [20000, 40000], [70000, 70000],'-k' ) % insert 20um sclaebar
daspect([1 1 1]); % make aspect ratio [1 1 1]
%set (gca,'Ydir','reverse');
set (gca,'XTick',[], 'YTick',[],'ZTick', [], 'Ydir','reverse');
view([-180,90]); % xy view


