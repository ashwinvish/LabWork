function treeVisualizer(tree,highlightedNodes,inducingNodes,specialNodes,newFigure,colorString,validNodes,pixelUnits)
relativeRes = [5 5 45]; % in nm
rndclr = colorString;
symCell={'o','s','v','x','d','*'};
synapseColor = {'r', 'g'};
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
    figure;hold;
   % h = gcf;
   % OuterPos = [83 , 37, 1023, 1538];
   % set(h, 'Position',OuterPos);
else
    hold on;
    %h = gcf;
    %OuterPos = [83 , 37, 1023, 1538];
    %set(h, 'Position',OuterPos);
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
            plot3(tempy,tempx,tempz,'color', [1,0.3,0]);
        else
            if ismember(kk,validNodes)
                %plot3(tempy,tempx,tempz,colorString);
                %plot3(tempy,tempx,tempz,'Visible','off');
                h = plot3(tempy,tempx,tempz,'color', rndclr{1});
                h.Color(4) = 0.5;
            end
        end
    end
end

for kk=1:numel(highlightedNodes)
    if pixelUnits
        plot3(tree{highlightedNodes(kk)}{3}(2)*relativeRes(2),tree{highlightedNodes(kk)}{3}(1)*relativeRes(1),-tree{highlightedNodes(kk)}{3}(3)* relativeRes(3),'Marker','o','MarkerSize',10 , 'MarkerFaceColor',rndclr{1},'MarkerEdgeColor','k' )
    else
        plot3(tree{highlightedNodes(kk)}{3}(2),tree{highlightedNodes(kk)}{3}(1),-tree{highlightedNodes(kk)}{3}(3),'Marker','o','MarkerSize' , 10,  'MarkerFaceColor',rndclr{1},'MarkerEdgeColor' , 'k' )
    end
end



for mm = 1: numel(specialNodes)
    for kk=1:size(specialNodes{mm}, 1)
        if pixelUnits
            plot3(specialNodes{mm}(kk,2)*relativeRes(2),specialNodes{mm}(kk,1)*relativeRes(1),-specialNodes{mm}(kk,3)*relativeRes(3),'Marker',symCell{rem(mm-1,numel(symCell))+1}, 'MarkerSize' , 4,  'MarkerFaceColor', synapseColor{mm} , 'MarkerEdgeColor' , 'none');
        else
            plot3(specialNodes{mm}(kk,2),specialNodes{mm}(kk,1),-specialNodes{mm}(kk,3),'Marker',symCell{rem(mm-1,numel(symCell))+1}, 'MarkerSize' , 4,  'MarkerFaceColor', synapseColor{mm} , 'MarkerEdgeColor' , 'none');
        end
    end
end

%daspect([1 1 1]);
axis([60000 250000 20000 140000 -60000 0]);
plot([240000, 240000], [120000, 140000], '-k' ) % insert 20um sclaebar
%view([90,90]); % xy view


