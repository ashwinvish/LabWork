function [ h4 ] = PlotViews( figureHandle )
%PlotViews Plots the three views xy,xz,yz in a single figure window
%   figureHandle is the handle for the figure, typically the input is the
%   XY view

% handle for the XY view
h1 = figureHandle;

%YZ view
h2 = figure(2);
copyobj(get(h1,'children'),h2);
%view([0,0]);camroll(90);
view([-90,0]); camroll(90);

% XZ view
h3 = figure(3);
copyobj(get(h1,'children'),h3);
view([-180, 0]);

h4 = figure(4);

newh1 = copyobj(get(h1,'children'),h4);
set(newh1,'Units','Normalized','Position',[0.3 0.3 0.3250 0.8150]);
pause(2);
newh2 = copyobj(get(h2,'children'),h4);
set(newh2,'Units','Normalized','Position',[0.5 0.3 0.3250 0.8150]);
newh3 = copyobj(get(h3,'children'),h4);
set(newh3,'Units','Normalized','Position',[0.3 -0.05 0.3250 0.8150]);
pbaspect([1 1 1]);


%close(h1);
%close(h2);
%close(h3);
end

