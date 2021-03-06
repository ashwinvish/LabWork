%clc;
clear all;

if ismac
     addpath(genpath('/Users/ashwin/Documents/'));
    df = readtable('/Users/ashwin/Google Drive/Zfish/SynapseDetector/04152019.csv');
else
    addpath(genpath('/usr/people/ashwinv/seungmount/research/Ashwin/Scripts'));
    df = readtable('/usr/people/ashwinv/seungmount/research/Ashwin/SynapseDetector/11252018.csv');
end

% Columns 1 through 8
%     'psd_segid'    'BBOX_bx'    'BBOX_by'    'BBOX_bz'    'BBOX_ex'    'BBOX_ey'    'BBOX_ez'    'postsyn_sz' 
% Columns 9 through 16 
%     'postsyn_wt'    'postsyn_x'    'postsyn_y'    'postsyn_z'    'presyn_sz'    'presyn_wt'    'presyn_x'    'presyn_y' 
% Columns 17 through 23
%     'presyn_z'    'size'    'postsyn_segid'    'presyn_segid'    'centroid_x'    'centroid_y'    'centroid_z'
 


%eg = (df.presyn_seg(df.postsyn_seg==76181));% find presynaptic partners of a cell

%eg = (df.postsyn_seg(df.presyn_seg==76181));% find postsynaptic partners of a cell

% find location of synapses that need to be seeded

% Fully chopped axon definitions
cellIDs_old = {'Int1_1','Int1_2', 'Int1_3','Int1_4', 'Int1_5' ,'Int1_6','Int1_7' ,'Int2_1' , 'Int2_2','Int2_3' , ...
   'Int2_4','Int2_5','Int2_6', 'Int2_7', 'Int2_8',  'Int2_9', 'Int3_1','Int3_2', 'Int3_3' 'Int3_4', 'Int3_5',  'Int3_6' };    % all CellIDS
functionalCellIDs_old = [76198, 76199, 76200, 76201, 76181, 76187, 76210, 76182, 76183, 76185, 76186, 76189, 76184, 76190, ...
    76191, 76192, 76193, 76194, 76195, 76196,76188, 76197];

cellIDsAlx = {'Int1_4','Int1_5','Int1_6','Int2_6','Int2_9','Int3_6'};                                               % all the alx cells, Int2_8 , Int 3_5, Int 1_7
cellIDsDbx = {'Int1_2','Int1_3','Int2_1','Int2_2','Int2_3','Int2_4','Int2_5', 'Int2_8', 'Int3_5'};                  % all bdx1b cells
cellIDsTrans = {'Int1_7'};                                                                                          % all cells with ipsi and contra projections            
cellIDsL = {'Int1_1', 'Int2_7', 'Int3_1', 'Int3_2', 'Int3_3', 'Int3_4'};                                            % all barhl1 cells
cellIDsAxon = {'Int1_4','Int1_5','Int1_6','Int1_7','Int2_6','Int2_9','Int3_5','Int3_6'};

% load Trees from manual tracing
if ismac
    fname = '/Users/ashwin/Documents/LabWork/CurrBio/AVTraces-Exported-01122016-Chopped14';
else
    fname = '/usr/people/ashwinv/seungmount/research/Ashwin/Scripts/CurrBio/AVTraces-Exported-01122016-Chopped14';
end
    
resolution = [5,5,45];

for kk = 1: numel(cellIDs_old)
    disp(fullfile(fname,[cellIDs_old{kk} , '_WithTags.swc']));
    [thisTree,rawLength,thisPreSynapse] = generateIrreducibleDoubleLinkedTree_WithDim(fullfile(fname,[cellIDs_old{kk} , '_WithTags.swc']),[-1:10],6, true, resolution);     % 6 = presynaptic sites
    [thisTree,rawLength,thisPostSynapse] = generateIrreducibleDoubleLinkedTree_WithDim(fullfile(fname,[cellIDs_old{kk} , '_WithTags.swc']),[-1:10],5, true, resolution);    % 5 = postsynapses
    [thisTree,rawLength,thisSpine] = generateIrreducibleDoubleLinkedTree_WithDim(fullfile(fname,[cellIDs_old{kk} , '_WithTags.swc']),[-1:10],9, true, resolution);          % 9 = spines
    [thisTree,rawLength,thisEnd] = generateIrreducibleDoubleLinkedTree_WithDim(fullfile(fname,[cellIDs_old{kk} , '_WithTags.swc']),[-1:10],4, true, resolution);            % 4 = exits thevolume
    allTrees{kk} = thisTree;
    allPreSynapse{kk} = thisPreSynapse; 
    allPostSynapse{kk} = thisPostSynapse;
    allSpine{kk} = thisSpine;
    allRawLength{kk} = rawLength;
    allEnds{kk} = thisEnd;
    allPost{kk} = vertcat(thisPostSynapse, thisSpine); % location of all postsynapses
end

% load synapses from AutoTracing

cellIDs_new = {'Int1_5' , 'Int2_1', 'Int2_2', 'Int2_6', 'Int2_3', 'Int2_4', 'Int1_6', 'Int3_5','Int2_5' ...
    ,'Int2_7', 'Int2_8', 'Int2_9', 'Int3_1', 'Int3_2', 'Int3_3', 'Int3_4', 'Int3_6', 'Int1_1', 'Int1_2', 'Int1_3', 'Int1_4', 'Int1_7'};

functionalCellIDs_new = [76181,76182,76183,76184,76185,76186,76187,76188,76189,76190,76191 ...
    ,76192,76193,76194,76195,76196,76197,76198,76199,76200,76201,76210];

% old and new cell sort order cellIds_new = cellIds_old(sortOrder)
sortOrder = [];
for i = 1:1:size(functionalCellIDs_new,2)
    sortOrder = [sortOrder, find(functionalCellIDs_old==functionalCellIDs_new(i))];
end

% Alx,Dbx,Barhl,T order

[a,b] = intersect(cellIDs_new,cellIDsAlx);
[c,d] = intersect(cellIDs_new,cellIDsDbx);
[e,f] = intersect(cellIDs_new,cellIDsL);
[g,h] = intersect(cellIDs_new,cellIDsTrans);

CellDisplayOrder = [b;d;f;h];
%map = cbrewer('seq','PuRd',255);
map = colorcet('L8');
%% compare synapses
colors = cbrewer('qual','Dark2',10);

for i = 1:size(functionalCellIDs_new,2)
    PrePartners{i} = SynapticPartners(functionalCellIDs_new(i),1,df);
    PostPartners{i} = SynapticPartners(functionalCellIDs_new(i),2,df);
end

synapses_old = cellfun(@length, allPost);
synapses_new = cellfun(@length, PrePartners);

%g = gramm('x', {1:22,1:22},'y',{synapses_old, synapses_new});

plot(1:22, synapses_old(sortOrder), 'o','MarkerFaceColor',colors(1,:));
hold on;
plot(1:22, synapses_new, 'o', 'MarkerFaceColor', colors(2,:));
% x = [0, max(synapses_old)+10];
% y = [0,max(synapses_old)+10];
% line(x,y,'color','black','LineStyle','--');

xlabel('CellNumber');
ylabel('Number of Synapses');
box off;
axis square;
set(gca,'XLim',[0,23],'FontName','Arial','FontSize',25);
legend({'old','new'})


%% Integrator partners only
temp = [];
intPartners = [PrePartners, PostPartners];
for i = 1: size(intPartners,2)
    temp = [temp ; intPartners{i}];
end
intPartners = unique(temp(temp<1e5)); 
IntConnMatrixPre = zeros(size(intPartners));
for i = 1:size(intPartners,1) % pre
        tempPrePartner =  df.presyn_segid(df.postsyn_segid==intPartners(i));
        tempPrePartner = tempPrePartner(tempPrePartner<1e5);
        if ~tempPrePartner == 0
            %i
            [N,edges] = histc(tempPrePartner, unique(tempPrePartner));
            [a,b] = intersect(intPartners,unique(tempPrePartner));
            for j = 1:size(b,1)
             IntConnMatrixPre(i,b(j)) = N(find(a(j)==unique(tempPrePartner)));
            end
        else 
            continue;
        end
        IntConnMatrixPre(i,i) = 0;
        clear tempPrePartner;      
end
save('IntPartners.mat','intPartners');
save('IntConnMatrixPre.mat','IntConnMatrixPre');

%% Construct connctivity matrix

AllCells = [];

AllCells = unique([unique(df.postsyn_segid(df.postsyn_segid<1e5)); unique(df.presyn_segid(df.presyn_segid<1e5))]);
AllCells = AllCells(AllCells<1e5);
ConnMatrixPre = zeros(length(AllCells));
clear tempPrePartner

[a,b,c] = intersect(functionalCellIDs_new(CellDisplayOrder),AllCells(1:22),'stable');
AllCellsSwap(1:size(c,1),:) = AllCells(c,:);
AllCellsSwap(size(c,1)+1:size(AllCells,1),:) = AllCells(size(c,1)+1:end,:);

AllCells = AllCellsSwap;

for i = 1:size(AllCells,1) % pre
        tempPrePartner =  df.presyn_segid(df.postsyn_segid==AllCells(i));
        tempPrePartner = tempPrePartner(tempPrePartner<1e5);
        if ~tempPrePartner == 0
            %i
            [N,edges] = histc(tempPrePartner, unique(tempPrePartner));
            [a,b] = intersect(AllCells,unique(tempPrePartner));
            for j = 1:size(b,1)
             ConnMatrixPre(i,b(j)) = N(j);
            end
        else 
            continue;
        end
        ConnMatrixPre(i,i) = 0;
        clear tempPrePartner;      
end



ConnMatrixPost = zeros(length(AllCells));
clear tempPostPartner
for i = 1:size(AllCells,1) % pre
        tempPostPartner =  df.postsyn_segid(df.presyn_segid==AllCells(i));
        tempPostPartner = tempPostPartner(tempPostPartner<1e5);
        if ~tempPostPartner == 0
            %i
            [N,edges] = histc(tempPostPartner, unique(tempPostPartner));
            [a,b] = intersect(AllCells,unique(tempPostPartner));
            for j = 1:size(b,1)
             ConnMatrixPost(i,b(j)) = N(j);
            end
        else 
            continue;
        end
        ConnMatrixPost(i,i) = 0;
        clear tempPostPartner;      
end 

save('AllCells.mat','AllCells');

%% Pre and Postsynaptic connectivity matrix


figure;
cspy(ConnMatrixPre,'Colormap',map,'Levels',255,'MarkerSize',15);
%imagesc(ConnMatrixPre); axis square;df.
xlabel('Presynaptic cell');
ylabel('Postsynaptic cell');
%map1 = cbrewer('seq','BuPu',255);
%colormap hot;
%colorcet('CBTL1');
set(gcf,'Color','white');
set(gca, 'FontName','Arial','FontSize',25);

figure;
cspy(ConnMatrixPost,'Colormap', map,'Levels',255,'MarkerSize',15)
%imagesc(ConnMatrixPost); axis square;
xlabel('Postsynaptic cell');
ylabel('Presynaptic cell');
%colormap hot;
%colorcet('CBTL1');
set(gcf,'Color','white');
set(gca, 'FontName','Aria','FontSize',25);

save('ConnMatrixPre.mat','ConnMatrixPre');


%% Conectivity for Intcells Only
[a,b] = ismember(functionalCellIDs_new(CellDisplayOrder),AllCells);

IntPreSynapseConn = ConnMatrixPre(b',:); 
IntPreSynSum = sum(IntPreSynapseConn,1);

figure;
subplot_tight(2,1,1,0.05);
bar(1:size(AllCells,1),IntPreSynSum,'FaceColor', colors(1,:));
set(gca, 'XLim', [1,size(AllCells,1)],'XTick',[],'XColor','none');
set(gca, 'FontName','Aria','FontSize',25);
ylabel('count');
box off;

subplot_tight(2,1,2,0.05);
cspy(IntPreSynapseConn,'Colormap',parula(255),'Levels',255,'MarkerSize',40);
set(gca,'XTick',1:size(AllCells,1),'XTickLabel',[],'XTickLabelRotation',45,'XAxisLocation','top');
box on;
%imagesc(IntPreSynapseConn);
%colorcet('CBTL2');
%colormap hot; colorbar
hold on;
line([0,size(AllCells,1)],[6.5,6.5],   'Color','w','LineWidth',4); % Alx block
line([0,size(AllCells,1)],[15.5,15.5], 'Color','w','LineWidth',4); % Dbx block
line([0,size(AllCells,1)],[21.5,21.5], 'Color','w','LineWidth',4); % Barhl block

xlabel('Presynaptic cell');
ylabel('Integrator cell');
%colorcet('CBTL1');
set(gcf,'Color','white');
set(gca, 'FontName','Aria','FontSize',25);


%% Combine Pre and Post into single matrix

IntConn = zeros(2*length(cellIDs_new),2*size(AllCells,1));
IntConn(1:length(cellIDs_new),1:size(AllCells,1)) = IntPreSynapseConn;
IntConn(length(cellIDs_new)+1:end,size(AllCells,1)+1:end) = IntPostSynapseConn;

cspy(IntConn,'Colormap',map,'Levels',255,'MarkerSize',25);
box on
%imagesc(IntConn);
%colorcet('CBTL2');
%colormap hot; 
%colorbar
hold on;
line([0,size(AllCells,1)],[6.5,6.5],   'Color','black'); % Alx block
line([0,size(AllCells,1)],[15.5,15.5], 'Color','black'); % Dbx block
line([0,size(AllCells,1)],[21.5,21.5], 'Color','black'); % Barhl block

line([0,2*size(AllCells,1)],[22.5, 22.5],  'Color','black'); % end of block
line([size(AllCells,1), size(AllCells,1) ],[0, 44.5],  'Color','black');

line([size(AllCells,1),2*size(AllCells,1)],[6.5+22,6.5+22],   'Color','black'); % Alx block
line([size(AllCells,1),2*size(AllCells,1)],[15.5+22,15.5+22], 'Color','black'); % Dbx block
line([size(AllCells,1),2*size(AllCells,1)],[21.5+22,21.5+22], 'Color','black'); % Barhl block


xlabel('Presynaptic cell / Postsynaptic cell');
ylabel('Integrator cell');
%colorcet('CBTL1');
set(gcf,'Color','white');
set(gca, 'FontName','Aria','FontSize',25);


%% post synaptic cells to IntCells
IntPostSynapseConn = zeros(length(cellIDs_new), size(AllCells,1));
[m,n] = ismember(functionalCellIDs_new(CellDisplayOrder),AllCells);

IntPostSynapseConn = ConnMatrixPost(n',:);
IntPostSynSum = sum(IntPostSynapseConn,1);

figure;
subplot_tight(2,1,1,0.05)
bar(1:size(AllCells,1),IntPostSynSum,'FaceColor',colors(1,:))
set(gca, 'XLim', [1,size(AllCells,1)],'YLim',[0,20],'XTick',[],'XColor','none');
set(gca, 'FontName','Aria','FontSize',25);
ylabel('count');
box off;

subplot_tight(2,1,2,0.05)
cspy(IntPostSynapseConn,'Colormap',map,'Levels',255,'MarkerSize',25);
box on;
%imagesc(IntPostSynapseConn);
%colormap hot; colorbar
%colorcet('CBTL1');
hold on;
line([0,size(AllCells,1)],[6.5,6.5],   'Color','black'); % Alx block
line([0,size(AllCells,1)],[15.5,15.5], 'Color','black'); % Dbx block
line([0,size(AllCells,1)],[21.5,21.5], 'Color','black'); % Barhl block

xlabel('Postsynaptic cell');
ylabel('Integrator cell');
%colorcet('CBTL2');
set(gcf,'Color','white');
set(gca, 'FontName','Aria','FontSize',25);

%% Distributions

histogram(IntPreSynapseConn(IntPreSynapseConn>0), 'LineWidth', 2 ...
    , 'FaceColor',colors(1,:));
hold on;
histogram(IntPostSynapseConn(IntPostSynapseConn>0), 'LineWidth', 2 ...
    ,'FaceColor',colors(2,:));
axis square
box off
xlabel('Number of synapses');
ylabel('count');
legend({'PreSynapses','PostSynapses'});
set(gca,'FontName','Arial','FontSize',25);

%%

WPre = zeros(length(cellIDs_new));

for i = 1:1:length(cellIDs_new)
    for j = 1:1:length(cellIDs_new)
        WPre(i,j) = IntPreSynapseConn(i,n(j));
    end
end
figure;
imagesc(WPre);
colormap hot; colorbar
colorcet('CBTL1');
hold on;
line([0,size(cellIDs_new,2)],[6.5,6.5],   'Color','white'); % Alx block
line([0,size(cellIDs_new,2)],[15.5,15.5], 'Color','white'); % Dbx block
line([0,size(cellIDs_new,2)],[21.5,21.5], 'Color','white'); % Barhl block

line([6.5,6.5],[0,size(cellIDs_new,2)],   'Color','white'); % Alx block
line([15.5,15.5],[0,size(cellIDs_new,2)], 'Color','white'); % Dbx block
line([21.5,21.5],[0,size(cellIDs_new,2)],  'Color','white'); % Barhl block


xlabel('Pre cell');
ylabel('Post cell');
%colorcet('CBTL2');
set(gcf,'Color','white');
set(gca, 'FontName','Aria','FontSize',25);
axis square;

WPost = zeros(length(cellIDs_new));

for i = 1:1:length(cellIDs_new)
    for j = 1:1:length(cellIDs_new)
        WPost(i,j) = IntPostSynapseConn(i,n(j));
    end
end

figure;
imagesc(WPost);
colormap hot; colorbar
colorcet('CBTL1');
hold on;

line([0,size(cellIDs_new,2)],[6.5,6.5],   'Color','white'); % Alx block
line([0,size(cellIDs_new,2)],[15.5,15.5], 'Color','white'); % Dbx block
line([0,size(cellIDs_new,2)],[21.5,21.5], 'Color','white'); % Barhl block

line([6.5,6.5],[0,size(cellIDs_new,2)],   'Color','white'); % Alx block
line([15.5,15.5],[0,size(cellIDs_new,2)], 'Color','white'); % Dbx block
line([21.5,21.5],[0,size(cellIDs_new,2)],  'Color','white'); % Barhl block

xlabel('Post cell');
ylabel('Pre cell');
%colorcet('CBTL2');
set(gcf,'Color','white');
set(gca, 'FontName','Aria','FontSize',25);
axis square;

%%



