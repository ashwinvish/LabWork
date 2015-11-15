%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%;
clc;
clear all;
resolution= [5,5,45];

% Axonal nodes for all the cells
Int1_1_axon = [];
Int1_2_axon = [15 5 6 40];
Int1_3_axon = [1,73];
Int1_4_axon = [30 42 43 44 45 47 48 41 52 49 50 51 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74];
Int1_5_axon = [56 106 107 109 110 114 115 116 159 160 161 162 163 164 165 166 167 169 170 171 168 77 80 81 82 84 86 87 85 92 99 100 103 115 117 123 125 126 137 139 150 151 153 127 124 128 131 132 134 140 143 152 154 155 156 157 158 136 138 141 142 144 106 110 114 116 133 135 145 146 147 148 149 ];
Int1_6_axon = [20 82  96 97 98 133 138 143 147 149 150 156 154 148 151 152 123  126 132 128 129 130 134 140 142 153 165 139 141 155 157 159 160 163 162 158 164 166 131 120 135 145 161 167 168 169 170 171 172 129 136 144 146 137 74 118 99];
Int1_7_axon = [9 30 75 81 99 115 120 128 110 114 116 117 121 122 123 129 76 107];
Int2_1_axon = [2 5 7 14 15 25 26 73];
Int2_2_axon = [5 84];
Int2_3_axon = [3 67];
Int2_4_axon = [2 101];
Int2_5_axon = [6 64 65 66 67 69];
Int2_6_axon = [];
Int2_7_axon = [];
Int2_8_axon = []; % 10 12 22 64
Int2_9_axon = [7 17 43 47 55 66 67 78 79 80 82 83 84 85 86 87 88 81];
Int3_1_axon = [];
Int3_2_axon = [];
Int3_3_axon = [];
Int3_4_axon = [];
Int3_5_axon = [34 52 80 82 83 86 88 90 103 104 107 111 112 113 114 115 92 81 84 87 89 95 96 98 99 100 93 94 101 102 105 106 108 109 110];
Int3_6_axon = [5 29 37 39 42 52 53 5];
MauthnerCell = [5*14474,5*49530,-45*448];                                                       % cartesian coordinates for the center of the Mauthner cell

cellIDs = {'Int1_1','Int1_2', 'Int1_3','Int1_4', 'Int1_5' ,'Int1_6','Int1_7' ,'Int2_1' , 'Int2_2','Int2_3' , 'Int2_4','Int2_5','Int2_6', 'Int2_7', 'Int2_8',  'Int2_9', 'Int3_1','Int3_2', 'Int3_3' 'Int3_4', 'Int3_5',  'Int3_6' };								      % all CellIDS
cellIDsAlx = {'Int1_4','Int1_5','Int1_6','Int1_7','Int2_6','Int2_8','Int2_9','Int3_5','Int3_6'};			% all the alx cells
cellIDsDbx = {'Int1_2','Int1_3','Int2_1','Int2_2','Int2_3','Int2_4','Int2_5'};			% all bdx1b cells
cellIDsL = {'Int1_1', 'Int2_7', 'Int3_1', 'Int3_2', 'Int3_3', 'Int3_4'};                        % all barhl1 cells

% Convert from .swc file to tree structre with presynapses,postsynapses,
for kk = 1: numel(cellIDs)
    disp([cellIDs{kk} , '_WithTags.swc']);
    [thisTree,rawLength,thisPreSynapse] = generateIrreducibleDoubleLinkedTree_WithDim([cellIDs{kk} , '_WithTags.swc'],[-1:10],6, true, resolution);
    [thisTree,rawLength,thisPostSynapse] = generateIrreducibleDoubleLinkedTree_WithDim([cellIDs{kk} , '_WithTags.swc'],[-1:10],5, true, resolution);
    [thisTree,rawLength,thisSpine] = generateIrreducibleDoubleLinkedTree_WithDim([cellIDs{kk} , '_WithTags.swc'],[-1:10],9, true, resolution);
    allTrees{kk} = thisTree; allPreSynapse{kk} = thisPreSynapse; allPostSynapse{kk} = thisPostSynapse;allSpine{kk} = thisSpine;
    allRawLength{kk} = rawLength; allPost{kk} = vertcat(thisPostSynapse, thisSpine);
end

% control cells
cellControl = {'C1','C2','C3','C4','C5','C6','C7'};
TreeC1 = generateIrreducibleDoubleLinkedTree_WithDim('C1 [treeline] #317070.swc',[-1:10],5,true, resolution);
TreeC2 = generateIrreducibleDoubleLinkedTree_WithDim('C2 [treeline] #317065.swc',[-1:10],5,true, resolution);
TreeC3 = generateIrreducibleDoubleLinkedTree_WithDim('C3 [treeline] #317067.swc',[-1:10],5,true, resolution);
TreeC4 = generateIrreducibleDoubleLinkedTree_WithDim('C4 [treeline] #317072.swc',[-1:10],5,true, resolution);
TreeC5 = generateIrreducibleDoubleLinkedTree_WithDim('C5 [treeline] #317049.swc',[-1:10],5,true, resolution);
TreeC6 = generateIrreducibleDoubleLinkedTree_WithDim('C6 [treeline] #317054.swc',[-1:10],5,true, resolution);
TreeC7 = generateIrreducibleDoubleLinkedTree_WithDim('C7 [treeline] #316984.swc',[-1:10],5,true, resolution);

% load anatomical features

% Commisure1 = load('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/SWCFiles/8kTraces/Commisure1.mat');
% Commisure1 = Commisure1.Commisure1;
% %AnatomicalSpline(5*Commisure1(:,1),5*Commisure1(:,2),-45*Commisure1(:,3), [0.5,0.5,0.5]);
% 
% Commisure2 = load('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/SWCFiles/8kTraces/Commisure2.mat');
% Commisure2 = Commisure2.Commisure2;
% %AnatomicalSpline(5*Commisure2(:,1),5*Commisure2(:,2),-45*Commisure2(:,3), [0.5,0.5,0.5]);
% 
% Commisure3 = load('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/SWCFiles/8kTraces/Commisure3.mat');
% Commisure3 = Commisure3.Commisure3;
% %AnatomicalSpline(5*Commisure3(:,1),5*Commisure3(:,2),-45*Commisure3(:,3), [0.5,0.5,0.5]);
% 
% Commisure4 = load('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/SWCFiles/8kTraces/Commisure4.mat');
% Commisure4 = Commisure4.Commisure4;
% %AnatomicalSpline(5*Commisure4(:,1),5*Commisure4(:,2),-45*Commisure4(:,3), [0.5,0.5,0.5]);
% 
% Commisure5 = load('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/SWCFiles/8kTraces/Commisure5.mat');
% Commisure5 = Commisure5.Commisure5;
% %AnatomicalSpline(5*Commisure5(:,1),5*Commisure5(:,2),-45*Commisure5(:,3), [0.5,0.5,0.5]);
% 
% Commisure6 = load('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/SWCFiles/8kTraces/Commisure6.mat');
% Commisure6 = Commisure6.Commisure6;
% %AnatomicalSpline(5*Commisure6(:,1),5*Commisure6(:,2),-45*Commisure6(:,3), [0.5,0.5,0.5]);
% 
% Commisure7 = load('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/SWCFiles/8kTraces/Commisure7.mat');
% Commisure7 = Commisure7.Commisure7;
% %AnatomicalSpline(5*Commisure7(:,1),5*Commisure7(:,2),-45*Commisure7(:,3), [0.5,0.5,0.5]);
% 
% MauthnerAxon = load('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/SWCFiles/8kTraces/MauthnerAxon.mat');
% MauthnerAxon = MauthnerAxon.MauthnerAxon;
% %AnatomicalSpline(5*MauthnerAxon(:,1),5*MauthnerAxon(:,2),-45*MauthnerAxon(:,3), [0.5,0.5,0.5]);


% box on;
% axis([ 20000 140000 60000 250000 -60000 0]);
% plot( [20000, 40000], [70000, 70000],'-k' ) % insert 20um sclaebar
% daspect([1 1 1]); % make aspect ratio [1 1 1]
% %set (gca,'Ydir','reverse');
% set (gca,'XTick',[], 'YTick',[],'ZTick', [], 'Ydir','reverse');
% view([-180,90]); % xy view

% Load Time constants

%time constants and persistence measure
tau = [7.16323454600000,6.23690945500000,9.88816223500000,17.7411357900000,5.48599643300000,10.1172859600000,7.51255049400000,18.1882559300000,57.0990589000000,100,38.7012773800000,11.6985985200000 ...
    11.0576103400000,1.44892223200000,19.5401836000000,100,16.5597699400000,10.3971051900000,26.7744003700000,7.92791676600000,6.10909614600000,5.21884052300000];
rho = [0.5207986709, 0.6661018539, 0.6648735491, 0.8190294252, 0.6370840437, 0.3520511176, 0.6180659652, 0.7770252108, 0.8445793478, 0.5048882377, 0.6169858888, 0.4693437009, 0.7319155263 ...
    0.2311146742, 0.5600555905, 0.69884915, 0.8217790613, 0.5105536878, 0.8749909893, 0.3470572669, 0.5251032841, 0.5768530374];

% Control cells details

ControlCellSoma(1,:) =  TreeC1{1,1}{1,3};
ControlCellSoma(2,:) =  TreeC2{1,1}{1,3};
ControlCellSoma(3,:) =  TreeC3{1,1}{1,3};
ControlCellSoma(4,:) =  TreeC4{1,1}{1,3};

stripe1 = [ControlCellSoma(1,:);ControlCellSoma(2,:);ControlCellSoma(3,:);ControlCellSoma(4,:)];	%stripe1, corresponds with alx transcription factor

ControlCellSoma(5,:) =  TreeC5{1,1}{1,3};
ControlCellSoma(6,:) =  TreeC6{1,1}{1,3};
ControlCellSoma(7,:) =  TreeC7{1,1}{1,3};

stripe2 = [ControlCellSoma(7,:);ControlCellSoma(5,:);ControlCellSoma(6,:)];                     	%stripe2, corresponds with dbx transcription factor



