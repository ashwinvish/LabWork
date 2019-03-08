% make_final_struct8.m
% ----------------------------
% modified from: dsred_posneg_trends_new.m
% run after: image_analysis_08_12.m
%           analyze_spont2.m
%           analyze_sine.m
%           plot_data.m

%           cell_classify3.m
%           cell_loc2.m/mod_cell_loc2.m
%           draw_stripe_boundaries.m
%           assign_stripes.m

% batch processing
% - modified from make_stripe_bar_plot.m global_trends_Delta_new.m,
% compare-a1a2_taus.m
% - collect (x,y,z) positions and tau_s, tau_k, rsq values from cells in
% s1,s2 and plot gradients

% 2/28/13: modified from stripe_trends_new.m
% separate cells into classes based on whether they are DsRed positive or
% negative

% 3/4/13: use data from new constrained fits (OKS4, SPONT4)

% 3/22/13: create a new structure, indexed by cell (include all VPNI, VSNI,
% OS cells)

% 5/24/13: exclude IO cells from VPNI, VSNI, OS cell classes; add in
% updated XFP identification (using xfp_sep3.m, xfp_ip3.m)

% 6/11/13: save correlation coefficients used in functional classification

% 7/1/13: save OKS6 data; includes fits of direct pathway + measurements of
% k, Eth

% 7/8/13: save RGratio for Dsred cells

% 7/24/13: track glyc cells

% 8/21/13: track DsRed neg cells, unknown cells

% 8/23/13: work with SPONT8, new fits on csaps filtered STAs

% 8/26/13: added rho

% 8/31/13: rho from method 1 fits, SPONT9, run after analyze_spont10.m

% 7/1/14: save rh6/7 border

% 7/10/14: assess stripe distribution using user-defined XFP truth set

% 7/11/14: new ML truth set, run after eval_ML_truth

% 7/22/14: added XFP id params

% 10/28/14: save age of fish (dpf)

% 11/19/14: save deconvolved fluor and fits

% 2/3/15: for unilateral pair gradient measurements, record left
% (L)/right(R) position of cells 

close all; clear all;

% mainfolder = 'C:\Users\mel2011\Documents\Imaging\dual_xfp_expts';
% folders = [strcat(mainfolder,'\101812');strcat(mainfolder,'\101712');strcat(mainfolder,'\101612')];

mainfolder = 'C:\Users\mel2011\Dropbox\Lee_ML\data_mat_files_for xfp';
folders = [strcat(mainfolder,'\101812');strcat(mainfolder,'\101712');strcat(mainfolder,'\101612');strcat(mainfolder,'\081612');strcat(mainfolder,'\080612')];


fishid = [];
vpni = [];
vsni = [];
os = [];
XYZ = [];
stripe = [];
vglut = [];
glyc = []; % added 7/24/13
taus = [];
taus1 = [];
tauk = [];
rsqs = [];
rsqs1 = [];
rsqk = [];
staf = []; % ipsi STA fluor
stad = []; % contra STA fluor
corrspont = [];
corrokr = [];
ksens = []; % sensitivity of r(E)
Eth = [];   % threshold eye position
r0 = [];    % FR threshold
rgratio = []; % only in dual fish
rho = [];
rho_fit = [];
rh67 = [];
XYZadj = [];   % X adjusted to be distance from rh6/7 border (=0 um)
roi_frac = [];
roi_ann = [];
roi_asym = [];
roi_mean_int = [];
age = [];
dstaf = []; % deconvolved fluorescence
ndstaf = []; % normalized deconvolved fluorescence
t_fit = []; % for fit to normalized deconv fluor
fit_res = []; % fit to normalized deconv fluor
LR = []; % left-right position of cells 

for ff=1:3%size(folders,1)
    cd(folders(ff,:));
    files = dir('*_fluor_data.mat');
    
    for fl=1:length(files)
        fl
        %if~(ff==1&&fl==5)    % skip this set (single xfp)
        % if~(ff==1&&(ismember(fl,[1,3,4,5,6]))) % skip these sets (dualxfp)
        %if(ismember(fl,[1,4,10]))
        file = files(fl).name;
        load(file);
        
        % added 3/28/13
        clear SPONT4; clear OKS4;
        %SPONT4 = SPONT6; OKS4 = OKS6;
        SPONT4 = SPONT9; OKS4 = OKS6; % 8/23/13
        
        % correlation for OKS data
        tau_cirf = 1.89;
        ker = exp(-(OKS2.tt-OKS2.tt(1))./tau_cirf);
        coks = zeros(2,length(OKS2.stim_data));
        temp = conv(ker,OKS2.stim_data);
        coks(1,:) = temp(1:length(OKS2.stim_data));
        temp = conv(ker,[0,diff(OKS2.stim_data)]);
        coks(2,:) = temp(1:length(OKS2.stim_data));
        temp = conv(ker,OKS2.bsig);
        coks(3,:) = temp(1:length(OKS2.bsig));
        
        % classify OKS cells by correlation to CIRF-convolved stim, eye variables
        corroks = zeros(3,FLUOR(1).num_rois);
        for i=1:FLUOR(1).num_rois
            for j=1:3
                corroks(j,i) = abs(corr(coks(j,:)',OKS2.finterp(:,i)));
            end
        end
        
        % classify based on spont data
        corrspont_ipsi = max(SPONT.corrc_orig([1,2,4,5],:),[],1);
        corrspont_contra = max(SPONT.corrc_orig([3,6],:),[],1);
        
        dt = mean(diff(SPONT2.T));
        ddfl = zeros(1,FLUOR(2).num_rois); % difference in fluor before and at sacc left/ipsi
        ddfr = zeros(1,FLUOR(2).num_rois);
        vars = zeros(1,FLUOR(2).num_rois);
        for i=1:FLUOR(2).num_rois
            ddfl(i) = mean(SPONT2.staLun(3/dt:8/dt,i))- mean(SPONT2.staLun(1:1/dt,i));
            ddfr(i) = mean(SPONT2.staRun(3/dt:8/dt,i))- mean(SPONT2.staRun(1:1/dt,i));
            vars(i) = var(SPONT.finterp(:,i));
        end
        
        % 6/11/13 classification
        area1_cell = find(((corrspont_ipsi > 0.4 & FLUOR(2).roi_y_pos > 256)|...
            (corrspont_contra > 0.4 & FLUOR(2).roi_y_pos < 256)) & ...
            max(corrspont_ipsi,corrspont_contra)./max(corroks(1,:),corroks(2,:))<3);
        area2_cell = find(max(corroks(1,:),corroks(2,:)) > 0.4 & ...
            max(corroks(1,:),corroks(2,:))./max(corrspont_ipsi,corrspont_contra) > 2 & vars < 50);
        os_cell = find(corrspont_ipsi > 0.3 | corrspont_contra > 0.3 | max(corroks(1,:),corroks(2,:)) > 0.3);
        
        curr_corrspont = zeros(1,FLUOR(1).num_rois);
        curr_corrspont(FLUOR(2).roi_y_pos > 256) = corrspont_ipsi(FLUOR(2).roi_y_pos > 256);
        curr_corrspont(FLUOR(2).roi_y_pos < 256) = corrspont_contra(FLUOR(2).roi_y_pos < 256);
        curr_corrokr = max(corroks(1,:),corroks(2,:));
        
        % added 5/24/13: exclude IO cells
        area1_cell = setdiff(area1_cell,FLUOR(2).io_cells);
        area2_cell = setdiff(area2_cell,FLUOR(2).io_cells);
        os_cell = setdiff(os_cell,FLUOR(2).io_cells);
        
        % added 4/23/14: exclude rh6 cells
        area1_cell = setdiff(area1_cell,FLUOR(2).rh6);
        area2_cell = setdiff(area2_cell,FLUOR(2).rh6);
        os_cell = setdiff(os_cell,FLUOR(2).rh6);
        
        roi_x_pos = zeros(1,FLUOR(1).num_rois);
        curr_start = 1;
        for i=1:FLUOR(1).num_rois
            curr_end = curr_start+FLUOR(1).len_patches(i)-1;
            xi = FLUOR(1).xi_all(curr_start:curr_end);
            yi = FLUOR(1).yi_all(curr_start:curr_end);
            curr_start = curr_end + 1;
            roi_x_pos(i) = mean(xi);
        end
        
        lenstaf = size(SPONT2.T,2);
        curr_staf = zeros(lenstaf,length(area1_cell));
        curr_stad = zeros(lenstaf,length(area1_cell));
        curr_dstaf = zeros(lenstaf,length(area1_cell));
        curr_ndstaf = zeros(lenstaf,length(area1_cell));
        for i=1:length(area1_cell)
            if(FLUOR(2).roi_y_pos(area1_cell(i)) > 256)
                curr_staf(:,i) = SPONT2.staLun(:,area1_cell(i));
                curr_stad(:,i) = SPONT2.staRun(:,area1_cell(i));
            else
                curr_staf(:,i) = SPONT2.staRun(:,area1_cell(i));
                curr_stad(:,i) = SPONT2.staLun(:,area1_cell(i));
            end
            curr_dstaf(:,i) = SPONT4.stafr(:,area1_cell(i));
            curr_ndstaf(:,i) = SPONT4.nstafr(:,area1_cell(i));
        end
        
        curr_vpni = zeros(1,FLUOR(1).num_rois);
        curr_vsni = zeros(1,FLUOR(1).num_rois);
        curr_os = zeros(1,FLUOR(1).num_rois);
        curr_XYZ = zeros(3,FLUOR(1).num_rois);
        curr_rh67 = FLUOR(1).r67_border*ones(1,FLUOR(1).num_rois); % 7/2/14
        
        curr_vpni(area1_cell) = 1;
        curr_vsni(area2_cell) = 1;
        curr_os(os_cell) = 1;
        curr_XYZ(1,:) = roi_x_pos;
        curr_XYZ(2,:) = FLUOR(1).roi_y_pos;
        curr_XYZ(3,:) = FLUOR(1).dv_loc*ones(1,FLUOR(1).num_rois);
        curr_XYZadj = curr_XYZ;         % 7/2/14
        if(FLUOR(1).r67_border < 0)
            curr_XYZadj(1,:) = curr_XYZ(1,:) + FLUOR(1).r67_border;
        else
            curr_XYZadj(1,:) = curr_XYZ(1,:) - FLUOR(1).r67_border;
        end
        currLR = FLUOR(1).roi_y_pos(1:FLUOR(1).num_rois) > mean(ONION.ym); % 1 = RHS, 0 = LHS
        
        curr_vglut = zeros(1,FLUOR(1).num_rois);
        curr_glyc = zeros(1,FLUOR(1).num_rois); % 7/24/13
        %             curr_vglut(FLUOR(1).pos_alg_vglut) = 1; % 7/4/13 - algorithm
        %             curr_vglut(FLUOR(1).neg_alg_vglut) = -1; % 8/21/13
        curr_glyc(FLUOR(1).keep_roi_gly) = 1; % 7/24/13
        
        %             curr_vglut(FLUOR(1).tp_vglut) = 1; % 7/10/14 - user-defined 'truth'
        %             curr_vglut(FLUOR(1).tn_vglut) = -1;
        
        %             curr_vglut(FLUOR(1).ML_pos) = 1; % 7/11/14 - user-defined 'truth'
        %             curr_vglut(FLUOR(1).ML_neg) = -1;
        
        %             curr_vglut(FLUOR(1).EA_pos2) = 1; % 7/18/14 - EA-defined 'truth'
        %             curr_vglut(FLUOR(1).EA_neg2) = -1;
        
        curr_vglut(FLUOR(1).EAML_pos) = 1; % 7/21/14 - EA and ML agreed upon 'truth'
        curr_vglut(FLUOR(1).EAML_neg) = -1;
        curr_vglut(FLUOR(1).EAML_unkn) = 2;
        curr_vglut(FLUOR(1).EAML_dis) = 3;
        
        % store stripe indices
        curr_stripe = zeros(4,FLUOR(1).num_rois);
        curr_stripe(1,MAP2.s1) = 1;
        curr_stripe(2,MAP2.s2) = 1;
        curr_stripe(3,MAP2.s3) = 1;
        curr_stripe(4,MAP2.s4) = 1;
        
        % store fish ID of fish -- dual xfp
        if(ff==1 && fl==2)
            fishnum = 1; disp('fish 1');
            curr_age = 9*ones(1,FLUOR(1).num_rois);
        elseif(ff==1 && ismember(fl,7:13))
            fishnum = 2; disp('fish 2');
            curr_age = 9*ones(1,FLUOR(1).num_rois);
        elseif(ff==2 && ismember(fl,1:7))
            fishnum = 3; disp('fish 3');
            curr_age = 7*ones(1,FLUOR(1).num_rois);
        elseif(ff==2 && ismember(fl,8:12))
            fishnum = 4; disp('fish 4');
            curr_age = 7*ones(1,FLUOR(1).num_rois);
        else
            fishnum = 5; disp('fish 5');
            curr_age = 6*ones(1,FLUOR(1).num_rois);
        end
        
        curr_fish = fishnum*ones(1,FLUOR(1).num_rois);
        
        curr_taus = zeros(1,FLUOR(1).num_rois);
        curr_taus1 = zeros(1,FLUOR(1).num_rois);
        curr_tauk = zeros(1,FLUOR(1).num_rois);
        curr_rsqs = zeros(1,FLUOR(1).num_rois);
        curr_rsqs1 = zeros(1,FLUOR(1).num_rois);
        curr_rsqk = zeros(1,FLUOR(1).num_rois);
        curr_rho = zeros(1,FLUOR(1).num_rois);
        curr_rho_fit = zeros(1,FLUOR(1).num_rois);
        curr_fit_res = zeros(length(SPONT4.t_fit),FLUOR(1).num_rois);
        curr_t_fit = repmat(SPONT4.t_fit',[1 FLUOR(1).num_rois]);
        %
        % store tau_s, tau_k, rsq values (may need to make this more
        % efficient)
        allrois = 1:FLUOR(1).num_rois;
        % use cnull that results in best fit
        [maxr,indm] = max(SPONT4.rsq1full);
        %[maxr,indm] = max(SPONT2.rsq(:,a1s1_rois));
        prs = [indm' allrois'];
        for k=1:FLUOR(1).num_rois
            curr_taus(k) = SPONT4.taup1(prs(k,1),prs(k,2));
            curr_rho_fit(k) = SPONT4.rho_fit(prs(k,1),prs(k,2));
            curr_rho(k) = SPONT4.rho(prs(k,1),prs(k,2));
            %addtaus(k) = SPONT2.taup(prs(k,1),prs(k,2));
            curr_fit_res(:,k) = squeeze(SPONT4.fit_res1(:,prs(k,1),prs(k,2)));
        end
        curr_rsqs = maxr;
        
        
        curr_tauk = OKS4.beta_all(1,:);
        curr_rsqk = OKS4.rsq;
        
        curr_ksens = OKS7.kfit;
        curr_Eth = OKS7.Eth;
        curr_r0 = OKS7.r0;
        
        curr_rgratio = FLUOR(1).RGratio;
        
        curr_roi_frac = FLUOR(1).ROI_frac;
        curr_roi_ann = FLUOR(1).ROI_xfp_mean_vglut./FLUOR(1).ROI_ann_mean_vglut;
        curr_roi_asym = FLUOR(1).ROI_asym_vglut;
        
        % calculate ROI mean int
        curr_roi_mean_int = FLUOR(1).ROI_xfp_mean_vglut;
        
        %             curr_roi_mean_int = zeros(1,FLUOR(1).num_rois);
        %             curr_start = 1;
        %             for k=1:FLUOR(1).num_rois
        %                 curr_end = curr_start+FLUOR(1).len_patches(k)-1;
        %                 xi = FLUOR(1).xi_all(curr_start:curr_end);
        %                 yi = FLUOR(1).yi_all(curr_start:curr_end);
        %                 curr_start = curr_end + 1;
        %                 BW = poly2mask(xi,yi,size(FLUOR(1).ref_red,1),size(FLUOR(1).ref_red,1));
        %                 curr_roi_mean_int(k) = mean(FLUOR(1).ref_red(BW));
        %             end
        
        % save data
        fishid = [fishid,curr_fish];
        vpni = [vpni,curr_vpni];
        vsni = [vsni,curr_vsni];
        os = [os,curr_os];
        XYZ = [XYZ,curr_XYZ];
        stripe = [stripe,curr_stripe];
        vglut = [vglut,curr_vglut];
        glyc = [glyc,curr_glyc]; % 7/24/13
        taus = [taus,curr_taus];
        %taus1 = [taus1,curr_taus1];
        tauk = [tauk,curr_tauk];
        rsqs = [rsqs,curr_rsqs];
        %rsqs1 = [rsqs1,curr_rsqs1];
        rsqk = [rsqk,curr_rsqk];
        staf = [staf curr_staf];
        stad = [stad curr_stad];
        corrspont = [corrspont curr_corrspont];
        corrokr = [corrokr curr_corrokr];
        ksens = [ksens curr_ksens];
        Eth = [Eth curr_Eth];
        r0 = [r0 curr_r0];
        rgratio = [rgratio curr_rgratio];
        rho = [rho curr_rho];
        rho_fit = [rho_fit curr_rho_fit];
        rh67 = [rh67 curr_rh67];        % 7/2/14
        XYZadj = [XYZadj curr_XYZadj];
        roi_frac = [roi_frac curr_roi_frac];
        roi_ann = [roi_ann curr_roi_ann];
        roi_asym = [roi_asym curr_roi_asym];
        roi_mean_int = [roi_mean_int curr_roi_mean_int];
        age = [age curr_age];
        dstaf = [dstaf curr_dstaf];
        ndstaf = [ndstaf curr_ndstaf];
        fit_res = [fit_res curr_fit_res];
        t_fit = [t_fit curr_t_fit];
        LR = [LR,currLR];
        
        % end  % if~....for excluding sets
    end
end

% mainfolder = 'C:\Users\mel2011\Documents\Imaging\xfp_expts';
% % folders = [strcat(mainfolder,'\032712');strcat(mainfolder,'\050112');strcat(mainfolder,'\032612');...
% %     strcat(mainfolder,'\081612');strcat(mainfolder,'\080612')];
% folders = [strcat(mainfolder,'\081612');strcat(mainfolder,'\080612')];

for ff=4:size(folders,1) %1:size(folders,1)
    cd(folders(ff,:));
    files = dir('*_fluor_data.mat');
    
    for fl=1:length(files)
        fl
        %  if~(ff==1&&fl==5)    % skip this set (single xfp)
        %if~(ff==1&&(ismember(fl,[1,3,4,5,6]))) % skip these sets (dualxfp)
        %if(ismember(fl,[1,4,10]))
        file = files(fl).name;
        load(file);
        
        % added 3/28/13
        clear SPONT4; clear OKS4;
        %SPONT4 = SPONT6; OKS4 = OKS6;
        SPONT4 = SPONT9; OKS4 = OKS6; % 8/23/13
        
        % correlation for OKS data
        tau_cirf = 1.89;
        ker = exp(-(OKS2.tt-OKS2.tt(1))./tau_cirf);
        coks = zeros(2,length(OKS2.stim_data));
        temp = conv(ker,OKS2.stim_data);
        coks(1,:) = temp(1:length(OKS2.stim_data));
        temp = conv(ker,[0,diff(OKS2.stim_data)]);
        coks(2,:) = temp(1:length(OKS2.stim_data));
        temp = conv(ker,OKS2.bsig);
        coks(3,:) = temp(1:length(OKS2.bsig));
        
        % classify OKS cells by correlation to CIRF-convolved stim, eye variables
        corroks = zeros(3,FLUOR(1).num_rois);
        for i=1:FLUOR(1).num_rois
            for j=1:3
                corroks(j,i) = abs(corr(coks(j,:)',OKS2.finterp(:,i)));
            end
        end
        
        % classify based on spont data
        corrspont_ipsi = max(SPONT.corrc_orig([1,2,4,5],:),[],1);
        corrspont_contra = max(SPONT.corrc_orig([3,6],:),[],1);
        
        dt = mean(diff(SPONT2.T));
        ddfl = zeros(1,FLUOR(2).num_rois); % difference in fluor before and at sacc left/ipsi
        ddfr = zeros(1,FLUOR(2).num_rois);
        vars = zeros(1,FLUOR(2).num_rois);
        for i=1:FLUOR(2).num_rois
            ddfl(i) = mean(SPONT2.staLun(3/dt:8/dt,i))- mean(SPONT2.staLun(1:1/dt,i));
            ddfr(i) = mean(SPONT2.staRun(3/dt:8/dt,i))- mean(SPONT2.staRun(1:1/dt,i));
            vars(i) = var(SPONT.finterp(:,i));
        end
        
        % 6/11/13 classification
        area1_cell = find(((corrspont_ipsi > 0.4 & FLUOR(2).roi_y_pos > 256)|...
            (corrspont_contra > 0.4 & FLUOR(2).roi_y_pos < 256)) & ...
            max(corrspont_ipsi,corrspont_contra)./max(corroks(1,:),corroks(2,:))<3);
        area2_cell = find(max(corroks(1,:),corroks(2,:)) > 0.4 & ...
            max(corroks(1,:),corroks(2,:))./max(corrspont_ipsi,corrspont_contra) > 2 & vars < 50);
        os_cell = find(corrspont_ipsi > 0.3 | corrspont_contra > 0.3 | max(corroks(1,:),corroks(2,:)) > 0.3);
        
        curr_corrspont = zeros(1,FLUOR(1).num_rois);
        curr_corrspont(FLUOR(2).roi_y_pos > 256) = corrspont_ipsi(FLUOR(2).roi_y_pos > 256);
        curr_corrspont(FLUOR(2).roi_y_pos < 256) = corrspont_contra(FLUOR(2).roi_y_pos < 256);
        curr_corrokr = max(corroks(1,:),corroks(2,:));
        
        % added 5/24/13: exclude IO cells
        area1_cell = setdiff(area1_cell,FLUOR(2).io_cells);
        area2_cell = setdiff(area2_cell,FLUOR(2).io_cells);
        os_cell = setdiff(os_cell,FLUOR(2).io_cells);
        
        %             % added 4/23/14: exclude rh6 cells
        %             area1_cell = setdiff(area1_cell,FLUOR(2).rh6);
        %             area2_cell = setdiff(area2_cell,FLUOR(2).rh6);
        %             os_cell = setdiff(os_cell,FLUOR(2).rh6);
        
        roi_x_pos = zeros(1,FLUOR(1).num_rois);
        curr_start = 1;
        for i=1:FLUOR(1).num_rois
            curr_end = curr_start+FLUOR(1).len_patches(i)-1;
            xi = FLUOR(1).xi_all(curr_start:curr_end);
            yi = FLUOR(1).yi_all(curr_start:curr_end);
            curr_start = curr_end + 1;
            roi_x_pos(i) = mean(xi);
        end
        
        lenstaf = size(SPONT2.T,2);
        curr_staf = zeros(lenstaf,length(area1_cell));
        curr_stad = zeros(lenstaf,length(area1_cell));
        curr_dstaf = zeros(lenstaf,length(area1_cell));
        curr_ndstaf = zeros(lenstaf,length(area1_cell));
        for i=1:length(area1_cell)
            if(FLUOR(2).roi_y_pos(area1_cell(i)) > 256)
                curr_staf(:,i) = SPONT2.staLun(:,area1_cell(i));
                curr_stad(:,i) = SPONT2.staRun(:,area1_cell(i));
            else
                curr_staf(:,i) = SPONT2.staRun(:,area1_cell(i));
                curr_stad(:,i) = SPONT2.staLun(:,area1_cell(i));
            end
            curr_dstaf(:,i) = SPONT4.stafr(:,area1_cell(i));
            curr_ndstaf(:,i) = SPONT4.nstafr(:,area1_cell(i));
        end
        
        curr_vpni = zeros(1,FLUOR(1).num_rois);
        curr_vsni = zeros(1,FLUOR(1).num_rois);
        curr_os = zeros(1,FLUOR(1).num_rois);
        curr_XYZ = zeros(3,FLUOR(1).num_rois);
        curr_rh67 = FLUOR(1).r67_border*ones(1,FLUOR(1).num_rois); % 7/2/14
        
        curr_vpni(area1_cell) = 1;
        curr_vsni(area2_cell) = 1;
        curr_os(os_cell) = 1;
        curr_XYZ(1,:) = roi_x_pos;
        curr_XYZ(2,:) = FLUOR(1).roi_y_pos(1:FLUOR(1).num_rois);
        curr_XYZ(3,:) = FLUOR(1).dv_loc*ones(1,FLUOR(1).num_rois);
        curr_XYZadj = curr_XYZ;         % 7/2/14
        if(FLUOR(1).r67_border < 0)
            curr_XYZadj(1,:) = curr_XYZ(1,:) + FLUOR(1).r67_border;
        else
            curr_XYZadj(1,:) = curr_XYZ(1,:) - FLUOR(1).r67_border;
        end
        currLR = FLUOR(1).roi_y_pos(1:FLUOR(1).num_rois) > mean(ONION.ym); % 1 = RHS, 0 = LHS
        
        curr_vglut = zeros(1,FLUOR(1).num_rois);
        curr_glyc = zeros(1,FLUOR(1).num_rois); % 7/24/13
        %             curr_vglut(FLUOR(1).pos_alg_vglut) = 1; % 7/4/13
        %             curr_vglut(FLUOR(1).neg_alg_vglut) = -1; % 8/21/13
        
        %              curr_vglut(FLUOR(1).tp_vglut) = 1; % 7/10/14 - user-defined 'truth'
        %             curr_vglut(FLUOR(1).tn_vglut) = -1;
        
        %             curr_vglut(FLUOR(1).ML_pos) = 1; % 7/11/14 - user-defined 'truth'
        %             curr_vglut(FLUOR(1).ML_neg) = -1;
        
        %             curr_vglut(FLUOR(1).EA_pos2) = 1; % 7/18/14 - EA-defined 'truth'
        %             curr_vglut(FLUOR(1).EA_neg2) = -1;
        
        curr_vglut(FLUOR(1).EAML_pos) = 1; % 7/21/14 - EA and ML agreed upon 'truth'
        curr_vglut(FLUOR(1).EAML_neg) = -1;
        curr_vglut(FLUOR(1).EAML_unkn) = 2;
        curr_vglut(FLUOR(1).EAML_dis) = 3;
        
        % store stripe indices
        curr_stripe = zeros(4,FLUOR(1).num_rois);
        curr_stripe(1,MAP2.s1) = 1;
        curr_stripe(2,MAP2.s2) = 1;
        curr_stripe(3,MAP2.s3) = 1;
        curr_stripe(4,MAP2.s4) = 1;
        
        % store fish ID of fish -- single xfp
        if(ff==1)
            fishnum = 6; disp('fish 6');
            curr_age = 9*ones(1,FLUOR(1).num_rois);
        else
            fishnum = 7; disp('fish 7');
            curr_age = 12*ones(1,FLUOR(1).num_rois);
        end
        
        curr_fish = fishnum*ones(1,FLUOR(1).num_rois);
        
        curr_taus = zeros(1,FLUOR(1).num_rois);
        curr_taus1 = zeros(1,FLUOR(1).num_rois);
        curr_tauk = zeros(1,FLUOR(1).num_rois);
        curr_rsqs = zeros(1,FLUOR(1).num_rois);
        curr_rsqs1 = zeros(1,FLUOR(1).num_rois);
        curr_rsqk = zeros(1,FLUOR(1).num_rois);
        curr_rho = zeros(1,FLUOR(1).num_rois);
        curr_rho_fit = zeros(1,FLUOR(1).num_rois);
        curr_fit_res = zeros(length(SPONT4.t_fit),FLUOR(1).num_rois);
        curr_t_fit = repmat(SPONT4.t_fit',[1 FLUOR(1).num_rois]);
        %
        % store tau_s, tau_k, rsq values (may need to make this more
        % efficient)
        allrois = 1:FLUOR(1).num_rois;
        % use cnull that results in best fit
        [maxr,indm] = max(SPONT4.rsq1full);
        %[maxr,indm] = max(SPONT2.rsq(:,a1s1_rois));
        prs = [indm' allrois'];
        for k=1:FLUOR(1).num_rois
            curr_taus(k) = SPONT4.taup1(prs(k,1),prs(k,2));
            curr_rho_fit(k) = SPONT4.rho_fit(prs(k,1),prs(k,2));
            curr_rho(k) = SPONT4.rho(prs(k,1),prs(k,2));
            %addtaus(k) = SPONT2.taup(prs(k,1),prs(k,2));
            curr_fit_res(:,k) = squeeze(SPONT4.fit_res1(:,prs(k,1),prs(k,2)));
        end
        curr_rsqs = maxr;
        
        
        curr_tauk = OKS4.beta_all(1,:);
        curr_rsqk = OKS4.rsq;
        
        curr_ksens = OKS7.kfit;
        curr_Eth = OKS7.Eth;
        curr_r0 = OKS7.r0;
        
        curr_roi_frac = FLUOR(1).ROI_frac;
        curr_roi_ann = FLUOR(1).ROI_xfp_mean_vglut./FLUOR(1).ROI_ann_mean_vglut;
        curr_roi_asym = FLUOR(1).ROI_asym_vglut;
        
        % calculate ROI mean int
        curr_roi_mean_int = FLUOR(1).ROI_xfp_mean_vglut;
        %             curr_roi_mean_int = zeros(1,FLUOR(1).num_rois);
        %             curr_start = 1;
        %             for k=1:FLUOR(1).num_rois
        %                 curr_end = curr_start+FLUOR(1).len_patches(k)-1;
        %                 xi = FLUOR(1).xi_all(curr_start:curr_end);
        %                 yi = FLUOR(1).yi_all(curr_start:curr_end);
        %                 curr_start = curr_end + 1;
        %                 BW = poly2mask(xi,yi,size(FLUOR(1).ref_red,1),size(FLUOR(1).ref_red,1));
        %                 curr_roi_mean_int(k) = mean(FLUOR(1).ref_red(BW));
        %             end
        
        % save data
        fishid = [fishid,curr_fish];
        vpni = [vpni,curr_vpni];
        vsni = [vsni,curr_vsni];
        os = [os,curr_os];
        XYZ = [XYZ,curr_XYZ];
        stripe = [stripe,curr_stripe];
        vglut = [vglut,curr_vglut];
        glyc = [glyc,curr_glyc]; % 7/24/13
        taus = [taus,curr_taus];
        %taus1 = [taus1,curr_taus1];
        tauk = [tauk,curr_tauk];
        rsqs = [rsqs,curr_rsqs];
        %rsqs1 = [rsqs1,curr_rsqs1];
        rsqk = [rsqk,curr_rsqk];
        staf = [staf curr_staf];
        stad = [stad curr_stad];
        corrspont = [corrspont curr_corrspont];
        corrokr = [corrokr curr_corrokr];
        ksens = [ksens curr_ksens];
        Eth = [Eth curr_Eth];
        r0 = [r0 curr_r0];
        rho = [rho curr_rho];
        rho_fit = [rho_fit curr_rho_fit];
        rh67 = [rh67 curr_rh67];        % 7/2/14
        XYZadj = [XYZadj curr_XYZadj];
        roi_frac = [roi_frac curr_roi_frac];
        roi_ann = [roi_ann curr_roi_ann];
        roi_asym = [roi_asym curr_roi_asym];
        roi_mean_int = [roi_mean_int curr_roi_mean_int];
        age = [age curr_age];
        dstaf = [dstaf curr_dstaf];
        ndstaf = [ndstaf curr_ndstaf];
        fit_res = [fit_res curr_fit_res];
        t_fit = [t_fit curr_t_fit];
        LR = [LR,currLR];
        
        %   end  % if~....for excluding sets
    end
end


FIN.fishid = fishid;
FIN.vpni = vpni;
FIN.vsni = vsni;
FIN.os = os;
FIN.XYZ = XYZ;
FIN.stripe = stripe;
FIN.vglut = vglut;
FIN.glyc = glyc;
FIN.taus = taus;
%FIN.taus1 = taus1;
FIN.tauk = tauk;
FIN.rsqs = rsqs;
%FIN.rsqs1 = rsqs1;
FIN.rsqk = rsqk;
FIN.staf = staf;
FIN.stad = stad;
FIN.corrspont = corrspont;
FIN.corrokr = corrokr;
FIN.ksens = ksens;
FIN.Eth = Eth;
FIN.r0 = r0;
FIN.rgratio = rgratio;
FIN.rho = rho;
FIN.rho_fit = rho_fit;
FIN.rh67 = rh67;
FIN.XYZadj = XYZadj;
FIN.roi_frac = roi_frac;
FIN.roi_ann = roi_ann;
FIN.roi_asym = roi_asym;
FIN.roi_mean_int = roi_mean_int;
FIN.age = age;
FIN.dstaf = dstaf;
FIN.ndstaf = ndstaf;
FIN.fit_res = fit_res;
FIN.t_fit = t_fit;
FIN.LR = LR;


mainfolder = 'C:\Users\mel2011\Documents\Imaging\dual_xfp_expts';
cd(mainfolder);

%sname = 'dual_single_xfp_final_with_STA_082113.mat';    % lsqcurvefit, SPONT6, OKS6 (fit res saved), updated xfp class, functional class
%sname = 'dual_single_xfp_final_with_STA_082313.mat';    % SPONT8, OKS6
%sname = 'dual_single_xfp_final_with_STA_082613.mat';    % SPONT8, OKS6
%sname = 'dual_single_xfp_final_with_STA_082713.mat';    % SPONT8, OKS6
%sname = 'dual_single_xfp_final_with_STA_083113.mat';    % SPONT9, OKS6
%sname = 'dual_single_xfp_final_with_STA_090113.mat';    % SPONT9, OKS6
%sname = 'dual_single_xfp_final_with_STA_042514.mat';    % exclude rh6 cells
%sname = 'dual_single_xfp_final_with_STA_042514b.mat';    % exclude rh6 cells + borderline ones
%sname = 'dual_single_xfp_final_with_STA_070214.mat';    % exclude rh6 cells + borderline ones; adjust RC position rel to rh67 border
%sname = 'dual_single_xfp_final_with_STA_071014.mat';    % for XFP id, use user-ID'd truth set
%sname = 'dual_single_xfp_final_with_STA_071114.mat';    % for XFP id, use ML-ID'd truth set
%sname = 'dual_single_xfp_final_with_STA_071814.mat';    % for XFP id, use EA-ID'd truth set
%sname = 'dual_single_xfp_final_with_STA_072114.mat';    % for XFP id, use EA and ML-ID'd truth set
%sname = 'dual_single_xfp_final_with_STA_102814.mat';    % for XFP id, use EA and ML-ID'd truth set; added age of fish
%sname = 'dual_single_xfp_final_with_STA_112014.mat';    % for XFP id, use EA and ML-ID'd truth set; added age, dstaf, fit_res
sname = 'dual_single_xfp_final_with_STA_020315.mat';    % for XFP id, use EA and ML-ID'd truth set; added LR


save(sname,'FIN');


% % calculate statistics for paper (8/27/13): RHO metric
% % ----------------------------------------
% close all; clear all;
% mainfolder = 'C:\Users\mel2011\Documents\Imaging\dual_xfp_expts';
% cd(mainfolder);
% %sname = 'dual_single_xfp_final_with_STA_090113.mat';    % SPONT9, OKS6
% %sname = 'dual_single_xfp_final_with_STA_042514.mat';    % exclude rh6 cell
% %sname = 'dual_single_xfp_final_with_STA_042514b.mat';    % exclude rh6 cells + borderline ones
% %sname = 'dual_single_xfp_final_with_STA_070214.mat';    % exclude rh6 cells + borderline ones; adjust RC position rel to rh67 border
% %sname = 'dual_single_xfp_final_with_STA_071114.mat';    % for XFP id, use ML-ID'd truth set
% %sname = 'dual_single_xfp_final_with_STA_071814.mat';    % for XFP id, use EA-ID'd truth set
% sname = 'dual_single_xfp_final_with_STA_072114.mat';    % for XFP id, use EA and ML-ID'd truth set
%
% load(sname);
%
% % %a1tau = FIN.rho(find(FIN.vpni & FIN.rsqs(4,:) > 0.97)); % 8/26/13
% % a1tau = FIN.taus(find(FIN.vpni & FIN.rsqs > 0.97)); % 8/27/13
% % a1tau = min(a1tau,100);
% % a1tau = max(a1tau,1);
% % %num_a1 = length(a1tau)
% % mean_a1tau = mean(a1tau)
% % med_a1tau = median(a1tau)
% % min_a1tau = min(a1tau)
% % max_a1tau = max(a1tau)
%
% a1rho = FIN.rho_fit(find(FIN.vpni & FIN.rsqs > 0.97));
% %a1rho = FIN.rho(find(FIN.vpni & FIN.rsqs > 0.97));
% %a1rho = min(a1rho,2);   % cap at 2
% a1rho(a1rho>2)=[];  % discard rho>2
% num_a1 = length(a1rho);
% mean_a1rho = mean(a1rho)
% med_a1rho = median(a1rho)
% % min_a1rho = min(a1rho)
% % max_a1rho = max(a1rho)
% len20 = round(.2*num_a1);
% a1rho = sort(a1rho);
% min_a1rho = mean(a1rho(1:len20))
% max_a1rho = mean(a1rho(num_a1-len20:end))
%
% % write to excel file
% %a1vect = [num_a1 mean_a1tau med_a1tau min_a1tau max_a1tau mean_a1rho med_a1rho min_a1rho max_a1rho]';
% a1vect = [num_a1 mean_a1rho med_a1rho min_a1rho max_a1rho]';
% %xlswrite('tau_rho.xls',a1vect,'A1:A9');
%
% %figure; plot(a1tau,a1rho,'.'); xlabel('\tau'); ylabel('\rho'); title('all VPNI');
%
% %figure; hist(a1tau); xlabel('\tau'); ylabel('count'); title('all VPNI');
% figure; hist(a1rho); xlabel('\rho'); ylabel('count'); title('all VPNI');
%
% % 10/30/13: for cumulative histogram
% [n,bins] = hist(a1rho);
% nelt = histc(a1rho,bins);
% celt = cumsum(nelt);
% %figure; bar(bins,celt); xlabel('\rho'); ylabel('cumulative distribution'); title('all VPNI');
%
% % DsRed positive group
% %a1co_tau = FIN.rho(find(FIN.vpni&FIN.vglut==1& FIN.rsqs(4,:) > 0.97)); % 8/26/13
% a1co_tau = FIN.taus(find(FIN.vpni&FIN.vglut==1& FIN.rsqs > 0.97)); % 8/27/13
% a1co_tau = min(a1co_tau,100);
% a1co_tau = max(a1co_tau,1);
% %num_a1co = length(a1co_tau)
% mean_a1tau = mean(a1co_tau)
% med_a1tau = median(a1co_tau)
% min_a1tau = min(a1co_tau)
% max_a1tau = max(a1co_tau)
%
% a1co_rho = FIN.rho_fit(find(FIN.vpni&FIN.vglut==1& FIN.rsqs > 0.97));
% %a1co_rho = FIN.rho(find(FIN.vpni&FIN.vglut==1& FIN.rsqs > 0.97));
% %a1co_rho = min(a1co_rho,2);   % cap at 2
% a1co_rho(a1co_rho>2)=[];  % discard rho>2
% num_a1co = length(a1co_rho)
% mean_a1rho = mean(a1co_rho)
% med_a1rho = median(a1co_rho)
% % min_a1rho = min(a1co_rho)
% % max_a1rho = max(a1co_rho)
% len20 = round(.2*num_a1co);
% a1co_rho = sort(a1co_rho);
% min_a1rho = mean(a1co_rho(1:len20))
% max_a1rho = mean(a1co_rho(num_a1co-len20:end))
%
% %a1covect = [num_a1co mean_a1tau med_a1tau min_a1tau max_a1tau mean_a1rho med_a1rho min_a1rho max_a1rho]';
% a1covect = [num_a1co mean_a1rho med_a1rho min_a1rho max_a1rho]';
%
% %figure; hist(a1co_tau); xlabel('\tau'); ylabel('count'); title('DsRed+ VPNI');
% figure; hist(a1co_rho); xlabel('\rho'); ylabel('count'); title('DsRed+ VPNI');
%
% % 10/30/13: for cumulative histogram
% [nco,binsco] = hist(a1co_rho);
% neltco = histc(a1co_rho,binsco);
% celtco = cumsum(neltco);
% %figure; bar(binsco,celtco); xlabel('\rho'); ylabel('cumulative distribution'); title('DsRed+ VPNI');
%
%
% % DsRed negative group
% a1non_tau = FIN.taus(find(FIN.vpni&FIN.vglut==-1& FIN.rsqs > 0.97)); % 8/27/13
% a1non_tau = min(a1non_tau,100);
% a1non_tau = max(a1non_tau,1);
% %num_a1non = length(a1non_tau)
% mean_a1tau = mean(a1non_tau)
% med_a1tau = median(a1non_tau)
% min_a1tau = min(a1non_tau)
% max_a1tau = max(a1non_tau)
%
% a1non_rho = FIN.rho_fit(find(FIN.vpni&FIN.vglut==-1& FIN.rsqs > 0.97));
% %a1non_rho = FIN.rho(find(FIN.vpni&FIN.vglut==-1& FIN.rsqs > 0.97));
% %a1non_rho = min(a1non_rho,2);   % cap at 2
% a1non_rho(a1non_rho>2) = [];    % discard if  >2
% num_a1non = length(a1non_rho)
% mean_a1rho = mean(a1non_rho)
% med_a1rho = median(a1non_rho)
% % min_a1rho = min(a1non_rho)
% % max_a1rho = max(a1non_rho)
% len20 = round(.2*num_a1non);
% a1non_rho = sort(a1non_rho);
% min_a1rho = mean(a1non_rho(1:len20))
% max_a1rho = mean(a1non_rho(num_a1non-len20:end))
%
% %a1nonvect = [num_a1non mean_a1tau med_a1tau min_a1tau max_a1tau mean_a1rho med_a1rho min_a1rho max_a1rho]';
% a1nonvect = [num_a1non mean_a1rho med_a1rho min_a1rho max_a1rho]';
%
% %figure; hist(a1co_tau); title('DsRed+ VPNI');
% % figure; plot(a1co_tau,a1co_rho,'r.',a1non_tau,a1non_rho,'g.'); xlabel('\tau'); ylabel('\rho');
% % legend('DsRed+ VPNI','DsRed- VPNI');
%
% %figure; hist(a1non_tau); xlabel('\tau'); ylabel('count'); title('DsRed- VPNI');
% figure; hist(a1non_rho); xlabel('\rho'); ylabel('count'); title('DsRed- VPNI');
%
% % 10/30/13: for cumulative histogram
% %[nnon,binsnon] = hist(a1non_rho);
% [nnon,binsnon] = hist(a1non_rho,100); % 6/23/14: use 100 bins
% neltnon = histc(a1non_rho,binsnon);
% celtnon = cumsum(neltnon);
% %figure; bar(binsnon,celtnon); xlabel('\rho'); ylabel('cumulative distribution'); title('DsRed- VPNI');
%
% %a1co_non_comp = [ht pt hr pr]';
%
% % glut stripe 1
% %a1s1_tau = FIN.rho(find(FIN.vpni & FIN.vglut==1 & FIN.stripe(1,:) & FIN.rsqs(4,:) > 0.97)); % 8/26/13
% a1s1_tau = FIN.taus(find(FIN.vpni & FIN.vglut==1 & FIN.stripe(1,:) & FIN.rsqs > 0.97)); % 8/27/13
% a1s1_tau = min(a1s1_tau,100);
% a1s1_tau = max(a1s1_tau,1);
% %num_a1s1 = length(a1s1_tau)
% mean_a1tau = mean(a1s1_tau)
% med_a1tau = median(a1s1_tau)
% min_a1tau = min(a1s1_tau)
% max_a1tau = max(a1s1_tau)
%
% a1s1_rho = FIN.rho_fit(find(FIN.vpni & FIN.vglut==1 & FIN.stripe(1,:) & FIN.rsqs > 0.97));
% %a1s1_rho = FIN.rho(find(FIN.vpni & FIN.vglut==1 & FIN.stripe(1,:) & FIN.rsqs > 0.97));
% %a1s1_rho = min(a1s1_rho,2);   % cap at 2
% a1s1_rho(a1s1_rho>2) =[];   % discard if >2
% num_a1s1 = length(a1s1_rho)
% mean_a1rho = mean(a1s1_rho)
% med_a1rho = median(a1s1_rho)
% % min_a1rho = min(a1s1_rho)
% % max_a1rho = max(a1s1_rho)
% len20 = round(.2*num_a1s1);
% a1s1_rho = sort(a1s1_rho);
% min_a1rho = mean(a1s1_rho(1:len20))
% max_a1rho = mean(a1s1_rho(num_a1s1-len20:end))
%
% %a1s1vect = [num_a1s1 mean_a1tau med_a1tau min_a1tau max_a1tau mean_a1rho med_a1rho min_a1rho max_a1rho]';
% a1s1vect = [num_a1s1 mean_a1rho med_a1rho min_a1rho max_a1rho]';
%
% % non-glut stripe 1
% a1s1n_tau = FIN.taus(find(FIN.vpni & FIN.vglut==-1 & FIN.stripe(1,:) & FIN.rsqs >0.97)); % 8/27/13
% a1s1n_tau = min(a1s1n_tau,100);
% a1s1n_tau = max(a1s1n_tau,1);
% %num_a1s1n = length(a1s1n_tau)
% mean_a1tau = mean(a1s1n_tau)
% med_a1tau = median(a1s1n_tau)
% min_a1tau = min(a1s1n_tau)
% max_a1tau = max(a1s1n_tau)
%
% a1s1n_rho = FIN.rho_fit(find(FIN.vpni & FIN.vglut==-1 & FIN.stripe(1,:) & FIN.rsqs > 0.97));
% a1s1n_rho(a1s1n_rho>2) = [];    % discard if >2
% num_a1s1n = length(a1s1n_rho)
% mean_a1rho = mean(a1s1n_rho)
% med_a1rho = median(a1s1n_rho)
% % min_a1rho = min(a1s2n_rho)
% % max_a1rho = max(a1s2n_rho)
% len20 = round(.2*num_a1s1n);
% a1s1n_rho = sort(a1s1n_rho);
% min_a1rho = mean(a1s1n_rho(1:len20))
% max_a1rho = mean(a1s1n_rho(num_a1s1n-len20:end))
%
% %a1s2nvect = [num_a1s2n mean_a1tau med_a1tau min_a1tau max_a1tau mean_a1rho med_a1rho min_a1rho max_a1rho]';
% a1s1nvect = [num_a1s1n mean_a1rho med_a1rho min_a1rho max_a1rho]';
%
% figure; hist(a1s1n_rho); xlabel('\rho'); ylabel('count'); title('DsRed- stripe 1');
%
% % 10/30/13: for cumulative histogram
% %[ns1n,binss1n] = hist(a1s1n_rho);
% [ns1n,binss1n] = hist(a1s1n_rho,100); % 6/23/14: use 100 bins
% nelts1n = histc(a1s1n_rho,binss1n);
% celts1n = cumsum(nelts1n);
% %figure; bar(binss1n,celts1n); xlabel('\rho'); ylabel('cumulative distribution'); title('DsRed- stripe 1');
%
% %[pt,ht] = ranksum(a1s1_tau,a1s1n_tau)   % 9/23/13
% [pr,hr] = ranksum(a1s1_rho,a1s1n_rho)
%
% %a1s1_s1n_comp = [ht pt hr pr]';
%
%
% % glut stripe 2
% %a1s2_tau = FIN.rho(find(FIN.vpni & FIN.vglut==1 & FIN.stripe(2,:) & FIN.rsqs(4,:) >0.97)); % 8/26/13
% a1s2_tau = FIN.taus(find(FIN.vpni & FIN.vglut==1 & FIN.stripe(2,:) & FIN.rsqs >0.97)); % 8/27/13
% a1s2_tau = min(a1s2_tau,100);
% a1s2_tau = max(a1s2_tau,1);
% %num_a1s2 = length(a1s2_tau)
% mean_a1tau = mean(a1s2_tau)
% med_a1tau = median(a1s2_tau)
% min_a1tau = min(a1s2_tau)
% max_a1tau = max(a1s2_tau)
%
% a1s2_rho = FIN.rho_fit(find(FIN.vpni & FIN.vglut==1 & FIN.stripe(2,:) & FIN.rsqs > 0.97));
% %a1s2_rho = FIN.rho(find(FIN.vpni & FIN.vglut==1 & FIN.stripe(2,:) & FIN.rsqs > 0.97));
% %a1s2_rho = min(a1s2_rho,2);   % cap at 2
% a1s2_rho(a1s2_rho>2) =[];   % discard if >2
% num_a1s2 = length(a1s2_rho)
% mean_a1rho = mean(a1s2_rho)
% med_a1rho = median(a1s2_rho)
% % min_a1rho = min(a1s2_rho)
% % max_a1rho = max(a1s2_rho)
% len20 = round(.2*num_a1s2);
% a1s2_rho = sort(a1s2_rho);
% min_a1rho = mean(a1s2_rho(1:len20))
% max_a1rho = mean(a1s2_rho(num_a1s2-len20:end))
%
% %a1s2vect = [num_a1s2 mean_a1tau med_a1tau min_a1tau max_a1tau mean_a1rho med_a1rho min_a1rho max_a1rho]';
% a1s2vect = [num_a1s2 mean_a1rho med_a1rho min_a1rho max_a1rho]';
%
% % figure; plot(a1s1_tau,a1s1_rho,'r.',a1s2_tau,a1s2_rho,'m.'); xlabel('\tau'); ylabel('\rho');
% % legend('stripe 1 glut','stripe 2 glut');
%
% %figure; hist(a1s1_tau); xlabel('\tau'); ylabel('count'); title('DsRed+ stripe 1');
% figure; hist(a1s1_rho); xlabel('\rho'); ylabel('count'); title('DsRed+ stripe 1');
%
% % 10/30/13: for cumulative histogram
% %[ns1,binss1] = hist(a1s1_rho);
% [ns1,binss1] = hist(a1s1_rho,100); % 6/23/14: use 100 bins
% nelts1 = histc(a1s1_rho,binss1);
% celts1 = cumsum(nelts1);
% %figure; bar(binss1,celts1); xlabel('\rho'); ylabel('cumulative distribution'); title('DsRed+ stripe 1');
%
% %[pt,ht] = ranksum(a1s1_tau,a1s2_tau) % 9/23/13
% [pr,hr] = ranksum(a1s1_rho,a1s2_rho)
% [pr,hr] = ranksum([a1s1_rho,a1s1n_rho],a1s2_rho) % 5/10/14: compare all s1 cells to s2 gluts
%
% %a1s1_s2_comp = [ht pt hr pr]';
%
% % non-glut stripe 2
% %a1s2n_tau = FIN.rho(find(FIN.vpni & FIN.vglut==-1 & FIN.stripe(2,:) & FIN.rsqs(4,:) >0.97)); % 8/26/13
% a1s2n_tau = FIN.taus(find(FIN.vpni & FIN.vglut==-1 & FIN.stripe(2,:) & FIN.rsqs >0.97)); % 8/27/13
% a1s2n_tau = min(a1s2n_tau,100);
% a1s2n_tau = max(a1s2n_tau,1);
% %num_a1s2n = length(a1s2n_tau)
% mean_a1tau = mean(a1s2n_tau)
% med_a1tau = median(a1s2n_tau)
% min_a1tau = min(a1s2n_tau)
% max_a1tau = max(a1s2n_tau)
%
% a1s2n_rho = FIN.rho_fit(find(FIN.vpni & FIN.vglut==-1 & FIN.stripe(2,:) & FIN.rsqs > 0.97));
% %a1s2n_rho = FIN.rho(find(FIN.vpni & FIN.vglut==-1 & FIN.stripe(2,:) & FIN.rsqs > 0.97));
% %a1s2n_rho = min(a1s2n_rho,2);   % cap at 2
% a1s2n_rho(a1s2n_rho>2) = [];    % discard if >2
% num_a1s2n = length(a1s2n_rho)
% mean_a1rho = mean(a1s2n_rho)
% med_a1rho = median(a1s2n_rho)
% % min_a1rho = min(a1s2n_rho)
% % max_a1rho = max(a1s2n_rho)
% len20 = round(.2*num_a1s2n);
% a1s2n_rho = sort(a1s2n_rho);
% min_a1rho = mean(a1s2n_rho(1:len20))
% max_a1rho = mean(a1s2n_rho(num_a1s2n-len20:end))
%
% %a1s2nvect = [num_a1s2n mean_a1tau med_a1tau min_a1tau max_a1tau mean_a1rho med_a1rho min_a1rho max_a1rho]';
% a1s2nvect = [num_a1s2n mean_a1rho med_a1rho min_a1rho max_a1rho]';
%
% % figure; plot(a1s2n_tau,a1s2n_rho,'g.',a1s2_tau,a1s2_rho,'m.'); xlabel('\tau'); ylabel('\rho');
% % legend('stripe 2 GABA','stripe 2 glut');
%
% %figure; hist(a1s2_tau); xlabel('\tau'); ylabel('count'); title('DsRed+ stripe 2');
% figure; hist(a1s2_rho); xlabel('\rho'); ylabel('count'); title('DsRed+ stripe 2');
%
% % 10/30/13: for cumulative histogram
% %[ns2,binss2] = hist(a1s2_rho);
% [ns2,binss2] = hist(a1s2_rho,100); % 6/23/14: use 100 bins
% nelts2 = histc(a1s2_rho,binss2);
% celts2 = cumsum(nelts2);
% %figure; bar(binss2,celts2); xlabel('\rho'); ylabel('cumulative distribution'); title('DsRed+ stripe 2');
%
%
% %figure; hist(a1s2n_tau); xlabel('\tau'); ylabel('count'); title('DsRed- stripe 2');
% figure; hist(a1s2n_rho); xlabel('\rho'); ylabel('count'); title('DsRed- stripe 2');
%
% % 10/30/13: for cumulative histogram
% %[ns2n,binss2n] = hist(a1s2n_rho);
% [ns2n,binss2n] = hist(a1s2n_rho,100); % 3/23/14: use 100 bins
% nelts2n = histc(a1s2n_rho,binss2n);
% celts2n = cumsum(nelts2n);
% % figure; bar(binss2n,celts2n); xlabel('\rho'); ylabel('cumulative distribution'); title('DsRed- stripe 2');
%
% % % all cumulative histograms
% % figure; stairs(bins,celt,'k','LineWidth',2); hold on;
% % stairs(binss1,celts1,'r','LineWidth',2); hold on;
% % stairs(binss1n,celts1n,'b','LineWidth',2); hold on;
% % stairs(binss2,celts2,'m','LineWidth',2); hold on;
% % stairs(binss2n,celts2n,'g','LineWidth',2);
% % xlabel('\rho'); ylabel('cumulative distribution');
% % legend('all','stripe 1 +','stripe 1 -','stripe 2 +','stripe 2 -');
%
% % % all cumulative histograms, excluding all VPNI popn
% % figure;
% % stairs(binss1,celts1,'r','LineWidth',2); hold on;
% % stairs(binss1n,celts1n,'b','LineWidth',2); hold on;
% % stairs(binss2,celts2,'m','LineWidth',2); hold on;
% % stairs(binss2n,celts2n,'g','LineWidth',2);
% % xlabel('\rho'); ylabel('cumulative distribution');
% % legend('stripe 1 +','stripe 1 -','stripe 2 +','stripe 2 -');
%
% % % all cumulative histograms, excluding all VPNI popn, normalized
% % figure;
% % stairs(binss1,celts1./max(celts1),'r','LineWidth',2); hold on;
% % stairs(binss1n,celts1n./max(celts1n),'b','LineWidth',2); hold on;
% % stairs(binss2,celts2./max(celts2),'m','LineWidth',2); hold on;
% % stairs(binss2n,celts2n./max(celts2n),'g','LineWidth',2);
% % xlabel('\rho'); ylabel('cumulative distribution');
% % legend('stripe 1 +','stripe 1 -','stripe 2 +','stripe 2 -');
%
% % all cumulative histograms, excluding all VPNI popn, normalized -- add 0
% % (6/23/14)
% figure;
% stairs([0,binss1],[0,celts1./max(celts1)],'r','LineWidth',2); hold on;
% stairs([0,binss1n],[0,celts1n./max(celts1n)],'b','LineWidth',2); hold on;
% stairs([0,binss2],[0,celts2./max(celts2)],'m','LineWidth',2); hold on;
% stairs([0,binss2n],[0,celts2n./max(celts2n)],'g','LineWidth',2);
% xlabel('\rho'); ylabel('cumulative distribution');
% legend('stripe 1 +','stripe 1 -','stripe 2 +','stripe 2 -');
%
% %[pt,ht] = ranksum(a1s2_tau,a1s2n_tau) % 9/23/13
% [pr,hr] = ranksum(a1s2_rho,a1s2n_rho)
%
% [p,h] = ranksum([a1s1_rho,a1s2n_rho],a1s2n_rho)  % 5/10/14: compare all s1 with s2 GABAs
%
% %a1s2_s2n_comp = [ht pt hr pr]';
%
%
%
% % generate histogram summary for all populations
% [a1rho_n,binc] = hist(a1rho);
% [a1s1rho_n] = hist(a1s1_rho,binc);
% [a1s1nrho_n] = hist(a1s1n_rho,binc);
% [a1s2rho_n] = hist(a1s2_rho,binc);
% [a1s2nrho_n] = hist(a1s2n_rho,binc);
% % figure; plot(binc,a1rho_n,'k.-',binc,a1s1rho_n,'r.-',binc,a1s1nrho_n,'b.-',binc,a1s2rho_n,'m.-',binc,a1s2nrho_n,'g.-','markersize',20);
% % xlabel('\rho','FontName','MyriadPro-Regular','FontSize',14);
% % ylabel('# cells','FontName','MyriadPro-Regular','FontSize',14);
% % legend('all VPNI','stripe 1 DsRed+','stripe 1 DsRed-','stripe 2 DsRed+','stripe 2 DsRed-');
% % set(gca,'box','off','TickDir','out','FontName','MyriadPro-Regular','FontSize',12);
% %
% % figure; plot(binc,a1s1rho_n,'r.-',binc,a1s1nrho_n,'b.-',binc,a1s2rho_n,'m.-',binc,a1s2nrho_n,'g.-','markersize',20);
% % xlabel('\rho','FontName','MyriadPro-Regular','FontSize',14);
% % ylabel('# cells','FontName','MyriadPro-Regular','FontSize',14);
% % legend('stripe 1 DsRed+','stripe 1 DsRed-','stripe 2 DsRed+','stripe 2 DsRed-');
% % set(gca,'box','off','TickDir','out','FontName','MyriadPro-Regular','FontSize',12);
%
% % % manhattan/skyscraper plot
% % figure; stairs(binc,a1rho_n,'k','LineWidth',2); hold on;
% % stairs(binc,a1s1rho_n,'r','LineWidth',2); hold on;
% % stairs(binc,a1s1nrho_n,'b','LineWidth',2); hold on;
% % stairs(binc,a1s2rho_n,'m','LineWidth',2); hold on;
% % stairs(binc,a1s2nrho_n,'g','LineWidth',2);
% % xlabel('\rho','FontName','MyriadPro-Regular','FontSize',14);
% % ylabel('# cells','FontName','MyriadPro-Regular','FontSize',14);
% % legend('all VPNI','stripe 1 DsRed+','stripe 1 DsRed-','stripe 2 DsRed+','stripe 2 DsRed-');
% % set(gca,'box','off','TickDir','out','FontName','MyriadPro-Regular','FontSize',12);
% %
% % figure; hold on;
% % stairs(binc,a1s1rho_n,'r','LineWidth',2); hold on;
% % stairs(binc,a1s1nrho_n,'b','LineWidth',2); hold on;
% % stairs(binc,a1s2rho_n,'m','LineWidth',2); hold on;
% % stairs(binc,a1s2nrho_n,'g','LineWidth',2);
% % xlabel('\rho','FontName','MyriadPro-Regular','FontSize',14);
% % ylabel('# cells','FontName','MyriadPro-Regular','FontSize',14);
% % legend('stripe 1 DsRed+','stripe 1 DsRed-','stripe 2 DsRed+','stripe 2 DsRed-');
% % set(gca,'box','off','TickDir','out','FontName','MyriadPro-Regular','FontSize',12);
%
% % manhattan plot for all VPNIcells
% [binc,a1rho_n] = stairs(binc,a1rho_n);
% % close off ends... (to find d*, enter: diff(*))
% dbinc = 0.1872;
% binc = [binc(1);binc;binc(end)+dbinc;binc(end)+dbinc];
% a1rho_n = [0;a1rho_n;a1rho_n(end);0];
% figure; plot(binc,a1rho_n,'k','LineWidth',2);
% xlabel('\rho','FontName','MyriadPro-Regular','FontSize',14);
% ylabel('# cells','FontName','MyriadPro-Regular','FontSize',14);
% set(gca,'box','off','TickDir','out','FontName','MyriadPro-Regular','FontSize',12);
%
% [a1rho_n,binc] = hist(a1rho);
% [a1s1rho_n] = hist(a1s1_rho,binc);
% [a1s1nrho_n] = hist(a1s1n_rho,binc);
% [a1s2rho_n] = hist(a1s2_rho,binc);
% [a1s2nrho_n] = hist(a1s2n_rho,binc);
%
%
% % make normalized manhattan plots
% figure; hold on;
% [a1s1_x,a1s1_y] = stairs(binc,a1s1rho_n);
% [a1s1n_x,a1s1n_y] = stairs(binc,a1s1nrho_n);
% % close off ends... (to find d*, enter: diff(*))
% dbinc = 0.1872;
% a1s1_x = [a1s1_x(1);a1s1_x;a1s1_x(end)+dbinc;a1s1_x(end)+dbinc];
% a1s1n_x = [a1s1n_x(1);a1s1n_x;a1s1n_x(end)+dbinc;a1s1n_x(end)+dbinc];
% a1s1_y = [0;a1s1_y;a1s1_y(end);0];
% a1s1n_y = [0;a1s1n_y;a1s1n_y(end);0];
% [ax,h1,h2] = plotyy(a1s1_x,a1s1_y,a1s1n_x,a1s1n_y,'plot');
% set(h1,'LineWidth',2,'Color','r'); %set(ax,'YLim',[0 20]);
% set(h2,'LineWidth',2,'Color','b');
% xlabel('\rho','FontName','MyriadPro-Regular','FontSize',14);
% ylabel('# cells','FontName','MyriadPro-Regular','FontSize',14);
% legend('stripe 1 DsRed+','stripe 1 DsRed-');
% set(gca,'box','off','TickDir','out','FontName','MyriadPro-Regular','FontSize',12);
%
% figure; hold on;
% [a1s2_x,a1s2_y] = stairs(binc,a1s2rho_n);
% [a1s2n_x,a1s2n_y] = stairs(binc,a1s2nrho_n);
% % close off ends... (to find d*, enter: diff(*))
% dbinc = 0.1872; % diff(binc)
% a1s2_x = [a1s2_x(1);a1s2_x;a1s2_x(end)+dbinc;a1s2_x(end)+dbinc];
% a1s2n_x = [a1s2n_x(1);a1s2n_x;a1s2n_x(end)+dbinc;a1s2n_x(end)+dbinc];
% a1s2_y = [0;a1s2_y;a1s2_y(end);0];
% a1s2n_y = [0;a1s2n_y;a1s2n_y(end);0];
% [ax,h1,h2] = plotyy(a1s2_x,a1s2_y,a1s2n_x,a1s2n_y,'plot');
% set(h1,'LineWidth',2,'Color','m');
% set(h2,'LineWidth',2,'Color','g');
% xlabel('\rho','FontName','MyriadPro-Regular','FontSize',14);
% ylabel('# cells','FontName','MyriadPro-Regular','FontSize',14);
% legend('stripe 2 DsRed+','stripe 2 DsRed-');
% set(gca,'box','off','TickDir','out','FontName','MyriadPro-Regular','FontSize',12);
%
%
% % % write data to excel file
% % % xlswrite('tau_rho.xls',a1vect,'B2:B10');
% % % xlswrite('tau_rho.xls',a1covect,'C2:C10');
% % % xlswrite('tau_rho.xls',a1nonvect,'D2:D10');
% % % xlswrite('tau_rho.xls',a1s1vect,'E2:E10');
% % % xlswrite('tau_rho.xls',a1s2vect,'F2:F10');
% % % xlswrite('tau_rho.xls',a1s2nvect,'G2:G10');
% % % xlswrite('tau_rho.xls',a1co_non_comp,'B15:B18');
% % % xlswrite('tau_rho.xls',a1s1_s2_comp,'C15:C18');
% % % xlswrite('tau_rho.xls',a1s2_s2n_comp,'D15:D18');
% %
% % xlswrite('tau_rho.xls',a1vect,3,'B2:B6');
% % xlswrite('tau_rho.xls',a1covect,3,'C2:C6');
% % xlswrite('tau_rho.xls',a1nonvect,3,'D2:D6');
% % xlswrite('tau_rho.xls',a1s1vect,3,'E2:E6');
% % xlswrite('tau_rho.xls',a1s1nvect,3,'F2:F6');
% % xlswrite('tau_rho.xls',a1s2vect,3,'G2:G6');
% % xlswrite('tau_rho.xls',a1s2nvect,3,'H2:H6');
% % xlswrite('tau_rho.xls',a1co_non_comp,3,'B15:B18');
% % xlswrite('tau_rho.xls',a1s1_s2_comp,3,'C15:C18');
% % xlswrite('tau_rho.xls',a1s2_s2n_comp,3,'D15:D18');
% % xlswrite('tau_rho.xls',a1s1_s1n_comp,3,'E15:E18');
%
%
