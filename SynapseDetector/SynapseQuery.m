clc;
clear all;

df = readtable('final.csv');

% Columns 1 through 12
%     'psd_segid'    'BBOX_bx'    'BBOX_by'    'BBOX_bz'    'BBOX_ex'    'BBOX_ey'    'BBOX_ez'    'COM_x'    'COM_y'    'COM_z'    'postsyn_seg'    'postsyn_sz'
%  Columns 13 through 23
%     'postsyn_wt'    'postsyn_x'    'postsyn_y'    'postsyn_z'    'presyn_seg'    'presyn_sz'    'presyn_wt'    'presyn_x'    'presyn_y'    'presyn_z'    'size'

eg = (df.presyn_seg(df.postsyn_seg==76181));% find presynaptic partners of a cell

eg = (df.postsyn_seg(df.presyn_seg==76181));% find postsynaptic partners of a cell

% find location of synapses that need to be seeded
PresynapticPartners = (df.presyn_seg(df.postsyn_seg==76184));
PrePartnerSynapseSegID76184 = PresynapticPartners(PresynapticPartners>1e5);
PrePartnerCoordinatesID76184 = PrePartnerCoordinates(PrePartnerSynapseSegID76184,df);
PostSynapticPartners = 
PostPartnerCoordinatesID76184 = PostPartnerCoordinates();

clear temp
temp = (df.presyn_seg(df.postsyn_seg==76487))

