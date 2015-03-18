%EXAMPLE_SCRIPT_3D - demonstrates usage of segmentation code with 3d inputs
%
% Authors: Alessandro Giusti and Dan Claudiu Cire?an
% Dalle Molle Institute for Artificial Intelligence (IDSIA)
% Lugano, Switzerland
%
% Contact: alessandrog@idsia.ch - dan@idsia.ch
%
% FOR RESEARCH USE ONLY, DO NOT REDISTRIBUTE
%
% Oct 2012; Last revision: 26-October-2012

addpath apply-dnn

% Read source stack
%in_fn=fullfile('data','input.tif');
in_fn = fullfile(filesep,'home', 'ashwin', 'av_data', '07122012','Exported-W007-CompleteTest', 'Normalizedimg001_x10_y1.tif' );
%in_fn=fullfile('data','Exported_Stack_Area2_05202012D2_1024x1024_10nm.tif');
for i=1:length(imfinfo(in_fn))
   ims(:,:,i)=imread(in_fn,i);
end
% note: im is expected to be uint8 grayscale (range [0-255]) or double
%ims(:,:,1)=imread(in_fn,1);


% grayscale (range [0-1]).

% Read net
fprintf('Reading net...\n');

% 2d input
net=readnet(fullfile('nets','2d','53x53-48C4-MP2-48C4-MP2-48C4-MP2-100N-1355827863-BVerr-V-10.770-T-11.157-epoch-14.nnt'));
% 3d input
net=readnet(fullfile('nets','7HL-37x37-48C4-MP2-48C4-MP2-48C4-MP2-100N-1355785981-BVerr-V-10.502-T-10.720-epoch-33.nnt'));

% Process stack
fprintf('Processing stack...\n');
t=tic;
tic;

% consider net inputs
nmaps=net.layers(1).l.nMaps;
assert(mod(nmaps,2)==1); %% must have odd number of maps
dz=(-(nmaps-1)/2):(+(nmaps-1)/2);
assert(length(dz)==nmaps);

% Prepare filtering kernels and function
dom1=fspecial('disk',3)>0.01; ord1=ceil(length(find(dom1))/2);
dom2=fspecial('disk',3)>0.03; ord2=ceil(length(find(dom2))/2);
ffunc=@(im) mean(cat(3,ordfilt2(im,ord1,dom1),ordfilt2(im,ord2,dom2)),3);

outdir=['out-net' net.name];
mkdir(outdir);

for i=1:size(ims,3) % changed to 2 from 1
    % decide which slices to sample for output slice i
    % (slices out of the stack are replicated by extending the closest one)
    zs=min(size(ims,3),max(1,i+dz));
    % build input image (one channel per input map)
    im=ims(:,:,zs);
    
    % actually process
    maxpixels=0.5e6; showimages=false;
    out=processImageWithNet(im,net,maxpixels,showimages);                       
    
    % make sure that output is as we expect
    assert(size(im,1)==size(out,1));
    assert(size(im,2)==size(out,2));
    sums=sum(out,3);
    assert(all(abs(sums(:)-1)<1e-5));
    membrane=out(:,:,2);
    imwrite(membrane,fullfile(outdir,sprintf('out-slice%05g.png',i)));
    
    % Apply filtering function
    membrane_filtered=ffunc(membrane);
    imwrite(membrane_filtered,fullfile(outdir,sprintf('out-slice%05g-filtered.png',i)));
end
toc;

