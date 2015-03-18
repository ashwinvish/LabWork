%
% Function to generate borders between objects in a color-labeled 3D tiff image
% 
% Parameters: 
%   numSlices - number of slices in the 3D tiff to be processed
%	numThreads - number of threads to use (maximum is 8)
%	inputFileName - name of the input 3D tiff
%	outputFileName - name of the output 3d tiff
%
function [output] = generateBorders(numSlices, numThreads, inputFileName, outputFileName )

for kk = 1:numSlices
    input(:,:,kk) = imread(inputFileName, 'Index',kk);
end

output = input;

% maximum number of threads is 8
numThreads = min( numThreads, 8 );

matlabpool( numThreads );

parfor slice = 1:numSlices
    
    currentImage = input(:,:, slice);
    thisImage = input(:,:, slice);
    
    for kk1 = 1:size(thisImage,1)
        for kk2 = 1:size(thisImage,2)
            if thisImage(kk1,kk2) ~= 0
                thisPatch = thisImage(max(1,kk1-1):min(size(thisImage,1),kk1+1),max(1,kk2-1):min(size(thisImage,1),kk2+1));
                if numel(setdiff(unique(thisPatch),0)) > 1
                    currentImage(kk1, kk2) = 0;
                end
            end
        end
    end
    
    output(:,:,slice) = currentImage;
    fprintf( 'processed slice %d...\n', slice );
end

imwrite(output(:,:,1), outputFileName,'Compression','none');
for kk=2:numSlices
    imwrite(output(:,:,kk),outputFileName,'WriteMode','append','Compression','none');
end

matlabpool close
