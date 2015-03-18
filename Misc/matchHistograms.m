% Match histogram of source with that one of target, optionally it saves
% the result on a file.
function [adjusted] = matchHistograms(sourceFileName, targetFileName, outputFileName)

saveFile = false;

if exist('outputFileName','var')
    saveFile = true;
end
    

% Read source stack
sourceInfo = imfinfo( sourceFileName );
for i=1:length(sourceInfo)
    ims(:,:,i) = imread( sourceFileName, i);
end

% Read target stack
targetInfo = imfinfo( targetFileName );
for i=1:length( targetInfo )
    imt(:,:,i) = imread( targetFileName, i);
end

ht = hist( double(imt(:)), 256 );

adjusted = histeq( ims(:), ht );

adjusted = reshape( adjusted, size( ims ) );

if saveFile
    imwrite(adjusted(:,:,1), outputFileName,'Compression','none');
    for kk=2:length( sourceInfo )
        imwrite(adjusted(:,:,kk),outputFileName,'WriteMode','append','Compression','none');
    end
end

end