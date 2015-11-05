function fnames =findswcfiles(varargin)
% find the names of all files in the current directory that 
% match the pattern
% words [words]words.swc
if isempty(varargin)
    fcont = dir(pwd);
else
      fcont = dir(varargin{1});
end
cnt = 1;    
for j=1:length(fcont)
    if ~fcont(j).isdir
        if ~isempty(regexp(fcont(j).name,'.swc$'))
            fnames{cnt} = fcont(j).name;
            cnt = cnt +1;
        end
    end
end