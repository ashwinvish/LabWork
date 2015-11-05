function [f1,f2,id1] =find_samecell_swcfiles(d1,d2)
% find the names of all files in two user inputed directories,d1,d2
% that match the pattern words number something or nothging number [words]words.swc
f10 = findswcfiles(d1);
f20 = findswcfiles(d2);

% this uses a hack to find the numbers.  it knows that
% one of the files will have a Int
N1 = length(f10);
N2 = length(f20);
if N1~=N2
    error('program needs the same number of swc files in both directories');
end


plne = zeros(N1,2);
cell = zeros(N1,2);
for j=1:N1
    plne(j,1) = str2num(f10{j}(4));
    plne(j,2) = str2num(f20{j}(4));
    
    if isempty(regexp(f10{j}(5),'\d'))
        % not a number
        cell(j,1) = str2num(f10{j}(6));
        cell(j,2) = str2num(f20{j}(5));
    else
        cell(j,1) = str2num(f10{j}(5));
        cell(j,2) = str2num(f20{j}(6));
    end
end

% first we sort the first index 
[id1,psind]=sort(10*plne(:,1) +cell(:,1));
for j=1:N1
    % sort first index
    f1{j} = [d1 f10{psind(j)}];
    for k=1:N1
        % look for k that matches sorted first index
        if id1(j) == 10*plne(k,2) +cell(k,2)
            f2{j} = [d2 f20{k}];
            break
        end
    end
end

% show sorted files
for j=1:N1
    [f1{j} ',' f2{j}]
end
