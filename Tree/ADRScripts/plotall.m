% script that makes a figure for 
% every integrator completed so far 
% adr 
% 4/24/14
% 

% show work completed so far
fnames = findswcfiles;
fh=figure;
% as seperate files 
for j=1:length(fnames)
  q=showskeleton(fnames{j}); 
  figure(fh);
  plot3(q(:,1),q(:,2),q(:,3),'b.');hold on;
  title(fnames{j})
end
%title(fnames{j})

if 0 
% together in the same file
 showskeleton(fnames{j},'k');
 for j=2:length(fnames)
     showskeleton(fnames{j},'k',0);
 end
end