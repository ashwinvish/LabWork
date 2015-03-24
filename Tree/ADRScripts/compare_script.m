%% compare 2 skeletons with each other 

% get file names
f = uigetfile(pwd,'select first file'); 
f2 = uigetfile(pwd,'select second file'); 

% get skeletons
q = showskeleton(f);
q2 = showskeleton(f2);	

%% insert
plotboth(q,q2,'au')
legend('first','second');

