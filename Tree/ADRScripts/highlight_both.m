function highlight_both(q,q2,d,d2,th,unit)
% plotboth(q,q2,'unit'); 
% plot the 3D coordinates stored in the matrices q and q2 next to each other 
% unit is a string telling the units of the plots
% use the designation au for plotting z-axis as planes with x-y coordinates 
% properly adjusted. 
%close all;
plot3(q(:,1),q(:,2),-q(:,3),'go','MarkerSize',4) ; hold on;
plot3(q2(:,1),q2(:,2),-q2(:,3),'ro','MarkerSize',4)
% plot approximate soma locations
plot3(q(1,1),q(1,2),-q(1,3),'go','MarkerSize',15)
plot3(q2(1,1),q2(1,2),-q2(1,3),'ro','MarkerSize',15)

% plot overlapping points 
plot3(q(d<=th,1),q(d<=th,2),-q(d<=th,3),'y.')
plot3(q2(d2<=th,1),q2(d2<=th,2),-q2(d2<=th,3),'y.')



ys = ['caudal-rostral axis (' unit ')' ];
xs = ['medial-lateral axis (' unit ')'];
zs = ['dorsal-ventral axis'] ;

if strcmp(unit,'au')
   zs = [zs '(plane)'];
end 

zlabel(zs)
xlabel(xs)
ylabel(ys)
