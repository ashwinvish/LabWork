function plotboth(q,q2,unit)
% plotboth(q,q2,'unit'); 
% plot the 3D coordinates stored in the matrices q and q2 next to each other 
% unit is a string telling the units of the plots
% use the designation au for plotting z-axis as planes with x-y coordinates 
% properly adjusted. 
close all;
plot3(q2(:,1),q2(:,2),q2(:,3),'b.')
hold on;
plot3(q(:,1),q(:,2),q(:,3),'r.')

xs = ['caudal-rostral axis (' unit ')' ];
ys = ['medial-lateral axis (' unit ')'];
zs = ['dorsal-ventral axis'] ;

if strcmp(unit,'au')
   zs = [zs '(plane)'];
end 

zlabel(zs)
xlabel(xs)
ylabel(ys)
