vol = h5read('/usr/people/ashwinv/seungmount/research/Ashwin/Scripts/NG_scripts/77605.h5','/main');
numSections = size(vol,4);

for i =  1:1:numSections 
    imVol(:,:,i) = reshape(vol(1,:,:,i),[size(vol,2) , size(vol,3)]);
end


figure();
hold on;

for kk=1:size(edges,1)
 plot3([nodes(edges(kk,1),2); nodes(edges(kk,2),2)],[nodes(edges(kk,1),1); nodes(edges(kk,2),1)],[nodes(edges(kk,1),3); nodes(edges(kk,2),3)],'r');
end