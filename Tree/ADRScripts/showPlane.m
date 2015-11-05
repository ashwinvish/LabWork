function showPlane(q,q2,z)
%showPlane(q,q2,z,unit) - Show points in a given plane

unit = 'pixels';

% find indicies of points in plane z 
zind = q(:,3) == z; 
zind2 = q2(:,3) == z; 

figure
plot(q(zind,1),q(zind,2),'go','Markersize',8); hold on;
plot(q2(zind2,1),q2(zind2,2),'r.','Markersize',15); 

ys = ['caudal-rostral axis (' unit ')' ];
xs = ['medial-lateral axis (' unit ')'];

ylabel(ys)
xlabel(xs)
legend('ashwin','ashleigh');
axis ij

end

