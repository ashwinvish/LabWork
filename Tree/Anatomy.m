% Plot Anatomy of the fish
hold on;

Commisure1 = load('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/SWCFiles/8kTraces/Commisure1.mat');
Commisure1 = Commisure1.Commisure1;
shp1 = alphaShape(5*Commisure1(:,1),5*Commisure1(:,2),-45*Commisure1(:,3));
h1 = plot(shp1);
lightangle(-155,30);
h1.FaceColor = [0.8 0.8 0.8];
h1.EdgeColor = 'none';
h1.FaceLighting = 'gouraud';
h1.AmbientStrength = 0.3;
h1.DiffuseStrength = 0.1;
h1.SpecularStrength = 0.9;
h1.SpecularExponent = 25;
h1.BackFaceLighting = 'unlit';

Commisure2 = load('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/SWCFiles/8kTraces/Commisure2.mat');
Commisure2 = Commisure2.Commisure2;
shp2 = alphaShape(5*Commisure2(:,1),5*Commisure2(:,2),-45*Commisure2(:,3));
h2 = plot(shp2);
lightangle(-155,30);
h2.FaceColor = [0.8 0.8 0.8];
h2.EdgeColor = 'none';
h2.FaceLighting = 'gouraud';
h2.AmbientStrength = 0.3;
h2.DiffuseStrength = 0.1;
h2.SpecularStrength = 0.4;
h2.SpecularExponent = 25;
h2.BackFaceLighting = 'unlit';

Commisure3 = load('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/SWCFiles/8kTraces/Commisure3.mat');
Commisure3 = Commisure3.Commisure3;
shp3 = alphaShape(5*Commisure3(:,1),5*Commisure3(:,2),-45*Commisure3(:,3));
h3 = plot(shp3);
lightangle(-155,30);
h3.FaceColor = [0.8 0.8 0.8];
h3.EdgeColor = 'none';
h3.FaceLighting = 'gouraud';
h3.AmbientStrength = 0.3;
h3.DiffuseStrength = 0.1;
h3.SpecularStrength = 0.4;
h3.SpecularExponent = 25;
h3.BackFaceLighting = 'unlit';

Commisure4 = load('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/SWCFiles/8kTraces/Commisure4.mat');
Commisure4 = Commisure4.Commisure4;
shp4 = alphaShape(5*Commisure4(:,1),5*Commisure4(:,2),-45*Commisure4(:,3));
h4 = plot(shp4);
lightangle(-155,30);
h4.FaceColor = [0.8 0.8 0.8];
h4.EdgeColor = 'none';
h4.FaceLighting = 'gouraud';
h4.AmbientStrength = 0.3;
h4.DiffuseStrength = 0.1;
h4.SpecularStrength = 0.4;
h4.SpecularExponent = 25;
h4.BackFaceLighting = 'unlit';

Commisure5 = load('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/SWCFiles/8kTraces/Commisure5.mat');
Commisure5 = Commisure5.Commisure5;
shp5 = alphaShape(5*Commisure5(:,1),5*Commisure5(:,2),-45*Commisure5(:,3));
h5 = plot(shp5);
lightangle(-155,30);
h5.FaceColor = [0.8 0.8 0.8];
h5.EdgeColor = 'none';
h5.FaceLighting = 'gouraud';
h5.AmbientStrength = 0.3;
h5.DiffuseStrength = 0.1;
h5.SpecularStrength = 0.9;
h5.SpecularExponent = 25;
h5.BackFaceLighting = 'unlit';

Commisure6 = load('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/SWCFiles/8kTraces/Commisure6.mat');
Commisure6 = Commisure6.Commisure6;
shp6 = alphaShape(5*Commisure6(:,1),5*Commisure6(:,2),-45*Commisure6(:,3));
h6 = plot(shp6);
lightangle(-155,30);
h6.FaceColor = [0.8 0.8 0.8];
h6.EdgeColor = 'none';
h6.FaceLighting = 'gouraud';
h6.AmbientStrength = 0.3;
h6.DiffuseStrength = 0.1;
h6.SpecularStrength = 0.9;
h6.SpecularExponent = 25;
h6.BackFaceLighting = 'unlit';

Commisure7 = load('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/SWCFiles/8kTraces/Commisure7.mat');
Commisure7 = Commisure7.Commisure7;
shp7 = alphaShape(5*Commisure7(:,1),5*Commisure7(:,2),-45*Commisure7(:,3));
h7 = plot(shp7);
lightangle(-155,30);
h7.FaceColor = [0.8 0.8 0.8];
h7.EdgeColor = 'none';
h7.FaceLighting = 'gouraud';
h7.AmbientStrength = 0.3;
h7.DiffuseStrength = 0.1;
h7.SpecularStrength = 0.9;
h7.SpecularExponent = 25;
h7.BackFaceLighting = 'unlit';

MauthnerAxon = load('/usr/people/ashwinv/seungmount/research/Ashwin/MIT/Emre_HindBrain/SWCFiles/8kTraces/MauthnerAxon.mat');
MauthnerAxon = MauthnerAxon.MauthnerAxon;
AnatomicalSpline(5*MauthnerAxon(:,1),5*MauthnerAxon(:,2),-45*MauthnerAxon(:,3), [0.9,0.9,0.9]);



% shp8 = alphaShape(5*MauthnerAxon(:,1),5*MauthnerAxon(:,2),-45*MauthnerAxon(:,3));
% h8 = plot(shp8);
% lightangle(-155,30);
% h8.FaceColor = [0.8 0.8 0.8];
% h8.EdgeColor = 'none';
% h8.FaceLighting = 'gouraud';
% h8.AmbientStrength = 0.3;
% h8.DiffuseStrength = 0.8;
% h8.SpecularStrength = 0.4;
% h8.SpecularExponent = 25;
% h.BackFaceLighting = 'unlit';

box on;
axis([ 20000 140000 60000 250000 -60000 0]);
plot( [20000, 40000], [70000, 70000],'-k' ) % insert 20um sclaebar
daspect([1 1 1]); % make aspect ratio [1 1 1]
%set (gca,'Ydir','reverse');
set (gca,'XTick',[], 'YTick',[],'ZTick', [], 'Ydir','reverse');
view([-180,90]); % xy view



