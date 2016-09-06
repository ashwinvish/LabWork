% To calculate the overlapping area between cells

ksize = 12000;
res = 1000;

AlxPreVol = HeatMapFish(ksize,res, AlxPre,AlxSoma,cellIDsAlx,false);
AlxPostVol = HeatMapFish(ksize,res, AlxPost,AlxSoma,cellIDsAlx,false);
TransPreVol = HeatMapFish(ksize,res, TransPre,TransSoma,cellIDsTrans,false);
TransPostVol = HeatMapFish(ksize,res, TransPost,TransSoma,cellIDsTrans,false);
DbxPostVol = HeatMapFish(ksize,res, DbxPost,DbxSoma,cellIDsDbx,false);
BarhlPostVol = HeatMapFish(ksize,res, BarhlPost,BarhlSoma,cellIDsL,false);

%Ipsi axons with Ipsi dendrites
figure();
AlxArea = dotVol(AlxPreVol,AlxPostVol,AlxSoma, AlxSoma, res);
suptitle('Ipis Axons with Ipsi dendrites');


% IpsiContra axons with IpsiContraDendrites
figure();
TransArea = dotVol(TransPreVol,TransPostVol,TransSoma, TransSoma, res);
suptitle('IpisContra Axons with IpsiContra dendrites');

% IpsiAxons with Contra dendrites
figure();
subplot(1,2,1);
AlxDbxArea = dotVol(AlxPreVol,DbxPostVol,AlxSoma, DbxSoma, res);
subplot(1,2,2);
TransDbxArea = dotVol(TransPreVol,DbxPostVol,TransSoma, DbxSoma, res);
suptitle('Ipis Axons with Contra dendrites');


% IpsiAxons with Unknown dendrites
figure();
subplot(1,2,1);
AlxBarhlArea = dotVol(AlxPreVol,BarhlPostVol,AlxSoma, BarhlSoma, res);
subplot(1,2,2);
TransBarhlArea = dotVol(TransPreVol,BarhlPostVol,TransSoma, BarhlSoma, res);
suptitle('Ipis Axons with Unkown dendrites');


figure();
bar([ AlxArea/1e6,TransArea/1e6,(AlxDbxArea+TransDbxArea)/1e6, (AlxBarhlArea+TransBarhlArea)/1e6] , 'Marker','o','MarkerFaceColor', calx,  'MarkerEdgeColor','k');




