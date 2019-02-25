load('/Users/ashwin/Documents/LabWork/SynapseDetector/f48870437u97212/dual_single_xfp_final_with_STA_112014.mat')

ind_vpni = find(logical(FIN.vpni)); % indices for staf, stad

rsthr = 0.9;
T = -2:0.05:7;

ind_alx = find(logical(FIN.stripe(1,:) & FIN.rsqs > rsthr));
ind_vglut_pos = find(logical(FIN.vglut==1 & FIN.rsqs > rsthr));
ind_vglut_neg = find(logical(FIN.vglut==-1 & FIN.rsqs > rsthr));
ind_dbx = find(logical(FIN.stripe(2,:) & FIN.rsqs > rsthr));
[~,alx_pop] = intersect(ind_vpni,ind_alx);
ind_alx_vglut = intersect(ind_alx,ind_vglut_pos);
[~,alx_pos_pop] = intersect(ind_vpni,ind_alx_vglut);
ind_alx_neg = intersect(ind_alx,ind_vglut_neg);
[~,alx_neg_pop] = intersect(ind_vpni,ind_alx_neg);
ind_dbx_vglut = intersect(ind_dbx,ind_vglut_pos);
[~,dbx_vglut_pop] = intersect(ind_vpni,ind_dbx_vglut);
ind_dbx_neg = intersect(ind_dbx,ind_vglut_neg);
[~,dbx_vglut_neg_pop] = intersect(ind_vpni,ind_dbx_neg);

subplot(2,2,1); plot(T,FIN.staf(:,alx_pos_pop),'r');  title('DsRed + alx'); ylabel('\DeltaF/F');
set(gca,'box','off','TickDir','out','FontName','MyriadPro-Regular','FontSize',12); set(gca,'ylim',[-5 50]); set(gca,'xlim',[-2 7]);
subplot(2,2,3); plot(T,FIN.staf(:,alx_neg_pop),'b'); title('DsRed - alx');
set(gca,'box','off','TickDir','out','FontName','MyriadPro-Regular','FontSize',12); set(gca,'ylim',[-5 50]); set(gca,'xlim',[-2 7]);
subplot(2,2,2); plot(T,FIN.staf(:,dbx_vglut_pop,'m'); title('DsRed + dbx1b');
set(gca,'box','off','TickDir','out','FontName','MyriadPro-Regular','FontSize',12); set(gca,'ylim',[-5 50]); set(gca,'xlim',[-2 7]);
subplot(2,2,4); plot(T,FIN.staf(:,dbx_vglut_neg_pop),'g');  title('DsRed - dbx1b');
set(gca,'box','off','TickDir','out','FontName','MyriadPro-Regular','FontSize',12); set(gca,'ylim',[-5 50]); set(gca,'xlim',[-2 7]);
xlabel('time (sec)'); 

figure;
hold on;
shadedErrorBar(T,mean(FIN.staf(:,alx_pos_pop)/100,2), std(FIN.staf(:,alx_pos_pop)/100,[],2),'lineprops',{'r'});
shadedErrorBar(T,mean(FIN.staf(:,alx_neg_pop)/100,2), std(FIN.staf(:,alx_neg_pop)/100,[],2),'lineprops',{'b'});
figure;
shadedErrorBar(T,mean(FIN.staf(:,dbx_vglut_pop)/100,2), std(FIN.staf(:,dbx_vglut_pop)/100,[],2),'lineprops',{'m'});
hold on;
shadedErrorBar(T,mean(FIN.staf(:,dbx_vglut_neg_pop)/100,2), std(FIN.staf(:,dbx_vglut_neg_pop)/100,[],2),'lineprops',{'g'})

