
%% 1.Corr of each conn with task stimuli, example EFT           
clc;clear;close all;
dFNC_mat = 'R5_dFNC_face/';
load(fullfile(dFNC_mat,'k4/FU2_face_ICA_dfnc.mat'))
load(fullfile(dFNC_mat,'/FU2_face_ICA_dfnc_sub_001_sess_001_results.mat'), 'FNCdyn');
 
Nwin = size(FNCdyn,1);                       % number of sliding windows
Nconn =size(FNCdyn,2);  clear FNCdyn         % number of connections = Ncomp*(Ncomp-1)/2
Nsub = length (dfncInfo.outputFiles);		 % number of subjects
Ncomp = length(dfncInfo.comps);              % number of ICA components
wsize = dfncInfo.wsize;                      % window length

%%%%%%%%%%%% Face SPM design %%%%%%%%%%%%
r_SPM_conn_corr_face = zeros(Nsub,Nconn,4);
p_SPM_conn_corr_face = zeros(Nsub,Nconn,4);

for sub=1:Nsub
    sub
    % Load TC, FU2_face_ICA_dfnc_sub_1160_sess_001_results.mat
    if sub<999
        sub_str = sprintf('%03d',sub);
    else
        sub_str = sprintf('%d',sub);
    end
    
    load(fullfile(dFNC_mat,['FU2_face_ICA_dfnc_sub_',sub_str,'_sess_001_results.mat']), 'FNCdyn')
    
    tmp_SPM= SPM_xX_face_incl{sub};
    tmp_SPM= tmp_SPM((wsize/2+1): (Nwin+wsize/2), :);  % num of TR to win
    
    tmp_SPM_angry = sum(tmp_SPM(:,1:4),2);
    tmp_SPM_neutral = sum(tmp_SPM(:,5:8),2);
    tmp_SPM_happy = sum(tmp_SPM(:,9:12),2);
    tmp_SPM_control = tmp_SPM(:,13);
    
    for i=1:Nconn        
        tmp_conn = FNCdyn(:,i);
        [r_SPM_conn_corr_face(sub,i,1),p_SPM_conn_corr_face(sub,i,1)]=corr(tmp_conn, tmp_SPM_angry);
        [r_SPM_conn_corr_face(sub,i,2),p_SPM_conn_corr_face(sub,i,2)]=corr(tmp_conn, tmp_SPM_happy);
        [r_SPM_conn_corr_face(sub,i,3),p_SPM_conn_corr_face(sub,i,3)]=corr(tmp_conn, tmp_SPM_neutral);
        [r_SPM_conn_corr_face(sub,i,4),p_SPM_conn_corr_face(sub,i,4)]=corr(tmp_conn, tmp_SPM_control);
    end
end


% demo data first 100 subjects
demo_sub = 100;
r_SPM_conn_corr_face = r_SPM_conn_corr_face(1:demo_sub,:,:);
r_SPM_conn_corr_mid = r_SPM_conn_corr_mid(1:demo_sub,:,:);
r_SPM_conn_corr_sst = r_SPM_conn_corr_sst(1:demo_sub,:,:);

% save('data\R_corr_SPM_regional_FNC.mat', 'r_SPM_conn_corr_*')

%% 2. Plot top Nfc that distinguish task conditions in different tasks
clc; clear; close all;
load('data\R_corr_SPM_regional_FNC.mat')
% r_SPM_conn_corr_*      % Corr of each conn with task stimuli,Nsub*Nconn*Nevents

[~,p_SPM_conn_corr_face,~,stats]=ttest(r_SPM_conn_corr_face(:,:,1),r_SPM_conn_corr_face(:,:,3));% angry vs neutral
t_SPM_conn_corr_face = stats.tstat;
[~, t_idx_face] = sort((t_SPM_conn_corr_face),'descend'); % positive contrast

[~,p_SPM_conn_corr_mid,~,stats]=ttest(r_SPM_conn_corr_mid(:,:,1),r_SPM_conn_corr_mid(:,:,3));   % largewin vs nowin
t_SPM_conn_corr_mid = stats.tstat;
[~, t_idx_mid] = sort((t_SPM_conn_corr_mid),'descend');   

[~,p_SPM_conn_corr_sst,~,stats]=ttest(r_SPM_conn_corr_sst(:,:,1),r_SPM_conn_corr_sst(:,:,2));   % stopsucc vs stopfail
t_SPM_conn_corr_sst = stats.tstat;
[~, t_idx_sst] = sort((t_SPM_conn_corr_sst),'descend');  


Nfc = [1:1:20];

cmap = [[0,197,0]; [251,72,74]; [4,101,240]]/255;

for i = 1:length(Nfc)
    % t stats of top Nfc conn indexed from EFT
    t_face_conn_face(i) = t_SPM_conn_corr_face(t_idx_face(Nfc(i)));
    t_mid_conn_face(i) = t_SPM_conn_corr_mid(t_idx_face(Nfc(i)));
    t_sst_conn_face(i) = t_SPM_conn_corr_sst(t_idx_face(Nfc(i)));
    
    % t stats of top Nfc conn indexed from MID
    t_face_conn_mid(i) = t_SPM_conn_corr_face(t_idx_mid(Nfc(i)));
    t_mid_conn_mid(i) = t_SPM_conn_corr_mid(t_idx_mid(Nfc(i)));    
    t_sst_conn_mid(i) = t_SPM_conn_corr_sst(t_idx_mid(Nfc(i)));    
    
    % t stats of top Nfc conn indexed from SST
    t_face_conn_sst(i)= t_SPM_conn_corr_face(t_idx_sst(Nfc(i)));
    t_mid_conn_sst(i) = t_SPM_conn_corr_mid(t_idx_sst(Nfc(i)));
    t_sst_conn_sst(i) = t_SPM_conn_corr_sst(t_idx_sst(Nfc(i)));

end

clc;close all;
x0=50; y0=0;
width=800; height=1000;

figure; 
subplot(311);plot(t_face_conn_face,'.-','Color',cmap(1,:),'MarkerSize',24,'LineWidth',1.5);hold on; %title('EFT specific FNC');
plot(t_mid_conn_face,'.-','Color',cmap(2,:),'MarkerSize',24,'LineWidth',1.5);
plot(t_sst_conn_face,'.-','Color',cmap(3,:),'MarkerSize',24,'LineWidth',1.5);

FNC_all = t_SPM_conn_corr_face;
xshade = [0.5,20.5,20.5,0.5];
yshade = [prctile(FNC_all,25), prctile(FNC_all,25), prctile(FNC_all,75), prctile(FNC_all,75)];
fill(xshade, yshade, 'k', 'FaceAlpha', 0.1,'EdgeColor','none')
line([0.5,20.5],[prctile(FNC_all,50), prctile(FNC_all,50)],'Color',[0.2,0.2,0.2],'LineWidth',1.5,'LineStyle','--'); hold off;
ylim([-12 15]); xlim([0 21]); box off


subplot(312);plot(t_face_conn_mid,'.-','Color',cmap(1,:),'MarkerSize',24,'LineWidth',1.5);hold on; %title('MID specific FNC');
plot(t_mid_conn_mid,'.-','Color',cmap(2,:),'MarkerSize',24,'LineWidth',1.5);
plot(t_sst_conn_mid,'.-','Color',cmap(3,:),'MarkerSize',24,'LineWidth',1.5);

FNC_all = t_SPM_conn_corr_mid;
xshade = [0.5,20.5,20.5,0.5];
yshade = [prctile(FNC_all,25), prctile(FNC_all,25), prctile(FNC_all,75), prctile(FNC_all,75)];
fill(xshade, yshade, 'k', 'FaceAlpha', 0.1,'EdgeColor','none')
line([0.5,20.5],[prctile(FNC_all,50), prctile(FNC_all,50)],'Color',[0.2,0.2,0.2],'LineWidth',1.5,'LineStyle','--'); hold off;
ylim([-12 15]); xlim([0 21]); box off



subplot(313);plot(t_face_conn_sst,'.-','Color',cmap(1,:),'MarkerSize',24,'LineWidth',1.5);hold on; %title('SST specific FNC');
plot(t_mid_conn_sst,'.-','Color',cmap(2,:),'MarkerSize',24,'LineWidth',1.5);
plot(t_sst_conn_sst,'.-','Color',cmap(3,:),'MarkerSize',24,'LineWidth',1.5);

FNC_all = t_SPM_conn_corr_sst;
xshade = [0.5,20.5,20.5,0.5];
yshade = [prctile(FNC_all,25), prctile(FNC_all,25), prctile(FNC_all,75), prctile(FNC_all,75)];
fill(xshade, yshade, 'k', 'FaceAlpha', 0.1,'EdgeColor','none')
line([0.5,20.5],[prctile(FNC_all,50), prctile(FNC_all,50)],'Color',[0.2,0.2,0.2],'LineWidth',1.5,'LineStyle','--'); hold off;
ylim([-12 15]); xlim([0 21]); box off


set(gcf,'position',[x0,y0,width,height])



