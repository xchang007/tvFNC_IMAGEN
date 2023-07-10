%% Average FC at each time point (FNCdyn_all, group) 

%% S23_conn_corr_rest.m
clc;clear;close all
cd G:\Projects\IMAGEN\FU2_tbfMRI\Scripts
load('R3_subject_outlier.mat','tf_*')

load('R5_dFNC_rest_filter_8TR/k4/FU2_rest_ICA_dfnc.mat','dfncInfo')
dFNC_mat = 'R5_dFNC_rest_filter_8TR';
load(fullfile(dFNC_mat,dfncInfo.outputFiles{1}))
incl = find(tf_FD_rest==1);% 991 
clear tf_*

Nwin = size(FNCdyn,1); 
Nconn =size(FNCdyn,2);  clear FNCdyn
Nsub = length (dfncInfo.outputFiles);		 % 991
Ncomp = length(dfncInfo.comps);
wsize = dfncInfo.wsize;                      % 8

FNCdyn_all = zeros(Nwin, Nconn);
for sub=1:Nsub
    sub
    % Load TC, FU2_face_ICA_dfnc_sub_1160_sess_001_results.mat
    if sub<999
        sub_str = sprintf('%03d',sub);
    else
        sub_str = sprintf('%d',sub);
    end    
    %%%%%%% change here according to sessions %%%%%%%
    load(fullfile(dFNC_mat,['FU2_rest_ICA_dfnc_sub_',sub_str,'_sess_001_results.mat']), 'FNCdyn')

    FNCdyn_all = FNCdyn_all + FNCdyn; clear FNCdyn
end

FNCdyn_all = FNCdyn_all/Nsub;

r_FCMatrix=corr(FNCdyn_all');
figure; imagesc(r_FCMatrix); caxis([0.98,1]); colorbar; 

% save('R23_conn_corr_rest_filter_8TR.mat', 'FNCdyn_all','N*','wsize','dFNC_mat')


%% 2.S23_conn_corr_face.m  %3
clc;clear;close all
load('R3_subject_outlier.mat','tf_*')
load('R5_dFNC_face_new_8TR_reg_all_con/k4/FU2_face_ICA_dfnc.mat','dfncInfo')
dFNC_mat = 'R5_dFNC_face_new_8TR_reg_all_con';

load(fullfile(dFNC_mat,dfncInfo.outputFiles{1}))

load('R10_corr_stimuli_oc_face.mat', 'mean_SPM*')
incl = find(tf_QC_task.*tf_SPM_face.*tf_FD_face.*tf_face.*tf_vol_face.*tf_demo_task==1);% 1263 
clear tf_*


Nwin = size(FNCdyn,1); 
Nconn =size(FNCdyn,2);  clear FNCdyn
Nsub = length (dfncInfo.outputFiles);		 % 1263
Ncomp = length(dfncInfo.comps);
wsize = dfncInfo.wsize;                      % 8


FNCdyn_all = zeros(Nwin, Nconn);
for sub=1:Nsub
    sub
    % Load TC, FU2_face_ICA_dfnc_sub_1160_sess_001_results.mat
    if sub<999
        sub_str = sprintf('%03d',sub);
    else
        sub_str = sprintf('%d',sub);
    end    
    %%%%%%% change here according to sessions %%%%%%%
    load(fullfile(dFNC_mat,['FU2_face_ICA_dfnc_sub_',sub_str,'_sess_001_results.mat']), 'FNCdyn')

    FNCdyn_all = FNCdyn_all + FNCdyn; clear FNCdyn
end

FNCdyn_all = FNCdyn_all/Nsub;

r_FCMatrix=corr(FNCdyn_all');
figure; imagesc(r_FCMatrix); caxis([0.94,1]); colorbar; 

% save('R23_conn_corr_task_face_8TR_reg_all_con.mat', 'FNCdyn_all','*SPM_conn_corr_face_grp','-append')


%% 3.S23_conn_corr_mid.m  %3
clc;clear;close all
% dFNC_mat = 'R5_dFNC_mid_new_4TR_reg_all/';
load(fullfile(dFNC_mat,'k4/FU2_mid_ICA_dfnc.mat'))
load(fullfile(dFNC_mat,'/FU2_mid_ICA_dfnc_sub_001_sess_001_results.mat'), 'FNCdyn');

load(fullfile(dFNC_mat,dfncInfo.outputFiles{1}))

load('R3_subject_outlier.mat','tf_*')
incl = find(tf_QC_task.*tf_SPM_mid.*tf_FD_mid.*tf_mid.*tf_vol_mid.*tf_demo_task==1);% 1221 
clear tf_*

Nwin = size(FNCdyn,1); 
Nconn =size(FNCdyn,2);  clear FNCdyn		 % 178
Nsub = length (dfncInfo.outputFiles);		 % 1221
Ncomp = length(dfncInfo.comps);
wsize = dfncInfo.wsize;                      % 8


FNCdyn_all = zeros(Nwin, Nconn);
for sub=1:Nsub
    sub
    % Load TC, FU2_face_ICA_dfnc_sub_1160_sess_001_results.mat
    if sub<999
        sub_str = sprintf('%03d',sub);
    else
        sub_str = sprintf('%d',sub);
    end    
    load(fullfile(dFNC_mat,['FU2_mid_ICA_dfnc_sub_',sub_str,'_sess_001_results.mat']), 'FNCdyn')
    FNCdyn_all = FNCdyn_all + FNCdyn; clear FNCdyn
end

FNCdyn_all = FNCdyn_all/Nsub;

r_FCMatrix=corr(FNCdyn_all');
figure; imagesc(r_FCMatrix); caxis([0.98,1]); colorbar; 

% save('R23_conn_corr_task_mid_4TR_reg_all.mat', 'FNCdyn_all','*SPM_conn_corr_mid_grp','dFNC_mat','dfncInfo','-append')

%% 4.S23_conn_corr_sst.m  %3
clc;clear;close all
load('R3_subject_outlier.mat','tf_*')

load('R5_dFNC_sst_new_8TR_reg_stop/k4/FU2_sst_ICA_dfnc.mat','dfncInfo')
dFNC_mat = 'R5_dFNC_sst_new_8TR_reg_stop';

load(fullfile(dFNC_mat,dfncInfo.outputFiles{1}))

incl = find(tf_FD_sst.*tf_sst.*tf_vol_sst.*tf_demo_task==1);% 1292 
clear tf_*

Nwin = size(FNCdyn,1); 
Nconn =size(FNCdyn,2);  clear FNCdyn
Nsub = length (dfncInfo.outputFiles);		 % 1330
Ncomp = length(dfncInfo.comps);
wsize = dfncInfo.wsize;                      % 8


FNCdyn_all = zeros(Nwin, Nconn);
for sub=1:Nsub
    sub
    % Load TC, FU2_face_ICA_dfnc_sub_1160_sess_001_results.mat
    if sub<999
        sub_str = sprintf('%03d',sub);
    else
        sub_str = sprintf('%d',sub);
    end    
    %%%%%%% change here according to sessions %%%%%%%
    load(fullfile(dFNC_mat,['FU2_sst_ICA_dfnc_sub_',sub_str,'_sess_001_results.mat']), 'FNCdyn')

    FNCdyn_all = FNCdyn_all + FNCdyn; clear FNCdyn
end

FNCdyn_all = FNCdyn_all/Nsub;

r_FCMatrix=corr(FNCdyn_all');
figure; imagesc(r_FCMatrix); caxis([0.98,1]); colorbar; 

% save('R23_conn_corr_task_sst_new_8TR_reg_stop.mat', 'FNCdyn_all','N*','wsize','dFNC_mat')
