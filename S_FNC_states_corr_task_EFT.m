% FNC state occurrences corr with task stimuli, group level, EFT
% S10_corr_stimuli_oc_face.m 1-3.1

%% 1. Face mean SPM_design 

clc;clear;close all
load('R3_subject_outlier.mat','tf_*')
load('R12_SPM_design.mat', '*_face')

incl = find(tf_QC_task.*tf_SPM_face.*tf_FD_face.*tf_face.*tf_vol_face.*tf_demo_task==1);% 1263 
clear tf_*

SPM_tbl_face_incl = SPM_tbl_face(incl);
SPM_xX_face_incl = SPM_xX_face(incl);


%%%%%%%%%%%% 2.plot mean SPM_design %%%%%%%%%%%%
nTR = min(cellfun(@length, SPM_xX_face_incl));  % min TR 192

mean_SPM_face_angry = zeros(nTR,1);
mean_SPM_face_neutral = zeros(nTR,1);
mean_SPM_face_happy = zeros(nTR,1);
mean_SPM_face_control = zeros(nTR,1);

for sub=1:size(SPM_xX_face_incl,1)
    tmp_xX = SPM_xX_face_incl{sub};
    mean_SPM_face_angry = mean_SPM_face_angry+sum(tmp_xX(1:nTR,1:4),2); 
    mean_SPM_face_neutral = mean_SPM_face_neutral+sum(tmp_xX(1:nTR,5:8),2); 
    mean_SPM_face_happy = mean_SPM_face_happy+sum(tmp_xX(1:nTR,9:12),2); 
    mean_SPM_face_control = mean_SPM_face_control+tmp_xX(1:nTR,end);
end

mean_SPM_face_angry = mean_SPM_face_angry/size(SPM_xX_face_incl,1);
mean_SPM_face_neutral = mean_SPM_face_neutral/size(SPM_xX_face_incl,1);
mean_SPM_face_happy = mean_SPM_face_happy/size(SPM_xX_face_incl,1);
mean_SPM_face_control = mean_SPM_face_control/size(SPM_xX_face_incl,1);

% save('R10_corr_stimuli_oc_face.mat','mean_SPM*')


%% 2.Face FNC state occurrences

clc;clear;close all;
cd 'G:\Projects\IMAGEN\FU2_tbfMRI\Scripts';
load('R10_corr_stimuli_oc_face.mat','mean_SPM_face*')

%%%%% get occurrance %%%%%
cd('R5_dFNC_face_new_8TR_reg_all_con/k4')
load('FU2_face_ICA_dfnc.mat')
load('FU2_face_ICA_dfnc_post_process.mat','clusterInfo')

cd 'G:\Projects\IMAGEN\FU2_tbfMRI\Scripts';

Nsub = length (dfncInfo.outputFiles);		 % 1263
Nwin = length(clusterInfo.IDXall) / Nsub;	 % 184 = 197-5-wsize
wsize = dfncInfo.wsize;                      % 
aIND = reshape(clusterInfo.IDXall, Nsub, Nwin); % 1263*184
aIND_tmp=[nan(size(aIND,1),ceil(wsize/2)),aIND,nan(size(aIND,1),floor(wsize/2))];
aIND_tmp = aIND_tmp';                                
Nk = dfncInfo.postprocess.num_clusters;


occurrence_face = zeros(Nwin+wsize,Nk);

for ii = 1:dfncInfo.postprocess.num_clusters % line 657
    
    %%%% mean occurrence %%%%
    oc = 100*mean(aIND_tmp == ii,2); oc([1:wsize/2, (end-wsize/2+1):end])=nan;
    occurrence_face(:,ii) = oc;
    
end


% save('R10_corr_stimuli_oc_face_reg_all_con.mat','occurrence_face','dfncInfo','N*','wsize','aIND_tmp','-append')


%% 3.correlation between stimuli and oc, group level
clc;clear;close all;

load('R10_corr_stimuli_oc_face_reg_all_con.mat')
wsize = dfncInfo.wsize;                       
Nk=dfncInfo.postprocess.num_clusters;

for ii = 1:Nk
    oc = occurrence_face(:,ii);
    
    [r_oc(ii,1),p_oc(ii,1)]=corr(oc,mean_SPM_face_angry,'rows','complete');
    [r_oc(ii,2),p_oc(ii,2)]=corr(oc,mean_SPM_face_happy,'rows','complete');
    [r_oc(ii,3),p_oc(ii,3)]=corr(oc,mean_SPM_face_neutral,'rows','complete');
    [r_oc(ii,4),p_oc(ii,4)]=corr(oc,mean_SPM_face_control,'rows','complete');
    
end

% save('R10_corr_stimuli_oc_face_reg_all_con.mat','r_oc','p_oc','-append')




%% 3-1 Plot for paper, r_oc, save to /R10_plot_corr_stimuli_oc/
clc;clear;close all;
load('R10_corr_stimuli_oc_face_reg_all_con.mat')
load('R_utility_reorder_index.mat')

output = 'R10_plot_corr_stimuli_oc/EFT_8TR_reg_all_con/';  % save to
reorder_idx = R5_dFNC_face_new_8TR_reg_all_con_k4_reorder_idx;  % R5_dFNC_face_new_8TR_reg_all_con
clear *_reorder_idx

% Re-ordered occurrence_face and r_oc
occurrence_face_reorder = occurrence_face(:,reorder_idx);
r_oc_reorder = r_oc(:,[1,3,2,4]);           % Note:since 2rd col of r_oc is happy, 3rd col is neutral, swap to match above plot
r_oc_reorder = r_oc_reorder(reorder_idx,:); % Reorder FC states
p_oc_reorder = p_oc(:,[1,3,2,4]);           % Note:since 2rd col of r_oc is happy, 3rd col is neutral, swap to match above plot
p_oc_reorder = p_oc_reorder(reorder_idx,:); % Reorder FC states

% save('R10_corr_stimuli_oc_face_reg_all_events.mat','*_reorder','-append')
% save([output,'r_oc_task_individual_reorder_compare.mat'],'*reorder*','-append')


cmap = [[0,68,27];[90,174,97];[118,42,131];[70,70,70]]/255;
linewidth = 1.5;
cmap_oc = [0.7,0.7,0.7];
close all;
for k = 1:Nk
    figure; ylim([-0.4,1.4])
    oc = occurrence_face_reorder(:,k);

    subplot(411);plot(oc,'Color',cmap_oc,'LineWidth',linewidth+0.5);yyaxis right; plot(mean_SPM_face_angry,'-.','Color',cmap(1,:),'LineWidth',linewidth);set(gca,'xtick',[]);box off
    subplot(412);plot(oc,'Color',cmap_oc,'LineWidth',linewidth+0.5);yyaxis right; plot(mean_SPM_face_neutral,'-.','Color',cmap(2,:),'LineWidth',linewidth);set(gca,'xtick',[]);box off
    subplot(413);plot(oc,'Color',cmap_oc,'LineWidth',linewidth+0.5);yyaxis right; plot(mean_SPM_face_happy,'-.','Color',cmap(3,:),'LineWidth',linewidth);set(gca,'xtick',[]);box off
    subplot(414);plot(oc,'Color',cmap_oc,'LineWidth',linewidth+0.5);yyaxis right; plot(mean_SPM_face_control,'-.','Color',cmap(4,:),'LineWidth',linewidth);box off

    saveas(figure(k),[output,'r_oc_task_group_reorder_state',num2str(k)],'tif')
end


%%%%% Plot for paper, barplot of correlation values %%%%%
x =[1,2,3,4];
figure; h=bar(x,r_oc_reorder,'EdgeColor','none'); 
h(1).FaceColor=cmap(1,:);h(2).FaceColor=cmap(2,:);
h(3).FaceColor=cmap(3,:);h(4).FaceColor=cmap(4,:); xticklabels([]); box off; ylim([-1,1])
saveas(gcf,[output,'r_oc_task_group_reorder_barplotR',num2str(k)],'tif')
