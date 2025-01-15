% FNC state occurrences corr with task stimuli, group level, EFT

%% 1. Face mean SPM_design 
clc;clear;close all
load('data\R_FNC_states_corr_task_EFT.mat')
% mean_SPM_face_*(nTR*1)     % time-course of task stimuli derived from SPM.mat, averaged across all subjects 



%% 2.Face FNC state occurrences
% dfncInfo                   % parameter file of time-varying FNC 
% clusterInfo                % results of k-means clustering
% aIND                       % cluster labels, Nsub*Nwin
% aIND_tmp                   % cluster labels, nTR*Nsub

Nsub = length (dfncInfo.outputFiles);		 % number of subjects
Nwin = length(clusterInfo.IDXall) / Nsub;	 % number of sliding windows
wsize = dfncInfo.wsize;                      % window length
aIND = reshape(clusterInfo.IDXall, Nsub, Nwin); 
aIND_tmp=[nan(size(aIND,1),ceil(wsize/2)),aIND,nan(size(aIND,1),floor(wsize/2))];
aIND_tmp = aIND_tmp';                                
Nk = dfncInfo.postprocess.num_clusters;		 % number of FNC states


occurrence_face = zeros(Nwin+wsize,Nk);

for ii = 1:dfncInfo.postprocess.num_clusters % line 657
    
    %%%% mean occurrence %%%%
    oc = 100*mean(aIND_tmp == ii,2); oc([1:wsize/2, (end-wsize/2+1):end])=nan;
    occurrence_face(:,ii) = oc;
    
end

occurrence_face = occurrence_face(:,reorder_idx); % FNC state sequence match with REST 



%% 3.correlation between stimuli and oc, group level
r_oc = nan(Nk,4);        % Nk*Ncondition
p_oc = nan(Nk,4); 
ci_oc_lower = nan(Nk,4); % 95% CI of corr
ci_oc_upper = nan(Nk,4); % 95% CI of corr
for ii = 1:Nk
    oc = occurrence_face(:,ii);

    %%%% add 95%CI %%%%
    [R,P,RL,RU] = corrcoef([oc,mean_SPM_face_angry, mean_SPM_face_neutral, mean_SPM_face_happy, mean_SPM_face_control],'rows','pairwise');
    r_oc(ii,:) = R(1,2:end);   
    p_oc(ii,:) = P(1,2:end);    
    ci_oc_lower(ii,:) = RL(1,2:end); 
    ci_oc_upper(ii,:) = RU(1,2:end); 
    clear R P RL RU
end

save('data\R_FNC_states_corr_task_EFT.mat')


%% 4. Plot correlation between stimuli and oc
clc;close all;

cmap = [[0,68,27];[90,174,97];[118,42,131];[70,70,70]]/255;
linewidth = 1.5;
cmap_oc = [0.7,0.7,0.7];
close all;
for k = 1:Nk
    figure; ylim([-0.4,1.4])
    oc = occurrence_face(:,k);

    subplot(411);plot(oc,'Color',cmap_oc,'LineWidth',linewidth+0.5);yyaxis right; plot(mean_SPM_face_angry,'-.','Color',cmap(1,:),'LineWidth',linewidth);set(gca,'xtick',[]);box off
    subplot(412);plot(oc,'Color',cmap_oc,'LineWidth',linewidth+0.5);yyaxis right; plot(mean_SPM_face_neutral,'-.','Color',cmap(2,:),'LineWidth',linewidth);set(gca,'xtick',[]);box off
    subplot(413);plot(oc,'Color',cmap_oc,'LineWidth',linewidth+0.5);yyaxis right; plot(mean_SPM_face_happy,'-.','Color',cmap(3,:),'LineWidth',linewidth);set(gca,'xtick',[]);box off
    subplot(414);plot(oc,'Color',cmap_oc,'LineWidth',linewidth+0.5);yyaxis right; plot(mean_SPM_face_control,'-.','Color',cmap(4,:),'LineWidth',linewidth);box off

    saveas(figure(k),[output,'r_oc_task_group_reorder_state',num2str(k)],'tif')
end


%%%%% Plot for paper, barplot of correlation values %%%%%
x =[1,2,3,4];
figure; h=bar(x,r_oc,'EdgeColor','none'); 
h(1).FaceColor=cmap(1,:);h(2).FaceColor=cmap(2,:);
h(3).FaceColor=cmap(3,:);h(4).FaceColor=cmap(4,:); xticklabels([]); box off; ylim([-1,1])
saveas(gcf,[output,'r_oc_task_group_reorder_barplotR',num2str(k)],'tif')
