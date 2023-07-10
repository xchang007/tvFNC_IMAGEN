% FNC state occurrences corr with task stimuli, individual level, SST
% S10_corr_stimuli_oc_sst.m 6,2-5

%% 1.individual FC pattern corr with task design, partial corr  

clc;clear;close all
cd G:\Projects\IMAGEN\FU2_tbfMRI\Scripts
load('R10_corr_stimuli_oc_sst_8TR_reg_stop.mat')

load('R12_SPM_design.mat', '*_sst')
load('R3_subject_outlier.mat','tf_*');
incl = find(tf_QC_task.*tf_SPM_sst.*tf_FD_sst.*tf_sst.*tf_vol_sst.*tf_demo_task==1);  % 1218
clear tf_*
SPM_tbl_sst_incl = SPM_tbl_sst(incl);
SPM_xX_sst_incl = SPM_xX_sst(incl);

%%%%% 3. individual FC pattern corr with task design %%%%%
r_oc_task_individual_partialcorr = nan(Nsub,Nk,4); % ss,sf,gl,gw
p_oc_task_individual_partialcorr = nan(Nsub,Nk,4); % ss,sf,gl,gw

for sub=1:Nsub
    sub
    tmp_aIND = aIND_tmp(:,sub);
    tmp_SPM = SPM_xX_sst_incl{sub}(1:Nscan,:);
    % Note: subjects have 197/202 timepoints, only take 1:197 (remove first 5)
    
    for k=1:max(tmp_aIND)
        tmp_aIND_k = double(tmp_aIND==k); 
        tmp_aIND_k(isnan(tmp_aIND))=nan;  % if tmp_aIND doesn't contain k, corr will be NaN
        
        %%%%%% stop success %%%%%%
        [r,p]=partialcorr(tmp_aIND_k,tmp_SPM(:,1),tmp_SPM(:,2:4),'rows','complete');
        r_oc_task_individual_partialcorr(sub,k,1) = r;
        p_oc_task_individual_partialcorr(sub,k,1) = p;
        
        %%%%%% stop failure %%%%%%
        [r,p]=partialcorr(tmp_aIND_k,tmp_SPM(:,2),tmp_SPM(:,[1,3,4]),'rows','complete');
        r_oc_task_individual_partialcorr(sub,k,2) = r;
        p_oc_task_individual_partialcorr(sub,k,2) = p;
        
        %%%%%% go late %%%%%%
        [r,p]=partialcorr(tmp_aIND_k,tmp_SPM(:,3),tmp_SPM(:,[1,2,4]),'rows','complete');
        r_oc_task_individual_partialcorr(sub,k,3) = r;
        p_oc_task_individual_partialcorr(sub,k,3) = p;
        
        %%%%%% go wrong %%%%%%
        [r,p]=partialcorr(tmp_aIND_k,tmp_SPM(:,4),tmp_SPM(:,1:3),'rows','complete');
        r_oc_task_individual_partialcorr(sub,k,4) = r;
        p_oc_task_individual_partialcorr(sub,k,4) = p;
    end
    
end

% save('R10_corr_stimuli_oc_sst_8TR_reg_stop.mat','*_oc_task_individual_partialcorr','-append')


%% 2. Violinplot individual FC pattern corr with task design, reordered
clc;clear;close all;

load('R10_corr_stimuli_oc_sst_8TR_reg_stop.mat')
load('R_utility_reorder_index.mat')
reorder_idx = R5_dFNC_sst_new_8TR_reg_stop_k4_reorder_idx; clear R5_dFNC_*  % R5_dFNC_sst_new_8TR_reg_stop

% output = 'R10_plot_corr_stimuli_oc/SST_8TR_reg_stop/';  % save to
% r_oc_task_individual_reorder = r_oc_task_individual(:,reorder_idx,:); % reordered,Nsub,Nk,Ncondition % ss,sf,gl,gw

output = 'R10_plot_corr_stimuli_oc/SST_8TR_reg_stop_partialcorr/';  % save to
r_oc_task_individual_reorder = r_oc_task_individual_partialcorr(:,reorder_idx,:); % reordered,Nsub,Nk,Ncondition % ss,sf,gl,gw

cmap = [[8,81,156];[66,146,198];[199,233,180];[237,248,177]]/255;

close all;
for k=1:Nk
    figure;h=violinplot(squeeze(r_oc_task_individual_reorder(:,k,:)),[],'Width',0.3,  'ViolinAlpha', 0.1); 
    h(1).ViolinColor = cmap(1,:); h(2).ViolinColor = cmap(2,:); h(3).ViolinColor = cmap(3,:); h(4).ViolinColor = cmap(4,:);
    set(gca, 'ylim',[-0.7,0.7]);
    set (gcf,'Position', [100 100 400 400])
    
    if ~exist(output, 'dir')
        mkdir(output)
    end
    saveas(figure(k),[output,'r_oc_task_individual_reorder_state',num2str(k)],'tif')

end


%% 3.Compare r_oc_task_individual with zeros in each state
Ncondition = 4;  % ss,sf,gl,gw

p_ttest_oc_task_individual_zeros = zeros(Nk,Ncondition); % paired-t between each two of 4 conditions
t_ttest_oc_task_individual_zeros = zeros(Nk,Ncondition); % paired-t between each two of 4 conditions
for k=1:Nk
    for cond = 1:Ncondition
        [~,p,~,stats] = ttest(squeeze(r_oc_task_individual_reorder(:,k,cond)));
        p_ttest_oc_task_individual_zeros(k,cond) = p; 
        t_ttest_oc_task_individual_zeros(k,cond) = stats.tstat;
    end
end

[h_fdr,p_thres]=fdr_bh(p_ttest_oc_task_individual_zeros,0.05,'pdep','yes')  % 8TR_reg_stop p_thres=0.0187

p_ttest_oc_task_individual_zeros=array2table(p_ttest_oc_task_individual_zeros,...
    'VariableNames',{'StopSucc', 'StopFail', 'GoLate', 'GoWrong'},'RowNames',{'State 1,reordered','State 2','State 3','State 4'});
t_ttest_oc_task_individual_zeros=array2table(t_ttest_oc_task_individual_zeros,...
    'VariableNames',{'StopSucc', 'StopFail', 'GoLate', 'GoWrong'},'RowNames',{'State 1,reordered','State 2','State 3','State 4'});


%% 4.Compare r_oc_task_individual between conditions in each state

cond_paires = nchoosek(1:Ncondition,2); % all possible combinations of the elements of vector v taken k at a time

p_anova_oc_task_individual_condition= zeros(Nk,1);   % anova  between conditions 
F_anova_oc_task_individual_condition= zeros(Nk,1);   % anova  between conditions 

p_ttest_oc_task_individual_condition = zeros(Nk,size(cond_paires,1)); % paired-t between each two of 4 conditions
t_ttest_oc_task_individual_condition = zeros(Nk,size(cond_paires,1)); % paired-t between each two of 4 conditions

for k=1:Nk
    % anova  between conditions 
    [p,tbl] =anova1(squeeze(r_oc_task_individual_reorder(:,k,:)),[],'off');
    p_anova_oc_task_individual_condition(k) = p;
    F_anova_oc_task_individual_condition(k) = cell2mat(tbl(2,5));
    
    % post-hoc t-test between each pair of conditions
    
    for i = 1:size(cond_paires,1)  
        [~,p,~,stats] = ttest(squeeze(r_oc_task_individual_reorder(:,k,cond_paires(i,1))),squeeze(r_oc_task_individual_reorder(:,k,cond_paires(i,2))));
        p_ttest_oc_task_individual_condition(k,i) = p; 
        t_ttest_oc_task_individual_condition(k,i) = stats.tstat;
    end
end

p_anova_oc_task_individual_condition = array2table(p_anova_oc_task_individual_condition,'RowNames',{'State 1','State 2','State 3','State 4'});
F_anova_oc_task_individual_condition = array2table(F_anova_oc_task_individual_condition,'RowNames',{'State 1','State 2','State 3','State 4'});

p_ttest_oc_task_individual_condition = array2table(p_ttest_oc_task_individual_condition,...
    'VariableNames',cellstr(num2str(cond_paires)),'RowNames',{'State 1','State 2','State 3','State 4'});
t_ttest_oc_task_individual_condition = array2table(t_ttest_oc_task_individual_condition,...
    'VariableNames',cellstr(num2str(cond_paires)),'RowNames',{'State 1','State 2','State 3','State 4'});

save([output,'r_oc_task_individual_reorder_compare.mat'],'p_*','t_*','F_*','r_oc_task_individual*','reorder*','N*','cond_paires')


%% 5 Plot for paper, r_oc, example sub save to /R10_plot_corr_stimuli_oc/

%%%%% Plot for paper, task design, oc, example oc %%%%%
clc;clear;close all;
load('R3_subject_outlier.mat','tf_*')
load('R12_SPM_design.mat', '*_sst')
incl = find(tf_QC_task.*tf_SPM_sst.*tf_FD_sst.*tf_sst.*tf_vol_sst.*tf_demo_task==1);% 1263 
clear tf_*

% load('R10_corr_stimuli_oc_sst_8TR_reg_stop.mat')
% reorder_idx = [2,3,1,4];  % R5_dFNC_sst_new_8TR_reg_stop
load('R10_corr_stimuli_oc_sst_8TR.mat')
reorder_idx = [2,3,4,1];  % R5_dFNC_sst_new_8TR

output = 'R10_plot_corr_stimuli_oc/SST_8TR/';  % save to

SPM_xX_sst_incl = SPM_xX_sst(incl);

% example_sub oc
% tmp_r = squeeze(r_oc_task_individual(:,3,:)); % r_oc_individual, state 3, in case example_sub doesn't have state 3
example_sub = 152;  %98,110,156
oc_sub = aIND_tmp(:,example_sub);
r_oc_sub = squeeze(r_oc_task_individual(example_sub,reorder_idx,:));  % 4 State*4 Condition, not reordered

% example_sub SPM
SPM_sub = SPM_xX_sst_incl{example_sub}(1:(Nwin+wsize),:);

cmap = [[8,81,156];[66,146,198];[199,233,180];[237,248,177]]/255;
linewidth = 1.5;
cmap_oc = [0.7,0.7,0.7];
close all;
for k = 1:Nk
    figure; ylim([-0.4,1.4])
    oc = oc_sub==reorder_idx(k);
    subplot(411);plot(oc,'Color',cmap_oc,'LineWidth',linewidth+0.5);yyaxis right; plot(SPM_sub(:,1),'-.','Color',cmap(1,:),'LineWidth',linewidth);set(gca,'xtick',[]);box off
    subplot(412);plot(oc,'Color',cmap_oc,'LineWidth',linewidth+0.5);yyaxis right; plot(SPM_sub(:,2),'-.','Color',cmap(2,:),'LineWidth',linewidth);set(gca,'xtick',[]);box off
    subplot(413);plot(oc,'Color',cmap_oc,'LineWidth',linewidth+0.5);yyaxis right; plot(SPM_sub(:,3),'-.','Color',cmap(3,:),'LineWidth',linewidth);set(gca,'xtick',[]);box off
    subplot(414);plot(oc,'Color',cmap_oc,'LineWidth',linewidth+0.5);yyaxis right; plot(SPM_sub(:,4),'-.','Color',cmap(4,:),'LineWidth',linewidth);box off

    saveas(figure(k),[output,'r_oc_task_example_sub',num2str(example_sub),'_state',num2str(k)],'tif')
end


