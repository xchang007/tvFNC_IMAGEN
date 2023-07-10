% pls_cca_toolkit\S_plot.m % 1-3
addpath(genpath('cca_pls_toolkit-master\cca_pls_toolkit-master'))

%% 1.Plot data projections for significant models
clc;clear; close all
cd 'G:\Projects\IMAGEN\FU2_tbfMRI\Scripts'

savefig = 'G:\Projects\IMAGEN\Manuscript\Figures\Figure_5_pls_cca_toolkit\spls_dwell_zeros_beh29_Nperm5000\';

p_thres = 0.05;
cmap=[[95,155,213];[255,127,0]]/255;
session = [{'R5_dFNC_rest_filter_8TR_zeros_beh29'},{'R5_dFNC_face_new_8TR_reg_all_con_zeros_beh29'},...
    {'R5_dFNC_mid_new_4TR_reg_all_zeros_beh29'},{'R5_dFNC_sst_new_8TR_reg_stop_zeros_beh29'}];

for k=1:4

    for ses = 1:length(session)
        % loop through all sessions
        res = load(['pls_cca_toolkit\',session{ses},'\state',num2str(k),'\framework\spls_holdout3-0.20_subsamp5-0.20\res\level1\res_1.mat']);
        res = res.res;
        model = load(fullfile(res.dir.res, 'model_1.mat'));
        [~,idx] = max(model.correl);        % split with highest correl (in test set)

        % Plot data projections if model pval<0.05
        if res.stat.pval(idx)< p_thres     
            plot_proj(res, {'X' 'Y'}, res.frwork.level, 'osplit',idx, 'training+test', ...
                '2d_group', 'gen.axes.FontSize', 15, ...
                'gen.legend.FontSize', 15, 'gen.legend.Location', 'best', ...
                'proj.scatter.SizeData', 45, 'proj.scatter.MarkerEdgeColor', 'none', ...
                'proj.xlabel', 'State dwell time latent variable',...
                'proj.scatter.MarkerFaceColor', cmap,'proj.lsline','on');

            saveas(gcf,[savefig,'state',num2str(k),'_',session{ses},'.tif']) % TIFF 24-bit (not compressed)
        end

    end
end


%% 2.1 Plot behaviour weights, dwell, heatmap
% S8_dFNC_Beh_pls_task_stability.m
clc;clear; close all
cd 'G:\Projects\IMAGEN\FU2_tbfMRI\Scripts'

Beh_label = readcell('S46_dFNC_Beh_cca_plot_Beh_label32.xlsx');
Beh_label = Beh_label(:,end);
Beh_label([3,5,10]) = []; 

savefig = 'G:\Projects\IMAGEN\Manuscript\Figures\Figure_5_pls_cca_toolkit\spls_dwell_zeros_beh29_Nperm5000\';

p_thres = 0.05;

for k=1:4
    
res_rest = load(['pls_cca_toolkit\R5_dFNC_rest_filter_8TR_zeros_beh29\state',num2str(k),'\framework\spls_holdout3-0.20_subsamp5-0.20\res\level1\res_1.mat']);
res_face = load(['pls_cca_toolkit\R5_dFNC_face_new_8TR_reg_all_con_zeros_beh29\state',num2str(k),'\framework\spls_holdout3-0.20_subsamp5-0.20\res\level1\res_1.mat']);
res_mid = load(['pls_cca_toolkit\R5_dFNC_mid_new_4TR_reg_all_zeros_beh29\state',num2str(k),'\framework\spls_holdout3-0.20_subsamp5-0.20\res\level1\res_1.mat']);
res_sst = load(['pls_cca_toolkit\R5_dFNC_sst_new_8TR_reg_stop_zeros_beh29\state',num2str(k),'\framework\spls_holdout3-0.20_subsamp5-0.20\res\level1\res_1.mat']);

model = load(fullfile(res_rest.res.dir.res, 'model_1.mat'));
[~,idx] = max(model.correl);        % split with highest correl (in test set)
beh_weight_rest = model.wY(idx,:)'; 
sig = res_rest.res.stat.pval(idx)< p_thres;

model = load(fullfile(res_face.res.dir.res, 'model_1.mat'));
[~,idx] = max(model.correl);        % split with highest correl
beh_weight_face = model.wY(idx,:)'; 
sig = [sig;res_face.res.stat.pval(idx)< p_thres];

model = load(fullfile(res_mid.res.dir.res, 'model_1.mat'));
[~,idx] = max(model.correl);        % split with highest correl
beh_weight_mid = model.wY(idx,:)'; 
sig = [sig;res_mid.res.stat.pval(idx)< p_thres];

model = load(fullfile(res_sst.res.dir.res, 'model_1.mat'));
[~,idx] = max(model.correl);        % split with highest correl
beh_weight_sst = model.wY(idx,:)'; 
sig = [sig;res_sst.res.stat.pval(idx)< p_thres];

clims=[-0.8,0.8];
cmap=flipud(cbrewer('div' ,'RdBu',128)); cmap(cmap>1)=1; cmap(cmap<0)=0;
x0=100;y0=100; width=800; height=300;

tmp_weight = [beh_weight_rest, beh_weight_face, beh_weight_mid, beh_weight_sst]';
% figure;plt=heatmap(tmp_weight, 'ColorLimits',clims,'CellLabelColor','none');
figure;plt=heatmap(tmp_weight.*sig, 'ColorLimits',clims,'CellLabelColor','none');
plt.Colormap=cmap; plt.XDisplayLabels=Beh_label;plt.YDisplayLabels={'Rest','EFT','MID','SST'};
set(gcf,'position',[x0,y0,width,height])
set(struct(plt).NodeChildren(3), 'XTickLabelRotation', 60);
title(['Behavioural weights'])
saveas(gcf,[savefig,'state',num2str(k),'.tif']) % TIFF 24-bit (not compressed)

end


%% 3. Load correl, pvals
clc;clear; close all
cd 'G:\Projects\IMAGEN\FU2_tbfMRI\Scripts'

% data_dir = 'G:\Projects\IMAGEN\FU2_tbfMRI\Scripts\pls_cca_toolkit\R5_dFNC_rest_filter_8TR_zeros_beh29\';
% data_dir = 'G:\Projects\IMAGEN\FU2_tbfMRI\Scripts\pls_cca_toolkit\R5_dFNC_face_new_8TR_reg_all_con_zeros_beh29\';
% data_dir = 'G:\Projects\IMAGEN\FU2_tbfMRI\Scripts\pls_cca_toolkit\R5_dFNC_mid_new_4TR_reg_all_zeros_beh29\';
data_dir = 'G:\Projects\IMAGEN\FU2_tbfMRI\Scripts\pls_cca_toolkit\R5_dFNC_sst_new_8TR_reg_stop_zeros_beh29\';

Nk=4;
pval = zeros(Nk,1);
correl = zeros(Nk,1);
trcorrel = zeros(Nk,1);
trexvary = zeros(Nk,1);

for k=1:Nk
    res = load([data_dir,'state',num2str(k),'\framework\spls_holdout3-0.20_subsamp5-0.20\res\level1\res_1.mat']);
    model = load([data_dir,'state',num2str(k),'\framework\spls_holdout3-0.20_subsamp5-0.20\res\level1\model_1.mat']);
    [~,idx]=max(model.correl);
    pval(k) = res.res.stat.pval(idx);
    correl(k) = model.correl(idx);
    trcorrel(k) = model.trcorrel(idx);
    trexvary(k) = model.trexvary(idx);
end

tmp=[trcorrel,correl,pval,trexvary]; open tmp
