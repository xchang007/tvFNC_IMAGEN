
% Plot re-ordered FC states from *ICA_dfnc_post_process.mat, 'clusterInfo'
% adopted from S14_pattern_modularity_distance.m %3

clc;clear;close all
addpath(genpath('E:\Matlab\Toolbox_Fudan\GroupICATv4.0b\'));

% cd /mnt/GroupUser/xchang/Projects/FU2_tbfMRI/Scripts/
% cd 'C:\Users\k1809663\OneDrive - King''s College London\Projects\IMAGEN\FU2_tbfMRI\Scripts'
cd G:\Projects\IMAGEN\FU2_tbfMRI\Scripts
load('R_utility_reorder_index.mat')

dfnc_folder = 'R5_dFNC_rest_filter_8TR';
state_num = 'k4';
reorder_idx = R5_dFNC_rest_filter_8TR_k4_reorder_idx; clear R5_dFNC_*


dfnc = strsplit(dfnc_folder,'_');
session = dfnc{3};
% TR = dfnc{end};
load(fullfile(dfnc_folder,state_num,['FU2_',session,'_ICA_dfnc.mat']), 'dfncInfo')
load(fullfile(dfnc_folder,state_num,['FU2_',session,'_ICA_dfnc_post_process.mat']), 'clusterInfo')
% load('R5_dFNC_rest_nocov_8TR/k4/FU2_rest_ICA_dfnc.mat','dfncInfo')
% load('R5_dFNC_rest_nocov_8TR/k4/FU2_rest_ICA_dfnc_post_process.mat','clusterInfo');

cmap = flipud(cbrewer('div' ,'RdBu',128)); cmap(cmap>1)=1; cmap(cmap<0)=0;
clims=[-0.4,0.4];
% [FC_pattern,FC_pattern_re] = plot_dfnc(dfncInfo, clusterInfo, reorder_idx, cmap, clims, 1, 1, 0);
[FC_pattern,FC_pattern_re] = plot_dfnc(dfncInfo, clusterInfo, reorder_idx, cmap, clims, 1, 1, 1);


function [FC_pattern,FC_pattern_re] =plot_dfnc(dfncInfo, clusterInfo, reorder_idx, cmap, clims, savepic, savemat, plotmodule)
close all

outfolder = dfncInfo.outputDir;
outfolder = strsplit(outfolder,{'/','\'},'CollapseDelimiters',true); %Split Character Vector with Multiple Delimiters
outfolder = outfolder{end};

if ~exist('savepic', 'var')
    savepic=0;
end

if ~exist('savemat', 'var')
    savemat=0;
end

if ~exist('plotmodule', 'var')
    plotmodule=0;  % Plot upper triangle as module mean
end
%%%%%% get FNC matrix from clusterInfo.Call %%%%%%
Nk= size(clusterInfo.Call,1);
Ncomp = 61;
FC_pattern = zeros(Ncomp,Ncomp,Nk);

% icatb_post_process_dfnc.m  line435
% Call is the clustering center of all timepoints 
% Cp is the center of representative timepoints
for k=1:Nk
    FC_pattern(:,:,k) = icatb_vec2mat(clusterInfo.Call(k,:));
end
FC_pattern_re = FC_pattern(:,:,reorder_idx);


%%%%%% plot FNC matrix %%%%%%
module_size = cellfun(@length,{dfncInfo.userInput.comp.value});
module_name = {dfncInfo.userInput.comp.name};

for k=1:Nk
    F = figure('Color','w');
    tmp = FC_pattern_re(:,:,k);
    
    if ~plotmodule
        h=icatb_plot_FNC(tmp,clims, [],1:61,F,'Connectivity',[],module_size,module_name);
    else  % Plot upper triangle as module mean
        tmp_module = S_utility_FNC_module(tmp,module_size);
        h=icatb_plot_FNC_new(tmp_module,clims, [],1:61,F,'Connectivity',[],module_size,module_name);
    end
    
    set(gca, 'YTick', [],'XTick', []); h.Colormap = cmap;
    set(gcf,'position',[100,100,750,650])
end


%%%%%% save picture? %%%%%%
if savepic
    if ~plotmodule
        if ~exist(['R_plot_FC_States/',outfolder], 'dir')
            mkdir(['R_plot_FC_States/',outfolder])
        end
        for k=1:Nk
            saveas(figure(k),['R_plot_FC_States/',outfolder,'/k',num2str(Nk),'_State',num2str(k)],'tif')
        end
        
    else
        if ~exist(['R_plot_FC_States_module/',outfolder], 'dir')
            mkdir(['R_plot_FC_States_module/',outfolder])
        end
        for k=1:Nk
            saveas(figure(k),['R_plot_FC_States_module/',outfolder,'/k',num2str(Nk),'_State',num2str(k)],'tif')
        end
    end
    
end

%%%%%% savemat? %%%%%%
if savemat
    if ~plotmodule
        save(['R_plot_FC_States/',outfolder,'/k',num2str(Nk),'.mat'],'FC_pattern*')
    else
        save(['R_plot_FC_States_module/',outfolder,'/k',num2str(Nk),'.mat'],'FC_pattern*')
    end
end

end
