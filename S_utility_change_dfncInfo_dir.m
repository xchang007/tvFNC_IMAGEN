function dfncInfo_new = S_utility_change_dfncInfo_dir(dfncInfo, newpath, changecompFiles, R4_folder)
% change R5_dFNC_rest_nocov_8TR\FU2_rest_ICA_dfnc.mat  dfncInfo folder info
%
% example:
% load('R5_dFNC_rest_nocov_8TR\FU2_rest_ICA_dfnc.mat');
% dfncInfo = S_utility_change_dfncInfo_dir(dfncInfo, 'G:\Projects\IMAGEN\FU2_tbfMRI\Scripts');
% dfncInfo = S_utility_change_dfncInfo_dir(dfncInfo, 'G:\Projects\IMAGEN\FU2_tbfMRI\Scripts', 1, 'G:\Projects\IMAGEN\FU2_tbfMRI\Scripts\R4_rest_ICA');

% if component files need to be changed
if (~exist('changecompFiles', 'var'))
    changecompFiles = 0;
end

dfncInfo_new = dfncInfo;

tmp_path = strsplit(dfncInfo.outputDir,'Scripts');
dfncInfo_new.outputDir = [newpath,tmp_path{end}];
tmp_path = strsplit(dfncInfo.userInput.outputDir,'Scripts');
dfncInfo_new.userInput.outputDir = [newpath,tmp_path{end}];
tmp_path = strsplit(dfncInfo.userInput.ica_param_file,'Scripts');
dfncInfo_new.userInput.ica_param_file = [newpath,tmp_path{end}];


%%%%% if component files need to be changed %%%%%
if changecompFiles
    compFiles_names = dfncInfo.userInput.compFiles;    % *_mean_component_ica_s_all_.nii,*
    
    for i=1:size(compFiles_names,1)
        [~, name, ext] = fileparts(compFiles_names(i,:));
        compFiles_names_new(i,:) = fullfile(R4_folder,[name, ext]);
    end
    dfncInfo_new.userInput.compFiles = compFiles_names_new;
end

end