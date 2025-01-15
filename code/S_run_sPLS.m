%% Input data, dwell time
clc;clear; close all
load('data\R_run_sPLS.mat')
addpath(genpath('...\cca_pls_toolkit-master\cca_pls_toolkit-master'))

% mean_dwell_time  % dwell_time of FNC states from k-means clustering, example resting-state, Nsub*Nk
% Beh              % 29 behavioural items, Nsub*Nbeh
% covar            % covariates, Nsub*10 (age,sex,site,meanFD)
% folder_dir       % folder directory

%%%%%% create data folder structure %%%%%%
if ~exist([folder_dir,'state1'])
    mkdir([folder_dir,'state1\']); mkdir([folder_dir,'state1\data']);
    mkdir([folder_dir,'state2\']); mkdir([folder_dir,'state2\data']);
    mkdir([folder_dir,'state3\']); mkdir([folder_dir,'state3\data']);
    mkdir([folder_dir,'state4\']); mkdir([folder_dir,'state4\data']);
else
    if ~exist([folder_dir,'state1\data'])
        mkdir([folder_dir,'state1\data']);
        mkdir([folder_dir,'state2\data']);
        mkdir([folder_dir,'state3\data']);
        mkdir([folder_dir,'state4\data']);
    end
end


%%%%%% write data to data_dir %%%%%%

for k = 1:size(mean_dwell_time,2)

    X = mean_dwell_time(:,k);
    Y = Beh;
    C = covar;
    sub_excl = (~any(isnan(X),2)) & (~any(isnan(Y),2));
    X = X(sub_excl,:);
    Y = Y(sub_excl,:);
    C = C(sub_excl,:);
    
    data_dir = [folder_dir,'state',num2str(k),'\data'];
    save(fullfile(data_dir, 'X.mat'), 'X');
    save(fullfile(data_dir, 'Y.mat'), 'Y');
    save(fullfile(data_dir, 'C.mat'), 'C');

end

%% Analysis, loop FNC states 

for k =1:4
    k
    clear cfg res
    data_dir = [folder_dir,'state',num2str(k),'\data'];
    
    set_path;
    % Project folder
    cfg.dir.project = fullfile(fileparts(data_dir));
    % cfg.machine.name = 'rcca';
    cfg.machine.name = 'spls';
    
    % Framework settings
    cfg.frwork.name = 'holdout';
    cfg.frwork.split.nout = 3;
    cfg.machine.metric = {'trcorrel' 'correl' 'simwx' 'simwy' ...
        'trexvarx' 'trexvary'};
    cfg.machine.param.crit = 'correl';
    cfg.machine.simw = 'correlation-Pearson';

    % Environment settings
    cfg.env.comp = 'local';
    % Statistical inference settings
    cfg.stat.nperm = 5000;
    
    % Update cfg with defaults
    cfg = cfg_defaults(cfg);
    
    % Run analysis
    main(cfg);
    
    % Clean up analysis files to save disc space
    cleanup_files(cfg);
end