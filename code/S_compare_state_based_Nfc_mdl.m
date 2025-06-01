% matlab: Statistics and Machine Learning Toolbox

%% 1. fitlme model_all, model_dFNC, model_sFNC compare, average across Nfc
clc;clear;close all;
load('data\R_compare_state_based_Nfc_mdl.mat')
% dFNC                 % time-varying FNC connectivity, Nsub*Nconn*Nk
% staticFNC            % static FNC connectivity,       Nsub*Nconn
% dFNC_reg             % time-varying FNC averaged for top Nfc task specific connectivity, regressed covariates, Nsub*Nk
% sFNC_reg             % static FNC averaged for top Nfc task specific connectivity, regressed covariates, Nsub*1
% Beh_reg              % 29 behavioural items, regressed covariates
% t_idx                % task-specificity connectivity index
% Nfc                  % number of top task-specificity connectivity

Nfc = 20;

%%%%% regress covariates from dFNC & sFNC %%%%%%
conn_idx = t_idx(1:Nfc);
Nk=size(dFNC,3);

dFNC_reg = squeeze(nanmean(dFNC(:,conn_idx,:),2));      % Nsub*1
dFNC_reg = my_out_zscore_covar(dFNC_reg,[],[],covar,'off');

sFNC_reg = nanmean(staticFNC(:,conn_idx),2);
sFNC_reg = my_out_zscore_covar(sFNC_reg,[],[],covar,'off');



mdl_fitting_results = cell(size(Beh_reg,2), 1);
mdl_compare_results = cell(size(Beh_reg,2), 1);

for i=1:size(Beh_reg,2)

    tmp_beh = Beh_reg(:,i);
    tbl =array2table([dFNC_reg, sFNC_reg, tmp_beh], 'VariableNames',{'state1','state2','state3','state4','sFNC','beh'});

    %%%% model estimation %%%%
    mdl_all_beh29 =  fitlme(tbl, 'beh~state1+state2+state3+state4+sFNC');
    mdl_dFNC_beh29 = fitlme(tbl, 'beh~state1+state2+state3+state4');
    mdl_sFNC_beh29 = fitlme(tbl, 'beh~sFNC');

    % mdl_fitting_results, rows: mdl_all, mdl_dFNC, mdl_sFNC; 
    % cols:R2, R2adj, MSE, numObs, numVariables, AIC, BIC, LogLikelihood, Deviance
    tmp_results = array2table([[mdl_all_beh29.Rsquared.Ordinary,mdl_all_beh29.Rsquared.Adjusted, mdl_all_beh29.MSE,...
        mdl_all_beh29.NumObservations,  mdl_all_beh29.NumVariables,  double(mdl_all_beh29.ModelCriterion)];...
        ...
        [mdl_dFNC_beh29.Rsquared.Ordinary,mdl_dFNC_beh29.Rsquared.Adjusted, mdl_dFNC_beh29.MSE,...
        mdl_dFNC_beh29.NumObservations,  mdl_dFNC_beh29.NumVariables,  double(mdl_dFNC_beh29.ModelCriterion)];...

        [mdl_sFNC_beh29.Rsquared.Ordinary,mdl_sFNC_beh29.Rsquared.Adjusted, mdl_sFNC_beh29.MSE,...
        mdl_sFNC_beh29.NumObservations,  mdl_sFNC_beh29.NumVariables,  double(mdl_sFNC_beh29.ModelCriterion)]],...
         'RowNames',{'mdl_all', 'mdl_dFNC', 'mdl_sFNC'}, 'VariableNames',{'R2', 'R2adj', 'MSE', 'numObs', 'numVar', 'AIC', 'BIC', 'LogLikelihood', 'Deviance'});
    mdl_fitting_results{i} = tmp_results;  clear tmp_results 

    %%%% model comparison %%%%
    % mdl_compare_results, rows: sFNC_all, dFNC_all;          
    % cols:LRStat, deltaDF, pValue
    mdl_compare_sFNC_all = compare(mdl_sFNC_beh29, mdl_all_beh29);
    mdl_compare_dFNC_all = compare(mdl_dFNC_beh29, mdl_all_beh29);

    tmp_results = array2table([[mdl_compare_sFNC_all.LRStat(2),mdl_compare_sFNC_all.deltaDF(2),mdl_compare_sFNC_all.pValue(2)];...
        [mdl_compare_dFNC_all.LRStat(2),mdl_compare_dFNC_all.deltaDF(2),mdl_compare_dFNC_all.pValue(2)]],...
        'RowNames',{'sFNC_all', 'dFNC_all'}, 'VariableNames',{'LRStat', 'deltaDF', 'pValue'});
    mdl_compare_results{i} = tmp_results; clear tmp_results 

end




