
%% Compare same and different FNC states
clc;clear;close all
load('data\R_FNC_states_network.mat')
% FC_pattern(Ncomponent*Ncomponent*Nk) % cluster centroids from time-varying FNC and k-means clustering
% Ncomp                                % 61 ICA components
% Nk                                   % number of FNC states


%%%%% all states %%%%%
FC_pattern_all_vec =[icatb_mat2vec(FC_pattern_rest_re(:,:,1)), icatb_mat2vec(FC_pattern_face_re(:,:,1)), icatb_mat2vec(FC_pattern_mid_re(:,:,1)), icatb_mat2vec(FC_pattern_sst_re(:,:,1)),...
    icatb_mat2vec(FC_pattern_rest_re(:,:,2)),icatb_mat2vec(FC_pattern_face_re(:,:,2)),icatb_mat2vec(FC_pattern_mid_re(:,:,2)),icatb_mat2vec(FC_pattern_sst_re(:,:,2)),...
    icatb_mat2vec(FC_pattern_rest_re(:,:,3)),icatb_mat2vec(FC_pattern_face_re(:,:,3)),icatb_mat2vec(FC_pattern_mid_re(:,:,3)),icatb_mat2vec(FC_pattern_sst_re(:,:,3)),...
    icatb_mat2vec(FC_pattern_rest_re(:,:,4)),icatb_mat2vec(FC_pattern_face_re(:,:,4)),icatb_mat2vec(FC_pattern_mid_re(:,:,4)),icatb_mat2vec(FC_pattern_sst_re(:,:,4))];

r = corr(FC_pattern_all_vec);
figure;imagesc(r,[0.4,1]);colorbar


%%%%% Statistic All states correlation %%%%%
tf_same_state = ones(size(r)); % FC matrices of different states =1
tf_same_state(1:4,1:4)=2;      % FC matrices of the same state =2
tf_same_state(5:8,5:8)=2;
tf_same_state(9:12,9:12)=2;
tf_same_state(13:16,13:16)=2;
tf_same_state(1:(size(r,1)+1):end)=0;
[h,p,ci,stats]=ttest2(r(triu(tf_same_state)==2),r(triu(tf_same_state)==1))
% p = 6.9296e-11, tstat: 7.1742,       df: 118   95CI: [0.1349 0.2378]
tmp = mes(r(triu(tf_same_state)==2),r(triu(tf_same_state)==1),'hedgesg');
cohend = tmp.hedgesg

r_within = r(triu(tf_same_state)==2);
r_between = r(triu(tf_same_state)==1);
[mean(r_within),mean(r_between)] % 0.8585    0.6722
[std(r_within),std(r_between)]   % 0.1658    0.0971


%% FNC network metrics
gamma = 1;
Nrand = 100;
for k=1:Nk
    k
    [net_global_rest(k), net_node_rest(k)] = S_utility_network_metric(FC_pattern_rest_re(:,:,k),gamma, Nrand);
    [net_global_face(k), net_node_face(k)] = S_utility_network_metric(FC_pattern_face_re(:,:,k),gamma, Nrand);
    [net_global_mid(k), net_node_mid(k)] = S_utility_network_metric(FC_pattern_mid_re(:,:,k),gamma, Nrand);
    [net_global_sst(k), net_node_sst(k)] = S_utility_network_metric(FC_pattern_sst_re(:,:,k),gamma, Nrand);
end


Ppos_all = [net_global_rest.Ppos_global; net_global_face.Ppos_global;...
    net_global_mid.Ppos_global; net_global_sst.Ppos_global];  % Nsession*Nk
Pneg_all = [net_global_rest.Pneg_global; net_global_face.Pneg_global;...
    net_global_mid.Pneg_global; net_global_sst.Pneg_global];  % Nsession*Nk
Q_all = [net_global_rest.Q; net_global_face.Q;...
    net_global_mid.Q; net_global_sst.Q]; % Nsession*Nk
Epos_all = [net_global_rest.Epos_global; net_global_face.Epos_global;...
    net_global_mid.Epos_global; net_global_sst.Epos_global];  % Nsession*Nk


[p,tbl] = anova1(Ppos_all); % p=3.5199e-05 ,F=22.1267
eta2=tbl{2,2}/tbl{4,2}      % SSeffect/SStotal
[p,tbl] = anova1(Pneg_all); % p=3.5843e-06 ,F=34.3838
eta2=tbl{2,2}/tbl{4,2} % SSeffect/SStotal
[p,tbl] = anova1(Q_all);    % p=2.1552e-06 ,F=37.8078
eta2=tbl{2,2}/tbl{4,2} % SSeffect/SStotal
[p,tbl] = anova1(Epos_all); % p=5.5857e-07 ,F=48.4374


% post-hoc t-test between each pair of conditions
Ngroup = 4;
cond_paires = nchoosek(1:Ngroup,2); % all possible combinations of the elements of vector v taken k at a time

for i = 1:size(cond_paires,1)
    [~,p,~,stats] = ttest(Ppos_all(:,cond_paires(i,1)),Ppos_all(:,cond_paires(i,2)));
    p_ttest_Ppos_all(i) = p;    t_ttest_Ppos_all(i) = stats.tstat;
    
    [~,p,~,stats] = ttest(Pneg_all(:,cond_paires(i,1)),Pneg_all(:,cond_paires(i,2)));
    p_ttest_Pneg_all(i) = p;    t_ttest_Pneg_all(i) = stats.tstat;
    
    [~,p,~,stats] = ttest(Q_all(:,cond_paires(i,1)),Q_all(:,cond_paires(i,2)));
    p_ttest_Q_all(i) = p;    t_ttest_Q_all(i) = stats.tstat;
    
    [~,p,~,stats] = ttest(Epos_all(:,cond_paires(i,1)),Epos_all(:,cond_paires(i,2)));
    p_ttest_Epos_all(i) = p;    t_ttest_Epos_all(i) = stats.tstat;
end

p_all = [p_ttest_Q_all;p_ttest_Ppos_all;p_ttest_Pneg_all];
[h, crit_p]=fdr_bh(p_all,0.05);         