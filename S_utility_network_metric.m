function [net_global, net_node] = S_utility_network_metric(FC,gamma, Nrandom)

if ~exist('gamma','var') 
    gamma = 1;
end

if ~exist('Nrandom','var') 
    Nrandom = 0;
end

% global strength, both pos and negative
[Spos_node,Sneg_node,Spos_global,Sneg_global] = strengths_und_sign(FC);

% community, participation, within-module degree z-score, clustering
[M,Q] = community_louvain(FC,gamma,[],'negative_sym');
[Ppos_node,Pneg_node]=participation_coef_sign(FC, M);
Ppos_global = mean(Ppos_node);  
Pneg_global = mean(Pneg_node);
Within_degree_node = module_degree_zscore(FC, M);
Within_degree_global = mean(Within_degree_node);

[Cpos_node,Cneg_node,Cpos_global,Cneg_global] = clustering_coef_wu_sign(FC);

% efficiency, only based on positive conn
FC_pos = FC.*(FC>0);
FC_pos = FC_pos./(Spos_global/nnz(FC>0));  % divided by global mean pos
Epos_global = efficiency_wei(FC_pos);
Epos_node = efficiency_wei(FC_pos,2);

% compare with randomized network, only calculate global measurement
if Nrandom ~= 0
    for n=1:Nrandom
        FC_rand=randmio_und_signed(FC,10);
        [Mrand,Qrand(n)] = community_louvain(FC_rand,gamma,[],'negative_sym');
        [Ppos_node_rand,Pneg_node_rand]=participation_coef_sign(FC_rand, Mrand);
        Ppos_rand(n) = mean(Ppos_node_rand);
        Pneg_rand(n) = mean(Pneg_node_rand);
        Within_degree_node_rand=module_degree_zscore(FC_rand, Mrand);
        Within_degree_rand(n)=mean(Within_degree_node_rand);
        
        [~,~,Cpos_rand(n),Cneg_rand(n)] = clustering_coef_wu_sign(FC_rand);
        
        % efficiency, only based on positive conn
        FC_rand_pos = FC_rand.*(FC_rand>0);
        FC_rand_pos = FC_rand_pos./(sum(sum(FC_rand_pos))/nnz(FC_rand_pos));  % divided by global mean pos
        Epos_rand(n) = efficiency_wei(FC_rand_pos);
    end
    Q = Q/mean(Qrand);
    Ppos_global = Ppos_global/mean(Ppos_rand);
    Pneg_global = Pneg_global/mean(Pneg_rand);
    Within_degree_global = Within_degree_global/mean(Within_degree_rand);
    Cpos_global = Cpos_global/mean(Cpos_rand);
    Cneg_global = Cneg_global/mean(Cneg_rand);
    Epos_global = Epos_global/mean(Epos_rand);
end

net_global.Q = Q; 
net_global.M = M; 
net_global.Spos_global = Spos_global; 
net_global.Sneg_global = Sneg_global; 
net_global.Ppos_global = Ppos_global; 
net_global.Pneg_global = Pneg_global; 
net_global.Within_degree_global = Within_degree_global; 
net_global.Epos_global = Epos_global; 
net_global.Cpos_global = Cpos_global; 
net_global.Cneg_global = Cneg_global; 

net_node.Spos_node = Spos_node'; 
net_node.Sneg_node = Sneg_node'; 
net_node.Ppos_node = Ppos_node; 
net_node.Pneg_node = Pneg_node; 
net_node.Within_degree_node = Within_degree_node; 
net_node.Cpos_node = Cpos_node; 
net_node.Cneg_node = Cneg_node; 
net_node.Epos_node = Epos_node; 

end 