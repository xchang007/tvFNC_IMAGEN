function FNC_module = S_utility_FNC_module(FNC, module_size, module_count)
% Plot module mean at lower triangle, FNC at upper triangle

if (~exist('module_size','var'))
     module_size = [5     6    13     8    17     5     7];
end

% module_count = 0, calculate non-zeros mean within module (default)
% module_count = 1, calculate non-zeros cells / total cells within module
if (~exist('module_count','var'))
     module_count = 0;
end

if isvector(FNC)
    FNC = icatb_vec2mat(FNC, 1);
end

node_index = cumsum(module_size);
node_index = [0,node_index];

FNC(1:(1+size(FNC,1)):end)=0;   % clear diagonal
FNC_upper = triu(FNC);          % Plot FNC at upper triangle

FNC_module_mean = zeros(size(FNC));

for module_i = 1:length(module_size)
    node_module_i = (node_index(module_i)+1):node_index(module_i+1);
    
    for module_j = module_i:length(module_size)
        node_module_j = (node_index(module_j)+1):node_index(module_j+1);
        if ~module_count % calculate non-zeros mean within module
            FNC_module_mean(node_module_i,node_module_j) = sum(sum(FNC_upper(node_module_i,node_module_j)))...
                ./nnz(FNC_upper(node_module_i,node_module_j));
        else             % calculate non-zeros cells / total cells within module
            FNC_module_mean(1:(size(FNC_module_mean,1)+1):end)=nan;   % make diagonal NaN, to exclude counting
            FNC_module_mean(node_module_i,node_module_j) = nnz(FNC_upper(node_module_i,node_module_j))...
                ./nnz(~isnan(FNC_upper(node_module_i,node_module_j)));
        end
    end
end
FNC_module_mean = triu(FNC_module_mean,1);
FNC_module = FNC_module_mean' + FNC_upper; % Plot FNC at upper triangle

end