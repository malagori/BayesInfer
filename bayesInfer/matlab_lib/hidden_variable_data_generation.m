function [ok, data_file, result_file] = hidden_variable_data_generation(n, p, seed, result_directory, data_directory, pathMatlabLib)
% n = number of samples
% p = parameter controlling conditional distributions, If p << 1, this encourages "deterministic" CPTs (one entry near 1, the rest near 0). If p = 1, each entry is drawn from U[0,1]. If p >> 1, the entries will all be near 1/k, where k is the arity of this node, i.e., each row will be nearly uniform. 
% seed = seed for the random nnumber generator
% result_directory = where results are stored
% data_directory = where data is stored

%addpath('toolboxes');
addpath(genpathKPM(pathMatlabLib))

ok = 0;

% Initialize random number generator
rng(seed, 'twister');

% generate a BN and learn statistics
[bnet, statistics] = generate_hidden_var_bn(p);

% sample from the BN
samples = cell(5, n);
for i=1:n
  samples(:,i) = sample_bnet(bnet);
end

data = cell2mat(samples)';

% Use constraint-based methods to find the hidden variable
data_new = data(:, 1:4);
% pdag_tmp = learn_struct_pdag_ic_star('cond_indep_chisquare', 4, 3, data_new');
% hidden_found_fci = pdag_tmp(2, 4) == 2;
% 
% pdag_with_hidden = learn_struct_pdag_ic_star('cond_indep_chisquare', 5, 4, data');

% Store results
file_suffix = [int2str(n) '_' num2str(p) '_' int2str(seed)];
result_file = [result_directory '/statistics_' file_suffix '.mat'];
%save(result_file, 'statistics', 'hidden_found_fci', 'pdag_tmp', 'pdag_with_hidden');
save(result_file, 'statistics');

data_file = [data_directory '/data_' file_suffix '.txt'];
%[rows, ia, ic] = unique(data_new, 'rows', 'sorted');
[rows, ia, ic] = unique(data_new, 'rows');

counts = zeros(length(rows(:, 1)), 1);
for i = 1:length(rows(:, 1))
    counts(i) = sum(ic == i);
end

output1 = [rows counts];
fid = fopen(data_file, 'w');
fprintf(fid, '%s\n', 'A,B,C,D,Counts');
for i = 1:length(output1(:, 1))
    fprintf(fid, '%d,%d,%d,%d,%d\n', output1(i, :));
end

fclose(fid);

ok = 1;
