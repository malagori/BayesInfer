seed = 2;
rand('state', seed);

tic;

p = 1;

n = 10000;
combined_statistics = zeros(n, 2);

for i = 1:n
    [bnet, statistics] = generate_hidden_var_bn(p);
    combined_statistics(i, :) = statistics;
end

subplot(2, 1, 1);
hist(combined_statistics(:, 1));
subplot(2, 1, 2);
hist(combined_statistics(:, 2));

toc;