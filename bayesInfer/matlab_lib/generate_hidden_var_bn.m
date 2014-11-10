function [bnet, statistics] = generate_hidden_var_bn(p)

statistics = zeros(2, 1);

N = 5;
dag = zeros(N, N);

A = 1;
B = 2;
C = 3;
D = 4;
H = 5;

dag(A, B) = 1;
dag(H, [B, C]) = 1;
dag(D, C) = 1;

discrete_nodes = 1:N;
node_sizes = 2*ones(1,N);
bnet = mk_bnet(dag, node_sizes, 'discrete', discrete_nodes);
CPTs = cell(N, 1);

for i=1:N
    k = node_sizes(i);
    ps = parents(dag, i);
    psz = prod(node_sizes(ps));
    CPT = dirichletSample(p*ones(1,k), psz);
    bnet.CPD{i} = tabular_CPD(bnet, i, 'CPT', CPT);
    CPTs{i} = CPT;
end

CPTs{B} = reshape(CPTs{B}, [2, 2, 2]);
CPTs{C} = reshape(CPTs{C}, [2, 2, 2]);

statistics(1) = CPTs{A}(1)*abs(CPTs{B}(1, 1, 1) - CPTs{B}(1, 1, 2) ) + CPTs{A}(2)*abs(CPTs{B}(2, 1, 1) - CPTs{B}(2, 1, 2) ) + CPTs{D}(1)*abs(CPTs{C}(1, 1, 1) - CPTs{C}(1, 1, 2) ) + CPTs{D}(2)*abs(CPTs{C}(2, 1, 1) - CPTs{C}(2, 1, 2) );
statistics(2) = CPTs{A}(1)*abs(CPTs{B}(1, 1, 1) - CPTs{B}(1, 1, 2) )^2 + CPTs{A}(2)*abs(CPTs{B}(2, 1, 1) - CPTs{B}(2, 1, 2) )^2 + CPTs{D}(1)*abs(CPTs{C}(1, 1, 1) - CPTs{C}(1, 1, 2) )^2 + CPTs{D}(2)*abs(CPTs{C}(2, 1, 1) - CPTs{C}(2, 1, 2) )^2;

