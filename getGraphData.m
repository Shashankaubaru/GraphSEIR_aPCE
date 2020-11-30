function [G] = getGraphData(param, G_old)
% generate adjacency matrix
n = param.n;

if isempty(G_old)
    A = sprand(n,n,0.1);
    A = A+A';
    A = A-diag(sparse(diag(A)));
    G = graph(A);
else
    Go = G_old.laplacian;
   %G = G_old;
     Go(Go<0.01*rand) = 0;
     A = Go+sprand(n,n,0.01);
     A = A+A';
     A =  A-diag(sparse(diag(A)));
end

% if param.f1show
%     figure(param.f1), spy(A)
%     title('Adjacency matrix')
% end

 G = graph(A);
% D = G.degree; % degree matrix;
% L = G.laplacian; % graph laplacian
% 
% % compute graph laplacian
% LL = diag(D)-A;
% % In = incidence(G);

% edge weight information will arrive from actual BT data
%G.Edges.Weight = rand(size(G.Edges,1),1);


% if param.f2show
%     figure(param.f2); p=plot(G,'EdgeLabel',G.Edges.Weight);
%     p.Marker = 's';
%     p.NodeColor = 'r';
%     title('Interaction Graph (from BT data)')
% end
% % edges = G.Edges;
% nodes = G.Nodes;
% weights = rand(size(edges,1),1);
% G = graph(edges,nodes,weights);
