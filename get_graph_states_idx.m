function [Lout,yt] = get_graph_states_idx(G, ys, Id)

for t = 1:length(G)
    L = (G{t}.laplacian);
    Lout{t} = L(Id,Id);
end
yt = ys(Id,:,:);
