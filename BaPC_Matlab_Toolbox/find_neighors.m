function [Out,Id] = find_neighors(G, ys, N)

n = length(ys);

L = (G.laplacian);

Out = zeros(N,n);
for i = 1: n
    
   [~,idx] = sort(L(:,i),'descend');
   Out(:,i) = ys(idx(1:N));
   if i==1
      Id = idx(1:N);
   end
   
end
