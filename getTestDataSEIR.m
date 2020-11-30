function [m, inds, P] = getTestDataSEIR(testType, G, param, inds) % early onset antigen 
% IGM early onset antigen, IGG recovered antigen 

n = param.n;
m = sparse(n,1);
L = (G.laplacian); % should remain sparse, but AD cannot handle sparse matrices... 

recovered = param.recovered;
if lower(testType) == 'm' %(M for IGM 'G' for IGG 
    nTests = round(rand(1)*n/100)+1;% number of tests in the time step
else
    nTests = round(rand(1)*n/50)+1;% number of tests in the time step
end
nonRecovered = find(~recovered); % identify those who are not known to be recovered
nonrecov = sparse(ones(n,1) - recovered);
if isempty(inds)
    inds = nonRecovered(randperm(nnz(~recovered),nTests)); % define those who got tested
end
nTests = length(inds);
    %temp = abs(L*nonrecov);
m(inds) = round(rand(nTests,1)); % assign either 0 or 1 to the tested individuals
%m(inds) = temp(inds)>max(temp)*rand;
% build a projector of the tested
if isempty(inds)
    P = 0;
else
    P = sparse(1:nTests,inds,ones(nTests,1),nTests,n);
end
% update recovered list outside the function