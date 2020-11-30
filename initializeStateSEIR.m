function [y0] = initializeStateSEIR(param)

%% initialize state
n = param.n;
ts = param.ts;
y0(:,:,1) = repmat([1,0,0,0],n,1); % corresponding to 4 dimensional probability [S E i r]
% y0 = repmat([1,0,0,0],n,1,ts); % corresponding to 4 dimensional probability [S E i r]

seeds = randperm(n,param.ny0);
y0(seeds,:,1) = repmat([0,0,1,0],param.ny0,1); % initiallize with few sick individuals
