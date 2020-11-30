clc
close all;clear all;
dbstop if error
addpath(genpath('cvx'))
addpath(genpath('BaPC_Matlab_Toolbox'))
%% set paramters
param.n = 1000; % for now assume graph nodes are fixed, later nodal dynamics can be incorporated
param.ts = 100; % timesteps

param.N = 5; % No. of neighbours
param.d = 3; %Order of expansion 

param.recovered = zeros(param.n,1); % recovered individuals
param.ny0 = 5; % number of initial infected individuals at the first time step
param.updateAllRecovered  = true; % flag for updating all recovered or just single step
param.dt = 0.1; % time stepping
param.lr = 2e-2; % learning rate
param.hyperR = 1e-2; % define hyper regularization parameter
% param.xi = 1e1; %hyper regularization parameter
param.iterNL = 5; %number of non-linear iterations

% define graph SEIR model parameters (can be learned later on)
param.kappaS = 0.1;
param.kappaE = 0.1;
param.kappaI = 0.25;

% reaction parameters
param.alpha = 0.02;
param.beta = 0.05;
param.gamma = 0.01;
param.mu = 0.05;
param.eps = 100;

%% initialize state
y = zeros(param.n, 4, param.ts);
[y(:,:,1)] = initializeStateSEIR(param);

d = sparse(param.n, param.ts); % IGM tests
r = sparse(param.n, param.ts); % IGG tests

%% get data
G_old = [];
ys = zeros(param.n, 4, param.ts);
[ys(:,:,1)] = initializeStateSEIR(param);
for ii=1:param.ts
    [G] = getGraphData(param,G_old); % get interaction data
    DG{ii} = G; % dynamic graph construct
    
    % get test data
    [d(:,ii), dinds{ii}, Pd{ii}] = getTestDataSEIR('M',G, param,[]); % get early onset antigen (infected)
    [r(:,ii), rinds{ii}, Pr{ii}] = getTestDataSEIR('G',G, param,[]); % get immunity antigen (recovered)
    param.recovered(logical(r(:,ii))) = r(logical(r(:,ii)),ii); % update the recovered list
     ys(:,:,ii) = updateStateTestingSEIR(ys(:,:,ii),d(:,ii),r(:,ii),dinds{ii}, rinds{ii}, G, param);
     [ys(:,:,ii+1)] = evolveGraphSEIRModel(G, ys(:,:,ii), param);
    G_old = G; 
end
%% Bayesian aPCE
for jj = 1:4
     
     [yN, Idx] = find_neighors(DG{1}, ys(:,jj,4),param.N);  %get neighbours ad parameters
     [L,yt] = get_graph_states_idx(DG, ys, Idx);
     [ResponseSurface{jj},Output{jj}] = BaPC_Covid(yN, ys(:,jj,:), yt,jj, L, param);
    
end
%% %% get stats

PosteriorOutputMean=mean(ResponseSurface{2});
PosteriorOutputVar=var(ResponseSurface{2});

%% optimal testing

w = optmial_testing_cvx(PosteriorOutputMean', PosteriorOutputVar', param);

% w(w<1) = 0;
 w(w>=0.0002) = 1;

 %% Plot results
%plot(err)
figure()
stem(w)
xlabel('w_t')
title('Optimal Test prescription')
plot2styles = {'-b';'-y'; '-r'; '-g'};
idx1 = ys(:,3,1)==1;

for j=1:4
    figure()
    x= 1:param.ts+1;
    y = squeeze(ys(:,j,:));
    %y= ResponseSurface{j}';
    shadedErrorBar(x,10*y,{@mean,@std},'lineprops',plot2styles{1}); 
    %plot(squeeze(mean(ys(idx1,j,:))))
    hold on
    x= 1:param.ts;
    y= ResponseSurface{j}';
    shadedErrorBar(x,y,{@mean,@std},'lineprops',plot2styles{3});
    legend('Pior','Posterior')
    xlabel('individuals')
    ylabel('mean probability')
end
