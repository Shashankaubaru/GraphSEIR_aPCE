clc
close all;%clear all;
dbstop if error
addpath(genpath('cvx'))
addpath(genpath('BaPC_Matlab_Toolbox'))
%% set paramters
param.n = 1000; % for now assume graph nodes are fixed, later nodal dynamics can be incorporated
param.ts = 10; % timesteps

param.N = 5; % No. of neighbours
param.d = 3; %Order of expansion

param.recovered = zeros(param.n,1); % recovered individuals
param.ny0 = 5; % number of initial infected individuals at the first time step
param.updateAllRecovered  = true; % flag for updating all recovered or just single step
param.dt = 0.1; % time stepping
param.lr = 2e-2; % learning rate
param.hyperR = 1e-4; % define hyper regularization parameter
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
ys2 = zeros(param.n, 4, param.ts);

[ys(:,:,1)] = initializeStateSEIR(param);
ys2(:,:,1)  = ys(:,:,1);
for ii=1:param.ts*10
    [G] = getGraphData(param,G_old); % get interaction data
    DG{ii} = G; % dynamic graph construct
    
    % get test data
    [d(:,ii), dinds{ii}, Pd{ii}] = getTestDataSEIR('M',G, param, []); % get early onset antigen (infected)
    [r(:,ii), rinds{ii}, Pr{ii}] = getTestDataSEIR('G',G, param,[]); % get immunity antigen (recovered)
    param.recovered(logical(r(:,ii))) = r(logical(r(:,ii)),ii); % update the recovered list
    ys2(:,:,ii) = updateStateTestingSEIR(ys2(:,:,ii),d(:,ii),r(:,ii),dinds{ii}, rinds{ii}, G, param);
    [ys2(:,:,ii+1)] = evolveGraphSEIRModel(G, ys2(:,:,ii), param);
    if mod(ii,10)==0
        [yN, Idx] = find_neighors(DG{1}, ys(:,2,mod(ii,10)+4),param.N);
        [L,yt] = get_graph_states_idx(DG, ys, Idx);
        [ResponseSurface,Output] = BaPC_Covid(yN, ys(:,2,:), yt,2, L, param);
        PosteriorOutputMean=mean(ResponseSurface);
        PosteriorOutputVar=var(ResponseSurface);
        w = optmial_testing_cvx(PosteriorOutputMean', PosteriorOutputVar', param);
        inds = find(w>0.0001);
        [d(:,ii), dinds{ii}, Pd{ii}] = getTestDataSEIR('M',G, param, inds); % get early onset antigen (infected)
        ys(:,:,ii) = updateStateTestingSEIR(ys(:,:,ii),d(:,ii),r(:,ii),dinds{ii}, rinds{ii}, G, param);
        [ys(:,:,ii+1)] = evolveGraphSEIRModel(G, ys(:,:,ii), param);
    else
         ys(:,:,ii) = updateStateTestingSEIR(ys(:,:,ii),d(:,ii),r(:,ii),dinds{ii}, rinds{ii}, G, param);
        [ys(:,:,ii+1)] = evolveGraphSEIRModel(G, ys(:,:,ii), param);  
    end
     G_old = G;
end

plot2styles = {'-b';'-y'; '-r'; '-g'};

for j=2
    figure()
    x= 1:param.ts*10;
    y = squeeze(ys3(:,j,1:end-1));
    %y= ResponseSurface{j}';
    shadedErrorBar(x,y,{@mean,@std},'lineprops',plot2styles{1});
    %plot(squeeze(mean(ys(idx1,j,:))))
    hold on
    x= 1:param.ts*10;
    y = squeeze(ys(:,j,1:end-1));
    shadedErrorBar(x,y,{@mean,@std},'lineprops',plot2styles{3});
    legend('Random testing','Optimal testing')
    xlabel('time -->')
    ylabel('mean probability')
end
