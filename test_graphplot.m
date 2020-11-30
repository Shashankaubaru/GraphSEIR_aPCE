close all
param.n = 100; % for now assume graph nodes are fixed, later nodal dynamics can be incorporated
param.ts = 20; % timesteps

param.N = 5; % No. of neighbours
param.d = 3; %Order of expansion

param.recovered = zeros(param.n,1); % recovered individuals
param.ny0 = 3; % number of initial infected individuals at the first time step
param.updateAllRecovered  = true; % flag for updating all recovered or just single step
param.dt = 0.1; % time stepping
param.lr = 2e-2; % learning rate
param.hyperR = 1e-3; % define hyper regularization parameter
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


G_old = [];
ys = zeros(param.n, 4, param.ts);
[ys(:,:,1)] = initializeStateSEIR(param);
    figure()
for ii=1:param.ts
    [G] = getGraphData(param,G_old); % get interaction data
    DG{ii} = G;
    G = DG{ii} ;
    [ys(:,:,ii+1)] = evolveGraphSEIRModel(G, ys(:,:,ii), param);
    title('State evolution')
    xlabel('time steps')
    
    subplot(4,1,1),imagesc(squeeze(ys(:,1,:))); ylabel('susceptive prob.'),colorbar
    subplot(4,1,2),imagesc(squeeze(ys(:,2,:))); ylabel('exposed prob.'),colorbar
    subplot(4,1,3),imagesc(squeeze(ys(:,3,:))); ylabel('infected prob.'),colorbar
    subplot(4,1,4),imagesc(squeeze(ys(:,4,:))); ylabel('recovered prob.'),colorbar
    
    G_old = G;
end
ys(ys<1e-4) = 0;
for ii=1:10:param.ts
    G= DG{ii};
    for jj= 1:param.n
        a{jj} = mat2str(ys(jj,:,ii),2);
    end
    b = ys(:,3,ii)>0.5;
    c= ys(:,3,ii)>0.03;
    e = ys(:,3,ii)>1e-3;
    figure()
    h= plot(G);
    h.LineWidth = 2;
    h.MarkerSize = 10;
    labelnode(h,1:param.n,a)
    highlight(h,e,'NodeColor','y')
    highlight(h,c,'NodeColor','m')
    highlight(h,b,'NodeColor','r')
end
