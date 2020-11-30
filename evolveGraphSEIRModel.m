function [yp] = evolveGraphSEIRModel(G, yt, param)
% extract model parameters (these should be learnt at a later stage(
kappaS = param.kappaS;
kappaE = param.kappaE;
kappaI = param.kappaI;
alpha = param.alpha;
beta = param.beta;
gamma = param.gamma;
mu = param.mu;
dt = param.dt;
% construct graph laplacian
L = (G.laplacian); % should remain sparse, but AD cannot handle sparse matrices... 

% extract state components
S = yt(:,1); % Susceptive
E = yt(:,2); % Exposed
I = yt(:,3); % Infected
R = yt(:,4); % Recovered

% update state
Sp = S + dt *(-kappaS * L * S - beta * E.*S - gamma * I .* S);
Ep = E + dt *(-kappaE * L * E + beta * E.*S + gamma * I .* S - alpha *E);
Ip = I + dt *(-kappaI * L * I + alpha * E - mu * I);
Rp = R + dt *(mu * I);

yp = [Sp, Ep, Ip, Rp];

yp(yp<0) = 0;
% normalize model 
yp = yp./sum(yp,2);

