function [yp] = evolveModel(L, xt, yt, jj, param)

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
%L = (G.laplacian); % should remain sparse, but AD cannot handle sparse matrices... 

% extract state components
S = yt(:,1); % Susceptive
E = yt(:,2); % Exposed
I = yt(:,3); % Infected
R = yt(:,4); % Recovered

% update state
if jj==1
    S= xt';
    yp = S + dt *(-kappaS * L * S - beta * E.*S - gamma * I .* S);
elseif jj==2
    E= xt';
    yp = E + dt *(-kappaE * L * E + beta * E.*S + gamma * I .* S - alpha *E);
elseif jj == 3
    I = xt';
    yp = I + dt *(-kappaI * L * I + alpha * E - mu * I);
elseif jj == 4
    R = xt';
    yp = R + dt *(mu * I);
end
%yp = [Sp, Ep, Ip, Rp];

% normalize model 
yp = yp./sum(yp);