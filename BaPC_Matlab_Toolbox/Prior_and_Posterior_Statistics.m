%% BaPC Matlab Toolbox
% Bayesian Data-driven Arbitrary Polynomial Chaos Expansion
% Author: Sergey Oladyshkin
% Developed year: 2010
% Stuttgart Center for Simulation Science
% Department of Stochastic Simulation and Safety Research for Hydrosystems,
% Institute for Modelling Hydraulic and Environmental Systems
% University of Stuttgart, Pfaffenwaldring 5a, 70569 Stuttgart
% E-mail: Sergey.Oladyshkin@iws.uni-stuttgart.de
% Phone: +49-711-685-60116
% Fax: +49-711-685-51073
% http://www.iws.uni-stuttgart.de

% The current BaPC Matlab Toolbox use the definition of aPC that is presented in the following manuscripts: 
% Oladyshkin S. and Nowak W. Data-driven uncertainty quantification using the arbitrary polynomial chaos expansion. Reliability Engineering & System Safety, Elsevier, V. 106, P. 179–190, 2012.
% Oladyshkin S. and Nowak W. Incomplete statistical information limits the utility of high-order polynomial chaos expansions. Reliability Engineering & System Safety, 169, 137-148, 2018.

% BaPC Matlab Toolbox is also use the iterative Bayesian framework that is presented in the following manuscripts: 
% Oladyshkin S., Class H. and Nowak W. Bayesian updating via bootstrap filtering combined with data-driven polynomial chaos expansions: methodology and application to history matching for carbon dioxide storage in geological formations. Computational Geosciences, 17(4), 671-687, 2013.
% Oladyshkin S., Schroeder P., Class H. and Nowak, W. Chaos expansion based Bootstrap filter to calibrate CO2 injection models. Energy Procedia, 40, 398-407. 2013

clear all;
load('BaPC.mat');

%% Prior output statistics: mean and variance
for j=1:1:NumberOfMeasurments;    
    PriorOutputMean(j)=Coeffs(1,j,end); %Analytical form of mean based on the last iteration 
    PriorOutputVar(j)=sum(Coeffs(2:P,j,end).^2); %Analytical form of variance based on the last iteration 
end

%% Posterior Estimation via Rejection Sampling
ii=0;
for i=1:MCsize
    if Weight(i)/max(Weight(end,:))>rand(1)
        ii=ii+1;
        Posterior_distribution(:,ii)=Prior_distribution(:,i);
        for NofM=1:NumberOfMeasurments; 
            Posterior_ResponseSurface(NofM,ii)=ResponseSurface(NofM,i); %Response Surface from the last iteration
        end
    end
end
% Posterior output statistics: mean and variance
for NofM=1:1:NumberOfMeasurments;    
    PosteriorOutputMean(NofM)=mean(Posterior_ResponseSurface(NofM,:)); %Posterior mean 
    PosteriorOutputVar(NofM)=var(Posterior_ResponseSurface(NofM,:)); %Posterior variance
end
