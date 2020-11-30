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


%% Visualization of Prior and Posterior: mean value and standard deviation of model output
figure(1)
hold on;
plot(Time,Observation,'kx'); 
plot(Time,ResponseSurface(:,index_ML),'k-');
plot(Time,PriorOutputMean,'b-'); 
plot(Time,PriorOutputMean+sqrt(PriorOutputVar),'b--');
plot(Time,PriorOutputMean-sqrt(PriorOutputVar),'b-.');
plot(Time,PosteriorOutputMean,'r-'); 
plot(Time,PosteriorOutputMean+sqrt(PosteriorOutputVar),'r--');
plot(Time,PosteriorOutputMean-sqrt(PosteriorOutputVar),'r-.');
legend('Measurements','Max. Likelihood','\mu_{Prior}','\mu_{Prior}+\sigma_{Prior}','\mu_{Prior}-\sigma_{Prior}','\mu_{Post}','\mu_{Post}+\sigma_{Post}','\mu_{Post}-\sigma_{Post}');
grid on;
box on;
xlabel('Time');
ylabel('Model output');
title('Model output statistics','FontWeight', 'bold');
set(gca,'Color',[1 1 1])
set(gcf,'Color',[1 1 1]);
set(gcf,'Position',[200 600 400 300]);


%% Visualization of Likelihood
figure(2)
scatter3(Prior_distribution(1,:)',Prior_distribution(2,:)', Weight(end,:),15,Weight(end,:),'filled');
colormap(jet); colorbar;
caxis([0 max(Weight(end,:))]);
h= findobj(gca,'Type','patch');
set(h,'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0. 0. 0.]);
grid on;
xlabel('Parameter 1');
ylabel('Parameter 2');
title('Likelihood','FontWeight', 'bold');
set(gca,'Color',[1 1 1])
set(gcf,'Color',[1 1 1]);
set(gcf,'Position',[800 600 400 300]);
box on;
view(0,90);
hold on
plot3(CollocationPoints(1:end-1,1),CollocationPoints(1:end-1,2),max(Weight(:))+0.*CollocationPoints(1:end-1,2),'gx' ,'LineWidth',2,'MarkerSize',6);
