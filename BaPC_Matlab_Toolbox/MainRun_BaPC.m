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
fprintf('\n=> BaPC Toolbox: Initialization ... \n');
tic

%% Initialization of polynomial chaos expansion
N=2; %Number of uncertain parameters
d=3; %Order of expansion 
P=factorial(N+d)/(factorial(N)*factorial(d)); %Total number of terms

%% Initialization of Prior distribution
MCsize=1000; %Sample size
Prior_distribution(1,:)=2*random('Beta',2,1,1,MCsize)-1; %This is an example, you can load your data, etc.
Prior_distribution(2,:)=2*random('Beta',2,0.5,1,MCsize)-1; %This is an example, you can load your data, etc.etc.

%% Initialization of Physical Space
Time=0:0.05:1;
PhysicalSpace=struct('index', Time);

%% Initialization of Measurments location in Physical Space
NumberOfMeasurments=length(Time); %Synthetic Measurment Time Series
MeasurmentSpace=struct('index', 1:NumberOfMeasurments);

%% Initialization of Synthetic Measurments 
SyntheticParameters = [0,0]; %Reference Parameters for Synthetic Measurments
Observation=PhysicalModel(PhysicalSpace,SyntheticParameters); %Synthetic Measurment Time Series

%% Initialization of Bayesian Framework
MeasurementError=0.5; %Initialization of Measurement Error  
CovarianceMatrix=MeasurementError.^2*eye(NumberOfMeasurments); %Initialization of Covariance Matrix
IterationLimit=5; %Initialization of Number of Bayesian Iteratoins 

%% Run the kernel BaPC code
BaPC

%% Computation of Prior and Posterior statistics: mean and variance
Prior_and_Posterior_Statistics

%% Visualization of Prior, Posterior and Likelihood
Visualization