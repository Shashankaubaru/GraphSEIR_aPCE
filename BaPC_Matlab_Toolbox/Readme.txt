BaPC Matlab Toolbox: Bayesian Data-driven Arbitrary Polynomial Chaos Expansion

AUTHOR: 
Sergey Oladyshkin

AFFILIATION: 
Stuttgart Center for Simulation Science,
Department of Stochastic Simulation and Safety Research for Hydrosystems, 
Institute for Modelling Hydraulic and Environmental Systems, 
University of Stuttgart, Pfaffenwaldring 5a, 70569 Stuttgart

CONTACT INFORMATION: 
E-mail: Sergey.Oladyshkin@iws.uni-stuttgart.de
Phone: +49-711-685-60116
Fax: +49-711-685-51073
Website: http://www.iws.uni-stuttgart.de

SCIENTIFIC LITERATURE:
Oladyshkin S. and Nowak W. Data-driven uncertainty quantification using the arbitrary polynomial chaos expansion. Reliability Engineering & System Safety, Elsevier, V. 106, P. 179–190, 2012.
Oladyshkin S. and Nowak W. Incomplete statistical information limits the utility of high-order polynomial chaos expansions. Reliability Engineering & System Safety, 169, 137-148, 2018.
Oladyshkin S., Class H. and Nowak W. Bayesian updating via bootstrap filtering combined with data-driven polynomial chaos expansions: methodology and application to history matching for carbon dioxide storage in geological formations. Computational Geosciences, 17(4), 671-687, 2013.
Oladyshkin S., Schroeder P., Class H. and Nowak W. Chaos expansion based Bootstrap filter to calibrate CO2 injection models. Energy Procedia, 40, 398-407, 2013

 
GENERAL INFORMATION:
BaPC Matlab Toolbox offers an advanced framework for stochastic model calibration and parameter inference based on the arbitrary polynomial chaos expansion (aPC) and strict Bayesian principles. BaPC framework follows the idea of aPC technique (Oladyshkin and Nowak, 2012) and can use arbitrary distributions for modelling parameters, which can be either discrete, continuous, or discretized continuous and can be specified either analytically, numerically as histogram or as raw data sets. BaPC Matlab Toolbox approximates the dependence of simulation model output on model parameters by expansion in an orthogonal polynomial basis and the resulting response surface can be seen as a reduced (surrogate) model. BaPC Matlab Toolbox employs an iterative Bayesian approach (Oladyshkin and Nowak, 2013) in incorporate the available measurement data and to construct the accurate reduced model in the relevant regions of high posterior probability.

SPECIFIC INFORMATION:
The main script MainRun_BaPC.m initialize the problem settings and then run the script BaPC.mat for the computations. User of the BaPC Matlab Toolbox are invited to initialize the following in the MainRun_BaPC.m script:  prior distribution (Prior_distribution), number of uncertain parameters (N), degree of expansion (d), Physical Space (PhysicalSpace), Measurments location in Physical Space (MeasurmentSpace), Measurments (Observation), Measurement Error (MeasurementError), Covariance Matrix (CovarianceMatrix) and number of Bayesian iteratoin (IterationLimit). Please, specify also the physical model in the file PhysicalModel.m or replace it by an external software.  The BaPC Matlab Toolbox call the original physical model as a Black Box model according to the initialized Physical Space. The current BaPC Matlab Toolbox offers a simple example of non-linear dynamic model determined in PhysicalModel.m. After successful computation, BaPC Matlab Toolbox computes Prior/Posterior statistics in Prior_and_Posterior_Statistics.m and provides visualization in Visualization.m.
