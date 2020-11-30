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

%% Construction of Data-driven Arbitrary Orthonormal Polynomial Basis
function OrthonormalBasis = aPC_OrthonormalBasis(Data, Degree)
% Input:
% Data - raw data array
% Degree - degree of the orthonormal polynomial basis
% Output:
% OrthonormalBasis - orthonormal polynomial basis 

%% Initialization
d=Degree; %Degree of polynomial expansion
dd=d+1; %Degree of polynomial defenition
L_norm=1; % L-norm for polynomial normalization
NumberOfDataPoints=length(Data);

%% Forward linear transformation
MeanOfData=mean(Data);
VarOfData=var(Data);
Data=Data/MeanOfData;

%% Raw Moments of the Input Data
for i=0:(2*dd+1)
    m(i+1)=sum(Data.^i)/length(Data); 
end

%% Polynomial up to degree dd
for degree=0:dd;
    %Defenition of Moments Matrix Mm
    for i=0:degree;
        for j=0:degree;                    
            if (i<degree) 
                Mm(i+1,j+1)=m(i+j+1); 
            end            
            if (i==degree) && (j<degree)
                Mm(i+1,j+1)=0;
            end
            if (i==degree) && (j==degree)
                Mm(i+1,j+1)=1;
            end        
        end
        Mm(i+1,:)=Mm(i+1,:)/max(abs(Mm(i+1,:))); %Matrix Normalization 
    end
    %Defenition of ortogonality conditions Vc
    for i=0:degree;
        if (i<degree) 
             Vc(i+1)=0; 
        end            
        if (i==degree)
            Vc(i+1)=1;
        end
    end
  
    %Coefficients of Non-Normal Orthogonal Polynomial: Vp
    inv_Mm=pinv(Mm);
    Vp=Mm\Vc';
    PolyCoeff_NonNorm(degree+1,1:degree+1)=Vp';   
       
    fprintf('=> aPC Toolbox: Inaccuracy for polynomial basis of degree %1i is %5f pourcents',degree, 100*abs(sum(abs(Mm*PolyCoeff_NonNorm(degree+1,1:degree+1)'))-sum(abs(Vc)))); 
    if 100*abs(sum(abs(Mm*PolyCoeff_NonNorm(degree+1,1:degree+1)'))-sum(abs(Vc)))>0.5
        fprintf('\n=> Warning: Computational error of the linear solver is too high.');           
    end
    fprintf('\n');
    

    %Normalization of polynomial coefficients
    P_norm=0;
    for i=1:NumberOfDataPoints;        
        Poly=0;
        for k=0:degree;
            Poly=Poly+PolyCoeff_NonNorm(degree+1,k+1)*Data(i)^k;     
        end
        P_norm=P_norm+Poly^2/NumberOfDataPoints;        
    end
    P_norm=sqrt(P_norm);
    for k=0:degree;
        Polynomial(degree+1,k+1)=PolyCoeff_NonNorm(degree+1,k+1)/P_norm;
    end
end

%% Backward linear transformation to the data space
Data=Data*MeanOfData;
for k=1:length(Polynomial);
    Polynomial(:,k)=Polynomial(:,k)/MeanOfData^(k-1);
end

%% Data-driven Arbitrary Orthonormal Polynomial Basis
OrthonormalBasis=Polynomial;
