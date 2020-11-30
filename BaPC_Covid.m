function [ResponseSurface,Output] = BaPC_Covid(Prior_distribution, Observation, yt,jj, L, param)

%%This function computes BaPC - Bayesian Data-driven Arbitrary Polynomial Chaos Expansion
%% Get paramters

[N, MCsize] = size(Prior_distribution);
d = param.d;
P=factorial(N+d)/(factorial(N)*factorial(d));

MeasurementError=0.5; %Initialization of Measurement Error  
CovarianceMatrix=MeasurementError.^2*eye(param.ts); %Initialization of Covariance Matrix

Observation = squeeze(Observation);
%% Construction of orthonormal aPC polynomial basis 
for i=1:N
    fprintf('=> BaPC Toolbox: Construction of data-driven polynomial basis for parameter %1i ... \n', i);
    OrthonormalBasis(:,:,i)=aPC_OrthonormalBasis(Prior_distribution(i,:), d);
end
%Construction of corresponding collocation points 
for i=1:N
    Polynomial=OrthonormalBasis(:,:,i);
    AvailableCollocationPoints(i,:)=roots(fliplr(Polynomial(d+2,:))); %Roots of d+1 polynomials
end

%% Construction of optimal integration points
%Unique combination of collocation points
[List{N:-1:1}]=ndgrid(1:d+1);
PointsVector=1:d+1;
UniqueCombinations=PointsVector(reshape(cat(N+1,List{:}),[],N));        
for i=1:1:length(UniqueCombinations) 
    DigitalPointsWeight(i)=sum(UniqueCombinations(i,:));                 
end
%Sorting of weights
[SortDigitalPointsWeight, index_SDPW]=sort(DigitalPointsWeight);
SortUniqueCombinations=UniqueCombinations(index_SDPW,:);
%Relative ranking via mean value 
for j=1:N    
    temp(j,:)=abs(AvailableCollocationPoints(j,:)-mean(Prior_distribution(j,:)));
end
[temp_sort, index_CP]=sort(temp,2);
for j=1:N    
    SortAvailableCollocationPoints(j,:)=AvailableCollocationPoints(j,index_CP(j,:));
end
%Mapping of unique combination to collocation points
for i=1:1:length(SortUniqueCombinations) 
    for j=1:N
        SortUniqueCombinations(i,j)=SortAvailableCollocationPoints(j,SortUniqueCombinations(i,j));
    end
end
%Optimal P-Collocation Points
CollocationPoints=SortUniqueCombinations(1:P,:);

%% Construction of Polynomial Degrees 
PosibleDegree=0:1:length(AvailableCollocationPoints(1,:))-1;
for i=2:1:N
    PosibleDegree=[PosibleDegree,0:1:length(AvailableCollocationPoints(i,:))-1];    
end
UniqueDegreeCombinations=unique(nchoosek(PosibleDegree,N),'rows');
%Possible degree computation
for i=1:1:length(UniqueDegreeCombinations) 
    DegreeWeight(i)=0;
    for j=1:1:N
        DegreeWeight(i)=DegreeWeight(i)+UniqueDegreeCombinations(i,j);
    end
end
%Sorting of possible degree
[SortDegreeWeight, i]=sort(DegreeWeight);
SortDegreeCombinations=UniqueDegreeCombinations(i,:);
%Set up of polynomail degrees 
PolynomialDegree=SortDegreeCombinations(1:P,:);

%% Initial evaluation of the original physical model and reading of corresponding model outputs
fprintf('=> BaPC Toolbox: Initial evaluation of the original physical model ... \n');
for i=1:1:P    
    Output(i,:)=evolveModel_time(L, CollocationPoints(i,:), yt, jj, param);
end

%% Pre-computation: Orthonormal Basis In all MC Points
OrthonormalBasisInMCPoints=zeros(MCsize,N,d);
for i=1:MCsize;  
    for ii=1:N;  
        for l=0:d
            OrthonormalBasisInMCPoints(i,ii,l+1)=polyval(OrthonormalBasis(l+1,length(Polynomial):-1:1,ii),Prior_distribution(ii,i));
        end
    end
end

%% Bayesain Iterations 
Weight=zeros(1,MCsize);
for index=1:param.iterNL+1
    
    %% Updating values of P_total    
    Size=size(CollocationPoints);
    P_total=Size(1); 
   
    %% Setting up of space/time independent matrix of arbitrary polynomials P*P_total
    clear Psi;
    Psi=zeros(P_total,P);
    for i=1:1:P;        
        for j=1:1:P_total;        
            Psi(j,i)=1;
            for ii=1:1:N;                                        
                Psi(j,i)=Psi(j,i)*polyval(OrthonormalBasis(PolynomialDegree(i,ii)+1,length(Polynomial):-1:1,ii),CollocationPoints(j,ii));
            end        
        end
    end

    %% Computation of the multi-dimensional expansion coefficients for arbitrary Polynomial Chaos
    clear Coeffs inv_PsiPsi temp ;       
    inv_PsiPsi=pinv(Psi'*Psi);
    for ip=1:param.ts
        for ii=1:1:P_total
            temp(ii)=Output(ii,ip);
        end                                         
        Coeffs(1:P,ip) = inv_PsiPsi*Psi'*temp';                                        
    end

    %% Respose Surfance in all MC points
    for ip=1:param.ts
        for i=1:1:MCsize
            ResponseSurface(ip,i)=0;
            for j=1:1:P;
                multi=1;
                for ii=1:1:N;          
                    multi=multi*OrthonormalBasisInMCPoints(i,ii,PolynomialDegree(j,ii)+1);
                end             
                ResponseSurface(ip,i)=ResponseSurface(ip,i)+Coeffs(j,ip)*multi;
            end 
        end
    end

    %% Bayesian Updating
    %Deviation of Model from Observations
    for im=1:param.ts
        Deviation(im,:)=Observation(:,im)'-ResponseSurface(im,:); 
    end
    %Computation of Weight/Likelihood according to the Deviation  
    for i=1:1:MCsize     
        Weight(i)=1/(sqrt(2*pi)*MeasurementError)^param.ts*exp(-0.5*Deviation(:,i)'*inv(CovarianceMatrix)*Deviation(:,i));
    end
    
    %% Selection of new Collocation Point as current Maximum Likelihood Point
    [Value,index_ML]=max(Weight);            
    CollocationPoints(P_total+1,:)=Prior_distribution(:,index_ML)';
        
    %% Evaluation of the original physical model in the new collocation point
    Output(P_total+1,:,:,:,:)=evolveModel_time(L, CollocationPoints(P_total+1,:), yt, jj, param);
    
    ResponseSurface = ResponseSurface./sum(ResponseSurface,1);
end


