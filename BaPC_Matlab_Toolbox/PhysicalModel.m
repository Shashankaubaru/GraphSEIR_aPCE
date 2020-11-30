%% Remark: you can replace it or call an external software
function ModelResponse = PhysicalModel(PhysicalSpace,P)
%Very simple example for non-linear dynamic system
t=PhysicalSpace.index;
ModelResponse=(P(1).^2+P(2)-1).^2+P(1).^2+0.1*P(1)*exp(P(2))-sqrt(0.5*t)*2*P(1)+1+sin(5*t);
