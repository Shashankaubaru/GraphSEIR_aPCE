
lam = [ 1e-0, 1e-2, 1e-4];

for i =1:3
param.hyperR = lam(i);

%w(i,:) = optmial_testing_cvx(PosteriorOutputMean', PosteriorOutputVar', param);

% normf = @(x) sum(x(:).^2);
% 
% y = covidAssimilation_cvx(ys, d, r, Pd, Pr, DG, param);
% 
% for ii=1:param.ts
%     err(ii) = normf(y(:,:,ii)- ys(:,:,ii));
% end
% w(w<1) = 0;
w1 = w(i,:);
if i ==1
    w1(w1>=0.001) = 1;
elseif i==2
    w1(w1>=0.0002) = 1;
else
    w1(w1>=0.00005) = 1;
end

 subplot(3,1,i)
 stem(w1)
 xlabel('w_t')
 tit = sprintf('Optimal Test prescription \\lambda = %.4f', lam(i));
 title(tit)
end

% obj = [1.079562e+01,5.479121e+00,1.991465e+00,3.501246e-01,5.361193e-02, 4.858157e-02, 5.489051e-02];
% 
% figure()
%  semilogx(lam, obj, 'LineWidth',2)
%  set(gca, 'XDir','reverse')
%  xlabel('\lambda')
%  ylabel('Objective function')
 
