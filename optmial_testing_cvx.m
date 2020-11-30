function w = optmial_testing_cvx(postmean, postvar, param)

%normf = @(x) sum(x(:).^2);
fun = @(w) (norm(w-postmean)+ norm(w-postvar)+param.hyperR*norm(w,1));

n = param.n;
w = rand(n);
cvx_begin
   variable w(n)
   minimize(fun(w))
      subject to
        w >= 0;
        ones(n,1) >= w; 
cvx_end