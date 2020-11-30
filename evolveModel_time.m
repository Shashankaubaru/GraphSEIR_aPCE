function [yp] = evolveModel_time(L, xt, yt, jj, param)

st = xt;
for t = 1:length(L)
    
     ytemp = evolveModel(L{t}, st, yt(:,:,t), jj, param);
     st = ytemp';
     yp(t) = ytemp(1);
end
     