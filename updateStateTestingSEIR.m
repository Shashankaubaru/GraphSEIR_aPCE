function [yt] = updateStateTestingSEIR(y,d, r, dinds, rinds, G, param)
% update state per antigen testing

yt = y;
% update infected
% yt(dinds,3) = d(dinds);
yt(dinds(logical(d(dinds))),:) = repmat([0,0,1,0],nnz(logical(d(dinds))),1); %update the infected prob.
yt(dinds(~logical(d(dinds))),3) = zeros(nnz(logical(~d(dinds))),1);% updated the non-infected prob.
yt(dinds(~logical(d(dinds))),:)./sum(yt(dinds(~logical(d(dinds))),:),2); % renormalize prob.

% update recovered
% yt(rinds,4) = r(rinds);
if param.updateAllRecovered
    recovered = param.recovered; %identify those who are not known to be recovered
else
    recovered = [];
    % yt(rinds(logical(r(rinds))),:) = repmat([0,0,0,1],nnz(logical(r(rinds))),1); %update the recovered prob. (holds for all future time steps)
    %recovered = rinds(logical(r(rinds)));
end
recovered = unique([rinds(logical(r(rinds))); find(recovered)]); % if update of recovered takes place outside function, then this line is not needed
%yt(recovered,:) = repmat([0,0,0,1],nnz(recovered),1); %update the recovered prob. (holds for all future time steps)
yt(recovered,:) = repmat([0,0,0,1],length(recovered),1); %update the recovered prob. (holds for all future time steps)

yt(rinds(~logical(r(rinds))),3) = zeros(nnz(logical(~r(rinds))),1); % update the non-recovered prob. (only known at a given times step)
yt(rinds(~logical(r(rinds))),:)./sum(yt(rinds(~logical(r(rinds))),:),2); % renormalize prob.


% figure(param.f4); p=plot(G); %,'EdgeLabel',G.Edges.Weight);
% gx0=1400; gy0=400; width=550; height=550;
% set(gcf,'position',[gx0,gy0,width,height])
% p.Marker = '.';
% hold on
% plot3(p.XData,p.YData, p.ZData, 'sb'); % blue for uncertain
% infected = find(yt(:,3)>0.5);
% recovered = find(yt(:,4)>0.5);
% plot3(p.XData(infected),p.YData(infected), p.ZData(infected),'dr', 'MarkerSize',10);  % red for highly probable sick
% plot3(p.XData(recovered),p.YData(recovered), p.ZData(recovered),'og','MarkerSize',10);% green for recovered
% title('Infected and Recovered Graph')
% hold off
