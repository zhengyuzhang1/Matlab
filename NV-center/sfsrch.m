function [rlt,gos] = sfsrch (omiga, n, azx, range, step, data, t0, deltat, tn)
m1=round(range(1)/step(1));
m2=round(range(2)/step(2));
domain=[repmat(step(1)*(-m1:m1),1,2*m2+1)',reshape(repmat(step(2)*(-m2:m2),2*m1+1,1),(2*m1+1)*(2*m2+1),1)]+repmat(azx,(2*m1+1)*(2*m2+1),1);
t = t0:deltat:tn;
[gos,minp]=min(sum(abs(diff(repmat(data',size(domain,1),1)./M(omiga,n,domain,t),2,2)>0.1),2));
rlt = domain(minp,:);
end     