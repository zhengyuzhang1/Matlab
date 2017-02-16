w = 5e5;
n = 2:2:80;
t = 1e-8:1e-8:10e-6;
m = zeros(length(n),length(t));
for i = 1:length(n)
    m(i,:) = prod(M(w,n(i),sample,t));
end
fid = zeros(size(n));
%% 
minm = min(m,[],1);
mint = find(diff(minm)>0,1);
[~,minn] = min(m(:,mint));
figure;plot(m(:,mint));
 
t1 = t(mint)+(-10e-8:5e-9:10e-8);
figure;plot(t1,M(w,n(minn),sample(1,:),t1))
figure;plot(t1,prod(M(w,n(minn),sample,t1)))
%% 
%minn = 20; 
fitfun = @(azx, t) M(w,n(minn),azx,t);
fitazx0 = sample(1,:) * 1.05 ;
fitt = t1;
fitdata = prod(M(w,n(minn),sample,t1));
[~,mid] = min(fitdata);
l = max(1,mid - 1);
r = min(length(fitt),mid + 1);
while l > 1 && fitdata(l) < fitdata(l-1)
    l = l - 1;
end
while r < length(fitt) && fitdata(r) < fitdata(r+1)
    r = r + 1;
end
fitt = fitt(l:r);
fitdata = fitdata(l:r);
azx = lsqcurvefit(fitfun,fitazx0,fitt,fitdata);
(azx - sample(1,:))./sample(1,:)
figure;plot(fitt,fitdata,fitt,M(w,n(minn),azx,fitt));
figure;plot(fitt,fitdata,fitt,M(w,n(minn),sample(1,:),fitt));