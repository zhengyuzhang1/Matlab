function f = srchhypcall(w, n, data, data0, t0, deltat, nuc)
tn=t0-deltat+deltat*length(data);
t=(t0:deltat:tn)';
[dps,locs]=findpeaks(-data);
dps=-dps;
alldips=[t(locs),dps];
mid=4;lp=4;rp=4;mindv=dps(2);
while (find(dps==mindv)>1 && 1-mindv<0.5*(1-dps(find(dps==mindv)-1))) || (find(dps==mindv)<size(dps,1) && (1-mindv)<0.5*(1-dps(find(dps==mindv)+1))) || (data(mid)>0 && (data(lp)<0.95 || data(rp)<0.95)) || mid-lp<2 || rp-mid<2 || abs(data(mid)-data(lp))/abs(data(mid)-data(rp))<0.5 || abs(data(mid)-data(lp))/abs(data(mid)-data(rp))>2 || (abs(data(mid)-data(mid-2))<0.015 && abs(data(mid)-data(mid+2))<0.015) || abs(data(mid)-data(mid+3))<0.02 || abs(data(mid)-data(mid-3))<0.02 || isnan(abs(data(mid)-data(lp))/abs(data(mid)-data(rp)))
    [mindv,mindloc]=min(alldips(:,2));
    mid=round(alldips(mindloc,1)/deltat);
    lp=mid;rp=mid;
    while lp>1 && data(lp-1)>data(lp)
        lp=lp-1;
    end
    while rp<length(data) && data(rp+1)>data(rp)
        rp=rp+1;
    end
    alldips(mindloc,:)=[];
end
dips=[t(locs(1-dps>0.3*(1-mindv))),dps(1-dps>0.3*(1-mindv))];
prd=srchprdall(w,find(dips(:,2)==mindv),dips(:,1),tn);
fitt=t(max(round((lp+mid)/2),mid-3):min(round((mid+rp)/2),mid+3));
if min([data(lp),data(rp)])>0.9
    fitdata=data(max(round((lp+mid)/2),mid-3):min(round((mid+rp)/2),mid+3))/min([data(lp),data(rp)]);
else
    fitdata=data(max(round((lp+mid)/2),mid-3):min(round((mid+rp)/2),mid+3));
end
fitfun= @(azx, t) M(w,n,azx,t);
fitazx0=[(2*prd-1)*pi/(t0+(mid-1)*deltat)-2*w,100000*pi];
fitlb(2)=4e4*pi;
fitub(2)=2e6*pi;
fitlb(1)=min(0.5*fitazx0(1),1.5*fitazx0(1));
fitub(1)=max(0.5*fitazx0(1),1.5*fitazx0(1));
azx=lsqcurvefit(fitfun,fitazx0,fitt,fitdata,fitlb,fitub,optimset('MaxFunEvals',400,'Display','Off'));
azx=min(azx,2*pi*[500e3,500e3]);
gos=1;
range=max(abs(azx/3),2*pi*30e3*[1,1]);
data1=data;
if azx(2)<5e5
    range=abs(azx).*[0.1,0.3];
    data1=data0;
    [~,gos]=fsrch(w,n,azx,range,range/10,data1,t0,deltat,tn);
end
if ~isempty(nuc) && gos~=0
    [minchnuc,loc]=min(max(abs((repmat(azx,size(nuc,1),1)-nuc)./nuc),[],2));
    if minchnuc<0.15
        azx=lsqcurvefit(fitfun,nuc(loc,:),fitt,fitdata,fitlb,fitub,optimset('MaxFunEvals',400,'Display','Off'));
        [minchnuc,loc]=min(max(abs((repmat(azx,size(nuc,1),1)-nuc)./nuc),[],2));
    end
else
    minchnuc=1;
end
if minchnuc<0.05
    range=3*abs(nuc(loc,:)-azx);
end
step=min(range/2,max(min(2000*pi*[1,1],range/30),10*2*pi*[1,1]));
while gos>0 && min(step)>2*pi
    [azx,gos]=fsrch(w,n,azx,range,step,data1,t0,deltat,tn);
    range=range/2;
    step=step/2;
end
%[azx,gofs]=sfsrch(w,n,azx,range,step,data,t0,deltat,tn);
gofs=0;
f=[azx;gos,gofs];    
end