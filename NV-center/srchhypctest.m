function f = srchhypctest(w, n, t, data, data0, nuc)
[dps,locs]=findpeaks(-data);
dps=-dps;
locspart=locs(dps<0.85);
dipq=zeros(length(locspart),2);
for i=1:length(locspart)
    dipq(i,:)=presrchprd(w,t,i,locspart);
end
[~,mindloc]=min(dipq(:,1));
mindv=data(locspart(mindloc));
dips=[locs(1-dps>0.3*(1-mindv)),dps(1-dps>0.3*(1-mindv))];
prds=srchprd(w,t,data,find(dips(:,2)==mindv),dips(:,1),5);
azx=fitdip(w,n,t,data,prds);
azx=min(azx,2*pi*[500e3,500e3]);
gos=1;
range=abs(azx./[10,3]);
data1=data0;
if azx(2)<5e5
    range=abs(azx).*[0.1,0.3];
    data1=data0;
    [~,gos]=fsrch(w,n,azx,range,range/10,data1,t);
end
if ~isempty(nuc) && gos~=0
    [minchnuc,loc]=min(max(abs((repmat(azx,size(nuc,1),1)-nuc)./nuc),[],2));
else
    minchnuc=1;
end
if minchnuc<0.05
    range=3*abs(nuc(loc,:)-azx);
end
step=min(range/2,max(min(2000*pi*[1,1],range/30),10*2*pi*[1,1]));
while gos>0 && min(step)>2*pi
    [azx,gos]=fsrch(w,n,azx,range,step,data1,t);
    range=range/2;
    step=step/2;
end
%[azx,gofs]=sfsrch(w,n,azx,range,step,data,t0,deltat,tn);
gofs=0;
f=[azx;gos,gofs];    
end