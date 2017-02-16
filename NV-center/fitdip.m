function f = fitdip(w, n, t, data, prds)
nod=size(prds,1);
azx=zeros(nod,2);
for i=1:nod
    mid=prds(i,1);
    lp=prds(i,2);
    rp=prds(i,3);
    fitt=t(max(round((lp+mid)/2),mid-3):min(round((mid+rp)/2),mid+3));
    fitdata=data(max(round((lp+mid)/2),mid-3):min(round((mid+rp)/2),mid+3));
    fitfun= @(azx, t) M(w,n,azx,t);
    fitazx0=[(2*prds(i,4)-1)*pi/t(mid)-2*w,100000*pi];
    fitlb(2)=4e4*pi;
    fitub(2)=2e6*pi;
    fitlb(1)=min(0.5*fitazx0(1),1.5*fitazx0(1));
    fitub(1)=max(0.5*fitazx0(1),1.5*fitazx0(1));
    [azx(i,1:2),azx(i,3)]=lsqcurvefit(fitfun,fitazx0,fitt,fitdata,fitlb,fitub,optimset('MaxFunEvals',400,'Display','Off'));
    azx(i,3)=azx(i,3)/length(fitt);
end
azmat=repmat(azx(:,1),1,nod);
axmat=repmat(azx(:,2),1,nod);
mp=abs((azmat-azmat')./azmat)+abs((axmat-axmat')./axmat);
mp=mp<0.02;
groupq=zeros(nod,1);
groupn=0;
stack=zeros(nod,2);
for i=1:nod
    if groupq(i)==0
        groupn=groupn+1;
        groupq(i)=groupn;
        top=1;
        stack(1,:)=[i,1];
        while top>0 
            while stack(top,2)<nod
                stack(top,2)=stack(top,2)+1;
                if groupq(stack(top,2))==0 && mp(stack(top,1),stack(top,2))==1
                    groupq(stack(top,2))=groupq(stack(top,1));
                    top=top+1;
                    stack(top,1)=stack(top-1,2);
                    stack(top,2)=1;
                end
            end
        top=top-1;
        end
    end
end
[maxgn,maxg]=max(sum(repmat(groupq,1,max(groupq))==repmat(1:max(groupq),nod,1)));
if maxgn>nod/2
    f=mean(azx(groupq==maxg,:),1);
else
    dipvs=data(prds(:,1));
    if sum(dipvs<0)>sum(dipvs>0)
        f=[median(azx(:,1)),max(azx(dipvs<0,2)),mean(azx(:,3))];
    else
        f=[median(azx(:,1)),min(azx(dipvs>0,2)),mean(azx(:,3))];
    end
end
end