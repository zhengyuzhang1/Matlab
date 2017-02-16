function [f,gtu]=srchprd(w, t, data, k, dipsloc, nod)
gprd=floor(0.5*(4*w*t(end)/pi+1));
a=NaN(gprd);
for i=1:gprd
    for j=1:min([gprd,floor(0.5*((2*i-1)*t(end)/t(dipsloc(k))+1))])
        a(j,i)=min(abs(t(dipsloc(k))*(2*j-1)/(2*i-1)-t(dipsloc)));
    end
end
a(a>repmat(nanmedian(a,1),j,1))=NaN;
[~,i]=min(nanmean(a,1)+(t(dipsloc(k))<0.25*(2*(1:gprd)-1)*pi/w | t(dipsloc(k))>(2*(1:gprd)-1)*pi/w));
j=find(~isnan(a(:,i)));
f=zeros(size(j,1),5);
for l=1:size(j,1)
    [~,f(l,1)]=min(abs(t(dipsloc(k))*(2*j(l)-1)/(2*i-1)-t(dipsloc)));
    f(l,1)=dipsloc(f(l,1));
    f(l,2)=f(l,1);
    f(l,3)=f(l,1);
    while f(l,2)>1 && data(f(l,2)-1)>data(f(l,2))
        f(l,2)=f(l,2)-1;
    end
    while f(l,3)<length(data) && data(f(l,3)+1)>data(f(l,3))
        f(l,3)=f(l,3)+1;
    end
end
f(:,4)=j;
f(:,5)=min([abs(data(f(:,1))-data(f(:,2))),abs(data(f(:,1))-data(f(:,3)))],[],2);
f=sortrows(f,-5);
gtu=sum(f(:,5)./(1-data(f(:,1)))<0.2)<=2 || median(f(f(:,5)./(1-data(f(:,1)))<0.2,1))<dipsloc(k);
if find(f(:,1)==dipsloc(k))>nod && find(f(:,1)==dipsloc(k))<=size(f,1)
    f=f(f(:,1)==dipsloc(k),1:4);
else
    f=f(1:min(nod,size(f,1)),1:4);
end
end