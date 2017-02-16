function f=presrchprd(w, t, k, dipsloc)
gprd=floor(0.5*(4*w*t(end)/pi+1));
a=NaN(gprd);
for i=1:gprd
    for j=1:min([gprd,floor(0.5*((2*i-1)*t(end)/t(dipsloc(k))+1))])
        a(j,i)=min(abs(t(dipsloc(k))*(2*j-1)/(2*i-1)-t(dipsloc)));
    end
end
a(a>repmat(nanmedian(a,1),j,1))=NaN;
[x,y]=min(nanmean(a,1)+(t(dipsloc(k))<0.25*(2*(1:gprd)-1)*pi/w | t(dipsloc(k))>(2*(1:gprd)-1)*pi/w));
f=[x,y];
end