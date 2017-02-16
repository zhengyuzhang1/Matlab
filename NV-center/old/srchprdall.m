function f=srchprdall(w, k, dipsloc, range)
gprd=floor(0.5*(4*w*range/pi+1));
a=NaN(gprd);
for i=1:gprd
    for j=1:min([gprd,floor(0.5*((2*i-1)*range/dipsloc(k)+1))])
        a(j,i)=min(abs(dipsloc(k)*(2*j-1)/(2*i-1)-dipsloc));
    end
end
a(a>repmat(nanmedian(a,1),j,1))=NaN;
[~,f]=min(nanmean(a,1)+(dipsloc(k)<0.25*(2*(1:gprd)-1)*pi/w | dipsloc(k)>(2*(1:gprd)-1)*pi/w));
end