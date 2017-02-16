function f=srchprd3(w, k, dips, range)
gprd=floor(0.5*(4*w*range/pi+1));
a=NaN(gprd);
b=zeros(gprd);
for i=1:gprd
    for j=1:min([floor(0.5*((2*i-1)*range/dips(k,1)+1)),gprd])
        [a(j,i),b(j,i)]=min(abs(dips(k,1)*(2*j-1)/(2*i-1)-dips(:,1)));
    end
    if (i>1 && ~isnan(a(i-1,i)) && dips(b(i,i),2)>dips(b(i-1,i),2)) || (i<gprd && ~isnan(a(i+1,i)) && dips(b(i,i),2)>dips(b(i+1,i),2)) 
        a(i,i)=1;
    end
end
a(a>repmat(nanmedian(a,1),j,1))=NaN;
[~,f]=min(nanmean(a,1)+(dips(k,1)<0.25*(2*(1:gprd)-1)*pi/w | dips(k,1)>(2*(1:gprd)-1)*pi/w));
end