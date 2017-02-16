function f = rmerrdata2(data)
f=data(2:end-1);
absdiff2=abs(diff(data,2));
while max(absdiff2)>0.3
    f=f.*(absdiff2<=0.3)+(absdiff2>0.3).*([data(1);f(1:end-1)]+[f(2:end);data(end)])/2;
    absdiff2=abs(diff([data(1);f;data(end)],2));  
end
f=[data(1);f;data(end)];
error=1;
while error>0
    dif=diff(f);
    absdiff=abs(dif);
    errorq=logical([0;0;(dif(2:end-2).*dif(3:end-1)<0) & (((absdiff(2:end-2)./absdiff(1:end-3)>3) & (absdiff(3:end-1)./absdiff(4:end)>3)) | max([absdiff(2:end-2)./absdiff(1:end-3),absdiff(3:end-1)./absdiff(4:end)],[],2)>5);0;0]);
    newdata=[f(1:2);(f(2:end-3)+f(4:end-1))/2;f(end-1:end)];
    f(errorq)=newdata(errorq);
    error=sum(errorq);
    f(f>1)=1;
end
end