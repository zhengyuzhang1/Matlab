function [rlt,gos] = fsrch (w, n, azx, range, step, data, t)
m1=round(range(1)/step(1));
m2=round(range(2)/step(2));
domain=[repmat(step(1)*(-m1:m1),1,2*m2+1)',reshape(repmat(step(2)*(-m2:m2),2*m1+1,1),(2*m1+1)*(2*m2+1),1)]+repmat(azx,(2*m1+1)*(2*m2+1),1);
domain(any(abs(domain)<4e4*pi,2),:)=[];
errpt=abs(bsxfun(@rdivide,data',M(w,n,domain,t')));
a=sum(errpt>1,2);
gos=min(a);
minp=find(a==gos);
rlt = domain(minp(round(end/2)),:);
end     