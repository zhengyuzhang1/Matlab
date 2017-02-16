function f=axes0(w,azx,t,axis)
%axis=1 gives nz, axis=2 gives nx
if size(t,1)==1
    t0=repmat(t,size(azx,1),1);
else
    t0=t;
end
f=zeros(size(t0));
for j=1:size(t0,2)
    for i=1:size(t0,1)
        logmat=logm(U0(w,azx(i,:),t0(i,j)))/rotang(w,azx(i,:),t0(i,j));
        f(i,j)=-2*imag(logmat(1,axis));
    end
end
        
        
