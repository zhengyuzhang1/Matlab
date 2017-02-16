function f=newaxes0(w,azx,t,axis)
if size(t,1)==1
    t0=repmat(t,size(azx,1),1);
else
    t0=t;
end
f=zeros(size(t0));
for j=1:size(t0,2)
    for i=1:size(t0,1)
        logmat=logm(newU0(w,azx(i,:),t0(i,j)))/(2*acos(0.5*trace(newU0(w,azx(i,:),t0(i,j)))));
        f(i,j)=-2*imag(logmat(1,axis));
    end
end
        
        
