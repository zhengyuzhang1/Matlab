function f = oldrmerrdata(data)
f=data;
for i=2:length(f)-1
    if abs(f(i+1)-f(i))>0.1 && abs(f(i-1)-f(i))>0.1 && abs(f(i+1)+f(i-1)-2*f(i))>0.3 
        f(i)=(f(i+1)+f(i-1))/2;
    end
end
f(f>1)=1;
end