function f = rmerrdata(data)
f=data(2:end-1);
diff2=diff(data,2);
while max(abs(diff2))>0.3 
    f=f.*(abs(diff2)<=0.3)+(abs(diff2)>0.3).*([data(1);f(1:end-1)]+[f(2:end);data(end)])/2;
    diff2=diff([data(1);f;data(end)],2);
end
f(f>1)=1;
f=[data(1);f;data(end)];
end