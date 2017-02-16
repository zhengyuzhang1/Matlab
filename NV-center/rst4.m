function f=rst4(w,azx,n)
w1=repmat(sqrt((azx(:,1)+w).^2+azx(:,2).^2),1,n);
w=repmat(w,size(azx,1),n);
t=repmat((1:n),size(azx,1),1)*pi./(w+w1);
b=repmat(azx(:,2),1,n);
f=t...
    +sin(w.*t)/2./(w+w1).*(b./w1).^2-sin(w.*t)/8.*(2*w1+2*w-(w1-w).*cos(w.*t)).*(b./w1).^4./(w+w1).^2;
end
