function f=rotang(w,azx,t)
w0=repmat(w,size(azx,1),1);
w1=realsqrt(realpow(azx(:,1)+w0,2)+realpow(azx(:,2),2));
mz=(azx(:,1)+w0)./w1;
f=2*acos(cos(bsxfun(@times,w1,t)).*cos(bsxfun(@times,w0,t))-bsxfun(@times,mz,sin(bsxfun(@times,w1,t)).*sin(bsxfun(@times,w0,t))));
end