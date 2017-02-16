function f=fidelity(w,azx,env,pulset,pulsen)
tmp=U0(w,azx,pulset(2))^pulsen(2)*U0(w,azx,pulset(1))^pulsen(1);
tmp1=abs(tmp(1));
tmp=U1(w,azx,pulset(2))^pulsen(2)*U1(w,azx,pulset(1))^pulsen(1);
tmp2=abs(tmp(2));
for k=1:size(env,1)
    tmp=U0(w,env(k,:),pulset(2))^pulsen(2)*U0(w,env(k,:),pulset(1))^pulsen(1);
    tmp1=tmp1*abs(tmp(1));
    tmp=U1(w,env(k,:),pulset(2))^pulsen(2)*U1(w,env(k,:),pulset(1))^pulsen(1);
    tmp1=tmp1*abs(tmp(1));
end
f=0.5*(tmp1+tmp2);
end
