function [pulse,phase,fid]=optimCnotn(w,azx,nop,tcap,tol,rst)

non=size(azx,1);

ubn=floor(tcap./(4*rst)); %absolute upper-bound for each resonant time
phi=rotang(w,azx,rst);
pulse=zeros(non,4);
phase=zeros(non,non,2);
fid=zeros(non,3);
for i=1:non
    for j=1:2:nop
        for l=2:2:nop
            n1=1:ubn(i,j);
            n1=n1(abs(abs(mod(phi(i,j)*n1./(0.5*pi),4)-2)-1)<tol);
            n2=1:ubn(i,l);
            n2=n2(abs(abs(mod(phi(i,l)*n2./(0.5*pi),4)-2)-1)<tol);
            for in1=n1
                for in2=n2
                    if 4*(rst(i,j)*in1+rst(i,l)*in2)>tcap
                        continue;
                    end
                    totalU0=U0(w,azx(i,:),rst(i,l))^in2*U0(w,azx(i,:),rst(i,j))^in1;
                    tmp1=abs(totalU0(1));
                    totalU1=U1(w,azx(i,:),rst(i,l))^in2*U1(w,azx(i,:),rst(i,j))^in1;
                    tmp2=abs(totalU1(2));
                    for k=[1:i-1,i+1:non]
                        tmp=U0(w,azx(k,:),rst(i,l))^in2*U0(w,azx(k,:),rst(i,j))^in1;
                        tmp1=tmp1*abs(tmp(1));
                        tmp=U1(w,azx(k,:),rst(i,l))^in2*U1(w,azx(k,:),rst(i,j))^in1;
                        tmp2=tmp2*abs(tmp(1));
                    end
                    if 0.5*(tmp1+tmp2)>fid(i,3)
                        pulse(i,:)=[rst(i,j),in1,rst(i,l),in2];
                        fid(i,:)=[tmp1,tmp2,0.5*(tmp1+tmp2)];
                    end
                end 
            end
        end
    end
end
for j=1:non
    for i=1:non
        tmp=(U0(w,azx(i,:),pulse(j,3))^pulse(j,4))*(U0(w,azx(i,:),pulse(j,1))^pulse(j,2));
        phase(i,j,1)=tmp(1)/abs(tmp(1));
        tmp=(U1(w,azx(i,:),pulse(j,3))^pulse(j,4))*(U1(w,azx(i,:),pulse(j,1))^pulse(j,2));
        phase(i,j,2)=tmp(2)/abs(tmp(2));
    end
end  
end