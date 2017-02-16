function [pulse,Ueachnuc,totalU]=optimpulse(w,azx,gate,tcap,rst,phi)
%nop is the max number of period, tcap the max total time

phi = rotang(w, azx, rst); %rotation angle
ubn1=floor(tcap./(4*rst));
ubn=repmat(permute(1:max(ubn1(:)),[3,1,2]),non,nop,1);
ubn(ubn>repmat(ubn1,1,1,max(ubn1(:))))=NaN;
n=zeros(non,nop,2);
[~,n(:,:,1)]=min(abs(mod(bsxfun(@times,ubn,phi),2*pi)-pi/2),[],3);
[~,n(:,:,2)]=min(abs(mod(bsxfun(@times,ubn,phi),2*pi)-3*pi/2),[],3);
pulse=zeros(non,4);
phase=zeros(non,non,2);
fid=zeros(non,3);
for i=1:non
    for j=1:2:nop
        for l=2:2:nop
            if sgn(i,j)==sgn(i,l)
                tmp=U0(w,azx(i,:),rst(i,l))^n(i,l,2)*U0(w,azx(i,:),rst(i,j))^n(i,j,1);
                tmp1=abs(tmp(1));
                tmp=U1(w,azx(i,:),rst(i,l))^n(i,l,2)*U1(w,azx(i,:),rst(i,j))^n(i,j,1);
                tmp2=abs(tmp(2));
                for k=[1:i-1,i+1:non]
                    tmp=U0(w,azx(k,:),rst(i,l))^n(i,l,2)*U0(w,azx(k,:),rst(i,j))^n(i,j,1);
                    tmp1=tmp1*abs(tmp(1));
                    tmp=U1(w,azx(k,:),rst(i,l))^n(i,l,2)*U1(w,azx(k,:),rst(i,j))^n(i,j,1);
                    tmp2=tmp2*abs(tmp(1));
                end
                a=[tmp1,tmp2,0.5*(tmp1+tmp2)];
                tmp=U0(w,azx(i,:),rst(i,l))^n(i,l,1)*U0(w,azx(i,:),rst(i,j))^n(i,j,2);
                tmp1=abs(tmp(1));
                tmp=U1(w,azx(i,:),rst(i,l))^n(i,l,1)*U1(w,azx(i,:),rst(i,j))^n(i,j,2);
                tmp2=abs(tmp(2));
                for k=[1:i-1,i+1:non]
                    tmp=U0(w,azx(k,:),rst(i,l))^n(i,l,1)*U0(w,azx(k,:),rst(i,j))^n(i,j,2);
                    tmp1=tmp1*abs(tmp(1));
                    tmp=U1(w,azx(k,:),rst(i,l))^n(i,l,1)*U1(w,azx(k,:),rst(i,j))^n(i,j,2);
                    tmp2=tmp2*abs(tmp(1));
                end
                b=[tmp1,tmp2,0.5*(tmp1+tmp2)];
                if a(3)>b(3) && a(3)>fid(i,3)
                    pulse(i,:)=[rst(i,j),n(i,j,1),rst(i,l),n(i,l,2)];
                    fid(i,:)=a;
                elseif b(3)>fid(i,3)
                    pulse(i,:)=[rst(i,j),n(i,j,2),rst(i,l),n(i,l,1)];
                    fid(i,:)=b;
                end
            else
                tmp=U0(w,azx(i,:),rst(i,l))^n(i,l,1)*U0(w,azx(i,:),rst(i,j))^n(i,j,1);
                tmp1=abs(tmp(1));
                tmp=U1(w,azx(i,:),rst(i,l))^n(i,l,1)*U1(w,azx(i,:),rst(i,j))^n(i,j,1);
                tmp2=abs(tmp(2));
                for k=[1:i-1,i+1:non]
                    tmp=U0(w,azx(k,:),rst(i,l))^n(i,l,1)*U0(w,azx(k,:),rst(i,j))^n(i,j,1);
                    tmp1=tmp1*abs(tmp(1));
                    tmp=U1(w,azx(k,:),rst(i,l))^n(i,l,1)*U1(w,azx(k,:),rst(i,j))^n(i,j,1);
                    tmp2=tmp2*abs(tmp(1));
                end
                a=[tmp1,tmp2,0.5*(tmp1+tmp2)];
                tmp=U0(w,azx(i,:),rst(i,l))^n(i,l,2)*U0(w,azx(i,:),rst(i,j))^n(i,j,2);
                tmp1=abs(tmp(1));
                tmp=U1(w,azx(i,:),rst(i,l))^n(i,l,2)*U1(w,azx(i,:),rst(i,j))^n(i,j,2);
                tmp2=abs(tmp(2));
                for k=[1:i-1,i+1:non]
                    tmp=U0(w,azx(k,:),rst(i,l))^n(i,l,2)*U0(w,azx(k,:),rst(i,j))^n(i,j,2);
                    tmp1=tmp1*abs(tmp(1));
                    tmp=U1(w,azx(k,:),rst(i,l))^n(i,l,2)*U1(w,azx(k,:),rst(i,j))^n(i,j,2);
                    tmp2=tmp2*abs(tmp(1));
                end
                b=[tmp1,tmp2,0.5*(tmp1+tmp2)];
                if a(3)>b(3) && a(3)>fid(i,3)
                    pulse(i,:)=[rst(i,j),n(i,j,1),rst(i,l),n(i,l,1)];
                    fid(i,:)=a;
                elseif b(3)>fid(i,3)
                    pulse(i,:)=[rst(i,j),n(i,j,2),rst(i,l),n(i,l,2)];
                    fid(i,:)=b;
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