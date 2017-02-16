function [pulse,phase,fid]=optimCnot_short(w,azx,nop)
t=((nop+0.5)*pi./(2*w+azx(:,1)))*(0:1e-4:1);
n0x=axes0(w,azx,t,2);
non=size(azx,1);
rst=zeros(non,nop);
sgn=zeros(size(rst));
for i=1:size(azx,1)
    [~,locsp]=findpeaks(n0x(i,:),'MINPEAKHEIGHT',0);
    [~,locsd]=findpeaks(-n0x(i,:),'MINPEAKHEIGHT',0);
    while length(locsp)+length(locsd)>nop
        if locsp(end)>locsd(end)
            locsp(end)=[];
        else locsd(end)=[];
        end
    end
    funp=@(t) -axes0(w,azx(i,:),t,2);
    fund=@(t) axes0(w,azx(i,:),t,2);
    for j=1:length(locsp)
        rst(i,j)=fminbnd(funp,t(i,locsp(j)-1),t(i,locsp(j)+1),optimset('TolX',1e-12));   
    end
    for j=1:length(locsd)
        rst(i,length(locsp)+j)=fminbnd(fund,t(i,locsd(j)-1),t(i,locsd(j)+1),optimset('TolX',1e-12));
    end
    [rst(i,:),idx]=sort(rst(i,:));
    sgn(i,:)=[ones(1,length(locsp)),-ones(1,length(locsd))];
    sgn(i,:)=sgn(i,idx);
end
phi=rotang(w,azx,rst);
n=ones(non,nop,2);
for j=1:nop
    for i=1:non
        tmp=2*pi;
        while abs(mod(n(i,j,1)*phi(i,j),2*pi)-pi/2)<tmp
            tmp=abs(mod(n(i,j,1)*phi(i,j),2*pi)-pi/2);
            n(i,j,1)=n(i,j,1)+1;
        end
        tmp=2*pi;
        while abs(mod(n(i,j,2)*phi(i,j),2*pi)-3*pi/2)<tmp
            tmp=abs(mod(n(i,j,2)*phi(i,j),2*pi)-3*pi/2);
            n(i,j,2)=n(i,j,2)+1;
        end
    end
end
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