function [pulse,phase,fid]=optimCnot(w,azx,nop)
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
ubn1=floor(1e-3./(4*rst));
ubn=repmat(permute(1:max(ubn1(:)),[3,1,2]),non,nop,1);
ubn(ubn>repmat(ubn1,1,1,max(ubn1(:))))=NaN;
rsd=zeros(non,nop,2);
n=zeros(non,nop,2);
[rsd(:,:,1),n(:,:,1)]=min(abs(mod(bsxfun(@times,ubn,phi),2*pi)-pi/2),[],3);
[rsd(:,:,2),n(:,:,2)]=min(abs(mod(bsxfun(@times,ubn,phi),2*pi)-3*pi/2),[],3);
totrotang=bsxfun(@times,rotang(w,azx,rst(:).'),reshape(n,1,non*nop,2));
n0xerr=reshape(sum(realpow(bsxfun(@times,axes0(w,azx,rst(:).',2),sin(totrotang/2)),2)),non,nop,2);
pulse=zeros(non,4);
phase=zeros(non,non,2);
fid=ones(non,6);
for i=1:non
    tmpminrsd=1000;
    for j=1:2:nop
        for l=2:2:nop
            if sgn(i,j)==sgn(i,l)
                a=rsd(i,j,1)+rsd(i,l,2)+n0xerr(i,j,1)+n0xerr(i,l,2);
                b=rsd(i,j,2)+rsd(i,l,1)+n0xerr(i,j,2)+n0xerr(i,l,1);
                if a<b && a<tmpminrsd
                    pulse(i,:)=[rst(i,j),n(i,j,1),rst(i,l),n(i,l,2)];
                elseif b<tmpminrsd
                    pulse(i,:)=[rst(i,j),n(i,j,2),rst(i,l),n(i,l,1)];
                end
            else
                a=rsd(i,j,1)+rsd(i,l,1)+n0xerr(i,j,1)+n0xerr(i,l,1);
                b=rsd(i,j,2)+rsd(i,l,2)+n0xerr(i,j,2)+n0xerr(i,l,2);
                if a<b && a<tmpminrsd
                    pulse(i,:)=[rst(i,j),n(i,j,1),rst(i,l),n(i,l,1)];
                elseif b<tmpminrsd
                    pulse(i,:)=[rst(i,j),n(i,j,2),rst(i,l),n(i,l,2)];
                end
            end
            tmpminrsd=min([tmpminrsd,a,b]);
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
for i=1:non
    tmp=U0(w,azx(i,:),pulse(i,3))^pulse(i,4)*U0(w,azx(i,:),pulse(i,1))^pulse(i,2);
    fid(i,1)=abs(tmp(1));
    tmp=U1(w,azx(i,:),pulse(i,3))^pulse(i,4)*U1(w,azx(i,:),pulse(i,1))^pulse(i,2);
    fid(i,2)=abs(tmp(2));
    fid(i,3)=(fid(i,1)+fid(i,2))/2;
    fid(i,4)=fid(i,1);
    fid(i,5)=fid(i,2);
    for j=[1:i-1,i+1:non]
        tmp=U0(w,azx(j,:),pulse(i,3))^pulse(i,4)*U0(w,azx(j,:),pulse(i,1))^pulse(i,2);
        fid(i,4)=fid(i,4)*abs(tmp(1));
        tmp=U1(w,azx(j,:),pulse(i,3))^pulse(i,4)*U1(w,azx(j,:),pulse(i,1))^pulse(i,2);
        fid(i,5)=fid(i,5)*abs(tmp(1));
    end
    fid(i,6)=(fid(i,4)+fid(i,5))/2;
end    
end