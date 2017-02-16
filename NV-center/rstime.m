function [rst, sgn] = rstime(w, azx, nop)
%resonant time for NV-center, nop is the max number of period

t=((nop+0.5)*pi./(2*w+azx(:,1)))*(0:1e-4:1);
n0x=axes0(w,azx,t,2);
non=size(azx,1); %number of nuclei
rst=zeros(non,nop); %resonant time
sgn=zeros(size(rst)); %rotating axis of each resonant: +1 for x, -1 for -x 
for i=1:non
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