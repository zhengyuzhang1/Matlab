function [gt,phi] = gatetau(w, azx, nop)
%gate time tau and rotation angle for NV-center, nop is the max number of period;
%gt and phi are non*nop*3 matrices, 3rd dimension correspond to z-rotation,
%unconditional x-rotation and conditional x-rotation, respectively.
%gt gives inter-pulse time tau and phi gives rotation angle

t = ((2 * nop + 0.5) * pi ./ (2 * w + azx(:,1))) * (0:1e-4:1);
n0x = axes0(w,azx,t,2);
n0z = axes0(w,azx,t,1);
non = size(azx,1); %number of nuclei
gt = zeros(non,nop,3); 
phi = zeros(non,nop,3);
 
for i = 1:non
    %find x-rotation tau
    [~,locsp] = findpeaks(n0x(i,:),'MINPEAKHEIGHT',0);
    [~,locsd] = findpeaks(-n0x(i,:),'MINPEAKHEIGHT',0);
    while length(locsp) + length(locsd) > 2 * nop
        if locsp(end) > locsd(end)
            locsp(end) = [];
        else locsd(end) = [];
        end
    end
    funp = @(t) -axes0(w,azx(i,:),t,2);
    fund = @(t) axes0(w,azx(i,:),t,2);
    xtau = zeros(1,2 * nop);
    for j = 1:length(locsp)
        xtau(j) = fminbnd(funp,t(i,locsp(j)-1),t(i,locsp(j)+1),optimset('TolX',1e-12));   
    end
    for j = 1:length(locsd)
        xtau(length(locsp)+j) = fminbnd(fund,t(i,locsd(j)-1),t(i,locsd(j)+1),optimset('TolX',1e-12));
    end
    [xtau,idx] = sort(xtau);
    xsgn = [ones(1,length(locsp)),-ones(1,length(locsd))];
    xsgn = xsgn(idx);
    gt(i,:,2) = xtau(2:2:end);
    phi(i,:,2) = rotang(w,azx(i,:),gt(i,:,2)) .* xsgn(2:2:end);
    gt(i,:,3) = xtau(1:2:end);
    phi(i,:,3) = rotang(w,azx(i,:),gt(i,:,3)) .* xsgn(1:2:end);
    %find z rotation tau
    [~,locsp] = findpeaks(n0z(i,:),'MINPEAKHEIGHT',0);
    [~,locsd] = findpeaks(-n0z(i,:),'MINPEAKHEIGHT',0);
    while length(locsp) + length(locsd) > nop
        if locsp(end) > locsd(end)
            locsp(end) = [];
        else locsd(end) = [];
        end
    end
    funp = @(t) -axes0(w,azx(i,:),t,1);
    fund = @(t) axes0(w,azx(i,:),t,1);
    ztau = zeros(1,nop);
    for j = 1:length(locsp)
        ztau(j) = fminbnd(funp,t(i,locsp(j)-1),t(i,locsp(j)+1),optimset('TolX',1e-12));   
    end
    for j = 1:length(locsd)
        ztau(length(locsp)+j) = fminbnd(fund,t(i,locsd(j)-1),t(i,locsd(j)+1),optimset('TolX',1e-12));
    end
    [ztau,idx] = sort(ztau);
    zsgn = [ones(1,length(locsp)),-ones(1,length(locsd))];
    zsgn = zsgn(idx);
    gt(i,:,1) = ztau;
    phi(i,:,1) = rotang(w,azx(i,:),gt(i,:,1)) .* zsgn;
end