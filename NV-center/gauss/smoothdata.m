%smooth data in delta
sop =  zeros(size(sopU));
for i = 1:size(sop,1)
    for j = 1:size(sop,2)
        for k = 1:size(sop,3)
            tmp = sopU{i,j,k};
            sop(i,j,k) = norm(tmp(2))^2;
        end
    end
end
sop2 = zeros(size(sop));
for i = 1:length(t)
    for j = 1:size(sigma,2)
        cc = reshape(sop(i,j,:),1,length(delta));
        [pv,locs] = findpeaks(cc);
        err = pv > 0.5;
        if isempty(pv) || all(err) 
            continue;
        else
            locs(err) = [];
            pv(err) = [];
        end
        err = diff(pv) > 0;
        pv(err) = [];
        locs(err) = [];
        err = diff(pv,2) < 0;
        pv(err) = [];
        locs(err) = [];
        if length(pv) == 1
            cc(cc - pv < 0) = pv;
            sop2(i,j,:)=cc;
            continue;
        end
        vq = interp1(delta(locs),pv,delta,'spline');
        l = 1;
        while vq(l) < sop(i,j,l)
            l = l + 1;
        end
        locs = [l-1,l,locs];
        pv = [sop(i,j,l-1),vq(l),pv];
        vq2 = interp1(delta(locs),pv,delta(l):delta(end),'pchip');
        sop2(i,j,:) = [reshape(sop(i,j,1:l-1),1,l-1),vq2];
    end
end
%%
%smooth data in sigma
sop3 = zeros(size(sop2));
for i = 1:length(t)
    for j = 1:length(delta)
        tmp = reshape(sop2(i,:,j),1,size(sigma,2));
        tmp = smooth(tmp);
        sop3(i,:,j) = tmp;
    end
end
%%
 deln = 101;
 tn = 16:20;
 [popt,sigmaopt] = min(sop3(tn,:,deln),[],2);
 sigmaT = sigma(tn,:)';
 figure;plot(t(tn)*delta(deln),popt,'-o');
 xlabel('T*\Delta');
 ylabel('minP');
 figure;plot(t(tn)*delta(deln),log(popt),'-o');
 xlabel('T*\Delta');
 ylabel('log(minP)');
 figure;plot(t(tn)*delta(deln),sigmaT(size(sigmaT,1)*(0:(size(sigmaT,2)-1))+sigmaopt')/delta(deln),'-o');
 xlabel('T*\Delta');
 ylabel('\deltam/\Delta');
 %%
 deln1 = 101;
 popt1 = zeros(length(t),1);
 sigmaopt1 = zeros(size(popt1));
 for i = 1:length(t)
     cc = sop3(i,:,deln1);
     [tmp1,tmp2] = findpeaks(-cc);
     popt1(i) = -tmp1(1);
     sigmaopt1(i) = tmp2(1);
 end
 sigmaT = sigma';
 figure;plot(t*delta(deln1),popt1,'-o');
 xlabel('T*\Delta');
 ylabel('minP');
 set(findall(gcf,'-property','FontSize'),'FontSize',14)
 figure;plot(t*delta(deln1),log(popt1),'-o');
 xlabel('T*\Delta');
 ylabel('log(minP)');
 set(findall(gcf,'-property','FontSize'),'FontSize',14)
 figure;plot(t*delta(deln1),sigmaT(size(sigma,2)*(0:(size(sigma,1)-1))+sigmaopt1')/delta(deln1),'-o');
 xlabel('T*\Delta');
 ylabel('\deltam/\Delta');
 set(findall(gcf,'-property','FontSize'),'FontSize',14)
 %%
 i = 5:8;
 j = 201;
 k = 30:75;
 figure;plot(sigma(i,k)',sop3(i,k,j)','-');
 %%
 i = 9;
 j = 18;
 k = 60:501;
 figure;plot(delta(k),reshape(sop2(i,j,k),length(j),length(k))');
 figure;plot(delta(k),reshape(sop(i,j,k),length(j),length(k))');
 %%
 i = 11;
 j = 6;
 cc = reshape(sop(i,j,:),1,length(delta));
 [pv,locs] = findpeaks(cc);
 err = pv > 0.5;
 locs(err) = [];
 pv(err) = [];
 err = diff(pv) > 0;
 pv(err) = [];
 locs(err) = [];
 err = diff(pv,2) < 0;
 pv(err) = [];
 locs(err) = [];
 locs = locs(2:end);
 pv = pv(2:end);
 fun1 = @(x,xdata) x(1) * exp(-x(2) * xdata.^2) + x(3);
 fun2 = @(x,xdata) x(1) * exp(-x(2) * xdata) + x(3);
 x01 = [1 0.001 0];
 x02 = [1 0.01 0];
 xdata = locs * t(i);
 ydata = pv;
 x1 = lsqcurvefit(fun1, x01, xdata, ydata);
 x2 = lsqcurvefit(fun2, x0, xdata, ydata);
 figure;plot(t(i)*delta(min(locs)-10:max(locs)+1),cc(min(locs)-10:max(locs)+1),xdata,fun1(x1,xdata),xdata,fun2(x2,xdata))
 xlabel('T*\Delta');
 ylabel('P');
 legend(sprintf('P(T=%0.2f,\\delta=%0.f)',t(i),sigma(i,j)),sprintf('%.2g*exp(-%.2g*x^2)+%.2g',x1(1),x1(2),x1(3)),sprintf('%.2g*exp(-%.2g*x)+%.2g',x2(1),x2(2),x2(3)));
 set(findall(gcf,'-property','FontSize'),'FontSize',14)
 figure;plot(xdata,log(ydata),'-o');
 xlabel('T*\Delta');
 ylabel('log(maxima of P)');
 legend(sprintf('log(maxima of P(T=%0.2f,\\delta=%0.f))',t(i),sigma(i,j)));
 set(findall(gcf,'-property','FontSize'),'FontSize',14)
 figure;plot(t(i)*delta(1:250),[cc(1:250);reshape(sop(i,j+5,1:250),1,250);reshape(sop(i,j+10,1:250),1,250)]);
 xlabel('T*\Delta');
 ylabel('P');
 legend(sprintf('P(T=%0.2f,\\delta=%0.f)',t(i),sigma(i,j)),sprintf('P(T=%0.2f,\\delta=%0.f)',t(i),sigma(i,j+5)),sprintf('P(T=%0.2f,\\delta=%0.f)',t(i),sigma(i,j+10)));
 set(findall(gcf,'-property','FontSize'),'FontSize',14)
 %%
 xdata = cell(size(sigma));
 ydata = cell(size(sigma));
 for i = 1:size(sigma,1)
     for j = 1:size(sigma,2)
         cc = reshape(sop(i,j,:),1,length(delta));
         [pv,locs] = findpeaks(cc);
         err = pv > 0.5;
         locs(err) = [];
         pv(err) = [];
         err = diff(pv) > 0;
         pv(err) = [];
         locs(err) = [];
         err = diff(pv,2) < 0;
         pv(err) = [];
         locs(err) = [];
         locs = locs(2:end);
         pv = pv(2:end);
         xdata{i,j} = locs * t(i);
         ydata{i,j} = pv;
     end
 end
 
 %%
 figure;
 i = 42;
 j = 81;
 hold on
 for i = 42:size(sigma,1)
     plot(xdata{i,81},log(ydata{i,81}),'-o');
 end
 hold off
 xlabel('T*\Delta');
 ylabel('log(maxima of P)');
 legend(sprintf('log(maxima of P(T=%0.2f,\\delta=%0.f))',t(i),sigma(i,j)));
 set(findall(gcf,'-property','FontSize'),'FontSize',14)