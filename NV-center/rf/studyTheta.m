sz = [1 0;0 -1];
sx = [0 1;1 0];
sy = [0 -1i;1i 0];

w = 2*pi*505*0.00107;
wrf = 2*pi*1.5;
rabi = 2*pi*30e-3;
[R,c]=AxelRot(0,[sin(pi/180*5),0,cos(pi/180*5)]);
b = [1,0,0]*R';
b = b/norm(b);

thetan = pi/180 * (0:20)';
maxtrf = 0;
for i = 1:length(thetan)
    theta = thetan(i);
    trfi = pi/2/rabi/sqrt(1-(dot(b,[sin(theta),0,cos(theta)]))^2);
    if trfi>maxtrf
        maxtrf = trfi;
    end
end

% tpn = linspace(0,5*rabi/wrf*maxtrf,100);
% tpn = linspace(1.5*rabi/wrf*maxtrf,2*rabi/wrf*maxtrf,100);

t = zeros(length(thetan),1);
tp = zeros(length(thetan),1);
trf = zeros(length(thetan),1);
infid = zeros(length(thetan),1);
normx = zeros(length(thetan),1);
for i = 1:length(thetan)
    theta = thetan(i);
    trf(i) = pi/2/rabi/norm(b-dot(b,[sin(theta),0,cos(theta)])*[sin(theta),0,cos(theta)]);
    sequencetrf = trf(i) * [1,1,0];
    w1 = [wrf,cos(theta),sin(theta)];
    sxp = (w1(2)*b(1)-w1(3)*b(3))*(w1(2)*sx-w1(3)*sz)+b(2)*sy;
    normx(i) = sqrt(abs(sxp(1))^2+abs(sxp(2))^2);
    t(i) = ceil(w*trf(i)/2/pi) * 2*pi/w;
    tp(i) = (2*pi-mod(wrf*t(i),2*pi))/(w+wrf);
    sequencet = t(i) * [1,1,0] + tp(i) * [0,1,1];
    [U0,U1] = generalRf(sequencet,sequencetrf,zeros(3,1),w,wrf,w1,rabi,b);
    infid(i) = 1-1/4*(trace(abs(expm(1i*pi/2/normx(i)*sxp/2)*U0))+trace(abs(expm(1i*pi/2/normx(i)*sxp/2)*U1)));
end
% figure;surface(infid);
% figure;plot(min(infid,[],2));
% figure;plot(log(sin(thetan)),log(min(infid,[],2)));
%% 
simuinfid = infid;
eps = rabi * normx / (wrf - w);
beta = pi / 4 ./ eps;
alpha = real(acos((cos(thetan)*b(1)-sin(thetan)*b(3))./normx));
infide2 = 1/8*eps.^2.*sin(beta).^2.*(3+cos(2*(beta-alpha-wrf*t)));
infideth = -1/4*eps.*thetan.*sin(beta).*sin(w*tp/2).*(cos(beta-w*tp/2)+cos(beta-2*(alpha+wrf*t)-w*tp/2)+2*cos(alpha).*cos(beta-alpha+w*tp/2));
infidth2 = 1/8*thetan.^2.*sin(w*tp/2).^2.*(2+cos(2*alpha-w*tp)+cos(2*alpha+w*tp+2*wrf*t));
calcinfid = infide2 + infideth + infidth2;
figure;subplot(1,3,1);plot(thetan,infide2);
xlabel('theta');ylabel('infidelity');title('infidelity ~(rabi/wrf)^2');
subplot(1,3,2);plot(thetan,infideth);
xlabel('theta');ylabel('infidelity');title('infidelity ~(rabi/wrf)*theta');
subplot(1,3,3);plot(thetan,infidth2);
xlabel('theta');ylabel('infidelity');title('infidelity ~theta^2');
figure;plot(thetan,[calcinfid,simuinfid]);
xlabel('theta');ylabel('infidelity');
legend('calculated infidelity to 2nd order','exact infidelity','Location','northwest');
figure;plot(log(thetan),log([calcinfid,simuinfid]));
xlabel('log(theta)');ylabel('log(infidelity)');
legend('calculated infidelity to 2nd order','exact infidelity','Location','northwest');
axis equal