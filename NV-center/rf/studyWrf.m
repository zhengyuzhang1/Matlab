sz = [1 0;0 -1];
sx = [0 1;1 0];
sy = [0 -1i;1i 0];

w = 2*pi*505*0.00107;
rabi = 2*pi*30e-3;
b = [1,1,1];
b = b/norm(b);
theta = pi/180 * (20)';

wrfn = 2*pi*(1:0.5:15)';

t = ceil(w*trf/2/pi) * 2*pi/w;
normx = norm(b-dot(b,[sin(theta),0,cos(theta)])*[sin(theta),0,cos(theta)]);
trf = pi/2/rabi/normx;
sequencetrf = trf * [1,1,0];
sxp = (w1(2)*b(1)-w1(3)*b(3))*(w1(2)*sx-w1(3)*sz)+b(2)*sy; 
infid = zeros(length(wrfn),1);
for i = 1:length(wrfn)
    wrf = wrfn(i);
    w1 = [wrf,cos(theta),sin(theta)];
    tp = (2*pi-mod(wrf*t,2*pi))/(w+wrf);
    sequencet = t * [1,1,0] + tp * [0,1,1];
    [U0,U1] = generalRf(sequencet,sequencetrf,w,wrf,w1,rabi,b);
    infid(i) = 1-1/4*(trace(abs(expm(1i*pi/2/normx*sxp/2)*U0))+trace(abs(expm(1i*pi/2/normx*sxp/2)*U1)));
end
% figure;surface(infid);
% figure;plot(min(infid,[],2));
% figure;plot(log(sin(thetan)),log(min(infid,[],2)));
%% 
simuinfid = infid;
eps = rabi .* normx ./ (wrfn - w);
beta = pi / 4 ./ eps;
alpha = real(acos((cos(theta)*b(1)-sin(theta)*b(3))./normx));
infide2 = 1/8*eps.^2.*sin(beta).^2.*(3+cos(2*(beta-alpha-wrfn*t)));
infideth = -1/4*eps.*theta.*sin(beta).*sin(w*tp/2).*(cos(beta-w*tp/2)+cos(beta-2*(alpha+wrfn*t)-w*tp/2)+2*cos(alpha).*cos(beta-alpha+w*tp/2));
infidth2 = 1/8*theta.^2.*sin(w*tp/2).^2.*(2+cos(2*alpha-w*tp)+cos(2*alpha+w*tp+2*wrfn*t));
calcinfid = infide2 + infideth + infidth2;
figure;subplot(1,3,1);plot(eps,infide2);
xlabel('eps');ylabel('infidelity');title('infidelity ~(rabi/wrf)^2');
subplot(1,3,2);plot(eps,infideth);
xlabel('eps');ylabel('infidelity');title('infidelity ~(rabi/wrf)*theta');
subplot(1,3,3);plot(eps,infidth2);
xlabel('eps');ylabel('infidelity');title('infidelity ~theta^2');
figure;plot(eps,[calcinfid,simuinfid]);
xlabel('eps');ylabel('infidelity');
legend('calculated infidelity to 2nd order','exact infidelity','Location','northwest');
figure;plot(log(eps),log([calcinfid,simuinfid]));
xlabel('log(eps)');ylabel('log(infidelity)');
axis equal
legend('calculated infidelity to 2nd order','exact infidelity','Location','northwest');
axis equal