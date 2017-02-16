sz = [1 0;0 -1];
sx = [0 1;1 0];
sy = [0 -1i;1i 0];

w = 2*pi*505*0.00107;
wrf = 2*pi*13.7;
rabi = 2*pi*30e-3;
b = [0,1,0];
b = b/norm(b);

thetan = pi/180 * (0:0.5:10)';
tpn = (0:10);

trf = zeros(length(thetan),1);
infid = zeros(length(thetan),length(tpn));
normx = zeros(length(thetan),1);
for i = 1:length(thetan)
    theta = thetan(i);
    w1 = [wrf,cos(theta),sin(theta)];
    sxp = (w1(2)*b(1)-w1(3)*b(3))*(w1(2)*sx-w1(3)*sz)+b(2)*sy;
    normx(i) = sqrt(abs(sxp(1))^2+abs(sxp(2))^2);
    trf(i) = pi/2/rabi/normx(i);
    sequencetrf = trf(i) * [1,1];
    for j =1:length(tpn)
        t = trf(i) + tpn(j);
        sequencet = t(i) * [1,1];
        sequencephi = [0, -mod(w*t,2*pi)];
        [U0,U1] = generalRf(sequencet,sequencetrf,w,wrf,w1,rabi,b);
        infid(i,j) = 1-1/4*(trace(abs(expm(1i*pi/2/normx(i)*sxp/2)*U0))+trace(abs(expm(1i*pi/2/normx(i)*sxp/2)*U1)));
    end
end
figure;surface(infid);
%%
tpmin = (2*pi-mod(wrf*t,2*pi))./(w+wrf);
alpha = real(acos((cos(thetan)*b(1)-sin(thetan)*b(3))./normx));
simuinfid = min(infid,[],2);
simuinfid = simuinfid - simuinfid(1);
calcinfid = 1/8*thetan.^2.*sin(w*tpmin/2).^2.*(2+cos(2*alpha-w*tpmin)+cos(2*alpha+w*tpmin+2*wrf*t));
figure;plot(thetan,[calcinfid,simuinfid]);
figure;plot(log(thetan),log([calcinfid,simuinfid]));
