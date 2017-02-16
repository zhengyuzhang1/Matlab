sz = [1 0;0 -1];
sx = [0 1;1 0];
sy = [0 -1i;1i 0];

w = 2*pi*505*0.00107;
wrf = 2*pi*13.5;
rabi = 2*pi*30e-3;
b0 = [1,0,0];
b0 = b0/norm(b0);

theta = pi/180 * (5)';
alphan = linspace(0,2*pi,10)';

t = ceil(w*trf/2/pi) * 2*pi/w;
tp = (2*pi-mod(wrf*t,2*pi))/(w+wrf);
sequencet = t * [1,1,0] + tp * [0,1,1];
normx = norm(b0-dot(b0,[sin(theta),0,cos(theta)])*[sin(theta),0,cos(theta)]);
trf = pi/2/rabi/normx;
sequencetrf = trf * [1,1,0];
w1 = [wrf,cos(theta),sin(theta)];
infid = zeros(length(alphan),1);
b = zeros(length(alphan),3);
for i = 1:length(alphan)
    alpha = alphan(i);
    [R,~]=AxelRot(alpha/pi*180,[sin(theta),0,cos(theta)]);%rotation angle in Degree
    b(i,:) = b0 * R';    
    sxp = (w1(2)*b(i,1)-w1(3)*b(i,3))*(w1(2)*sx-w1(3)*sz)+b(i,2)*sy;
    [U0,U1] = generalRf(sequencet,sequencetrf,w,wrf,w1,rabi,b(i,:));
    infid(i) = 1-1/4*(trace(abs(expm(1i*pi/2/normx*sxp/2)*U0))+trace(abs(expm(1i*pi/2/normx*sxp/2)*U1)));
end
% figure;surface(infid);
% figure;plot(min(infid,[],2));
% figure;plot(log(sin(thetan)),log(min(infid,[],2)));
%% 
simuinfid = infid;
eps = rabi * normx / (wrf - w);
beta = pi / 4 ./ eps;
% alphan = real(acos((cos(theta)*b(:,1)-sin(theta)*b(:,3))./normx));
infide2 = 1/8*eps.^2.*sin(beta).^2.*(3+cos(2*(beta-alphan-wrf*t)));
infideth = -1/4*eps.*theta.*sin(beta).*sin(w*tp/2).*(cos(beta-w*tp/2)+cos(beta-2*(alphan+wrf*t)-w*tp/2)+2*cos(alpha).*cos(beta-alpha+w*tp/2));
infidth2 = 1/8*theta.^2.*sin(w*tp/2).^2.*(2+cos(2*alphan-w*tp)+cos(2*alphan+w*tp+2*wrf*t));
calcinfid = infide2 + infideth + infidth2;
figure;subplot(1,3,1);plot(alphan,infide2);
xlabel('theta');ylabel('infidelity');title('infidelity ~(rabi/wrf)^2');
subplot(1,3,2);plot(alphan,infideth);
xlabel('theta');ylabel('infidelity');title('infidelity ~(rabi/wrf)*theta');
subplot(1,3,3);plot(alphan,infidth2);
xlabel('theta');ylabel('infidelity');title('infidelity ~theta^2');
figure;plot(alphan,[calcinfid,simuinfid]);
xlabel('theta');ylabel('infidelity');
legend('calculated infidelity to 2nd order','exact infidelity','Location','northwest');
figure;plot(log(alphan),log([calcinfid,simuinfid]));
xlabel('log(theta)');ylabel('log(infidelity)');
legend('calculated infidelity to 2nd order','exact infidelity','Location','northwest');
axis equal