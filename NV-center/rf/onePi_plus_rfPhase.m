w = 2*pi*500*0.00107;
wrf = 2*pi*13.7;
theta = pi/180 * 7;
w1 = [wrf,cos(theta),sin(theta)];
rabi = 2*pi*40e-3;
b = [1,0,0];
b = b/norm(b);

t = 7;
trf = pi/2/rabi/norm(b-dot(b,[sin(theta),0,cos(theta)])*[sin(theta),0,cos(theta)]);
phi = (w - wrf)*t;

a = zeros(10,1);
U0 = zeros(2,2,10);
U1 = zeros(2,2,10);
U0(:,:,1) = eye(2);
U1(:,:,1) = eye(2);
for i = 2:10
    [tempU0, tempU1] = generalRf([t t], [trf trf], phi*[i-2 i-1]+wrf*t*[2*(i-2) 2*(i-2)+1], w, wrf, w1, rabi, b);
    if mod(i,2) == 0
        U0(:,:,i) = tempU0 * U0(:,:,i-1);
        U1(:,:,i) = tempU1 * U1(:,:,i-1);
    else
        U0(:,:,i) = tempU1 * U0(:,:,i-1);
        U1(:,:,i) = tempU0 * U1(:,:,i-1);
    end
end
%lab frame to rf frame
for i = 2:10
    Ur = rSpin(-wrf*2*t*(i-1)*[sin(theta),0,cos(theta)]); 
    U0(:,:,i) = Ur * U0(:,:,i) * Ur';
    U1(:,:,i) = Ur * U1(:,:,i) * Ur';
end
figure;plot(reshape(abs(U0(1,1,:)).^2,10,1));
%%
tpmin = (2*pi-mod(wrf*t,2*pi))./(w+wrf);
alpha = real(acos((cos(thetan)*b(1)-sin(thetan)*b(3))./normx));
simuinfid = min(infid,[],2);
simuinfid = simuinfid - simuinfid(1);
calcinfid = 1/8*thetan.^2.*sin(w*tpmin/2).^2.*(2+cos(2*alpha-w*tpmin)+cos(2*alpha+w*tpmin+2*wrf*t));
figure;plot(thetan,[calcinfid,simuinfid]);
figure;plot(log(thetan),log([calcinfid,simuinfid]));
