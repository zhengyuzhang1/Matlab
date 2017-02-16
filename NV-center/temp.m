sx = [0 1; 1 0];
sy = [0 -1i; 1i 0];
sz = [1 0; 0 -1];
r = rand(2,1);
theta = rand(2,1) * pi;
phi = rand(2,1) * 2*pi;
p = r(1) * [sin(theta) .* cos(phi), sin(theta) .* sin(phi), cos(theta)];
rou1 = 0.5 * (eye(2) + p(1,1) * sx + p(1,2) * sy + p(1,3) * sz);
rou2 = 0.5 * (eye(2) + p(2,1) * sx + p(2,2) * sy + p(2,3) * sz);
x1 = sqrtm(sqrtm(rou1) * rou2 * sqrtm(rou1));
x2 = sqrtm(rou1) * sqrtm(rou2);