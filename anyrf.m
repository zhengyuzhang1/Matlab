function U = anyrf(t, w0, w1, w, phi, n0, n1)
%H0 = w0*(n0*sigma/2), H1(tao) = w1*sin(w*tao+phi)*(n1*sigma), U = U(t)
%under rotating wave approximation in lab frame

%normalization of n0,n1
n0 = n0 / norm(n0);
n1 = n1 / norm(n1);

sx = [0 1;1 0];
sy = [0 -1i;1i 0];
sz = [1 0;0 -1];


nx1 = n1-dot(n0,n1)*n0;
sintheta = norm(nx1);
nx1 = nx1 / sintheta;
ny1 = cross(n0,nx1);

sx1 = nx1(1)*sx + nx1(2)*sy + nx1(3)*sz;
sy1 = ny1(1)*sx + ny1(2)*sy + ny1(3)*sz;
sz1 = n0(1)*sx + n0(2)*sy + n0(3)*sz;
H1 = (w0 - w) * sz1/2 + w1 * sintheta * (cos(phi) * sx1 + sin(phi) * sy1) / 2;

U = expm(-1i * t * w * sz1 / 2) * expm(-1i * t * H1);
end


