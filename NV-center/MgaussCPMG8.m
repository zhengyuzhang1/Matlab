function M = MgaussCPMG8(w,n,azx,t,pt)
%gauss pi-pulse CPMG, w is the Larmor frequency of C13, n the number of
%CPMG-8 pulses, azx is n*2 matrix of coupling coefficients, t is vector of free evolution time, and
%pt is the Pi pulse duration time
sz = [1,0;0,-1];
sx = [0,1;1,0];
sy = [0,-1i;1i,0];
M = zeros(size(azx,1),length(t));

s0 = kron(0.5 * ones(2), 0.5 * eye(2)); %initial density matrix
sigma = pt / 2 / 3; %standard deviation of gauss pulse
a = sqrt(pi / 2) / sigma / erf(pt / 2 / sqrt(2) / sigma); %maximal amplitude of pi pulse

nt = 100;
dt = pt / nt;
tn = linspace(-pt / 2, pt / 2, nt + 1);

for nuc = 1:size(azx,1)
    H0 = w / 2 * sz;
    H1 = (azx(nuc,1) + w) / 2 * sz + azx(nuc,2) / 2 * sx;
%     H = blkdiag(H0, H1);
%     xflipU = eye(4);
%     yflipU = eye(4);
%     for i = 1:nt
% %         xflipU = expm(-1i * dt * (H + a * exp(- tn(i)^2 / 2 / sigma^2) / 2 * kron(sx,eye(2)))) * xflipU;
% %         yflipU = expm(-1i * dt * (H + a * exp(- tn(i)^2 / 2 / sigma^2) / 2 * kron(sy,eye(2)))) * yflipU;
% %         xflipU = expm(-1i * dt * (H + pi / pt / 2 * kron(sx,eye(2)))) * xflipU;
% %         yflipU = expm(-1i * dt * (H + pi / pt / 2 * kron(sy,eye(2)))) * yflipU;
%     end
    for i = 1:length(t)
%         freeU = expm(-1i * t(i) * H);
%         U = (freeU * xflipU * freeU^2 * yflipU * freeU)^2 * (freeU * yflipU * freeU^2 * xflipU * freeU)^2;
%         U = U^n;
% %         U = (freeU * xflipU * freeU) ^ (8 * n);
%         sn = U * s0 * U';
%         M(nuc,i) = real(sn(3) + sn(8) + sn(9) + sn(14));
        U0 = (expm(-1i * t(i) * H0) * expm(-2i * t(i) * H1) * expm(-1i * t(i) * H0))^(n / 2);
        U1 = (expm(-1i * t(i) * H1) * expm(-2i * t(i) * H0) * expm(-1i * t(i) * H1))^(n / 2);
        M(nuc,i) = 1 / 2 * real(trace(U0 * 0.5 * ones(2) * U1'));
    end
end

