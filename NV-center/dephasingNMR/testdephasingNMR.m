deltan = 2*pi*([(0.1:0.05:0.5),(0.6:0.2:3)]); %unit = Mhz
gamma2n = 2./ (50:50:500).^2; %unit = 1/ms

gaussP = zeros(3,length(deltan),length(gamma2n));
mininfidGauss = ones(length(deltan),length(gamma2n));
% pSquare = zeros(2,length(deltan),length(gamma2n));
% mininfidSquare = ones(length(deltan),length(gamma2n));

gaussA = pi;
gaussTn = [1,5,10]; %unit = ms
for deltai = 1:length(deltan)
    delta = deltan(deltai);
    for gamma2i = 1:length(gamma2n)
        display([deltai,gamma2i]);
        gamma2 = gamma2n(gamma2i);
        funGauss = @(x) optimGauss(x,delta,gamma2);
%         funSquare = @(x) optimSquare(x,delta,gamma2);
        for gaussTi = 1:length(gaussTn)
            gaussT = gaussTn(gaussTi) / delta;
            gaussSigma = gaussT/4;
            [temp1,temp2] = fminsearch(funGauss,[gaussA;gaussSigma;gaussT]);
            if temp2 < mininfidGauss(deltai,gamma2i)
                gaussP(:,deltai,gamma2i) = temp1;
                mininfidGauss(deltai,gamma2i) = temp2;
            end
%             [temp1,temp2] = fminsearch(funSquare,[pi/gaussT;gaussT]);
%             if temp2 < mininfidSquare(deltai,gamma2i)
%                 pSquare(:,deltai,gamma2i) = temp1;
%                 mininfidSquare(deltai,gamma2i) = temp2;
%             end            
        end
    end
end
%% infidelity
[X,Y] = meshgrid(deltan(1:12)/2/pi,(50:50:200));
figure;surface(X,Y,mininfidGauss(1:12,1:4)');
xlabel('\Delta(Mhz)');ylabel('T_{2}^{*}(\mus)');zlabel('infidelity');
%% T
i = length(deltan);j = length(gamma2n);
figure;surface(X,Y,reshape(gaussP(3,1:12,1:4),12,4)');
xlabel('\Delta(Mhz)');ylabel('T_{2}^{*}(\mus)');zlabel('T(\mus)');
%% cut-off ratio
% figure;surface(X,Y,reshape(gaussP(1,:,:).*erf(gaussP(3,:,:)./gaussP(2,:,:)/2/sqrt(2)),i,j)');
% figure;surface(X,Y,reshape(gaussP(1,:,:)/sqrt(2*pi)./gaussP(2,:,:).*exp(-(gaussP(3,:,:)./gaussP(2,:,:)/2/sqrt(2)).^2),i,j)');
% figure;surface(X,Y,reshape(gaussP(1,:,:)/sqrt(2*pi)./gaussP(2,:,:),i,j)');
figure;surface(X,Y,reshape(exp(-(gaussP(3,1:12,1:4)./gaussP(2,1:12,1:4)/2/sqrt(2)).^2),12,4)');
xlabel('\Delta(Mhz)');ylabel('T_{2}^{*}(\mus)');zlabel('cut-off ratio');
%% 
gaussP = gaussP(:,end,end);
w = @(t) gaussP(1) / sqrt(2*pi) / gaussP(2) * exp(-(t-gaussP(3)/2).^2/2/gaussP(2)^2);
figure;plot(linspace(0,gaussP(3),100),w(linspace(0,gaussP(3),100)))
