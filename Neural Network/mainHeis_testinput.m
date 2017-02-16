%restricted Boltzmann machine for 1D Anti-Ferro Heisenberg model

N = 20; %visible layer nodes
M = 1; %hidden layer nodes

%initial values of vareiables of RBM
K = 1+M+M*N;

%X = 1/N*(rand(K,1) + 1i * rand(K,1))-1/2/N*(1+1i)*ones(K,1);
% load('heis20.mat','regX');
%X = regX(:,end-1);

% read or write to file for fortran debug
% fileID = fopen('data.txt','w');
% for i = 1:K
%     fprintf(fileID, '%f %f\n',real(X(i)),imag(X(i)));
% end
fileID = fopen('SYK10_10.txt','r');
for i = 1:626
    fgetl(fileID);
end
tline = fgetl(fileID);
C = strsplit(tline, {',', '(',')',') (',' '});
C = C(3:end);
X = zeros(K,1);
for i = 1:K
    X(i) = str2double(C{2*i-1})+str2double(C{2*i})*1i;
end

a = X(1);
b = X(2:1+M);
W = reshape(X(1+M+1:end),M,N);
% a = [0;0];
% b = [0];
% W = [1.3366i -0.0985i];


%Monte Carlo parameters
maxIte = 3; %max iteration number
regX = zeros(K,maxIte);
Nmont = 10000; %Monte Carlo sweep number
equilpoint = 0; %calculating average after equilpoint samples, minimal value is 0
Nmonte = Nmont - equilpoint; %effective Monte Carlo sample number
gamma0 = 0.05;
lambda0 = 0;
bb = 0.9;
bbg = 0.95;
lambdamin = 0.000;
gammamin = 0.001;

E = zeros(maxIte,1);
flipnumber = zeros(maxIte,1);
ite = 1;
while ite < maxIte + 1
    %initial sample
    walker = 2*randi([0 1], N, 1) - 1;
%    C = strsplit(fgetl(fileID));
%     for i = 1:N
%         walker(i) = str2double(C{i});
%     end
%     fprintf(fileID,'\n');
    alpha = sum(walker)*a;
    theta = repmat(b,1,N) + W * TransInv(walker);
    phi = Phi(alpha,theta);
    meanEloc = 0;
    meanElocTOC = zeros(K,1);
    conjO = zeros(K,Nmont-equilpoint);
    flipnumber(ite) = 0;
    for imont = 1:Nmont
         fs = randi(N); %fs is site of the spin to be flipped
%         C = strsplit(fgetl(fileID));
%         fs = str2double(C{1});
        newalpha = alpha - 2 * walker(fs) * a;
        newtheta = theta - 2 * walker(fs) * W * circshift(eye(N),fs-1);
        newphi = Phi(newalpha,newtheta);
        p = min([1 abs(newphi/phi)^2]);
        rtemp = rand;
%         fprintf(fileID,'%d %f\n', fs, rtemp);
%         rtemp = str2double(C{2});
        if rtemp < p
            flipnumber(ite) = flipnumber(ite)+1;
            alpha = newalpha;
            theta = newtheta;
            phi = newphi;
            walker(fs) = -walker(fs);
        end
        if imont > equilpoint
            Eloc = 0;
            for j = 1:N
                jn = mod(j,N)+1;
                if walker(j)~=walker(jn)
                    Eloc = Eloc + 2/phi ...
                        *Phi(alpha-2*(walker(j)+walker(jn))*a,theta-2*W*(walker(j)*circshift(eye(N),j-1)+walker(jn)*circshift(eye(N),jn-1)));
                end
            end
            Eloc = Eloc + walker'* circshift(walker,1);
            conjO(:,imont-equilpoint) = conj([sum(walker) ; sum(tanh(theta),2) ; reshape(tanh(theta)*TransInv(walker)',M*N,1)]);
            meanEloc = meanEloc + Eloc/Nmonte;
            meanElocTOC = meanElocTOC + Eloc * conjO(:,imont-equilpoint)/Nmonte;
        end
    end
    E(ite) = meanEloc/N;
%    if isnan(meanEloc) || ( ite>2 && E(ite)-E(ite-1)>0.2*abs(E(ite-1)) )
%        a = regX(1,ite-2);
%        b = regX(2:M+1,ite-2);
%        W = reshape(regX(M+2:end,ite-2),M,N);
%        ite = ite - 1;
%        continue;
%    end
    meanconjO = mean(conjO,2);
    F = meanElocTOC - meanEloc * meanconjO;
    lambda =  max(lambda0*bb^(ite),lambdamin);
    A = @(x)S(x,conjO,meanconjO,lambda);
%     gamma = max(gamma0*bbg^(ite),gammamin);
    gamma = 0.01;
    [dX,flag] = minresQLP(A,-gamma*F,[],10*N);
%    if flag ~= 1
%        continue;
%    end
    if ite > 1
        regX(:,ite) = regX(:,ite-1) + dX;
    else
        regX(:,ite) = X + dX;
    end
    a = regX(1,ite);
    b = regX(2:M+1,ite);
    W = reshape(regX(M+2:end,ite),M,N);
    ite = ite + 1;
%    if mod(ite,10) == 0
%        save(mfilename);
%    end
end

% fclose(fileID);
display(E);
%%
% sx = [0 1;1 0];sy = [0 -1i;1i 0];sz = [1 0;0 -1];
% H = kron(kron(sx, sx),eye(2)) + kron(kron(sy, sy),eye(2)) + kron(kron(sz, sz),eye(2)) ...
%     + kron(eye(2),kron(sx, sx)) + kron(eye(2),kron(sy, sy)) + kron(eye(2),kron(sz, sz)) ...
%     + kron(sx,kron(eye(2), sx)) + kron(sy,kron(eye(2), sy)) + kron(sz,kron(eye(2), sz));
% [V,D] = eig(H);
% phi = [0.6533;0.2706;0.2706;0.6533];
% H*phi./phi
% finalE = phi'*H*phi/(phi'*phi)/N;

% newphi = zeros(N,1);
% for i = 1:N
%     neww = walker;
%     neww(i) = - neww(i);
%     neww = 2*randi([0 1], N, 1) - 1;
%     alpha = sum(neww)*a;
%     theta = repmat(b,1,N) + W * TransInv(neww);
%     newphi(i) = Phi(alpha,theta);
% end