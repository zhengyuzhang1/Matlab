%restricted Boltzmann machine for 1D Anti-Ferro Heisenberg model using
%translational symmetry

N = 10; %number of visible layer nodes
M = 1; %number of hidden layer nodes = M*N. Translational symmetry => M = 1

%initial values of vareiables of RBM
K = 1+M+M*N;
X = 1/N*((rand(K,1)-0.5));% + 1/N*1i*(rand(K,1)-0.5);
% load('heis20.mat','regX');
%X = regX(:,end-1);

% read or write to file for fortran debug
% fileID = fopen('data.txt','w');
% for i = 1:K
%     fprintf(fileID, '%f %f\n',real(X(i)),imag(X(i)));
% end
%fileID = fopen('data.txt','r');
% for i = 1:K
%     tline = fgetl(fileID);
%     C = strsplit(tline);
%     X(i) = str2double(C{1})+str2double(C{2})*1i;
% end

a = X(1); %parameters on visible nodes
b = X(2:1+M);%parameters on hidden nodes
W = reshape(X(1+M+1:end),M,N);%links between visible and hidden nodes



%Monte Carlo parameters
maxIte = 100; %max iteration number
regX = zeros(K,maxIte);
Nmont = 1000; %Monte Carlo sweep number
equilpoint = 0; %calculating average after equilpoint samples, minimal value is 0
Nmonte = Nmont - equilpoint; %effective Monte Carlo sample number
%gamma is step length of gradient descend
gamma0 = 0.05;
gammarate = 0.9;
gammamin = 0.01;
%lambda is the same as in Carleo's paper, here for quicker converge I set
%to 0
lambda0 = 0;
lambdarate = 0.9;
lambdamin = 0.000;

E = zeros(maxIte,1);
%flipnumber = zeros(maxIte,1);
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
    meanElocTOC = zeros(K,1); %mean of Eloc Times Conjugate(O)
    conjO = zeros(K,Nmont-equilpoint);
%    flipnumber(ite) = 0;
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
        if rtemp < p %flip
%             flipnumber(ite) = flipnumber(ite)+1;
            alpha = newalpha;
            theta = newtheta;
            phi = newphi;
            walker(fs) = -walker(fs);
        end
        if imont > equilpoint %run sampling
            Eloc = 0;
            for j = 1:N
                jn = mod(j,N)+1;
                if walker(j)~=walker(jn)
                    Eloc = Eloc - 2/phi ...
                        *Phi(alpha-2*(walker(j)+walker(jn))*a,theta-2*W*(walker(j)*circshift(eye(N),j-1)+walker(jn)*circshift(eye(N),jn-1)));
                end
            end
            Eloc = Eloc + walker'* circshift(walker,1);
            conjO(:,imont-equilpoint) = conj([sum(walker) ; sum(tanh(theta),2) ; reshape(tanh(theta)*TransInv(walker)',M*N,1)]);
            meanEloc = meanEloc + Eloc/Nmonte;
            meanElocTOC = meanElocTOC + Eloc * conjO(:,imont-equilpoint)/Nmonte;
        end
    end
    E(ite) = meanEloc/N;%energy per site
%    if isnan(meanEloc) || ( ite>2 && E(ite)-E(ite-1)>0.2*abs(E(ite-1)) )
%        a = regX(1,ite-2);
%        b = regX(2:M+1,ite-2);
%        W = reshape(regX(M+2:end,ite-2),M,N);
%        ite = ite - 1;
%        continue;
%    end
    meanconjO = mean(conjO,2);
    F = meanElocTOC - meanEloc * meanconjO;
    lambda =  max(lambda0*lambdarate^(ite),lambdamin);
    A = @(x)S(x,conjO,meanconjO,lambda);
    gamma = max(gamma0*gammarate^(ite),gammamin);
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
figure;plot([abs(E), real(E)]);
