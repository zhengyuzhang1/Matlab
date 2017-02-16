%restricted Boltzmann machine for Heisenberg model

N = 3; %visible layer nodes
M = 1; %hidden layer nodes
hnondiag = -2;

%initial values of vareiables of RBM
K = 1+M+M*N;
% X =(-0.02:0.01:0.03)';
X = rand(K,1)-0.5;% + 1i * rand(K,1)-1/2*(1+1i)*ones(K,1);
% load('heis20.mat','regX');
% X = regX(:,end-1);
a = X(1);
b = X(1+1:1+M);
W = reshape(X(1+M+1:end),M,N);
maxIte = 100; %max iteration number
regX = zeros(K,maxIte);
dof = 2^N; %state numbers

lambda0 = 0;
lambdarate = 0.9;
lambdamin = 0.0000;
gamma0 = 0.05;
gammarate = 0.9;
gammamin = 0.001;

E = zeros(maxIte,1);
allp = zeros(dof,1);
meanconjO = zeros(K,1);
meanElocTOC = zeros(K,1);
conjO = zeros(K,dof);
ite = 1;
while ite < maxIte + 1
    %initial sample
    meanEloc = 0;
    meanconjO(:) = 0;
    meanElocTOC(:) = 0;
    conjO(:) = 0;
    for imont = 1:dof
        walker = de2bi(imont-1,N)';
        %walker(walker == 0 ) = -1;
        alpha = sum(walker)*a;
        theta = repmat(b,1,N) + W * TransInv(walker);
        philog = Phi_log_bin(alpha,theta);
        allp(imont) = abs(exp(philog))^2;
        if isnan(allp(imont)) || isinf(allp(imont))
            display('allp');
            pause();
        end
        Eloc = 0;
        for j = 1:N
            jn = mod(j,N)+1;
            if walker(j)~=walker(jn)
                thetanew = theta-2*W*(walker(j)*circshift(eye(N),j-1)+walker(jn)*circshift(eye(N),jn-1));
                Etemp = hnondiag * exp(Phi_log_bin(alpha-2*(walker(j)+walker(jn))*a,thetanew) + philog');
                if isnan(Etemp) || isinf(Etemp)
                    display('Etemp');
                    pause();
                end
                Eloc = Eloc + Etemp;
            end
        end
        walker1 = walker;
        walker1(walker1==0)=-1;
        hdiag = walker1'* circshift(walker1,1);
        Eloc = Eloc + hdiag * allp(imont);
        conjO(:,imont) = conj([sum(walker) ; sum(exp(theta)./(1+exp(theta)),2) ; reshape((exp(theta)./(1+exp(theta)))*TransInv(walker)',M*N,1)]);
        meanEloc = meanEloc + Eloc;
        meanconjO = meanconjO + conjO(:,imont) * allp(imont);
        meanElocTOC = meanElocTOC + Eloc * conjO(:,imont);
%        meanconjOTO = meanconjOTO + conjO(:,imont) * conjO(:,imont)'*allp(imont);
    end
    meanEloc = meanEloc / sum(allp);
    meanElocTOC = meanElocTOC / sum(allp);
    E(ite) = meanEloc/N;
%    if isnan(meanEloc) || ( ite>2 && E(ite)-E(ite-1)>0.2*abs(E(ite-1)) )
%        a = regX(1,ite-2);
%        b = regX(2:M+1,ite-2);
%        W = reshape(regX(M+2:end,ite-2),M,N);
%        ite = ite - 1;
%        continue;
%    end
    meanconjO = meanconjO / sum(allp);
%    meanconjOTO = meanconjOTO / sum(allp);
    F = meanElocTOC - meanEloc * meanconjO;
    lambda =  max(lambda0*lambdarate^(ite),lambdamin);
    allp = allp / sum(allp);
    A = @(x)S(x,conjO.*repmat(sqrt(dof*allp'),K,1),meanconjO,lambda);
%    A = meanconjOTO - meanconjO * meanconjO';
%    A = A + lambda*diag(diag(A));
    [dX,flag] = minresQLP(A,-max(gamma0*gammarate^(ite),gammamin)*F,[],10*N);
%    dX = -max(gamma0*gammarate^(ite),gammamin)*pinv(A)*F;
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
% figure;plot(eps);
figure;plot([abs(E), real(E)]);


%%
% 
% varE = zeros(16,100);
% varX = linspace(-0.5, 0.5, 100);
% finalstate = zeros(dof,1);
% for i = 1:16
%     for k = 1:100
%         X = finalX;
%         if mod(i,2) == 1
%             X(floor((i+1)/2)) = X(floor((i+1)/2)) + varX(k);
%         else
%             X(floor((i+1)/2)) = X(floor((i+1)/2)) + varX(k)*1i;
%         end
%         a = X(1);
%         b = X(2:1+M);
%         W = reshape(X(1+M+1:end),M,N);
%         meanEloc = 0;
%         meanconjO = zeros(K,1);
%         meanElocTOC = zeros(K,1);
%         meanconjOTO = zeros(K,K);
%         conjO = zeros(K,dof);
%         for imont = 1:dof
%             walker = de2bi(imont-1,N)';
%             walker(walker == 0 ) = -1;
%             alpha = sum(walker)*a;
%             theta = repmat(b,1,N) + W * TransInv(walker);
%             phi = Phi(alpha,theta);
%             finalstate(imont) = phi;
%             allp(imont) = abs(phi)^2;
%             Eloc = 0;
%             for j = 1:N
%                 jn = mod(j,N)+1;
%                 if walker(j)~=walker(jn)
%                     Eloc = Eloc + 2/phi ...
%                         *Phi(alpha-2*(walker(j)+walker(jn))*a,theta-2*W*(walker(j)*circshift(eye(N),j-1)+walker(jn)*circshift(eye(N),jn-1)));
%                 end
%             end
%             Eloc = Eloc + walker'* circshift(walker,1);
%             conjO(:,imont) = conj([sum(walker) ; sum(tanh(theta),2) ; reshape(tanh(theta)*TransInv(walker)',M*N,1)]);
%             meanEloc = meanEloc + Eloc*allp(imont);
%             meanconjO = meanconjO + conjO(:,imont)*allp(imont);
%             meanElocTOC = meanElocTOC + Eloc * conjO(:,imont)*allp(imont);
%             meanconjOTO = meanconjOTO + conjO(:,imont) * conjO(:,imont)'*allp(imont);
%         end
%         meanEloc = meanEloc / sum(allp);
%         meanElocTOC = meanElocTOC / sum(allp);
%         varE(i,k) = meanEloc/N;
%     end
% end
% finalstate = finalstate / sqrt(sum(abs(finalstate).^2));
% %%
% figure;plot(real(varE(11:12,:)'))