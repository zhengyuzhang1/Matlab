rng('shuffle');
%restricted Boltzmann machine for Heisenberg model
ntry = 80;
sr = 0;
N = 6; %visible layer nodes
M = 1; %hidden layer nodes
hnondiag = 2;

%initial values of vareiables of RBM
K = 1+M+M*N;
maxIte = 100; %max iteration number
dof = 2^N; %state numbers

lambda0 = 1;
lambdarate = 0.9;
lambdamin = 0.0001;
gamma0 = 0.05;
gammarate = 0.9;
gammamin = 0.001;

finalE = zeros(ntry,1);
finalX = zeros(K,ntry);
initialX = zeros(K,ntry);
parfor itry = 1:ntry
    X_re = rand(K,1) - 0.5;
    X_im = 1i*(rand(K,1) - 0.5);
    X = 0.1*(X_re + X_im);
%    X(1) = 0;
    initialX(:,itry) = X;
    regX = zeros(K,maxIte);

    E = zeros(maxIte,1);
    allp = zeros(1,dof);
    flipnumber = zeros(maxIte,1);
    meanconjO = zeros(K,1);
    meanElocTOC = zeros(K,1);
    conjO = zeros(K,dof);
    ElocTO = zeros(K,1);
    phi_log = zeros(1,dof);
    ite = 1;
    while ite < maxIte + 1
        a = X(1);
        b = X(1+1:1+M);
        W = reshape(X(1+M+1:end),M,N);
        %initial sample
        meanEloc = 0;
        meanconjO(:) = 0;
        meanElocTOC(:) = 0;
        conjO(:) = 0;
        for imont = 1:dof
            walker = de2bi(dof-imont,N,'left-msb')';
            walker(walker == 0 ) = -1;
            alpha = sum(walker)*a;
            theta = repmat(b,1,N) + W * TransInv(walker);
            phi_log(imont) = Phi_log(alpha,theta);
        end
        maxphi_log = max(real(phi_log));
        phi_log = phi_log - maxphi_log + 20;
        phi_log(real(phi_log)<0) = 0;
        allp = exp(2*real(phi_log));
        allp = allp / sum(allp);
        for imont = 1:dof
            walker = de2bi(dof-imont,N,'left-msb')';
            walker(walker == 0 ) = -1;
            alpha = sum(walker)*a;
            theta = repmat(b,1,N) + W * TransInv(walker);
            Eloc = 0;
            for j = 1:N
                jn = mod(j,N)+1;
                if walker(j)~=walker(jn)
                    imontn = imont + walker(j)*2^(N-j) + walker(jn)*2^(N-jn);
                    Etemp = hnondiag * exp( phi_log(imontn) - phi_log(imont) ) ;
                    Eloc = Eloc + Etemp;
                end
            end
            hdiag = walker'* circshift(walker,1);
            Eloc = Eloc + hdiag;
            conjO(:,imont) = conj([sum(walker) ; sum(tanh(theta),2) ; reshape(tanh(theta)*TransInv(walker)',M*N,1)]);
            meanEloc = meanEloc + Eloc * allp(imont);
            meanconjO = meanconjO + conjO(:,imont) * allp(imont);
            meanElocTOC = meanElocTOC + Eloc * conjO(:,imont) * allp(imont);
    %        meanconjOTO = meanconjOTO + conjO(:,imont) * conjO(:,imont)'*allp(imont);
        end
        E(ite) = meanEloc/N;
    %    if isnan(meanEloc) || ( ite>2 && E(ite)-E(ite-1)>0.2*abs(E(ite-1)) )
    %        a = regX(1,ite-2);
    %        b = regX(2:M+1,ite-2);
    %        W = reshape(regX(M+2:end,ite-2),M,N);
    %        ite = ite - 1;
    %        continue;
    %    end
        F = meanElocTOC - meanEloc * meanconjO;
        gamma = max(gamma0*gammarate^(ite),gammamin);
        if sr == 1 
            lambda =  max(lambda0*lambdarate^(ite),lambdamin);
            A = @(x)S(x,conjO.*repmat(sqrt(dof*allp),K,1),meanconjO,lambda);
            [dX,flag] = minresQLP(A,-gamma*F,[],10*N);
        else
            dX = - gamma * F;
        end
    %    if flag ~= 1
    %        continue;
    %    end
        X = X + dX;
        regX(:,ite) = X;
        ite = ite + 1;
    %    if mod(ite,10) == 0
    %        save(mfilename);
    %    end
    end
    %figure;plot([abs(E), real(E)]);
    finalE(itry) = real(E(end));
    finalX(:,itry) = X;
    display(itry);
end
figure;plot(finalE);


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