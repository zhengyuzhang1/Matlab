%restricted Boltzmann machine for 1D Anti-Fero Heisenberg model

N = 20; %visible layer nodes
M = 20; %hidden layer nodes

%initial values of variables of RBM
K = N*M+N+M;
% load('lasttime.mat','regX');
X = 0.1*(rand(K,1) + 1i * rand(K,1))-0.05*(1+1i)*ones(K,1);
a = X(1:N);
b = X(N+1:N+M);
W = reshape(X(N+M+1:end),M,N);
% a = [0;0];
% b = [0];
% W = [1.3366i -0.0985i];


%Monte Carlo parameters
maxIte = 1000; %max iteration number
regX = zeros(K,maxIte);
Nmont = 10000; %Monte Carlo sweep number
equilpoint = floor(0.5*Nmont); %calculating average after equilpoint samples, minimal value is 0
gamma0 = 0.5;
lambda0 = 100;
bb = 0.9;
bbg = 0.9;
lambdamin = 1e-4;
gammamin = 1e-2;

E = zeros(maxIte,1);
for ite = 1:maxIte
    %initial sample
    walker = 2*randi([0 1], N, 1) - 1;
    alpha = walker'*a;
    theta = b + W * walker;
    phi = Phi(alpha,theta);
    totalEloc = 0;
    totalconjO = zeros(K,1);
    totalElocTOC = zeros(K,1);
    
    for imont = 1:Nmont
        fs = randi(N); %fs is site of the spin to be flipped
        newalpha = alpha - 2 * walker(fs) * a(fs);
        newtheta = theta - 2 * walker(fs) * W(:,fs);
        newphi = Phi(newalpha,newtheta);
        p = min([1 abs(newphi/phi)^2]);
        if rand < p
            alpha = newalpha;
            theta = newtheta;
            phi = newphi;
            walker(fs) = -walker(fs);
        end
        if imont > equilpoint
            Eloc = 0;
            for j = 1:N
                jn = mod(j,N)+1;
                Eloc = Eloc + (1-walker(j)*walker(jn))/phi ...
                    *Phi(alpha-2*walker(j)*a(j)-2*walker(jn)*a(jn),theta-2*walker(j)*W(:,j)-2*walker(jn)*W(:,jn));
            end
            Eloc = Eloc + walker'* circshift(walker,1);
            conjO = conj([walker ; tanh(theta) ; kron(theta, walker)]);
            totalEloc = totalEloc + Eloc;
            totalconjO = totalconjO + conjO;
            totalElocTOC = totalElocTOC + Eloc * conjO;
        end
    end
    meanconjO = totalconjO/(Nmont-equilpoint);
    meanEloc = totalEloc/(Nmont-equilpoint);
    meanElocTOC = totalElocTOC/(Nmont-equilpoint);
    F = meanElocTOC - meanEloc * meanconjO;
    S =  totalOCTO/(Nmont-equilpoint)...
        - meanconjO * meanconjO';
    S = S + max(lambda0*bb^ite,lambdamin)*diag(diag(S));
    X = minresQLP(S,-max(gamma0*bbg^ite,gammamin)*F) + X;
    if any(isnan(X))
        pause;
    end
    regX(:,ite) = X;
    a = X(1:N);
    b = X(N+1:N+M);
    W = reshape(X(N+M+1:end),M,N);
    E(ite) = meanEloc/N;
end
% figure;plot(eps);
figure;plot([abs(E), real(E)]);
%%
% H = 2*[1 0 0 0; 0 -1 2 0; 0 2 -1 0; 0 0 0 1];
% [V,D] = eig(H);
% phi = [0.6533;0.2706;0.2706;0.6533];
% H*phi./phi
% finalE = phi'*H*phi/(phi'*phi)/N;