%restricted Boltzmann machine for 1D transverse-field Ising (TFI) model

N = 40; %visible layer nodes
M = 1; %hidden layer nodes
h = 2; %Transverse field strength

%initial values of variables of RBM
K = N*M+N+M;
X = X100;
% X = 0.1*(rand(K,1) + 1i * rand(K,1));
a = X(1:N);
b = X(N+1:N+M);
W = reshape(X(N+M+1:end),M,N);
% a = [0;0];
% b = [0];
% W = [1.3366i -0.0985i];


%Monte Carlo parameters
maxIte = 100; %max iteration number
Nmont = 100000; %Monte Carlo sweep number
equilpoint = floor(0.5*Nmont); %calculating average after equilpoint samples, minimal value is 0
gamma0 = 0.1;
lambda0 = 100;
bb = 0.9;
lambdamin = 1e-4;
gammamin = 1e-3;

regwalker = zeros(N,Nmont);
O = zeros(length(X),Nmont-equilpoint);
E = zeros(maxIte,1);
for ite = 1:maxIte
    %initial sample
    walker = 2*randi([0 1], N, 1) - 1;
    alpha = walker'*a;
    theta = b + W * walker;
    phi = Phi(walker'*a,theta);
    Eloc = zeros(1,Nmont-equilpoint);
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
        regwalker(:,imont) = walker;
        if imont > equilpoint
            i = imont - equilpoint;
            for j = 1:N
                Eloc(i) = Eloc(i)-h * Phi(alpha-2*walker(j)*a(j),theta-2*walker(j)*W(:,j))/phi;
            end
            Eloc(i) = Eloc(i) - walker'* circshift(walker,1);
            O(1:N,i) = walker;
            O(N+1:N+M,i) = tanh(theta);
            O(N+M+1:end,i) = kron(theta, walker);
        end
    end
    meanconjO = mean(conj(O),2);
    F = mean(bsxfun(@times,conj(O),Eloc),2) - mean(Eloc) * meanconjO;
    S = mean(bsxfun(@times,reshape(conj(O),K,1,Nmont-equilpoint),reshape(O,1,K,Nmont-equilpoint)),3) ...
        - meanconjO * meanconjO';
    S = S + max(lambda0*bb^ite,lambdamin)*diag(diag(S));
    X = minresQLP(S,-max(gamma0*bb^ite,gammamin)*F) + X;
    a = X(1:N);
    b = X(N+1:N+M);
    W = reshape(X(N+M+1:end),M,N);
    E(ite) = mean(Eloc)/N;
end
Ereal = 2/pi*(1+1/h)*ellipticE(pi/2,2*sqrt(1/h)/(1+1/h));
eps = abs(abs(E)-Ereal)/Ereal;
% figure;plot(eps);
figure;plot([abs(E), real(E)]);
%%
% H = [-2 -1 -1 0; -1 2 0 -1; -1 0 2 -1; 0 -1 -1 -2];
% phi = [0.6533;0.2706;0.2706;0.6533];
% H*phi./phi
% finalE = phi'*H*phi/(phi'*phi)/N;