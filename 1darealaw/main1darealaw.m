rng('shuffle');
N = 10;
mode = 3; % 1: no minus sign; 2: random minus sign before different i,j term; 3:  
boundary = 'p'; % 'o' for open boundary; 'p' for periodic boundary
%% Periodic boundary condition

alpha_list = [ (0.1:0.1:1), 2, 3 ];
entropy = zeros(floor(N/2)-1,length(alpha_list));
parfor alpha_idx = 1:length(alpha_list)
    co = zeros(N,N,3);
    alpha = alpha_list(alpha_idx);
    for i = 1:N-1
        for j = i+1:N
            if mode == 1
                coef = [1,1,1];
            elseif mode == 2
                coef = sign([rand * ones(2,1); rand]-0.5);
            elseif mode == 3
                coef = sign(rand(3,1)-0.5);
            end
            if boundary == 'p'
                co(i,j,:) = coef / min(j-i,i+N-j)^alpha;
            elseif boundary == 'o'
                co(i,j,:) = coef / (j-i)^alpha;
            end                
        end
    end
    maxHnorm = sum(sum(triu(sum(abs(co),3))));
    H = Sparse_Heisenlike_alltoall(N,co);
    [v, d] = eigs(H - maxHnorm*speye(2^N),1);
    d = d + maxHnorm;
    ent = zeros(floor(N/2)-1,1);
    for i = 2:floor(N/2)
       density_red = partrace_pure_large(v,i+1:N,2*ones(N,1));
       ent(i-1) = Entropy(density_red);
    end
    entropy(:,alpha_idx) = ent;
end
save(['arealaw_',boundary,num2str(mode),'.mat'],'entropy');