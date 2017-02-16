function nucs = scan(sample,Aall,n,d,x)
%scan first n frequencies in Aall
nucs = zeros(1,n);
for i = 1:n
    if isempty(sample)
        break;
    end
    if abs(sample(1)-Aall(i)) < 0.1
        nucs(i) = 0.5;  
        while ~isempty(sample) && abs(sample(1)-Aall(i)) < 0.1
            sample(1) = [];
        end
        continue;
    end
    frq = Aall(i);
    t = 0.1;
    M = 1;
    for j = 1:length(sample)
        if t * (frq - sample(j)) < d(end)
            [~,ind] = min(abs(d - t * (frq - sample(j))));
            M = M / 2 * trace(x(:,:,ind)');
        end
    end
    nucs(i) = (1 - real(M)) / 2;   
end
figure;plot(nucs,'o');
end
