fmin7 = ones(1,5);
taomin7 = zeros(5,7);
for i = 1:10
    load(['search7_',num2str(i)]);
    for j = 1:5
        if fmin{1,1}(j)<fmin7(j)
            fmin7(j) = fmin{1,1}(j);
            taomin7(j,:) = taomin{1,1}(j,:);
        end
    end
end