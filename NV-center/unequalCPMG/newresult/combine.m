nof = 30;
fminall = cell(3,3);
taominall = cell(3,3,nof);
for i = 1:3
    for j =1:3
        fminall{i,j} = zeros(2*i+3,nof);
        for k = 1:nof
            taominall{i,j,k} = zeros(2*i+3,2*j+3);
        end
    end
end
for i = 1:nof
    load(['search_',num2str(i),'.mat'],'fmin','taomin');
    if i == 1
        for j = 1:3
            for k = 1:3
                fminall{j,k}(:,1) = fmin{j,k};
                taominall{j,k,1} = taomin{j,k};
            end
        end
    else
        for j = 1:3
            for k = 1:3
                for l = 1:2*j+3
                    if fmin{j,k}(l)<fminall{j,k}(l,i-1)
                        fminall{j,k}(l,i) = fmin{j,k}(l);
                        taominall{j,k,i}(l,:) = taomin{j,k}(l,:);
                    else
                        fminall{j,k}(l,i) = fminall{j,k}(l,i-1);
                        taominall{j,k,i}(l,:) = taominall{j,k,i-1}(l,:);
                    end
                end
            end
        end
    end                        
end
%% 
for i = 1:3
    for j = 1:3
        for k = 1:nof
            taominall{i,j,k}(:,2*j+4) = (sum(taominall{i,j,k}(:,2:2:end),2)*2).*taominall{i,j,k}(:,end);
        end
    end
end
%%
for i = 1:1
    for j = 1:1
        figure;plot((100:100:3000),log10(fminall{i,j}'));
    end
end
xlabel('optimizing times');
ylabel('log10(error)');
title('optimization for one sample of 5 nuclei');
fig2pdf2('optimizing times');
%% combine random
nor = 10;
fr = zeros(5*nor,3);
for i = 1:nor
    load(['searchrandom_',num2str(i),'.mat'],'w1','fmin');
    fr(5*i-4:5*i,1) = fmin{1,1};
    fr(5*i-4:5*i,2) = w1{1,1}(:,3);
    fr(5*i-4:5*i,3) = w1{1,1}(:,1);
end
figure;plot(fr(:,1));xlabel('nucleus');ylabel('error');
fig2pdf2('random_case1');
figure;plot(fr(:,2));xlabel('nucleus');ylabel('Azx/w1');
fig2pdf2('random_case2');
%% segments
seg = cell(3,1);
for i = 1:3
    seg{i,1} = zeros(2*i+3,3);
    seg{i,1} = [fminall{i,1}(:,end),fminall{i,2}(:,end),fminall{i,3}(:,end)];
end
figure;subplot(1,3,1);plot([5,7,9],seg{1,1}');xlabel('number of segments');ylabel('error');title('error on one sample of 5 nuclei');
subplot(1,3,2);plot([5,7,9],seg{2,1}');xlabel('number of segments');ylabel('error');title('error on one sample of 7 nuclei');
subplot(1,3,3);plot([5,7,9],seg{3,1}');xlabel('number of segments');ylabel('error');title('error on one sample of 9 nuclei');
fig2pdf(gcf,[10 5],'segment');
%% number of nuclei
nucn = zeros(3);
for i = 1:3
    for j = 1:3
        nucn(i,j) = mean(fminall{i,j}(:,end));
    end
end
figure;plot([5 7 9],nucn');xlabel('number of nuclei');ylabel('average error');
legend('5 segments','7 segments','9 segments');
fig2pdf2('number of nuclei');
