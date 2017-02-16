%% This test script test the function EnergyGapSq

%% Test 1
tic;
g=1; gamma=0; gamma_z=0.2; delta_z=0; mu=0.5; T=10^-20; N=50;
orderParam=[-1.025462167212862,2.065145482098000,0,2.065145482098000,1.025462167212862,0];
param=[g gamma gamma_z delta_z mu T N];
EnergyGapSqFunc=@(x)EnergyGapSq(x,param,orderParam);
options = optimoptions(@fmincon,'Algorithm','trust-region-reflective','GradObj','on','Display','off','TolX',10^-7,'TolFun',10^-7);
bklower=[0,0,0]; bkupper=[2*pi,2*pi,2*pi]+0.1;

problem = createOptimProblem('fmincon','objective',...
    EnergyGapSqFunc,'x0',rand(3,1)*2*pi,'lb',bklower,'ub',bkupper,'options',options);
gs = GlobalSearch;
[bk,enGapSq] = run(gs,problem);
enGap=sqrt(enGapSq);
toc;
%[bk,enGapSq]=fmincon(EnergyGapSqFunc,rand(3,1)*2*pi,[],[],[],[],bklower,bkupper,[],options);

