% This script performs the mean field calculation for a 3D p-wave Feshbach resonance system.
% It uses the function fminunc to find the minimum for the thermodynamic potential
% See the note OptLatPwave-Note for detailed calculations.
% use an index phase to indicate different phases:
% phase=0, normal gas; phase=1, pz wave;
% phase=2, px+ipy wave; phase=3, px+ i alpha py, unequal magnitude
% phase=4, over condensation;
% ------ This script explicitly includes the constraint u.^2+v.^2<=1 for the condensation, and u.*v==0 for the gauge. Use the fmincon function

g=1; % coupling strength
gamma=0; % Feshbach molecule chemical potential
gamma_z=1; % z direction anisotropy
delta_z=0; % z direction molecule chemical potential
mu=1; % chemical potential
T=10^-20; % temperature
N=200;
cr1= 10^-3; % normal gas criterion: |u| <=cr1 && |v| <= cr1, then phase=0
% if either |u| <=cr1 or |v| <=cr1, but NOT both, then one vector is zero, phase = 1 (Gauge is u.*v =0)
cr2=0.5; % over condensation criterion: sqrt(u^2+v^2)>cr2, then phase = 4
cr3=10^-3; % criterion for isotropic px+ ipy: |uMag-vMag|/uMag < cr3, then phase 2

%[varAll1, varAll2]=ndgrid([0.5:0.05:0.8, 0.825:0.025:1.2, 1.25:0.05:2], 2*mu+(1:0.05:10)); % two parameters are changed here; see below.
[varAll1, varAll2]=ndgrid([1,1.1], 2*mu+[1,2]); % two parameters are changed here; see below.

phaseAll=nan(size(varAll1)); OmegaAll=nan(size(varAll1)); uvMagAll=nan(size(varAll1)); enGapAll=nan(size(varAll1));
uvAll=nan(size(varAll1,1),size(varAll1,2),6);

%----------- Discretization: put them outside of the loop and ThermoPot function to speed up the calculation -----------% 
dk=2*pi/N; %discretization in k 
[kx, ky, kz]=ndgrid(0:dk:2*pi,0:dk:2*pi,0:dk:2*pi); % momentum space discretization
sin_kx=sin(kx); sin_ky=sin(ky); sin_kz=sin(kz);
cos_kx=cos(kx); cos_ky=cos(ky); cos_kz=cos(kz);
%----------- Discretization: put them outside of the loop and ThermoPot function to speed up the calculation -----------% 

for Ind1=1:size(varAll1,1);
    tic;
    for Ind2=1:size(varAll1,2);
        gamma_z=varAll1(Ind1,Ind2);
        gamma=varAll2(Ind1,Ind2);
        
        param=[g gamma gamma_z delta_z mu T N];
        
        epsilon= 2*(cos_kx + cos_ky +gamma_z*cos_kz ) - mu; % see note 
        epsilonSq= epsilon.^2; clear epsilon
        
        OmegaFunc=@(uv)ThermoPot(uv,param,sin_kx,sin_ky,sin_kz,epsilonSq); %input with (ux,uy,uz,vx,vy,vz)
        options = optimoptions(@fmincon,'Algorithm','interior-point','GradObj','on','Display','off','GradConstr','on','TolX',10^-7,'TolFun',10^-7,'TolCon',10^-10);
        
        % Sweeping along varAll1, i.e., use the orderParam for the previous run as the next initial condition
        if Ind2==1 
            if Ind1==1
                % The first point is obtained from random initial points
                % enforce the nonlinear constraint u^2+v^2<=1, and the gauge choice u.v==0.
                uvInit=orth(rand(3,2)); % generate two random orthgonal vectors with unit length.
                uvInit=[uvInit(:,1)*sqrt(0.5)*rand; uvInit(:,2)*sqrt(0.5)*rand]; % make the lengths of the vectors random.
            else
                uvInit=uvAll(Ind1-1,Ind2,:); uvInit=uvInit(:); % The first column, (i,1) obtained from (i-1,1)
            end
        else
            uvInit=uvAll(Ind1,Ind2-1,:); uvInit=uvInit(:); % other than the first column, (i,j) obtained from (i,j-1), sweeping along i
        end
        
        [orderParam,Omega]=fmincon(OmegaFunc,uvInit,[],[],[],[],[],[],@uvConstraint,options);
        
        u=orderParam(1:3); v=orderParam(4:6);
        uMag=sqrt(sum(u.^2)); vMag=sqrt(sum(v.^2)); uvMag=sqrt(sum(u.^2)+sum(v.^2));
        uD=u/uMag; vD=v/vMag;
                
        % Note the order of the different lines
        while 1
            if uvMag>=cr2, phase=4; break; end % full condensate criterion
            if uMag<cr1 && vMag<cr1, phase=0; break; end % normal free gas criterion
            if uMag<cr1 || vMag<cr1, phase=1; break; end % pz wave
            %if the above three criteria are not met, it is a px+ipy superfluid
            if abs(uMag-vMag)/uMag<cr3, phase=2; break; end % isotropic px+i py wave
            phase= 3; break; % anisotropic px + i alpha py
        end
        
        OmegaAll(Ind1, Ind2)= Omega; phaseAll(Ind1, Ind2)= phase; 
        uvMagAll(Ind1, Ind2)= uvMag; uvAll(Ind1,Ind2,:)= orderParam;
        
        %%%%% ----- Find the minimum band gap from the mean field parameters ----- %%%%%
%         EnergyGapSqFunc=@(x)EnergyGapSq(x,param,orderParam);
%         options2 = optimoptions(@fmincon,'Algorithm','trust-region-reflective','GradObj','on','Display','off','TolX',10^-7,'TolFun',10^-7);
%         bklower=[0,0,0]; bkupper=[2*pi,2*pi,2*pi]+0.1;
%         
%         problem = createOptimProblem('fmincon','objective',...
%             EnergyGapSqFunc,'x0',rand(3,1)*2*pi,'lb',bklower,'ub',bkupper,'options',options2);
%         gs = GlobalSearch('Display','off');
%         [~,enGapSq] = run(gs,problem);
%         enGap=sqrt(enGapSq);
%         enGapAll(Ind1, Ind2)= enGap;
        %%%%% ----- Find the minimum band gap from the mean field parameters ----- %%%%%
    end
    toc;
end

clear kx ky kz sin_kx sin_ky sin_kz cos_kx cos_ky cos_kz epsilonSq

FileName=sprintf('pwave_mu%3.1f_deltaz%3.1f_T%3.1f.mat',mu,delta_z,T);
save(FileName);
