function [ enGapSq, enGapSqGrad] = EnergyGapSq(bk,param,orderParam)
% EnergyGapSq finds the square of the minimum energy gap given the mean field order parameter u and v
% bk is the momentum (kx,ky,kz), orderParam= (u1,u2,u3,v1,v2,v3)
% param include all the input parameters: param=[g gamma gamma_z delta_z mu T N];
% output is the energy gap square enGapSq and the gradient of enGapSq

g=param(1); gamma_z=param(3); mu=param(5); 

u=orderParam(1:3); v=orderParam(4:6);

kx=bk(1); ky=bk(2); kz=bk(3); 

epsilon= 2*(cos(kx) + cos(ky) +gamma_z*cos(kz) ) - mu; % see note 
epsilonSq= epsilon.^2; 

u_sin_k= sin(kx)*u(1) + sin(ky)*u(2)+ gamma_z*sin(kz)*u(3);
v_sin_k= sin(kx)*v(1) + sin(ky)*v(2)+ gamma_z*sin(kz)*v(3);
Lsq= 16*g^2*( (u_sin_k).^2  + (v_sin_k).^2 );

enGapSq=epsilonSq + Lsq; % quasiparticle excitation spectrum square

%----------- Gradient of enGapSq: enGapSqGrad -----------% 
enGapSqGrad=zeros(3,1);
enGapSqGrad(1)= -4*sin(kx)*epsilon + 32*g^2*cos(kx)*(u(1)*u_sin_k+ v(1)*v_sin_k );
enGapSqGrad(2)= -4*sin(ky)*epsilon + 32*g^2*cos(ky)*(u(2)*u_sin_k+ v(2)*v_sin_k );
enGapSqGrad(3)= -4*gamma_z*sin(kz)*epsilon + 32*gamma_z*g^2*cos(kz)*(u(3)*u_sin_k+ v(3)*v_sin_k );
%----------- Gradient of enGapSq: enGapSqGrad -----------% 
