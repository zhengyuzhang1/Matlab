function [ Om, OmGrad] = ThermoPot(x,param,sin_kx,sin_ky,sin_kz,epsilonSq)
% ThermoPot find the thermodynamic potential in terms of the mean field order parameter u and v
% x is the order parameter (ux,uy,uz,vx,vy,vz)
% param include all the input parameters: param=[g gamma gamma_z delta_z mu T N];
% sincosFcn includes the functions (sin_kx,sin_ky,sin_kz,epsilonSq), which are put outside for faster speed
% output is the thermodynamic potential Om and the gradient of Om.
% integral makes use of the Simpson's rule, accurate to O(dk^4). See Wikipedia etc

g=param(1);  gamma=param(2); gamma_z=param(3);
delta_z=param(4); mu=param(5); T=param(6); N=param(7);

u= x(1:3); v=x(4:6);

%----------- thermodynamic Potential Om -----------% 
dk=2*pi/N; %discretization in k 

u_sin_k= sin_kx*u(1) + sin_ky*u(2)+ gamma_z*sin_kz*u(3);
v_sin_k= sin_kx*v(1) + sin_ky*v(2)+ gamma_z*sin_kz*v(3);

Lsq= 16*g^2*( (u_sin_k).^2  + (v_sin_k).^2 );
En= sqrt(epsilonSq + Lsq); % quasiparticle excitation spectrum
Om_Int=1/2 *En+T*log( 1+ exp(-En/T) ); % integrand in the integration   
Om= -1/(2*pi)^3 * Simp3D(Om_Int,N,dk)  + (gamma-2*mu)*sum(u.*u +v.*v) + delta_z*(u(3)^2+v(3)^2);  % thermodynamic potential
%----------- thermodynamic Potential Om -----------% 

%----------- Gradient of Om: OmGrad -----------% 
OmGrad=zeros(6,1);
EnFactor= (1-2./(exp(En/T) +1) )./En; % refer to the note
uFactor=u_sin_k.*EnFactor; vFactor=v_sin_k.*EnFactor;

OmGrad(1)= -8*g^2/(2*pi)^3* Simp3D(sin_kx.*uFactor,N,dk)  + 2*(gamma-2*mu)*u(1);
OmGrad(2)= -8*g^2/(2*pi)^3* Simp3D(sin_ky.*uFactor,N,dk) + 2*(gamma-2*mu)*u(2);
OmGrad(3)= -8*g^2/(2*pi)^3*gamma_z* Simp3D(sin_kz.*uFactor,N,dk) + 2*(gamma-2*mu + delta_z)*u(3);

OmGrad(4)= -8*g^2/(2*pi)^3*Simp3D(sin_kx.*vFactor,N,dk)  + 2*(gamma-2*mu)*v(1);
OmGrad(5)= -8*g^2/(2*pi)^3*Simp3D(sin_ky.*vFactor,N,dk)  + 2*(gamma-2*mu)*v(2);
OmGrad(6)= -8*g^2/(2*pi)^3*gamma_z* Simp3D(sin_kz.*vFactor,N,dk)  + 2*(gamma-2*mu + delta_z)*v(3);
%----------- Gradient of Om: OmGrad -----------% 


