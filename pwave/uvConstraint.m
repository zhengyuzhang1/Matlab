function [c, ceq, gradc, gradceq]=uvConstraint(uv)
% This function contains the constraints for the order parameters u and v
% c: u^2+v^2-1<=0;  ceq: u.*v==0; 
% gradc and gradceq are the gradients of the constraints c and ceq

u=uv(1:3); v=uv(4:6);

c=sum(u.^2 + v.^2) -1; % nonlinear inequality
ceq= sum(u.*v); % nonlinear equality, gauge choice

gradc=2*uv;

gradceq=zeros(length(uv),1);
gradceq(1:3)=uv(4:6); gradceq(4:6)=uv(1:3);

end