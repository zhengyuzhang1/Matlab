function [axis,angle] = axisangle(U)

globalprefactor = det(U)^(1/4);
U = U/globalprefactor;
if real(U(1))<0
    U = -U;
end
angle = acos(abs(trace(U))/2)*2;

axis = zeros(3,1);
axis(1) = -imag(U(2))/sin(angle/2);
axis(2) = -real(U(2))/sin(angle/2);
axis(3) = -imag(U(1))/sin(angle/2);

end
