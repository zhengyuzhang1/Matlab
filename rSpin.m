function [ U ] = rSpin( p )
% rotation expm(-1i*p*sigma)

a = norm(p);
if a ~= 0
    p = p/a;
end
c = cos(a/2); s = sin(a/2);
U = [c-1i*p(3)*s (-1i*p(1)-p(2))*s;(-1i*p(1)+p(2))*s c+1i*p(3)*s];

end

