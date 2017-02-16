function dy = rigid(t,y,a,b,c)
dy = zeros(2,1);
dy(1) = a * y(1) + b * exp(c * t^2) * y(2);
dy(2) = b *exp(c * t^2) * y(1) - a * y(1);
end