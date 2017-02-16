function f = Msingleecho( omiga, azx, t )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if isempty(azx)
    f=ones(1,length(t));
    return;
else
omiga1 = realsqrt(realpow(omiga+azx(:,1),2)+realpow(azx(:,2),2));
f = 1-realpow((azx(:,2)./(2*omiga1)),2)*(1-cos(omiga*t)).*(1-cos(omiga1*t));
end
end
