function f = M(w, n, azx, t )
% w is the Larmor frequency of C13 nucleus, n is number of Pi pulse, azx is
% a n*2 matrix of coupling constant in the form of [Az1,Ax1;Az2,Ax2;...],t
% is a row vector of time

if isempty(azx)
    f=ones(1,length(t));
    return;
else
    w1 = realsqrt(realpow(w+azx(:,1),2)+realpow(azx(:,2),2));
    coshalffi=bsxfun(@times,cos(w1*t),cos(w*t))-((azx(:,1)+w)./w1)*sin(w*t).*sin(w1*t);
    f = 1-realpow((azx(:,2)./w1),2)*(1-cos(w*t)).*(1-cos(w1*t))./(1+coshalffi).*realpow(sin(n/2*acos(coshalffi)),2);
end
end
