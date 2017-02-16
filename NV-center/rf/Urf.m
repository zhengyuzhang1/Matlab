function U = Urf(t, trf, phi, w, wrf, w1, rabi, b, electronState)
%in lab frame;
%b is the direction of rf magnetic field

b = b/norm(b);

if electronState == 0
    r = sqrt(1-b(3)^2);
    cth = b(1)/r;
    sth = b(2)/r;
    bx = cth * cos(phi) - sth * sin(phi);
    by = sth * cos(phi) + cth * sin(phi);
    U = rSpin([0 0 1]*((t-trf)*w+trf*wrf))...
        *rSpin(trf*[rabi*r*bx,rabi*r*by,w-wrf]);
else
    szp = [w1(3);0;w1(2)];
    [R,~] = AxelRot(phi/pi*180,[w1(3),0,w1(2)]);
    sxp = R * [(w1(2)*b(1)-w1(3)*b(3))*w1(2); b(2); -(w1(2)*b(1)-w1(3)*b(3))*w1(3)]; 
    U = rSpin(szp*((t-trf)*w1(1)+trf*wrf))...
        *rSpin(trf*((w1(1)-wrf)*szp+rabi*sxp));

end
    