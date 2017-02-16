S_x = 1/2 * [0 1; 1 0];
S_y = 1/2 * [0 -1i; 1i 0];
S_z = [0 0; 0 -1];
I_x = 1/2 * [0 1; 1 0];
I_y = 1/2 * [0 -1i; 1i 0];
I = [1 0; 0 1];

Initial_state = [1,0,0,0]';
E_Rabi = 300;Bz = 505; 
 T_final1 = 3.125;   tau1 = 9.2507;
 T_final2 = 0.96;     tau2 = 7.0803;%7.0803;
 T_final3 = 2.08;     tau3 = 7.2308;
A1= 2*pi*(-2.162); I1_z = [1 0; 0 0];         gamma1 = 2*pi*(-0.000308);  Omega1 = 2*pi*0.04; 
A2= 2*pi*13.7;      I2_z = 1/2*[1 0; 0 -1]; gamma2 = 2*pi*(-0.00107);    Omega2 = 2*pi*0.13; 
A3= 2*pi*(-6.5);     I3_z = 1/2*[1 0; 0 -1]; gamma3 = 2*pi*(-0.00107);    Omega3 = 2*pi*0.06;

A = A1; I_z = I1_z; Omega = Omega1; T_final = T_final1; tau = tau1;
H_with_RF =   A*kron(I,I_z) +A*kron(S_z,I_z) + Omega*kron(I,I_x);% -A*kron(S_z,I) +
H_without_RF =  A*kron(I,I_z) +A*kron(S_z,I_z);

RotationX = expm(-1i*(H_without_RF+2*pi*E_Rabi*kron(S_x,I))*0.5/E_Rabi);
Hadamard_E = expm(-1i*(H_without_RF+2*pi*E_Rabi*kron(S_x,I))*0.25/E_Rabi);
Final_stateH = RotationX*Initial_state;
 
Evolution1 = expm(-1i*H_without_RF*(tau-T_final));
Evolution2 = expm(-1i*H_without_RF*2*(tau-T_final));

U = Evolution1*expm(-1i*H_with_RF*T_final)*RotationX*Evolution2*...
    expm(-1i*H_with_RF*2*T_final)*RotationX*Evolution1*expm(-1i*H_with_RF*T_final)...
    *Evolution1*expm(-1i*H_with_RF*T_final)*RotationX*Evolution2*...
    expm(-1i*H_with_RF*2*T_final)*RotationX*Evolution1*expm(-1i*H_with_RF*T_final);
S1 = expm(-1i*H_with_RF*T_final)*Initial_state;
S2 = RotationX*Evolution1*S1;
S3 = expm(-1i*H_with_RF*2*T_final)*S2;
S4 = RotationX*Evolution2*S3;
S5 = expm(-1i*H_with_RF*T_final)*S4;
S6 = U*Initial_state;
ss = dot(S6,S6);

U_target = kron(I,2*S_x);
Utria = @(x)kron( conj(x) , x );
Fidelity_N = real( trace( ctranspose(Utria(U)) * Utria(U_target) ) )/16 

 