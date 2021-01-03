%% It works correctly, only when QBS_Bell_GHZ.m ran.

syms p;
assume(p>0&p<1);

%% noise
%% phase flipping noise
E0=sqrt(p)*[1 0;0 1];
E1=sqrt(1-p)*[1 0;0 -1];
e1 =mykron(eye(2),eye(2),E0,eye(2),E0,E0,E0);e2 =mykron(eye(2),eye(2),E0,eye(2),E0,E0,E1);
e3 =mykron(eye(2),eye(2),E0,eye(2),E0,E1,E0);e4 =mykron(eye(2),eye(2),E0,eye(2),E0,E1,E1);
e5 =mykron(eye(2),eye(2),E0,eye(2),E1,E0,E0);e6 =mykron(eye(2),eye(2),E0,eye(2),E1,E0,E1);
e7 =mykron(eye(2),eye(2),E0,eye(2),E1,E1,E0);e8 =mykron(eye(2),eye(2),E0,eye(2),E1,E1,E1);
e9 =mykron(eye(2),eye(2),E1,eye(2),E0,E0,E0);e10=mykron(eye(2),eye(2),E1,eye(2),E0,E0,E1);
e11=mykron(eye(2),eye(2),E1,eye(2),E0,E1,E0);e12=mykron(eye(2),eye(2),E1,eye(2),E0,E1,E1);
e13=mykron(eye(2),eye(2),E1,eye(2),E1,E0,E0);e14=mykron(eye(2),eye(2),E1,eye(2),E1,E0,E1);
e15=mykron(eye(2),eye(2),E1,eye(2),E1,E1,E0);e16=mykron(eye(2),eye(2),E1,eye(2),E1,E1,E1);

% output with phase-flipping noise
rho_source_phaseflip=simplify(e1*rho_source*e1'+e2*rho_source*e2'+e3*rho_source*e3'+e4*rho_source*e4'+e5*rho_source*e5'+e6*rho_source*e6'+...
    e7*rho_source*e7'+e8*rho_source*e8'+e9*rho_source*e9'+e10*rho_source*e10'+e11*rho_source*e11'+e12*rho_source*e12'+e13*rho_source*e13'+...
    e14*rho_source*e14'+e15*rho_source*e15'+e16*rho_source*e16');
rho_all_phaseflip=U*rho_source_phaseflip*U';
t=(1:4);
rho_out_phaseflip=rho_all_phaseflip(t,t);
while(t(end)<=2^7)
    t=t+4;
    if t(end)<=2^7
    rho_out_phaseflip=collect(rho_out_phaseflip+rho_all_phaseflip(t,t));
    end
end
disp(rho_out_phaseflip);
    
% calculate the fidelity
f_out_phaseflip=collect(phi_in'*rho_out_phaseflip*phi_in);disp('f_out_phaseflip =');disp(f_out_phaseflip);
% 4*p^3 - 6*p^2 + 3*p
