close all;clear;clc;
genhao2=sym(sqrt(0.5));
%% 
L0=[1;0];
L1=[0;1];

L00=kron(L0,L0);
L01=kron(L0,L1);
L10=kron(L1,L0);
L11=kron(L1,L1);

L000=kron(L0,L00);
L001=kron(L0,L01);
L010=kron(L0,L10);
L011=kron(L0,L11);
L100=kron(L1,L00);
L101=kron(L1,L01);
L110=kron(L1,L10);
L111=kron(L1,L11);

%% initialization
phi_in=L10;
zeta=genhao2*(L000+L111);
bell=genhao2*(L00+L11);
rho_in=phi_in*phi_in';
% rho_source=kron(rho_in,kron(zeta*zeta',bell*bell'));
rho_source=kron(phi_in,kron(bell,zeta))*kron(phi_in,kron(bell,zeta))';

%% quantum gate

%CNOT_BE
CNOT_BE=eye(size(rho_source));
CNOT_BE([9,10,11,12,13,14,15,16,25,26,27,28,29,30,31,32],[9,10,11,12,13,14,15,16,25,26,27,28,29,30,31,32])=kron(eye(8),[0 1;1 0]);
t=[9,10,11,12,13,14,15,16,25,26,27,28,29,30,31,32]+32*ones(1,16);
CNOT_BE(t,t)=kron(eye(8),[0 1;1 0]);
t=t+64*ones(1,16);
CNOT_BE(t,t)=kron(eye(8),[0 1;1 0]);
t=t-32*ones(1,16);
CNOT_BE(t,t)=kron(eye(8),[0 1;1 0]);

%CNOT_AB
CNOT_AB=eye(size(rho_source));
CNOT_AB([17,25,18,26,19,27,20,28,21,29,22,30,23,31,24,32],[17,25,18,26,19,27,20,28,21,29,22,30,23,31,24,32])=kron(eye(8),[0 1;1 0]);
t=[17,25,18,26,19,27,20,28,21,29,22,30,23,31,24,32]+32*ones(1,16);
CNOT_AB(t,t)=kron(eye(8),[0 1;1 0]);
t=t+64*ones(1,16);
CNOT_AB(t,t)=kron(eye(8),[0 1;1 0]);
t=t-32*ones(1,16);
CNOT_AB(t,t)=kron(eye(8),[0 1;1 0]);

%H_A
H=genhao2*[1 1;1 -1];
H_A=eye(size(rho_source));
t1=(0:15)+1;
t2=t1+16;
t3=t1+32;
t4=t2+32;
t5=t1+64;
t6=t2+64;
t7=t1+96;
t8=t2+96;
tt=[t1,t3,t5,t7;t2,t4,t6,t8];
t=reshape(tt,1,128);
H_A(t,t)=kron(eye(64),H);

%H_E
H_E=eye(size(rho_source));
t1=(0:63)*2+1;
t2=t1+1;
tt=[t1;t2];
t=reshape(tt,1,128);
H_E(t,t)=kron(eye(64),H);

%CNOT_ED
CNOT_ED=eye(size(rho_source));
CNOT_ED([2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32],[2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32])=kron(eye(8),[0 1;1 0]);
t=[2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32]+32*ones(1,16);
CNOT_ED(t,t)=kron(eye(8),[0 1;1 0]);
t=t+64*ones(1,16);
CNOT_ED(t,t)=kron(eye(8),[0 1;1 0]);
t=t-32*ones(1,16);
CNOT_ED(t,t)=kron(eye(8),[0 1;1 0]);

%H_B
H_B=eye(size(rho_source));
t1=(0:7)+1;
    % t2=t1+8;
t3=(0:7)*2^4;
t1=[t1+t3(1),t1+t3(2),t1+t3(3),t1+t3(4),t1+t3(5),t1+t3(6),t1+t3(7),t1+t3(8)];
t2=t1+8;
tt=[t1;t2];
t=reshape(tt,1,128);
H_B(t,t)=kron(eye(64),H);

%CNOT_BC
CNOT_BC=eye(size(rho_source));
CNOT_BC([9,13,10,14,11,15,12,16,25,29,26,30,27,31,28,32],[9,13,10,14,11,15,12,16,25,29,26,30,27,31,28,32])=kron(eye(8),[0 1;1 0]);
t=[9,13,10,14,11,15,12,16,25,29,26,30,27,31,28,32]+32*ones(1,16);
CNOT_BC(t,t)=kron(eye(8),[0 1;1 0]);
t=t+64*ones(1,16);
CNOT_BC(t,t)=kron(eye(8),[0 1;1 0]);
t=t-32*ones(1,16);
CNOT_BC(t,t)=kron(eye(8),[0 1;1 0]);

%CNOT_1E
CNOT_1E=eye(size(rho_source));
t1=(0:31)*2+64+1;
t2=t1+1;
tt=[t1;t2];
t=reshape(tt,1,64);
CNOT_1E(t,t)=kron(eye(32),[0 1;1 0]);

%CNOT_1A
CNOT_1A=eye(size(rho_source));
t1=(0:15)+64+1;
t2=t1+16;
t3=t1+32;
t4=t2+32;
tt=[t1,t3;t2,t4];
t=reshape(tt,1,64);
CNOT_1A(t,t)=kron(eye(32),[0 1;1 0]);

%CZ_2A
CZ_2A=eye(size(rho_source));
t1=(0:15)+32+1;
t2=t1+16;
t3=t1+64;
t4=t2+64;
tt=[t1,t3;t2,t4];
t=reshape(tt,1,64);
CZ_2A(t,t)=kron(eye(32),[1 0;0 -1]);

%% output
U=CNOT_BE*CNOT_AB*H_A*H_E*CNOT_ED*H_B*CNOT_BC*CNOT_1E*CNOT_1A*CZ_2A;
rho_source_noise=rho_source;
rho_all=U*rho_source_noise*U';
t=(1:4);
rho_out=rho_all(t,t);
while(t(end)<=2^7)
    t=t+4;
    if t(end)<=2^7
    rho_out=rho_out+rho_all(t,t);
    end
end
disp('rho_in =');disp(rho_in);
disp('rho_out =');disp(rho_out);
f=phi_in'*rho_in*phi_in;disp('f =');disp(f);
f_out=phi_in'*rho_out*phi_in;disp('f_out =');disp(f_out);