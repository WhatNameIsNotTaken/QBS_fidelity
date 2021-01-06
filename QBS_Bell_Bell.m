clear;
% genhao2=sym(sqrt(0.5));

L0=[1;0];
L1=[0;1];

L00=kron(L0,L0);
L01=kron(L0,L1);
L10=kron(L1,L0);
L11=kron(L1,L1);


%% initialization
phi_in=L00;
bell=sqrt(0.5)*(L00+L11);
rho_in=phi_in*phi_in';
% rho_source=kron(rho_in,kron(zeta*zeta',bell*bell'));
rho_source=kron(phi_in,kron(bell,bell))*kron(phi_in,kron(bell,bell))';

%% quantum gate

%CNOT_BD
CNOT_BD=eye(size(rho_source));
CNOT_BD([5,6,7,8,13,14,15,16],[5,6,7,8,13,14,15,16])=kron(eye(4),[0 1;1 0]);
t=[5,6,7,8,13,14,15,16]+16*ones(1,8);
CNOT_BD(t,t)=kron(eye(4),[0 1;1 0]);
t=t+32*ones(1,8);
CNOT_BD(t,t)=kron(eye(4),[0 1;1 0]);
t=t-16*ones(1,8);
CNOT_BD(t,t)=kron(eye(4),[0 1;1 0]);

%CNOT_AB
CNOT_AB=eye(size(rho_source));
CNOT_AB([9,13,10,14,11,15,12,16],[9,13,10,14,11,15,12,16])=kron(eye(4),[0 1;1 0]);
t=[9,13,10,14,11,15,12,16]+16*ones(1,8);
CNOT_AB(t,t)=kron(eye(4),[0 1;1 0]);
t=t+32*ones(1,8);
CNOT_AB(t,t)=kron(eye(4),[0 1;1 0]);
t=t-16*ones(1,8);
CNOT_AB(t,t)=kron(eye(4),[0 1;1 0]);

%H_A
H=sqrt(0.5)*[1 1;1 -1];
H_A=eye(size(rho_source));
t1=(1:8);
t2=(0:3)*2^4;
t1=[t1+t2(1),t1+t2(2),t1+t2(3),t1+t2(4)];
tt=[t1;t1+8];
t=reshape(tt,1,64);
H_A(t,t)=kron(eye(32),H);

%H_D
H_D=eye(size(rho_source));
t1=(0:31)*2+1;
t2=t1+1;
tt=[t1;t2];
t=reshape(tt,1,64);
H_D(t,t)=kron(eye(32),H);

%CNOT_DC
CNOT_DC=eye(size(rho_source));
CNOT_DC([2,4,6,8,10,12,14,16],[2,4,6,8,10,12,14,16])=kron(eye(4),[0 1;1 0]);
t=[2,4,6,8,10,12,14,16]+16*ones(1,8);
CNOT_DC(t,t)=kron(eye(4),[0 1;1 0]);
t=t+32*ones(1,8);
CNOT_DC(t,t)=kron(eye(4),[0 1;1 0]);
t=t-16*ones(1,8);
CNOT_DC(t,t)=kron(eye(4),[0 1;1 0]);

%H_B
H_B=eye(size(rho_source));
t1=(0:3)+1;
    % t2=t1+8;
t3=(0:7)*2^3;
t1=[t1+t3(1),t1+t3(2),t1+t3(3),t1+t3(4),t1+t3(5),t1+t3(6),t1+t3(7),t1+t3(8)];
t2=t1+4;
tt=[t1;t2];
t=reshape(tt,1,64);
H_B(t,t)=kron(eye(32),H);

%CNOT_1D
CNOT_1D=eye(size(rho_source));
t1=(0:15)*2+32+1;
t2=t1+1;
tt=[t1;t2];
t=reshape(tt,1,32);
CNOT_1D(t,t)=kron(eye(16),[0 1;1 0]);

%CNOT_1A
CNOT_1A=eye(size(rho_source));
t1=(0:7)+32+1;
t2=t1+8;
t3=t1+16;
t4=t3+8;
tt=[t1,t3;t2,t4];
t=reshape(tt,1,32);
CNOT_1A(t,t)=kron(eye(16),[0 1;1 0]);

%CZ_2A
CZ_2A=eye(size(rho_source));
t1=(0:7)+16+1;
t2=t1+8;
t3=t1+32;
t4=t3+8;
tt=[t1,t3;t2,t4];
t=reshape(tt,1,32);
CZ_2A(t,t)=kron(eye(16),[1 0;0 -1]);

%% output
U=CNOT_BD*CNOT_AB*H_A*H_D*CNOT_DC*H_B*CNOT_1D*CNOT_1A*CZ_2A;
rho_source_noise=rho_source;
rho_all=U*rho_source_noise*U';
t=(1:4);
rho_out=rho_all(t,t);
while(t(end)<=2^6)
    t=t+4;
    if t(end)<=2^6
    rho_out=rho_out+rho_all(t,t);
    end
end
disp('rho_in =');disp(rho_in);
disp('rho_out =');disp(rho_out);
f=phi_in'*rho_in*phi_in;disp('f =');disp(f);
f_out=phi_in'*rho_out*phi_in;disp('f_out =');disp(f_out);