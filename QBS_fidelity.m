%% It works correctly, only when QBS_Bell_GHZ.m ran.
% rho_source, phi_in calculated in the file QBS_Bell_GHZ.m
% toc= 292.4903, 4min 52s

close all;
tic;
p=0:0.05:1;

%% noise
%% phase flipping noise
E0=zeros(2,2,size(p,2));
E1=zeros(2,2,size(p,2));
for i=1:size(p,2)
    E0(:,:,i)=sqrt(1-p(i))*[1 0;0 1];
    E1(:,:,i)=sqrt(p(i))*[1 0;0 -1];
    e1=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E0(:,:,i),E0(:,:,i),E0(:,:,i));e2=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E0(:,:,i),E0(:,:,i),E1(:,:,i));
    e3=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E0(:,:,i),E1(:,:,i),E0(:,:,i));e4=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E0(:,:,i),E1(:,:,i),E1(:,:,i));
    e5=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E1(:,:,i),E0(:,:,i),E0(:,:,i));e6=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E1(:,:,i),E0(:,:,i),E1(:,:,i));
    e7=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E1(:,:,i),E1(:,:,i),E0(:,:,i));e8=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E1(:,:,i),E1(:,:,i),E1(:,:,i));
    e9=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E0(:,:,i),E0(:,:,i),E0(:,:,i));e10=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E0(:,:,i),E0(:,:,i),E1(:,:,i));
    e11=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E0(:,:,i),E1(:,:,i),E0(:,:,i));e12=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E0(:,:,i),E1(:,:,i),E1(:,:,i));
    e13=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E1(:,:,i),E0(:,:,i),E0(:,:,i));e14=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E1(:,:,i),E0(:,:,i),E1(:,:,i));
    e15=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E1(:,:,i),E1(:,:,i),E0(:,:,i));e16=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E1(:,:,i),E1(:,:,i),E1(:,:,i));

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
    rho_out_phaseflip=rho_out_phaseflip+rho_all_phaseflip(t,t);
    end
end
f_out_phaseflip(i)=phi_in'*rho_out_phaseflip*phi_in;
end

%% amplitude damping noise
E0=zeros(2,2,size(p,2));
E1=zeros(2,2,size(p,2));
for i=1:size(p,2)
    E0(:,:,i)=[1 0;0 sqrt(1-p(i))];
    E1(:,:,i)=[0 sqrt(p(i));0 0];
    e1=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E0(:,:,i),E0(:,:,i),E0(:,:,i));e2=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E0(:,:,i),E0(:,:,i),E1(:,:,i));
    e3=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E0(:,:,i),E1(:,:,i),E0(:,:,i));e4=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E0(:,:,i),E1(:,:,i),E1(:,:,i));
    e5=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E1(:,:,i),E0(:,:,i),E0(:,:,i));e6=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E1(:,:,i),E0(:,:,i),E1(:,:,i));
    e7=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E1(:,:,i),E1(:,:,i),E0(:,:,i));e8=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E1(:,:,i),E1(:,:,i),E1(:,:,i));
    e9=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E0(:,:,i),E0(:,:,i),E0(:,:,i));e10=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E0(:,:,i),E0(:,:,i),E1(:,:,i));
    e11=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E0(:,:,i),E1(:,:,i),E0(:,:,i));e12=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E0(:,:,i),E1(:,:,i),E1(:,:,i));
    e13=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E1(:,:,i),E0(:,:,i),E0(:,:,i));e14=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E1(:,:,i),E0(:,:,i),E1(:,:,i));
    e15=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E1(:,:,i),E1(:,:,i),E0(:,:,i));e16=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E1(:,:,i),E1(:,:,i),E1(:,:,i));

% output with amplitude damping noise
rho_source_ampdamp=simplify(e1*rho_source*e1'+e2*rho_source*e2'+e3*rho_source*e3'+e4*rho_source*e4'+e5*rho_source*e5'+e6*rho_source*e6'+...
    e7*rho_source*e7'+e8*rho_source*e8'+e9*rho_source*e9'+e10*rho_source*e10'+e11*rho_source*e11'+e12*rho_source*e12'+e13*rho_source*e13'+...
    e14*rho_source*e14'+e15*rho_source*e15'+e16*rho_source*e16');
rho_all_ampdamp=U*rho_source_ampdamp*U';
t=(1:4);
rho_out_ampdamp=rho_all_ampdamp(t,t);
while(t(end)<=2^7)
    t=t+4;
    if t(end)<=2^7
    rho_out_ampdamp=rho_out_ampdamp+rho_all_ampdamp(t,t);
    end
end
f_out_ampdamp(i)=phi_in'*rho_out_ampdamp*phi_in;
end

%% bit flipping noise
E0=zeros(2,2,size(p,2));
E1=zeros(2,2,size(p,2));
for i=1:size(p,2)
    E0(:,:,i)=sqrt(1-p(i))*[1 0;0 1];
    E1(:,:,i)=sqrt(p(i))*[0 1;1 0 ];
    e1 =mykron(eye(2),eye(2),E0(:,:,i),eye(2),E0(:,:,i),E0(:,:,i),E0(:,:,i));e2 =mykron(eye(2),eye(2),E0(:,:,i),eye(2),E0(:,:,i),E0(:,:,i),E1(:,:,i));
    e3 =mykron(eye(2),eye(2),E0(:,:,i),eye(2),E0(:,:,i),E1(:,:,i),E0(:,:,i));e4 =mykron(eye(2),eye(2),E0(:,:,i),eye(2),E0(:,:,i),E1(:,:,i),E1(:,:,i));
    e5 =mykron(eye(2),eye(2),E0(:,:,i),eye(2),E1(:,:,i),E0(:,:,i),E0(:,:,i));e6 =mykron(eye(2),eye(2),E0(:,:,i),eye(2),E1(:,:,i),E0(:,:,i),E1(:,:,i));
    e7 =mykron(eye(2),eye(2),E0(:,:,i),eye(2),E1(:,:,i),E1(:,:,i),E0(:,:,i));e8 =mykron(eye(2),eye(2),E0(:,:,i),eye(2),E1(:,:,i),E1(:,:,i),E1(:,:,i));
    e9 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E0(:,:,i),E0(:,:,i),E0(:,:,i));e10=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E0(:,:,i),E0(:,:,i),E1(:,:,i));
    e11=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E0(:,:,i),E1(:,:,i),E0(:,:,i));e12=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E0(:,:,i),E1(:,:,i),E1(:,:,i));
    e13=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E1(:,:,i),E0(:,:,i),E0(:,:,i));e14=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E1(:,:,i),E0(:,:,i),E1(:,:,i));
    e15=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E1(:,:,i),E1(:,:,i),E0(:,:,i));e16=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E1(:,:,i),E1(:,:,i),E1(:,:,i));

% output with bit-flipping noise
rho_source_bitflip=simplify(e1*rho_source*e1'+e2*rho_source*e2'+e3*rho_source*e3'+e4*rho_source*e4'+e5*rho_source*e5'+e6*rho_source*e6'+...
    e7*rho_source*e7'+e8*rho_source*e8'+e9*rho_source*e9'+e10*rho_source*e10'+e11*rho_source*e11'+e12*rho_source*e12'+e13*rho_source*e13'+...
    e14*rho_source*e14'+e15*rho_source*e15'+e16*rho_source*e16');
rho_all_bitflip=U*rho_source_bitflip*U';
t=(1:4);
rho_out_bitflip=rho_all_bitflip(t,t);
while(t(end)<=2^7)
    t=t+4;
    if t(end)<=2^7
    rho_out_bitflip=rho_out_bitflip+rho_all_bitflip(t,t);
    end
end
f_out_bitflip(i)=phi_in'*rho_out_bitflip*phi_in;
end
disp(toc);

%% plot
plot(p,f_out_phaseflip)
hold on;
[f_out_phaseflip_min,f_out_phaseflip_I]=min(f_out_phaseflip);                
plot(p(f_out_phaseflip_I),f_out_phaseflip(f_out_phaseflip_I),'.','MarkerSize',10)     
str = ['(' num2str(p(f_out_phaseflip_I)) ',' num2str(double(f_out_phaseflip(f_out_phaseflip_I))) ')'];
text(p(f_out_phaseflip_I),f_out_phaseflip(f_out_phaseflip_I),str)                       
hold on;

plot(p,f_out_ampdamp)
hold on;
[f_out_ampdamp_min,f_out_ampdamp_I]=min(f_out_ampdamp); 
plot(p(f_out_ampdamp_I),f_out_ampdamp(f_out_ampdamp_I),'.','MarkerSize',10)     
str = ['(' num2str(p(f_out_ampdamp_I)) ',' num2str(double(f_out_ampdamp(f_out_ampdamp_I))) ')'];
text(p(f_out_ampdamp_I),f_out_ampdamp(f_out_ampdamp_I),str)                       
hold on;

plot(p,f_out_bitflip)
hold on;
[f_out_bitflip_min,f_out_bitflip_I]=min(f_out_bitflip); 
plot(p(f_out_bitflip_I),f_out_bitflip(f_out_bitflip_I),'.','MarkerSize',10)     
str = ['(' num2str(p(f_out_bitflip_I)) ',' num2str(double(f_out_bitflip(f_out_bitflip_I))) ')'];
text(p(f_out_bitflip_I),f_out_bitflip(f_out_bitflip_I),str)                       
hold on;