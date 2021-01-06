%% It works correctly, only when QBS_Bell_GHZ.m ran.
% It is not necessary to divide the file into two parts. It is the 'sym' parameter caused the slow course.
% rho_source, phi_in calculated in the file QBS_Bell_GHZ.m
% toc= 3.7178

tic;
p=0:0.05:1;

%% noise
%% x noise
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

% output with x noise
rho_source_x=e1*rho_source*e1'+e2*rho_source*e2'+e3*rho_source*e3'+e4*rho_source*e4'+e5*rho_source*e5'+e6*rho_source*e6'+...
    e7*rho_source*e7'+e8*rho_source*e8'+e9*rho_source*e9'+e10*rho_source*e10'+e11*rho_source*e11'+e12*rho_source*e12'+e13*rho_source*e13'+...
    e14*rho_source*e14'+e15*rho_source*e15'+e16*rho_source*e16';
rho_all_x=U*rho_source_x*U';
t=(1:4);
rho_out_x=rho_all_x(t,t);
while(t(end)<=2^7)
    t=t+4;
    if t(end)<=2^7
    rho_out_x=rho_out_x+rho_all_x(t,t);
    end
end
f_out_x(i)=phi_in'*rho_out_x*phi_in;
end

%% y noise
E0=zeros(2,2,size(p,2));
E1=zeros(2,2,size(p,2));
for i=1:size(p,2)
    E0(:,:,i)=sqrt(1-p(i))*[1 0;0 1];
    E1(:,:,i)=sqrt(p(i))*[0 -1;1 0];
    e1=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E0(:,:,i),E0(:,:,i),E0(:,:,i));e2=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E0(:,:,i),E0(:,:,i),E1(:,:,i));
    e3=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E0(:,:,i),E1(:,:,i),E0(:,:,i));e4=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E0(:,:,i),E1(:,:,i),E1(:,:,i));
    e5=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E1(:,:,i),E0(:,:,i),E0(:,:,i));e6=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E1(:,:,i),E0(:,:,i),E1(:,:,i));
    e7=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E1(:,:,i),E1(:,:,i),E0(:,:,i));e8=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E1(:,:,i),E1(:,:,i),E1(:,:,i));
    e9=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E0(:,:,i),E0(:,:,i),E0(:,:,i));e10=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E0(:,:,i),E0(:,:,i),E1(:,:,i));
    e11=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E0(:,:,i),E1(:,:,i),E0(:,:,i));e12=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E0(:,:,i),E1(:,:,i),E1(:,:,i));
    e13=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E1(:,:,i),E0(:,:,i),E0(:,:,i));e14=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E1(:,:,i),E0(:,:,i),E1(:,:,i));
    e15=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E1(:,:,i),E1(:,:,i),E0(:,:,i));e16=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E1(:,:,i),E1(:,:,i),E1(:,:,i));

% output with y noise
rho_source_y=e1*rho_source*e1'+e2*rho_source*e2'+e3*rho_source*e3'+e4*rho_source*e4'+e5*rho_source*e5'+e6*rho_source*e6'+...
    e7*rho_source*e7'+e8*rho_source*e8'+e9*rho_source*e9'+e10*rho_source*e10'+e11*rho_source*e11'+e12*rho_source*e12'+e13*rho_source*e13'+...
    e14*rho_source*e14'+e15*rho_source*e15'+e16*rho_source*e16';
rho_all_y=U*rho_source_y*U';
t=(1:4);
rho_out_y=rho_all_y(t,t);
while(t(end)<=2^7)
    t=t+4;
    if t(end)<=2^7
    rho_out_y=rho_out_y+rho_all_y(t,t);
    end
end
f_out_y(i)=phi_in'*rho_out_y*phi_in;
end

%% z noise
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

% output with z noise
rho_source_z=e1*rho_source*e1'+e2*rho_source*e2'+e3*rho_source*e3'+e4*rho_source*e4'+e5*rho_source*e5'+e6*rho_source*e6'+...
    e7*rho_source*e7'+e8*rho_source*e8'+e9*rho_source*e9'+e10*rho_source*e10'+e11*rho_source*e11'+e12*rho_source*e12'+e13*rho_source*e13'+...
    e14*rho_source*e14'+e15*rho_source*e15'+e16*rho_source*e16';
rho_all_z=U*rho_source_z*U';
t=(1:4);
rho_out_z=rho_all_z(t,t);
while(t(end)<=2^7)
    t=t+4;
    if t(end)<=2^7
    rho_out_z=rho_out_z+rho_all_z(t,t);
    end
end
f_out_z(i)=phi_in'*rho_out_z*phi_in;
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
rho_source_ampdamp=e1*rho_source*e1'+e2*rho_source*e2'+e3*rho_source*e3'+e4*rho_source*e4'+e5*rho_source*e5'+e6*rho_source*e6'+...
    e7*rho_source*e7'+e8*rho_source*e8'+e9*rho_source*e9'+e10*rho_source*e10'+e11*rho_source*e11'+e12*rho_source*e12'+e13*rho_source*e13'+...
    e14*rho_source*e14'+e15*rho_source*e15'+e16*rho_source*e16';
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

%% depolarizing noise
E0=zeros(2,2,size(p,2));
E1=zeros(2,2,size(p,2));
E2=zeros(2,2,size(p,2));
E3=zeros(2,2,size(p,2));
for i=1:size(p,2)
E0(:,:,i)=sqrt(1-p(i))*[1 0;0 1];
E1(:,:,i)=sqrt(p(i)/3)*[0 1;1 0];
E2(:,:,i)=sqrt(p(i)/3)*[0 -1i;1i 0];
E3(:,:,i)=sqrt(p(i)/3)*[1 0;0 -1];

e1 =mykron(eye(2),eye(2),E0(:,:,i),eye(2),E0(:,:,i),E0(:,:,i),E0(:,:,i)); e65 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E0(:,:,i),E0(:,:,i),E0(:,:,i)); e129=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E0(:,:,i),E0(:,:,i),E0(:,:,i)); e193=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E0(:,:,i),E0(:,:,i),E0(:,:,i));
e2 =mykron(eye(2),eye(2),E0(:,:,i),eye(2),E0(:,:,i),E0(:,:,i),E1(:,:,i)); e66 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E0(:,:,i),E0(:,:,i),E1(:,:,i)); e130=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E0(:,:,i),E0(:,:,i),E1(:,:,i)); e194=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E0(:,:,i),E0(:,:,i),E1(:,:,i));
e3 =mykron(eye(2),eye(2),E0(:,:,i),eye(2),E0(:,:,i),E0(:,:,i),E2(:,:,i)); e67 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E0(:,:,i),E0(:,:,i),E2(:,:,i)); e131=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E0(:,:,i),E0(:,:,i),E2(:,:,i)); e195=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E0(:,:,i),E0(:,:,i),E2(:,:,i));
e4 =mykron(eye(2),eye(2),E0(:,:,i),eye(2),E0(:,:,i),E0(:,:,i),E3(:,:,i)); e68 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E0(:,:,i),E0(:,:,i),E3(:,:,i)); e132=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E0(:,:,i),E0(:,:,i),E3(:,:,i)); e196=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E0(:,:,i),E0(:,:,i),E3(:,:,i));
e5 =mykron(eye(2),eye(2),E0(:,:,i),eye(2),E0(:,:,i),E1(:,:,i),E0(:,:,i)); e69 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E0(:,:,i),E1(:,:,i),E0(:,:,i)); e133=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E0(:,:,i),E1(:,:,i),E0(:,:,i)); e197=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E0(:,:,i),E1(:,:,i),E0(:,:,i));
e6 =mykron(eye(2),eye(2),E0(:,:,i),eye(2),E0(:,:,i),E1(:,:,i),E1(:,:,i)); e70 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E0(:,:,i),E1(:,:,i),E1(:,:,i)); e134=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E0(:,:,i),E1(:,:,i),E1(:,:,i)); e198=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E0(:,:,i),E1(:,:,i),E1(:,:,i));
e7 =mykron(eye(2),eye(2),E0(:,:,i),eye(2),E0(:,:,i),E1(:,:,i),E2(:,:,i)); e71 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E0(:,:,i),E1(:,:,i),E2(:,:,i)); e135=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E0(:,:,i),E1(:,:,i),E2(:,:,i)); e199=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E0(:,:,i),E1(:,:,i),E2(:,:,i));
e8 =mykron(eye(2),eye(2),E0(:,:,i),eye(2),E0(:,:,i),E1(:,:,i),E3(:,:,i)); e72 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E0(:,:,i),E1(:,:,i),E3(:,:,i)); e136=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E0(:,:,i),E1(:,:,i),E3(:,:,i)); e200=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E0(:,:,i),E1(:,:,i),E3(:,:,i));
e9 =mykron(eye(2),eye(2),E0(:,:,i),eye(2),E0(:,:,i),E2(:,:,i),E0(:,:,i)); e73 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E0(:,:,i),E2(:,:,i),E0(:,:,i)); e137=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E0(:,:,i),E2(:,:,i),E0(:,:,i)); e201=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E0(:,:,i),E2(:,:,i),E0(:,:,i));
e10=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E0(:,:,i),E2(:,:,i),E1(:,:,i)); e74 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E0(:,:,i),E2(:,:,i),E1(:,:,i)); e138=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E0(:,:,i),E2(:,:,i),E1(:,:,i)); e202=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E0(:,:,i),E2(:,:,i),E1(:,:,i));
e11=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E0(:,:,i),E2(:,:,i),E2(:,:,i)); e75 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E0(:,:,i),E2(:,:,i),E2(:,:,i)); e139=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E0(:,:,i),E2(:,:,i),E2(:,:,i)); e203=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E0(:,:,i),E2(:,:,i),E2(:,:,i));
e12=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E0(:,:,i),E2(:,:,i),E3(:,:,i)); e76 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E0(:,:,i),E2(:,:,i),E3(:,:,i)); e140=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E0(:,:,i),E2(:,:,i),E3(:,:,i)); e204=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E0(:,:,i),E2(:,:,i),E3(:,:,i));
e13=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E0(:,:,i),E3(:,:,i),E0(:,:,i)); e77 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E0(:,:,i),E3(:,:,i),E0(:,:,i)); e141=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E0(:,:,i),E3(:,:,i),E0(:,:,i)); e205=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E0(:,:,i),E3(:,:,i),E0(:,:,i));
e14=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E0(:,:,i),E3(:,:,i),E1(:,:,i)); e78 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E0(:,:,i),E3(:,:,i),E1(:,:,i)); e142=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E0(:,:,i),E3(:,:,i),E1(:,:,i)); e206=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E0(:,:,i),E3(:,:,i),E1(:,:,i));
e15=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E0(:,:,i),E3(:,:,i),E2(:,:,i)); e79 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E0(:,:,i),E3(:,:,i),E2(:,:,i)); e143=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E0(:,:,i),E3(:,:,i),E2(:,:,i)); e207=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E0(:,:,i),E3(:,:,i),E2(:,:,i));
e16=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E0(:,:,i),E3(:,:,i),E3(:,:,i)); e80 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E0(:,:,i),E3(:,:,i),E3(:,:,i)); e144=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E0(:,:,i),E3(:,:,i),E3(:,:,i)); e208=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E0(:,:,i),E3(:,:,i),E3(:,:,i));
e17=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E1(:,:,i),E0(:,:,i),E0(:,:,i)); e81 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E1(:,:,i),E0(:,:,i),E0(:,:,i)); e145=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E1(:,:,i),E0(:,:,i),E0(:,:,i)); e209=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E1(:,:,i),E0(:,:,i),E0(:,:,i));
e18=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E1(:,:,i),E0(:,:,i),E1(:,:,i)); e82 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E1(:,:,i),E0(:,:,i),E1(:,:,i)); e146=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E1(:,:,i),E0(:,:,i),E1(:,:,i)); e210=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E1(:,:,i),E0(:,:,i),E1(:,:,i));
e19=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E1(:,:,i),E0(:,:,i),E2(:,:,i)); e83 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E1(:,:,i),E0(:,:,i),E2(:,:,i)); e147=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E1(:,:,i),E0(:,:,i),E2(:,:,i)); e211=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E1(:,:,i),E0(:,:,i),E2(:,:,i));
e20=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E1(:,:,i),E0(:,:,i),E3(:,:,i)); e84 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E1(:,:,i),E0(:,:,i),E3(:,:,i)); e148=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E1(:,:,i),E0(:,:,i),E3(:,:,i)); e212=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E1(:,:,i),E0(:,:,i),E3(:,:,i));
e21=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E1(:,:,i),E1(:,:,i),E0(:,:,i)); e85 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E1(:,:,i),E1(:,:,i),E0(:,:,i)); e149=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E1(:,:,i),E1(:,:,i),E0(:,:,i)); e213=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E1(:,:,i),E1(:,:,i),E0(:,:,i));
e22=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E1(:,:,i),E1(:,:,i),E1(:,:,i)); e86 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E1(:,:,i),E1(:,:,i),E1(:,:,i)); e150=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E1(:,:,i),E1(:,:,i),E1(:,:,i)); e214=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E1(:,:,i),E1(:,:,i),E1(:,:,i));
e23=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E1(:,:,i),E1(:,:,i),E2(:,:,i)); e87 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E1(:,:,i),E1(:,:,i),E2(:,:,i)); e151=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E1(:,:,i),E1(:,:,i),E2(:,:,i)); e215=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E1(:,:,i),E1(:,:,i),E2(:,:,i));
e24=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E1(:,:,i),E1(:,:,i),E3(:,:,i)); e88 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E1(:,:,i),E1(:,:,i),E3(:,:,i)); e152=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E1(:,:,i),E1(:,:,i),E3(:,:,i)); e216=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E1(:,:,i),E1(:,:,i),E3(:,:,i));
e25=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E1(:,:,i),E2(:,:,i),E0(:,:,i)); e89 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E1(:,:,i),E2(:,:,i),E0(:,:,i)); e153=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E1(:,:,i),E2(:,:,i),E0(:,:,i)); e217=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E1(:,:,i),E2(:,:,i),E0(:,:,i));
e26=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E1(:,:,i),E2(:,:,i),E1(:,:,i)); e90 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E1(:,:,i),E2(:,:,i),E1(:,:,i)); e154=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E1(:,:,i),E2(:,:,i),E1(:,:,i)); e218=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E1(:,:,i),E2(:,:,i),E1(:,:,i));
e27=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E1(:,:,i),E2(:,:,i),E2(:,:,i)); e91 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E1(:,:,i),E2(:,:,i),E2(:,:,i)); e155=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E1(:,:,i),E2(:,:,i),E2(:,:,i)); e219=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E1(:,:,i),E2(:,:,i),E2(:,:,i));
e28=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E1(:,:,i),E2(:,:,i),E3(:,:,i)); e92 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E1(:,:,i),E2(:,:,i),E3(:,:,i)); e156=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E1(:,:,i),E2(:,:,i),E3(:,:,i)); e220=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E1(:,:,i),E2(:,:,i),E3(:,:,i));
e29=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E1(:,:,i),E3(:,:,i),E0(:,:,i)); e93 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E1(:,:,i),E3(:,:,i),E0(:,:,i)); e157=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E1(:,:,i),E3(:,:,i),E0(:,:,i)); e221=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E1(:,:,i),E3(:,:,i),E0(:,:,i));
e30=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E1(:,:,i),E3(:,:,i),E1(:,:,i)); e94 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E1(:,:,i),E3(:,:,i),E1(:,:,i)); e158=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E1(:,:,i),E3(:,:,i),E1(:,:,i)); e222=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E1(:,:,i),E3(:,:,i),E1(:,:,i));
e31=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E1(:,:,i),E3(:,:,i),E2(:,:,i)); e95 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E1(:,:,i),E3(:,:,i),E2(:,:,i)); e159=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E1(:,:,i),E3(:,:,i),E2(:,:,i)); e223=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E1(:,:,i),E3(:,:,i),E2(:,:,i));
e32=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E1(:,:,i),E3(:,:,i),E3(:,:,i)); e96 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E1(:,:,i),E3(:,:,i),E3(:,:,i)); e160=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E1(:,:,i),E3(:,:,i),E3(:,:,i)); e224=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E1(:,:,i),E3(:,:,i),E3(:,:,i));
e33=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E2(:,:,i),E0(:,:,i),E0(:,:,i)); e97 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E2(:,:,i),E0(:,:,i),E0(:,:,i)); e161=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E2(:,:,i),E0(:,:,i),E0(:,:,i)); e225=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E2(:,:,i),E0(:,:,i),E0(:,:,i));
e34=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E2(:,:,i),E0(:,:,i),E1(:,:,i)); e98 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E2(:,:,i),E0(:,:,i),E1(:,:,i)); e162=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E2(:,:,i),E0(:,:,i),E1(:,:,i)); e226=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E2(:,:,i),E0(:,:,i),E1(:,:,i));
e35=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E2(:,:,i),E0(:,:,i),E2(:,:,i)); e99 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E2(:,:,i),E0(:,:,i),E2(:,:,i)); e163=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E2(:,:,i),E0(:,:,i),E2(:,:,i)); e227=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E2(:,:,i),E0(:,:,i),E2(:,:,i));
e36=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E2(:,:,i),E0(:,:,i),E3(:,:,i)); e100=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E2(:,:,i),E0(:,:,i),E3(:,:,i)); e164=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E2(:,:,i),E0(:,:,i),E3(:,:,i)); e228=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E2(:,:,i),E0(:,:,i),E3(:,:,i));
e37=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E2(:,:,i),E1(:,:,i),E0(:,:,i)); e101=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E2(:,:,i),E1(:,:,i),E0(:,:,i)); e165=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E2(:,:,i),E1(:,:,i),E0(:,:,i)); e229=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E2(:,:,i),E1(:,:,i),E0(:,:,i));
e38=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E2(:,:,i),E1(:,:,i),E1(:,:,i)); e102=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E2(:,:,i),E1(:,:,i),E1(:,:,i)); e166=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E2(:,:,i),E1(:,:,i),E1(:,:,i)); e230=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E2(:,:,i),E1(:,:,i),E1(:,:,i));
e39=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E2(:,:,i),E1(:,:,i),E2(:,:,i)); e103=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E2(:,:,i),E1(:,:,i),E2(:,:,i)); e167=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E2(:,:,i),E1(:,:,i),E2(:,:,i)); e231=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E2(:,:,i),E1(:,:,i),E2(:,:,i));
e40=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E2(:,:,i),E1(:,:,i),E3(:,:,i)); e104=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E2(:,:,i),E1(:,:,i),E3(:,:,i)); e168=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E2(:,:,i),E1(:,:,i),E3(:,:,i)); e232=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E2(:,:,i),E1(:,:,i),E3(:,:,i));
e41=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E2(:,:,i),E2(:,:,i),E0(:,:,i)); e105=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E2(:,:,i),E2(:,:,i),E0(:,:,i)); e169=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E2(:,:,i),E2(:,:,i),E0(:,:,i)); e233=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E2(:,:,i),E2(:,:,i),E0(:,:,i));
e42=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E2(:,:,i),E2(:,:,i),E1(:,:,i)); e106=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E2(:,:,i),E2(:,:,i),E1(:,:,i)); e170=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E2(:,:,i),E2(:,:,i),E1(:,:,i)); e234=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E2(:,:,i),E2(:,:,i),E1(:,:,i));
e43=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E2(:,:,i),E2(:,:,i),E2(:,:,i)); e107=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E2(:,:,i),E2(:,:,i),E2(:,:,i)); e171=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E2(:,:,i),E2(:,:,i),E2(:,:,i)); e235=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E2(:,:,i),E2(:,:,i),E2(:,:,i));
e44=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E2(:,:,i),E2(:,:,i),E3(:,:,i)); e108=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E2(:,:,i),E2(:,:,i),E3(:,:,i)); e172=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E2(:,:,i),E2(:,:,i),E3(:,:,i)); e236=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E2(:,:,i),E2(:,:,i),E3(:,:,i));
e45=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E2(:,:,i),E3(:,:,i),E0(:,:,i)); e109=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E2(:,:,i),E3(:,:,i),E0(:,:,i)); e173=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E2(:,:,i),E3(:,:,i),E0(:,:,i)); e237=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E2(:,:,i),E3(:,:,i),E0(:,:,i));
e46=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E2(:,:,i),E3(:,:,i),E1(:,:,i)); e110=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E2(:,:,i),E3(:,:,i),E1(:,:,i)); e174=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E2(:,:,i),E3(:,:,i),E1(:,:,i)); e238=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E2(:,:,i),E3(:,:,i),E1(:,:,i));
e47=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E2(:,:,i),E3(:,:,i),E2(:,:,i)); e111=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E2(:,:,i),E3(:,:,i),E2(:,:,i)); e175=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E2(:,:,i),E3(:,:,i),E2(:,:,i)); e239=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E2(:,:,i),E3(:,:,i),E2(:,:,i));
e48=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E2(:,:,i),E3(:,:,i),E3(:,:,i)); e112=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E2(:,:,i),E3(:,:,i),E3(:,:,i)); e176=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E2(:,:,i),E3(:,:,i),E3(:,:,i)); e240=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E2(:,:,i),E3(:,:,i),E3(:,:,i));
e49=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E3(:,:,i),E0(:,:,i),E0(:,:,i)); e113=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E3(:,:,i),E0(:,:,i),E0(:,:,i)); e177=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E3(:,:,i),E0(:,:,i),E0(:,:,i)); e241=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E3(:,:,i),E0(:,:,i),E0(:,:,i));
e50=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E3(:,:,i),E0(:,:,i),E1(:,:,i)); e114=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E3(:,:,i),E0(:,:,i),E1(:,:,i)); e178=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E3(:,:,i),E0(:,:,i),E1(:,:,i)); e242=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E3(:,:,i),E0(:,:,i),E1(:,:,i));
e51=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E3(:,:,i),E0(:,:,i),E2(:,:,i)); e115=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E3(:,:,i),E0(:,:,i),E2(:,:,i)); e179=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E3(:,:,i),E0(:,:,i),E2(:,:,i)); e243=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E3(:,:,i),E0(:,:,i),E2(:,:,i));
e52=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E3(:,:,i),E0(:,:,i),E3(:,:,i)); e116=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E3(:,:,i),E0(:,:,i),E3(:,:,i)); e180=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E3(:,:,i),E0(:,:,i),E3(:,:,i)); e244=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E3(:,:,i),E0(:,:,i),E3(:,:,i));
e53=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E3(:,:,i),E1(:,:,i),E0(:,:,i)); e117=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E3(:,:,i),E1(:,:,i),E0(:,:,i)); e181=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E3(:,:,i),E1(:,:,i),E0(:,:,i)); e245=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E3(:,:,i),E1(:,:,i),E0(:,:,i));
e54=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E3(:,:,i),E1(:,:,i),E1(:,:,i)); e118=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E3(:,:,i),E1(:,:,i),E1(:,:,i)); e182=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E3(:,:,i),E1(:,:,i),E1(:,:,i)); e246=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E3(:,:,i),E1(:,:,i),E1(:,:,i));
e55=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E3(:,:,i),E1(:,:,i),E2(:,:,i)); e119=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E3(:,:,i),E1(:,:,i),E2(:,:,i)); e183=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E3(:,:,i),E1(:,:,i),E2(:,:,i)); e247=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E3(:,:,i),E1(:,:,i),E2(:,:,i));
e56=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E3(:,:,i),E1(:,:,i),E3(:,:,i)); e120=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E3(:,:,i),E1(:,:,i),E3(:,:,i)); e184=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E3(:,:,i),E1(:,:,i),E3(:,:,i)); e248=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E3(:,:,i),E1(:,:,i),E3(:,:,i));
e57=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E3(:,:,i),E2(:,:,i),E0(:,:,i)); e121=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E3(:,:,i),E2(:,:,i),E0(:,:,i)); e185=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E3(:,:,i),E2(:,:,i),E0(:,:,i)); e249=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E3(:,:,i),E2(:,:,i),E0(:,:,i));
e58=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E3(:,:,i),E2(:,:,i),E1(:,:,i)); e122=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E3(:,:,i),E2(:,:,i),E1(:,:,i)); e186=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E3(:,:,i),E2(:,:,i),E1(:,:,i)); e250=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E3(:,:,i),E2(:,:,i),E1(:,:,i));
e59=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E3(:,:,i),E2(:,:,i),E2(:,:,i)); e123=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E3(:,:,i),E2(:,:,i),E2(:,:,i)); e187=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E3(:,:,i),E2(:,:,i),E2(:,:,i)); e251=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E3(:,:,i),E2(:,:,i),E2(:,:,i));
e60=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E3(:,:,i),E2(:,:,i),E3(:,:,i)); e124=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E3(:,:,i),E2(:,:,i),E3(:,:,i)); e188=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E3(:,:,i),E2(:,:,i),E3(:,:,i)); e252=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E3(:,:,i),E2(:,:,i),E3(:,:,i));
e61=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E3(:,:,i),E3(:,:,i),E0(:,:,i)); e125=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E3(:,:,i),E3(:,:,i),E0(:,:,i)); e189=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E3(:,:,i),E3(:,:,i),E0(:,:,i)); e253=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E3(:,:,i),E3(:,:,i),E0(:,:,i));
e62=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E3(:,:,i),E3(:,:,i),E1(:,:,i)); e126=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E3(:,:,i),E3(:,:,i),E1(:,:,i)); e190=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E3(:,:,i),E3(:,:,i),E1(:,:,i)); e254=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E3(:,:,i),E3(:,:,i),E1(:,:,i));
e63=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E3(:,:,i),E3(:,:,i),E2(:,:,i)); e127=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E3(:,:,i),E3(:,:,i),E2(:,:,i)); e191=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E3(:,:,i),E3(:,:,i),E2(:,:,i)); e255=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E3(:,:,i),E3(:,:,i),E2(:,:,i));
e64=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E3(:,:,i),E3(:,:,i),E3(:,:,i)); e128=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E3(:,:,i),E3(:,:,i),E3(:,:,i)); e192=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E3(:,:,i),E3(:,:,i),E3(:,:,i)); e256=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E3(:,:,i),E3(:,:,i),E3(:,:,i));

rho_source_depolarze=e1*rho_source*e1'+e2*rho_source*e2'+e3*rho_source*e3'+e4*rho_source*e4'+e5*rho_source*e5'+e6*rho_source*e6'+e7*rho_source*e7'+e8*rho_source*e8'+e9*rho_source*e9'+e10*rho_source*e10'+e11*rho_source*e11'+e12*rho_source*e12'+e13*rho_source*e13'+e14*rho_source*e14'+e15*rho_source*e15'+e16*rho_source*e16'+e17*rho_source*e17'+e18*rho_source*e18'+e19*rho_source*e19'+e20*rho_source*e20'+e21*rho_source*e21'+e22*rho_source*e22'+e23*rho_source*e23'+e24*rho_source*e24'+e25*rho_source*e25'+e26*rho_source*e26'+e27*rho_source*e27'+e28*rho_source*e28'+e29*rho_source*e29'+e30*rho_source*e30'+e31*rho_source*e31'+e32*rho_source*e32'+e33*rho_source*e33'+e34*rho_source*e34'+e35*rho_source*e35'+e36*rho_source*e36'+e37*rho_source*e37'+e38*rho_source*e38'+e39*rho_source*e39'+e40*rho_source*e40'+e41*rho_source*e41'+e42*rho_source*e42'+e43*rho_source*e43'+e44*rho_source*e44'+e45*rho_source*e45'+e46*rho_source*e46'+e47*rho_source*e47'+e48*rho_source*e48'+e49*rho_source*e49'+e50*rho_source*e50'+e51*rho_source*e51'+e52*rho_source*e52'+e53*rho_source*e53'+e54*rho_source*e54'+e55*rho_source*e55'+e56*rho_source*e56'+e57*rho_source*e57'+e58*rho_source*e58'+e59*rho_source*e59'+e60*rho_source*e60'+e61*rho_source*e61'+e62*rho_source*e62'+e63*rho_source*e63'+e64*rho_source*e64'+e65*rho_source*e65'+e66*rho_source*e66'+e67*rho_source*e67'+e68*rho_source*e68'+e69*rho_source*e69'+e70*rho_source*e70'+e71*rho_source*e71'+e72*rho_source*e72'+e73*rho_source*e73'+e74*rho_source*e74'+e75*rho_source*e75'+e76*rho_source*e76'+e77*rho_source*e77'+e78*rho_source*e78'+e79*rho_source*e79'+e80*rho_source*e80'+e81*rho_source*e81'+e82*rho_source*e82'+e83*rho_source*e83'+e84*rho_source*e84'+e85*rho_source*e85'+e86*rho_source*e86'+e87*rho_source*e87'+e88*rho_source*e88'+e89*rho_source*e89'+e90*rho_source*e90'+e91*rho_source*e91'+e92*rho_source*e92'+e93*rho_source*e93'+e94*rho_source*e94'+e95*rho_source*e95'+e96*rho_source*e96'+e97*rho_source*e97'+e98*rho_source*e98'+e99*rho_source*e99'+e100*rho_source*e100'+...
	e101*rho_source*e101'+e102*rho_source*e102'+e103*rho_source*e103'+e104*rho_source*e104'+e105*rho_source*e105'+e106*rho_source*e106'+e107*rho_source*e107'+e108*rho_source*e108'+e109*rho_source*e109'+e110*rho_source*e110'+e111*rho_source*e111'+e112*rho_source*e112'+e113*rho_source*e113'+e114*rho_source*e114'+e115*rho_source*e115'+e116*rho_source*e116'+e117*rho_source*e117'+e118*rho_source*e118'+e119*rho_source*e119'+e120*rho_source*e120'+e121*rho_source*e121'+e122*rho_source*e122'+e123*rho_source*e123'+e124*rho_source*e124'+e125*rho_source*e125'+e126*rho_source*e126'+e127*rho_source*e127'+e128*rho_source*e128'+e129*rho_source*e129'+e130*rho_source*e130'+e131*rho_source*e131'+e132*rho_source*e132'+e133*rho_source*e133'+e134*rho_source*e134'+e135*rho_source*e135'+e136*rho_source*e136'+e137*rho_source*e137'+e138*rho_source*e138'+e139*rho_source*e139'+e140*rho_source*e140'+e141*rho_source*e141'+e142*rho_source*e142'+e143*rho_source*e143'+e144*rho_source*e144'+e145*rho_source*e145'+e146*rho_source*e146'+e147*rho_source*e147'+e148*rho_source*e148'+e149*rho_source*e149'+e150*rho_source*e150'+e151*rho_source*e151'+e152*rho_source*e152'+e153*rho_source*e153'+e154*rho_source*e154'+e155*rho_source*e155'+e156*rho_source*e156'+e157*rho_source*e157'+e158*rho_source*e158'+e159*rho_source*e159'+e160*rho_source*e160'+e161*rho_source*e161'+e162*rho_source*e162'+e163*rho_source*e163'+e164*rho_source*e164'+e165*rho_source*e165'+e166*rho_source*e166'+e167*rho_source*e167'+e168*rho_source*e168'+e169*rho_source*e169'+e170*rho_source*e170'+e171*rho_source*e171'+e172*rho_source*e172'+e173*rho_source*e173'+e174*rho_source*e174'+e175*rho_source*e175'+e176*rho_source*e176'+e177*rho_source*e177'+e178*rho_source*e178'+e179*rho_source*e179'+e180*rho_source*e180'+e181*rho_source*e181'+e182*rho_source*e182'+e183*rho_source*e183'+e184*rho_source*e184'+e185*rho_source*e185'+e186*rho_source*e186'+e187*rho_source*e187'+e188*rho_source*e188'+e189*rho_source*e189'+e190*rho_source*e190'+e191*rho_source*e191'+e192*rho_source*e192'+e193*rho_source*e193'+e194*rho_source*e194'+e195*rho_source*e195'+e196*rho_source*e196'+e197*rho_source*e197'+e198*rho_source*e198'+e199*rho_source*e199'+e200*rho_source*e200'+...
	e201*rho_source*e201'+e202*rho_source*e202'+e203*rho_source*e203'+e204*rho_source*e204'+e205*rho_source*e205'+e206*rho_source*e206'+e207*rho_source*e207'+e208*rho_source*e208'+e209*rho_source*e209'+e210*rho_source*e210'+e211*rho_source*e211'+e212*rho_source*e212'+e213*rho_source*e213'+e214*rho_source*e214'+e215*rho_source*e215'+e216*rho_source*e216'+e217*rho_source*e217'+e218*rho_source*e218'+e219*rho_source*e219'+e220*rho_source*e220'+e221*rho_source*e221'+e222*rho_source*e222'+e223*rho_source*e223'+e224*rho_source*e224'+e225*rho_source*e225'+e226*rho_source*e226'+e227*rho_source*e227'+e228*rho_source*e228'+e229*rho_source*e229'+e230*rho_source*e230'+e231*rho_source*e231'+e232*rho_source*e232'+e233*rho_source*e233'+e234*rho_source*e234'+e235*rho_source*e235'+e236*rho_source*e236'+e237*rho_source*e237'+e238*rho_source*e238'+e239*rho_source*e239'+e240*rho_source*e240'+e241*rho_source*e241'+e242*rho_source*e242'+e243*rho_source*e243'+e244*rho_source*e244'+e245*rho_source*e245'+e246*rho_source*e246'+e247*rho_source*e247'+e248*rho_source*e248'+e249*rho_source*e249'+e250*rho_source*e250'+e251*rho_source*e251'+e252*rho_source*e252'+e253*rho_source*e253'+e254*rho_source*e254'+e255*rho_source*e255'+e256*rho_source*e256';

rho_all_depolarize=U*rho_source_depolarze*U';
t=(1:4);
rho_out_depolarize=rho_all_depolarize(t,t);
while(t(end)<=2^7)
    t=t+4;
    if t(end)<=2^7
    rho_out_depolarize=rho_out_depolarize+rho_all_depolarize(t,t);
    end
end
f_out_depolarize(i)=phi_in'*rho_out_depolarize*phi_in;
end

disp(toc);

%% draw the figure, and show the legend of all these 3 noise with its minimum.
plot(p,f_out_x)
hold on;
[f_out_x_min,f_out_x_I]=min(f_out_x); 
plot(p(f_out_x_I),f_out_x(f_out_x_I),'.','MarkerSize',10)     
str = ['(' num2str(p(f_out_x_I)) ',' num2str(double(f_out_x(f_out_x_I))) ')'];
text(p(f_out_x_I),f_out_x(f_out_x_I),str)                       
hold on;

plot(p,f_out_y)
hold on;
[f_out_y_min,f_out_y_I]=min(f_out_y);                
plot(p(f_out_y_I),f_out_y(f_out_y_I),'.','MarkerSize',10)     
str = ['(' num2str(p(f_out_y_I)) ',' num2str(double(f_out_y(f_out_y_I))) ')'];
text(p(f_out_y_I),f_out_y(f_out_y_I),str)                       
hold on;

plot(p,f_out_z)
hold on;
[f_out_z_min,f_out_z_I]=min(f_out_z); 
plot(p(f_out_z_I),f_out_z(f_out_z_I),'.','MarkerSize',10)     
str = ['(' num2str(p(f_out_z_I)) ',' num2str(double(f_out_z(f_out_z_I))) ')'];
text(p(f_out_z_I),f_out_z(f_out_z_I),str)                       
hold on;

plot(p,f_out_ampdamp)
hold on;
[f_out_ampdamp_min,f_out_ampdamp_I]=min(f_out_ampdamp); 
plot(p(f_out_ampdamp_I),f_out_ampdamp(f_out_ampdamp_I),'.','MarkerSize',10)     
str = ['(' num2str(p(f_out_ampdamp_I)) ',' num2str(double(f_out_ampdamp(f_out_ampdamp_I))) ')'];
text(p(f_out_ampdamp_I),f_out_ampdamp(f_out_ampdamp_I),str)                       
hold on;

plot(p,f_out_depolarize)
hold on;
[f_out_depolarize_min,f_out_depolarize_I]=min(f_out_depolarize); 
plot(p(f_out_depolarize_I),f_out_depolarize(f_out_depolarize_I),'.','MarkerSize',10)     
str = ['(' num2str(p(f_out_depolarize_I)) ',' num2str(double(f_out_depolarize(f_out_depolarize_I))) ')'];
text(p(f_out_depolarize_I),f_out_depolarize(f_out_depolarize_I),str)                       
hold on;

ff=[f_out_x;f_out_y;f_out_z;f_out_ampdamp;f_out_depolarize]';

legend('x noise','x noise min','y noise','y noise min','z noise','z noise min','amplitude damping noise','amplitude damping noise min','depolarizing noise','depolarizing noise min');