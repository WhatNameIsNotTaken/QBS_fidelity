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
    e1 =mykron(eye(2),eye(2),E0(:,:,i),eye(2),E0(:,:,i),E0(:,:,i));e2 =mykron(eye(2),eye(2),E0(:,:,i),eye(2),E0(:,:,i),E1(:,:,i));
    e3 =mykron(eye(2),eye(2),E0(:,:,i),eye(2),E1(:,:,i),E0(:,:,i));e4 =mykron(eye(2),eye(2),E0(:,:,i),eye(2),E1(:,:,i),E1(:,:,i));
    e5 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E0(:,:,i),E0(:,:,i));e6 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E0(:,:,i),E1(:,:,i));
    e7 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E1(:,:,i),E0(:,:,i));e8 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E1(:,:,i),E1(:,:,i));

% output with x noise
rho_source_x=e1*rho_source*e1'+e2*rho_source*e2'+e3*rho_source*e3'+e4*rho_source*e4'+e5*rho_source*e5'+e6*rho_source*e6'+...
    e7*rho_source*e7'+e8*rho_source*e8';
rho_all_x=U*rho_source_x*U';
t=(1:4);
rho_out_x=rho_all_x(t,t);
while(t(end)<=2^6)
    t=t+4;
    if t(end)<=2^6
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
    e1 =mykron(eye(2),eye(2),E0(:,:,i),eye(2),E0(:,:,i),E0(:,:,i));e2 =mykron(eye(2),eye(2),E0(:,:,i),eye(2),E0(:,:,i),E1(:,:,i));
    e3 =mykron(eye(2),eye(2),E0(:,:,i),eye(2),E1(:,:,i),E0(:,:,i));e4 =mykron(eye(2),eye(2),E0(:,:,i),eye(2),E1(:,:,i),E1(:,:,i));
    e5 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E0(:,:,i),E0(:,:,i));e6 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E0(:,:,i),E1(:,:,i));
    e7 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E1(:,:,i),E0(:,:,i));e8 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E1(:,:,i),E1(:,:,i));

% output with y noise
rho_source_y=e1*rho_source*e1'+e2*rho_source*e2'+e3*rho_source*e3'+e4*rho_source*e4'+e5*rho_source*e5'+e6*rho_source*e6'+...
    e7*rho_source*e7'+e8*rho_source*e8';
rho_all_y=U*rho_source_y*U';
t=(1:4);
rho_out_y=rho_all_y(t,t);
while(t(end)<=2^6)
    t=t+4;
    if t(end)<=2^6
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
    e1 =mykron(eye(2),eye(2),E0(:,:,i),eye(2),E0(:,:,i),E0(:,:,i));e2 =mykron(eye(2),eye(2),E0(:,:,i),eye(2),E0(:,:,i),E1(:,:,i));
    e3 =mykron(eye(2),eye(2),E0(:,:,i),eye(2),E1(:,:,i),E0(:,:,i));e4 =mykron(eye(2),eye(2),E0(:,:,i),eye(2),E1(:,:,i),E1(:,:,i));
    e5 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E0(:,:,i),E0(:,:,i));e6 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E0(:,:,i),E1(:,:,i));
    e7 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E1(:,:,i),E0(:,:,i));e8 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E1(:,:,i),E1(:,:,i));

% output with z noise
rho_source_z=e1*rho_source*e1'+e2*rho_source*e2'+e3*rho_source*e3'+e4*rho_source*e4'+e5*rho_source*e5'+e6*rho_source*e6'+...
    e7*rho_source*e7'+e8*rho_source*e8';
rho_all_z=U*rho_source_z*U';
t=(1:4);
rho_out_z=rho_all_z(t,t);
while(t(end)<=2^6)
    t=t+4;
    if t(end)<=2^6
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
    e1 =mykron(eye(2),eye(2),E0(:,:,i),eye(2),E0(:,:,i),E0(:,:,i));e2 =mykron(eye(2),eye(2),E0(:,:,i),eye(2),E0(:,:,i),E1(:,:,i));
    e3 =mykron(eye(2),eye(2),E0(:,:,i),eye(2),E1(:,:,i),E0(:,:,i));e4 =mykron(eye(2),eye(2),E0(:,:,i),eye(2),E1(:,:,i),E1(:,:,i));
    e5 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E0(:,:,i),E0(:,:,i));e6 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E0(:,:,i),E1(:,:,i));
    e7 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E1(:,:,i),E0(:,:,i));e8 =mykron(eye(2),eye(2),E1(:,:,i),eye(2),E1(:,:,i),E1(:,:,i));

% output with amplitude damping noise
rho_source_ampdamp=e1*rho_source*e1'+e2*rho_source*e2'+e3*rho_source*e3'+e4*rho_source*e4'+e5*rho_source*e5'+e6*rho_source*e6'+...
    e7*rho_source*e7'+e8*rho_source*e8';
rho_all_ampdamp=U*rho_source_ampdamp*U';
t=(1:4);
rho_out_ampdamp=rho_all_ampdamp(t,t);
while(t(end)<=2^6)
    t=t+4;
    if t(end)<=2^6
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

e1 =mykron(eye(2),eye(2),E0(:,:,i),eye(2),E0(:,:,i),E0(:,:,i)); e33=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E0(:,:,i),E0(:,:,i));
e2 =mykron(eye(2),eye(2),E0(:,:,i),eye(2),E0(:,:,i),E1(:,:,i)); e34=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E0(:,:,i),E1(:,:,i));
e3 =mykron(eye(2),eye(2),E0(:,:,i),eye(2),E0(:,:,i),E2(:,:,i)); e35=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E0(:,:,i),E2(:,:,i));
e4 =mykron(eye(2),eye(2),E0(:,:,i),eye(2),E0(:,:,i),E3(:,:,i)); e36=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E0(:,:,i),E3(:,:,i));
e5 =mykron(eye(2),eye(2),E0(:,:,i),eye(2),E1(:,:,i),E0(:,:,i)); e37=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E1(:,:,i),E0(:,:,i));
e6 =mykron(eye(2),eye(2),E0(:,:,i),eye(2),E1(:,:,i),E1(:,:,i)); e38=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E1(:,:,i),E1(:,:,i));
e7 =mykron(eye(2),eye(2),E0(:,:,i),eye(2),E1(:,:,i),E2(:,:,i)); e39=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E1(:,:,i),E2(:,:,i));
e8 =mykron(eye(2),eye(2),E0(:,:,i),eye(2),E1(:,:,i),E3(:,:,i)); e40=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E1(:,:,i),E3(:,:,i));
e9 =mykron(eye(2),eye(2),E0(:,:,i),eye(2),E2(:,:,i),E0(:,:,i)); e41=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E2(:,:,i),E0(:,:,i));
e10=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E2(:,:,i),E1(:,:,i)); e42=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E2(:,:,i),E1(:,:,i));
e11=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E2(:,:,i),E2(:,:,i)); e43=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E2(:,:,i),E2(:,:,i));
e12=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E2(:,:,i),E3(:,:,i)); e44=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E2(:,:,i),E3(:,:,i));
e13=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E3(:,:,i),E0(:,:,i)); e45=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E3(:,:,i),E0(:,:,i));
e14=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E3(:,:,i),E1(:,:,i)); e46=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E3(:,:,i),E1(:,:,i));
e15=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E3(:,:,i),E2(:,:,i)); e47=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E3(:,:,i),E2(:,:,i));
e16=mykron(eye(2),eye(2),E0(:,:,i),eye(2),E3(:,:,i),E3(:,:,i)); e48=mykron(eye(2),eye(2),E2(:,:,i),eye(2),E3(:,:,i),E3(:,:,i));
e17=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E0(:,:,i),E0(:,:,i)); e49=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E0(:,:,i),E0(:,:,i));
e18=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E0(:,:,i),E1(:,:,i)); e50=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E0(:,:,i),E1(:,:,i));
e19=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E0(:,:,i),E2(:,:,i)); e51=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E0(:,:,i),E2(:,:,i));
e20=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E0(:,:,i),E3(:,:,i)); e52=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E0(:,:,i),E3(:,:,i));
e21=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E1(:,:,i),E0(:,:,i)); e53=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E1(:,:,i),E0(:,:,i));
e22=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E1(:,:,i),E1(:,:,i)); e54=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E1(:,:,i),E1(:,:,i));
e23=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E1(:,:,i),E2(:,:,i)); e55=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E1(:,:,i),E2(:,:,i));
e24=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E1(:,:,i),E3(:,:,i)); e56=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E1(:,:,i),E3(:,:,i));
e25=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E2(:,:,i),E0(:,:,i)); e57=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E2(:,:,i),E0(:,:,i));
e26=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E2(:,:,i),E1(:,:,i)); e58=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E2(:,:,i),E1(:,:,i));
e27=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E2(:,:,i),E2(:,:,i)); e59=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E2(:,:,i),E2(:,:,i));
e28=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E2(:,:,i),E3(:,:,i)); e60=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E2(:,:,i),E3(:,:,i));
e29=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E3(:,:,i),E0(:,:,i)); e61=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E3(:,:,i),E0(:,:,i));
e30=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E3(:,:,i),E1(:,:,i)); e62=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E3(:,:,i),E1(:,:,i));
e31=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E3(:,:,i),E2(:,:,i)); e64=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E3(:,:,i),E3(:,:,i));
e32=mykron(eye(2),eye(2),E1(:,:,i),eye(2),E3(:,:,i),E3(:,:,i)); e63=mykron(eye(2),eye(2),E3(:,:,i),eye(2),E3(:,:,i),E2(:,:,i));

rho_source_depolarze=e1*rho_source*e1'+e2*rho_source*e2'+e3*rho_source*e3'+e4*rho_source*e4'+e5*rho_source*e5'+e6*rho_source*e6'+e7*rho_source*e7'+e8*rho_source*e8'+e9*rho_source*e9'+e10*rho_source*e10'+e11*rho_source*e11'+e12*rho_source*e12'+e13*rho_source*e13'+e14*rho_source*e14'+e15*rho_source*e15'+e16*rho_source*e16'+e17*rho_source*e17'+e18*rho_source*e18'+e19*rho_source*e19'+e20*rho_source*e20'+e21*rho_source*e21'+e22*rho_source*e22'+e23*rho_source*e23'+e24*rho_source*e24'+e25*rho_source*e25'+e26*rho_source*e26'+e27*rho_source*e27'+e28*rho_source*e28'+e29*rho_source*e29'+e30*rho_source*e30'+e31*rho_source*e31'+e32*rho_source*e32'+e33*rho_source*e33'+e34*rho_source*e34'+e35*rho_source*e35'+e36*rho_source*e36'+e37*rho_source*e37'+e38*rho_source*e38'+e39*rho_source*e39'+e40*rho_source*e40'+e41*rho_source*e41'+e42*rho_source*e42'+e43*rho_source*e43'+e44*rho_source*e44'+e45*rho_source*e45'+e46*rho_source*e46'+e47*rho_source*e47'+e48*rho_source*e48'+e49*rho_source*e49'+e50*rho_source*e50'+e51*rho_source*e51'+e52*rho_source*e52'+e53*rho_source*e53'+e54*rho_source*e54'+e55*rho_source*e55'+e56*rho_source*e56'+e57*rho_source*e57'+e58*rho_source*e58'+e59*rho_source*e59'+e60*rho_source*e60'+e61*rho_source*e61'+e62*rho_source*e62'+e63*rho_source*e63'+e64*rho_source*e64';

rho_all_depolarize=U*rho_source_depolarze*U';
t=(1:4);
rho_out_depolarize=rho_all_depolarize(t,t);
while(t(end)<=2^6)
    t=t+4;
    if t(end)<=2^6
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

legend('x noise','x noise min','y noise','y noise min','z noise','z noise min','amplitude damping noise','amplitude damping noise min','depolarizing noise','depolarizing noise min'...
   , 'x noise 1','x noise min 1','y noise 1','y noise min 1','z noise 1','z noise min 1','amplitude damping noise 1','amplitude damping noise min 1','depolarizing noise 1','depolarizing noise min 1');