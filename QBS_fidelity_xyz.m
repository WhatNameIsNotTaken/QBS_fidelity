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

legend('x noise','x noise min','y noise','y noise min','z noise','z noise min','amplitude damping noise','amplitude damping noise min');