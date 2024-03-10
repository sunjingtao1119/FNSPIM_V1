%Forth-order explicit symplectic Runge-Kutta method 
clear;close all;
w=50;            % frequency
dt=0.05;        % time steps
N=1000000;        %iterations
%initialize
p1=zeros(1,N);    % displacement
q1=zeros(1,N);    % velocity
p1(1)=0;
q1(1)=1;
% Runge-Kutta coefficient
e1=0.13883725894365473;
e2=0.46958619250378464;
e3=0.751399209882663;
e4=-0.3598226613301023;
d1=0.3726518368174738;
d2=0.41264784985125225;
d3=-0.04864313400799411;
d4=0.26334344733926796;
H=zeros(N,1);
E0=w^2*1^2;
H(1)=w^2*q1(1)^2+p1(1)^2; % Energy
for i=2:N
    qq1=q1(i-1)+(e1*dt)*(p1(i-1));
    pp1=p1(i-1)+(d1*dt)*(-w^2*qq1);
    qq1=qq1+(e2*dt)*(pp1);
    pp1=pp1+(d2*dt)*(-w^2*qq1);
    qq1=qq1+(e3*dt)*(pp1);
    pp1=pp1+(d3*dt)*(-w^2*qq1);
    q1(i)=qq1+dt*(e4)*(pp1);
    p1(i)=pp1+d4*dt*(-w^2*q1(i));
    H(i)=w^2*q1(i)^2+p1(i)^2; 
end
t=dt.*[0:N-1];
pe=cos(w.*t);   % exact solution
E=abs((H-E0))./E0;
max(E)
subplot(2,1,1)
plot(t,E)
ylim([-1e-7,5e-7])
grid on
title('SPRK4')
xlabel('Time/(s)')
ylabel('\it E_{H}','FontName','Helvetica')
set(gca,'FontName','Helvetica ','FontSize',12);
Ru=norm((q1-pe),2)/norm(pe,2)
subplot(2,1,2)
plot(t,q1-pe)
grid on
xlabel('Time/(s)')
ylabel('\it R_{u}','FontName','Helvetica')
grid on
set(gca,'FontName','Helvetica ','FontSize',12);
set(gcf,'Position',[100 100 500 350])
