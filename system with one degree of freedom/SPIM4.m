% Forth-order nearly symplectic precise integration method 
clear;close all;
clc;
w=50;               % frequency 
A=[0,1;-w^2,0];     % system matrix 
I=eye(size(A));
dt=100;            %  time step
%% Parameter N estimation
b=norm(A,1);
mb=b*dt;
u=log2(mb);
L=5;
w1=log2(10)*L;
n=ceil(u+(u+w1-log2(12))/4);
%% symplectic precise integration method 
m=2^n;
tau=dt/m;
Ta1=(A*tau/2)+(A*tau)^2/12;
Ta2=-(A*tau/2)+(A*tau)^2/12;
for i=1:n
    Ta1=2*Ta1+Ta1*Ta1;
    Ta2=2*Ta2+Ta2*Ta2;
end
T1=I+Ta1;
T2=I+Ta2;
T=T1*inv(T2);  % transfer matrix
v=[1;0];       % Initialize displacement and velocity
N=1000000;      % interations

t=dt.*[0:N-1];
V=zeros(2,N);  
V(:,1)=v;
H=zeros(N,1);
E0=w^2;
H(1)=w^2*v(1)^2+v(2)^2;  % Energy
% Time extrapolation
for i=2:N
    v=T*v;
    V(:,i)=v;
    H(i)=w^2*v(1)^2+v(2)^2;
end
E=abs((H-E0))./E0;
maxE=max(E)
pe=cos(w.*t); % exact solution
subplot(2,1,1)
plot(t,E,'k')
grid on
title(' SPIM4')
xlabel('t(s)')
ylabel('\it{E_{H}}')
set(gca,'FontName','Helvetica ','FontSize',12);
Ru=norm(pe-V(1,:),2)/norm(pe,2)
subplot(2,1,2)
E=pe-V(1,:);
plot(t,E,'k')
grid on
xlabel('t(s)')
ylabel('\it{E_{u}}')
grid on
set(gca,'FontSize',12);
set(gcf,'Position',[100 100 500 350])
