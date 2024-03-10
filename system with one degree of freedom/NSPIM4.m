% Forth-order nearly symplectic precise integration method 
clear;close all;
clc;
w=50;
A=[0,1;-w^2,0];
I=eye(size(A));
dt=100;  % time step
b=norm(A,1);
mb=b*dt;
% if n1<0 | n1==0
% n1=2;
% end
u=log2(mb);
L=5;
w1=log2(10)*L;
n=ceil(u+(u+w1-log2(12))/4);

m=2^n;
tau=dt/m;
Ta=(A*tau/2)+(A*tau)^2/12;
% nearly Páde approximation 
B=((A*tau/2)-(A*tau)^2/12);
Ta1=B+B^2+B^3+B^4+B^5+B^6;
Ta=Ta*Ta1+Ta+Ta1;
% precise integration
for i=1:n
    Ta=2*Ta+Ta*Ta;
end
T=I+Ta;    % transfer matrix
v=[1;0];   % Initialize displacement and velocity
N=1000000;  % interations
t=dt.*[0:N-1];
V=zeros(2,N);
V(:,1)=v;
H=zeros(N,1);
E0=w^2;
H(1)=w^2*v(1)^2+v(2)^2;  % energy
for i=2:N
    v=T*v;
    V(:,i)=v;
    H(i)=w^2*v(1)^2+v(2)^2;
end
E=abs((H-E0))./E0;  
max(E)
pe=cos(w.*t);  % exact solution
Ru=norm(pe-V(1,:),2)/norm(pe,2)
subplot(2,1,1)
plot(t,E,'k')
% ylim([-0.001,0.003])
grid on
title(' NSPIM4')
xlabel('t(s)')
ylabel('\it E_{H}')
set(gca,'FontSize',12);
subplot(2,1,2)
E=pe-V(1,:);
plot(t,E,'k')
%  ylim([-4,4])
grid on
xlabel('t(s)')
ylabel('\it E_{u}')
set(gca,'FontSize',12);
set(gcf,'Position',[100 100 500 350])
