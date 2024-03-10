%%
% one-dimension fast nearly symplectic precise integration method 
%%
clear;close all;
% Initialization of setup
nt    = 7500 ;         % number of time steps
xmax  = 8000;      % Length of domain [m]
N     = 4;              % Order of Lagrange polynomials 
ne    = 150 ;         % Number of elements
Tdom  =0.4 ;         % Dominant period of Ricker source wavelet
iplot = 30  ;          % Plotting each iplot snapshot
vs    = 1500;          % S velocity [m/s]
rho   = 1 ;       % Density [kg/m^
% variables for elemental matrices
Me =zeros(N+1,1);
Ke = zeros(N+1,N+1);
%Initialization of GLL points integration weights
[xi,w]= gll(N) ;  % xi, N+1 coordinates [-1 1] of GLL points
                           %w Integration weights at GLL locations
%Space domain
le = xmax/ne;        % Length of elements
%Vector with GLL points  
k = 1;
xg = zeros((N*ne)+1,1) ;
xg(1) = 0;
for i =1:ne
    for j =1:N
        k = k+1;
        xg(k) = (i-1)*le + .5*(xi(j+1)+1)*le;
    end
end

%%
dxmin = min(diff(xg)) ; 
ep = 2;              % Courant value
dt = ep*dxmin/vs ;     % time step
% Mapping - Jacobian
J = le/2 ;
Ji = 1/J ;   % Inverse Jacobian

%1st derivative of Lagrange polynomials
l1d = lagrange1st(N) ;  % Array with GLL as columns for each N+1 polynomial

%%
vs  = vs * ones(N*ne +1,1);
rho = rho * ones(N*ne +1,1);

% vs(b-a/2:b+a/2) = max(vs) * percent;
mu  = rho .* vs.^2  ;                % Shear modulus mu
figure(1)
plot(xg, vs, 'linewidth',1.5);
title('Velocity model','fontsize',16);
ylabel('Velocity (m/s)','fontsize',16);
xlabel('space (m)','fontsize',16);
xlim([0,xmax])
ylim([0,max(vs)+50])
axis("equal")
%%
k = 0;
m = 0;
ng = (ne-1)*N + N + 1;
M = zeros(ng,1) ;
for i =1:ne
%     # ------------------------------------
%     # Elemental Mass matrix
%     # ------------------------------------
    for l =0:N
        m=m+ 1;
        l=l+1;
        Me(l)= rho(m) * w(l) * J ;   %stored as a vector since it's diagonal
    end
    m =0;
    
    for j =0:N
        m=j+1;
        k = k + 1;
        if i>1
            if j==0
                k = k - 1;             
            end
        end
        M(k) = M(k) + Me(m);
    end
end
  % Computational Seismology A Practical Introduction p196
% # Inverse matrix of M 
% # --------------------------------------------------------------- 
Minv = eye(ng);
for i =0:ng-1
    k=i+1;
    Minv(k,k)= 1./M(k);
end
% # --------------------------------------------------------------- 
K = zeros(ng,ng);
xe = 0 ;

for e = 1:ne
    i0 = (e - 1)*N + 1;
    j0 = i0;
%     # ------------------------------------
%     # Elemental Stiffness Matrix
%     # ------------------------------------
    for i =0:N 
        for j = 0:N
            sum = 0;
            for k =0:N              
                sum = sum + mu(k+1+xe) * w(k+1) * Ji^2 * J * l1d(i+1,k+1) * l1d(j+1,k+1);
                Ke(i+1, j+1)= sum ;
            end
        end
    end
    xe=xe+N;
   for i =0:N
        for j = 0:N
            K(i0+i, j0+j) =  K(i0+i, j0+j) +Ke(i+1, j+1);
        end
    end
end
K(1,1)=1;K(end,end)=1;
K(1,2:end)=0;
K(end,1:end-1)=0;
Minv(1,1)=1; Minv(end,end)=1;
u = zeros(ng,1);
uold = u;
unew = u;
f = u ;
C=zeros(ng,ng);
% A=[ zeros(ng,ng),eye(ng,ng);-Minv*K,zeros(ng,ng)];
A=[ C,eye(ng,ng);-Minv*K,C];
LI=Minv*K;
I=eye(size(A));
% dt=0.1;
b=norm(A,2);
mb=b*dt;
% N=ceil(log2(mb));
%%
u=log2(mb);
L=5;
w=log2(10)*L;
N=ceil(u+(u+w-log2(720))/4);
m=2^N;
tau=dt/m;
Ta=(A*tau/2)+(A*tau)^2/12;
B=((A*tau/2)-(A*tau)^2/12);
% norm(B,inf)
Ta1=B+B^2+B^3+B^4+B^5+B^6;
Ta=Ta*Ta1+Ta+Ta1;
% cond(Ta) 
clear Ta1 B;
for i=1:N
    Ta=2*Ta+Ta*Ta;
    [m,n,src]=find(Ta);
    k=abs(src)<1e-25;   % Sparsification 
    src(k)=0;
    Ta=sparse(m,n,src,2*ng,2*ng);
end
clear m n k src;
T=sparse(I+Ta);  % Sparse storage
clear Ta;
%  Newton-Cotes integral coefficients
c=[7/90;16/45;2/15;16/45;7/90];
v=zeros(2*ng,1);
u0=exp(-(xg-4000).^2/200^2);  % initial displacement
up=-(xg-4000)/100^2.*u0;      % initial velocity
v(1:ng)=u0;
v(ng+1:end)=up;
r=v;
r2=zeros(size(r));
r1=r2;
r3=r2;
r4=r3;
r0=r1;
% V=zeros(4,101);
% V(:,1)=v;
t0=0;
rr=zeros(nt,1);

H=zeros(1,200000);
H(1)=1/2*(v(ng+1:end)'*v(ng+1:end)+v(1:ng)'*LI*v(1:ng));

for it=1:200000
v=T*v;
v(ng+1)=0;
v(end)=0;
v(1)=0;
v(ng)=0;
H(it)=1/2*(v(ng+1:end)'*v(ng+1:end)+v(1:ng)'*LI*v(1:ng));
plot(xg,v(1:ng),'linewidth',1.5);
set(gcf,'color','w');            
ylim([-1, 1]);
drawnow()
end
% max(abs((H-H(1)))./H(1))
% plot(xg,v(1:ng),'linewidth',1.5);
% set(gcf,'color','w');            
% ylim([-1, 1]);







