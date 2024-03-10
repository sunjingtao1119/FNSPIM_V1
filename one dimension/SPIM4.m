%%
% one-dimension  nearly symplectic precise integration method 
%%
clear;close all;
% Initialization of setup
nt    = 7500 ;     % number of time steps
xmax  = 8000;      % Length of domain [m]
N     = 4;         % Order of Lagrange polynomials 
ne    = 150 ;      % Number of elements
Tdom  =0.4 ;       % Dominant period of Ricker source wavelet
vs    = 1500;      % S velocity [m/s]
rho   = 1 ;         % Density [kg/m^
% variables for elemental matrices
Me =zeros(N+1,1);
Ke = zeros(N+1,N+1);
%Initialization of GLL points integration weights
[xi,w]= gll(N) ;  % xi, N+1 coordinates [-1 1] of GLL points
                           %w Integration weights at GLL locations
%Space domain
le = xmax/ne;        % Length of elements
%Vector with GLL points  坐标变换
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
alpha = 10;          % Courant value
dt = alpha*dxmin/vs ; % Global time step
% u0=exp(-(xg-4000).^2/8000^2)-exp(-(xg+4000).^2/8000^2);
% u0=exp(-(xg-2000).^2/20^2);
% Mapping - Jacobian
J = le/2 ;
Ji = 1/J ;   % Inverse Jacobian

%1st derivative of Lagrange polynomials
l1d = lagrange1st(N) ;  % Array with GLL as columns for each N+1 polynomial

%%
% el_span=25;
% percent=0.4;
% a=el_span*N+1;
% b=floor((N*ne+1)/2);
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
% treatment of boundary condition , fixed boundary
K(1,1)=1;K(end,end)=1;
K(1,2:end)=0;
K(end,1:end-1)=0;
Minv(1,1)=1; Minv(end,end)=1;
% # initialize source time function and force vector f
src  = ricker(dt,Tdom);
% isrc = int(np.floor(ng/2))   # Source location
isrc = floor(ng/2);
% # Initialization of solution vectors
u = zeros(ng,1);
uold = u;
unew = u;
f = u ;
LI=Minv*K;
%
C=zeros(ng,ng);
A=[ zeros(ng,ng),eye(ng,ng);-Minv*K,zeros(ng,ng)];
I=eye(size(A));
% dt=0.1;
b1=norm(A,1);
b2=norm(A,2);
b=max(b1,b2);
mb=b*dt;
% N=ceil(log2(mb));
%%
u=log2(mb);
L=5;
w=log2(10)*L;
N=ceil(u+(u+w-log2(720))/4);
% N=45;
m=1;
tau=dt/m;
Ta=(A*tau/2)+(A*tau)^2/12;
%利用幂级数近似
B=(A*tau/2)-(A*tau)^2/12;
T=(I+Ta)*inv(I-B);
% T(find(abs(T)<1e-15))=0;
% T=sparse(T);
spy(T,'red');
% 
c=[7/90;16/45;2/15;16/45;7/90];
v=zeros(2*ng,1);u0=exp(-(xg-4000).^2/200^2);
up=-(xg-4000)./100^2.*u0;
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
% T=sparse(T);
tic
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
toc











