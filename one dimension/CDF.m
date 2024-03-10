%  central finite difference (CDF) method. 
clear;close all;
% Initialization of setup
nt    = 200000;         % number of time steps
xmax  = 8000;      % Length of domain [m]
N     = 4;              % Order of Lagrange polynomials 
ne    = 150 ;         % Number of elements
Tdom  =0.4;         % Dominant period of Ricker source wavelet
iplot = 30  ;          % Plotting each iplot snapshot
vs    = 4000;          % S velocity [m/s]
rho   = 2000  ;       % Density [kg/m^
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
ep = 0.8;                   % Courant value
dt = ep*dxmin/vs ;          % Time step
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
M = zeros(2*ng,1) ;
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
 Minv=sparse(Minv);

% # Display inverse mass matrix inv(M)
% # --------------------------------------------------------------
figure(2)
image(Minv)
title('Mass Matrix M')
axis("off")
axis("equal")
%
% # Global Stiffness Matrix
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

t = (1:nt)*dt;

u=exp(-(xg-4000).^2/200^2);
v=(-(xg-4000)/100^2.*u);
uold=u-v*dt;

% # ---------------------------------------------------------------
% # Time extrapolation
% # ---------------------------------------------------------------
xt= [ ];
%固定边界条件
K(1,1)=1;K(end,end)=1;
K(1,2:end)=0;
K(end,1:end-1)=0;
K=sparse(K);
% 边界条件处理
L=Minv*K;
for it =1:nt
    unew = dt^2. * Minv * (- K *u) + 2 * u - uold;
    uold=u;
    u = unew;
plot(xg,u,'linewidth',1.5);
set(gcf,'color','w'); 
ylim([-1, 1])
snt=string(num2str(nt));
str4='Time Step nt ='+string(it);
 title( str4) 
drawnow();    

end



