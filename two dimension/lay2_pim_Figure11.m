clear;clc
%% Generating mesh
seiscolor=load('colorbar.txt');
lex=6000;                  % Length in x direction
ley=6000;                  % Length in y direction
nex=200;                   % The number of units in x direction
ney=200;                   % The number of units in x direction
x=linspace(0,lex,nex);
y=linspace(0,ley,ney);
[x,y]=meshgrid(x,y);
dex=lex/nex;
dey=ley/ney;
X=x(:);
Y=y(:);
Point=[X,Y];
DT = delaunayTriangulation(Point);
Top=DT.ConnectivityList;
p=DT.Points;
ppx=p(:,1);
V=3000.*ones(size(ppx));
index=find(ppx<3150);
V(index)=1500;
F = freeBoundary(DT);
b=unique(F(:)); % boundary node
%% parameter setting
nt = 1500;  % total number of iteration times
v=2000;
ept=dex/v;
dt=0.4*ept;     % time step
T = (1:4*nt)*dt/4;
number_of_notes = length(ppx); % the total number of node
fmain =40;
%% Source
s_t = (1-2*pi^2*fmain*(T-0.2).^2).*exp(-fmain*pi^2*(T-0.2).^2);%ricker wavelet
%% Assembly stiffness matrix
dictK= TStiffness(p,Top,V);    % Stiffness matrix -double key value dictionary

K=dic2sp2(dictK,number_of_notes,number_of_notes);
%% boundary condition treatment
K(b,:)=0;
K(:,b)=0;
K(b,b)=1;
%%
clear dictK;
dictM= TMass(p,Top);              % Mass matrix -double key value dictionary
Minv=dic2sp2(dictM,number_of_notes,number_of_notes,1);
Minv(b,:)=0;
Minv(:,b)=0;
% Minv(b,b)=1;
clear dictM;
rec=10; y=linspace(1000,5000,100);
s_x=6000;s_y=6000;
[source_x,m] = find(abs(p(:,1))>=s_x/2-dex&abs(p(:,1))<=s_x/2+dex&abs(p(:,2))>=s_y/2-dey&abs(p(:,2))<=s_y/2+dey);% Source location
isrc=source_x(1);
rx=2500;
for i=1:100
   [rec,m] = find(abs(p(:,1))>=rx-dex&abs(p(:,1))<=rx+dex&abs(p(:,2))>=y(i)-dey&abs(p(:,2))<=y(i)+dey);
    indrec(i,:)=rec; % reciver location
end


%% Transfer matrix assembly
I1=speye(number_of_notes,number_of_notes);
L=Minv*K;
mm=Minv(isrc,isrc);
[m1,n1,src1]=find(I1);
n1=n1+number_of_notes;
[m2,n2,src2]=find(L);
m2=m2+number_of_notes;
src2=-1.*src2;
m=[m1;m2];
n=[n1;n2];
src=[src1;src2];
A=sparse(m,n,src,2*number_of_notes,2*number_of_notes);  % Generate sparse matrix
I=speye(size(A));
N=10;
m1=2^N;
tau=dt/m1;
Ta=(A*tau/2)+(A*tau)^2/12;
% Ta=(A*tau)+(A*tau)^2/2+(A*tau)^3/6+(A*tau)^4/24;
Ta=sparse(Ta);
% nearly pade approximation
B=(A*tau/2)-(A*tau)^2/12;
B=sparse(B);
Ta1=B+B^2+B^3+B^4+B^5+B^6;
Ta=Ta*Ta1+Ta+Ta1;
clear B Ta1 A Minv C;

% fast nearly symplectic precise integration algorithm
for i=1:N
    tic
    Ta=2*Ta+Ta*Ta;
    [m,n,src]=find(Ta);
    k=find(abs(src)<1e-15); % Matrix sparseness
    src(k)=0;
    Ta=sparse(m,n,src,2*number_of_notes,2*number_of_notes);  
    toc
     if i==N-2
        T3=Ta;
    end
end
T=I+Ta;        % 1
T3=I+T3;       % 1/4
clear Ta I;
c=[7/90;16/45;2/15;16/45;7/90];
numN=number_of_notes;
v=zeros(2*number_of_notes,1);
r=v;
r2=zeros(size(r));
r1=r2;
r3=r2;
r4=r3;
r0=r1;

tic

% disp=zeros(nt,100);
% disp=zeros(nt,1);
% isrc=ne/2;
for it =1:500   
   k=(it-1)*4+1;
  if  k<length(s_t)-4
    r0(numN+isrc)=mm*s_t (k);
    r1(numN+isrc)=mm*s_t (k+1);
    r2(numN+isrc)=mm*s_t (k+2);
    r3(numN+isrc)=mm*s_t (k+3);
    r4(numN+isrc)=mm*s_t (k+4);
    vf=dt*(c(1)*r4+T3*(c(2)*r3+T3*(c(3)*r2+T3*(c(4)*r1+T3*c(5)*r0))));
    v=T*v+vf;
  else
   v=T*v;

  end
 U=v(1:numN);
 if it==100
 save U100.mat U
 end
 if it==500
 save U500.mat U
 end

    trisurf(Top, p(:,1), p(:,2), U)
     colorbar
     clim([-1e-8,1e-8])
     tsr=num2str(it);
     title(tsr)
    shading interp
    colormap(gca,[seiscolor])
    view([90,90])
    drawnow() 
end


