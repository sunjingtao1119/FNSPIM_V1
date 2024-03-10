clear;clc
%% color spectrum of seismic 
seiscolor=load('colorbar.txt');

%% Generate triangular mesh
lex=4600;   % Length in x direction
ley=2500;   % Length in y direction
nex=460;   % The number of units in x direction
ney=250;   % The number of units in x direction
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
%% velocity model
seismic=read_segy_file('MODEL_P-WAVE_VELOCITY_1.25m.segy');
data=seismic.traces;
[m,n]=size(data);
x=linspace(1,n,nex);
x=fix(x);
y=linspace(1,m,ney);
y=fix(y);
V=data(y,x);
V=V(:);
F = freeBoundary(DT);
b=unique(F(:));    % boundary node
%% Parameter Setting
nt = 500;
v=max(V);
de=min([dex,dey]);
ept=dex/v;
dt=0.4*ept;
T = (1:4*nt)*dt/4;
number_of_notes = length(ppx);
fmain =40;
%% Setting the source
s_t = (1-2*pi^2*fmain*(T-0.2).^2).*exp(-fmain*pi^2*(T-0.2).^2);%雷克子波
fricker = @(x) (1-2*pi^2*fmain*(x-0.2).^2).*exp(-fmain*pi^2*(x-0.2).^2);

% 顶部边界
[rec,m] = find(abs(p(:,2))==0);
b=setdiff(b,rec);
%% Assembly stiffness matrix
dictK= TStiffness(p,Top,V); % 
dictM= TMass(p,Top); % 
K=dic2sp2(dictK,number_of_notes,number_of_notes);
K(b,:)=0;
K(:,b)=0;
K(b,b)=1;
clear dictK;
Minv=dic2sp2(dictM,number_of_notes,number_of_notes,1);
Minv(b,:)=0;
Minv(:,b)=0;
% K(b,b)=1;
clear dictM;
Fr = freeBoundary(DT);
boundarynodes=unique(Fr(:)); % 边界节点
U = zeros(number_of_notes,nt);
% source_x= find(abs(p(:,1))>=lex/2-de&abs(p(:,1))<=lex/2+de&abs(p(:,2))>=ley/2-de&abs(p(:,2))<=ley/2+de);
[source_x,m] = find(abs(p(:,2))==0&abs(p(:,1))>=lex/2-dey&abs(p(:,1))<=lex/2+dey);
[rec,m] = find(abs(p(:,1))>=1000-dex&abs(p(:,1))<=1000+dex&abs(p(:,2))>=50-dey&abs(p(:,2))<=50+dey);

isrc=source_x(1);
%% 利用递推关系求波场值
L=Minv*K;
mm=Minv(isrc,isrc);
I1=speye(number_of_notes,number_of_notes);
[m1,n1,src1]=find(I1);
n1=n1+number_of_notes;
[m2,n2,src2]=find(L);
m2=m2+number_of_notes;
src2=-1.*src2;
% L=Minv*K;
% I1=speye(number_of_notes,number_of_notes);
% [m1,n1,src1]=find(I1);
% n1=n1+number_of_notes;
% [m2,n2,src2]=find(L);
% m2=m2+number_of_notes;
% src2=-1.*src2;
m=[m1;m2];
n=[n1;n2];
src=[src1;src2];
A=sparse(m,n,src,2*number_of_notes,2*number_of_notes);  
I=speye(size(A));
N=10;
m1=2^N;
tau=dt/m1;
Ta=(A*tau/2)+(A*tau)^2/12;
Ta=sparse(Ta);
%利用幂级数近似
B=(A*tau/2)-(A*tau)^2/12;
B=sparse(B);
Ta1=B+B^2+B^3+B^4+B^5+B^6;
Ta=Ta*Ta1+Ta+Ta1;
clear B Ta1 A Minv C I1 L;
% Ta= gpuArray(Ta);
for i=1:N
    tic
    Ta=2*Ta+Ta*Ta;
    [m,n,src]=find(Ta);
    k=find(abs(src)<1e-15);
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
disp=zeros(nt,1);
mm=full(mm);
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
%     v(boundarynodes) = 0; 
    else
    v=T*v;
%     v(boundarynodes) = 0;
  end
    U=v(1:numN);
%     disp(it)=U(rec(1));
%     trisurf(Top, p(:,1), p(:,2), U)
%     xlim([0,lex])
%     ylim([0,ley])
%     colorbar
%     clim([-1e-8,1e-8])
%     colormap(gca,[seiscolor])
%    shading interp
%     view([0,-90]);
%     dx=lex*5/1000;
%     dy=ley*5/1000;
%     set(gcf,'unit','centimeters','position',[10 5 dx dy]);
%     drawnow() 
end
save MSpimu500.mat U
% save mm.mat mm

