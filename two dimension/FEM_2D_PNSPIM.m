clear;clc
%% 生成单元
seiscolor=load('colorbar.txt');
le=2000;  % size 2000*2000
ne=200;   % 
[x,y]=meshgrid(linspace(0,le,ne));
de=le/ne;
X=x(:);
Y=y(:);
P=[X,Y];
DT = delaunayTriangulation(P);
F = freeBoundary(DT);
b=unique(F(:)); % 边界节点
%生成结构化三角网格
p=DT.Points;
% p=p';
Top=DT.ConnectivityList;
% t=t';
num_nodes = length(p);
num_elements = length(Top);
t_s = sprintf('nodes=%d  elements=%d',num_nodes,num_elements);
title(t_s)
%% 参数设置
nt = 1500;
v = 1500; 
ept=de/v;
dt=0.5*ept;
T = (1:4*nt)*dt/4;
number_of_notes = length(p);
fmain =40;
%% 设置震源
s_t = (1-2*pi^2*fmain*(T-0.2).^2).*exp(-fmain*pi^2*(T-0.2).^2);%雷克子波
fricker = @(x) (1-2*pi^2*fmain*(x-0.2).^2).*exp(-fmain*pi^2*(x-0.2).^2);

%% 调用函数计算矩阵
% FEM = FEM_2D_func();
% S = FEM.assemble_matrix_2D(p, t);%刚度矩阵
% Top=t';
% pp=p';

dictK= TStiffness(p,Top); % 

dictM= TMass(p,Top); % 

K=dic2sp2(dictK,number_of_notes,number_of_notes);

clear dictK;

Minv=dic2sp2(dictM,number_of_notes,number_of_notes,1);

clear dictM;
Fr = freeBoundary(DT);
boundarynodes=unique(Fr(:)); % 边界节点
U = zeros(number_of_notes,nt);
source_x= find(abs(p(:,1))>=le/2-de&abs(p(:,1))<=le/2+de&abs(p(:,2))>=le/2-de&abs(p(:,2))<=le/2+de);

% 精细积分法
% L =Minv*K.*v^2;
% [m1,n1,src1]=find(Minv);
% n1=n1+num_nodes;
% I1=speye(num_nodes,num_nodes);
% [m2,n2,src2]=find(K);
% m2=m2+num_nodes;
% src2=src2.*v^2.*-1;
L=Minv*K.*v^2;
I1=speye(number_of_notes,number_of_notes);
[m1,n1,src1]=find(I1);
n1=n1+number_of_notes;
[m2,n2,src2]=find(L);
m2=m2+number_of_notes;
src2=-1.*src2;
m=[m1;m2];
n=[n1;n2];
src=[src1;src2];
A=sparse(m,n,src,2*num_nodes,2*num_nodes);  
% A=[ -Minv*C./2,Minv;C*Minv*C./4-v^2*S,-C*Minv./2];
%   [m,n,src]=find(A1);
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

isrc=source_x(1);
ms=Minv(isrc,isrc);
clear B Ta1 A Minv C;
% Ta=full(Ta);
for i=1:N
    tic
    Ta=2*Ta+Ta*Ta;
    [m,n,src]=find(Ta);
    k=find(abs(src)<1e-15);
    src(k)=0;
    Ta=sparse(m,n,src);  
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
isrc=source_x(1);
tic

for it =1:100   
   k=(it-1)*4+1;
  if  k<length(s_t)-4
    r0(numN+isrc)=ms*s_t (k);
    r1(numN+isrc)=ms*s_t (k+1);
    r2(numN+isrc)=ms*s_t (k+2);
    r3(numN+isrc)=ms*s_t (k+3);
    r4(numN+isrc)=ms*s_t (k+4);
    vf=dt*(c(1)*r4+T3*(c(2)*r3+T3*(c(3)*r2+T3*(c(4)*r1+T3*c(5)*r0))));
    v=T*v+vf;
%     v(boundarynodes) = 0;
  else
   v=T*v;
%    v(boundarynodes) = 0;
  end
% U=v(1:numN);
% trisurf(Top, p(:,1), p(:,2), U)
%      colorbar
% clim([-1e-8,1e-8])
%   shading interp
%     colormap(gca,[seiscolor])
%     view([0,90])
%     drawnow() 
end
toc





