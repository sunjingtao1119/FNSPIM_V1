% Assembly stiffness matrix based on sparse dictionary
function A= dic2sp2(dict,m,n,inv)
% input
% dic: 
% m:row
% n:column
% inv: 1 Inverse moment of mass matrix
% output
% A：sparse matrix
if nargin<4
    inv=0;
end
names=dict.keys;
if inv==1
   data=1./dict.values;
else
   data=dict.values;
end
len=length(names);
index=ones(len,2);
for i=1:len
    tindex=names(i);
    tindex=strsplit(tindex,',');
    for j=1:2
    index(i,j)=str2num(char(tindex(j)));
    end
end
A=sparse(index(:,1),index(:,2),data,m,n);
end

