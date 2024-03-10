function dictK= TStiffness(p,Top,v)
% input:
% p：The  total node 
% Top：The element matirx
% output
% dicK：The dictionary of stiffness matrix;
if nargin<3
    v=1;
end
[Ne,~]=size(Top); %Ne : the number of element , nn the number of element
if v==1
   v=ones(Ne,3);
end

for  kk=1:Ne
    id=Top(kk,:); 
    Ve=sum(v(id))/3;
    names=Top2Keys(id);  
    point=p(id,:);
    Ke=StifnessT3(point);
    Ke=Ke.*Ve^2;
    Ke=Ke';
    values=Ke(:);
    values=values';
    if kk==1
        dictK=dictionary(names,values);  
    else
       dictK=sumdict(dictK,names,values);
    end

end

end

