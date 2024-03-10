function dictM= TMass(p,Top)
% input:
% p：The  total node 
% Top：The element matirx
% output
% dicK：The dictionary of stiffness matrix;
[Ne,~]=size(Top); %Ne : the number of element , nn the number of element

for  kk=1:Ne
    id=Top(kk,:);         
    names=Top2KeysM(id);  
    point=p(id,:);
    Me=MassT3(point);
    Me=sum(Me,"all")/sum(diag(Me))*diag(Me);% 集中质量矩阵
    values=Me';
    if kk==1
        dictM=dictionary(names,values);  
    else
       dictM=sumdict(dictM,names,values);
    end

end

end

