function Keys = Top2Keys(top)
% Matlab does not support double key value dictionary
% The Top matrix of  element is transformed into dictionary keys of 
% the global sparse stiffness matrix .
n=length(top);
Keys=strings(1,n*n);
k=1;
for i=1:n
    for j=1:n
        Keys(k)=string(top(i))+","+string(top(j));
        k=k+1;
    end
end
end

