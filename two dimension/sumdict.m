function  dict=sumdict(dict,key,value)
% If the key value exists, the values are added
% If the key value does not exist, add the key and value
tf=isKey(dict,key);
n=length(tf);
for i=1:n
    k=key(i);
    if tf(i)
        dict(k)=dict(k)+value(i);
    else
        dict(k)=value(i); 
    end
end
end