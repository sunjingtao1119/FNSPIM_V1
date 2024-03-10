function Keys = Top2KeysM(top)
% The index of mass matrix
n=length(top);
Keys=strings(1,n);
  for i=1:n
    Keys(i)=string(top(i))+","+string(top(i));
  end
end

