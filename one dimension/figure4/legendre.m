function p= legendre(N,x)
%Returns the value of Legendre Polynomial P_N(x) at position x[-1, 1]. 
% N  is the order of legendre polynomials
% x  is the value in [-1, 1]
% p  is the value of Legendre Polynomial P_N(x) at position x[-1, 1]. 
 P = zeros(2 * N,1);
  if N == 0
        P(1) = 1;
    elseif N == 1
        P(2) = x;
    else
        P(1) = 1;
        P(2)= x;
    for i= 2:N
        k=i+1;
        P(k) = (1.0 / i) * ((2 * i - 1) * x * P(k-1) - (i - 1) * P(k-2));
    end
   end
 p=P(N+1);
end
% Computational Seismology A Practical Introduction p196
