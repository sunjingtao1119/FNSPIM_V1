function fac= lagrange(N,i,x)
%     """
%     Function to calculate  Lagrange polynomial for order N and polynomial
%     i[0, N] at location x.
%     """
[xi, ~] = gll(N);
    fac = 1;
    k=i+1;
    for j =0:N      
        m=j+1;
        if j ~= i
            fac = fac .* ((x - xi(m)) ./ (xi(k)- xi(m)));
        end
    end
end

