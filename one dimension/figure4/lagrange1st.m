function out = lagrange1st(N)
%     """
%     # Calculation of 1st derivatives of Lagrange polynomials
%     # at GLL collocation points
%     # out = legendre1st(N)
%     # out is a matrix with columns -> GLL nodes
%     #                        rows -> order
%     """
out=zeros(N+1,  N+1);
[xi, ~] = gll(N);
% # initialize dij matrix (see Funaro 1993)
    d =zeros([N + 1, N + 1]);

    for i =0:N
        k=i+1;
        for j = 0:N
            m=j+1;
            if i ~= j
                d(k,m)= legendre(N, xi(k)) /legendre(N, xi(m)) *1.0 / (xi(k)-xi(m));
            end
            if i == 0
                if j == 0
                    d(k, m)= -1.0 / 4.0 * N * (N + 1);
                end
            end
            if i == N
                if j == N
                    d(k, m)= 1.0 / 4.0 * N * (N + 1);
                end
            end
        end
     end
  % Computational Seismology A Practical Introduction p196
%   Calculate matrix with 1st derivatives of Lagrange polynomials
    for n=0:N
        l=n+1;
        for i=0:N
            k=i+1;
            sum = 0;
            for j =0:N
                m=j+1;
                    sum = sum + d(k, m) * lagrange(N, n, xi(m));
            end
            out(l,k)= sum;
        end
    end
end

