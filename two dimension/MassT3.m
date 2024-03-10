function Me=MassT3(point)
% Stiffness matrix of linear triangular element
%% input
%  point : The node coordinates of  element
%% output
%  Ke : Element stiffness matrix

        N{1} = @(xi, eta) 1 - xi - eta;  % The shapefunction of node 1 
        N1_xi = -1; N1_eta = -1;       % The derivtive of shapefunction
       
        N{2} = @(xi, eta) xi;            % The shapefunction of node  2 
        N2_xi = 1; N2_eta = 0;         % The derivtive of shapefunction
      
        N{3} = @(xi, eta) eta;           % The shapefunction of node  3
        N3_xi = 0; N3_eta = 1;         % The derivtive of shapefunction
     
%         ymax = @(xi) 1 - xi;
%       xx = point(:,1); yy = point(:,2);
       JN=[N1_xi,N2_xi,N3_xi;
           N1_eta,N2_eta,N3_eta];
       % The Jacobin matrix
       J=JN*point;
       detJ = abs(det(J));
       % The inverse of the Jacobin matrix
       Me = zeros(3, 3);

        for i = 1: 3
             for j = 1: 3
                    func=N{j}(1/3,1/3)*N{i}(1/3,1/3);
                    Me(i,j)= func*detJ*0.5;
 
             end
        end

end