function Ke=StifnessT3(point)
% Stiffness matrix of linear triangular element
%% input
%  point : The node coordinates of  element
%% output
%  Ke : Element stiffness matrix
        DN=zeros(2,3);
%       N1 = @(xi, eta) 1 - xi - eta;  % The shapefunction of node 1 
        N1_xi = -1; N1_eta = -1;       % The derivtive of shapefunction
        DN(:,1)=[N1_xi;N1_eta];
%       N2 = @(xi, eta) xi;            % The shapefunction of node  2 
        N2_xi = 1; N2_eta = 0;         % The derivtive of shapefunction
        DN(:,2)=[N2_xi;N2_eta];
%       N3 = @(xi, eta) eta;           % The shapefunction of node  3
        N3_xi = 0; N3_eta = 1;         % The derivtive of shapefunction
        DN(:,3)=[N3_xi;N3_eta];
%         ymax = @(xi) 1 - xi;
%       xx = point(:,1); yy = point(:,2);
       JN=[N1_xi,N2_xi,N3_xi;
           N1_eta,N2_eta,N3_eta];
       % The Jacobin matrix
       J=JN*point;
       detJ = abs(det(J));
       % The inverse of the Jacobin matrix
       J_inv_T =(1/detJ)* [J(2,2),-J(1,2);-J(2,1),J(1,1)];
       Ke = zeros(3, 3);
        for i = 1: 3
             for j = 1: 3
                    func=(J_inv_T*DN(:,j))'*(J_inv_T*DN(:,i));
                    %Ke(i,j)= func*detJ*0.5;
                    Ke(i,j)= func*detJ;
             end
        end

end