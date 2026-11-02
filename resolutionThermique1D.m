function [Tet, k, dx]=resolutionThermique1D(Nx, k, kc, ki, Lx)

    Tetc = 1; %Température haut
    Tetf = 0; %Température bas
    dx=Lx/(Nx-1); %Pas x

    k = ki + k*(kc-ki); % On remplit le vecteur avec les valeurs de conductivité
    
    % Construction matrice A
    
    A = spalloc(Nx,Nx,5*Nx);
    B = zeros(Nx,1);

    for i = 2:Nx-1
    A(i,i-1) = k(i)/dx^2;
    A(i,i) = -(k(i)+k(i+1))/dx^2;
    A(i,i+1) = k(i+1)/dx^2;

    end
 
    % CL Dirichlet

    % Paroi gauche ( Tet = Tetc)
    A(1,1) = 1;
    B(1)  = Tetc;

    % Paroi droite ( Tet = Tetf)
    A(Nx,Nx) = 1;   
    B(Nx) = Tetf;   

    %Résolution

    Tet = A\B;

end




   



