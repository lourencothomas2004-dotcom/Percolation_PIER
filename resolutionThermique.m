function [X, k, dx, dy]=resolutionThermique(Nx, Ny, k, kc, ki, Lx, Ly)
    Tetc = 1; %Température haut
    Tetf = 0; %Température bas
    N = Nx*Ny;

    dx=Lx/(Nx-1); %Pas x
    dy=Ly/(Ny-1); %Pas y
    
    
    k = ki + k*(kc-ki); % On remplit la matrice avec les valeurs de conductivité
    
    % Construction matrice A
    
    I = zeros(1, 2*N);
    J = zeros(1, 2*N);
    V = zeros(1, 2*N);
    B = zeros(N,1);
    
    %Points internes
    v = 0;
    for i = 2:Ny-1
        for j = 2:Nx-1
            i0 = i+ (j-1)*Ny;  
            
            % Diagonale principale (Correcte)
            v = v+1;
            I(v) = i0;
            J(v) = i0;
            V(v) = -(k(i,j+1)/dx^2 + k(i,j)/dx^2 + k(i+1,j)/dy^2 + k(i,j)/dy^2);
            
            % Voisin HAUT (i-1) -> Déplacement Y (on utilise dy et la conductivité du haut/centre)
            v = v+1;
            I(v) = i0;
            J(v) = i0-1;
            V(v) = k(i,j)/dy^2; 
            
            % Voisin BAS (i+1) -> Déplacement Y (on utilise dy et la conductivité du bas)
            v = v+1;
            I(v) = i0;
            J(v) = i0+1;
            V(v) = k(i+1,j)/dy^2;
            
            % Voisin DROITE (j+1) -> Déplacement X (on utilise dx et la conductivité de droite)
            v = v+1;
            I(v) = i0;
            J(v) = i0+Ny;
            V(v) = k(i,j+1)/dx^2;
            
            % Voisin GAUCHE (j-1) -> Déplacement X (on utilise dx et la conductivité de gauche/centre)
            v = v+1;
            I(v) = i0;
            J(v) = i0-Ny;
            V(v) = k(i,j)/dx^2;
        end
    end
    
    %CONDITIONS LIMITES
    
    % CONDITION ADIABATIQUE dT/dy=0 ( paroi gauche et droite)
    
    % ORDRE 1 1*T(2:Nx-1,1)-1*T(2:Nx-1,2)=0
    
    % Paroi gauche ( T(i,1) - T(i,2) = 0 )
    
    for i = 2:Ny-1
        i0 = i;

        v = v+1;
        I(v) = i0;
        J(v) = i0;
        V(v) = 1;

        v = v+1;
        I(v) = i0;
        J(v) = i0+Ny;
        V(v) = -1;

        B(i0) = 0;
    end
    
    % Paroi droite ( T(i,Nx) - T(i,Nx-1) = 0 )
    
    for i = 2:Ny-1
        i0 = i + (Nx-1)*Ny;

        v = v+1;
        I(v) = i0;
        J(v) = i0;
        V(v) = 1;
        
        v = v+1;
        I(v) = i0;
        J(v) = i0-Ny;
        V(v) = -1;

        B(i0) = 0;
    end
    
    % CONDITION LIMITE DIRCHLET (température imposée)
    
    % Paroi haut

    A = sparse(I(1:v), J(1:v), V(1:v), N, N);
    for j = 1:Nx
        i0 = 1 + (j-1)*Ny;
        
        A(i0,:) = 0;
        A(i0,i0) = 1;
        B(i0) = 1; % Température Tc imposée paroi haut
    end
    
    % Paroi bas
    
    for j = 1:Nx
        i0 = Ny + (j-1)*Ny;
        A(i0,:) = 0;
        A(i0,i0) = 1;
        B(i0) = 0; % Température Tf imposée paroi bas
    end
    
    
    %Résolution
    X=A\B;

end