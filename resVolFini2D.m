function [X, k, dx, dy] = resVolFini2D(Nx, Ny, k, kc, ki, Lx, Ly)
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
    for i = 1:Ny
        for j = 1:Nx
            i0 = i + (j-1)*Ny;  
            
            % 1. CONDITIONS DE DIRICHLET (Températures imposées Haut et Bas)
            if i == 1
                v = v+1;
                I(v) = i0;
                J(v) = i0;
                V(v) = 1; % Diagonale à 1, aucun voisin ajouté
                B(i0) = Tetc; 
                
            elseif i == Ny
                v = v+1;
                I(v) = i0;
                J(v) = i0;
                V(v) = 1; % Diagonale à 1, aucun voisin ajouté
                B(i0) = Tetf; 
                
            % 2. POINTS INTERNES ET PAROIS ADIABATIQUES (Gauche et Droite)
            else
                % i est forcément entre 2 et Ny-1, les voisins verticaux existent toujours
                Rim = dy/(2*dx) * (1/k(i-1, j) + 1/k(i, j));
                Rip = dy/(2*dx) * (1/k(i+1, j) + 1/k(i, j));
                
                % Les parois adiabatiques latérales restent gérées par les inf
                if j > 1
                    Rjm = dx/(2*dy) * (1/k(i, j-1) + 1/k(i, j));
                else
                    Rjm = inf;
                end
                
                if j < Nx
                    Rjp = dx/(2*dy) * (1/k(i, j+1) + 1/k(i, j));
                else
                    Rjp = inf;
                end
                
                % Diagonale principale
                v = v+1;
                I(v) = i0;
                J(v) = i0;
                V(v) = -1 * (1/Rim + 1/Rip + 1/Rjm + 1/Rjp);
                
                % Voisins verticaux (toujours présents)
                v = v+1; I(v) = i0; J(v) = i0-1; V(v) = 1/Rim;
                v = v+1; I(v) = i0; J(v) = i0+1; V(v) = 1/Rip;
                
                % Voisins horizontaux (ajoutés seulement si pas sur un bord)
                if j > 1
                    v = v+1; I(v) = i0; J(v) = i0-Ny; V(v) = 1/Rjm;
                end
                if j < Nx
                    v = v+1; I(v) = i0; J(v) = i0+Ny; V(v) = 1/Rjp;
                end
            end
        end
    end
    
    % CONDITION LIMITE DIRCHLET (température imposée)
    % Paroi haut

    A = sparse(I(1:v), J(1:v), V(1:v), N, N);    
    
    %Résolution
    X=A\B;
end