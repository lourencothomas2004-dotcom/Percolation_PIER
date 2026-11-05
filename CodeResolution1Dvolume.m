%----------- Conduction thermique 1D, cas hétérogène ------------------

%------------------- Méthode des volumes finis -------------------------


clear all;
close all;

Tetc = 1; % Température normalisée paroi gauche (Tet = 1)
Tetf = 0; % Température normalisée paroi droite (Tet = 0)
Lx = 1; % Longueur
Nx = 512; % Nb de points
dx = Lx/(Nx-1); % Pas
r = 6; % Rayon des charges (en points)
fs = 0.2; % Fraction de charges
kc = 1000; % Conductivité charges
ki = 1; % Conductivité isolant

[S, Fs] = PercolationGenerator1D(Nx, r, fs, false);
k = ki + double(S(:))*(kc-ki);    % On remplit le vecteur des conductivités


% Coefficient alpha_i

alpha = k(1:Nx-1).*k(2:Nx)./(k(1:Nx-1)+k(2:Nx));   

% Initialisation matrice de coefficients A et second membre B

A = spalloc(Nx,Nx,3*Nx);
B = zeros(Nx,1);

% Points internes 
for i = 2:Nx-1
    A(i,i-1) = alpha(i-1);
    A(i,i) = -alpha(i-1) - alpha(i);
    A(i,i+1) = alpha(i);
end

% CL Dirichlet

% Paroi gauche : Tet = Tetc

A(1,1) = 1;     
B(1) = Tetc; 

% Paroi doite : Tet = Tetf

A(Nx,Nx) = 1;   
B(Nx) = Tetf;

% Résolution
Tet = A\B;

% Tracé
x = 0:dx:Lx;
figure;
subplot(2,1,1)
plot(x, Tet, 'LineWidth', 1.5)
title('Profil de température normalisé')
ylabel('T^*'); grid on
subplot(2,1,2)
plot(x, k, 'LineWidth', 1.5)
title('Conductivité'); xlabel('x (en m)'); ylabel('k'); grid on

% Flux aux interfaces i+1/2 

Phi = 2*dx*alpha.*(Tet(1:Nx-1)-Tet(2:Nx));

figure;
plot(x(1:Nx-1)+dx/2, Phi)
ylim([0.9*mean(Phi) 1.1*mean(Phi)])
xlabel('x (m)'); ylabel('Flux'); title('Flux aux interfaces')

% Conductivité effective et vérification 

k_eff = (mean(Phi)*Lx)/(dx^2*(Tetc-Tetf));
k_eff_th = Lx/(dx*(sum(1./k) - 0.5/k(1) - 0.5/k(Nx)));
fprintf('Fs = %.3f , k_eff = %.6f , k_eff théorique = %.6f\n', Fs, k_eff, k_eff_th)

