function [f, keq] = Reso1Cas2D(Lx, Ly, r, fs, kc, ki, Nx, Ny)

    [kl, f, filename]= PercolationGenerator2D(Nx, Ny, Lx, Ly, r, fs, false, false);
    [X, k, dx, dy] = resVolFini2D(Nx, Ny, kl, kc, ki, Lx, Ly);

    Tet=reshape(X,Ny,Nx);
    Tet = flipud(Tet); 
    k = flipud(k);

    Phi = mean(sum((1./(1./k(1:Ny-1,:) + 1./k(2:Ny,:))).*(Tet(2:Ny,:)-Tet(1:Ny-1,:))*(2*dx/dy),2));

    keq = Phi;
end