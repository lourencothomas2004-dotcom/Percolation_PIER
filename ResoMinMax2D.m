function [fmin, keqMin, fmax, keqMax] = ResoMinMax2D(Lx, Ly, r, fs, kc, ki, Nx, Ny)

    [kMin, fmin, filename]= PercolationGenerator2Dmin(Nx, Ny, Lx, Ly, r, fs, false, false);

    [Xmin, kmin, dx, dy] = resVolFini2D(Nx, Ny, kMin, kc, ki, Lx, Ly);

    Tet=reshape(Xmin,Ny,Nx);
    TetMin = flipud(Tet); 
    kmin = flipud(kmin);

    Phi = mean(sum((1./(1./kmin(1:Ny-1,:) + 1./kmin(2:Ny,:))).*(TetMin(2:Ny,:)-TetMin(1:Ny-1,:))*(2*dx/dy),2));

    keqMin = Phi * (Nx-1) / Ny;

    [kMax, fmax, filename]= PercolationGenerator2Dmax(Nx, Ny, Lx, Ly, r, fs, false, false);
    [Xmax, kmax, dx, dy] = resVolFini2D(Nx, Ny, kMax, kc, ki, Lx, Ly);

    Tet=reshape(Xmax,Ny,Nx);
    TetMax = flipud(Tet); 
    kmax = flipud(kmax);

    Phi = mean(sum((1./(1./kmax(1:Ny-1,:) + 1./kmax(2:Ny,:))).*(TetMax(2:Ny,:)-TetMax(1:Ny-1,:))*(2*dx/dy),2));

    keqMax = Phi * (Nx-1) / Ny;
end