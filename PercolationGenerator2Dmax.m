% CAS LE PLUS CONDUCTEUR : lignes de charges verticales (pleine hauteur)
% Nx and Ny (int) are the matrix dimension
% Lx and Ly (float) are the domain lengths
% r (float) is the half-width of a line (a line is 2*r wide)
% fs (float) is the surface ratio of percolation
% OverL (bool) is we want overlap or not

function [S, tempFs, filename]=PercolationGenerator2Dmax(Nx, Ny, Lx, Ly, r, fs, OverL)
    maxAttempt = 1000;
    S = false(Nx, Ny);
    X = linspace(0, Ly, Ny);
    Y = linspace(0, Lx, Nx);
    [X, Y] = meshgrid(X, Y);
    tempFs = 0;
    surface = (Nx*Ny);
    attempt = 0;
while tempFs < fs - 0.5/Ny && attempt<maxAttempt
        x = rand()*Ly;
        rr = min(r, (fs-tempFs)*Ly/2);   % la dernière ligne est raccourcie pour tomber sur fs
        TempLine = abs(X-x) <= rr;       % bande verticale sur toute la hauteur
        Overlap = S & TempLine;
if (~any(Overlap(:)) || OverL) && any(TempLine(:))
            S = S |TempLine;
            tempFs = tempFs + (sum(TempLine(:)) - sum(Overlap(:)))/surface;
            attempt=0;
else
            attempt = attempt + 1;
end
end
    S = sparse(S);
    filename = sprintf('Matrices/max_%d_%d_%d_0%.0f.h5', Nx, Ny, r, int16(tempFs*100));
    SaveMatrix(S,filename);
end