% Nx and Ny (int) are the matrix dimension
% r (int) is the percolation load's radius
% fs (float) is the surface ratio of percolation
% OverL (bool) is we want overlap or not

function [S, tempFs, filename]=PercolationGenerator(Nx, Ny, r, fs, OverL)
    maxAttempt = 1000;
    S = false(Nx, Ny);

    [X, Y] = meshgrid(1:Ny, 1:Nx);
    tempFs = 0;
    surface = (Nx*Ny);
    
    attempt = 0;
    while tempFs<fs && attempt<maxAttempt
        x = randi([1,Ny]);
        y = randi([1,Nx]);

        TempCircle = ((X-x).^2 + (Y-y).^2) <= r^2;

        Overlap = S & TempCircle;
        if ~any(Overlap(:)) || OverL
            S = S |TempCircle;
            tempFs = tempFs + (sum(TempCircle(:)) - sum(Overlap(:)))/surface;

            attempt=0;
        else
            attempt = attempt + 1;
        end
    end
    
    S = sparse(S);
    filename = sprintf('Matrices/%d_%d_%d_0%.0f.h5', Nx, Ny, r, int16(tempFs*100));

    SaveMatrix(S,filename);
end
