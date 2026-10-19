% Nx and Ny are the matrix dimension
% r is the percolation load's radius
% fs is the surface ratio of percolation

function [S, tempFs, filename]=GeneratePercolation(Nx, Ny, r, fs)
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
        if ~any(Overlap(:))
            S = S |TempCircle;
            tempFs = tempFs + sum(TempCircle(:))/surface;

            attempt=0;
        else
            attempt = attempt + 1;
        end
    end
    
    S = sparse(S);
    filename = sprintf("%d_%d_%d_%.0f.h5", Nx, Ny, r, int16(tempFs*100));

    SaveMatrix(S,filename);
end
