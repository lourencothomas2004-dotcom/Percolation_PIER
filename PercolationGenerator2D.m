% Nx and Ny (int) are the matrix dimension
% r (int) is the percolation load's radius
% fs (float) is the surface ratio of percolation
% OverL (bool) is we want overlap or not

function [S, tempFs, filename]=PercolationGenerator(Nx, Ny, Lx, Ly, r, fs, OverL, save)
    maxAttempt = 1000;
    S = false(Nx, Ny);
    
    x_vec = linspace(0, Ly, Ny);
    y_vec = linspace(0, Lx, Nx);
    [X, Y] = meshgrid(x_vec, y_vec);
    tempFs = 0;
    surface = (Nx*Ny);
    
    attempt = 0;
    while tempFs<fs && attempt<maxAttempt
        x = rand()*Ly;
        y = rand()*Lx;

        col_idx = find(x_vec >= (x - r) & x_vec <= (x + r));
        row_idx = find(y_vec >= (y - r) & y_vec <= (y + r));
        X_sub = X(row_idx, col_idx);
        Y_sub = Y(row_idx, col_idx);

        TempCircle = ((X_sub-x).^2 + (Y_sub-y).^2) <= r^2; %((X-x).^2<=r^2 & (Y-y).^2<=r^2);

        Overlap = S(row_idx, col_idx) & TempCircle;
        if ~any(Overlap(:)) || OverL
            S(row_idx, col_idx) = S(row_idx, col_idx) |TempCircle;
            tempFs = tempFs + (sum(TempCircle(:)) - sum(Overlap(:)))/surface;

            attempt=0;
        else
            attempt = attempt + 1;
        end
    end

    S = sparse(S);   
    filename = sprintf('Matrices/%d_%d_%d_0%.0f.h5', Nx, Ny, r, int16(tempFs*100));
    if save    
        SaveMatrix(S,filename);
    end
end
