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

% 1. Define your parameters
Nx = 512;       % Height
Ny = 512;       % Width
r = 7;         % Radius of inclusions
fs = 0.25;      % 40% surface filling

% 2. Generate the matrix
% (Make sure your GeneratePercolation and SaveMatrix functions are in the same folder)
[S, f, filename] = GeneratePercolation(Nx, Ny, r, fs);

% 3. Visualize it
figure;
imagesc(full(S));
colormap(hot); % Gives a nice thermal color theme (black to white/yellow)
axis equal;
axis tight;
title(['Percolation Grid: fs = ', num2str(f)]);