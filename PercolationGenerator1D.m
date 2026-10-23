% Nx and Ny (int) are the matrix dimension
% r (int) is the percolation load's radius
% fs (float) is the surface ratio of percolation
% OverL (bool) is we want overlap or not

function [S, tempFs]=PercolationGeneratorD(Nx, r, fs, OverL)
    maxAttempt = 1000;
    S = zeros(1, Nx);

    Y = 1:Nx;
    tempFs = 0;
    surface = (Nx);
    
    attempt = 0;
    while tempFs<fs && attempt<maxAttempt
        y = randi([1,Nx]);

        TempCircle = ((Y-y).^2) <= r^2;

        Overlap = S & TempCircle;
        if ~any(Overlap(:)) || OverL
            S = S |TempCircle;
            tempFs = tempFs + (sum(TempCircle(:)) - sum(Overlap(:)))/surface;

            attempt=0;
        else
            attempt = attempt + 1;
        end
    end
end


S = PercolationGeneratorD(512, 1, 0.2, false);
figure;
imagesc(S);
colormap(hot)