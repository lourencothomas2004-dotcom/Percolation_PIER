function S=LoadMatrix(filename)
    row = h5read(filename, '/rows') + 1;
    col = h5read(filename, '/cols') + 1;
    gridSize = h5read(filename, '/size');

    S=sparse(row, col, ones(length(row), 1), gridSize(1), gridSize(2));
end