function SaveMatrix(S, filename)
    if(isfile(filename))
        delete(filename);
    end
    
    [row, col, val] = find(S);
    gridSize = size(S);

    h5create(filename, '/rows', length(row), 'Datatype', 'int32');
    h5create(filename, '/cols', length(col), 'Datatype', 'int32');
    h5create(filename, '/size', length(gridSize), 'Datatype', 'int32');

    h5write(filename, '/rows', int32(row-1));
    h5write(filename, '/cols', int32(col-1));
    h5write(filename, '/size', int32(gridSize));

end