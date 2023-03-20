function [ind] = sub2ind(row, col, n)
    % Row-wise transformation
    ind = (row-1)*n+col;
end

