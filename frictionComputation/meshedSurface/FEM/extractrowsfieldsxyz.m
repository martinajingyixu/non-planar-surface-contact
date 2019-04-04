function [coordOutput] = extractrowsfieldsxyz(coord,idxRow)
%applytofieldsxyz apply the certain function to all x y z fields of the
%input struct coord and return the struct coordOutput with same fields
%   Detailed explanation goes here
coordOutput.x = extractrows(coord.x,idxRow);
coordOutput.y = extractrows(coord.y,idxRow);
coordOutput.z = extractrows(coord.z,idxRow);

end

function [extractedRows] = extractrows(matrix,idxRow)
%extractelements Extract rows of a matrix with given row idx
    extractedRows = matrix(idxRow,:);
end

