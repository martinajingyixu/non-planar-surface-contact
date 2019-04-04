function [ rotation_matrix ] = getCoordTransMat( x_old, y_old, z_old, x_new,y_new,z_new )
%compute transformation matrix given old x,y,z and new x,y,z
%   size of x,y,z old and new are 3x1
    if size(x_old,1) ~= 3 || size(x_old,2) ~= 1 || size(x_new,1) ~= 3 || size(x_new,2) ~= 1
        error('old or new x axis wrong dimention!')
    end         
        x_old = x_old./norm(x_old);
        y_old = y_old./norm(y_old);
        z_old = z_old./norm(z_old);
        x_new = x_new./norm(x_new);
        y_new = y_new./norm(y_new);
        z_new = z_new./norm(z_new);

        new_coord = [x_new,y_new,z_new];
        old_coord = [x_old,y_old,z_old];
        
        rotation_matrix = getRotationMatrix(old_coord,new_coord);

end

