function [] = createdir(dir)
%createdir Create directory is does not exist
    if ~exist(dir, 'dir')
        mkdir(dir); 
    end

end

