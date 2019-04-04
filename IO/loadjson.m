function [valueStruct] = loadjson(filename)
fid = fopen(filename); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
valueStruct = jsondecode(str);
end

