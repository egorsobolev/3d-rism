function x=matrix_read(fname, type)

if nargin < 1
    error('MATLAB::badopt', 'Error using');
    return;
elseif nargin < 2
    type='double';
end    

fid=fopen(fname, 'r');

[S, cnt]=fread(fid, 2, 'uint');
[x, cnt]=fread(fid, [S(1), S(2)], type);

fclose(fid);
return;