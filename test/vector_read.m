function x=matrix_read(fname, type)

if nargin < 1
    error('MATLAB::badopt', 'Error using');
    return;
elseif nargin < 2
    type='double';
end    

fid=fopen(fname, 'r');

[S, cnt]=fread(fid, 1, 'uint');
[x, cnt]=fread(fid, [S(1), 1], type);

fclose(fid);
return;