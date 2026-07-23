function M = readmatrix(fname, varargin)
% READMATRIX  Octave shim for MATLAB's readmatrix (numeric CSV).
%   M = readmatrix(FNAME) reads a comma-separated file that has one header
%   row of text followed by numeric data, and returns the numeric matrix.
%
%   This mirrors how the analysis scripts call readmatrix on the
%   data/*_interpolated.csv files. It exists because Octave (as of 8.x)
%   has no built-in readmatrix. Any extra arguments are ignored.
%
%   See also: dlmread, octave_compat/README.md
  M = dlmread(fname, ',', 1, 0);   % skip 1 header row, start at column 0
end
