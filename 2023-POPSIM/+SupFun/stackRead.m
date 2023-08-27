function [stack, stackinfo] = stackRead(stackpath)
%
% [stack, stackinfo] = stackRead(stackpath)
%
% Wrapper function for Fran�ois Nedelec's tiffread.m
% This works on STK files as well as multipage TIFFs.
% Outputs:
%
%   stack     : 3D data array
%   stackinfo : Any other information included in original file
%
% Francois Aguet, 01/2010

if (nargin == 0 || isempty(stackpath))
    stackinfo = tiffRead();
else
    stackinfo = SupFun.tiffRead(stackpath);
end

stack = stackinfo;%cat(3,stackinfo(:).data);
%stackinfo = rmfield(stackinfo, 'data');