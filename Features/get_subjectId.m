function [ Subject ] = get_subjectId_rerc( n,type )
% *************************************************************************
% get_subjectId.m
%
% This function takes a participant number and type (control or cva) and 
% converts it into the subject identifier unique to the study.
%
%     Input:   n       = subject number (integer)
%              type    = subject type ('controls' or 'cva')
%     Output:  Subject = subject ID from study (string)
%
% Examples: Healthy subject #1 = HC01, Stroke subject #20 = CVA20
%
% Megan O'Brien, 2018
% *************************************************************************

if isnumeric(n)
    if floor(n) ~= n
       error('Error. Input must be int.')
    end
    
    str = num2str( n );
    if length( str ) == 1 
        if strcmp(type,'cva')
            Subject=['CVA0' str];
        elseif strcmp(type,'controls')
            Subject=['HC0' str];
        else
            error('Subject type not recognized. Please use inputs ''controls'' or ''cva''\n');
        end
    else
        if strcmp(type,'cva')
            Subject=['CVA' str];
        elseif strcmp(type,'controls')
            Subject=['HC' str];
        else
            error('Subject type not recognized. Please use inputs ''controls'' or ''cva''\n');
        end
    end
    
else
    error('MyComponent:incorrectType',...
       'Error. \nInput must be int, not %s.',class(n))
end

end

