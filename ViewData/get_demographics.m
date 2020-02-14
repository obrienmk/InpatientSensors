function [ Dem_Out ] = get_demographics( Subject_type,dirDem )
%**************************************************************************
% get_demographics.m
%
% Written by Megan O'Brien
% version 9.18.2018
%
% Get demographic info for all subjects from Excel file
%
% Inputs:
%       Subject_type = 'cva' or 'controls'
%       dirDem = location of Excel file
%**************************************************************************

if strcmp(Subject_type,'cva')
    filename = 'CVA_Clinical Outcome Measure Scores Master Sheet.xlsx';
elseif strcmp(Subject_type,'controls')
    filename = 'HC_Clinical Outcome Measure Scores Master Sheet.xlsx';
else
    fprintf('Subject type not recognized. \n')
end

Demographics=table2cell(readtable([dirDem filename],'ReadRowNames',true));

Dem_Out = Demographics;

end
