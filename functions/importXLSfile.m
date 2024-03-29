%filetoRead1 = '/home/posfailab/Documents/emergencyMeeting/Gata_Spreadsheet1.xls'; %full path to the data file, including the filename of the data.xlsx file. 
%importXLSssfile(filetoRead1)
function [newData1] = importXLSfile(fileToRead1)
%IMPORTFILE(FILETOREAD1)
%  Imports data from the specified file
%  FILETOREAD1:  file to read

%  Auto-generated by MATLAB on 17-Sep-2019 21:16:49

% Import the file
%[~, sheetName] = xlsfinfo(fileToRead1);
sheetName = sheetnames(fileToRead1);
newData1 =containers.Map;
for kay = 1:numel(sheetName)
    %[numbers, strings] = xlsread(fileToRead1, sheetName(kay));
    numbers = table2array(readtable(fileToRead1,"Sheet",sheetName(kay)));
    if ~isempty(numbers)
        newData1(sheetName{kay}) =  numbers;
    end
    if ~isempty(strings)
        %newData1(sheetName{kay}) =  strings;
        disp(strings);
    end
end
end

