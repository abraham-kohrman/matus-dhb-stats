function [newData1] = importXLSfile(fileToRead1)
%IMPORTFILE(FILETOREAD1)
%  Imports data from the specified file
%  FILETOREAD1:  file to read

%  Auto-generated by MATLAB on 17-Sep-2019 21:16:49

% Import the file
[~, sheetName] = xlsfinfo(fileToRead1);
newData1 =containers.Map;
for kay = 1:numel(sheetName)
    [numbers, strings] = xlsread(fileToRead1, sheetName{kay});
    if ~isempty(numbers)
        newData1(sheetName{kay}) =  numbers;
    end
    if ~isempty(strings)
        %newData1(sheetName{kay}) =  strings;
        disp(strings);
    end
end
end
