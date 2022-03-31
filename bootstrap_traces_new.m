close all;

%% read the stuff in
tic
replace = 0;
allData = importXLSfile('/home/posfailab/Documents/emergencyMeeting/Nanog_Spreadsheet.xls'); %full path to the data file, including the filename of the data.xlsx file. 

basedir = '/home/posfailab/Documents/emergencyMeeting/'; %path to the folder where the comparison_table.xlsx is located. Also the directory that will contain the outputs!

time_stamp = datestr(datetime('now'));
time_stamp = strrep(time_stamp,':','.'); 
savedir = [basedir time_stamp '/'];
mkdir(savedir);
%% start here if resetting to avoid reloading large datasets into memory.  
[comparison_code,headers] = xlsread([ basedir 'diagonal_board.xls']); %filename of the comparison_table.xlsx file. You can name it whatever you want. 
keylist = allData.keys()

stretchDataset1 = 0; %do not use this mode
stretchDataset2 = 0; %REALLY do not use this mode. 

deaths = 0; %are we synchronized to first or last timepoint in the file? 

windowSmoothing = 0; %run moving average
smoothingWindow = 3; %size of smoothing window in time points
mirroredWindow = 1; %do you want to mirror at the first and last time points? 

abstrue = 1; %absolute value of the area? (yes)

delta_function_islist = 1; %which delta function to use (list =1, use this mode)

% dataset1LenRealTime = 97;
% dataset2LenRealTime = size(dataset2,1);

num_reps = 1000000; %number of simulations 
maxLenRealTime = 179; %clip length
maxLen = 179; %clip length, equal unless stretching. Which you shouldn't be.


%% Script begins. 
results_matrix = zeros(size(headers))-1;
value_matrix = {};
figure;


if windowSmoothing
    windowSize = smoothingWindow;
    beeWindow = (1/windowSize)*ones(1,windowSize);
    ehWindow = 1;
end
total_comparisons = sum(sum(comparison_code));
current_comparison = 0;

for aye = 1:size(headers,1)-1
    for bee = 1:size(headers,2)-1
        if comparison_code(aye,bee)
            
            disp(['current comparison: ' num2str(current_comparison+1) ' of ' num2str(total_comparisons)])
            disp(['roughly ' num2str(current_comparison/total_comparisons*100) '% complete'])
            current_comparison = current_comparison+1;
            
            dataset1code = headers(aye+1);
            dataset1code = dataset1code{1};
            dataset2code = headers(bee+1);
            dataset2code = dataset2code{1};
            disp(['comparing ' dataset1code ' and ' dataset2code '.'])
            
            dataset1 = allData(dataset1code);
            dataset2 = allData(dataset2code);
            if deaths
                dataset1 = flipud(dataset1);
                dataset2 = flipud(dataset2);
            end
            
            %get the true means
            
            
            % if stretchDataset1 || stretchDataset2
            %     error('You can only stretch one dataset at a time-- set the other stretch question to 0');
            % end
            % if stretchDataset1
            % %implement this later.
            % end
            
            
            
%             if stretchDataset2
%                 %
%                 %% remove negative numbers and NaNs.
%                 
%                 %first dataset
%                 dataset1(dataset1<=0) = nan;
%                 xdata=(1:length(dataset1))';
%                 dataset1prime = zeros(size(dataset1,1),1);
%                 county = 1;
%                 for aye = 1:size(dataset1,2)
%                     temp = dataset1(:,aye);
%                     if sum(~isnan(temp)) %this rather naughty thing says 'if there are non-nan elements'
%                         dataset1prime(:,county)=interp1(xdata(~isnan(temp)),temp(~isnan(temp)),xdata,'linear');
%                         county = county+1;
%                     end
%                 end
%                 dataset1 = dataset1prime;
%                 
%                 
%                 %second dataset
%                 dataset2(dataset2<=0) = nan;
%                 xdata=(1:length(dataset2))';
%                 dataset1prime = zeros(size(dataset2,1),1);
%                 county = 1;
%                 for aye = 1:size(dataset2,2)
%                     temp = dataset2(:,aye);
%                     if sum(~isnan(temp)) %this rather naughty thing says 'if there are non-nan elements'
%                         dataset1prime(:,county)=interp1(xdata(~isnan(temp)),temp(~isnan(temp)),xdata,'linear');
%                         county = county+1;
%                     end
%                 end
%                 dataset2 = dataset1prime;
%                 clear dataset1prime;
%                 
%                 %% stretch the unstretched data.
%                 maxLen = floor(maxLenRealTime/dataset1LenRealTime*10000); % replace 10000 with measured max from other trace.
%                 dataset2interp = zeros(maxLen,1);
%                 county=1;
%                 for aye = 1:size(dataset2,2)
%                     temp = dataset2(:,aye);
%                     if sum(~isnan(temp)) >4
%                         
%                         no_nan_prime = 1/4*(temp(1:end-3)+temp(2:end-2)+temp(3:end-1)+temp(4:end));
%                         
%                         cliplen = find(isnan(no_nan_prime),1)-1; %last nonNaN element.
%                         if cliplen>maxLenRealTime
%                             dataset2interp(:,county) = interp1(linspace(1,maxLen,maxLenRealTime)', no_nan_prime(1:maxLenRealTime), linspace(1,maxLen,maxLen),'spline');
%                             county = county+1;
%                         else
%                             ourlen = floor(cliplen/maxLenRealTime*maxLen); %determine the stretch length when there aren't the real-time-length number of values available.
%                             dataset2interp(1:ourlen,county) = interp1(linspace(1,maxLen,cliplen)', no_nan_prime(1:cliplen), linspace(1,ourlen,ourlen),'spline');
%                             county = county+1;
%                         end
%                     end
%                     dataset2 = dataset2interp;
%                 end
%                 
%                 %linspace(first number, last number, number of points to specify
%                 %in-between*) *in this case we want this to be the real-time length of
%                 %the original vector.
%             end
            
            %% Bootstrapping
            disp('Start Par')
            [truePair1 , truePair2] = makeMeans_JustData(dataset1,dataset2, maxLen);
            
            %Smoothing by 1d moving average.
            if windowSmoothing
                if mirroredWindow
                    %instantiate temporary variables.
                    temp1 = zeros(2*length(truePair1),1);
                    temp2 = zeros(2*length(truePair2),1);
                    
                    %fill temporary variables, right hand side
                    temp1(length(truePair1)+1:end) = truePair1;
                    temp2(length(truePair2)+1:end) = truePair2;
                    
                    %fill temporary variables left hand side,
                    temp1(1:length(truePair1)) = flipud(truePair1);
                    temp2(1:length(truePair2)) = flipud(truePair2);
                    
                    %filter on the temps
                    temp1 = filter(beeWindow, ehWindow, temp1);
                    temp2 = filter(beeWindow, ehWindow, temp2);
                    
                    %pass back filtered data.
                    truePair1 = temp1(length(truePair1)+1:end);
                    truePair2 = temp2(length(truePair2)+1:end);
                    
                else
                    %simply moving average filter mean as passed to system.
                    truePair1 = filter(beeWindow, ehWindow, truePair1);
                    truePair2 = filter(beeWindow, ehWindow, truePair2);
                end
                
                
            end
            
            
            clf;
            ax = gca;
            ax.YLim = [0.1,1.2];
            plot(truePair1,'c', 'LineWidth',2)
            hold on;
            plot(truePair2,'g','LineWidth',2)
            hold off;
            if delta_function_islist
                value = deltaFunctionLists(truePair1, truePair2);
            elseif abstrue
                value = abs(deltaFunction(truePair1, truePair2));
            else
                value = deltaFunction(truePair1, truePair2);
            end
            
            saveas(gcf,[savedir dataset1code '_' dataset2code '_true.png'])
            %bootstrap
            full_set = clippy_cat(dataset1,dataset2,maxLen);
            catchbucket = zeros(num_reps,1);
            len = size(dataset1,2);
            %% parfor!
            parfor aye = 1:num_reps
                [fake1, fake2] = randomizeTraces_new(full_set, len, replace);
                [mean1, mean2] = fakeMeans(fake1,fake2);
                
                
                %Smoothing by 1d moving average.
                if windowSmoothing
                    if mirroredWindow
                        %instantiate temporary variables.
                        temp1 = zeros(2*length(mean1),1);
                        temp2 = zeros(2*length(mean2),1);
                        
                        %fill temporary variables, right hand side
                        temp1(length(mean1)+1:end) = mean1;
                        temp2(length(mean2)+1:end) = mean2;
                        
                        %fill temporary variables left hand side,
                        temp1(1:length(mean1)) = flipud(mean1);
                        temp2(1:length(mean2)) = flipud(mean2);
                        
                        %filter on the temps
                        temp1 = filter(beeWindow, ehWindow, temp1);
                        temp2 = filter(beeWindow, ehWindow, temp2);
                        
                        %pass back filtered data.
                        mean1 = temp1(length(mean1)+1:end);
                        mean2 = temp2(length(mean2)+1:end);
                        
                    else
                        mean1 = filter(beeWindow, ehWindow, mean1);
                        mean2 = filter(beeWindow, ehWindow, mean2);
                        
                    end
                end
                %                     if mod(aye,1000) == 0
                %                         figure;
                %                         hold on;
                %                         ax = gca;
                %                         ax.YLim = [0,1];
                %                         plot(mean1,'b', 'LineWidth',2)
                %                         plot(mean2,'m','LineWidth',2)
                %                         hold off;
                %                         saveas(gcf,['C:\Users\data\Documents\dump\' num2str(aye) '.png'])
                %                     end
                if delta_function_islist
                    catchbucket(aye) = deltaFunctionLists(mean1,mean2);
                elseif abstrue
                    catchbucket(aye) = abs(deltaFunction(mean1,mean2));
                else
                    catchbucket(aye) = deltaFunction(mean1,mean2);
                end
                
            end
            
            if value < 0
                county = catchbucket>value;
            else
                county = catchbucket<value;
            end
            
            countysum = sum(county);
            format long
            pval = 1-(countysum/num_reps)
            
            hold off;
            histogram(catchbucket)
            hold on;
            xline(value,'color', 'r', 'LineWidth', 2);
            axxes = gca;
            xmin = axxes.XLim;
            ymax = axxes.YLim;
            text(xmin(1),ymax(2)-100,['P = ' num2str(pval)])
            saveas(gcf,  [savedir dataset1code '_' dataset2code '_noreplace_histo.png'])
            results_matrix(aye,bee) = pval;
            value_matrix{aye,bee} = catchbucket;
            toc
        end
    end
end
resultsTable = array2table(string(horzcat({'XX'},keylist))');
for aye = 1:length(keylist)
    resultsTable.(char(keylist(aye))) = results_matrix(:,aye); 
end
writetable( resultsTable,[savedir 'human_readable_results.xlsx'])
csvwrite([savedir 'results.csv'],results_matrix);
save([savedir 'distro_mat.mat'], 'value_matrix');


toc
