close all;
clear all; 
%read the stuff in
tic
basedir = 'C:\Users\data\Desktop\20200901_Stats_Package'; %path to the folder containing both the data.xlsx and comparison_table.xlsx files
filename = 'demo_points_data'; %basename for your data.xlsx file, also used in output naming. 
[allDataMatrix,allTextMatrix] = xlsread([basedir filename '.xlsx']); % data xls. 
[comparison_code,headers] = xlsread([ basedir 'demo_points_comparison_table.xlsx']); %comparison table .xls

mode = 'Mean'; %options implemented: 'Mean' 'Median' any other single parameter can be easily implemented by adding the option in deltaFunction_points. 
absmode = 1; %take the absolute value of the difference in parameter? 1=true, 0=5
num_reps = 100000000;

allData=containers.Map; 
allTextMatrix= allTextMatrix(1,:);
for aye = 1:length(allTextMatrix)
    currentCode = allTextMatrix(aye);
    allData(currentCode{1}) = allDataMatrix(:,aye);
end

%prevent overwriting
time_stamp = datestr(datetime('now'));
time_stamp = strrep(time_stamp,':','.'); 
savedir = [basedir time_stamp '\'];
mkdir(savedir);

keylist = allData.keys();

results_matrix = zeros(size(headers))-1;
value_matrix = {};

total_comparisons = sum(sum(comparison_code));
current_comparison = 0;

for aye = 1:size(headers,1)-1
    for bee = 1:size(headers,2)-1
        tic
        if comparison_code(aye,bee)
            disp(['current comparison: ' num2str(current_comparison+1) ' of ' num2str(total_comparisons)])
            disp(['roughly ' num2str(current_comparison/total_comparisons*100) '% complete'])
            current_comparison = current_comparison+1;
            
            dataset1 = headers(aye+1);
            dataset1 = dataset1{1};
            dataset2 = headers(bee+1);
            dataset2 = dataset2{1};
            disp(['comparing ' dataset1 ' and ' dataset2 '.'])
            
            
            
            %maxLen = 30;
            %get the true means

            truePair1 = allData(dataset1);
            truePair2 = allData(dataset2);
            
            %nan censorship
            truePair1 = rmmissing(truePair1);
            truePair2 = rmmissing(truePair2);
            if isempty(truePair1) || isempty(truePair2)
                disp(['Skipping ' dataset1 ' X ' dataset2 ' Because one of the inputs is empty.'])
            else
            % figure;
            % hold on;
            % ax = gca;
            % ax.YLim = [0,2];
            % plot(truePair1,'c', 'LineWidth',2)
            % plot(truePair2,'g','LineWidth',2)
            % hold off;
            
            if absmode
                value = abs(deltaFunction_points(truePair1', truePair2',mode));
            else
                value = deltaFunction_points(truePair1', truePair2',mode);
            end

            % saveas(gcf,['C:\Users\data\Documents\dump\true.png'])
            %bootstrap

            catchbucket = zeros(num_reps,1);

            parfor jaye = 1:num_reps
                [fake1, fake2] = randomizePoints(truePair1', truePair2');
                %     if mod(aye,1000) == 0
                %         figure;
                %         hold on;
                %         ax = gca;
                %         ax.YLim = [0,1];
                %         plot(mean1,'b', 'LineWidth',2)
                %         plot(mean2,'m','LineWidth',2)
                %         hold off;
                %         saveas(gcf,['C:\Users\data\Documents\dump\' num2str(aye) '.png'])
                %     end
                if absmode
                    catchbucket(jaye) = abs(deltaFunction_points(fake1,fake2,mode));
                else
                    catchbucket(jaye) = deltaFunction_points(fake1,fake2,mode);
                end
            end

            histogram(catchbucket,'Normalization','probability')
            

            county = catchbucket >= value;
            countysum = sum(county);
            disp(['Comparison between: ' dataset1 ' and ' dataset2 ':'])
            format long
            pval = (countysum/num_reps)
            hold on;
            ax = gca;
            xline(value,'color', 'r', 'LineWidth', 2); 
            axxes = gca;
            xmin = axxes.XLim;
            ymax = axxes.YLim;
            text(xmin(1),ymax(2)-100,['P = ' num2str(pval)])
            hold off;
            saveas(gcf,  [savedir filename '_histo_' dataset1 '_' dataset2 '.png'])
            toc
            results_matrix(aye,bee) = pval;
            value_matrix{aye,bee} = catchbucket;
            end
            toc
        end
    end
end

%make metadata mat

meta_table = containers.Map();
meta_table('basedir') = basedir;
meta_table('filename') = filename;
meta_table('headers') = headers;
meta_table('num_reps') = num_reps;
meta_table('mode') = mode;
meta_table('absmode') = absmode;



%saving
csvwrite([savedir filename '_results.csv'],results_matrix);
save([savedir filename '_distro_mat.mat'], 'value_matrix');
save([savedir filename '_meta_mat.mat'], 'meta_table');

