%% Pull Excel Data into Matlab
% Will take at least 500 seconds to run on a good computer
tic %grabs start time
xlfile = '5_REMOVETRIALS_Old_New_Experiment_NoMemory_2Halves_ViewThenMemory_Strong_15sec_PROCESSED_SETH.xlsx';
[status,sheets] = xlsfinfo(xlfile);

subjectdata = {}; %will hold reorganized subject data
for subject = 3:length(sheets) %1st and 2nd sheets are data summary
    
    [~,~,raw] = xlsread(xlfile,subject,'A1:AC2500');
    studyimagescol = raw(2:end,2); %column B
    testimagescol = raw(2:end,16); %column P
    
    studyimages = NaN(1,length(studyimagescol)); %image #'s, allocate space to speed up
    testimages = NaN(1,length(testimagescol)); %image #'s, allocate space to speed up
    for i = 1:length(studyimagescol);
        if ~isnan(studyimagescol{i}) %if the excel cell is not blank
            if ~strcmpi('cross.bmp',studyimagescol{i}) %if the excel cell does not contain a cross
                studyimages(i) = str2num(studyimagescol{i}(1:end-3));
            end
        end
    end
    for i = 1:length(testimagescol);
        if ~isnan(testimagescol{i}) %if the excel cell is not blank
            if ~strcmpi('cross.bmp',testimagescol{i}) %if the excel cell does not contain a cross
                testimages(i) = str2num(testimagescol{i}(1:end-3));
            end
        end
    end
    studyimages(isnan(studyimages)) = []; %remove NaN's
    testimages(isnan(testimages)) = []; %remove NaN's
    
    uniquestudyimages = unique(studyimages);
    uniquetestimages = unique(testimages);
    
    study_trials = NaN(1,length(uniquestudyimages));
    study_rowstart_end = NaN(2,length(uniquestudyimages));
    study_xy = cell(1,length(uniquestudyimages));
    for i = 1:length(uniquestudyimages);
        rows = find(studyimages == uniquestudyimages(i)); %find which cells in excel sheet are for that image
        trial = raw{rows(1)+1,1}; %determine the trial for that image
        study_trials(i) = trial;
        study_rowstart_end(:,i) = [rows(1) rows(end)]; %save which rows in excel sheet
        xs = NaN(1,length(rows));
        ys = NaN(1,length(rows));
        for r = 1:length(rows);
            xs(r) = raw{rows(r)+1,11};
            ys(r) = raw{rows(r)+1,12};
        end
        xs(isnan(xs)) = []; %if adjusted fixation positions are missing
        ys(isnan(ys)) = []; %if adjusted fixation positions are missing
        study_xy{i} = [xs;ys];
    end
    
    test_xy = cell(1,length(uniquetestimages));
    test_trials = NaN(1,length(uniquetestimages));
    test_rowstart_end = NaN(2,length(uniquetestimages));
    for i = 1:length(uniquetestimages); %find which cells in excel sheet are for that image
        rows = find(testimages == uniquetestimages(i));
        d = find(diff(rows)>1); %find when there is a break in row index
        if ~isempty(d)
            rows = rows(1:d); %only want the first presentation for the test phase
        end
        trial = raw{rows(1)+1,15}; %determine the trial for that image
        test_trials(i) = trial;
        test_rowstart_end(:,i) = [rows(1) rows(end)]; %save which rows in excel sheet
        xs = NaN(1,length(rows));
        ys = NaN(1,length(rows));
        for r = 1:length(rows);
            xs(r) = raw{rows(r)+1,25};
            ys(r) = raw{rows(r)+1,26};
        end
        xs(isnan(xs)) = []; %if adjusted fixation positions are missing
        ys(isnan(ys)) = []; %if adjusted fixation positions are missing
        test_xy{i} = [xs;ys];
    end
    
    subjectdata{subject-2}.subject = sheets{subject}; %subject name based on excel sheet name
    
    [~,studyorder] = sort(study_trials); %sort everything by trial order
    %data arranged so column represents same trial/image across variables
    subjectdata{subject-2}.study_images = uniquestudyimages(studyorder); %images
    subjectdata{subject-2}.study_trials = study_trials(studyorder); %trial numbers
    subjectdata{subject-2}.study_rowstart_end = study_rowstart_end(:,studyorder); %rows of data for images in excel spread sheet
    subjectdata{subject-2}.study_xy = study_xy(studyorder); % fixation locations in adjusted x (row 1) and y(row 2) coordinates
    
    [~,testorder] = sort(test_trials); %sort everything by trial order
    %data arranged so column represents same trial/image across variables
    subjectdata{subject-2}.test_images = uniquetestimages(testorder); %images
    subjectdata{subject-2}.test_trials = test_trials(testorder); %trial numbers
    subjectdata{subject-2}.test_rowstart_end = test_rowstart_end(:,testorder); %rows of data for images in excel spread sheet
    subjectdata{subject-2}.test_xy = test_xy(testorder); % fixation locations in adjusted x (row 1) and y(row 2) coordinates
    
end
save('SmithSquireData_35','subjectdata') %save data since takes so long to ru
clearvars -except subjectdata
toc %tells you how long code took to run
%% Calculate KL Divergence
% run after you have saved SquireData then won't have to run the above section again
% roughly 5,000+ seconds to run on a good computer
tic
load('SmithSquireData.mat');

f = fspecial('gaussian',[256,256],24); %~2 dva 2D gaussian filter
binsize = 32; %some integer divisor of horizontal and vertical size of image
imageX = 1024; %horizontal pixel size of images
imageY = 768; %vertical pixel size of images

overlap = cell(1,length(subjectdata)); % Overlap measured with KL Divergence
Chanceoverlap = cell(1,length(subjectdata)); % Expected overlap by Chance
trialdifference = cell(1,length(subjectdata)); %Difference in image presentation order
imagepairs = cell(1,length(subjectdata));
for subject = 1:length(subjectdata)
    
    studyimages = subjectdata{subject}.study_images;
    testimages = subjectdata{subject}.test_images;
    max_study_trial = length(subjectdata{subject}.study_trials); %The start the order of image presention during test phase
    
    imagepair = []; %row 1 study, row 2 test
    studyorphans = []; %pictures shown during study phase but not test phase
    for img = 1:length(studyimages);
        imgind = find(studyimages(img) == testimages);
        if ~isempty(imgind)
            imagepair = [imagepair [studyimages(img); testimages(imgind)]];
        else
            studyorphans = [studyorphans  studyimages(img)];
        end
    end
    [~,ia,~] = intersect(testimages,imagepair(2,:));
    testorphans = testimages(~ia); %pictures shown during test phase but not study phase
    
    for pairs = 1:size(imagepair,2);
        
        %pre-allocate space for fixations pdfs (maps)
        num_fix_groups = 8; %this is the number of groups of fixations, 1 is for all fixations
        novelmap = cell(1,num_fix_groups); %Fixation matrix for a study phase image
        familiarmap = cell(1,num_fix_groups); %Fixation matrix for the same image but during test phase image
        for i = 1:num_fix_groups
            novelmap{i} = zeros(imageY,imageX);
            familiarmap{i} = zeros(imageY,imageX);
        end
        
        novelfixations = subjectdata{subject}.study_xy{studyimages == imagepair(1,pairs)};
        novelfixations(1,:) = round(novelfixations(1,:)*1024); %convert normalized x to pixel coordinate
        novelfixations(2,:) = round(novelfixations(2,:)*768); %convert normalized y to pixel coordinate
        novelfixations = novelfixations(:,2:end); %remove first fixation since this is on the cross hair
        
        familiarfixations = subjectdata{subject}.test_xy{testimages == imagepair(2,pairs)};
        familiarfixations(1,:) = round(familiarfixations(1,:)*1024); %convert normalized x to pixel coordinate
        familiarfixations(2,:) = round(familiarfixations(2,:)*768); %convert normalized y to pixel coordinate
        familiarfixations = familiarfixations(:,2:end); %remove first fixation since this is on the cross hair
        
        if size(familiarfixations,2) >= 5 && size(novelfixations,2) >= 5
            maxfixes = min(size(familiarfixations,2),size(novelfixations,2)); %so same number of fixations comparison
            
            %break fixations into groups of 5 or so we can get order and location information
            for fix = 1:maxfixes
                fixx = novelfixations(1,fix);
                fixx(fixx < 1) = 1; fixx(fixx > imageX) = imageX;
                fixy = imageY - novelfixations(2,fix);
                fixy(fixy < 1) = 1; fixy(fixy > imageY) = imageY;
                %cateogrizes fixations into groups of 5
                if fix <=  5
                    novelmap{1}(fixy,fixx) = novelmap{1}(fixy,fixx) + 1;
                elseif fix <= 10 && maxfixes >= 10
                    novelmap{2}(fixy,fixx) = novelmap{2}(fixy,fixx) + 1;
                elseif fix <= 15 && maxfixes >= 15
                    novelmap{3}(fixy,fixx) = novelmap{3}(fixy,fixx) + 1;
                elseif fix <= 20 && maxfixes >= 20
                    novelmap{4}(fixy,fixx) = novelmap{4}(fixy,fixx) + 1;
                elseif fix <= 25 && maxfixes >= 25
                    novelmap{5}(fixy,fixx) = novelmap{5}(fixy,fixx) + 1;
                elseif fix <= 30 && maxfixes >= 30
                    novelmap{6}(fixy,fixx) = novelmap{6}(fixy,fixx) + 1;
                elseif fix <= 35 && maxfixes >= 35
                    novelmap{7}(fixy,fixx) = novelmap{7}(fixy,fixx) + 1;
                end
                novelmap{num_fix_groups}(fixy,fixx) = novelmap{num_fix_groups}(fixy,fixx) + 1; % all fixations
            end
            
            for fix = 1:maxfixes;
                fixx = familiarfixations(1,fix);
                fixx(fixx < 1) = 1; fixx(fixx > imageX) = imageX;
                fixy = imageY - familiarfixations(2,fix);
                fixy(fixy < 1) = 1; fixy(fixy > imageY) = imageY;
                %cateogrizes fixations into groups of 5
                if fix <=  5
                    familiarmap{1}(fixy,fixx) = familiarmap{1}(fixy,fixx) + 1;
                elseif fix <= 10 && maxfixes >= 10
                    familiarmap{2}(fixy,fixx) = familiarmap{2}(fixy,fixx) + 1;
                elseif fix <= 15 && maxfixes >= 15
                    familiarmap{3}(fixy,fixx) = familiarmap{3}(fixy,fixx) + 1;
                elseif fix <= 20 && maxfixes >= 20
                    familiarmap{4}(fixy,fixx) = familiarmap{4}(fixy,fixx) + 1;
                elseif fix <= 25 && maxfixes >= 25
                    familiarmap{5}(fixy,fixx) = familiarmap{5}(fixy,fixx) + 1;
                elseif fix <= 30 && maxfixes >= 30
                    familiarmap{6}(fixy,fixx) = familiarmap{6}(fixy,fixx) + 1;
                elseif fix <= 35 && maxfixes >= 35
                    familiarmap{7}(fixy,fixx) = familiarmap{7}(fixy,fixx) + 1;
                end
                familiarmap{num_fix_groups}(fixy,fixx) = familiarmap{num_fix_groups}(fixy,fixx) + 1; % all fixations
            end
            
            % For loop below...
            % 1. convoles fixation matrix with 24 pixel gaussian filter
            % to account for variability in fixation locations and eye tracking error
            % 2. Bins (the verb) convolved matrix
            % 3. removes 0's and replaces with minimum defined value in matlab eps (2 ^-52)
            % 4. creates fixation pdf by dividing matrix by the sum of the matrix
            for i = 1:size(novelmap,2);
                novelmap{i} = imfilter(novelmap{i},f);
                novelmap{i} = bin2(novelmap{i},binsize,binsize);
                novelmap{i}(novelmap{i} == 0) = eps;
                novelmap{i} = novelmap{i}./sum(sum(novelmap{i}));
                
                familiarmap{i} = imfilter(familiarmap{i},f);
                familiarmap{i} = bin2(familiarmap{i},binsize,binsize);
                familiarmap{i}(familiarmap{i} == 0) = eps;
                familiarmap{i} = familiarmap{i}./sum(sum(familiarmap{i}));
            end
            
            %calculate overlap using KL divergence. Close/more similar PDFs have smaller KL divergence (distance in bits)
            for i = 1:size(novelmap,2)
                if all(all(novelmap{i} == mean(mean(novelmap{i})))); %if no fixations aka flat PDF
                    overlap{subject}(pairs,i) = NaN;
                else
                    overlap{subject}(pairs,i) = ...
                        sum(sum(log2(novelmap{i}./familiarmap{i}).*novelmap{i})) + ...
                        sum(sum(log2(familiarmap{i}./novelmap{i}).*familiarmap{i})); %KL divergence
                end
            end
            trialdifference{subject}(pairs) = (find(testimages == imagepair(2,pairs))+max_study_trial)...
                - find(studyimages == imagepair(1,pairs)); %Difference in image presentation order
        else %if not enough fixatiosn keep structure by placing NaNs
            overlap{subject}(pairs,1:size(novelmap,2)) = NaN(1,size(novelmap,2));
            trialdifference{subject}(pairs) = NaN;
        end
    end
    imagepairs{subject} = imagepair;
    
    % calculate Chance levels for overlap using random pairings of study stage image
    % using subject as control if have differences say in central bias
    Chancepairs = tril(ones(length(studyimages))-eye(length(studyimages))); %all possible pairs
    [Chance_row, Chance_col] = find(Chancepairs);
    randomization = randperm(length(Chance_row));
    Chance_row = Chance_row(randomization); %randomize order of pairs
    Chance_col = Chance_col(randomization); %randomize order of pairs
    total = 0; %total of successfully calculated Chance overlap values
    count = 0; %index since some pairs might not have enough fixations
    while total < 100 %want 100 pairs to calculate Chance levels
        count = count + 1;
        
        %all data is from study phase but going to call 1st one novel and 2nd familiar
        %this will allow us to use the same code as before
        novelfixations = subjectdata{subject}.study_xy{studyimages == studyimages(Chance_row(count))};
        novelfixations(1,:) = round(novelfixations(1,:)*1024);
        novelfixations(2,:) = round(novelfixations(2,:)*768);
        novelfixations = novelfixations(:,2:end); %remove first fixation since this is on the cross hair
        familiarfixations = subjectdata{subject}.study_xy{studyimages == studyimages(Chance_col(count))};
        familiarfixations(1,:) = round(familiarfixations(1,:)*1024);
        familiarfixations(2,:) = round(familiarfixations(2,:)*768);
        familiarfixations = familiarfixations(:,2:end); %remove first fixation since this is on the cross hair
        
        if size(familiarfixations,2) >= 5 && size(novelfixations,2) >= 5
            %pre-allocate space for fixations pdfs (maps)
            maxfixes = min(size(familiarfixations,2),size(novelfixations,2));
            novelmap = cell(1,num_fix_groups); %Fixation matrix for a study phase image
            familiarmap = cell(1,num_fix_groups); %Fixation matrix for the same image but during test phase image
            for i = 1:num_fix_groups
                novelmap{i} = zeros(imageY,imageX);
                familiarmap{i} = zeros(imageY,imageX);
            end
            
            %break fixations into groups of 5 or so we can get order and location information
            for fix = 1:maxfixes;
                fixx = novelfixations(1,fix);
                fixx(fixx < 1) = 1; fixx(fixx > imageX) = imageX;
                fixy = imageY - novelfixations(2,fix);
                fixy(fixy < 1) = 1; fixy(fixy > imageY) = imageY;
                %cateogrizes fixations into groups of 5
                if fix <=  5
                    novelmap{1}(fixy,fixx) = novelmap{1}(fixy,fixx) + 1;
                elseif fix <= 10 && maxfixes >= 10
                    novelmap{2}(fixy,fixx) = novelmap{2}(fixy,fixx) + 1;
                elseif fix <= 15 && maxfixes >= 15
                    novelmap{3}(fixy,fixx) = novelmap{3}(fixy,fixx) + 1;
                elseif fix <= 20 && maxfixes >= 20
                    novelmap{4}(fixy,fixx) = novelmap{4}(fixy,fixx) + 1;
                elseif fix <= 25 && maxfixes >= 25
                    novelmap{5}(fixy,fixx) = novelmap{5}(fixy,fixx) + 1;
                elseif fix <= 30 && maxfixes >= 30
                    novelmap{6}(fixy,fixx) = novelmap{6}(fixy,fixx) + 1;
                elseif fix <= 35 && maxfixes >= 35
                    novelmap{7}(fixy,fixx) = novelmap{7}(fixy,fixx) + 1;
                end
                novelmap{num_fix_groups}(fixy,fixx) = novelmap{num_fix_groups}(fixy,fixx) + 1; % all fixations
            end
            
            for fix = 1:maxfixes;
                fixx = familiarfixations(1,fix);
                fixx(fixx < 1) = 1; fixx(fixx > imageX) = imageX;
                fixy = imageY - familiarfixations(2,fix);
                fixy(fixy < 1) = 1; fixy(fixy > imageY) = imageY;
                %cateogrizes fixations into groups of 5
                if fix <=  5
                    familiarmap{1}(fixy,fixx) = familiarmap{1}(fixy,fixx) + 1;
                elseif fix <= 10 && maxfixes >= 10
                    familiarmap{2}(fixy,fixx) = familiarmap{2}(fixy,fixx) + 1;
                elseif fix <= 15 && maxfixes >= 15
                    familiarmap{3}(fixy,fixx) = familiarmap{3}(fixy,fixx) + 1;
                elseif fix <= 20 && maxfixes >= 20
                    familiarmap{4}(fixy,fixx) = familiarmap{4}(fixy,fixx) + 1;
                elseif fix <= 25 && maxfixes >= 25
                    familiarmap{5}(fixy,fixx) = familiarmap{5}(fixy,fixx) + 1;
                elseif fix <= 30 && maxfixes >= 30
                    familiarmap{6}(fixy,fixx) = familiarmap{6}(fixy,fixx) + 1;
                elseif fix <= 35 && maxfixes >= 35
                    familiarmap{7}(fixy,fixx) = familiarmap{7}(fixy,fixx) + 1;
                end
                familiarmap{num_fix_groups}(fixy,fixx) = familiarmap{num_fix_groups}(fixy,fixx) + 1; % all fixations
            end
            
            % For loop below...
            % 1. convoles fixation matrix with 1 dva (24 pixel) gaussian filter
            % to account for variability in fixation locations and eye tracking error
            % 2. Bins convolved matrix
            % 3. removes 0's and replaces with minimum defined value in matlab eps (2 ^-52)
            % 4. creates fixation pdf by dividing matrix by the sum of the matrix
            for i = 1:size(novelmap,2);
                novelmap{i} = imfilter(novelmap{i},f);
                novelmap{i} = bin2(novelmap{i},binsize,binsize);
                novelmap{i}(novelmap{i} == 0) = eps;
                novelmap{i} = novelmap{i}./sum(sum(novelmap{i}));
                
                familiarmap{i} = imfilter(familiarmap{i},f);
                familiarmap{i} = bin2(familiarmap{i},binsize,binsize);
                familiarmap{i}(familiarmap{i} == 0) = eps;
                familiarmap{i} = familiarmap{i}./sum(sum(familiarmap{i}));
            end
            
            %calculate overlap using KL divergence. Close more similar PDFs have smaller KL divergence (distance in bits)
            for i = 1:size(novelmap,2)
                if all(all(novelmap{i} == mean(mean(novelmap{i})))); %if no fixations aka flat PDF
                    Chanceoverlap{subject}(total+1,i) = NaN;
                else
                    Chanceoverlap{subject}(total+1,i) = ...
                        sum(sum(log2(novelmap{i}./familiarmap{i}).*novelmap{i})) + ...
                        sum(sum(log2(familiarmap{i}./novelmap{i}).*familiarmap{i})); % KL divergence
                end
            end
            total = total + 1; %successfully got a piar of images to use for data
        end
    end
end
save('KLDivergenceData_25','overlap','Chanceoverlap','trialdifference','imagepairs')
%save raw data since takes an extremely long time to run
%overlap/Chance overlap column format: fixations 1-5, fixations 6-10,
%11-15,... all fixations 1-35
toc
%% Plot Number of images in between presentations vs overlap in fixations
% probably not very helpful
for subject = 1:length(subjectdata);
    figure
    hold on
    for i = 1:size(overlap{subject},1);
        plot(trialdifference{subject}(i),overlap{subject}(i,3),'*')
    end
    xlabel('Number of images in between novel and repeated presentation')
    ylabel('Overlap/KL Divergence (Distance in bits)')
end

%% Plot Median Chance Overlap (blue) & Median overlap (red) by Subject
% probably not very helpful
figure
hold on
for subject = 1:length(subjectdata);
    plot(subject,median(Chanceoverlap{subject}(:,3)),'*')
end
for subject = 1:length(subjectdata);
    plot(subject,nanmedian(overlap{subject}(:,3)),'r*')
end
xlabel('Subject "number"')
ylabel('Overlap/KL Divergence (Distance in bits)')
%% Convert Data back to Excel
filename = 'Subject_and_KL_Data_15sec.xls';
load('SmithSquireData_35.mat')
load('KLDivergenceData_35.mat')

for subject = 1:length(subjectdata);
    sheet = subjectdata{subject}.subject; %sheet number
    
    %raw will hold all the data that's written to a sheet
    raw = {'Study Trial','Image','Adjusted X','Adjusted Y','Test Trial','Image',...
        'Adjusted X','Adjusted Y','Study Image','Test Image','KL Divergence 1-5',...
        'KL Divergence 6-10','KL Divergence 11-15','KL Divergence 16-20','KL Divergence 21-25',...
        'KL Divergence 26-30','KL Divergence 31-35','KL Divergence All','Chance 1-5',...
        'Chance 6-10','Chance 11-15','Chance 16-20','Chance 21-25',...
        'Chance 26-30','Chance 31-35','Chance All'}; %1st line, title line
    
    % write study subject data to raw cell array
    previousline = 1;
    for trial = 1:length(subjectdata{subject}.study_xy);
        for lines = 1:size(subjectdata{subject}.study_xy{trial},2);
            raw{previousline+lines,1} = subjectdata{subject}.study_trials(trial);
            raw{previousline+lines,2} = subjectdata{subject}.study_images(trial);
            raw{previousline+lines,3} = subjectdata{subject}.study_xy{trial}(1,lines);
            raw{previousline+lines,4} = subjectdata{subject}.study_xy{trial}(2,lines);
        end
        previousline = previousline+size(subjectdata{subject}.study_xy{trial},2);
    end
    
    % write test subject data to raw cell array
    previousline = 1;
    for trial = 1:length(subjectdata{subject}.test_xy);
        for lines = 1:size(subjectdata{subject}.test_xy{trial},2);
            raw{previousline+lines,5} = subjectdata{subject}.test_trials(trial);
            raw{previousline+lines,6} = subjectdata{subject}.test_images(trial);
            raw{previousline+lines,7} = subjectdata{subject}.test_xy{trial}(1,lines);
            raw{previousline+lines,8} = subjectdata{subject}.test_xy{trial}(2,lines);
        end
        previousline = previousline+size(subjectdata{subject}.test_xy{trial},2);
    end
    
    % write KL divergence data to raw cell array
    for trial = 1:size(overlap{subject},1);
        raw{trial+1,9} = imagepairs{subject}(1,trial); %study image
        raw{trial+1,10} = imagepairs{subject}(2,trial); %test image
        raw{trial+1,11} = overlap{subject}(trial,1); %KL divergence fixations 1-5
        raw{trial+1,12} = overlap{subject}(trial,2); %KL divergence fixations 6-10
        raw{trial+1,13} = overlap{subject}(trial,3); %KL divergence fixations 11-15
        raw{trial+1,14} = overlap{subject}(trial,4); %KL divergence fixations 16-20
        raw{trial+1,15} = overlap{subject}(trial,5); %KL divergence fixations 21-25
        raw{trial+1,16} = overlap{subject}(trial,6); %KL divergence fixations 26-30
        raw{trial+1,17} = overlap{subject}(trial,7); %KL divergence fixations 31-35
        raw{trial+1,18} = overlap{subject}(trial,8); %KL divergence all fixations;
    end
    
    % write Chance KL divergence data to raw cell array
    for trial = 1:size(chanceoverlap{subject});
        raw{trial+1,19} = chanceoverlap{subject}(trial,1); %chance KL divergence fixations 1-5
        raw{trial+1,20} = chanceoverlap{subject}(trial,2); %chance KL divergence fixations 6-10
        raw{trial+1,21} = chanceoverlap{subject}(trial,3); %chance KL divergence fixations 11-15
        raw{trial+1,22} = chanceoverlap{subject}(trial,4); %chance KL divergence fixations 16-20
        raw{trial+1,23} = chanceoverlap{subject}(trial,5); %chance KL divergence fixations 21-25
        raw{trial+1,24} = chanceoverlap{subject}(trial,6); %chance KL divergence fixations 26-30
        raw{trial+1,25} = chanceoverlap{subject}(trial,7); %chance KL divergence fixations 31-35
        raw{trial+1,26} = chanceoverlap{subject}(trial,8); %chance KL divergence all fixations;
    end
    xlswrite(filename,raw,sheet)
end