%% Pull Excel Data into Matlab
% Will take at least 500 seconds to run on a good computer
tic %grabs start time
xlfile = '0317.xlsx';
[status,sheets] = xlsfinfo(xlfile);

%Call directory for where the salience heat maps are located
%addpath('P:/Christine/Experiments/EyeTracking/Old_New_NoMemory/2Halves_Experiment/scenes for heatmaps')
%for k=1:80
    %matfilename=sprintf('%03d-saliencemap.mat',k);
    %saliencemaps{k}=load(matfilename);
%end;

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
    study_xy_pixel = cell(1,length(uniquetestimages));
    study_salience_values = cell(1,length(uniquetestimages));
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
        
        %x_pixel=round(xs*1024); %convert coordinate onto 1024 scale, column #
       % y_pixel=round(ys*768); %convert coordinate onto 768 scale, row #
        
       
        if  round(xs*1024)<1 %convert coordinate onto 1024 scale, column #. If this rounds to zero force it to be a 1.
            x_pixel = 1;
        elseif round(xs*1024)>0
            x_pixel=round(xs*1024);
        end
        
        if  round(ys*768)<1 %convert coordinate onto 768 scale, row #. If this rounds to zero force it to be a 1.
            y_pixel = 1;
        elseif round(ys*768)>0
            y_pixel=round(ys*768);
        end   
        
        study_xy_pixel{i} = [x_pixel;y_pixel];
        
        current_image=uniquestudyimages(i);
        addpath('C:\Users\seth.koenig\Documents\MATLAB\Smith and Squire\ALL-saliencemap\');
        matfilename=sprintf('%03d-saliencemap.mat',current_image);
        current_saliencemap_temp=load(matfilename);
        current_saliencemap=current_saliencemap_temp.saliencemap;
        for s=1:length(rows);
            salience{s}=[current_saliencemap((study_xy_pixel{1,i}(2,s)),(study_xy_pixel{1,i}(1,s)))];
            study_salience_values{s,i}=[salience{s}];
        end;
    end
    
    test_xy = cell(1,length(uniquetestimages));
    test_xy_pixel = cell(1,length(uniquetestimages));
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
       % x_pixel=round(xs*1024); %convert coordinate onto 1024 scale, column #. 
       % y_pixel=round(ys*768); %convert coordinate onto 768 scale, row #
        if  round(xs*1024)<1 %convert coordinate onto 1024 scale, column #. If this rounds to zero force it to be a 1.
            x_pixel = 1;
        elseif round(xs*1024)>0
            x_pixel=round(xs*1024);
        end
        
        if  round(ys*768)<1 %convert coordinate onto 768 scale, row #. If this rounds to zero force it to be a 1.
            y_pixel=1;
        elseif round(ys*768)>0
            y_pixel=round(ys*768);
        end
        
        
        test_xy_pixel{i} = [x_pixel;y_pixel];
        
        current_image=uniquetestimages(i);
        addpath('C:\Users\seth.koenig\Documents\MATLAB\Smith and Squire\ALL-saliencemap\');
        matfilename=sprintf('%03d-saliencemap.mat',current_image);
        current_saliencemap_temp=load(matfilename);
        current_saliencemap=current_saliencemap_temp.saliencemap;
        for s=1:length(rows);
            
            salience{s}=[current_saliencemap((test_xy_pixel{1,i}(2,s)),(test_xy_pixel{1,i}(1,s)))];
            test_salience_values{s,i}=[salience{s}];
        end;
    end
    
    subjectdata{subject-2}.subject = sheets{subject}; %subject name based on excel sheet name
    
    [~,studyorder] = sort(study_trials); %sort everything by trial order
    %data arranged so column represents same trial/image across variables
    subjectdata{subject-2}.study_images = uniquestudyimages(studyorder); %images
    subjectdata{subject-2}.study_trials = study_trials(studyorder); %trial numbers
    subjectdata{subject-2}.study_rowstart_end = study_rowstart_end(:,studyorder); %rows of data for images in excel spread sheet
    subjectdata{subject-2}.study_xy = study_xy(studyorder); % fixation locations in adjusted x (row 1) and y(row 2) coordinates
    subjectdata{subject-2}.study_xy_pixel = study_xy_pixel(studyorder); %converted x and y coordinates
    subjectdata{subject-2}.study_salience_values = study_salience_values(studyorder);
    
    [~,testorder] = sort(test_trials); %sort everything by trial order
    %data arranged so column represents same trial/image across variables
    subjectdata{subject-2}.test_images = uniquetestimages(testorder); %images
    subjectdata{subject-2}.test_trials = test_trials(testorder); %trial numbers
    subjectdata{subject-2}.test_rowstart_end = test_rowstart_end(:,testorder); %rows of data for images in excel spread sheet
    subjectdata{subject-2}.test_xy = test_xy(testorder); % fixation locations in adjusted x (row 1) and y(row 2) coordinates
    subjectdata{subject-2}.test_xy_pixel = test_xy_pixel(testorder);  %converted x and y coordinates
    subjectdata{subject-2}.test_salience_values = test_salience_values(testorder); 
end

save('SmithSquireData','subjectdata') %save data since takes so long to run
%clearvars -except subjectdata saliencemaps
toc %tells you how long code took to run

salience_study = NaN(1,length(subjectdata));
%AA = cell(1,length(subjectdata));
%salience_study = cell2mat(AA);
for s=1:length(subjectdata)
    %Cycles through all 40 images
    for im=1:length(subjectdata{1,s}.study_images)
        
        current_image=subjectdata{1,s}.study_images(im);
        addpath('C:\Users\seth.koenig\Documents\MATLAB\Smith and Squire\ALL-saliencemap\');
        matfilename=sprintf('%03d-saliencemap.mat',current_image);
        current_saliencemap_temp=load(matfilename);
        current_saliencemap=current_saliencemap_temp.saliencemap;
        
        for f=1:length(subjectdata{1,s}.study_xy_pixel{im})
            salience_study(im,f)=current_saliencemap(subjectdata{1,s}.study_xy_pixel{1,im}(2,f),subjectdata{1,s}.study_xy_pixel{1,im}(1,f));
        end;
    end;
    subjectdata{s}.salience_study=salience_study;
end;

salience_test = NaN(1,length(subjectdata));
%AB = cell(1,length(subjectdata));
%salience_test = cell2mat(AB);
for s=1:length(subjectdata)
    %Cycles through all 40 images
    for im=1:length(subjectdata{1,s}.test_images)
        
        current_image=subjectdata{1,s}.test_images(im);
        addpath('C:\Users\seth.koenig\Documents\MATLAB\Smith and Squire\ALL-saliencemap\');
        matfilename=sprintf('%03d-saliencemap.mat',current_image);
        current_saliencemap_temp=load(matfilename);
        current_saliencemap=current_saliencemap_temp.saliencemap;
        
        for f=1:length(subjectdata{1,s}.test_xy_pixel{im})
            salience_test(im,f)=current_saliencemap(subjectdata{1,s}.test_xy_pixel{1,im}(2,f),subjectdata{1,s}.test_xy_pixel{1,im}(1,f));
        end;
    end;
    subjectdata{s}.salience_test=salience_test;
end;

save('SalienceData','subjectdata')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%Creates output file
filename='0317_Salience_TestingGround.xls';
load('SalienceData.mat')
for subject=1:length(subjectdata)
   sheet=subjectdata{subject}.subject; %putting number of subject to sheet
    
   raw = {'Image','Fixation 1','Fixation 2','Fixation 3','Fixation 4','Fixation 5','Fixation 6','Fixation 7','Fixation 8','Fixation 9','Fixation 10',...
    'Fixation 11','Fixation 12','Fixation 13','Fixation 14','Fixation 15','Fixation 16','Fixation 17','Fixation 18','Fixation 19','Fixation 20',...
    'Fixation 21','Fixation 22','Fixation 23','Fixation 24','Fixation 25'};
    
    subjectdata{subject}.study_images = subjectdata{subject}.study_images';
    subjectdata{subject}.test_images = subjectdata{subject}.test_images';
 
    
       for image=1:40;
        raw{image+1,1}=subjectdata{subject}.study_images(image,1);
        raw{image+41,1}=subjectdata{subject}.test_images(image,1);
        
        end;
    
    subjectdata{subject}.study_images = subjectdata{subject}.study_images';
    subjectdata{subject}.test_images = subjectdata{subject}.test_images';
    
    
    salrows = size(subjectdata{subject}.salience_study,1); %number of rows in your data
    salcols = size(subjectdata{subject}.salience_study,2); %number of columns in your data
    raw(2:salrows+1,2:salcols+1) = num2cell(subjectdata{subject}.salience_study);
    
    salrows = size(subjectdata{subject}.salience_test,1); %number of rows in your data
    salcols = size(subjectdata{subject}.salience_test,2); %number of columns in your data
    raw(42:salrows+41,2:salcols+1) = num2cell(subjectdata{subject}.salience_test);
    
    %subjectdata{subject}.salience_study(1:length(subjectdata{subject}.study_images),1:length(subjectdata{subject}.salience_study));
    
    %for image=1:length(subjectdata{subject}.study_xy); 
       % raw{1+image,1}=subjectdata{subject}.study_images(image,1);
    %end;
    
    %numel(subjectdata{subject}.study_xy); 
    %numel(subjectdata{subject}.salience_study);
    %for fixation=1:length(subjectdata{subject}.salience_study);
        %raw{2:length(subjectdata{subject}.study_xy),2:length(subjectdata{subject}.salience_study)}=subjectdata{subject}.salience_study{1:length(study_xy),1:length(salience_study)};
    %end;
     
    %a = subjectdata{subject}.salience_study(1,:)
%    raw{2,2} = cell2mat(a)
    %nums = [raw{2,2}]

    
   
    
       %for image=1:40 
        %raw{image+1,1}=subjectdata{subject}.study_images(image,1);
        %raw{image+41,1}=subjectdata{subject}.test_images(image,1);
        %for fixation=1:end(salience_test)
            %raw{image+1,fixation+1}=subjectdata{subject}.salience_study{image,fixation};
            %raw{image+41,fixation+1}=subjectdata{subject}.salience_test{image,fixation};
        %end;
    %end;
    
    xlswrite(filename,raw,sheet)

end;

