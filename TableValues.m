% Define the deceptive groups for blocks 3 & 4 and blocks 6 & 7
deceptiveInBlocks34 = [1, 2, 3, 7, 21, 22, 23, 24, 27, 29];
deceptiveInBlocks67 = setdiff(1:29, deceptiveInBlocks34);
deceptiveInBlocks67(deceptiveInBlocks67 == 3) = [];  % Exclude subject number 3

% Dummy data for gender, replace with actual data loading if available
genderData = struct();
for k = 1:29
    genderData(k).Subject = k;
    genderData(k).Gender = 'Unknown';  % Default to unknown
end
genderData(1).Gender = 'Male';
genderData(2).Gender = 'Female';

function filteredData = applyLowPassFilter(data, sampleRate, cutoffFreq)
    % Butterworth Low-Pass Filter
    [b, a] = butter(3, cutoffFreq/(sampleRate/2), 'low');
    filteredData = filtfilt(b, a, data);
end

function filteredData = applyBandPassFilter(data, sampleRate, lowCutoff, highCutoff)
    % Butterworth Band-Pass Filter
    [b, a] = butter(3, [lowCutoff highCutoff]/(sampleRate/2), 'bandpass');
    filteredData = filtfilt(b, a, data);
end


% Folder selection and file loading
dataFolder = uigetdir('', 'Select the folder containing all subject files');
if dataFolder == 0
    error('No folder selected. Exiting script.');
end
fileList = dir(fullfile(dataFolder, '*.mat'));

% Initialize results structure
allResults = repmat(struct(), length(fileList), 1);

% Physiological signal calculation functions
function SCL = calculateSCL(EDA)
    SCL = mean(EDA);
end

function SCR = calculateSCR(EDA, sampleRate)
    % Apply a low-pass filter to reduce high-frequency noise
    filteredEDA = applyLowPassFilter(EDA, sampleRate, 1); % 1 Hz cutoff frequency
    [peaks, ~] = findpeaks(filteredEDA, 'MinPeakProminence', 0.01, 'MinPeakDistance', 50);
    SCR = length(peaks);
end

function HRV = calculateHRV(ECG, sampleRate)
    % Apply a band-pass filter to isolate the QRS complex
    filteredECG = applyBandPassFilter(ECG, sampleRate, 0.5, 40); % Bandpass between 0.5 Hz and 40 Hz
    [qrs_peaks, ~] = findpeaks(filteredECG, 'MinPeakHeight', max(filteredECG)*0.4, 'MinPeakDistance', 50);
    RR_intervals = diff(qrs_peaks);
    HRV = std(RR_intervals);
end

function gen = getGender(subNum, genderData)
    idx = find([genderData.Subject] == subNum);
    if ~isempty(idx)
        gen = genderData(idx).Gender;
    else
        gen = 'Unknown'; % Default case if no gender data is available
    end
end

% Main data loading and processing loop
sampleRate = 1000; % Define the sample rate of your data
for i = 1:length(fileList)
    fileName = fullfile(dataFolder, fileList(i).name);
    subjectNum = str2double(regexp(fileName, '\d+', 'match', 'once'));
    if ~isempty(subjectNum) && (ismember(subjectNum, deceptiveInBlocks34) || ismember(subjectNum, deceptiveInBlocks67))
        data = load(fileName);
        if isfield(data, 'bp_data')
            bp_data = data.bp_data;
            labels = cellstr(bp_data.labels);
            EDA_index = find(contains(lower(labels), 'eda'));
            ECG_index = find(contains(lower(labels), 'ecg'));
            if isempty(EDA_index) || isempty(ECG_index)
                continue;
            end
            gender = getGender(subjectNum, genderData);
            blockStarts = bp_data.BLOCK_START;
            blocksOfInterest = [3, 4, 6, 7];
            for block = blocksOfInterest
                if block <= length(blockStarts)
                    startIdx = round(blockStarts(block));
                    if block == length(blockStarts)
                        endIdx = size(bp_data.data, 1);
                    else
                        endIdx = round(blockStarts(block + 1)) - 1;
                    end

                    EDA_block = bp_data.data(startIdx:endIdx, EDA_index);
                    ECG_block = bp_data.data(startIdx:endIdx, ECG_index);
                    allResults(i).Gender = gender;

                    allResults(i).(strcat('SCL_Block', num2str(block))) = calculateSCL(EDA_block);
                    allResults(i).(strcat('SCR_Block', num2str(block))) = calculateSCR(EDA_block, sampleRate);
                    allResults(i).(strcat('HRV_Block', num2str(block))) = calculateHRV(ECG_block, sampleRate);
                end
            end
        end
    end
end

% Prepare data for ANOVA
numSubjects = length(fileList);
dataMatrix = NaN(numSubjects, 12);  % SCL, SCR, HRV for two blocks of interest per subject
genderVector = cell(numSubjects, 1);

for i = 1:numSubjects
    if ~isempty(allResults(i).Gender)
        genderVector{i} = allResults(i).Gender;
    end
    if isfield(allResults(i), 'SCL_Block3') && isfield(allResults(i), 'SCL_Block4') && isfield(allResults(i), 'SCL_Block6') && isfield(allResults(i), 'SCL_Block7')
        dataMatrix(i, 1:4) = [allResults(i).SCL_Block3, allResults(i).SCL_Block4, allResults(i).SCL_Block6, allResults(i).SCL_Block7];
    end
    if isfield(allResults(i), 'SCR_Block3') && isfield(allResults(i), 'SCR_Block4') && isfield(allResults(i), 'SCR_Block6') && isfield(allResults(i), 'SCR_Block7')
        dataMatrix(i, 5:8) = [allResults(i).SCR_Block3, allResults(i).SCR_Block4, allResults(i).SCR_Block6, allResults(i).SCR_Block7];
    end
    if isfield(allResults(i), 'HRV_Block3') && isfield(allResults(i), 'HRV_Block4') && isfield(allResults(i), 'HRV_Block6') && isfield(allResults(i), 'HRV_Block7')
        dataMatrix(i, 9:12) = [allResults(i).HRV_Block3, allResults(i).HRV_Block4, allResults(i).HRV_Block6, allResults(i).HRV_Block7];
    end
end

% Assume dataMatrix and genderVector are correctly populated
% Convert genderVector to categorical and clean data
genderVector = categorical(genderVector);
nanRows = any(isnan(dataMatrix), 2) | isundefined(genderVector);

% Remove rows with missing or incomplete data
dataMatrix(nanRows, :) = [];
genderVector(nanRows) = [];

% Check that all columns now have the same number of rows
if length(genderVector) == size(dataMatrix, 1)
    varNames = {'Gender', 'SCL3', 'SCL4', 'SCL6', 'SCL7', 'SCR3', 'SCR4', 'SCR6', 'SCR7', 'HRV3', 'HRV4', 'HRV6', 'HRV7'};
    rmTable = table(genderVector, dataMatrix(:,1), dataMatrix(:,2), dataMatrix(:,3), dataMatrix(:,4), ...
                    dataMatrix(:,5), dataMatrix(:,6), dataMatrix(:,7), dataMatrix(:,8), ...
                    dataMatrix(:,9), dataMatrix(:,10), dataMatrix(:,11), dataMatrix(:,12), ...
                    'VariableNames', varNames);
    
    TIME = [3, 4, 6, 7];  % Define blocks as time points for within-subject factor
    
    % Create separate repeated measures models for each physiological measure
    Rm_SCL = fitrm(rmTable, 'SCL3-SCL7 ~ Gender', 'WithinDesign', TIME);
    Rm_SCR = fitrm(rmTable, 'SCR3-SCR7 ~ Gender', 'WithinDesign', TIME);
    Rm_HRV = fitrm(rmTable, 'HRV3-HRV7 ~ Gender', 'WithinDesign', TIME);

    % Run and display ANOVA for each model
    ranovatbl_SCL = ranova(Rm_SCL);
    ranovatbl_SCR = ranova(Rm_SCR);
    ranovatbl_HRV = ranova(Rm_HRV);

    disp('ANOVA Results for SCL:');
    disp(ranovatbl_SCL);
    disp('ANOVA Results for SCR:');
    disp(ranovatbl_SCR);
    disp('ANOVA Results for HRV:');
    disp(ranovatbl_HRV);
else
    error('Mismatch in the number of rows between genderVector and dataMatrix');
end
