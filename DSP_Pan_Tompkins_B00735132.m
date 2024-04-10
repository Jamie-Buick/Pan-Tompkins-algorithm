%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                          %
%                           EEE826 Coursework                              %
%                                                                          %
%                         Jamie Buick - B00735132                          %
%                                                                          %
% This program contains a functioning Pan-Tompkins Algorithm which runs    %
% though the MIT-BIH database. The program returns the accuracy of the     %
% algorthm in comparison to the hand annotations placed on the ECG by a    %
% professional.                                                            %
%                                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Clear workspace, clear command window and close all open figures
close all;
clc;
clear;


%***************************************************************************%
%                                                                           %
%                           Initial Set-up                                  %
%                                                                           %
%***************************************************************************%

% All files in the directory which are .mat are saved into a 48x1 struct called
% files. The length of this struct is saved in the variable n.
files = dir('*.mat');
n = length(files);

for fullDataBase=1:1:n

% The patient files is loaded in by changing the number in files(x) from 1
% to 48 - this is done using for loop. The Fieldnames are retrieved and
% then tranferred into a char array. Each component of the patient .mat file is now loaded in.
patient = load(files(fullDataBase).name);
y = fieldnames(patient);
patientFileNumber = char(y);
patientFileNumberstr = convertCharsToStrings(patientFileNumber);
ecgSignal = patient.(patientFileNumber).sig;
ecgAnn = patient.(patientFileNumber).ann;
ecgAnnType = patient.(patientFileNumber).annType;

% This section of code removes the unwanted annotations from the ecgAnnType
% array. This is essential as the unwanted annotations (anything which is
% not a beat) will influence the final accuracy calculation.
remove = [];
t = 1;

for o = 1:length(ecgAnnType)
    if ecgAnnType(o) == "+"
        remove(t) = 0;
    elseif ecgAnnType(o) == "["
        remove(t) = 0;
    elseif ecgAnnType(o) == "]"
        remove(t) = 0;
    elseif ecgAnnType(o) == ")"
        remove(t) = 0;
    elseif ecgAnnType(o) == "("
        remove(t) = 0;
    elseif ecgAnnType(o) == "!"
        remove(t) = 0;
    elseif ecgAnnType(o) == "x"
        remove(t) = 0;
    elseif ecgAnnType(o) == "~"
        remove(t) = 0;
    elseif ecgAnnType(o) == "r"
        remove(t) = 0;
    elseif ecgAnnType(o) == "?"
        remove(t) = 0;
    elseif ecgAnnType(o) == "|"
        remove(t) = 0;
    elseif ecgAnnType(o) == '"'
        remove(t) = 0;
    elseif ecgAnnType(o) == "p"
        remove(t) = 0;
    elseif ecgAnnType(o) == "t"
        remove(t) = 0;
    elseif ecgAnnType(o) == "u"
        remove(t) = 0;
    elseif ecgAnnType(o) == "@"
        remove(t) = 0;
    elseif ecgAnnType(o) == "*"
        remove(t) = 0;
    elseif ecgAnnType(o) == "D"
        remove(t) = 0;
    elseif ecgAnnType(o) == "s"
        remove(t) = 0;
    elseif ecgAnnType(o) == "'"
        remove(t) = 0;
    elseif ecgAnnType(o) == "^"
        remove(t) = 0;
    elseif ecgAnnType(o) == "`"
        remove(t) = 0;
    elseif ecgAnnType(o) == "="
        remove(t) = 0;
    else
        remove(t) = 1;
    end
    t=t+1;
end

% FixedAnn is created by multipying the remove array with the original
% ecgAnn array. This sets the unwanted annotation to 0. The ann variable
% is then corrected by removing the 0 values.
fixedAnn = remove'.*ecgAnn;
ann = fixedAnn(fixedAnn~=0);

% Table of the original annotations side by side with the annotation types.
origAnnotation = table(ecgAnn,ecgAnnType);

% The sampling frequency is 360Hz. I have used this sampling frequency to
% convert into time domain for seconds and minutes.
fs = 360;
signalSamples = length(ecgSignal);
signalTimeSecs = (signalSamples/fs);
signalTimeMins = signalTimeSecs/60;



%***************************************************************************%
%                                                                           %
%                           Preprocessing                                   %
%                                                                           %
%***************************************************************************%

% A bandpass filter is used to remove noise between 5 and 23Hz. Butter is
% the preferred method as the filter response does not contain passband ripples.
% Filtfilt is used here as this removes any phase shift between the raw ECG
% and the filtered ECG.
[loaded,b] = butter(1, [5 23] /(fs/2)); % second order filter, cut off at 5-23Hz
filtered_ecg = filtfilt(loaded,b,ecgSignal); % filtfilt used to remove phase shift


% The derivative of the signal is computer using the diff function. This
% part of the pre-processing is used to calculate the slope of the QRS.
diff_ecg = diff(filtered_ecg);


% All elements of the derivative signal are now squared. This enhances the
% dominant QRS peaks to avoid T-waves being wrongly detected in the
% algorithm.
Square_ecg = diff_ecg.^2;


% The signal is now integrated in a moving window to find
% information about the length of the QRS.
moveMean = movmean(Square_ecg,54); % 54 used to match pan-tompkins paper 150ms window


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%***************************************************************************%
%                                                                           %
%                       Thresholdings & Decision                            %
%                                                                           %
%***************************************************************************%


% FIDUCIAL MARK - Finding the peaks in the Integrated signal and the
% filtered signal to be used for initalizing the thresholding.
[pksi,locsi] =  findpeaks(moveMean,'MINPEAKDISTANCE',72);
[pksf,locsf] =  findpeaks(filtered_ecg,'MINPEAKDISTANCE',72);

% Empty Arrays for storing thresholding values.
ThresI1 = [];
ThresF1 = [];
ThresI2 = [];
learningWindow = 5;

% This for loop runs through both the Integrated signal and the
% Filtered signal. It takes a mean and max value every two seconds for
% the length of the signal. The adaptive thresholding lines are
% calculated using variants of the equations from the Pan-Tompkins
% original paper.
for x=1:1:length(moveMean)-(fs*learningWindow)
    NPKI1 = 1/2*(mean(moveMean(x:x+(fs*learningWindow))));
    SPKI1 = 1/3*(max(moveMean(x:x+(fs*learningWindow))));
    ThresI1(x) = NPKI1+0.25*(SPKI1-NPKI1);
    ThresI2(x) = ThresI1(x)*0.5;
    
    
    NPKF1 = 1/2*(mean(filtered_ecg(x:x+(fs*learningWindow))));
    SPKF1 = 1/3*(max(filtered_ecg(x:x+(fs*learningWindow))));
    ThresF1(x) = NPKF1+0.25*(SPKF1-NPKF1);
    ThresF2(x) = ThresF1(x)*0.5;
end


% The thesholding arrays have to be padded - repeating the final value
% of the array until it matches the length of the ecgSignal array. This
% is to avoid any sizing errors.
sampDifferenceI1 = length(ecgSignal)-length(ThresI1); % finds the sample difference
meanThresI1 = mean(ThresI1); % finds the mean of the theshold 1
meanArrayI1 =  meanThresI1.*ones(sampDifferenceI1,1); % creates array mean values
padThresI1 = vertcat(ThresI1',meanArrayI1); %  fully padded new threshold
padThresI2 = padThresI1*0.5;

sampDifferenceF1 = length(ecgSignal)-length(ThresF1); % finds the sample difference
meanThresF1 = mean(ThresF1); % finds the mean of the theshold 1
meanArrayF1 =  meanThresF1.*ones(sampDifferenceF1,1); % creates array mean values
padThresF1 = vertcat(ThresF1',meanArrayF1); %  fully padded new threshold
padThresF2 = padThresF1*0.5;

% Empty arrays initialized to store the peak dectections from the filtered
% signal.
pksiNew = [];
pksiFinal = [];

%For and If statements to logically compute whether the peak is greater
%than the adaptive threshold, for the integrated signal, at the correct
%location. If the peak is greater than the threshold a value of 1 is passed
%to an array and else a 0 is passed to the array.
for i=1:length(locsi)
    if pksi(i) > padThresI1(locsi(i))
        pksiNew = 1;
        pksiFinal = vertcat(pksiFinal, pksiNew);
    else
        pksiNew = 0;
        pksiFinal = vertcat(pksiFinal, pksiNew);
    end
end

% Multiplying the logical peak array and the original peak locations, to
% cause false detections to become 0. These 0 values are then removed.
locRemovali = [pksiFinal.*locsi];
newiLocs = locRemovali(locRemovali~=0);

% Multiplying the logical peak array vs the original peak detections to
% cause false detections to become 0. These 0 values are then removed.
peakRemovali = [pksiFinal.*pksi];
newiPeaks = peakRemovali(peakRemovali~=0);

%Empty arrays initialized to store the peak dectections from the filtered
%signal.
pksfNew = [];
pksfFinal = [];

%For and If statements to logically compute whether the peak is greater
%than the adaptive threshold, for the integrated signal, at the correct
%location. If the peak is greater than the threshold a value of 1 is passed
%to an array and else a 0 is passed to the array.
for i=1:length(locsf)
    if pksf(i) > padThresF1(locsf(i))
        pksfNew = 1;
        pksfFinal = vertcat(pksfFinal, pksfNew);
    else
        pksfNew = 0;
        pksfFinal = vertcat(pksfFinal, pksfNew);
    end
end

% Multiplying the logical peak array vs the original peak locations, to
% cause false detections to become 0. These 0 values are then removed.
locRemovalf = [pksfFinal.*locsf];
newfLocs = locRemovalf(locRemovalf~=0);

% Multiplying the logical peak array vs the original peak detections to
% cause false detections to become 0. These 0 values are then removed.
peakRemovalf = [pksfFinal.*pksf];
newfPeaks = peakRemovalf(peakRemovalf~=0);

% Attempting to implement RR-Avergaing
peaktime = newiLocs./360; % Converting the locations from sample domain to time domain
for i=1:length(peaktime)
    diffRR = diff(peaktime(end-8:end)); % calculate RR interval of the first 8 locations
    mean_RR = mean(diffRR); % calculate the mean of 8 previous R waves interval
end

% Ismember function is used with a tolerance of 20 (reasoning behind
% tolerance has been explained). This function provides a logical array
% according to like values within the tolerance.
[logic] = ismembertol(newiLocs, newfLocs, 20, 'DataScale',1,'ByRows', true, 'OutputAllIndices', true);

% Final value of the locations
finalLocRemoval = [logic.*newiLocs];
finalLocations = finalLocRemoval(finalLocRemoval~=0);

% Final value of the peaks
finalPeakRemoval = [logic.*newiPeaks];
finalPeaks = finalPeakRemoval(finalPeakRemoval~=0);


% Printing the new peaks and their locations
PeaksLocsi = [finalLocations,finalPeaks];
PeaksLocsf = [newfLocs,newfPeaks];


% Calculating the beatheats of each of the patients that have been ran
% through the Algorithm
beats = length(finalLocations);
time = length(ecgSignal)/fs;
bpm = (beats*60)/time;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%***************************************************************************%
%                                                                           %
%                              Performance                                  %
%                                                                           %
%***************************************************************************%


% Ismember function is used with a tolerance of 20 (reasoning behind
% tolerance has been explained). This function provides a logical array
% according to like values within the tolerance..
[compareLogic] = ismembertol(finalLocations,ann,20, 'DataScale', 1, 'ByRows', true, 'OutputAllIndices', true);


% Performance equations carried out to determine the overal performance of
% the algorithm. These results are stored in a structure.
TP = sum(compareLogic ==1); % My detection and clinician are within tolerance
FP = sum(compareLogic ==0); % When my algorithm has detected a peak but clinician did not
FN = abs(length(ann)-length(finalLocations)); % When my algorithm has missed a peak detection from the clinician

% Sensitivity Calculation as a percentage
sensitivity = TP/(TP+FN)*100;

% Calculating the Positive predictive value - which is the probability that
% following a positive test result the result is truly positive.
ppv = TP/(TP+FP);
ppvPercentage = ppv*100;

% Percentage of detections vs the clinicians annotations.
overallPercentage = (TP/length(ann))*100;

%Storing the algorithms results in a 1x1 structure with 9 fields.
structArray1(fullDataBase)= struct('PatientNumber',patientFileNumber,'TP',TP,'FP',FP,...
    'FN',FN,'Sensitivity',sensitivity,'ppv',...
    ppvPercentage,'Percentage_compare_to_ann',overallPercentage,'BPM',bpm);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%***************************************************************************%
%                                                                           %
%                       Plotting Section                                    %
%                                                                           %
%***************************************************************************%



%ECG print test
% figure
% plot(ecgSignal)
% hold on
% plot(ecgAnn,ecgSignal(ecgAnn),'ro')


% This is a plot of the raw ECG signal with the hand annotations printed on
% top. The other plot shows the signal after the bandpass filter.
% figure
% plot(ecgSignal)
% hold on
% grid on
% plot(ann,ecgSignal(ann),'ro')
% title('Original ECG with Annotations')
% xlim([0,5000])
% 
% figure
% sgtitle('Pre-Processing')
% subplot(2,2,1)
% plot(filtered_ecg)
% xlim([0,5000]) % A limit of 5000 samples is set for visual purposes
% title('Filtered ECG - Butterworth')
% 
% 
% subplot(2,2,2)
% plot(diff_ecg)
% xlim([0,5000])
% title('Differiential of the Signal')
% 
% subplot(2,2,3)
% plot(Square_ecg)
% xlim([0,5000])
% title('Square of the signal')

%Moving the mean wth annotations on top
% subplot(2,2,4)
% plot(moveMean)
% title('Moving mean window')
% xlim([0,5000])
% hold on
% plot(ecgAnn,moveMean(ecgAnn),'ro')
% xlim([0,5000])

% figure
% sgtitle('Locating Peaks with findpeaks')
% subplot(2,1,1)
% plot(moveMean);
% hold on
% plot(locsi,pksi,'bx')
% %xlim([0,5000]);
% title('Peaks in the Integated signal')
% subplot(2,1,2)
% plot(filtered_ecg);
% hold on
% plot(locsf,pksf,'rx')
% xlim([0,5000]);
% title('Peaks in the filtered signal')

% Printing the dynamic Thresholds
% figure
% subplot(2,1,1)
% title('Threshold on Integrated Signal')
% plot(moveMean)
% xlim([0,5000]);
% hold on
% plot(ThresI1,'r--')
% xlim([0,5000]);
% subplot(2,1,2)
% title('Threshold on Filtered Signal')
% plot(filtered_ecg)
% xlim([0,5000]);
% hold on
% plot(ThresF1,'b--')
% xlim([0,5000]);
% sgtitle('Dynamic Thresholding')

% Showing Phase delay - will not work without changes at the filtering
%     figure
%     sgtitle('Phase Shift')
%     subplot(2,1,1)
%     plot(ecgSignal)
%     hold on
%     plot(filtered_ecg)
%     xlim([0,5000]);
%     legend('Original ECG','Filtered with filter')
%     title('Filtering with filter')
%     subplot(2,1,2)
%     plot(ecgSignal)
%     hold on
%     plot(filt_ecg)
%     xlim([0,5000]);
%     legend('Original ECG','Filtered with filtfilt')
%     title('Filtering with filtfilt')


% figure
% plot(filtered_ecg);
% hold on
% plot(newfLocs,newfPeaks,'bx');
% hold on
% plot(ThresF1,'--r');
% title('Find peaks in the filtered signal')

% Prints ingrated signal with algorithms peak detections along with the
% clinicians annotations. Adaptive thresholding is also included in this
% plot.
% figure
% plot(moveMean);
% hold on
% plot(finalLocations,finalPeaks,'bx')
% hold on
% plot(ann,moveMean(ann),'ro')
% hold on
% plot(padThresI1,'r--')
% hold on
% plot(ThresI2,'b--')
% legend('moveMean','Detected Peaks','Hand Annotations','ThresI1','ThresI2')
% title('Removing peaks using findpeaks and thresholds')

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%***************************************************************************%
%                                                                           %
%                              Excel                                        %
%                                                                           %
%***************************************************************************%


% Displaying the structure as a table in the command window
table = struct2table(structArray1)

% This function writes the results structure to an excel file for extra
% analysis.
writetable(struct2table(structArray1), 'results.xlsx')

% This section calculates the mean values of the main fields within the
% structure.
avgSensitivity = mean([structArray1.Sensitivity])
avgppv = mean([structArray1.ppv])
avgOverallPercentage = mean([structArray1.Percentage_compare_to_ann])



