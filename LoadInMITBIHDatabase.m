% add the data URl to the code:

% will download the zip file of the package
wfdb_url='https://physionet.org/physiotools/matlab/wfdb-app-matlab/wfdb-app-toolbox-0-10-0.zip';
% will save the web URL in the directory (you can change the path to it)
outfilename = websave('wfdb-app-toolbox-0-10-0.zip',wfdb_url);
% will unzip the zip file
unzip('wfdb-app-toolbox-0-10-0.zip');
% will add the path of teh unzipped folder to the code
addpath ('mcode')

% download the Records file and find the names in the file: 
urlwrite('https://physionet.org/files/mitdb/1.0.0/RECORDS', 'Records.txt');
Records = load ('Records.txt', '-ascii'); % these would be the name of all the files
% this makes a new folder in the directory to save the data in
mkdir('data')
for i = 1:length(Records) 
    %each file name
    file = num2str(Records(i));
    %load the data
     [signal,~,~]=rdsamp(['mitdb/' file]);
     signal = signal(: , 1); % save one channel of the signal
     %load the annotations
       [ann , anntype] = rdann(['mitdb/' file], 'atr');
       %create a structured dataset for each patient
       eval(['pat_', file, '.sig = signal' ]) ;
     eval(['pat_', file, '.ann = ann' ]) ;
     eval(['pat_', file, '.annType = anntype' ]) ;
     % save the data in data folder made in line 16
     save(['data/pat_',  file, '.mat'], ['pat_', file])

end
