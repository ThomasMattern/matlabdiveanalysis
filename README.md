# matlabdiveanalysis
Series of Matlab scripts to perform automated penguin dive analysis recorded with GPS data loggers.

The scripts are principally designed to analyise data recorded with TechnoSmart AxyTrek GPS dive loggers, 
although the dive analyis itself can be used on any data set that conforms to the file formatting of
AxyTrek data files.

The scripts can accept GPS data as part of the AxyTrek accelerometer file which also contains GPS fixes,
or as a separate GPS file if, for example, AxyDepth TDRs were used in conjunction with other GPS loggers
such as iGotU. Dive analysis can also be performed without any GPS data. However, this requires an empty, 
dummy GPS file.

## 1. Reduce accelerometer data

After conversion of the raw logger data using the X Manager,  the raw accelerometer data needs to be 
'reduced' to only data rows that contain actualy dive data. For example, in files that contain accelerometer 
data recorded at 25Hz, only every 25th row contains the required information for dive analyis.

To run, call the script 'reduce_axy' in Matlab and select the raw accelerometer file. After a few minutes
the reduced data file will be written into the folder containing the accelerometer file adding with the 
suffix '_reduced.txt' to the file name, e.g. '2_S1.csv' will be reduced to '2_S1_reduced.txt'.
