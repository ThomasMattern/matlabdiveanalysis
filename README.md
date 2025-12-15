# matlabdiveanalysis
Series of Matlab scripts to perform automated penguin dive analysis recorded with GPS data loggers.

The scripts are principally designed to analyise data recorded with TechnoSmart AxyTrek GPS dive loggers, 
although the dive analyis itself can be used on any data set that conforms to the file formatting of
AxyTrek data files.

The scripts can accept GPS data as part of the AxyTrek accelerometer file which also contains GPS fixes,
or as a separate GPS file if, for example, AxyDepth TDRs were used in conjunction with other GPS loggers
such as iGotU. Dive analysis can also be performed without any GPS data. However, this requires an empty, 
dummy GPS file.

## 1. Reduce accelerometer data (reduce_axy.m)

After conversion of the raw logger data using the X Manager,  the raw accelerometer data needs to be 
'reduced' to only data rows that contain actualy dive data. For example, in files that contain accelerometer 
data recorded at 25Hz, only every 25th row contains the required information for dive analyis.

To run, call the script `reduce_axy.m` in Matlab and select the raw accelerometer file. After a few minutes
the reduced data file will be written into the folder containing the accelerometer file adding with the 
suffix '_reduced.txt' to the file name, e.g. '2_S1.csv' will be reduced to '2_S1_reduced.txt'.

Note: The script assumes that data is exported in X Manager without any Metadata, pressure data coverted 
to metres (rather than mbar), with date in the format dd/mm/yyy and 24 hour time.

## 2. Dive Analysis (read_data.m & dive_analysis.m)

To perform a dive analysis of recorded GPS logger data, some data curation is necessary. Data from each 
deployment must be stored in a folder named after the deployment site. Each site folder should contain 
subfolders for individual deployments, named using a unique identifier (typically a bird ID).

Within each bird ID subfolder, the script expects to only find the reduced dive data file produced with
the reduce_axy script. The file must be renamed to 'dive.txt'.

For example, two deployments performed at Urenui Beach (New Zealand) would be stored as follows:

- Urenui
  - female_40195_GPS05
    - dive.txt
  - male_40549_GPS03
    - dive.txt
  
Run the `read_data.m` script and select the site folder (e.g., "Urenui" in the example above). The script will:

1. Iterate through all subfolders.
2. Perform a dive analysis on the contained data files.
3. Prompt for a name and location to save the dive analysis data in CSV format.

Note: this scripts requires helper scripts to work - baselinetracking.m, bottomtrap.m, and pos2dist.m. These 
files should all be present in the toolbox folder.

## 3. Extract raw GPS Data based of Dive Analysis (trip_fixes.m)

This auxiliary script extracts raw GPS data for each foraging trip identified during the dive analysis. The 
dive analysis discards any GPS fixes not associated with diving behavior (e.g., while drifting at the surface).

The script requires:

- A completed, unaltered dive analysis file.
- Access to the original dive and GPS data.

Run `trip_fixes.m` and:

1. Select the dive analysis CSV file generated in Step 2.
2. Select the site folder containing the raw dive and GPS data.

The script will output a new CSV file containing raw GPS fixes, named `{SITENAME}_gps.txt`

## 4. Visualise dive data (visdiv.m)

Another auxilary script that produces a graph that allows dive data to be explored visually.
