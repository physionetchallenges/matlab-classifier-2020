# Example prediction code for MATLAB for the PhysioNet/CinC Challenge 2020

## Contents

This prediction code uses three scripts:

* `get_12ECG_score.m` makes predictions on clinical time-series data.  Add your prediction code to the `get_12ECG_score` function. It calls all the functions inside the `functions` folder 
* `load_12ECG_model.m` loads model weights, etc. for making predictions.  To reduce your code's run time, add any code to the `load_12ECG_model` function that you only need to run once, such as loading weights for your model.
* `driver.m` calls `load_12ECG_model` once and `get_12ECG_score` many times. It also performs all file input and output.  **Do not** edit this script -- or we will be unable to evaluate your submission.

Check the code in these files for the input and output formats for the `load_12ECG_model` and `get_12ECG_score` functions.

## Running

You can run this prediction code by starting MATLAB and running

    driver(input_directory, output_directory)

where `input_directory` is a directory for input data files and `output_directory` is a directory for output prediction files.  The PhysioNet/CinC 2020 webpage provides a training database with data files and a description of the contents and structure of these files.

## Submission

The driver.py, get_12ECG_score.py, and get_12ECG_features.py scripts to be in the root path of the Github repository. If they are inside a folder, then the submission will fail.

## Details

The code uses three main toolboxes:
- HRV toolbox to compute the RR intervals. https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox.git
  "An Open Source Benchmarked Toolbox for Cardiovascular Waveform and Interval Analysis", 
   Physiological measurement 39, no. 10 (2018): 105004. DOI:10.5281/zenodo.1243111; 2018. 
 - ECGkit to find the ECG fiducial points: https://github.com/marianux/ecg-kit.git
  Demski AJ, Llamedo Soria M. "ecg-kit a Matlab Toolbox for Cardiovascular Signal Processing". 
  Journal of Open Research Software. 2016;4(1):e8. DOI: http://doi.org/10.5334/jors.86
- GEH parameter extraction and origin point: https://github.com/Tereshchenkolab/Global-Electrical-Heterogeneity.git and https://github.com/Tereshchenkolab/Origin.git
  Erick Perez-Alday,et al; "Importance of the Heart Vector Origin Point Definition for an ECG analysis: 
  The Atherosclerosis Risk in Communities (ARIC) study". Comp Biol Med, Volume 104, January 2019, 
  pages 127-138. https://doi.org/10.1016/j.compbiomed.2018.11.013
  Waks JW, et al. "Global Electric Heterogeneity Risk Score for Prediction of Sudden Cardiac Death in the General Population: 
  The Atherosclerosis Risk in Communities (ARIC) and Cardiovascular Health (CHS) Studies". Circulation. 2016;133:2222-2234.

 


