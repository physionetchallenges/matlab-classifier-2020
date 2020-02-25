# Example MATLAB classifier for the PhysioNet/CinC Challenge 2020

## Contents

This classifier uses three scripts:

* `run_12ECG_classifier.m` makes classifications on 12-Leads ECG data.  Add your prediction code to the `run_12ECG_classifier` function. `load_12ECG_model.m` loads model weights, etc. for making classifications.  To reduce your code's run time, add any code to the `load_12ECG_model` function that you only need to run once, such as loading weights for your model.
* `get_12ECG_features.py` extract the features from the clinical time-series data. This script and function are optional, but we have included it as an example. It calls all the functions inside the `Tools` folder 
* `driver.m` calls `load_12ECG_model` once and `run_12ECG_classifier` many times. It also performs all file input and output.  **Do not** edit this script -- or we will be unable to evaluate your submission.

Check the code in these files for the input and output formats for the `load_12ECG_model` and `run_12ECG_classifier` functions.

## Running

You can run this classifier code by starting MATLAB and running

    driver(input_directory, output_directory)

where `input_directory` is a directory for input data files and `output_directory` is a directory for output classification files.  The PhysioNet/CinC 2020 webpage provides a training database with data files and a description of the contents and structure of these files.

## Submission

The `driver.m`, `get_12ECG_score.m`, and `get_12ECG_features.m` scripts need to be in the base or root path of the Github repository. If they are inside a subfolder, then the submission will fail.

## Details
“The baseline classifiers are simple logistic regression models. They use global electrical heterogeneity (GEH) computed from the WFDB signal file (the `.mat` file) with the [PhysioNet Cardiovascular Signal Toolbox] and demographic data taken directly from the WFDB header file (the `.hea` file) as predictors. 

The code uses three main toolboxes:
- HRV toolbox to compute the RR intervals: https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox.git. 
  "An Open Source Benchmarked Toolbox for Cardiovascular Waveform and Interval Analysis", 
   Physiological measurement 39, no. 10 (2018): 105004. DOI:10.5281/zenodo.1243111; 2018. 
 - ECG-kit to find the ECG fiducial points: https://github.com/marianux/ecg-kit.git
  Demski AJ, Llamedo Soria M. "ecg-kit: a Matlab Toolbox for Cardiovascular Signal Processing".  
  Journal of Open Research Software. 2016;4(1):e8. DOI: http://doi.org/10.5334/jors.86
- GEH parameter extraction and origin point: https://github.com/Tereshchenkolab/Global-Electrical-Heterogeneity.git and https://github.com/Tereshchenkolab/Origin.git. 
  Perez-Alday, et al. "Importance of the Heart Vector Origin Point Definition for an ECG analysis: 
  The Atherosclerosis Risk in Communities (ARIC) study". Comp Biol Med, Volume 104, January 2019, 
  pages 127-138. https://doi.org/10.1016/j.compbiomed.2018.11.013
  Waks JW, et al. "Global Electric Heterogeneity Risk Score for Prediction of Sudden Cardiac Death in the General Population: 
  The Atherosclerosis Risk in Communities (ARIC) and Cardiovascular Health (CHS) Studies". Circulation. 2016;133:2222-2234.
