# Example prediction code for MATLAB for the PhysioNet/CinC Challenge 2020

## Contents

This code uses two main scripts to train the model and classify the data:

* `train_model.m` Train your model. Add your model code to the `train_12ECG_model` function. It also performs all file input and output. **Do not** edit this script or we will be unable to evaluate your submission.
* `driver.m` is the classifier which calls the output from your train_model. It also performs all file input and output. **Do not** edit this script or we will be unable to evaluate your submission.

Check the code in these files for the input and output formats for the `train_model` and `driver` scripts.

To create and save your model, you should edit `train_12ECG_classifier.m` script. Note that you should not change the input arguments of the `train_12ECG_classifier` function or add output arguments. The needed models and parameters should be saved in a separated file. In the sample code, an additional script, `get_12ECG_features.m`, is used to extract hand-crafted features. 

To run your classifier, you should edit the `run_12ECG_classifier.m` script, which takes a single recording as input and outputs the predicted classes and probabilities. Please, keep the formats of both outputs as they are shown in the example. You should not change the inputs and outputs of `run_12ECG_classifier` function. The needed models and parameters can be loaded in the `load_12ECG_model` function.

Check the code in these files for the input and output formats for the `load_12ECG_model` and `run_12ECG_classifier` functions.

## Running

You can run this prediction code by starting MATLAB and running

    train_model(training_data, model)
    driver(model, test_data, test_outputs)

where `training_data` is a directory of training data files, `model` is a directory of files for the model, `test_data` is the directory of test data files, and `test_outputs` is a directory of classifier outputs.  The [PhysioNet/CinC 2020 webpage](https://physionetchallenges.github.io/2020/) provides a training database with data files and a description of the contents and structure of these files.

## Submission

The `driver.m`, `get_12ECG_score.m`, and `get_12ECG_features.m` scripts must be in the root path of your repository. If they are inside a folder, then the submission will be unsuccessful.

## Details

The code uses three main toolboxes:
- HRV toolbox to compute the RR intervals. https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox.git. 
  "An Open Source Benchmarked Toolbox for Cardiovascular Waveform and Interval Analysis", 
   Physiological measurement 39, no. 10 (2018): 105004. DOI:10.5281/zenodo.1243111; 2018. 
 - ECGkit to find the ECG fiducial points: https://github.com/marianux/ecg-kit.git
  Demski AJ, Llamedo Soria M. "ecg-kit a Matlab Toolbox for Cardiovascular Signal Processing".  
  Journal of Open Research Software. 2016;4(1):e8. DOI: http://doi.org/10.5334/jors.86
- GEH parameter extraction and origin point: https://github.com/Tereshchenkolab/Global-Electrical-Heterogeneity.git and https://github.com/Tereshchenkolab/Origin.git. 
  Perez-Alday,et al; "Importance of the Heart Vector Origin Point Definition for an ECG analysis: 
  The Atherosclerosis Risk in Communities (ARIC) study". Comp Biol Med, Volume 104, January 2019, 
  pages 127-138. https://doi.org/10.1016/j.compbiomed.2018.11.013
  Waks JW, et al. "Global Electric Heterogeneity Risk Score for Prediction of Sudden Cardiac Death in the General Population: 
  The Atherosclerosis Risk in Communities (ARIC) and Cardiovascular Health (CHS) Studies". Circulation. 2016;133:2222-2234.
