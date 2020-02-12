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

## Details

See the PhysioNet webpage for more details, including instructions for the other files in this repository.
