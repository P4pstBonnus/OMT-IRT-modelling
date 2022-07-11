# OMT-IRT-model
This repository contains the code for the masterthesis regarding the Item-Response-Theory-models for the operant motive test.

## Utils
The python transformation function can be found in transformation_function.py.
The needed instructions to install pandas can be found at https://pypi.org/project/pandas/.
The transformation can easily be extended to work for motivelevel data as well, but since there are many ways to achieve the needed data structure, only the functionality to transform motive data is provided.

The R-utility functions can be found in glmm_model_fitting_utils.R.
The needed installs can be done via\
install.packages("glmmTMB", type="source")\
install.packages("mvtnorm", type="source")\
in a R-environment.

## Examples
There are two examples included.
The first example is a model iteration for the omt motive-coded data and can be found in example_main_motives.R.
The second example is a model iteration for the omt motivelevel-coded data and can be found in example_motive_level.R.
For further insights and instructions pls refer to the masterthesis (contact repository owner).

## Data
Three data sources are provided for ten imaginary persons to get the idea of the needed data structure.
original_omt_data_motive.csv contains the input for the transformation via the python file.
example_omt_data.csv contains the transformed data input for the example_main_motives.R file.
example_omt_data_motivelevel.csv contains the transformed data input for the example_motive_level.R file.
