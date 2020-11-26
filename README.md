# NNUAL
Neural Network with Unknown Activation and Link Function
## Description
We develop a supervised functional principal component method that extract the latent scores based on the relationship between functional covariates and response.
## Usage
NNUAL(ZZnchu,TTnchu,YYnchu,ntrain)
## Arguments
**ZZnchu**  &nbsp;&nbsp;  an list containing sampled values of curves.  
**TTnchu**  &nbsp;&nbsp;  an list containing sampled time points.  
**YYnchu**  &nbsp;&nbsp;  an list containing response.  
**ntrain**  &nbsp;&nbsp;&nbsp;&nbsp;  size of training sample.   
## Value
An object of the GAMiFM class containing:  

**ypre** &nbsp;&nbsp;&nbsp;&nbsp; the prediction error.  
**B**    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; the matrix of multi-index coefficient.  
**scores**  &nbsp; the matrix of scores.   
**gf**  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; a functional data object for link function.  
**psaif**  &nbsp;&nbsp;&nbsp;&nbsp; a functional data object for additive component functions.
