function run_strfRTA(resp, rawStim)
w = size(rawStim,1);
T = size(rawStim,3);
rawStim = reshape(rawStim, [w w T]);

params.RTAC = [1 0];
[stimtrain,params] = preprocRTAC(rawStim, params);

% We create a Generalized Linear Model (GLM) with the linInit function. 
% The first argument is the number of parameters in the fitting, and the second argument is 
% the number of delays to consider between stimulus and response
strf = linInit(size(stimtrain,2), [0:8]);

% Set the bias term of the GLM to the mean response. This generally
% speeds up the fitting process
strf.b1 = mean(resp);
% preprocessing parameters
strf.params = params;

% STRFlab uses a global variable to share data between functions.
global globDat

% strfData puts the stimulus and response that will be used to fit the
% model into the globab variable
strfData(stimtrain,resp)

% To train this model we are going to use the scaled conjugate gradient
% technique in the function trnSCG. Calling it with no arguments returns
% the default set of options, which can then be edited. We will change it
% to graphically display every five iterations and only perform at most
% 100 iterations
options=trnSCG;
options.display=-5;
options.maxIter = 100;

% Create an index of the samples used to train the model. 
trainingIdx = [1:globDat.nSample];

% strfOpt is a general routine that takes any strf structure and set of
% options and uses the appropriate error and fitting functions. It returns
% a new structure with the model parameters fit to the data. Because of
% the display option above, this will open a window showing model 
% parameter values and error values
figure
strfTrained=strfOpt(strf,trainingIdx,options);

% Now we can use the visualization routine that corresponds to the preproc
% routine to visualize the STRF
preprocRTAC_vis(strfTrained);
