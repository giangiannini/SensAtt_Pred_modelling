# SensAtt_Pred_modelling
This repository contains some code that was used to model and fit the behavioural responses of the task described in [Giannini et al., 2025](https://www.nature.com/articles/s41598-025-87244-9). 

The function ___ allows to model the updating of the beliefs of an agent undergoing a probabilistic learning task in which a stimulation (or not) is given based on an underlying hidden state over series of n trials. Each observation (stimulation or not) updates a Beta-Bernoulli distribution. The mean of the distribution mimics the belief of the participant, while its variance mimics the confidence of the agent.  

In the experiment, participants were required to learn the probabilistic state over sequences of 25 trials and to give an explicit estimate of their beliefs at the end of each sequence. Therefore, we obtained an estimate of the beliefs of each participant at the end of each sequence, through which we could fit our model parameters to. 

More in detail, the model is equipped with a learning parameter (k) that regulates how much each new observation weights on the posterior probability at each trial. The learning parameter is assumed to be fixed across the experiment. The learning parameter that best explains the behavioural responses obtained from participants at the end of each sequence is found through minimisation of negative log likelihood, implemented in MATLAB via the fminsearch function. 
