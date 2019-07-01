==================================================================================
THIS README FILE ADDRESSES THE FOLLOWING WARNING MESSAGE:

Warning! Negative diagonal elements of (I+dA); the estimate may not be meaningful.

OR ERROR MESSAGE:

 Error in sample.int(length(x), size, replace, prob) : 
  NA in probability vector 

WHICH CAN OCCUR WHEN USING THE Markov AND SemiMarkov FUNCTIONS. ITS OCCURRENCE CAN RESULT IN NEGATIVE STATE OCCUPANCY PROBABILITIES OR PROBABILITIES THAT ARE NOT ABLE TO BE
CALCULATED 
================================================================================== 


The Markov or SemiMarkov function, as appropriate, require that the increment of the cumulative hazard between consecutive timepoints is below one. The above warning message occurs when this is not the case for every pair of consecutive time points.

The first step to addressing this issue is to use the function MarkovHaz or SemiMarkovHaz, as appropriate, to determine which timepoints are causing the problem. The function should initially be used with the same timeseq and timeseq_ext arguments as used with the Markov or SemiMarkov function. This will show the cumulative hazards, for each transition, at each time in the same sequence of timepoints that was used by the Markov/SemiMarkov function, and allow identification of the problematic timepoints.

Running the following commands after using the function MarkovHaz or SemiMarkovHaz can help with identification of problems:

## assuming the result of MarkovHaz or SemiMarkovHaz is called results

## identify timepoints where the difference in consecutive cumulative hazards >1
which(diff(results$Haz)>1 & diff(results$trans)==0)

## identify timepoints where the difference in consecutive cumulative hazards <0
which(diff(results$Haz)<0 & diff(results$trans)==0)


## show the timepoints of interest to help identify where timeseq will need to expand
## for example,

results[2265:2295,]
results[4510:4530,]





