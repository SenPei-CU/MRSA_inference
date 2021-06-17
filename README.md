# MRSA_inference

Code for running a synthetic test in Pei S., Liljeros F., Shaman J. Identifying asymptomatic spreaders of antimicrobial-resistant pathogens in hospital settings.

To run the code, first compile the C programs in MATLAB using "mex agentbasedprob.c" and "mex agentbasedsimulation". Then run the function "SILIinference()". Explanations on data, variables, and algorithm details can be found in the comments in the function SILIinference.

The folder "net" contains example contact networks in 52 weeks for synthetic testing. Each mat file contains the network structure for one week. "deg" (degree) records the numbers of connections for all patients in that week. "nl" (neighbor list) and "part" (partition) jointly specify the network structure. The first column of "nl" lists all neighbors of patients; columns 2 to 8 are binary variables (0 or 1) that indicate whether this neighbor is connected with the focal patient. "part" records the rows of neighbors corresponding to each patient. For instance, the neighbors of patient i in the network is nl(part(i):part(i+1)-1,1). If nl(part(i)+x, 1+t)==1, where 0<=x<=part(i+1)-part(i)-1, this neighbor nl(part(i)+x,1) is connected to patient i on day t of this week; otherwise nl(part(i)+x,1) is not connected with patient i on day t of this week.
