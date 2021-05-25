# FairOT
This code is the implementation of fair distributed optimal transport using ADMoM.
We have two algoriths shown with ex3 and ex4. Ex3 is the standard algorithm, while ex4 is the online algorithm, meaning the parameters can be updated mid calculation
and the optimal solution can still be found. 

OTPRIMALIncomplete_fair.m is the central planner function (i.e. it does not use the distributed algorithm to calculate optimal transport solution)
We incorporate this to compare the solution of the optimal solution with that of a central planner. 

OTRelaxPrimalADMoMPrivateFormal_Fair.m is the distributed algorithm.
OTRelaxPrimalADMoMPrivateFormalOnline_Fair.m is the online version if the distributed algorithm. 

Ex3 sets up the use of the standard distributed algorithm.
Ex4 sets up the use of the oniline distributed algorithm.
