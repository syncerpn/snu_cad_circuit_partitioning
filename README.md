# SNU Introduction to Computer-Aided Design
Fall 2018

Term Project: Circuit Partitioning

Team 5:
- Nguyen Duy Thanh  : 2017-31828
- Nguyen Tuan Nghia : 2018-21525

== IMPLEMENTATION ==
We have implemented 3 partitioning algorithms: FM, KL and SA.
The code can be found in 'src/partition_alg.c' and 'src/partition_alg.h'.
We divide the file 'src/partition_alg.c' into sections for different algorithms.
For each algorithm, any related functions fall into the same section.
We also add comments for your reading convenience.
Visually, the file 'src/partition_alg.c' looks like this:


//======================  FM  ==============================
/* CODE FOR FM ALGORITHM*/

//======================  KL  ==============================
/* CODE FOR KL ALGORITHM*/

//======================  SA  ==============================
/* CODE FOR SA ALGORITHM*/


There are some additional helpers for the partitioning algorithms.
These helpers can be at the bottom of 'src/util.c' and at the beginning of 'src/partition_alg.c'.

== COMPILING GUIDE ==
We have included everything needed for compilation in makefile.
You can simply run 'make' in 'src' folder.
We sucessfully compiled and tested our code on the following environment:
-- OS: Ubuntu 18.04 LTS
-- Compiler: gcc version 7.3.0 (Ubuntu 7.3.0-27ubuntu1~18.04)
-- CPU: Intel(R) Core(TM) i5-4670 @ 3.4GHz
-- RAM: 8GBs

== OUTPUT RESULT ==
For each algorithm, the test is repeated 10 times.
Both the min and average cut sizes are reported, along with area ratio.
Cut sizes are computed on nets, thus they can be compared with each other.
