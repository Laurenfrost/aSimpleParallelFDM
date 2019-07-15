# aSimpleParallelFDM
Parallelize some code for some postgraduate bros in DUT  

## Software Enviroment
CentOS 7.6 1810 version  
Intel(R) C++ Compiler 2019 update 3  
Intel(R) MKL(R) 2019 update 3  
Intel(R) MPI 2019 update 3  

## Hardware Enviroment
Intel(R) Xeon(R) CPU E5-2660 v3 @ 2.60GHz  

# What is FDM
Acturally I'm not very clear about FDM, the shorter name of "Finite Differential Method".  
As far as I know, the FDM bases on the following forma:  
$$\frac{\ddot{a} T}{\ddot{a} t}=(\frac{\sqrt \ddot{a} T}{\ddot{a}\sqrt x}+\frac{\sqrt \ddot{a} T}{\ddot{a} \sqrt t})$$
