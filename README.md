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

<a href="https://www.codecogs.com/eqnedit.php?latex=$$\frac{\partial&space;T}{\partial&space;t}=&space;\alpha&space;(\frac{\partial^2&space;T}{\partial&space;x^2}&space;&plus;&space;\frac{\partial^2&space;T}{\partial&space;y^2})$$" target="_blank"><img src="https://latex.codecogs.com/gif.latex?$$\frac{\partial&space;T}{\partial&space;t}=&space;\alpha&space;(\frac{\partial^2&space;T}{\partial&space;x^2}&space;&plus;&space;\frac{\partial^2&space;T}{\partial&space;y^2})$$" title="$$\frac{\partial T}{\partial t}= \alpha (\frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial y^2})$$" /></a>  
>$$\frac{\partial T}{\partial t}= \alpha (\frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial y^2})$$  

Further infomation is beyond my capacibility.  
