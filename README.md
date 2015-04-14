###################################################################
#                                                                 #
#                           Data Flow                             #
#      Christoph Vogel, Konrad Schindler and Stefan Roth          #
#                          GCPR 2013                              #
#                                                                 #
#         Copyright 2013-2015 ETH Zurich (Christoph Vogel)        #
#                                                                 #
###################################################################



ABOUT:
This software implements our approach to optical flow estimation [1] 
with several data cost functions. 


The additional and optional library
 - Eigen
is not included.

To download that package follow the link:
http://eigen.tuxfamily.org/index.php?title=Main_Page
and read the licensing information provided there.


==========================================================================
DISCLAIMER:
This demo software has been rewritten for the sake of simplifying the
implementation. Therefore, the results produced by the code may differ
from those presented in the papers [1].
In fact the results should be better on the KITTI dataset:
http://www.cvlibs.net/datasets/kitti/.

==========================================================================


IMPORTANT:
If you use this software you should cite the following in any resulting publication:

    [1] An Evaluation of Data Costs for Optical Flow
        C. Vogel, S. Roth and K. Schindler
        In GCPR, Saarbruecken, Germany, September 2013


INSTALLING & RUNNING

1.	(Optional) Download and install eigen from 
	http://eigen.tuxfamily.org/index.php?title=Main_Page 
	and place it into the folder ./Source.
	Alternatively one can change the switch in the file compileMex
	to compile a slightly slower standalone version.
	
2.	Start MATLAB and run compileMex.m to build the utilities binaries.
	(This step can be omitted if you are using Windows 64 bit or Unix 64 bit
	And do not want ot use OpenMP.)
	Adjust the compiler flags accordingly for your purposes 
	(defaults should work in most cases). 
	
3.	From folder DataFlow run calltest( xx ) - 
	example given as comment in the code.
	This will execute a KITTI example.

4.	Otherwise load images I1, I2 and run :
	flow = Data_flow(1.25/255, 12.333, 3, 0.9, I1, I2, 10, 2, 1, 0, 0.5, 16, 0, 1)
	'flow' contains the computed 2d flow as usual.
	Here I1, I2 are the input images (gray-level only so far).
	All parameter are explained in the script, Data_flow.

Note that the code should perform slightly better as published.	
	
CHANGES
	1.0		April 19, 2014	Initial public release
