 Multi-contrast segmentation

 Author Information: 
 Veronica Corona
 
 Department of Applied Mathematics and Theoretical Physics
 University of Cambridge, UK

 Contact: vc324@cam.ac.uk


 Date: 2018-07-12

 This code is implemented in MATLAB to perform multi-contrast segmentation as 
 described in the preprint of the article at https://arxiv.org/abs/1807.10757.
 If you use this code in your research, please cite this work. 


 The main script is "thalamus_segmentation". 
 Before using this software, you need to install the following packages:
 
 1.prtools at http://prtools.org/software/
   Alternative: replace "nn3=knnc(data,3);" with the matlab  implementation of
   k-Nearest neighbors as "nn3 = fitcknn(double(data),labels,'NumNeighbors',3);".
   Parzen classifier implementation is not available in MATLAB. 

 2.NIfTI_20140122 (outdated). This opens "nii" files. The latest versions of 
   MATLAB have the "niftiread","niftiinfo" and "niftiwrite" built-in
   functions that you can use instead. 


 In case you experience any problems, please contact me at vc324@cam.ac.uk.


