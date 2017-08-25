
# ChromAberFix (ICCV 2017)

This repository contains the code for the ICCV 2017 paper "Revisiting Cross-channel Information Transfer for Chromatic Aberration Correction" authored by Tiancheng Sun, Yifan Peng, and Wolfgang Heidrich. [Project page](http://vccimaging.org/Publications/Sun2017RCI/)


## Getting Started

We use Matlab 2016 to implement the algorithm. In the code, we use the open source package [mmx](https://www.mathworks.com/matlabcentral/fileexchange/37515-mmx-multithreaded-matrix-operations-on-n-d-matrices) and a slightly modified version of [im2col](https://github.com/yjxiong/im2col). You can find these packages in the folder `util/`. Please run the code `build.m` inside to initialize these packages.


## Running the tests


We provide an example image taken from diffractive lens. You can use the code `test.m` to test the algorithm. The core concept is to transfer the information from green channel (which is less burry) to other channels. 

We also provide a blind deconvolution version of the green channel (deconvolved using [Normalized Sparsity Measure](http://cs.nyu.edu/~dilip/research/blind-deconvolution/)). You can toggle the variable `use_blind` in `test.m` to deblur the images.


# License


This code and data is released under the Creative Commons Attribution-NonCommercial 4.0 International license (CC BY-NC.) In a nutshell:

- The license is only for non-commercial use (commercial licenses can be obtained from KAUST, please contact the authors for details)
- The material is provided as-is, with no warranties whatsoever.
- If you publish any code, data, or scientific work based on this, cite our work. 

# Acknowledgements

This work was supported by KAUST baseline funding, as well as a UBC 4YF Doctoral Fellowship. The authors thank Tao Yue, Qiang Fu, and Felix Heide for the help on synthetic results.

If you find any bugs or have comments/questions, please contact 

	Tiancheng Sun      [kevin.kingo0627@gmail.com]
	
	Yifan (Evan) Peng  [evanpeng@cs.ubc.ca]
	
	Wolfgang Heidrich  [wolfgang.heidrich@kaust.edu.sa]

Aug.2017
