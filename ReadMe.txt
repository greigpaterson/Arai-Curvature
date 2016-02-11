This function calculate the curvature of Arai plot data following Paterson [2011, JGR].

The routine using best-fit circle routines from:
[1] Taubin, G. (1991), Estimation of planar curves, surfaces, and nonplanar space curves defned by implicit equations with applications to edge and range image segmentation, IEEE Trans. Pattern Anal. Mach. Intell., 13, 1115-1138, doi:10.1109/34.103273.
[2] Chernov, N., and C. Lesort (2005), Least squares fitting of circles, J. Math. Imaging Vision, 23, 239-252, doi:10.1007/s10851-005-0482-8.

Usage:
Download all three MATLAB scripts into the same folder. Start MATLAB and navigate to this folder.
The main function is "AraiCurvature.m" and it can be called using:

[parameters]=AraiCurvature(TRM, NRM)

The TRM and NRM inputs do not have to be normalized, the function does this.
The output, "parameters", is a four element vector containing:

1) The signed curvature
2) The x-coordinate of the circle centre
3) The y-coordinate of the circle centre
4) The sum of the squares of the errors (SSE, Paterson, 2011, equation 1)


If any problems are encountered please e-mail me (greig.paterson@mail.iggcas.ac.cn)
