# Catmull-s-Hidden-Surface-Algorithm

In this project, part of the paper titled "A hidden Surface Algorithm with 
anti-aliasing" authored by Edwin Catmull\cite{0} has been implemented. In this 
paper, an anti aliasing technique addressing hidden surface in image has been 
introduced to reduce symptoms of aliasing effectively when data for the pixel is 
complicated. To display an image which will always have much higher resolution than 
corresponding display, it need to be sampled at discrete points corresponding to 
pixels. These image space objects have sharpness and details (higher frequencies) 
that cannot be possibly reproduced due to lower sampling rate than required 
(Nyquist-shannon theorem). It is the attempt to sample that detail at discrete 
points in the image that causes the problem. Hence, the image should be filtered to 
filter these fine details before it can mess up attainable lower frequencies as 
aliases after which sampling can be done. 

One simple filter that is easy to implement analytically is two dimensional Fourier(box) window. Convolving with this filter in frequency domain is equivalent to integrating i.e. taking average visible intensities over the area of each pixel in spatial domain. To correctly integrate intensities of only visible objects at a pixel,a hidden surface algorithm at every pixel is also required. Therefore, what we need is an analytic continuous solution to both hidden surface algorithm and the filter convolution and this paper provides an algorithm for that.

How To Run:
1. Go to the Source folder.
2. Run the .exe file.
3. Type any number between 1 to 5 to select one of the following options to see results presented in the report in the Report folder.
	1. Figure 1(b)
	2. Figure 2(b)
	3. Figure 3(b)
	4. Figure 4(b)
	5. Figure 4(c)
