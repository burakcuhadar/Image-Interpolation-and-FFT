# Image Interpolation and Discrete Fourier Transformation

In this project I have implemented several interpolation algorithms and inverse fast fourier transform.

**How to test the code:** There is a class called ImageViewer which has a main method. This class provides a graphical
interface to load your own images and test the resizing of the image with different interpolation methods. The methods I 
implemented can also be tested numerically by writing tests like those in Test.java and Test_Interpolation.java files.

# Example outputs:
* The image in bilder/test.bmp before resizing, opened with ImageViewer:
![alt text](https://github.com/burakcuhadar/Image-Interpolation-and-FFT/blob/master/example/before_resizing.png)

* The image is scaled up by a factor of 2 using nearest neighbor method:
![alt text](https://github.com/burakcuhadar/Image-Interpolation-and-FFT/blob/master/example/after_resizing_with_nearest.png)

* The image is scaled up by a factor of 2 using cubic splines method:
![alt text](https://github.com/burakcuhadar/Image-Interpolation-and-FFT/blob/master/example/after_resizing_with_cubic.png)



The following methods are implemented by myself:
* **evaluate** method in **LinearInterpolation** class
* **computeCoefficients, addSamplingPoint** and **evaluate** methods in **NewtonPolynom** class
* **computeDerivatives** and **evaluate** methods in **CubicSpline** class
* **add, sub, mul** and **fromPolar** methods in **Complex** class which can be found in dft folder
* **ifft** method in **IFFT** class which can be found in dft folder
