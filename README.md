# contrast-openmp
Given an image in PNM (P5 or P6) format. The program increases the contrast of the image by stretching the range of used colors, according to user-chosen coefficient. For better performance, OpenMP 2.0 was used (a standard C++ module for parallel programming).
## For running
You will need CLion 2021.3.3 with g++9 and CMake 3.21.1.

1) Open the project (image-processor) as a C++ project in CLion.

2) Put a PNM file in the folder cmake-build-debug.

3) Find a file main.cpp. Now you should specify program arguments as follows:

<number_of_threads> <input_file_name> <output_file_name> <coeff>, where coeff is a float that belongs to [0; 0.5)

You can find more information about program in report.pdf (available in Russian only).
