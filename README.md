# EE6403-KLT-Matlab
EEE-EE6403 assignment, the KLT matlab code

KLT, the Karhunen loeve transform, is a method that makes transformation more effective by replacing the coordinates. By calculating the covariance matrix, the eigen value, which is like a shadow dropped to different directions, are obtained. With the eigenvalue and eigenvector, the most out-standing coordinate can be calculated, and the information amount along each direction can be measured. Such process makes it possible to compress the image with the highest efficience.

The matlab code usage is simple.

Details are in the `.m` file. The most simple usage is: `avg_err = KLT('./lena.png', 4095, 64)`, which partition the 512*512 lena image into blocks of size 64*64, and select the biggest 4095 eigen values out of 4096 total values.

Furthur enquiry can PM me. 
