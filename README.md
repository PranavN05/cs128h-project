Group Name - Rustaceans

Pranav Narayanan(pn21), Aarav Agarwal(aarav3), Arha Gatram(agatram2)

Project Introduction - Fast Fourier Transform

This project implements the commonly used Fast Fourier Transform, a computationally efficient version of the Discrete Fourier Transform. This Fourier transform turns a function of time or distance into a function of frequency, which describes the most present frequencies in the initial dataset. The DFT does this for discrete datasets. The goal of this project is to read in data from an input file, and output the Fourier transformed data in a new file. We chose this project because we noticed that it makes use of linear algebra, and because it is relevant in many fields and we thought it would be useful to get familiar with it.

Technical Overview: 
This project makes use of the num_complex crate as the Fourier Transform involves complex numbers

  By Checkpoint 1:
  
  Implement the bit-reversal method for the input dataset. This step works for the algorithm operating on datasets of a length that is a power of 2. In essence, all the indices are converted to binary, their digits are reversed, and they are reordered according to the reversed binary. This sets up the dataset to be transformed.
    
  Implement the butterfly method that takes two pairs of DFTs and computes the combined DFT for the total of that dataset. In the FFT, this is known as the butterfly step. It recomputes the transform of each input value by computing a linear combination of one value from each input DFT, using the relevant root of unity to scale these values.
    
  By Checkpoint 2:
  
  Set up a main method that calls the bit-reversal and combiner functions in the necessary order based on the input dataset. Each call to the butterfly will be done in a new thread in order to make use of parallelism.
    
  Implement methods to read data in from a file and store in a vector, and to output a vector of data into a file. These files will be whitespace delimited sets of data. Alongside these methods will be a method that ensures that the dataset is a power of two, in order for this algorithm to be applicable.
    
  By Deadline:
  
  Will implement Chirrp-Z Transform (CZT), which computes DFT for more general datasets. This requires implementing the inverse FFT, which makes use of the same functions as FFT. It will also require a method that zero-pads the dataset such that it becomes a power of two. The CZT method that we will write combines the FFT methods and the iFFT to compute DFT for non "power-of-2" datasets.
    
  Will refactor the main method and file reads to allow for datasets that are not base-2.

Possible Challenges

  We may have to re-implement the butterfly method 
