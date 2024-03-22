Group Name - Rust Aceans

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

  * We may have to re-implement the butterfly method
  * We will struggle to understand and have the patience to go through our long research papers we want to use
  * We will have to understand computer architecture to effectively utilize the full performance of the CPU

References:

  * [Construction of a High-Performance FFT](https://edp.org/work/Construction.pdf)
  * [Bluestein Design Document](https://rocm.docs.amd.com/projects/rocFFT/en/latest/design/bluestein.html)
  * [FFTs: A Tutorial Review](https://pdf.sciencedirectassets.com/271605/1-s2.0-S0165168400X01685/1-s2.0-016516849090158U/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjEDUaCXVzLWVhc3QtMSJHMEUCIQDb4owOzQw%2FICXzk8kCA9%2BsfAXR%2FBOvhhawAtUzDNZHmwIgAt%2Fd6VC4AvDEQRhQEVPp%2BvuBdzCh2hymZzcVoGEOc5wqswUITRAFGgwwNTkwMDM1NDY4NjUiDHjdUjKrxvmSxwftmCqQBULFlY5eGf73kAN9JI0WvCMyzg7gTTiX8b6bBpbrC89VPFj3WXu0ZAt%2F%2FT6pEG%2B4bAH5cgoglerXa3rXTx8TFoY5%2BqJX40PSSS5Pp2JGddKCGM%2BphcQspqJ6Nrhe8yo%2FcWqzVxmOh8UAExd0WeeN4B8kGH6F9KyPf8qpU%2BOJs9pBhoyUYYhCjBOEeX9lPOElafrGD3%2BtoWhwuwrCmApXI17d0RLtRvU%2BCrqh9EGNjaJbI0nee8hffpjzycJcyPL2XGSyqEzgco1Eg%2FbkasyQiEaFBwuVFErpe7%2FYOSRMTGLNDDtWsyNWePbLE%2B8aazpYTEzJr4R%2FkpPuTEK783WGFoeMTrr9idS%2B0r0xDUD4OO8lfJTuZVegya4DxUjLyOd6z5RADvgf641pWTnJUfunvRCtjLp3EJjeFnowR%2Fw6Zl6Fj91wdpF4Bep2SR9sjCj1iE0NukEht3MjX16BPM8z%2FlQgJer3a7Afy0IchCNNIcsbSKfWX0kCrbWTNTCqoij2hntH%2FKFaK%2FoU5rVZpnO9RljeCpeyCkZv0GA5nzb2FIVgykWOK%2FuZu8eIOx67PilQWxWPhtQBUta70aDEutPSLYmmUIEFzYkRGmZQ%2BHyLg8HV5uoMGR%2FlzAc9XFW90JSHh1VpBbmwo%2F8%2FoKwqrVHLTsXB0Gex6AN955vT%2FQXUEmizlGcYlBr7uy50XXEQSD23yuvGvGR1%2FIC0kxX3kL2%2FacndzfMC5LkMvfI6ZZK%2FGdGmDja34uFaPRA7Cj%2BjfX6%2BS6dVYdOUIjf6V2nl51l9ZVIa7hnfiwKs4y6U00MLXqFSprUAA8o4G4uFB4mXS49NSwETLn%2Fvht%2BXGnD1EuAua%2B9ZXDfin5mfFtxdXmCdLiPFMLjL968GOrEBelYk9P7fi8dB%2Bk%2FJ8jG4zxUpYN6hPXZyefmNx70efNx32eHkO51PlpRLDWrITvtmJPLnQjeBY9bgp4Vm7Mvbcks%2B1xXUERKUEToG%2Bbr3W38ZsEQySgwh24Em2FTmDSWwV17I4TlDEK8Lq8ocTUcD29KCatovSMduuwACN6gHYGqWR7f7hqTO0Z1%2BIQ9oLu%2Bnx2G1d2blPAuNvsAyejXguwLwoqvGYLv7sGo9kbGMOYM1&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20240322T212148Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTYZHNUN4VU%2F20240322%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=25d56f591e4757c5a9c07d41b14eb73545795dff49560f2bbf512921281eaac3&hash=0963d1df160689aed69e7e71aa1264adf7692f2a8a7634e279a201d04407554a&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=016516849090158U&tid=spdf-8298f423-c4aa-4d8f-ba32-1a95f359adce&sid=72574d30490f30451a984d312b7792ba61c0gxrqa&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=0f115c585c5751070053&rr=86893bc55ee8e178&cc=us)
  * [CZT vs FFT:Flexibility vs Speed](https://www.osti.gov/servlets/purl/816417)
