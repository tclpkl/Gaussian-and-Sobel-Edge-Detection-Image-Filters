#include <iostream>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <cstdlib>
#include "bmplib.h"

using namespace std;

//============================Add function prototypes here======================

void dummy(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB]);

void convolve(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB], 
        int N, double kernel[][11]);

void sobel(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB]);

void gaussian(double k[][11], int N, double sigma);

void gaussian_filter(unsigned char out[][SIZE][3], unsigned char 
     in[][SIZE][3], int N, double sigma);

void unsharp(unsigned char out[][SIZE][3], unsigned char in[][SIZE][3], 
     int N, double sigma, double alpha);

//============================Do not change code in main()======================

#ifndef AUTOTEST

int main(int argc, char* argv[])
{
   //First check argc
  if(argc < 3)
    {
      //we need at least ./filter <input file> <filter name> to continue
      cout << "usage: ./filter <input file> <filter name> <filter parameters>";
      cout << " <output file name>" << endl;
      return -1;
    }
   //then check to see if we can open the input file
   unsigned char input[SIZE][SIZE][RGB];
   unsigned char output[SIZE][SIZE][RGB];
   char* outfile;
   int N;
   double sigma, alpha;
   //double kernel[11][11];

   // read file contents into input array
   int status = readRGBBMP(argv[1], input); 
   if(status != 0)
   {
      cout << "unable to open " << argv[1] << " for input." << endl;
      return -1;
   }
   //Input file is good, now look at next argument
   if( strcmp("sobel", argv[2]) == 0)
   {
     sobel(output, input);
     outfile = argv[3];
   }
   else if( strcmp("blur", argv[2]) == 0)
   {
     if(argc < 6)
       {
	 cout << "not enough arguments for blur." << endl;
	 return -1;
       }
     N = atoi(argv[3]);
     sigma = atof(argv[4]);
     outfile = argv[5];
     gaussian_filter(output, input, N, sigma);
   }
   else if( strcmp("unsharp", argv[2]) == 0)
   {
     if(argc < 7)
       {
	 cout << "not enough arguments for unsharp." << endl;
	 return -1;
       }
     N = atoi(argv[3]);
     sigma = atof(argv[4]);
     alpha = atof(argv[5]);
     outfile = argv[6];
     unsharp(output, input, N, sigma, alpha);
   }
   else if( strcmp("dummy", argv[2]) == 0)
   {
     //do dummy
     dummy(output, input);
     outfile = argv[3];
   }
   else
   {
      cout << "unknown filter type." << endl;
      return -1;
   }

   if(writeRGBBMP(outfile, output) != 0)
   {
      cout << "error writing file " << argv[3] << endl;
   }

}   

#endif 

//=========================End Do not change code in main()=====================


// Creates an identity kernel (dummy kernel) that will simply
// copy input to output image and applies it via convolve()
//
// ** This function is complete and need not be changed.
// Use this as an example of how to create a kernel array, fill it in
// appropriately and then use it in a call to convolve.
void dummy(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB])
{
   double k[11][11];
   for (int i = 0; i < 3; i++)
   {
      for(int j = 0; j < 3; j++)
      {
         k[i][j] = 0;
      }
   }
   k[1][1] = 1;
   convolve(out, in, 3, k);
}


// Convolves an input image with an NxN kernel to produce the output kernel
// You will need to complete this function by following the 
//  instructions in the comments
void convolve(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB], 
	      int N, double kernel[][11])
{
 
   int padded[SIZE+10][SIZE+10][RGB];  // Use for input image with appropriate 
                                       // padding
   int temp[SIZE][SIZE][RGB];          // Use for the unclamped output pixel 
                                       // values then copy from temp to out, 
                                       // applying clamping 

   //first set all of padded to 0 (black)
  for (int i = 0; i < SIZE + 10; i += 1) {
    for (int j = 0; j < SIZE + 10; j += 1) {
      for (int c = 0; c < 3; c += 1) {
        padded[i][j][c] = 0;
      }
    }
  }

   //now copy input into padding to appropriate locations
  for (int i = 0; i < SIZE; i += 1) {
    for (int j = 0; j < SIZE; j += 1) {
      for (int c = 0; c < 3; c += 1) {
        padded[i + N/2][j + N/2][c] = in[i][j][c];
      }
    }
  }

   //initialize temp pixels to 0 (black)
  for (int i = 0; i < SIZE; i += 1) {
    for (int j = 0; j < SIZE; j += 1){
      for (int c = 0; c < 3; c += 1) {
        temp[i][j][c] = 0;
      }
    }
  }

  //now perform convolve (using convolution equation on each pixel of the 
  // actual image) placing the results in temp (i.e. unclamped result)
  //Here we give you the structure of the convolve for-loops, you need
  //to figure out the loop limits
  for (int y = 0; y < SIZE ; y += 1) {
    for (int x = 0; x < SIZE; x += 1) {
      for (int k = 0; k < RGB; k += 1) {
         for (int i = -N/2; i <= N/2; i += 1) {
            for (int j = -N/2; j <= N/2; j += 1) {
                  temp[x][y][k] += padded[x + N/2 + i][y + N/2 + j][k]
                    *kernel[N/2 + i][N/2 + j];
            }
         }
      }
    }
  }
  
   //now clamp and copy to output
   // You may need to cast to avoid warnings from the compiler:
   // (i.e. out[i][j][k] = (unsigned char) temp[i][j][k];)
  for (int i = 0; i < SIZE; i += 1) {
    for (int j = 0; j < SIZE; j += 1) {
      for (int c = 0; c < 3; c += 1) {
        if (temp[i][j][c] > 255) {
          temp[i][j][c] = 255;
        }
        if (temp[i][j][c] < 0) {
          temp[i][j][c] = 0;
        }
        out[i][j][c] = (unsigned char) temp[i][j][c];
      }
    }
  }

}

// You will need to complete this function by following the 
//  instructions in the comments
void sobel(unsigned char out[][SIZE][RGB], unsigned char in[][SIZE][RGB])
{
   double k[11][11];
   double s_h1[3][3] = { {-1, 0, 1}, 
                         {-2, 0, 2}, 
                         {-1, 0, 1} };
   double s_h2[3][3] = { {1, 0, -1}, 
                         {2, 0, -2}, 
                         {1, 0, -1} };
   
   unsigned char h1_sobel[SIZE][SIZE][RGB]; //hold intemediate images
   unsigned char h2_sobel[SIZE][SIZE][RGB]; 

   for (int i = 0; i < 11; i++)
   {
      for(int j = 0; j < 11; j++)
      {
         k[i][j] = 0;
      }
   }


   // Copy in 1st 3x3 horizontal sobel kernel (i.e. copy s_h1 into k)
  for (int i = 0; i < 3; i += 1) {
    for (int j = 0; j < 3; j += 1) {
      k[i][j] = s_h1[i][j];
    }
  }

   // Call convolve to apply horizontal sobel kernel with result in h1_sobel
  convolve(h1_sobel, in, 3, k);


   // Copy in 2nd 3x3 horizontal sobel kernel (i.e. copy s_h2 into k)
  for (int i = 0; i < 3; i += 1) {
    for (int j = 0; j < 3; j += 1) {
      k[i][j] = s_h2[i][j];
    }
  } 

   // Call convolve to apply horizontal sobel kernel with result in h2_sobel
  convolve(h2_sobel, in, 3, k);


   // Add the two results (applying clamping) to produce the final output 
  for (int i = 0; i < SIZE; i += 1) {
    for (int j = 0; j < SIZE; j += 1) {
      for (int c = 0; c < 3; c += 1) {
        out[i][j][c] = h1_sobel[i][j][c] + h2_sobel[i][j][c];
        if (out[i][j][c] > 255) {
          out[i][j][c] = 255;
        }
        if (out[i][j][c] < 0) {
          out[i][j][c] = 0;
        }
        out[i][j][c] = (unsigned char)out[i][j][c];
      }
    }
  }

}


// Add the rest of your functions here (i.e. gaussian, gaussian_filter, unsharp)

void gaussian(double k[][11], int N, double sigma) {
  
  //setup array for y, x value of each pixel in gaussian kernel
  double y_of_gaussian[11][11];
  double x_of_gaussian[11][11];
  
  //storing total of kernel values for normalization later
  double kernel_total = 0;
  
  //alter value of each y,x index to correspond to gaussian indexing values 
  for (int i = 0; i < N; i += 1) {
    for (int j = 0; j < N; j += 1) {
      y_of_gaussian[i][j] = i - N/2;
      x_of_gaussian[i][j] = j - N/2;
    }
  }
  
  //applying formula to calculate raw values, then copying to kernel
  //adding each calculated raw value to the kernel total 
  for (int i = 0; i < N; i += 1) {
    for (int j = 0; j < N; j += 1) {
      k[i][j] = exp(-((pow(y_of_gaussian[i][j], 2) / (2 * pow(sigma, 2))) 
                      + (pow(x_of_gaussian[i][j], 2) / (2 * pow(sigma, 2)))));
      kernel_total += k[i][j];
    }
  }
  
  //dividing each pixel of kernel by raw sum for normalized result
  for (int i = 0; i < N; i += 1) {
    for (int j = 0; j < N; j += 1) {
      k[i][j] = k[i][j] / kernel_total;
    }
  }
  
}

void gaussian_filter(unsigned char out[][SIZE][3], unsigned char 
     in[][SIZE][3], int N, double sigma) {
     
  //creating an 11 by 11 kernel to use
  double k[11][11];
     
  //initializing all the values in the kernel to 0
  for (int i = 0; i < 11; i += 1) {
    for (int j = 0; j < 11; j += 1) {
      k[i][j] = 0;
    }
  }
  
  //creating kernel values based off N and sigma
  gaussian(k, N, sigma);
  
  //applying gausian filter to input image and creating output image
  convolve(out, in, N, k);

}

void unsharp(unsigned char out[][SIZE][3], unsigned char in[][SIZE][3], 
     int N, double sigma, double alpha) {
  
  //creating a temporary array to hold blurred pixels
  unsigned char blurred[SIZE][SIZE][RGB];
  
  //creating a temporary array to hold detailed pixels
  int detailed[SIZE][SIZE][RGB];
  
  //blurring the input array and storing in new blurred array
  gaussian_filter(blurred, in, N, sigma);
  
  //iterating through each pixel creating detailed map 
  for (int i = 0; i < SIZE; i += 1) {
    for (int j = 0; j < SIZE; j += 1) {
      for (int c = 0; c < 3; c += 1) {
        detailed[i][j][c] = in[i][j][c] - blurred[i][j][c];
      }
    }
  }
  
  //iterating through each pixel calculating final sharpen pixel
  //value, clamping the value, and then copying the pixel value to 
  //the final output pixel
  for (int i = 0; i < SIZE; i += 1) {
    for (int j = 0; j < SIZE; j += 1) {
      for (int c = 0; c < 3; c += 1) {
        if ((in[i][j][c] + alpha * detailed[i][j][c]) > 255) {
          out[i][j][c] = 255;
        }
        else if ((in[i][j][c] + alpha * detailed[i][j][c]) < 0) {
          out[i][j][c] = 0;
        }
        else {
          out[i][j][c] = (unsigned char)(in[i][j][c] + alpha * 
            detailed[i][j][c]);
        }
      }
    }
  }
  
}