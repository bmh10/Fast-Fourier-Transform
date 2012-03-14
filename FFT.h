/*
 * FFT.h
 *
 *  Created on: Feb 18, 2011
 *      Author: dr
 */

#ifndef FFT_H_
#define FFT_H_

// Global includes
#include <string>
#include <math.h>
#include <vector>
#include <complex>

// Local includes
#include "Image.h"

using namespace std;

class FFT
{

private:

  // Run 2D FFT in the direction specified
  void Run_2D(Image<complexd> &image, int dir);

  // Run 1D FFT in the direction specified (uses a recursive implementation)
  vector<complexd> Run_1D_Recursive(const vector<complexd> &v, int n, int dir);

  // Run 1D FFT in the direction specified (uses an iterative implementation)
  vector<complexd> Run_1D_Iterative(const vector<complexd> &v, int n, int dir);

  // Compute magnitude of complex image
  Image<double> Magnitude(Image<complexd> c);
  
  // Computes distances from corner to a point in image
  double DistFromCorner(int cx, int cy, int px, int py);
  
  // Reverses bits in a vector
  vector<complexd> BitReverseCopy(const vector<complexd> &v, int n);
  
  // Applies filter to individual pixel and returns true if pixel should
  // be set to zero depending on filter type being applied
  bool pixel_filter(int i, int j, int rows, int cols, double maxcutd, 
					   double mincutd, int type);

public:

  // Constructor
  FFT();

  // Filter an image using low pass, high pass and band pass filters. The
  // parameters max and min specify the frequencies for the filters
  void Filter(Image<complexd> &image, int type, double max, double min);
  
  // Run forward 2D FFT on a 2D image of type double. Returns an image of
  // complex doubles
  Image<complexd> ForwardTransform(const Image<double> image);

  // Run reverse 2D FFT on a 2D image of type complex double. Returns an image
  // doubles
  Image<double> ReverseTransform(const Image<complexd> image);
  
};

// Constructor with not much to do
inline FFT::FFT()
{}

#endif /* FFT_H_ */
