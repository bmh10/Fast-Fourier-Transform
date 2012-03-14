/*
 * FFT.cxx - bmh10
 */

// Global includes
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <complex>
#include <math.h>

// Local includes
#include "FFT.h"

using namespace std;

/* 
 * Computes the magnitude of a complex image and required to
 * convert an image of complex doubles back into an image of doubles.
 */
Image<double> FFT::Magnitude(Image<complexd> input)
{
	unsigned int i, j, rows, cols;
	rows = input.numRows();
	cols = input.numCols();
	Image<double> ret = Image<double>(rows, cols);
	
	complexd curr;
	
	for (i = 0; i < rows; ++i) {
  	    for (j = 0; j < cols; ++j)
  		{
  		    // Get data from complex image pixel, calculate magnitude, then
  		    // set pixel in double image.
  			curr = input.getData(i, j);			 
			ret.setData(i, j, (double) sqrt( pow( real(curr), 2) 
											     + pow( imag(curr),2 )));
		}
	}
	   
	return ret;
}

double FFT::DistFromCorner(int cx, int cy, int px, int py)
{
    double x = (double) abs(px-cx);
    double y = (double) abs(py-cy);
    return sqrt(pow(x, 2) + pow(y, 2));
}

/* 
 * Applies filter to individual pixel and returns true if pixel should
 * be set to zero depending on filter type being applied.
 */
bool FFT::pixel_filter(int i, int j, int rows, int cols, double maxcutd, 
					   double mincutd, int type)
{
    double dist1, dist2, dist3, dist4;
	dist1 = DistFromCorner (0, 0, i, j);
    dist2 = DistFromCorner (rows, 0, i, j);
    dist3 = DistFromCorner (0, cols, i, j);
  	dist4 = DistFromCorner (rows, cols, i, j);
  	
   bool ret = false;
   
  	if (type != 2)
  	{
		// Set frequencies towards centre to zero.
		ret |= (dist1 > maxcutd
			&& dist2 > maxcutd
			&& dist3 > maxcutd
			&& dist4 > maxcutd);
	}
	if (type != 1)
	{
		// Set frequencies towards corners to zero.
		ret |= (dist1 < mincutd
			 || dist2 < mincutd
			 || dist3 < mincutd
			 || dist4 < mincutd);
	}
	return ret;
}



void FFT::Filter(Image<complexd> &image, int type, double maxcut, double mincut)
{
  unsigned int i, j, rows, cols;
  rows = image.numRows();
  cols = image.numCols();

  double maxdist = DistFromCorner(0, 0, rows/2, cols/2);
  double maxcutd = maxdist * maxcut;
  double mincutd = maxdist * mincut;
  
  if (type == 1)
  {
  	 cout << " low-pass ..." << endl;
  	
  	 // Sets all frequencies of the image above a certain frequency to 0.
  	 for (i = 0; i < rows; ++i) {
  		for (j = 0; j < cols; ++j)
  		{
  			if (pixel_filter(i, j, rows, cols, maxcutd, 0, type))
                image.setData(i, j, 0);
 		}
 	 }
  }
  
  else if (type == 2)
  {
  	 cout << " high-pass ..." << endl;
  	 
  	 // Sets all frequencies of the image below a certain frequency to 0.
	 for (i = 0; i < rows; ++i) {
  		for (j = 0; j < cols; ++j)
  		{
			 if (pixel_filter(i, j, rows, cols, 0, mincutd, type))
                 image.setData(i, j, 0);
 		}
 	 }
  }
  
  else if (type == 3)
  {
  	 cout << " band-pass ..." << endl;

     // Sets all freqs of the image above and below certain frequencies to 0.
	 for (i = 0; i < rows; ++i) {
  		for (j = 0; j < cols; ++j)
  		{
  			 if (pixel_filter(i, j, rows, cols, maxcutd, mincutd, type))
                image.setData(i, j, 0);
 		}
 	 }  
  }
}

Image<complexd> FFT::ForwardTransform(Image<double> image)
{
    int i, j;
    int rows = image.numRows();
    int cols = image.numCols();
    
    Image<complexd> res = Image<complexd>(rows, cols);
    
    for (i = 0; i < rows; ++i) {
        for (j = 0; j < cols; ++j)
            res.setData(i, j, complex<double>(image.getData(i, j), 0));
    }

    Run_2D(res, 1);
    return res;
}

Image<double> FFT::ReverseTransform(Image<complexd> image)
{
    Run_2D(image, -1);
    return Magnitude(image);
}

void FFT::Run_2D(Image<complexd> &image, int dir)
{
  // dir = +1 for forward transform
  // dir = -1 for reverse transform
  
  unsigned int rows, cols, i;
  rows = image.numRows();
  cols = image.numCols();

  // Apply 1D FFT to all rows.
  for (i = 0; i < rows; ++i)
      image.fillRow(i, Run_1D_Recursive(image.extractRow(i), cols, dir));

 // Apply 1D FFT to all columns.
 for (i = 0; i < cols; ++i)
     image.fillCol(i, Run_1D_Recursive(image.extractCol(i), rows, dir));

}

vector<complexd> FFT::Run_1D_Recursive(const vector<complexd> &v, int n, int dir)
{
  if (n == 1)
    return v;
    
  vector<complexd> a0, a1, y0, y1;
  complexd w, wn, tmp, tmp2, t;

  wn = complexd(cos(2 * M_PI / (double) n), dir * sin(2 * M_PI / (double) n));

  w = complexd(1, 0);
  
  for (int i = 0; i < n/2; ++i)
  {
    a0.push_back(v.at(2 * i));
    a1.push_back(v.at(2 * i+1));
  }
  
  y0 = Run_1D_Recursive(a0, n/2, dir);
  y1 = Run_1D_Recursive(a1, n/2, dir);
  
  vector<complexd> y (n);
  for (int k = 0; k < n/2; ++k)
  {
  	t = w * y1[k];
    tmp = y0[k] + t;
    tmp2 = y0[k] - t;
    
    if (dir == 1)
	{
    	y[k] = tmp;
    	y[k+n/2] = tmp2;
    }
    else
    {
    	y[k] = complexd (real(tmp)/2, imag(tmp)/2);
    	y[k+n/2] = complexd (real(tmp2)/2, imag(tmp2)/2);
    }
    	
    w *= wn;
  }
  
  return y;
}

/* Optional - Iterative 1D FFT implementation. */

vector<complexd> FFT::BitReverseCopy(const vector<complexd> &v, int n)
{
	vector<complexd> ret(n);
	int rev, val;
	
    for (int i = 0; i < n; ++i)
    {
    	rev = 0;
    	val = i;
    	for (int j = log(n)/log(2); j > 0; --j)
    	{
    	    rev <<= 1;
    	    rev |= val & 1;
    	    val >>= 1;
    	
    	}
    	ret[rev] = v[i];
    }
    
    return ret;
}

vector<complexd> FFT::Run_1D_Iterative(const vector<complexd> &a, int n, int dir)
{
    vector<complexd> A = BitReverseCopy(a, n);
	
	  int m;
	  complexd w, wm, t, u;

	  for (int s = 1; s <= log(n)/log(2); ++s)
	  {
		    m = pow (2, s);
		    wm = complexd(cos(2 * M_PI / (double) m), dir * sin(2 * M_PI / (double) m));
		
		    for (int k = 0; k < n; k += m)
		    {
		        w = complexd (1, 0);

		        for (int j = 0; j < m/2; ++j)
		        {
		            t = w * A[k+j+m/2];
		            u = A[k+j];
		            A[k+j] = u + t;
		            A[k+j+m/2] = u - t;
		            w = w * wm;
		        }
	      }
    }
    
    return A;
}
