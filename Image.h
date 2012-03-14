/*
 * Image.h
 *
 *  Created on: Feb 18, 2011
 *      Author: dr
 */

#ifndef IMAGE_H_

#define IMAGE_H_

#include <string>
#include <sstream>
#include <math.h>
#include <vector>
#include <fstream>
#include <time.h>
#include <iostream>
#include <complex>

using namespace std;

const int MAX_PGM_LINE_NUMBER = 1024;

typedef complex<double> complexd;

template <class PixelType> class Image : public vector<PixelType>
{

protected:

	// Number of rows and columns in image
  unsigned int m_cols, m_rows;

  // Maximum value in image
  PixelType m_max;

public:

  // Default constructor
  Image();

  // Copy constructor
  Image(const Image& ref);

  // Construct image with no. of rows and columns
  Image(int rows, int cols);

  // Destructor for image
  virtual ~Image();

	// Resize image to no. of rows and columns
  void init (int rows, int cols);

  // Read image from filename
  int read(string fileName);

  // Write image to filename
  int write(string fileName);

  // Extract column from image
  vector <PixelType> extractCol(int colNum) const;

  // Extract row from image
  vector <PixelType> extractRow(int rowNum) const;

  // Fill column in image
  void fillCol(int colNum, vector <PixelType> in);

  // Fill row in image
  void fillRow(int rowNum, vector <PixelType> in);

  // Return number of rows
  unsigned int numRows() const;

  // Return number of columns
  unsigned int numCols() const;

  // Return data at row, col
  PixelType getData(unsigned int row, unsigned int col) const;

  // Set data at row, col
  void setData(unsigned int row, unsigned int col, PixelType val);

};

template <class PixelType> inline unsigned int Image<PixelType>::numRows() const
{
  return m_rows;
}
template <class PixelType> inline unsigned int Image<PixelType>::numCols() const
{
  return m_cols;
}


#endif /* IMAGE_H_ */
