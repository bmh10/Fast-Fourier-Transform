/*
 * Image.cxx
 *
 *  Created on: Feb 18, 2011
 *      Author: dr
 */

#include "Image.h"
#include <stdlib.h>

template <class PixelType> Image<PixelType>::Image()
{
  m_max  = 0;
  m_cols = m_rows = 0;
}

template <class PixelType> Image<PixelType>::Image(const Image& ref):vector<PixelType>(ref)
{
  m_max  = ref.m_max;
  m_cols = ref.m_cols;
  m_rows = ref.m_rows;
}

template <class PixelType> Image<PixelType>::Image(int rows, int cols)
{
  init(rows, cols);
}

template <class PixelType> Image<PixelType>::~Image()
{}

template <class PixelType> void Image<PixelType>::init (int rows, int cols)
{
  m_max = 0;
  m_cols = cols;
  m_rows = rows;
  vector<PixelType>::resize(rows*cols);
}

template <> int Image<complexd>::read(string fileName)
{
  return 0;
}

template <class PixelType> int Image<PixelType>::read(string fileName)
{
  ifstream in (fileName.c_str());
  if(!in.is_open()) return 0;

  char buffer[MAX_PGM_LINE_NUMBER];
  in.getline(buffer, MAX_PGM_LINE_NUMBER);

  string magic = string(buffer);
  if (magic == "P5") { //Normal PGM

    // Get rid of comments
    for (in.getline(buffer, MAX_PGM_LINE_NUMBER); buffer[0] == '#'; in.getline(buffer, MAX_PGM_LINE_NUMBER));

    // Read column and row
    stringstream sizeStream(buffer);
    sizeStream >> m_cols;
    sizeStream >> m_rows;
    init(m_rows, m_cols);

    // Get rid of comments
    for (in.getline(buffer, MAX_PGM_LINE_NUMBER); buffer[0] == '#'; in.getline(buffer, MAX_PGM_LINE_NUMBER));

    // Read max value
    stringstream maxStream(buffer);
    maxStream >> m_max;
    int pixByte = (m_max < 256 ) ? 1 : 2;

    unsigned int i, pixel;
    unsigned char byte[1];

    for (i = 0; i < vector<PixelType>::size(); i++) {
      pixel =0;
      if (pixByte == 2) {
        in.read((char*)byte,1);
        pixel += byte[0]<<8;
      }
      in.read((char*)byte,1);
      pixel += byte[0];

      (*this)[i] = pixel;
    }
  } else if (magic == "P2") { //plain PGM
    cout << "Raw PGM reading has not been implemented." <<endl;
    return 0;

  } else {
    cout << "Incorrect PGM magic number" <<endl;
    return 0;
  }
  return 1;
}

template <> int Image<complexd>::write(string fileName)
{
  return 0;
}

template <class PixelType> int Image<PixelType>::write(string fileName)
{
  double min = 0;
  unsigned int pixel;
  unsigned char byte[1];


  ofstream out(fileName.c_str());
  if(!out.is_open()) return 0;

  // Write header
  out<<"P5"<<endl;
  out<<"# File created at ";
  char *hostName = getenv("HOST");
  if (hostName) out << hostName;
  char *user = getenv("USER");
  if (user) out <<" by "<< user;
  time_t rawtime;
  tm * timeinfo;
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  out << " at " << asctime (timeinfo);

  // Write rows and columns
  out << m_cols << " " << m_rows << endl;

  // Normalize image intensities
  m_max = 0;
  unsigned int i; for(i = 0; i < vector<PixelType>::size(); i++) {
    if(m_max < (*this)[i]) m_max = (*this)[i];
    if(min > (*this)[i]) min = (*this)[i];
  }
  if (min < 0) {
    for(i = 0; i < vector<PixelType>::size(); i++) (*this)[i] -= min;
    m_max -=min;
    min =0;
  }
  if (m_max > 255) {
    double scale = 255.0/m_max;
    for(i = 0; i < vector<PixelType>::size(); i++) {
      (*this)[i] *= scale;
    }
    m_max = 255;
  }

  // Write maximum value
  out << int(m_max) << endl;

  // Write pixels
  for (i = 0; i < vector<PixelType>::size(); i++) {
    pixel = int(round(min+(*this)[i]));
    byte[0] = (unsigned char)(pixel);
    out.write((char*)byte, 1);
  }

  out.close();
  return 1;
}

template <class PixelType> vector <PixelType> Image<PixelType>::extractCol(int colNum)const
{
  unsigned int i;
  vector <PixelType> ret (m_rows);

  for(i=0; i<m_rows; i++) ret[i] = (*this)[i*m_cols + colNum];
  return ret;
}

template <class PixelType> vector <PixelType> Image<PixelType>::extractRow(int rowNum)const
{
  unsigned int i;
  vector <PixelType> ret (m_cols);

  for(i=0; i<m_cols; i++) ret[i] = (*this)[rowNum * m_cols + i];
  return ret;
}

template <class PixelType> void Image<PixelType>::fillCol(int colNum, vector <PixelType> in)
{
  unsigned int i;
  for(i=0; i<m_rows; i++) (*this)[i* m_cols + colNum]= in[i];
}

template <class PixelType> void Image<PixelType>::fillRow(int rowNum, vector <PixelType> in)
{
  unsigned int i;
  for(i=0; i<m_cols; i++) (*this)[rowNum * m_cols + i]= in[i];
}

template <class PixelType> PixelType Image<PixelType>::getData(unsigned int row, unsigned int col) const
{
  return (*this)[row * m_cols + col];
}

template <class PixelType> void Image<PixelType>::setData(unsigned int row, unsigned int col, PixelType val)
{
  if (row<m_rows && col<m_cols && row >=0 && col>=0) (*this)[row*m_cols + col] = val;
}

// Force compiler to instantiate code for Image<double>
template class Image<double>;

// Force compiler to instantiate code for Image<complexd>
template class Image<complexd>;

