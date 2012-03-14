/*
 * main.cxx
 *
 *  Created on: Feb 18, 2011
 *      Author: dr
 */

// Global includes
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <vector>
#include <string.h>

// Local includes
#include "Image.h"
#include "FFT.h"

using namespace std;

template <class type>type stringTo (const char in[])
{
  std::istringstream tmp(in);
  type out;
  tmp >> out;
  return out;
}

template <class type>type stringTo (string& in)
{
  std::istringstream tmp(in);
  type out;
  tmp >> out;
  return out;
}

void usage()
{
  cout << "Usage: FastFourierTransform fileIn fileOut [options]" << endl;
  cout << "Option -f [1/2/3]: filter frequencies with:" << endl;
  cout << "1: low pass filter" << endl;
  cout << "2: high pass filter" << endl;
  cout << "3: band pass filter" << endl;
  cout << "Option -pmax: percentage of maximum frequency above which is ignored in LP and BP filters" << endl;
  cout << "Option -pmin: percentage of maximum frequency below which is ignored in HP and BP filters" << endl;
}

int main (int argc, char *argv[])
{
  if (argc < 3) {
  	usage();
    exit(-1);
  }

  // Default parameters
  int filter  = 0;
  double maxp = 0.1;
  double minp = 0.01;

  // Input filename
  string input_filename = argv[1];
  argv++;
  argc--;

  // Output filename
  string output_filename = argv[1];
  argv++;
  argc--;

  // Parse remaining options
  bool ok;
  while (argc > 1) {
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-f") == 0)) {
      argc--;
      argv++;
      filter = stringTo<int>(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-pmax") == 0)) {
      argc--;
      argv++;
      maxp = stringTo<double>(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-pmin") == 0)) {
      argc--;
      argv++;
      minp = stringTo<double>(argv[1]);
      argc--;
      argv++;
      ok = true;
    }
    if (ok == false) {
      cerr << "Can not parse argument " << argv[1] << endl;
      usage();
      exit(-1);
    }
  }

  Image<double> input, output;
  Image<complexd> tmp;

  // Read image
  cout << "Reading image ... " << endl;
  if (!input.read(input_filename)) {
    cout << "Cannot read image from file " << argv[1] << endl;
    exit(1);
  }

  cout << "Running FFT (forward) ..." << endl;
  FFT encoder;

  tmp = encoder.ForwardTransform(input);

  if (filter) {
    cout << "Filtering ..."; cout.flush();
    encoder.Filter(tmp,filter,maxp,minp);
  }

  cout << "Running FFT (reverse) ..." << endl;
  output = encoder.ReverseTransform(tmp);

  // Write output
  cout << "Writing image..." << endl;
  if (!output.write(output_filename)) {
    cout << "Can't write to file " << output_filename << endl;
    exit(-1);
  }

  return 0;
}

