// This file is part of VisionTools.
//
// Copyright 2011 Hauke Strasdat (Imperial College London)
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights  to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.

#ifndef SAMPLE_H_
#define SAMPLE_H_

#include <tr1/random>

namespace vk {

class Sample
{
public:
  static int
  uniform                    (int from, int to);
  static double
  uniform                    ();
  static double
  gaussian                   (double sigma);

  static tr1::ranlux_base_01 gen_real;
  static tr1::mt19937 gen_int;
};

tr1::ranlux_base_01 Sample::gen_real;
tr1::mt19937 Sample::gen_int;

int Sample
::uniform(int from, int to)
{
  tr1::uniform_int<int> unif(from, to);
  int sample = unif(gen_int);
  return sample;
}

double Sample
::uniform()
{
  tr1::uniform_real<double> unif(0.0, 1.0);
  double sample = unif(gen_real);
  return sample;
}

double Sample
::gaussian(double sigma)
{
  tr1::normal_distribution<double> gauss(0.0, sigma);
  double sample = gauss(gen_real);
  return sample;
}

}

#endif /* SAMPLE_H_ */
