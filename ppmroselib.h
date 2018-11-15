// ppmroselib.h
// Copyright (C) 2013 Michael Rose

// ============================================================================
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// ============================================================================

// Include files; no other libraries needed than 'libgpp':

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>
#include <math.h>
#include <new>

#define VERSION "1.4"

// ----------------------------------------------------------------------------

// Example to compile successfully under Visual C++ with
// cl -GX ppmroselib.cpp <prgname>.cpp

#ifdef _WIN32
  typedef unsigned __int8   u_int8_t;
  typedef unsigned __int16  u_int16_t;
  typedef          __int32    int32_t;
  typedef unsigned __int32  u_int32_t;
  inline double round(double z) { return floor(z+0.5); }
  #include <io.h>
  #include <fcntl.h>
  #define binmode(fh) _setmode(_fileno(fh),_O_BINARY)
#else
  #define binmode(fh)
#endif

// ----------------------------------------------------------------------------

// Global types:

typedef unsigned long ulong;
typedef u_int16_t     pixel_t;   // 16 bit RGB-values for pixel.
typedef u_int16_t     color_t;   // Value type for color       grid.
typedef int32_t       deform_t;  // Value type for deformation grid.

// ============================================================================
// General common functions.
// ============================================================================

// Error messages and normal messages (both on error channel):

void Error(const char *format, ...);
void Print(const char *format, ...);

// ============================================================================
// Unified handling of parameters/options.
// ============================================================================

// Options are defined via 'AddFlag', 'AddString', 'AddInt' and
// 'AddDouble' in the overloaded method 'Define()'. Additionally
// some concluding lines for the usage message can be supplied
// by 'ExtraUsage()'. The defining 'mode' has the following meaning:
//
// 0:  Initialize the option.
// 1:  Print option for debugging purpose.
// 2:  Calculate formatting tabulators for option list of usage method.
// 3:  Print option line for usage method.
// 4:  Read option from current command line argument.

class Parameter
{
  // Definition mode and suppression of normal messages:
  int mode, quiet;

  // Internal variables for option parsing:
  int  argi, argii, argc;
  char buf[256], *arg, **argv;

  void PrintOption();
  int  GetOption(const char *option, int extra);

protected:

  virtual void Define() { }
  virtual void Check()  { }

  int  AddFlag(const char *option, int& param,
               const char *show, const char *meaning);
  int  AddString(const char *option, const char *& param, const char *value,
                 const char *show, const char *meaning);
  int  AddInt(const char *option, int& param, int radix,
              int value, int low, int high,
              const char *show, const char *meaning);
  int  AddDouble(const char *option, double& param,
                 double value, double low, double high,
                 const char *show, const char *meaning);
  void ExtraUsage(const char *msg);

public:

  // Program name:
  const char *prgName;

  // Additional lonely input name:
  const char *inpName;

  Parameter(const char *prgName) { this->prgName = prgName; }

  void Debug();
  void Usage(int exitCode);
  void operator()(int argcPrg, char *argvPrg[]);
};

// ============================================================================
// Basic handling for raw PPM-images.
// ============================================================================

struct Image
{
  long     width,   // Image width.
           height,  // Image height.
           depth,   // Maximum color value per channel.
           size;    // 3 * width * height, each pixel has RGB-tripel.
  pixel_t *pixel;   // 3 * width * height RGB values (16 bit).

  Image()  { pixel = 0;  Reset(); }
  ~Image() {             Reset(); }

  void Reset(long width=0, long height=0, long depth=255);
  void Read(const char *filename, FILE *fh=0);
  void Write(const char *filename, FILE *fh, const char *creator);
};

// ============================================================================
// Constant/linear/cubic interpolation.
// ============================================================================

// Basic usage for interpolation 'value = f(x,y)' in 2D-array 'raw'
// with incremental offsets 'ox' and 'oy' in the corresponding dimensions
// ('c0' is a flag to force linear interpolation):
//
// Interpolation2 inter(c0, mx, ox, my, oy);
// <loop>
//   inter.cx.Set(x);  with 0.0 <= x <= mx
//   inter.cy.Set(y);  with 0.0 <= y <= my
//   value = inter(raw);

// 'f[4]' contains up to four interpolation coefficients, depending on
// the fraction between the aequidistant grid points. 'm' is the maximum
// possible grid index (between 0 and 'm'), 'o' the row/column offset
// to adjust the grid pointer for adjacent rows/columns and 'n' denotes
// the valid coefficients 'f[0 ... n]'.
// 'o2' is another offset suitable to create system matrices.
// The cubic interpolation is valid for 'm >= 2'.
// The interpolation changes automatically to linear order, if 'm < 2'
// and to constant order if 'm == 0'.

struct Interpolation1 {
  double f[4];
  long   order, m, mm, o, o2, i, n;
  void SetSize(int c0, long m) {
    this->m = m;
    order   = m < 2 ? m : (c0 ? 1 : 3);
    mm      = m ? m : 1; }
  long Offset() { return i * o; }
  Interpolation1(int c0, long m, long o) { SetSize(c0, m);  this->o = o; }
  Interpolation1(int c0, long m, long o, long o2) {
    SetSize(c0, m);  this->o = o;  this->o2 = o2; }
  void Set(double t1);
};

// ----------------------------------------------------------------------------

// Twodimensional constant/linear/cubic interpolation:

template<class ityp>
struct Interpolation2 {
  Interpolation1 cx,cy;  // Onedimensional interpolation along x- and y-axis.
  long           oo;
  long Offset() { return cx.Offset() + cy.Offset(); }
  Interpolation2(int c0, long mx, long ox, long my, long oy) :
    cx(c0, mx, ox), cy(c0, my, oy) { }
  Interpolation2(int c0, long mx, long ox, long my, long oy, long oo) :
    cx(c0, mx, ox, oo*ox), cy(c0, my, oy, oo*oy) { this->oo = oo; }
  double operator()(ityp *raw);
  void AssembleRhs(double scale, double *rhs);
  void AssembleMat(double *mat);
};

// ============================================================================
// Solver for symmetric positive definite systems with limited bandwidth.
// ============================================================================

// System modul to solve linear equation with band matrix.
// Matrix is of size 'n x n' with lower and upper bandwidth 'b'.
// The offset between valid entries on adjacent columns on the
// same row is 'o = 2*b', the adjacent diagonal elements have offset 'o + 1':

class BandedSystem
{
  long    n, b, o;
  double *sol, *rhs, *mat;

public:

  double *GetSol() { return sol; }
  double *GetRhs() { return rhs; }
  double *GetMat() { return mat; }

  BandedSystem()  { sol = 0;  Reset(); }
  ~BandedSystem() {           Reset(); }

  void Reset(long n=0, long b=0);
  void ClearRhs();
  void ClearMat();
  void Cholesky();
  void Solve();
};

// ============================================================================
// Main routine.
// ============================================================================

// A global 'Main'-routine is required:

extern void Main();
extern int  Handler(int argc, char *argv[], Parameter& param);
