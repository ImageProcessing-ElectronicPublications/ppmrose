// ppmunwarp.cc
// Copyright (C) 2013 Michael Rose

/* ============================================================================
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
============================================================================ */

// This program can be used to unwarp pictures from digital cameras
// with information obtained from a suitable calibration grid.
// The primary goal is to unwarp pictures from book pages or other material
// using a digital camera instead of a scanner. If the positions of the
// camera and the book pages as well as the camera settings (e.g. zoom,
// resolution) are not changed, this program can be used to automate the
// unwarp procedure.
// All pictures read and written by this program must be in raw PPM format.
// The program is testet with Linux but should also run on other platforms
// with only minor modifications. No special libraries are used.
// To compile use something like:
//
// $ g++ -o ppmunwarp ppmunwarp.cc
//
// A brief help text is printed, when wrong command syntax is used or simply
// by typing:
//
// $ g++ -h
//
// The user creates a regular grid with red points (the color can be
// changed with the '-pc' option) on a sheet of paper. These points should
// cover the picture area of interest. Then she puts this calibration grid
// at the same position as the pages to be unwarped and takes a digital photo.
// By using conversion tools like 'jpegtopnm' she creates a file,
// e.g. 'calibration.ppm'. From this file, the command
//
// $ ppmunwarp calibration.ppm > calibration.bin
//
// extracts the information needed to unwarp the book sheets. If some
// error occurs, probably the calibration points are not correctly detected.
// The option '-m check.ppm' writes an instructive picture about the
// calibration points, the program claims to detect. To adjust the detection
// process, the user can use the following options:
//
// -pc ff0000    Changes the hue value (via RGB) of the calibration points.
// -ph 85        Changes the allowable range around the hue value in percent.
// -ps 50        Changes the minimum saturation value in percent.
// -pv 50        Changes the minimum luminate value in percent.
//
// Increasing the percent values reduces the amount of detected points.
// The program accepts, if some few calibration points are not detected and
// also ignores some erroneous points. If in doubt always inspect the
// created 'check.ppm' file. All green points connected by blue lines
// are used as calibration input.
// The detected grid points define the area of the images to be unwarped
// in the following unwarping operation. A suitable resolution of the unwarped
// images is automatically determined. A small amount of extrapolation
// is also possible. If e.g. '-gm 2' is given, an additional boundary of the
// distance between two calibration points is added, but in that area
// there is probably more distortion in the final images. The option '-gs'
// determines the resolution of the cubic deformation grid (exported in
// 'calibration.bin') with respect to the calibration grid. With the default
// value '-gs 2' the deformation grid has twice the number of points as
// the calibration grid in each axis and should be fine in most circumstances.
// Technically the deformation grid is calculated by a moving least squares
// (MLS) algorithm from the calibration points and the later unwarping
// operation is done by cubic interpolation of the warped (x,y) positions
// on this deformation grid.
// If 'page01.ppm' is one of the images of the camera to be unwarped,
// this can be achieved by:
//
// $ ppmunwarp -d calibration.bin page01.ppm > page01_unwarped.ppm
//
// This program can handle 8 bit and 16 bit PPM files with channel depth
// values (maximal color value) from 1 to 65535. To specify an output color
// depth different from the input picture, use the option '-od'.
//
// Usually it is not crucial if the book page is not aligned with
// the calibration grid, because this leads only to a simple rotation
// in 'page01_unwarped.ppm', which can be corrected automatically by
// suitable programs like 'scantailor'.
//
// There is one drawback with this approach of unwarping:
// If a book is digitized page after page, the distance between camera
// and book pages is not constant. This continuously increases the errors
// made by the unwarping procedure. I have two principal ideas to prevent
// this annoying behaviour:
//
// * Keep the distant between book pages and camera constant during
//   the scan process by using a firmly mounted glass plate to which the
//   pages are attached one after one.
// * Take new calibration pictures after every 10 to 20 book pages and
//   organize Yourself well in the following postprocessing step.
//
// With the options '-gx 5.08' and '-gy 5.08' the number of calibration
// points per inch can be specified in x-direction and y-direction.
// The first options sets both values together for convenience, so
// the order of this pair of options is relevant. The default values
// are choosen for calibration points, which are seperated by 1/2 cm.
// These values are only used to calculate the PPI numbers of the
// unwarped pictures, which are displayed for information purpose.
// It is also possible to grab these values by suitable wrapper scripts.

// ============================================================================

// Include files; no other libraries needed than 'libgpp':

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>
#include <math.h>
#include <new>

#define VERSION "1.3"

// ----------------------------------------------------------------------------

// Example to compile successfully under Visual C++ 6.0 with
// cl -GX ppmunwarp.cpp

#ifdef _WIN32
  typedef unsigned __int8   u_int8_t;
  typedef unsigned __int16  u_int16_t;
  typedef __int32           int32_t;
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

typedef u_int16_t pixel_t;   // Type for RGB components of a pixel.
typedef int32_t   deform_t;  // Type for (x,y)-components of deformation grid.

typedef unsigned long ulong;

// ----------------------------------------------------------------------------

// Global variables:

const  int  nGlbBuf = 4096;
static char glbBuf[nGlbBuf];  // Used for error messages.

// ============================================================================

// Global parameters/options:

struct Parameter {

  // Program name:
  const char *prgName;

  // Suppression of normal messages:
  int quiet;

  // File names:
  const char *calibPicName;     // Input image with calibration points.
  const char *calibCheckName;   // Control image for detection process.
  const char *calibTextName;    // Calibration points in text format.
  const char *calibDeformName;  // Final deformation grid.
  const char *sourcePicName;    // Skewed input image.
  const char *destPicName;      // Unwarped output image.

  // Amount of additional area added around calibration grid mesh,
  // measured in units of distance between calibration points:
  int gridMargin;

  // Scale factor of deformation grid with respect to resolution
  // of the calibration grid:
  int gridScale;

  // Grid points per inch in x-direction:
  double gridPerInchX;

  // Grid points per inch in y-direction:
  double gridPerInchY;

  // The hue value of the calibration points (coded in 24 bit RGB):
  int pointColor;

  // The percentual width of the infeasible hue area
  // for calibration point detection:
  int pointRange;

  // The minimal saturation for calibration points:
  int pointSaturation;

  // The minimal luminance for calibration points:
  int pointThreshold;

  // Output color depth per channel:
  int outputDepth;

  Parameter();
  void Usage(int exitCode);
  void operator()(int argc, char *argv[]);
};

// ----------------------------------------------------------------------------

static Parameter param;  // Program parameters are global for ease of use.

// ----------------------------------------------------------------------------

// Standard values for global parameters:

Parameter::Parameter()
{
  prgName         = "ppmunwarp";
  quiet           = 0;
  calibPicName    = 0;
  calibCheckName  = 0;
  calibTextName   = 0;
  calibDeformName = 0;
  sourcePicName   = 0;
  destPicName     = 0;
  gridMargin      = 1;
  gridScale       = 2;
  gridPerInchX    = 5.08;
  gridPerInchY    = 5.08;
  pointColor      = 0xff0000;
  pointRange      = 85;
  pointSaturation = 50;
  pointThreshold  = 50;
  outputDepth     = 0;
}

// ----------------------------------------------------------------------------

// Usage/help information:

void Parameter::Usage(int exitCode)
{
  if (exitCode < 0) {
    fprintf(stdout,
            "%s " VERSION "\n\n"
            "Copyright (C) 2013 Michael Rose\n"
            "License GPLv3+: GNU GPL version 3 or later"
            " <http://gnu.org/licenses/gpl.html>\n"
            "This is free software: you are free to change"
            " and redistribute it.\n"
            "There is NO WARRANTY, to the extent permitted by law.\n\n",
            param.prgName);
    exit(0); }
  fprintf(stderr,
    "Usage: %s [options] [--] [inpname or stdin]\n\n"
    "Options:\n"
    "  --version     Print program version.\n"
    "  -h            Print program usage.\n"
    "  -q            Suppress normal program messages.\n"
    "  -c            Enforce calibration mode.\n"
    "  -cc (inpname) Set input PPM picture with calibration points.\n"
    "  -m  <name>    Set output PPM picture with checked calibration points.\n"
    "  -cp <name>    Set file name for textual calibration point list.\n"
    "  -d  (stdout)  Set file name for binary deformation grid.\n"
    "  -i  (inpname) Set input file name.\n"
    "  -o  (stdout)  Set ouput file name.\n"
    "  -gm (1)       Set additional margin around calibration points.\n"
    "  -gs (2)       Set scale factor for deformation grid.\n"
    "  -gx (5.08)    Set grid points per inch in x- and y-direction.\n"
    "  -gy (5.08)    Set grid points per inch in y-direction.\n"
    "  -pc (ff0000)  Point color for calibration point detection.\n"
    "  -ph (85)      Infeasible hue range  [%%] (point detection).\n"
    "  -ps (50)      Minimal    saturation [%%] (point detection).\n"
    "  -pv (50)      Minimal    luminance  [%%] (point detection).\n"
    "  -od (0)       Output color depth, zero means same as input depth.\n\n"
    "Simple calibration:  ppmunwarp [-m check.ppm] calib.ppm > deform.bin\n"
    "Simple unwarping:    ppmunwarp -d deform.bin skewed.ppm > unwarped.ppm\n"
    ,prgName);
  exit(exitCode);
}

// ----------------------------------------------------------------------------

// Sets global parameters/options from the command line:

void Parameter::operator()(int argc,char *argv[])
{
  const char  *inpName = 0;
  const char **doName = 0;
  int          doOpt = 1, doCalib = 0, *doNumber = 0, argi;
  double      *doReal = 0;
  char        *end;
  for (argi = 1; argi < argc; ++argi) {
    char *arg = argv[argi];
    if (!*arg) Usage(1);
    if (doName) { *doName = arg;  doName = 0; }
    else if (doNumber) {
      long  z = strtol(arg, &end, (doNumber == &pointColor ? 16 : 10));
      if (z < (doNumber == &gridMargin      ? 0        :
               doNumber == &outputDepth     ? 0        : 1)  ||
          z > (doNumber == &gridMargin      ? 4        :
               doNumber == &gridScale       ? 10       :
               doNumber == &pointColor      ? 0xffffff :
               doNumber == &outputDepth     ? 65535    : 100))  Usage(1);
      *doNumber = int(z);
      doNumber  = 0; }
    else if (doReal) {
      double z = strtod(arg, &end);
      if (z < 0.01 || z > 100.0)  Usage(1);
      if (doReal == &gridPerInchX)  gridPerInchY = z;
      *doReal = z;
      doReal  = 0; }
    else if (doOpt && *arg == '-') {
      if      (!strcmp(arg, "--version"))  Usage(-1);
      else if (!strcmp(arg, "--"))         doOpt = 0;
      else if (!strcmp(arg, "-h"))         Usage(0);
      else if (!strcmp(arg, "-q"))         quiet    = 1;
      else if (!strcmp(arg, "-c"))         doCalib  = 1;
      else if (!strcmp(arg, "-cc"))        doName   = &calibPicName;
      else if (!strcmp(arg, "-m"))         doName   = &calibCheckName;
      else if (!strcmp(arg, "-cp"))        doName   = &calibTextName;
      else if (!strcmp(arg, "-d"))         doName   = &calibDeformName;
      else if (!strcmp(arg, "-i"))         doName   = &sourcePicName;
      else if (!strcmp(arg, "-o"))         doName   = &destPicName;
      else if (!strcmp(arg, "-gm"))        doNumber = &gridMargin;
      else if (!strcmp(arg, "-gs"))        doNumber = &gridScale;
      else if (!strcmp(arg, "-gx"))        doReal   = &gridPerInchX;
      else if (!strcmp(arg, "-gy"))        doReal   = &gridPerInchY;
      else if (!strcmp(arg, "-pc"))        doNumber = &pointColor;
      else if (!strcmp(arg, "-ph"))        doNumber = &pointRange;
      else if (!strcmp(arg, "-ps"))        doNumber = &pointSaturation;
      else if (!strcmp(arg, "-pv"))        doNumber = &pointThreshold;
      else if (!strcmp(arg, "-od"))        doNumber = &outputDepth;
      else Usage(1); }
    else {
      if (argi+1 != argc || !*arg)  Usage(1);
      inpName = arg; } }
  if (doName || doNumber || doReal)  Usage(1);
  // Check options and set missing file names for calibration/unwarp mode:
  argi = 0;
  if (calibPicName)     ++argi;
  if (calibTextName)    ++argi;
  if (calibDeformName)  ++argi;
  if (argi != 1 || calibCheckName)  doCalib = 1;
  if (doCalib) {
    if (sourcePicName || destPicName) Usage(1);
    if (calibPicName) { if (inpName) Usage(1);  }
    else if (inpName) { calibPicName = inpName; }
    if (!calibPicName && (calibCheckName || !calibTextName)) calibPicName = "";
    if (!calibDeformName && (!calibPicName || !calibTextName))
      calibDeformName = ""; }
  else {
    if (!sourcePicName) sourcePicName = inpName ? inpName : "";
    else if (inpName)   Usage(1);
    if (!destPicName)   destPicName = ""; }
}

// ============================================================================
// General common functions:
// ----------------------------------------------------------------------------

// Error message with program termination:

static void Error(const char *format, ...)
{
  va_list ap;
  va_start(ap, format);
  vsprintf(glbBuf, format, ap);  // Hopefully this fits (long filenames)!
  va_end(ap);
  throw glbBuf;
}

// Normal messages, which can be suppressed:

static void Print(const char *format, ...)
{
  if (!param.quiet) {
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap); }
}

// ----------------------------------------------------------------------------

// Conversion of RGB-color components into HSV-components. This function
// is coded adhoc without respecting any norms. Its only purpose is
// based on calibration point detection:

static void convertRGBtoHSV(int r, int g, int b, int mcol,
                            double& h,double& s,double& v)
{
  double cfac = 1.0 / double(mcol);
  int    i;
  // Gray values have arbitrary hue:
  if (r == g && r == b) { h = s = 0.0;  v = cfac * double(r);  return; }
  // Rotate color components to have 'r,g > b' and adjust base of hue:
  if      (r < g) { if (r < b) { h = 1.0;  i = r;  r = g;  g = b;  b = i; }
                    else         h = 0.0; }
  else if (g < b) {              h = 2.0;  i = r;  r = b;  b = g;  g = i; }
  else                           h = 0.0;
  // r = v*(1-s)+v*s*(1-h),  g = v*(1-s)+v*s*h,  b = v*(1-s)
  r -= b;  g -= b;  r += g;  b += r;
  // r = v*s,  g = v*s*h,  b = v
  h = (h + double(g) / double(r)) / 3.0;
  if (h >= 1.0) h = 0.0;
  v = cfac * double(b);
  s = double(r) / double(b);
}

// ============================================================================

// Onedimensional cubic interpolation:
//
// 'f[4]' contains up to four interpolation coefficients, depending on
// the fraction between the aequidistant grid points. 'm' is the maximum
// possible grid index (between 0 and 'm'), 'o' the row/column offset
// to adjust the grid pointer to adjacent rows/columns. 'r = 1' indicates
// that the weighting process should start one step before the current
// row/column and 'n' denotes the valid coefficients 'f[0 ... n]'.
// The interpolation is correct if there are at least four grid
// points (that means 'm >= 3').

struct Cubic1
{
  double f[4];
  long   m, o, n, r;
  Cubic1(long m, long o) { this->m = m;  this->o = o; }
  long operator()(double t1);
};

// ----------------------------------------------------------------------------

// Initializes cubic interpolation at location '0.0 <= t1 <= m'.
// The interpolation is C1-continuous.

long Cubic1::operator()(double t1)
{
  double t0 = floor(t1);
  long   i  = long(t0);
  t1 -= t0;
  if (i >= m) { --i;  t1 += 1.0; }  // Border case.
  t0 = 1.0 - t1;
  if (i) {
    if (i == m - 1) {  // Highest grid interval:
      f[    0] = -0.5 * t0 * t1;
      f[r = 1] =       t0 * (2.0 - t0);
      f[n = 2] = 0.5 * t1 * (2.0 - t0); }
    else {  // General grid interval:
      f[    0] = -0.5 * t0 * t0 * t1;
      f[r = 1] = -t0 * (((1.5 * t1) - 1.0) * t1 - 1.0);
      f[    2] = -t1 * (((1.5 * t0) - 1.0) * t0 - 1.0);
      f[n = 3] = -0.5 * t1 * t1 * t0; } }
  else {  // Lowest grid interval:
    f[r = 0] = 0.5 * t0 * (2.0 - t1);
    f[    1] =       t1 * (2.0 - t1);
    f[n = 2] = -0.5 * t0 * t1; }
  return i;
}

// ----------------------------------------------------------------------------

// Twodimensional cubic interpolation:

template<class ityp>
struct Cubic2
{
  Cubic1 cx,cy;  // Onedimensional interpolation along x-axis and y-axis.
  Cubic2(long mx, long ox, long my, long oy) : cx(mx, ox), cy(my, oy) { }
  double operator()(ityp *raw);
};

// ----------------------------------------------------------------------------

// Returns interpolated value for previously specified grid point
// in the cell, whose lower left pointer is given by 'raw':
//
// The cubic interpolation on a rectangular grid has
// the tensor product property.

template<class ityp>
double Cubic2<ityp>::operator()(ityp *raw)
{
  double  value = 0.0;
  long    ix, iy;
  ityp   *r;
  if (cx.r) raw -= cx.o;
  if (cy.r) raw -= cy.o;
  for (ix=0; ix <= cx.n; ++ix,raw+=cx.o) {
    for (iy=0,r=raw; iy <= cy.n; ++iy,r+=cy.o) {
      value += double(*r) * (cx.f[ix] * cy.f[iy]); } }
  return value;
}

// ============================================================================

// Type for grid of calibration points:
//
// The detected points have indices '0 ... nAll-1'.
// The connected and used points have indices '0 ... n-1' with 'n <= nAll'.
// 'width' and 'height' is the size of the underlying image.
// 'x[]' and 'y[]' contain the coordinates of the detected points.
// 'ix[]' and 'iy[]' are integer values of their regular indices
// in the rectangular mesh, which connectivity is determined heuristically.

struct Grid
{
  long n, nAll, width, height, *x, *y, *ix, *iy;
  Grid()  { x = 0;  Reset(); }
  ~Grid() {         Reset(); }

  void Reset(long n=0, long width=0, long height=0);
  void Set(long i, long x, long y, long ix=0, long iy=0);
  void Read(const char *filename, FILE *fh=0);
  void Write(const char *filename, FILE *fh=0);
  void Limit(long& nx,long& ny);
};

// ----------------------------------------------------------------------------

// Initializes the internal arrays:

void Grid::Reset(long n, long width, long height)
{
  this->n      = nAll = n;
  this->width  = width;
  this->height = height;
  if (x) delete x;
  x = y = ix = iy = 0;
  if (n) {
    x  = new long[4*n];
    y  = x  + n;
    ix = y  + n;
    iy = ix + n; }
  for (long i=0; i < 4*n; ++i)  x[i] = 0;
}

// ----------------------------------------------------------------------------

// Sets data for one calibration point:

void Grid::Set(long i, long x, long y, long ix, long iy)
{
  if (i < 0 || i >= n)  throw "Fatal error in 'Grid::Set'";
  this->x[i]  = x;
  this->y[i]  = y;
  this->ix[i] = ix;
  this->iy[i] = iy;
}

// ----------------------------------------------------------------------------

// Reads the calibration point list:

void Grid::Read(const char *filename, FILE *fh)
{
  long k, kk, ixx, iyy, xx, yy;
  if (fh != stdin)  fh = fopen(filename, "r");
  if (!fh)  Error("Couldn't read calibration point file '%s'", filename);
  try {
    if (fscanf(fh,"Calibrationpoints %ld %ld %ld\n",&n,&width,&height) != 3 ||
        n <= 0 || width <= 0 || height <= 0 ||
        width > 25000 || height > 25000) throw 0;
    Reset(n, width, height);
    for (k=0; k < n; ++k) {
      if (fscanf(fh,"%ld\t%ld\t%ld\t%ld\t%ld\n",&kk,&ixx,&iyy,&xx,&yy) != 5 ||
          kk < 0 || ixx < 0 || iyy < 0 ||
          xx < 0 || xx >= width || yy < 0 || yy >= height)  throw 1;
      Set(k, xx, yy, ixx, iyy); }
    if (fh != stdin)  fclose(fh); }
  catch (int nmr) {
    const char *msg = nmr ? "has wrong point syntax" : "has wrong preamble";
    if (fh != stdin)  fclose(fh);
    Error("Calibration point file '%s' %s",filename,msg); }
}

// ----------------------------------------------------------------------------

// Writes the calibration point list:

void Grid::Write(const char *filename, FILE *fh)
{
  if (fh != stdout)  fh = fopen(filename, "w");
  if (!fh)  Error("Couldn't write calibration point file '%s'", filename);
  fprintf(fh, "Calibrationpoints %ld %ld %ld\n", n, width, height);
  // We write all detected points for possible human inspection:
  for (long k=0; k < nAll; ++k)
    fprintf(fh, "%ld\t%ld\t%ld\t%ld\t%ld\n", k, ix[k], iy[k], x[k], y[k]);
  if (fh != stdout)  fclose(fh);
}

// ----------------------------------------------------------------------------

// Calculates limits for grid indices:

void Grid::Limit(long& nx,long& ny)
{
  for (long i=nx=ny=0; i < n; ++i) {
    if (nx <= ix[i])  nx = ix[i] + 1;
    if (ny <= iy[i])  ny = iy[i] + 1; }
}

// ============================================================================

// Type for connectivity mesh beteen calibration points:

struct Mesh
{
  Grid& grid;
  long      nAll, n, *x, *y, *ix, *iy;
  long      *htab, *hnxt, *active;
  long      ixa, ixb, iya, iyb;
  u_int32_t nh, *hash;

  static u_int32_t HashValue(long xi, long yi);
  long             Find(long xi, long yi);
  void             Insert();
  long             Right(long k)  { return Find(ix[k] + 1, iy[k]);     }
  long             Left(long k)   { return Find(ix[k] - 1, iy[k]);     }
  long             Top(long k)    { return Find(ix[k],     iy[k] + 1); }
  long             Bottom(long k) { return Find(ix[k],     iy[k] - 1); }

  Mesh(Grid& g) : grid(g) { htab = 0; Reset(); }
  ~Mesh()                 {           Reset(); }

  void   Reset(int init=0);
  double CellError(long ind[9]);
  void   Seed();
  void   FindIndices();
  void   SortAndTranslate();
};

// ----------------------------------------------------------------------------

// Calculates 32bit hash value from grid indices:

u_int32_t Mesh::HashValue(long xi, long yi)
{
  u_int32_t xx = u_int32_t(xi), yy = u_int32_t(yi);
  for (int i=0; i < 5; ++i) {
    xx = xx ^ (xx << 5)  ^ (xx >> 13);
    yy = yy ^ (yy << 17) ^ (yy >> 3); }
  return xx ^ yy;
}

// ----------------------------------------------------------------------------

// Searches for a grid point in the hash table:

long Mesh::Find(long xi, long yi)
{
  u_int32_t h = HashValue(xi, yi);
  long      k = htab[h % nh];
  while (k >= 0) {
    if (ix[k] == xi && iy[k] == yi) return k;
    k = hnxt[k]; }
  return -1;
}

// ----------------------------------------------------------------------------

// Insert a new grid point with index 'n' into the hash table:

void Mesh::Insert()
{
  long      xx = ix[n], yy = iy[n];
  u_int32_t hi;
  // Track extreme grid indices:
  if (ixa > xx)  ixa = xx;  if (ixb < xx)  ixb = xx;
  if (iya > yy)  iya = yy;  if (iyb < yy)  iyb = yy;
  hi       = (hash[n] = HashValue(xx, yy)) % nh;
  hnxt[n]  = htab[hi];
  htab[hi] = n++;
  if (u_int32_t(n) >= nh && nh <= 0x10000000) {  // Enlarge hash table:
    delete htab;
    htab = 0;  // If new throws an error.
    htab = new long[nh <<= 1];
    for (hi=0; hi < nh; ++hi)  htab[hi] = -1;
    for (long k=0; k < n; ++k) {
      hi       = hash[k] % nh;
      hnxt[k]  = htab[hi];
      htab[hi] = k; } }
}

// ----------------------------------------------------------------------------

// Initializes internal arrays:

void Mesh::Reset(int init)
{
  if (htab) { delete htab;  delete hnxt;  delete hash; }
  nAll = n = ixa = ixb = iya = iyb = 0;
  nh = 0;
  x = y = ix = iy = htab = hnxt = active = 0;
  hash = 0;
  if (init) {
    n      = nAll = grid.nAll;
    x      = grid.x;
    y      = grid.y;
    ix     = grid.ix;
    iy     = grid.iy;
    htab   = new long[(nh = 16)];
    hnxt   = new long[2*n];
    active = hnxt + n;
    hash   = new u_int32_t[n];
    for (u_int32_t i=0; i < nh; ++i)  htab[i] = -1; }
}

// ----------------------------------------------------------------------------

// Calculates the cell error for the nine choosen points. Used in 'Seed()':
//
// This method optimizes 'rx,ry,ax,ay,bx,by' under the orthogonality
// constraint 'ax*bx+ay*by == 0' with respect to the squared distance
// error of the nine points with indices 'ind[0..8]'
//
//   P6  P3  P5
//   P2  P0  P1  ----> (ax,by)
//   P7  P4  P8
//
// and the rectangular cell with points
// '(rx+i*ax+k*bx,ry+i*ay+k*by)', 'i,k=-1,0,1'.
// It turns out, that '(rx,ry)' is just the mean value of the nine points
// '(x(0..8),y(0..8))'. The squared error '12 * psi + constant'
// is subsequently given by
//
// psi = (ax^2 + ay^2 + bx^2 + by^2)/2 - (ax*ux+ay*uy+bx*vx+by*vy)
//
// ux = (x1 + x5 + x8 - x2 - x6 - x7) / 6
// uy = (y1 + y5 + y8 - y2 - y6 - y7) / 6
// vx = (x3 + x5 + x6 - x4 - x7 - x8) / 6
// vy = (y3 + y5 + y6 - y4 - y7 - y8) / 6
//
// Using '(ax,ay) = alpha*(c,s)' and '(bx,by) = beta*(-s,c)' with
// 'c^2 + s^2 == 1' for the orthogonality restriction, this leads to
//
// psi = (alpha^2 + beta^2) / 2 - (c*ux+s*uy)*alpha - (c*vy-s*vx)*beta.
//
// Therefore neccessarily for a minimum value of 'psi' we have
//
// alpha = c * ux + s * uy,
// beta  = c * vy - s * vx,
//
// and substituting this into 'psi' gives
//
// 2*psi = 2*(vx*vy-ux*uy)*c*s - (ux^2+vy^2)*c^2 - (uy^2+vx^2)*s^2.
//
// We change the trigonometric terms 'c=cos(phi), s=sin(phi)' to
// 'c2=cos(2*phi), s2=sin(2*phi)' by
//
// c^2 = (1 + c2)/2, s^2 = (1 - c2)/2, 2*c*s = s2.
//
// This leads to
//
// -4*psi = (ux^2-uy^2+vy^2-vx^2)*c2 + 2*(ux*uy-vx*vy)*s2 + constant.
//
// This can be solved for 'phi' and the optimum is found.
// Note that 'c2,s2' cannot be found for the two cases
// 'u = v' or 'u = -v', but these cases are of degenerated nature
// and are discarded.

double Mesh::CellError(long ind[9])
{
  // Mapping of the indices 0..8 to the positions in the regular cell:
  double cellA[] = {  0.0,  1.0, -1.0,  0.0,  0.0,  1.0, -1.0, -1.0,  1.0 },
         cellB[] = {  0.0,  0.0,  0.0,  1.0, -1.0,  1.0,  1.0, -1.0, -1.0 };
  double rx, ry, ux, uy, vx, vy, c, s, phi, alpha, beta,
         ax, ay, bx, by, dx, dy, err;
  int    i;
  // Mean value of grid points is optimal centrum:
  for (i=0,rx=ry=0.0; i < 9; ++i) {
    rx += double(x[ind[i]]);
    ry += double(y[ind[i]]); }
  rx /= 9.0;  ry /= 9.0;
  // Calculate optimal orthogonal rotation 'c=cos(2*phi), s=sin(2*phi)':
  ux = double(x[ind[1]]+x[ind[5]]+x[ind[8]]
              -x[ind[2]]-x[ind[6]]-x[ind[7]]) / 6.0;
  uy = double(y[ind[1]]+y[ind[5]]+y[ind[8]]
              -y[ind[2]]-y[ind[6]]-y[ind[7]]) / 6.0;
  vx = double(x[ind[3]]+x[ind[5]]+x[ind[6]]
              -x[ind[4]]-x[ind[7]]-x[ind[8]]) / 6.0;
  vy = double(y[ind[3]]+y[ind[5]]+y[ind[6]]
              -y[ind[4]]-y[ind[7]]-y[ind[8]]) / 6.0;
  c  = vy * vy - vx * vx + ux * ux - uy * uy;
  s  = 2.0 * (ux * uy - vx * vy);
  if (fabs(c) + fabs(s) <= 1e-8)  return -1.0;  // Obviously degenerated case.
  phi = 0.5 * atan2(s,c);
  // Calculate regular grid axis '(ax,ay)' and '(bx,by)':
  c     = cos(phi);
  s     = sin(phi);
  alpha = c  * ux + s * uy;
  beta  = c  * vy - s * vx;
  ax    = c  * alpha;
  ay    = s  * alpha;
  bx    = -s * beta;
  by    = c  * beta;
  // Calculate summed squared distances between exact and real points:
  for (i=0,err=s=0.0; i < 9; ++i) {
    dx   = rx - double(x[ind[i]]);
    dy   = ry - double(y[ind[i]]);
    s   += dx * dx + dy * dy;
    dx  += cellA[i] * ax + cellB[i] * bx;
    dy  += cellA[i] * ay + cellB[i] * by;
    err += dx * dx + dy * dy; }
  // To be scale invariant, the error is divided by the point variances:
  return err / s;
}

// ----------------------------------------------------------------------------

// Creates 3x3 size start cell:
//
// The connectivity construction is quite heuristic. Skewed meshes are
// distorted, but hopefully not too much to find an almost rectangular
// 3x3 block of calibration points. This seed cell is searched for here.

void Mesh::Seed()
{
  long v = 0, cx, cy, dx, dy, i, j, k, m;
  long iBest[9], ind[10], value[10], xv[9], yv[9];
  double best = -1.0, bb;
  n = nAll;  // All detected points should be considered.
  if (n < 9) throw "At least nine calibration points are needed";
  // Check for each point, if a suitable cell can be found
  // with that point at the center (index 0):
  for (i=0; i < n; ++i) {
    // The point coordinates are saved in (cx,cy):
    cx = x[i];  cy = y[i];  m = 0;
    // Search the nine nearest points to (cx,cy), currently 'm' are found:
    for (j=0; j < n; ++j) {
      dx = x[j] - cx;  dy = y[j] - cy;  v = dx * dx + dy * dy;
      for (k=m; k > 0; --k) {
        if (v >= value[k-1])  break;
        ind[k] = ind[k-1];  value[k] = value[k-1]; }
      ind[k] = j;  value[k] = v;
      // Up to nine nearest points are saved:
      if (m < 9) ++m; }
    // Now nine nearest points with indices 'ind[0..8]' are found, including
    // the original center point. Arrange them in the best possible way
    // to make a rectangular cell:
    for (j=0; j < 9; ++j) { xv[j] = x[ind[j]] - cx;  yv[j] = y[ind[j]] - cy; }
    // The points have relative coordinates 'xv[0..8]' and 'yv[0..8]'.
    // Point 0 and 1 are already found according to the sorting by distance.
    for (j=2; j < 8; ++j) {
      // Set the best place (cx,cy) for the next cell point:
      switch (j) {
        case 2:   cx = -xv[1];           cy = -yv[1];           break;
        case 3:   cx = (yv[2]-yv[1])/2;  cy = (xv[1]-xv[2])/2;  break;
        case 4:   cx = -xv[3];           cy = -yv[3];           break;
        case 5:   cx = xv[1] + xv[3];    cy = yv[1] + yv[3];    break;
        case 6:   cx = xv[3] + xv[2];    cy = yv[3] + yv[2];    break;
        default:  cx = xv[2] + xv[4];    cy = yv[2] + yv[4]; }
      // Search in the available points for best fit to (cx,cy):
      for (k=j,m=0; k < 9; ++k) {
        dx = xv[k] - cx;  dy = yv[k] - cy;  dx = dx * dx + dy * dy;
        if (!m || v > dx) { m = k;  v = dx; } }
      // Swap cell points 'j' and 'm' to get best fit 'm' at position 'j':
      v = ind[j];  ind[j] = ind[m];  ind[m] = v;
      v = xv[j];   xv[j]  = xv[m];   xv[m]  = v;
      v = yv[j];   yv[j]  = yv[m];   yv[m]  = v; }
    // Now the best match of the nearest points must be rated:
    bb = CellError(ind);
    // If a better fit is found, save the objective value in 'best' and
    // the corresponding point indices in 'iBest[0..9]':
    if (bb >= 0.0 && (best < 0.0 || best > bb)) {
      for (j=0,best=bb; j < 9; ++j) iBest[j] = ind[j]; } }
  // This shouldn't happen too often:
  if (best <= 0.0) throw "Calibration points are too distorted";
  // Now 'iBest' contains the indices of the detected seed cell.
  // Rotate cell by 90 degrees if neccessary:
  dx = x[iBest[1]] - x[iBest[2]];  if (dx < 0) dx = -dx;
  dy = y[iBest[1]] - y[iBest[2]];  if (dy < 0) dy = -dy;
  if (dx < dy) {
    m = iBest[1];  iBest[1] = iBest[3];  iBest[3] = iBest[2];
                   iBest[2] = iBest[4];  iBest[4] = m;
    m = iBest[5];  iBest[5] = iBest[6];  iBest[6] = iBest[7];
                   iBest[7] = iBest[8];  iBest[8] = m; }
  // Rotate cell by 180 degrees if neccessary:
  dx = x[iBest[1]] - x[iBest[2]];
  if (dx < 0) {
    m = iBest[1];  iBest[1] = iBest[2];  iBest[2] = m;
    m = iBest[3];  iBest[3] = iBest[4];  iBest[4] = m;
    m = iBest[5];  iBest[5] = iBest[7];  iBest[7] = m;
    m = iBest[6];  iBest[6] = iBest[8];  iBest[8] = m; }
  // Moves points of cell to the beginning of the point arrays:
  for (j=0; j < 9; ++j) {
    v = x[j];  x[j] = x[iBest[j]];  x[iBest[j]] = v;
    v = y[j];  y[j] = y[iBest[j]];  y[iBest[j]] = v;
    for (k=j+1; k < 9; ++k)
      if (iBest[k] == j) { iBest[k] = iBest[j]; break; } }
}

// ----------------------------------------------------------------------------

// Find recursively neighbouring points from seed cell:

void Mesh::FindIndices()
{
  long ix_init[]  = { 0,  1, -1,  0,  0,  1, -1, -1,  1 },
       iy_init[]  = { 0,  0,  0,  1, -1,  1,  1, -1, -1 };
  long actBeg     = 0,  // Marks the range of extension
       actEnd     = 8,  // point indices in 'active'.
       actSync    = 8, r = 0, i, j, k, kk, m, dx, dy, rc;
  int  changed    = 0;  // To prevent infinite looping.
  // grid indices of seed cell are set to '-1,0,1':
  for (i=0; i < 9; ++i) {
    ix[i] = ix_init[i];
    iy[i] = iy_init[i];
    // The outer 8 points are extension points:
    active[i] = i + 1; }
  // Hash table init:
  for (u_int32_t hi=0; hi < nh; ++hi)  htab[hi] = -1;
  for (n=0; n < 9;)  Insert();
  // Initialize markers for lowest and largest grid indices found so far:
  ixa = iya = -1;
  ixb = iyb =  1;
  // Loop until all points are inserted or no extension points are left:
  while (actBeg != actEnd && n < nAll) {
    int redo = 0;  // Redo if there is currently no data for extrapolation.
    int ka   = active[actBeg];
    // If the mesh is not enlarged since the last synchronization
    // marker, there is no point to continue the loop forever:
    if (actBeg == actSync) {
      if (!changed)  break;
      actSync = actEnd;
      changed = 0; }
    if (++actBeg >= nAll)  actBeg = 0;
    // Try to continue point in all four directions:
    for (i=0; i < 4 && n < nAll; ++i) {
      // The work array contains known         ?     work[2]     ?
      // mesh points relative to one        work[1]   <new>   work[0]
      // <new>==-1 unknown positions:          ?     work[3]     ?
      long work[] = {-1,-1,-1,-1};
      switch (i) {
        case 0:   if (Right(ka)  < 0)  work[1] = ka;  break;
        case 1:   if (Left(ka)   < 0)  work[0] = ka;  break;
        case 2:   if (Top(ka)    < 0)  work[3] = ka;  break;
        default:  if (Bottom(ka) < 0)  work[2] = ka;  }
      for (j=0; j < 4; ++j)  if (work[j] >= 0)  break;
      // If any of the four neighbours is not in the mesh, this neighbour
      // position will be connected to a new point if possible.
      // If all neighbours are already in the mesh, the current point
      // will be discarded from the active list:
      if (j < 4) {
        long xc = 0, yc = 0, nc = 0;
        // Try to set other neighbour positions to get a better estimation
        // of a valid (x,y) position for the point to be joined to the mesh:
        for (j=0; j < 3; ++j) {
          if ((k = work[0]) >= 0) {
            if ((m = Top(k))    >= 0 && (m = Left(m)) >= 0)  work[2] = m;
            if ((m = Bottom(k)) >= 0 && (m = Left(m)) >= 0)  work[3] = m; }
          if ((k = work[1]) >= 0) {
            if ((m = Top(k))    >= 0 && (m = Right(m)) >= 0)  work[2] = m;
            if ((m = Bottom(k)) >= 0 && (m = Right(m)) >= 0)  work[3] = m; }
          if ((k = work[2]) >= 0) {
            if ((m = Left(k))  >= 0 && (m = Bottom(m)) >= 0)  work[1] = m;
            if ((m = Right(k)) >= 0 && (m = Bottom(m)) >= 0)  work[0] = m; }
          if ((k = work[3]) >= 0) {
            if ((m = Left(k))  >= 0 && (m = Top(m)) >= 0)  work[1] = m;
            if ((m = Right(k)) >= 0 && (m = Top(m)) >= 0)  work[0] = m; } }
        // Find all 'nc' possible mesh continuations by
        // extrapolation from inner mesh points:
        for (j=0; j < 4; ++j) {
          if ((k = work[j]) >= 0) {
            // Horizontal or vertical extrapolation:
            switch (j) {
              case 0:   m = Right(k);   break;
              case 1:   m = Left(k);    break;
              case 2:   m = Top(k);     break;
              default:  m = Bottom(k); }
            if (m >= 0) { ++nc;  xc += 2*x[k] - x[m];  yc += 2*y[k] - y[m]; }
            // Diagonal extrapolation:
            switch (j) {
              case 0:   kk = work[2];  m = Top(k);     break;
              case 1:   kk = work[3];  m = Bottom(k);  break;
              case 2:   kk = work[1];  m = Left(k);    break;
              default:  kk = work[0];  m = Right(k);   }
            if (kk >= 0 && m >= 0) {
              ++nc;
              xc += x[k] + x[kk] - x[m];
              yc += y[k] + y[kk] - y[m]; } } }
        // If at least one extrapolation is found, a new point for the
        // open position will be searched for, otherwise there might be
        // a possible extrapolation later when the mesh is grown at other
        // boundary points. Therefore in the latter case the active point
        // is reinserted into the active list:
        if (nc) {
          // Calculate the mean of all extrapolations found:
          xc /= nc;
          yc /= nc;
          // Calculate a reference squared distance 'rc'
          // to the known surrounding mesh points:
          for (j=0,nc=rc=0; j < 4; ++j) {
            if ((k = work[j]) >= 0) {
              ++nc;
              dx  = x[k] - xc;
              dy  = y[k] - yc;
              rc += dx * dx + dy * dy; } }
          rc /= nc;
          // From all available points the nearest point to the dedicated
          // position '(xc,yc)' has index 'k' and squared distance 'r':
          for (j=n,k=-1; j < nAll; ++j) {
            dx = x[j] - xc;
            dy = y[j] - yc;
            dx = dx * dx + dy * dy;
            if (k < 0 || r > dx) { k = j;  r = dx; } }
          // If the distance from the choosen point to the optimal position
          // '(xc,yc)' is about 1/3 to the reference distance, the point will
          // be inserted into the mesh, otherwise this mesh will not be
          // extended from the current boundary point in the current direction:
          if (r <= rc / 10) {
            changed = 1;  // The mesh will be extended now.
            // Put selected point 'k' in insertion position 'n':
            r = x[n];  x[n] = x[k];  x[k] = r;
            r = y[n];  y[n] = y[k],  y[k] = r;
            // Find correct grid '(ix[n],iy[n])' indices for the new point:
            for (j=0; j < 4; ++j) {
              if ((k = work[j]) >= 0) {
                switch (j) {
                  case 0:   ix[n]  = ix[k] - 1;  iy[n]  = iy[k];  break;
                  case 1:   ix[n]  = ix[k] + 1;  iy[n]  = iy[k];  break;
                  case 2:   iy[n]  = iy[k] - 1;  ix[n]  = ix[k];  break;
                  default:  iy[n]  = iy[k] + 1;  ix[n]  = ix[k]; }
                break; } }
            // The new point becomes active as possible new extension point:
            active[actEnd++] = n;
            if (actEnd >= nAll)  actEnd = 0;
            if (actBeg == actEnd)
              throw "Fatal error in Mesh::FindIndices";
            // Insert point 'n' into the hash table:
            Insert(); } }
        // Signal reinsertion of current point in active list.
        else redo = 1; } }
    // Reinsert current point in active list on demand:
    if (redo) {
      active[actEnd++] = ka;
      if (actEnd >= nAll)  actEnd = 0;
      if (actBeg == actEnd)
        throw "Fatal error in Mesh::FindIndices"; } }
  // Now 'n' from the 'nAll' detected points are connected by the mesh.
  // These points are located at the index positions '0 .. n-1'.
  if (n < nAll)  Print("Only %ld detected points used for calibration!\n",n);
  grid.n = n;  // Set correct number of points in coupled grid structure.
}

// ----------------------------------------------------------------------------

// Sort the grid points according to the underlying grid and translate
// the grid indices to nonnegative values:
//
// Even the translation of the indices will destroy the hash table,
// therefore this method has to be executed at the end of mesh handling.

void Mesh::SortAndTranslate()
{
  long nn = 0, xi,yi, k,ak,hk;
  // Walk through the grid in sort order using the hash table:
  for (yi=iya; yi <= iyb; ++yi) {
    for (xi=ixa; xi <= ixb; ++xi) {
      if ((k = Find(xi,yi)) >= 0)  active[nn++] = k; } }
  if (n != nn) throw "Fatal error in Mesh::SortAndTranslate";
  // Now the hash table is not needed anymore.
  // Build the inverse table for resorting:
  for (k=0; k < n; ++k) hnxt[active[k]] = k;
  // Now sort the points:
  for (k=0; k < n; ++k) {
    // Adjust sort information:
    ak         = active[k];
    hk         = hnxt[k];
    active[hk] = ak;
    hnxt[ak]   = hk;
    // Swap the entries:
    hk = x[k];   x[k]  = x[ak];   x[ak]  = hk;
    hk = y[k];   y[k]  = y[ak];   y[ak]  = hk;
    hk = ix[k];  ix[k] = ix[ak];  ix[ak] = hk;
    hk = iy[k];  iy[k] = iy[ak];  iy[ak] = hk; }
  // Translate grid indices to nonnegative values:
  for (k=0; k < n; ++k) { ix[k] -= ixa;  iy[k] -= iya; }
}

// ============================================================================

// Type for image data:

struct Picture {
  long     width,   // Image width
           height,  // Image height
           depth,   // Maximum color value per channel.
           size;    // 3 * width * height, each pixel has RGB-tripel.
  pixel_t *pixel;   // 3 * width * height RGB values (16 bit).

  Picture()  { pixel = 0;  Reset(); }
  ~Picture() {             Reset(); }

  void Reset(long width=0, long height=0, long depth=255);
  void Read(const char *filename, FILE *fh=0);
  void Write(const char *filename, FILE *fh=0);
  void FilterRedPoints();
  int  SwapRedGreenPoint(long x, long y, int canal=1);
  long CountRedPoints();
  void CollectPoints(Grid &grid);
  void Dot(long x, long y, int white,int swap=0);
  void Bresenham(long x0, long y0, long x1, long y1);
  void MakeWhite(Grid &grid);
  void CreateMesh(Grid& grid);
  void Calibrate(Grid& grid);
};

// ----------------------------------------------------------------------------

// Resizes the image array by given size parameters:

void Picture::Reset(long width, long height, long depth)
{
  if (pixel)  delete pixel;
  this->width  = width;
  this->height = height;
  this->depth  = depth;
  size         = 3 * width * height;
  pixel        = 0;
  if (size)  pixel = new pixel_t[size];
}

// ----------------------------------------------------------------------------

// Read raw PPM image file:

void Picture::Read(const char *filename, FILE *fh)
{
  int      mode = 0, ch = 0;
  long     z = 0, n, m;
  u_int8_t buf[4096], *src;
  pixel_t  val,       *dst;
  if (fh != stdin)  fh = fopen(filename, "rb");
  else              binmode(fh);  // Windows.
  if (!fh)  Error("Couldn't read '%s' as raw PPM image file", filename);
  try {
    // Read preamble of PPM file:
    while (mode != 5) {
      // Do we need a next character?
      if      (mode & 64)                mode &= 63;
      else if ((ch = fgetc(fh)) == EOF)  throw 1;
      // Do we need to read space with comments?
      if (mode & 8) {
        if (isspace(ch)) {
          mode |= 16;  // We have seen at least one space/comment.
          if (ch == '\n')  mode &= 31;
          continue; }
        else if (ch == '#') { mode |= 48;  continue; }
        else {
          if (mode & 32)  continue;   // Ignore comment characters in line.
          if (!(mode & 16)) throw 1;  // Was there a separation?
          // Switch to number reading:
          mode &= 23; } }
      // Do we need to read a number?
      if (mode & 16) {
        if (!(mode & 32))  z = 0;  // Initialize at beginning.
        if (ch >= '0' && ch <= '9') {
          mode |= 32;  // We have seen at least one digit.
          z     = 10 * z + (ch - '0');
          if (z > 65535)  throw 1;
          continue; }
        if (!(mode & 32))  throw 1;  // We need at least one digit.
        // Switch to global parsing mode:
        mode &= 7; }
      // Main parsing of preamble:
      switch (mode) {
        case 0:   if (ch == 'P') { mode =  1;  continue; }  throw 1;
        case 1:   if (ch == '6') { mode = 10;  continue; }  throw 1;
        case 2:   width  = z;  mode = 75;  continue;
        case 3:   height = z;  mode = 76;  continue;
        default:  depth  = z;  mode =  5; } }
    // Now the preamble has been consumed, check for final space character:
    if (width > 25000 || height > 25000 || !depth || !isspace(ch))  throw 1;
    Reset(width, height, depth);
    // Read binary data; either 8bit or 16bit RGB channel values:
    n   = size;
    dst = pixel;
    if (depth > 255) {
      while (n) {
        n -= m = n < 2048 ? n : 2048;
        if (fread((src = buf), sizeof(u_int8_t) << 1, m, fh) != ulong(m))
          throw 2;
        while (m--) {
          val    = pixel_t(*src++);
          *dst++ = val = (val << 8) | pixel_t(*src++);
          if (val > depth)  throw 4; } } }
    else {
      while (n) {
        n -= m = n < 4096 ? n : 4096;
        if (fread((src = buf), sizeof(u_int8_t), m, fh) != ulong(m))
          throw 2;
        while (m--) {
          *dst++ = val = pixel_t(*src++);
          if (val > depth)  throw 4; } } }
    if (fgetc(fh) != EOF)  throw 3;
    if (fh != stdin)  fclose(fh); }
  catch (int nmr) {
    const char *msg;
    if (fh != stdin)  fclose(fh);
    switch (nmr) {
      case 1:   msg = "has wrong preamble";     break;
      case 2:   msg = "is too small";           break;
      case 3:   msg = "is too big";             break;
      default:  msg = "has wrong color depth";  }
    Error("PPM image file '%s' %s", filename, msg); }
}

// ----------------------------------------------------------------------------

// Write raw PPM image file:

void Picture::Write(const char *filename, FILE *fh)
{
  long     n, m, k;
  u_int8_t buf[4096], *dst;
  pixel_t  val,       *src;
  if (fh != stdout)  fh = fopen(filename, "wb");
  else               binmode(fh);  // Windows.
  if (!fh) Error("Couldn't write '%s' as raw PPM image file", filename);
  fprintf(fh, "P6\n# CREATOR: %s\n%ld %ld\n%ld\n",
          param.prgName, width, height, depth);
  n   = size;
  src = pixel;
  if (depth > 255) {
    while (n) {
      n   -= m = k = n < 2048 ? n : 2048;
      dst  = buf;
      while (k--) {
        val    = *src++;
        *dst++ = u_int8_t(val >> 8);
        *dst++ = u_int8_t(val & 255); }
      fwrite(buf, sizeof(u_int8_t) << 1, m, fh); } }
  else {
    while (n) {
      n   -= m = k = n < 4096 ? n : 4096;
      dst  = buf;
      while (k--)  *dst++ = u_int8_t(*src++);
      fwrite(buf, sizeof(u_int8_t), m, fh); } }
  if (fh != stdout)  fclose(fh);
}

// ----------------------------------------------------------------------------

// Point detection routine. Detected points are in red color, all
// other pixels are reset to black color:
//
// mode == 0:  Calculate maximum saturation value in 'm' of all
//             points within hue range.
// mode == 1:  Make all points black, which are not in hue range
//             or less saturated than 'pSaturation * m'.
// mode == 2:  Find maximum value 'm' of 's * v'.
//             Filtering with 'saturation * volume' seems to work
//             better than solely filtering with 'volume'.
// mode == 3:  Make all points black with 's * v < pThreshold * m'.
//             All other points will be red with red canal value
//             proportional to 's * v'.

void Picture::FilterRedPoints()
{
  int      pRed        = (param.pointColor >> 16) & 0xff,
           pGreen      = (param.pointColor >>  8) & 0xff,
           pBlue       =  param.pointColor        & 0xff;
  double   pRange      = double(param.pointRange)      / 100.0,
           pSaturation = double(param.pointSaturation) / 100.0,
           pThreshold  = double(param.pointThreshold)  / 100.0;
  double   pHue = 0.0, m = 0.0, mz = 0.0, h, s, v, z;
  int      mode;
  pixel_t *raw;
  long     nRaw;
  if (pRed == pGreen && pGreen == pBlue) {
    mode = 2;
    Print("Hue specification '-ph' ignored for gray value points!\n"); }
  else  mode = 0;
  for (; mode < 4; ++mode) {
    switch (mode) {
      case 0:  convertRGBtoHSV(pRed, pGreen, pBlue, 255, pHue, s, v);
               m = 1e-8;
               break;
      case 1:  m = pSaturation * m;
               break;
      case 2:  mz = 0.0;
               break;
      default: if (mz > 0.0)  m = pThreshold * mz;
               else           throw "No calibration points detected"; }
    for (nRaw=size,raw=pixel; nRaw; nRaw-=3,raw+=3) {
      convertRGBtoHSV(raw[0], raw[1], raw[2], depth, h, s, v);
      if (mode < 2) {
        z = 2.0 * fabs(fmod(h + 1.0 - pHue, 1.0) - 0.5);
        if (mode) { if (z < pRange || s < m) raw[0] = raw[1] = raw[2] = 0; }
        else      { if (z >= pRange && m < s)  m = s; } }
      else {
        z = s * v;  // Filtering seems to work better including saturation.
        if (mode >= 3) {
          raw[0] = z >= m ? pixel_t((255.0 * z) / mz) : 0;
          raw[1] = raw[2] = 0; }
        else if (mz < z) mz = z; } } }
  // Control picture will have standard color depth:
  depth = 255;
}

// ----------------------------------------------------------------------------

// Determines center of points and changes color between red and green:
//
// The detected area must have a complete rectangular boundary of two points
// with zero color canal ('canal = 0,1'). The detected color in canal
// 'canal' is swaped with the component of canal '1-canal'.
// Returns 1, if point was found, otherwise 0.

int Picture::SwapRedGreenPoint(long x,long y,int canal)
{
  long     width3 = 3*width, limit[4];
  pixel_t *base   = pixel + width3*y + 3*x, *raw;
  long     xx, yy, dx, dy, i, k, n, iRaw, kRaw, kLow, kHig;
  double   weight, xvalue, yvalue;
  int      border[4], value, dir, empty;
  if ((value = base[canal])) {
    base[canal]   = 0;
    base[1-canal] = value;
    weight        = double(value);
    xvalue        = yvalue = 0.0;
    for (i=0; i < 4; ++i)  limit[i] = border[i] = 0;
    for (;;) {
      for (dir=0; dir < 4; ++dir)
        if (0 <= border[dir] && border[dir] < 2) break;
      if (dir >= 4) break;
      // Go one vertical/horizontal pixel line further in direction 'dir':
      i = limit[dir] + 1;
      if (dir < 2) {
        if (dir) { i = -i;  if (x + i < 0)      dir += 4; }
        else                if (x + i >= width) dir += 4;
        iRaw  = 3;  kRaw  = width3;     kLow  = 2;  kHig  = 3;
        xx    = i;  yy    = -limit[2];  dx    = 0;  dy    = 1; }
      else {
        if (dir == 2) { i = -i;  if (y + i < 0)       dir += 4; }
        else                     if (y + i >= height) dir += 4;
        iRaw  = width3;  kRaw  = 3;          kLow  = 1;  kHig = 0;
        yy    = i;       xx    = -limit[1];  dx    = 1;  dy   = 0; }
      if (dir < 4) {  // No picture boundary yet:
        n      = limit[kLow];
        raw    = base + (iRaw * i - kRaw * n);
        n     += limit[kHig] + 1;
        empty  = 1;
        for (k=0; k < n; ++k,raw+=kRaw,xx+=dx,yy+=dy) {
          if ((value = raw[canal])) {
            raw[canal]   = empty = 0;
            raw[1-canal] = value;
            weight  += double(value);
            xvalue  += double(value * xx);
            yvalue  += double(value * yy);
            if      (k == 0   && border[kLow] > 0)  border[kLow] = 0;
            else if (k == n-1 && border[kHig] > 0)  border[kHig] = 0;
            else if (k == 1   && border[kLow] > 1)  border[kLow] = 1;
            else if (k == n-2 && border[kHig] > 1)  border[kHig] = 1; } }
        if (empty) border[dir]++;
        else       border[dir] = 0;
        limit[dir]++; }
      // Picture boundary reached:
      else border[dir-4] = -1; }
    // Only mark center point if real detection red => green is done:
    if (!canal) {
      xvalue = round(xvalue / weight);
      yvalue = round(yvalue / weight);
      raw    = base + (width3 * long(yvalue) + 3 * long(xvalue));
      raw[2] = 255; }
    return 1; }
  return 0;
}

// ----------------------------------------------------------------------------

// Count all detected calibration points:

long Picture::CountRedPoints()
{
  long count = 0, x, y;
  for (y=0; y < height; ++y)
    for (x=0; x < width; ++x)  count += SwapRedGreenPoint(x, y, 0);
  return count;
}

// ----------------------------------------------------------------------------

// Collect all detected calibration points in the point grid.
// The detected points are characterized by a blue color component of 255:

void Picture::CollectPoints(Grid &grid)
{
  pixel_t *raw;
  long     i = 0, x, y;
  for (y=0,raw=pixel; y < height; ++y)
    for (x=0; x < width; ++x,raw+=3)
      if (raw[2] == 255)  grid.Set(i++, x, y);
  if (i != grid.n)  throw "Fatal error in Picture::CollectPoints";
}

// ----------------------------------------------------------------------------

// Draws a white/blue star at (x,y) or (y,x) if 'swap != 0':

void Picture::Dot(long x, long y, int white, int swap)
{
  long ix, iy, iy0, iy1;
  if (swap) { ix = x;  x = y;  y = ix; }  // Handy for 'Bresenham'.
  for (ix=-2; ix <= 2; ++ix) {
    iy0 = (ix<0?-ix:ix);
    iy1 = 2-iy0;
    iy0 = iy0-2;
    for (iy=iy0; iy <= iy1; ++iy) {
      long xx = x + ix, yy = y + iy;
      if (0 <= xx && xx < width && 0 <= yy && yy < height) {
        pixel_t *raw = pixel + 3*(width * yy + xx);
        if (white)  raw[0] = raw[1] = raw[2] = 255;
        else { if (!raw[0] && !raw[1])  raw[2] = 255; } } } }
}

// ----------------------------------------------------------------------------

// Draws a blue line to display the mesh connectivity:

void Picture::Bresenham(long x0, long y0, long x1, long y1)
{
  long dx = 1, dy = 1, rx = x1 - x0, ry = y1 - y0, z;
  int  swap = 0;
  if (rx < 0) { dx = -1;  rx = -rx; }
  if (ry < 0) { dy = -1;  ry = -ry; }
  if (rx < ry) {
    z = x0;  x0 = y0;  y0 = z;
    z = x1;  x1 = y1;  y1 = z;
    z = dx;  dx = dy;  dy = z;
    z = rx;  rx = ry;  ry = z;
    swap = 1; }
  z = rx >> 1;
  Dot(x0, y0, 0, swap);
  while (x0 != x1) {
    x0 += dx;
    Dot(x0, y0, 0, swap);
    z += ry;
    if (z >= rx) {
      y0 += dy;
      Dot(x0, y0, 0, swap);
      z -= rx; } }
}

// ----------------------------------------------------------------------------

// Display all calibration points with white stars:

void Picture::MakeWhite(Grid &grid)
{
  for (long i=0; i < grid.nAll; ++i)  Dot(grid.x[i], grid.y[i], 1);
}

// ----------------------------------------------------------------------------

// Creates and displays the mesh of calibration points:

void Picture::CreateMesh(Grid& grid)
{
  long i, k;
  Mesh mesh(grid);
  mesh.Reset(1);
  mesh.Seed();
  mesh.FindIndices();
  for (i=grid.n; i < grid.nAll; ++i)
    SwapRedGreenPoint(grid.x[i], grid.y[i], 1);
  for (i=0; i < mesh.n; ++i) {
    if ((k = mesh.Right(i)) >= 0)
      Bresenham(mesh.x[i], mesh.y[i], mesh.x[k], mesh.y[k]);
    if ((k = mesh.Top(i)) >= 0)
      Bresenham(mesh.x[i], mesh.y[i], mesh.x[k], mesh.y[k]); }
  mesh.SortAndTranslate();
}

// ----------------------------------------------------------------------------

// Determines the calibration grid points and mesh:

void Picture::Calibrate(Grid& grid)
{
  long nPoints;
  FilterRedPoints();
  nPoints = CountRedPoints();
  Print("Number of detected points: %ld\n",nPoints);
  grid.Reset(nPoints, width, height);
  CollectPoints(grid);
  CreateMesh(grid);
  MakeWhite(grid);
}

// ============================================================================

// Type for the deformation mesh:

struct Deformation
{
  long      nx,      // Point indices in x-direction: 0 ... nx-1
            ny,      // Point indices in y-direction: 0 ... ny-1
            size;    // 2 * nx * ny
  double    width,   // Calculated relative width  of unwarped image
            height;  // Calculated relative height of unwarped image
  deform_t *xyMap;   // Pixel coordinates in skewed image

  Deformation()  { xyMap = 0;  Reset(); }
  ~Deformation() {             Reset(); }

  void Reset(long nx=0, long ny=0);
  void Read(const char *filename, FILE *fh=0);
  void Write(const char *filename, FILE *fh=0);
  void SetFlattenedPictureSize(Grid& grid);
  void Calculate(Grid& grid, int gridMargin, int gridScale,
                 double gridPerInchX, double gridPerInchY);
  void Flatten(Picture& srcPic, Picture& dstPic);
};

// ----------------------------------------------------------------------------

// Reset the internal grid memory:

void Deformation::Reset(long nx, long ny)
{
  if (nx > 1000 || ny > 1000)
    Error("Deformation grid size %ldx%ld is too high", nx, ny);
  if (xyMap)  delete xyMap;
  this->nx = nx;
  this->ny = ny;
  size     = 2 * nx * ny;
  width    = height = 0.0;
  xyMap    = 0;
  if (size) {
    xyMap = new deform_t[size];
    for (long i=0; i < size; ++i)  xyMap[i] = 0; }
}

// ----------------------------------------------------------------------------

// Read a deformation map file:

void Deformation::Read(const char *filename, FILE *fh)
{
  u_int32_t magic;
  if (fh != stdin)  fh = fopen(filename, "rb");
  else              binmode(fh);  // Windows.
  if (!fh ||
      fread(&magic, sizeof(u_int32_t), 1, fh) != 1 ||
      magic != 0x7f3319bc ||
      fread(&width, sizeof(double), 1, fh) != 1 ||
      width < 0.0 || width > 4.0 ||
      fread(&height, sizeof(double), 1, fh) != 1 ||
      height < 0.0 || height > 4.0 ||
      fread(&nx, sizeof(long), 1, fh) != 1 ||
      nx < 4 || nx > 1000 ||
      fread(&ny, sizeof(long), 1, fh) != 1 ||
      ny < 4 || ny > 1000)  magic = 0;
  if (magic) {
    size  = 2 * nx * ny;
    xyMap = new deform_t[size];
    if (fread(xyMap,sizeof(deform_t),size,fh) != ulong(size))  magic = 0; }
  if (fh && fh != stdin)  fclose(fh);
  if (!magic)  Error("Couldn't read deformation map file '%s'", filename);
}

// ----------------------------------------------------------------------------

// Write a deformation map file:

void Deformation::Write(const char *filename, FILE *fh)
{
  u_int32_t magic = 0x7f3319bc;
  if (fh != stdout)  fh = fopen(filename,"wb");
  else               binmode(fh);  // Windows.
  if (!fh ||
      fwrite(&magic, sizeof(u_int32_t), 1, fh) != 1 ||
      fwrite(&width, sizeof(double), 1, fh) != 1 ||
      fwrite(&height, sizeof(double), 1, fh) != 1 ||
      fwrite(&nx, sizeof(long), 1, fh) != 1 ||
      fwrite(&ny, sizeof(long), 1, fh) != 1 ||
      fwrite(xyMap, sizeof(deform_t), size, fh) != ulong(size))  magic = 0;
  if (fh && fh != stdout)  fclose(fh);
  if (!magic) Error("Couldn't write deformation map file '%s'", filename);
}

// ----------------------------------------------------------------------------

// Calculate a suitable image size 'width * height' for the unwarped image:

void Deformation::SetFlattenedPictureSize(Grid& grid)
{
  double    sx   = double(nx - 1) / double(grid.width  - 1),
            sy   = double(ny - 1) / double(grid.height - 1),
            area = 0.0;
  deform_t *map  = xyMap;
  long      x=0, y=0, x1=map[0], y1=map[1], x0,y0;
  int    mode = 0;
  while (mode < 4) {
    x0 = x1;  y0 = y1;
    switch (mode) {
      case 0:   map += 2;     if (++x >= nx-1) ++mode;  break;
      case 1:   map += 2*nx;  if (++y >= ny-1) ++mode;  break;
      case 2:   map -= 2;     if (--x <= 0)    ++mode;  break;
      default:  map -= 2*nx;  if (--y <= 0)    ++mode; }
    x1 = map[0];  y1 = map[1];
    area += 0.5 * double(x0 + x1) * double(y1 - y0); }
  area   = fabs(area / (sx * sy));
  area   = sqrt(area) / double(0x10000000);
  width  = area * sx;
  height = area * sy;
  if (width < 1e-6 || height < 1e-6)  throw "Unwarping area is too small";
  if (width > 4.0  || height > 4.0)   throw "Unwarping area is too big";
  Print("Deskewed picture size: %ld x %ld   (%.2lf%% x %.2lf%%)\n",
        long(double(grid.width  - 1)  * width),
        long(double(grid.height - 1) * height),
        100.0 * width, 100.0 * height);
}

// ----------------------------------------------------------------------------

// Calculate a complete deformation map grid of suitable resolution
// with interpolated and extrapolated mesh points (by MLS):

void Deformation::Calculate(Grid& grid, int gridMargin, int gridScale,
                            double gridPerInchX, double gridPerInchY)
{
  long      gnx,gny, ix,iy, j,k;
  double    sx,sy;
  deform_t *map;
  if (grid.width <=1 || grid.height <= 1)
    throw "Calibration grid has too few pixel";
  sx = 1.0 / double(grid.width  - 1);
  sy = 1.0 / double(grid.height - 1);
  grid.Limit(gnx, gny);
  Reset(gridScale*(gnx-1+2*gridMargin)+1, gridScale*(gny-1+2*gridMargin)+1);
  if (nx < 4)  throw "Too less grid points in x-direction";
  if (ny < 4)  throw "Too less grid points in y-direction";
  // Moving least square Verfahren (this also smoothes the data a little bit):
  for (iy=0,map=xyMap; iy < ny; ++iy) {
    double y = double(iy) / gridScale - gridMargin;
    for (ix=0; ix < nx; ++ix) {
      double x = double(ix) / gridScale - gridMargin;
      double a[6]  = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
             bx[3] = {0.0, 0.0, 0.0},
             by[3] = {0.0, 0.0, 0.0}, alpha[3],detA,xx,yy;
      for (j=0; j < grid.n; ++j) {
        double dx      = x - grid.ix[j],
               dy      = y - grid.iy[j],
               omega   = exp(-0.5 * (dx * dx + dy * dy)),
               pol[3]  = { 1.0, dx, dy },
               wpol[3] = { omega, omega * dx, omega * dy },
               fx      = grid.x[j],
               fy      = grid.y[j];
        for (k=0; k < 3; ++k) {
          a[k]   += wpol[k] * pol[k];
          a[k+3] += wpol[k] * pol[(k+1)%3];
          bx[k]  += fx * wpol[k];
          by[k]  += fy * wpol[k]; } }
      // Solve A*alpha = [detA;0;0], A = {{a0,a3,a5},{a3,a1,a4},{a5,a4,a2}}:
      alpha[0]  = a[1] * a[2] - a[4] * a[4];
      alpha[1]  = -a[2] * a[3];
      detA      = a[3] * a[4];
      alpha[2]  = detA - a[1] * a[5];
      detA      = a[0] * alpha[0] + a[3] * alpha[1] + a[5] * (detA + alpha[2]);
      alpha[1] += a[4] * a[5];
      if (fabs(detA) <= 1e-12)  throw "Distance too large for extrapolation";
      for (k=0,xx=yy=0.0; k < 3; ++k) {
        xx += bx[k] * alpha[k];
        yy += by[k] * alpha[k]; }
      xx *= sx / detA;
      yy *= sy / detA;
      if (fabs(xx) > 2.0 || fabs(yy) > 2.0)  throw "Distortion is too large";
      *map++ = deform_t(round(double(0x10000000) * xx));
      *map++ = deform_t(round(double(0x10000000) * yy)); } }
  SetFlattenedPictureSize(grid);
  // Calculate PPI information:
  sx  = double(grid.width  - 1) * width;
  sy  = double(grid.height - 1) * height;
  sx *= gridPerInchX * double(gridScale) / double(nx - 1);
  sy *= gridPerInchY * double(gridScale) / double(ny - 1);
  Print("PPI in xy-direction:  %ld x %ld\n",
        long(round(sx)), long(round(sy)));
  Print("PPI: %ld\n", long(round(sqrt(sx * sy))));
}

// ----------------------------------------------------------------------------

// Unwarp image 'picSrc' onto 'picDst':

void Deformation::Flatten(Picture& srcPic, Picture& dstPic)
{
  long      srcMX  = srcPic.width  - 1,
            srcMY  = srcPic.height - 1, srcOY = srcPic.width * 3,
            defMX  = nx            - 1,
            defMY  = ny            - 1, defOY = nx * 2,
            dstMX  = long(width  * double(srcMX)), dstX,
            dstMY  = long(height * double(srcMY)), dstY, i;
  double    dstSX = double(defMX) / double(dstMX),
            dstSY = double(defMY) / double(dstMY),
            srcSX = double(srcMX) / double(0x10000000),
            srcSY = double(srcMY) / double(0x10000000), tf;
  long      srcDepth = srcPic.depth,
            dstDepth = param.outputDepth ? param.outputDepth : srcDepth;
  double    outDepth = double(dstDepth),
            facDepth = outDepth / double(srcDepth);
  pixel_t  *srcRaw   = srcPic.pixel, *dstRaw, *rawX, *rawY;
  deform_t *mapX, *mapY;
  Cubic2<deform_t> defCub(defMX, 2, defMY, defOY);
  Cubic2<pixel_t>  srcCub(srcMX, 3, srcMY, srcOY);
  dstPic.Reset(dstMX+1, dstMY+1, dstDepth);
  // Empty all pixel values with color white as default value:
  i             = dstPic.size;
  dstRaw = rawX = dstPic.pixel;
  while (i--)  *rawX++ = dstDepth;
  for (dstY=0; dstY <= dstMY; ++dstY) {
    mapY = xyMap + defOY * defCub.cy(dstSY * dstY);
    for (dstX=0; dstX <= dstMX; ++dstX) {
      mapX = mapY + 2 * defCub.cx(dstSX * dstX);
      tf   = srcSX * defCub(mapX);
      if (tf >= 0.0 && tf <= double(srcMX)) {
        rawX = srcRaw + 3 * srcCub.cx(tf);
        tf   = srcSY * defCub(mapX+1);
        if (tf >= 0.0 && tf <= double(srcMY)) {
          rawY = rawX + srcOY * srcCub.cy(tf);
          for (i=0; i < 3; ++i) {
            tf = round(facDepth * srcCub(rawY++));
            if      (tf < 0.0)       tf = 0.0;
            else if (tf > outDepth)  tf = outDepth;
            *dstRaw++ = pixel_t(tf); } } } } }
}

// ============================================================================

// Main routine:

void Main()
{
  int         hasData  = 0;
  Deformation deform;
  Grid        grid;
  if (param.calibPicName) {
    int     hasCheck = (param.calibCheckName != 0);
    Picture pic;
    if (*param.calibPicName)  pic.Read(param.calibPicName);
    else                      pic.Read("stdin",stdin);
    try {
      pic.Calibrate(grid);
      if (hasCheck) {
        hasCheck = 0;
        pic.Write(param.calibCheckName); } }
    catch(...) {
      if (hasCheck)  pic.Write(param.calibCheckName);
      throw; }
    if (param.calibTextName)  grid.Write(param.calibTextName);
    hasData = 1; }
  if (!hasData && param.calibTextName) {
    grid.Read(param.calibTextName);
    hasData = 1; }
  if (param.calibDeformName || param.sourcePicName) {
    if (hasData) {
      deform.Calculate(grid, param.gridMargin, param.gridScale,
                       param.gridPerInchX, param.gridPerInchY);
      if (param.calibDeformName) {
        if (*param.calibDeformName) deform.Write(param.calibDeformName);
        else                        deform.Write("stdout", stdout); } }
    else if (param.calibDeformName && *param.calibDeformName) {
      deform.Read(param.calibDeformName);
      hasData = 1; }
    if (hasData && param.sourcePicName && param.destPicName) {
      Picture src, dst;
      if (*param.sourcePicName) src.Read(param.sourcePicName);
      else                      src.Read("stdin", stdin);
      deform.Flatten(src, dst);
      if (*param.destPicName)  dst.Write(param.destPicName);
      else                     dst.Write("stdout", stdout); } }
}

// ============================================================================

// Main program to supply command line args, catch errors and
// call the main routine:

int main(int argc,char *argv[]) {
  try {
    try {
      param(argc, argv);
      Main(); }
    // Catch exception from 'new' operator:
    catch (std::bad_alloc& ba) { throw "Out of memory"; } }
  catch (const char *msg) {
    fprintf(stderr,"!!! Error in %s:\n!!! %s!\n", param.prgName, msg);
    return 1; }
  return 0;
}
