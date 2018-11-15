// ppmwhitebalance.cc
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

// This program can be used to correct the white balance of an image
// based on a calculated grid of correction values, which is derived
// from a calibrating photo made from a gray card.
// The primary goal is to normalize the colors of pictures from book pages
// or other material using a digital camera instead of a scanner.
// The setup for several images using the same calibration shot should
// remain as fixed as possible (lights, shadows, camera settings).
// All pictures read and written by this program must be in raw PPM format.
// The program is testet with Linux but should also run on other platforms
// with only minor modifications. No special libraries are used.
// To compile use something like:
//
// $ g++ -o ppmwhitebalance ppmwhitebalance.cc
//
// A brief help text is printed, when wrong command syntax is used or simply
// by typing:
//
// $ ppmwhitebalance -h
//
// The user uses a gray card fitting the whole area of the image, e.g.
// the gray card should cover the picture area of interest.
// The user puts this gray card at the same position as the pages to be
// color corrected and takes a digital photo.
// By using conversion tools like 'jpegtopnm' she creates a file,
// e.g. 'calibration.ppm'. From this file, the command
//
// $ ppmwhitebalance calibration.ppm > calibration.bin
//
// extracts the information needed to normalize the color of the book sheets.
// To adjust the detection process, the user can use the following options:
//
// -gb        Correct only brightness values not individual RGB-channels.
// -g1        Maximal interpolation order is linear.
// -gc  80    Changes the desired RGB gray value of the gray card to 80.
//            Take mean value of calibration picture if zero.
// -gm 4.0    Maximal Multiplicator for the RGB channels.
// -gd 4.0    Maximal Divisor       for the RGB channels.
// -nx 10     Number of grid points in x-direction (default 1 if nx=ny=0).
// -ny 10     Number of grid points in y-direction (default 1 if nx=ny=0).
//
// If one of the last two options is zero, that value is choosen
// approximately considering the aspect ratio of the image.
//
// If 'page01.ppm' is one of the images of the camera to be color
// adjusted, this can be achieved by:
//
// $ ppmwhitebalance -d calibration.bin page01.ppm > page01_enhanced.ppm
//
// This program can handle 8 bit and 16 bit PPM files with channel depth
// values (maximal color value) from 1 to 65535. To specify an output color
// depth different from the input picture, use the option '-od'.
//
// If a book is digitized page after page, and the environmental conditions
// change, a new calibration picture shot should be taken from time to time.
//
// If the border colors of the calibration image 'calibration.ppm' changes
// rapidly, the interpolation in the neighborhood can be disturbed.
// Either change the calibration image manually by a graphic editor,
// cut the problematic errors in all pictures before color calibration or
// avoid such a setup completely.

// ============================================================================

// Include files; no other libraries needed than 'libgpp' and 'ppmroselib.o':

#include "ppmroselib.h"

// ============================================================================

// We limit the number of grid points, otherwise the
// band matrix solver will become too expensive:

#define MAX_DOF  100000

// ============================================================================

// Howto obtain gray values:

#define PTR_TO_GRAY(value,ptr) \
  value += 0.299 * double(*ptr++); \
  value += 0.587 * double(*ptr++); \
  value += 0.114 * double(*ptr++)

// ============================================================================
// Global parameters/options.
// ============================================================================

struct ParameterWhiteBalance : Parameter
{
  // Calibration mode enforced:
  int doCalib;

  // File names:
  const char *calibPicName;    // Input image with gray card.
  const char *calibTextName;   // Calibration grid in text format.
  const char *calibColorName;  // Final color grid.
  const char *sourcePicName;   // Skewed input image.
  const char *destPicName;     // Unwarped output image.

  // Should adjustment only act on brightness:
  int colorBrightness;

  // Maximal linear interpolation:
  int gridLinear;

  // Desired gray RGB value of gray card (16 <= grayCard <= 240):
  int grayCard;

  // Maximal multiplicator for RGB values:
  double colorFactor;

  // Maximal divisor for RGB values:
  double colorDivisor;

  // Number of grid points in x-direction:
  int gridNX;

  // Number of grid points in y-direction:
  int gridNY;

  // Output color depth per channel:
  int outputDepth;

  ParameterWhiteBalance() : Parameter("ppmwhitebalance") { }

  void Define();
  void Check();
};

// ----------------------------------------------------------------------------

// Program parameters are global for ease of use:

static ParameterWhiteBalance param;

int main(int argc, char *argv[]) { return Handler(argc, argv, param); }

// ----------------------------------------------------------------------------

void ParameterWhiteBalance::Define()
{
  AddFlag("-c",  doCalib, 0,
          "Enforce calibration mode");
  AddString("-cc", calibPicName, 0, "(inpname)",
            "Set input PPM picture with gray card calibration image");
  AddString("-cp", calibTextName, 0, "<name>",
            "Set file name for textual color grid");
  AddString("-d", calibColorName, 0, "(stdout)",
            "Set file name for binary color grid");
  AddString("-i", sourcePicName, 0, "(inpname)",
            "Set input file name");
  AddString("-o", destPicName, 0, "(stdout)",
            "Set output file name");
  AddFlag("-gb", colorBrightness, 0,
          "Only change brightness of color");
  AddFlag("-g1", gridLinear, 0,
          "Maximum order of interpolation is linear");
  AddInt("-gc", grayCard, 10, 128, 0, 240, 0,
         "Set the desired RGB gray card value.\n"
         "A zero value indicates picture mean gray value");
  AddDouble("-gm", colorFactor, 1.5, 1.0, 16.0, 0,
            "Set maximal multiplicator for RGB values");
  AddDouble("-gd", colorDivisor, 1.5, 1.0, 16.0, 0,
            "Set maximal divisor       for RGB values");
  AddInt("-nx", gridNX, 10, 0, 1, 1000, 0,
         "Number of grid points in x-direction (1 if nx=ny=0)");
  AddInt("-ny", gridNY, 10, 0, 1, 1000, 0,
         "Number of grid points in y-direction (1 if nx=ny=0).\n"
         "If nx or ny is zero, the image aspect ratio is used");
  AddInt("-od", outputDepth, 10, 0, 0, 65535, 0,
         "Output color depth, zero means same as input depth");
  ExtraUsage(
    "Simple calibration:  ppwhitebalance calib.ppm > color.bin\n"
    "Simple usage:        ppwhitebalance -d color.bin src.ppm > dst.ppm");
}

// ----------------------------------------------------------------------------

// Special argument checking:

void ParameterWhiteBalance::Check()
{
  int n = 0;
  if (grayCard && grayCard < 16)  Usage(1);
  if (calibPicName)    ++n;
  if (calibTextName)   ++n;
  if (calibColorName)  ++n;
  if (n != 1)  doCalib = 1;
  if (doCalib) {
    if (sourcePicName || destPicName)  Usage(1);
    if (calibPicName) { if (inpName)  Usage(1); }
    else if (inpName) { calibPicName = inpName; }
    if (!calibPicName && !calibTextName) calibPicName = "";
    if (!calibColorName && (!calibPicName || !calibTextName))
      calibColorName = ""; }
  else {
    if (!sourcePicName)  sourcePicName = inpName ? inpName : "";
    else if (inpName)    Usage(1);
    if (!destPicName)    destPicName = ""; }
  if (!gridNX && !gridNY)  gridNX = gridNY = 1;
}

// ============================================================================

// Type for image data:

struct Picture : Image {
  void Write(const char *filename, FILE *fh=0) {
    Image::Write(filename, fh, param.prgName); }

  int  GetMeanValue();
};

// ----------------------------------------------------------------------------

int Picture::GetMeanValue() {
  pixel_t *ptr   = pixel;
  long     nn    = width * height, n = nn;
  double   value = 0.0;
  double   fac   = 255.0 / double(depth);
  while (n--) { PTR_TO_GRAY(value,ptr); }
  value = round(fac * (value / double(nn)));
  if      (value <  16.0)  value =  16.0;
  else if (value > 240.0)  value = 240.0;
  return int(value);
}

// ============================================================================

// Conversion modul:

struct Converter {
  Picture      *pic;
  double       *map;
  BandedSystem  sys;
  long          pn1, po1, pn2, po2;
  long          gn1, gn2, go2, gb;
  int           swap, linear;

  Converter()  { map = 0;  Reset(); }
  ~Converter() {           Reset(); }

  void Reset(Picture* pic=0, long nx=0, long ny=0, int linear=0);
  void SetMap(int channel = -1);
  void BuildtAndSolveSystem(int doMat=0);
  void ReadSolution(color_t *rgbMap);
  void GetRgbMap(color_t *rgbMap);

};

// ----------------------------------------------------------------------------

void Converter::Reset(Picture *pic, long nx, long ny, int linear) {
  long size;
  if (map)  delete map;
  this->pic    = pic;
  this->linear = linear;
  map          = 0;
  pn1          = po1 = pn2 = po2 = gn1 = gn2 = go2 = gb = 0;
  swap         = 0;
  if (pic) {
    pn1 = pic->width;
    po1 = 1;
    pn2 = pic->height;
    po2 = pn1;
    gn1 = nx;
    gn2 = ny;
    go2 = nx;
    if (gn1 > gn2) {
      po1 = pn1;  pn1 = pn2;  pn2 = po1; po1 = po2;  po2 = 1;
      go2 = gn1;  gn1 = gn2;  gn2 = go2; go2 = gn1;
      swap = 1; }
    gb   = (linear ? 1 : 3) * (go2 + 1);
    size = pn1 * pn2;
    map  = new double[size];
    for (long i=0; i < size; ++i)  map[i] = 0.0;
    sys.Reset(gn1 * gn2, gb); }
  else  sys.Reset();
}

// ----------------------------------------------------------------------------

// Creates the map for the RGB-channel (0,1,2) or for gray values (-1):

void Converter::SetMap(int channel) {
  pixel_t *pixel  = pic->pixel;
  double  *ptr    = map;
  long     n      = pn1 * pn2;
  double   fac    = 1.0 / double(pic->depth);
  double   orig;
  while (n--) {
    if (channel < 0) { orig = 0.0;  PTR_TO_GRAY(orig,pixel);       }
    else             { orig = double(pixel[channel]);  pixel += 3; }
    *ptr++ = fac * orig; }
}

// ----------------------------------------------------------------------------

// Creates the band matrix and right hand side for the selected channel:

void Converter::BuildtAndSolveSystem(int doMat) {
  long    gm1    = gn1 - 1;
  long    gm2    = gn2 - 1;
  double  scale1 = pn1 > 1 ? double(gm1) / double(pn1 - 1) : 1.0;
  double  scale2 = pn2 > 1 ? double(gm2) / double(pn2 - 1) : 1.0;
  double *mat    = sys.GetMat();
  double *rhs    = sys.GetRhs();
  double *src1, *src2;
  long    i1, i2;
  Interpolation2<double>  inter(linear, gm1, 1, gm2, go2, 2*gb);
  if (doMat)  sys.ClearMat();
  else        sys.ClearRhs();
  for (i1=0, src1=map; i1 < pn1; ++i1, src1 += po1) {
    inter.cx.Set(scale1 * double(i1));
    for (i2=0, src2=src1; i2 < pn2; ++i2, src2 += po2) {
      inter.cy.Set(scale2 * double(i2));
      if (doMat)  inter.AssembleMat(mat);
      else        inter.AssembleRhs(*src2, rhs); } }
  if (doMat) {
    Print("Cholesky factorization ...\n");
    sys.Cholesky(); }
  else sys.Solve();
}

// ----------------------------------------------------------------------------

// Solve the band matrix system:

void Converter::ReadSolution(color_t *rgbMap) {
  int     brightness = param.colorBrightness;
  long    gnx, gox, gny, goy, ix, iy;
  double *px, *py, z;
  if (swap) { gnx = gn2;  gox = go2;  gny = gn1;  goy = 1;   }
  else      { gnx = gn1;  gox = 1;    gny = gn2;  goy = go2; }
  for (iy=0, py=sys.GetSol(); iy < gny; ++iy, py += goy) {
    for (ix=0, px=py; ix < gnx; ++ix, px += gox) {
      z = *px;
      z = z < 0.0 ? 0.0 : z > 1.0 ? 1.0 : z;
      z = 65535.0 * z + 0.5;
      *rgbMap = color_t(z);
      if (brightness)  rgbMap[2] = rgbMap[1] = rgbMap[0];
      rgbMap += 3; } }
}

// ----------------------------------------------------------------------------

// Combines all actions to obtain color grid:

void Converter::GetRgbMap(color_t *rgbMap) {
  int channel;
  Print("Build system matrix ...\n");
  BuildtAndSolveSystem(1);
  if (param.colorBrightness) {
    // Set factors for brightness only:
    Print("Calibrate brightness ...\n");
    SetMap();
    BuildtAndSolveSystem();
    ReadSolution(rgbMap); }
  else {
    for (channel = 0; channel < 3; ++channel) {
      // Set factors for each pixel of selected channel:
      Print("Calibrate channel %d ...\n", channel);
      SetMap(channel);
      BuildtAndSolveSystem();
      ReadSolution(rgbMap++); } }
}

// ============================================================================

// Type for grid of color points:

struct Grid {
  long       nx,      // Point indices in x-direction: 0 ... nx-1
             ny,      // Point indices in y-direction: 0 ... ny-1
             size;    // 3 * nx * ny
  int        linear;  // Maximum interpolation order is linear
  int        target;  // Target RGB-value for gray map.
  u_int32_t  smin;    // Minimal scaling factor.
  u_int32_t  smax;    // Maximal scaling factor.
  color_t   *rgbMap;  // normalized RGB values.

  Grid()  { rgbMap = 0;  Reset(); }
  ~Grid() {              Reset(); }

  void Reset(long nx=0, long ny=0);
  void TextRead(const char *filename, FILE *fh=0);
  void TextWrite(const char *filename, FILE *fh=0);
  void Read(const char *filename, FILE *fh=0);
  void Write(const char *filename, FILE *fh=0);
  void Calibrate(Picture& pic);
  void Convert(Picture& src, Picture& dst);
};

// ----------------------------------------------------------------------------

void Grid::Reset(long nx, long ny) {
  size = 3 * nx * ny;
  if (nx > 1000 || ny > 1000 || size > (3*MAX_DOF))
    Error("Color grid size %ldx%ld is too high", nx, ny);
  if (rgbMap)  delete rgbMap;
  this->nx = nx;
  this->ny = ny;
  rgbMap   = 0;
  linear   = 0;
  target   = 0;
  smin     = smax = 0.0;
  if (size) {
    rgbMap = new color_t[size];
    for (long i=0; i < size; ++i)  rgbMap[i] = 0; }
}

// ----------------------------------------------------------------------------

// Reads the color calibration grid from text file:

void Grid::TextRead(const char *filename, FILE *fh) {
  int      z, zz, r, g, b;
  long     k, n, ix, iy;
  color_t *rgb;
  if (fh != stdin)  fh = fopen(filename, "r");
  if (!fh)  Error("Couldn't read color calibration file '%s'", filename);
  try {
    if (fscanf(fh, "Colorcalibrationgrid %ld %ld\n", &ix, &iy) != 2 ||
        ix <= 0 || iy <= 0 || ix > 1000 || iy > 1000 ||
        ix * iy > MAX_DOF)  throw 0;
    Reset(ix, iy);
    if (fscanf(fh, "Onlylinear %d\n", &linear) != 1 ||
        linear < 0 || linear > 1)  throw 0;
    if (fscanf(fh, "Target %d\n", &target) != 1 || target < 16 || target > 240)
      throw 0;
    if (fscanf(fh, "Factor %x %x ", &z, &zz) != 2 ||
        z < 4096 || z > 1048576 || zz < 4096 || zz > 1048576)  throw 0;
    smin = u_int32_t(z);
    smax = u_int32_t(zz);
    while (!feof(fh))  if (fgetc(fh) == '\n')  break;
    rgb = rgbMap;
    n   = nx * ny;
    for (k=0; k < n; ++k) {
      if (fscanf(fh, "%ld\t%ld\t%x\t%x\t%x ", &ix, &iy, &r, &g, &b) != 5 ||
          ix < 0 || ix >= nx || iy < 0 || iy >= ny ||
          r < 0 || r > 65535 ||
          g < 0 || g > 65535 ||
          b < 0 || b > 65535)  throw 1;
      *rgb++ = color_t(r);
      *rgb++ = color_t(g);
      *rgb++ = color_t(b);
      while (!feof(fh))  if (fgetc(fh) == '\n')  break; }
    if (fh != stdin)  fclose(fh); }
  catch (int nmr) {
    const char *msg = nmr ? "has wrong point syntax" : "has wrong preamble";
    if (fh != stdin)  fclose(fh);
    Error("Calibration color file '%s' %s", filename, msg); }
}

// ----------------------------------------------------------------------------

// Writes the color calibration grid to text file:

void Grid::TextWrite(const char *filename, FILE *fh) {
  color_t *rgb = rgbMap;
  if (fh != stdout)  fh = fopen(filename, "w");
  if (!fh)  Error("Couldn't write color calibration file '%s'", filename);
  fprintf(fh, "Colorcalibrationgrid %ld %ld\n", nx, ny);
  fprintf(fh, "Onlylinear %d\n", linear);
  fprintf(fh, "Target %d\n", target);
  fprintf(fh, "Factor %06x %06x   (%2.2lf  %2.2lf)\n",
          smin, smax, double(smin) / 65536.0, double(smax) / 65536.0);
  for (long iy=0; iy < ny; ++iy) {
    for (long ix=0; ix < nx; ++ix, rgb += 3) {
      fprintf(fh,"%ld\t%ld\t%04x\t%04x\t%04x   (%2.2lf  %2.2lf  %2.2lf)\n",
              ix, iy, rgb[0], rgb[1], rgb[2],
              double(rgb[0]) / 65535.0,
              double(rgb[1]) / 65535.0,
              double(rgb[2]) / 65535.0); } }
  if (fh != stdout)  fclose(fh);
}

// ----------------------------------------------------------------------------

// Reads the color calibration grid from binary file:

void Grid::Read(const char *filename, FILE *fh) {
  u_int32_t magic;
  if (fh != stdin)  fh = fopen(filename, "rb");
  else              binmode(fh);  // Windows.
  if (!fh ||
      fread(&magic, sizeof(u_int32_t), 1, fh) != 1 ||
      magic != 0x7f4418bd ||
      fread(&nx, sizeof(long), 1, fh) != 1 ||
      nx < 1 || nx > 1000 ||
      fread(&ny, sizeof(long), 1, fh) != 1 ||
      ny < 1 || ny > 1000 ||
      fread(&linear, sizeof(int), 1, fh) != 1 ||
      linear < 0 || linear > 1 ||
      fread(&target, sizeof(int), 1, fh) != 1 ||
      target < 16 || target > 240 ||
      fread(&smin, sizeof(u_int32_t), 1, fh) != 1 ||
      smin < 4096 || smin > 1048576 ||
      fread(&smax, sizeof(u_int32_t), 1, fh) != 1 ||
      smax < 4096 || smax > 1048576 ||
      smin > smax)  magic = 0;
  if (magic) {
    size = 3 * nx * ny;
    if (size <= (3*MAX_DOF)) {
      rgbMap = new color_t[size];
      if (fread(rgbMap, sizeof(color_t), size, fh) != ulong(size))
        magic = 0; }
    else magic = 0; }
  if (fh && fh != stdin)  fclose(fh);
  if (!magic) Error("Couldn't read color grid file '%s'", filename);
}

// ----------------------------------------------------------------------------

// Writes the color calibration grid to binary file:

void Grid::Write(const char *filename, FILE *fh) {
  u_int32_t magic = 0x7f4418bd;
  if (fh != stdout)  fh = fopen(filename, "wb");
  else               binmode(fh);  // Windows.
  if (!fh ||
      fwrite(&magic,  sizeof(u_int32_t), 1, fh) != 1 ||
      fwrite(&nx,     sizeof(long),      1, fh) != 1 ||
      fwrite(&ny,     sizeof(long),      1, fh) != 1 ||
      fwrite(&linear, sizeof(int),       1, fh) != 1 ||
      fwrite(&target, sizeof(int),       1, fh) != 1 ||
      fwrite(&smin,   sizeof(u_int32_t), 1, fh) != 1 ||
      fwrite(&smax,   sizeof(u_int32_t), 1, fh) != 1 ||
      fwrite(rgbMap,  sizeof(color_t), size, fh) != ulong(size))  magic = 0;
  if (fh && fh != stdout)  fclose(fh);
  if (!magic) Error("Couldn't write color grid file '%s'", filename);
}

// ----------------------------------------------------------------------------

// Create color calibration from gray card image:

void Grid::Calibrate(Picture& pic) {
  long nx = param.gridNX, ny = param.gridNY;
  Converter converter;
  if      (nx <= 0)  nx = (ny * pic.width + (pic.height >> 1)) / pic.height;
  else if (ny <= 0)  ny = (nx * pic.height + (pic.width >> 1)) / pic.width;
  if (nx <= 0)  nx = 1;
  if (ny <= 0)  ny = 1;
  if (nx > pic.width)   nx = pic.width;
  if (ny > pic.height)  ny = pic.height;
  Print("Calibration grid: %ld x %ld\n", nx, ny);
  Reset(nx, ny);
  if (!param.grayCard) {
    param.grayCard = pic.GetMeanValue();
    Print("Automatic graycard value: %d\n", param.grayCard); }
  // Save important parameters for later conversion:
  linear = param.gridLinear;
  target = param.grayCard;
  smin   = u_int32_t(65536.0 / param.colorDivisor + 0.5);
  smax   = u_int32_t(65536.0 * param.colorFactor  + 0.5);
  // Calculate the best grid:
  converter.Reset(&pic, nx, ny, linear);
  converter.GetRgbMap(rgbMap);
}

// ----------------------------------------------------------------------------

// Applies the color enhancement to image:

void Grid::Convert(Picture& src, Picture& dst) {
  long     width    = src.width;
  long     height   = src.height;
  long     mx       = nx - 1;
  long     my       = ny - 1;
  long     oy       = 3 * nx, ix, iy;
  double   scalex   = width  > 1 ? double(mx) / double(width  - 1) : 1.0;
  double   scaley   = height > 1 ? double(my) / double(height - 1) : 1.0;
  double   gray     = double(target) /   255.0;
  double   div      = double(smin)   / 65536.0;
  double   fac      = double(smax)   / 65536.0;
  long     srcDepth = src.depth;
  long     dstDepth = param.outputDepth ? param.outputDepth : srcDepth;
  double   outDepth = double(dstDepth);
  double   facDepth = outDepth / double(srcDepth), z;
  int      channel;
  pixel_t *ps, *pd;
  Interpolation2<color_t> inter(linear, mx, 3, my, oy);
  dst.Reset(width, height, dstDepth);
  for (iy=0, ps=src.pixel, pd=dst.pixel; iy < height; ++iy) {
    inter.cy.Set(scaley * double(iy));
    for (ix=0; ix < width; ++ix) {
      inter.cx.Set(scalex * double(ix));
      for (channel=0; channel < 3; ++channel) {
        z  = inter(rgbMap + channel) / 65535.0;
        z  = fac * z <= gray ? fac : div * z >= gray ? div : gray / z;
        z *= facDepth * double(*ps++);
        if      (z < 0.0)       z = 0.0;
        else if (z > outDepth)  z = outDepth;
        *pd++ = pixel_t(z); } } }
}

// ============================================================================

// Main routine:

void Main() {
  int     hasData = 0;
  Grid    grid;
  if (param.calibPicName) {
    Picture pic;
    if (*param.calibPicName)  pic.Read(param.calibPicName);
    else                      pic.Read("stdin", stdin);
    grid.Calibrate(pic);
    if (param.calibTextName)  grid.TextWrite(param.calibTextName);
    hasData = 1; }
  if (!hasData && param.calibTextName) {
    grid.TextRead(param.calibTextName);
    hasData = 1; }
  if (param.calibColorName || param.sourcePicName) {
    if (hasData) {
      if (param.calibColorName) {
        if (*param.calibColorName)  grid.Write(param.calibColorName);
        else                        grid.Write("stdout", stdout); } }
    else if (param.calibColorName && *param.calibColorName) {
      grid.Read(param.calibColorName);
      hasData = 1; }
    if (hasData && param.sourcePicName && param.destPicName) {
      Picture src, dst;
      if (*param.sourcePicName)  src.Read(param.sourcePicName);
      else                       src.Read("stdin", stdin);
      grid.Convert(src, dst);
      if (*param.destPicName)  dst.Write(param.destPicName);
      else                     dst.Write("stdout", stdout); } }
}
