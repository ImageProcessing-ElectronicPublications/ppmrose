// ppmroselib.cc
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

#include "ppmroselib.h"

// ============================================================================
// General common functions.
// ============================================================================

static int printFlag = 1;

// Error message with program termination:

#define ELEN 1024
void Error(const char *format, ...)
{
  static char buf[ELEN];  // Probably enough for error messages.
  buf[ELEN-5] = 'e';
  va_list ap;
  va_start(ap, format);
  vsnprintf(buf, ELEN, format, ap);
  va_end(ap);
  if (buf[ELEN-5])  strcpy(buf+(ELEN-5), " ...");
  throw buf;
}
#undef ELEN

// ----------------------------------------------------------------------------

// Normal messages, which can be suppressed:

void Print(const char *format, ...)
{
  if (printFlag) {
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
} }

// ============================================================================
// Unified handling of global parameters/options.
// ============================================================================

// 'mode == 2' calculates column width of the formatted option list,
// 'mode == 3' prints the formatted option list:

void Parameter::PrintOption()
{
  char *ptr = buf;
  if (mode == 2) {
    if (argi < int(strlen(ptr)))  argi  = strlen(ptr);
    ptr += strlen(ptr) + 1;
    if (argii < int(strlen(ptr)))  argii = strlen(ptr);
    return; }
  // For 'mode == 3':
  fprintf(stderr, "%-*s ", argi,  ptr);  ptr += strlen(ptr) + 1;
  fprintf(stderr, "%-*s ", argii, ptr);  ptr += strlen(ptr) + 1;
  for (;;) {
    while (*ptr && *ptr != '\n')  fputc(*ptr++, stderr);
    if (!*ptr++)  break;
    fprintf(stderr, "\n%-*s %-*s ", argi, "", argii, ""); }
  fprintf(stderr, ".\n");
}

// ----------------------------------------------------------------------------

// Check if current option is equal to 'option' and get next argument
// if 'extra' is set. On success return 1:

int Parameter::GetOption(const char *option, int extra)
{
  if (argi != argii || strcmp(arg, option))  return 0;
  ++argi;
  if (extra) {
    if (argi >= argc)  Usage(1);
    arg = argv[argi++]; }
  return 1;
}

// ----------------------------------------------------------------------------

int Parameter::AddFlag(const char *option, int& param,
                       const char *show, const char *meaning)
{
  switch (mode) {
    case 0:   param = 0;                          break;
    case 1:   Print("%s = %d\n", option, param);  break;
    case 2:
    case 3:   sprintf(buf, "%s%c%s%c%s",
                      option, 0, (show?show:""), 0, meaning);
              PrintOption();  break;
    default:  if (GetOption(option, 0))  return (param = 1); }
  return 0;
}

// ----------------------------------------------------------------------------

int Parameter::AddString(const char *option, const char *& param,
                         const char *value,
                         const char *show, const char *meaning)
{
  switch (mode) {
    case 0:   param = value;                          break;
    case 1:   Print("%s = \"%s\"\n", option, param);  break;
    case 2:
    case 3:   sprintf(buf, "%s%c%s%c%s", 
                      option, 0, (show?show:value), 0, meaning);
              PrintOption();  break;
    default:  if (GetOption(option, 1)) { param = arg;  return 1; } }
  return 0;
}

// ----------------------------------------------------------------------------

int Parameter::AddInt(const char *option, int& param, int radix,
                      int value, int low, int high,
                      const char *show, const char *meaning)
{
  char *end;
  switch (mode) {
    case 0:   param = value;                      break;
    case 1:   Print("%s = %d\n", option, param);  break;
    case 2:
    case 3:   if (show)
                sprintf(buf, "%s%c%s%c%s", option, 0, show, 0, meaning);
              else if (radix == 16)
                sprintf(buf, "%s%c(%x)%c%s", option, 0, value, 0, meaning);
              else
                sprintf(buf, "%s%c(%d)%c%s", option, 0, value, 0, meaning);
              PrintOption();  break;
    default:  if (GetOption(option, 1)) {
                param = strtol(arg, &end, radix);
                if (*end || param < low || param > high)  Usage(1);
                return 1; } }
  return 0;
}

// ----------------------------------------------------------------------------

int Parameter::AddDouble(const char *option, double& param,
                         double value, double low, double high,
                         const char *show, const char *meaning)
{
  char *end;
  switch (mode) {
    case 0:   param = value;                          break;
    case 1:   Print("%s = %.16lg\n", option, param);  break;
    case 2:
    case 3:   if (show)
                sprintf(buf, "%s%c%s%c%s", option, 0, show, 0, meaning);
              else
                sprintf(buf, "%s%c(%.4lg)%c%s", option, 0, value, 0, meaning);
              PrintOption();  break;
    default:  if (GetOption(option, 1)) {
                param = strtod(arg, &end);
                // Accept rounding errors at boundary:
                if (param < low  && param >= low  - 1e-14)  param = low;
                if (param > high && param <= high + 1e-14)  param = high;
                if (*end || param < low || param > high)  Usage(1);
                return 1; } }
  return 0;
}

// ----------------------------------------------------------------------------

void Parameter::ExtraUsage(const char *msg)
{
  if (mode == 3)  fprintf(stderr, "\n%s\n", msg);
}

// ----------------------------------------------------------------------------

// Displays all options for debugging purpose:

void Parameter::Debug()
{
  Print("BEG Debugging output of options:\n");
  Print("Program name = \"%s\"\n", prgName);
  Print("Quiet flag = %d\n", quiet);
  Print("Input name = \"%s\"\n", inpName);
  mode = 1;
  Define();
  Print("END Debugging output of options.\n");
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
            prgName);
    exit(0); }
  fprintf(stderr,
          "Usage: %s [options] [--] [inpname or stdin]\n\nOptions:\n",
          prgName);
  for (mode = 2, argi=argii=0; mode < 4; ++mode) {
    sprintf(buf, "--version%c%cPrint program version", 0, 0);
    PrintOption();
    sprintf(buf, "-h%c%cPrint program usage", 0, 0);
    PrintOption();
    sprintf(buf, "-q%c%cSuppress normal program messages", 0, 0);
    PrintOption();
    Define(); }
  exit(exitCode);
}

// ----------------------------------------------------------------------------

// Sets global parameters/options from the command line:

void Parameter::operator()(int argcPrg, char *argvPrg[])
{
  int doOpt = 1;
  inpName = 0;
  argi    = 1;
  argc    = argcPrg;
  argv    = argvPrg;
  mode    = quiet = 0;
  // Set default values:
  Define();
  while (argi < argc) {
    arg = argv[argi];
    if (!*arg)  Usage(1);
    if (doOpt && *arg == '-') {
      if      (!strcmp(arg, "--"))         { doOpt = 0;  ++argi; }
      else if (!strcmp(arg, "-q"))         { quiet = 1;  ++argi; }
      else if (!strcmp(arg, "--version"))  Usage(-1);
      else if (!strcmp(arg, "-h"))         Usage(0);
      else {
        mode  = 4;
        argii = argi;
        Define();
        if (argi == argii)  Usage(1); } }
    else {
      if (++argi != argc || !*arg)  Usage(1);
      inpName = arg; } }
  printFlag = !quiet;
  // Additional checks and settings:
  Check();
}

// ============================================================================
// Basic handling for raw PPM-images.
// ============================================================================

// Resizes the image array by given size parameters:

void Image::Reset(long width, long height, long depth)
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

void Image::Read(const char *filename, FILE *fh)
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
          if (mode & 32)     continue;  // Ignore comment characters in line.
          if (!(mode & 16))  throw 1;   // No separation found?
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
    if (!isspace(ch) || !width || !height || !depth)  throw 1;
    if (width > 25000 || height > 25000 || width*height > 100000000)  throw 4;
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
          if (val > depth)  throw 5; } } }
    else {
      while (n) {
        n -= m = n < 4096 ? n : 4096;
        if (fread((src = buf), sizeof(u_int8_t), m, fh) != ulong(m))
          throw 2;
        while (m--) {
          *dst++ = val = pixel_t(*src++);
          if (val > depth)  throw 5; } } }
    if (fgetc(fh) != EOF)  throw 3;
    if (fh != stdin)  fclose(fh); }
  catch (int nmr) {
    const char *msg;
    if (fh != stdin)  fclose(fh);
    switch (nmr) {
      case 1:   msg = "has wrong preamble";           break;
      case 2:   msg = "is too small";                 break;
      case 3:   msg = "is too big";                   break;
      case 4:   msg = "has too large dimensions";     break;
      default:  msg = "has inconsistent color depth"; }
    Error("PPM image file '%s' %s", filename, msg); }
}

// ----------------------------------------------------------------------------

// Write raw PPM image file:

void Image::Write(const char *filename, FILE *fh, const char *creator)
{
  long     n, m, k;
  u_int8_t buf[4096], *dst;
  pixel_t  val,       *src;
  if (fh != stdout)  fh = fopen(filename, "wb");
  else               binmode(fh);  // Windows.
  if (!fh)  Error("Couldn't write '%s' as raw PPM image file", filename);
  fprintf(fh, "P6\n");
  if (creator)  fprintf(fh, "# CREATOR: %s\n", creator);
  fprintf(fh, "%ld %ld\n%ld\n", width, height, depth);
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

// ============================================================================
// Constant/linear/cubic interpolation.
// ============================================================================

// Initializes interpolation at location '0.0 <= t1 <= m'.
// The interpolation is C1-continuous if 'order == 3':

void Interpolation1::Set(double t1) {
  double t0 = floor(t1);
  long   r;
  i   = long(t0);
  t1 -= t0;
  if (i >= mm) { --i;  t1 += 1.0; }  // Border case.
  switch (order) {
    case 0:   f[n = r = 0] = 1.0;                   break;
    case 1:   f[r = 0] = 1.0 - t1;  f[n = 1] = t1;  break;
    default:
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
        f[n = 2] = -0.5 * t0 * t1; } }
  i -= r;
}

// ----------------------------------------------------------------------------

// Returns interpolated value for previously specified grid point
// on array 'raw':
//
// The interpolation on a rectangular grid has
// the tensor product property.

template<class ityp>
double Interpolation2<ityp>::operator()(ityp *raw) {
  double  value = 0.0;
  long    ix,iy;
  ityp   *p;
  raw += Offset();
  for (ix=0; ix <= cx.n; ++ix, raw += cx.o) {
    for (iy=0, p=raw; iy <= cy.n; ++iy, p += cy.o) {
      value += double(*p) * (cx.f[ix] * cy.f[iy]); } }
  return value;
}

// ----------------------------------------------------------------------------

// Assembles right hand side vector "G' * rhs":

template <class ityp>
void Interpolation2<ityp>::AssembleRhs(double scale, double *rhs) {
  long   ix, iy;
  double *p, v;
  rhs += Offset();
  for (ix=0; ix <= cx.n; ++ix, rhs += cx.o) {
    v = scale * cx.f[ix];
    for (iy=0, p=rhs; iy <= cy.n; ++iy, p += cy.o) {
      *p += v * cy.f[iy]; } }
}

// ----------------------------------------------------------------------------

// Assembles full or band system matrix "G' * G":

template <class ityp>
void Interpolation2<ityp>::AssembleMat(double *mat) {
  long   ix, iy, ix2, iy2;
  double *p1, *p2, *p3, v1, v2, v3;
  mat += Offset() * (oo + 1);
  for (ix=0; ix <= cx.n; ++ix, mat += cx.o) {
    v1 = cx.f[ix];
    for (iy=0, p1=mat; iy <= cy.n; ++iy, p1 += cy.o) {
      v2   = v1 * cy.f[iy];
      for (ix2=0, p2=p1; ix2 <= cx.n; ++ix2, p2 += cx.o2) {
        v3 = v2 * cx.f[ix2];
        for (iy2=0, p3=p2; iy2 <= cy.n; ++iy2, p3 += cy.o2) {
          *p3 += v3 * cy.f[iy2]; } } } }
}

// ============================================================================
// Solver for symmetric positive definite systems with limited bandwidth.
// ============================================================================

void BandedSystem::Reset(long n, long b) {
  long    size;
  double *arr;
  if (sol)  delete sol;
  this->n = n;
  this->b = b;
  o       = b << 1;
  size    = (o + 3) * n;
  sol     = rhs = mat = 0;
  if (size) {
    sol = arr = new double[size];
    rhs = sol + n;
    mat = rhs + n;
    while (size--)  *arr++ = 0.0; }
}

// ----------------------------------------------------------------------------

void BandedSystem::ClearRhs() {
  long    s = n;
  double *p = rhs;
  while (s--)  *p++ = 0.0; }

// ----------------------------------------------------------------------------

void BandedSystem::ClearMat() {
  long    s = (o + 1) * n;
  double *p = mat;
  while (s--)  *p++ = 0.0; }

// ----------------------------------------------------------------------------

// Cholesky "C' * C = M" decomposition of band matrix (band is preserved):

void BandedSystem::Cholesky() {
  double *pi, *pj, *pji, *pkj, *pki, s;
  long    i, j, k, m, mo;
  for (i=0, mo=-o, pi=mat; i < n; ++i, pi += o + 1) {
    if (i <= b) { m = i;  mo += o; } else m = b;
    for (j = 0, pji=pi-m, pj=pji-mo; j <= m; ++j, ++pji, pj += o+1) {
      for (k=0, pkj=pj, pki=pji, s=*pji; k < j; ++k)  s -= *--pkj * *--pki;
      if (j == m) {
        if (s < 1e-14)  Error("System matrix not positive definite");
        *pji = sqrt(s); }
      else  *pji = s / *pj; } }
}

// ----------------------------------------------------------------------------

// Solve "C' * C * sol = rhs":

void BandedSystem::Solve() {
  double *ps, *pi, *pji, *pjs, s;
  long    i, j, m;
  for (i=0; i < n; ++i) sol[i] = rhs[i];
  // Solve for transposed cholesky factor:
  for (i=0, pi=mat, ps=sol; i < n; ++i, pi += o + 1) {
    m = i <= b ? i : b;
    for (j=0, pji=pi, pjs=ps, s=*ps; j < m; ++j)  s -= *--pji * *--pjs;
    *ps++ = s / *pi; }
  // Solve for cholesky factor:
  for (i=0, pi=mat+(n-1)*(o+1), ps=sol+(n-1); i < n; ++i, pi -= o + 1) {
    m = i <= b ? i : b;
    for (j=0, pji=pi+o, pjs=ps, s=*ps; j < m; ++j, pji += o)
      s -= *pji * *++pjs;
    *ps-- = s / *pi; }
}

// ============================================================================
// Main routine.
// ============================================================================

// Main program to supply command line args, catch errors and
// call the 'Main' routine:

int Handler(int argc, char *argv[], Parameter& param)
{
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

// ============================================================================
// Explicit template instantiation.
// ============================================================================

template class Interpolation2<double>;
template class Interpolation2<pixel_t>;
template class Interpolation2<deform_t>;
