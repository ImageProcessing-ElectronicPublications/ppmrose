# PPMROSE

```
    Author: mhr
    Posts: 37
    Joined: 07 May 2012, 10:12
    E-book readers owned: onyx-boox-m92 sony-trs-t1
    Number of books owned: 500
    Country: Germany
```

Dewarping with a regular calibration point grid

#### Post by mhr » 28 May 2012, 06:38

Hello all,

I refer to the post viewtopic.php?f=17&t=113&hilit=calibration, in which
a calibration point grid is suggested to do a calibrated dewarping if the
camera has a fixed position with respect to the scanned page. There are
other approaches like viewtopic.php?f=13&t=784 to automatically recognize
calibration information for dekeystoning. The latter seems appropriate
for setups where the relative position is changing during scanning.

I have written a small command line program, which is capable to correct
keystoning, page warping and camera lense distortion automatically by
using a calibration point grid as suggested in the first referenced post
from above. In the following posts I will discuss the approach taken
and finally will supply the program, which is GPLed.

The general setup is the following:

Your scanning setup guarantees that the camera is in fixed geometrical
position with respect to the platen. In contrast to all other setups,
the platen need not be flat! A non-flat platen seems strange at first
sight, but can be advantageous to avoid blur. I have estimated an absolute
depth of focus of about 1 cm for a camera like mine (a Powershot A490),
which seems to be roughly independent from the distance of the camera
and the platen. Therefore, with the pythagorean theorem you can show,
that noticeable blurring occurs quite fast for usual camera setups. If
dR denotes the maximal depth of focus, R the vertical camera to platen
distance and r the radius of a circle on a page with focused content,
then the formula

```
r = sqrt(dR * (2*R + dR))
```

gives the relationship between these values, providing that the camera
is orthogonal to the platen. Some examples for `dR = 1 cm`:

```
R = 44 cm ==> r = 9.4 cm
R = 60 cm ==> r = 11 cm
R = 80 cm ==> r = 12.6 cm
R = 100 cm ==> r = 14.1 cm
```

If the platen is flat, `2*r` is the maximal diameter of a page fully in
focus, but if the platen is part of a cylindrical surface, `2*r` is the
maximal height of a page fully in focus!  Therefore I plan to construct
a scanner with a cylindrical platen to enhance sharpness.  The distortion
will be no problem for me as the following posts should demonstrate.
