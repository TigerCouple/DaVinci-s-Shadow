# DaVinci-s-Shadow
comes from http://www.beneaththewaves.net/Software/DaVincis_Shadow.html

HERE I JUST SAVE THE 2.9.2 VERSION CODE

DaVinci's Shadow is an experimental fork I created of DaVinci, a command-line data-processing tool written by Arizona State University for NASA's Jet Propulsion Laboratory. The changes I made are minor, but very necessary for some of the core functionality I wanted to use in The Mirror's Surface Breaks.

Alterations From The Mainline Version

DaVinci's Shadow 2.9.2 is based on version 2.09 of DaVinci (as of version 1.3 of The Mirror's Surface Breaks - prior releases were based on version 2.05 of DaVinci), and has had these modifications performed:

libpng 1.2.3 has been replaced with libpng 1.2.46. This required a minor modification to Makefile.in for the iomedley subcomponent.

The standard DaVinci code has a long-standing bug related to saving 16-bit-per-channel (AKA "48-bit") TIFFs. The data is represented internally as 16-bit signed integers, but TIFFs use a 16-bit signed format. DaVinci shifts the data down by 32768 when TIFFs are read, but does not shift it up by 32768 when saving TIFFs, resulting in corrupted output images. DaVinci's Shadow corrects this.

The standard DaVinci code has a similar bug (also present for some time) related to 16-bit-per-channel PNG files, and in addition, this type of PNG file is not read correctly, because the 16-bit values are read in big-endian order instead of little-endian order. DaVinci's Shadow corrects this as well, although the code has not been tested on a big-endian platform to verify that it hasn't broken the same thing there.

The standard DaVinci code introduced a new bug related to TIFFs sometime after 2.05, in which 16- and 32-bit-per-channel images with three or four channels were stored in a way that caused other software to interpret them as a single greyscale channel and multiple alpha channels, instead of an RGB colour image or an ARGB colour image. DaVinci' Shadow corrects this as well.

The TIFF output code has been modified to use deflate/zip compression instead of LZW, because it results in much smaller files. [ This change was incorporated into the main-line DaVinci. One line of source code I wrote has a non-zero probability of being used for processing images from space probes! ]

Built-in functions for converting values in the RGB colourspace to and from values in the YUV, HSL, HSY, and IHSL colourspaces have been added.

rgb2yuv()

yuv2rgb()

rgb2hsv()

hsv2rgb()

rgb2hsl()

hsl2rgb()

rgb2hsy()

hsy2rgb()

rgb2ihsl()

ihsl2rgb()

rgb2ynuv() - same as rgb2yuv(), except returns U and V values in the range 0.0 - 1.0

ynuv2rgb() - same as yuv2rgb(), except accepts U and V values in the range 0.0 - 1.0

rgb2pdfhsy() - same as rgb2hsy(), except uses the slightly different conversion coefficients specified in Adobe's PDF specification.

pdfhsy2rgb() - same as hsy2rgb(), except uses the slightly different conversion coefficients specified in Adobe's PDF specification.

A built-in function for calculating the range (range()), median (med()), and an extended version of the "moment" function (momentext()) which returns the median and range in addition to the standard moments of the original function have been added.
The THEMIS-related code has been disabled in the Windows® build, because I couldn't get the module to load. It is still enabled in the Linux/Unix-like build, although I have not tested it.

Numerous custom (and customized) user-defined functions have been added in the dshadow.dvrc library file.
 
General-purpose data-manipulation functions:  

anormalize() - performs an aggressive normalization on the input data - this is discussed in detail in TMSB XML Schema Part 6: Transformation Profiles.

clipdata() - clips the upper and/or lower bound of the input data.

iiaz() - "invert independently about zero" - inverts the relationship between large and small values, while retaining the distinction between positive and negative values.

naav() - normalizes the data about an arbitrary value - this is discussed in detail in TMSB XML Schema Part 6: Transformation Profiles.

todouble() - a shortcut function to convert input data of arbitrary numeric type to double-precision floating point in the specified range (default: 0.0 - 1.0).

tofloat() - a shortcut function to convert input data of arbitrary numeric type to single-precision floating point in the specified range (default: 0.0 - 1.0).

tointeger() - a shortcut function to convert input data of arbitrary numeric type to byte (8-bit unsigned), short (16-bit signed), or int (32-bit signed) integer format.

tobyte() - a wrapper for tointeger() which always converts to byte (8-bit unsigned) format.

toint() - a wrapper for tointeger() which always converts to int (32-bit signed) format

toshort() - a wrapper for tointeger() which always converts to short (16-bit signed) format
 
General-purpose image-related functions:  

applypalette() - applies a palette to a greyscale image. Designed for use specifically with gradients generated by buildgradient().

blendcolourhsy() - performs a luminance/colour blend similar to the method used in Photoshop® (see Secrets of Photoshop's Colour Blend Mode Revealed (Sort Of))

blendcolourrgb() - a wrapper for blendcolourhsy() which accepts RGB input images instead of HSY

buildgradient() - generates a colour gradient based on user-defined control points/colours.

momentcube() - uses the one-dimensional output of DaVinci's moment() function to generate a derived three-dimensional image cube which contains statistically-calculated greyscale images (see Calculated Greyscale Images).
 
Specialized image-related functions. For more information on the DCS functionality, see Decorrelation Stretch Images. For more information on the vegetation index/satellite imagery-related functions, see Calculated Greyscale Images.  

dcsf() - a modified version of Noel Gorelick's dcs() function which uses floating-point values instead of 8-bit integers (for higher image fidelity).

dcss() - the same as dscf(), except it uses 16-bit signed integers.

rdcsf() - a modified version of C. Edwards' rdcs() function which uses floating-point values instead of 8-bit integers (for higher image fidelity).

rdcss() - the same as rdscf(), except it uses 16-bit signed integers.

arvi() - generates Atmospherically Resistant Vegetation Index images.

asvi() - generates Atmospheric and Soil Vegetation Index images.

avi() - generates Angular Vegetation Index images.

evi() - generates Enhanced Vegetation Index images. Note: I added a few optional parameters to allow the output to be tweaked if so desired.

gemi() - generates Global Environmental Monitoring Index images.

msavi() - generates Modified Soil-Adjusted Vegetation Index images.

msavi2() - generates Modified Soil-Adjusted Vegetation Index 2 images.

ndi() - a shortcut function for generating Normalized Difference images of various types (NDVI, et cetera).

slavi() - generates Specific Leaf Area Vegetation Index images.

kt1() - a shortcut function for generating KT1 ("Tasseled Cap") images.

kt2() - a shortcut function for generating KT2 ("Tasseled Cap") images

kt3() - a shortcut function for generating KT3 ("Tasseled Cap") images

Known Issues

This application uses a fairly large amount of memory. The use of the 64-bit version is highly recommended. See the TMSB Troubleshooting and Optimization for more information.

32-bit-per-channel RGB and ARGB TIFFs are written in a slightly odd format. Some software (for example, CinePaint) is perfectly happy to read them. Adobe Photoshop® will read the file, but instead of assigning the data to the red, green, and blue channels, it will treat the red channel as a greyscale channel, and everything else as a bunch of alpha channels. I have not been able to find a good fix for this, although it's easy enough to move the channels around in Photoshop® after the fact if necessary. I can't imagine many people are using images with this fidelity, so hopefully it's not a big problem.

The statistical calculations for an image cube are very slow compared to other operations. This becomes extremely noticeable with large window sizes. For example, while the calculations with a window size of 1 may take several minutes, the same calculations with a window size of 81 may take several days. While this is partly due to the large volume of data, it is also significantly exacerbated by DaVinci's relatively slow execution of for loops in user-defined functions. I do not know enough C to create a built-in function equivalent, which is the only reasonable way of addressing this limitation.

Gradient-mapping is also slow compared to other operations, for a similar reason (it involves the use of a for loop in a user-defined function).

DaVinci's Shadow does not include the optional add-ons (GNUPlot, ImageJ, etc.) that the full version does. This is because I didn't want to be obligated to provide source code for them. If you need these add-ons, install the regular version of DaVinci from ASU's site, then copy those add-ons to the DaVinci's Shadow directory.
