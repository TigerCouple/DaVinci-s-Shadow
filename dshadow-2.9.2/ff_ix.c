#include "parser.h"
#include "func.h"

/**
 ** histogram()
 **
 ** rgb2hsi()
 ** hsi2rgb()
 ** threshold()
 ** saturate()
 ** scale()
 ** stretch()
 ** color_index()
 **/

/**
 ** compute histogram of an image
 **/

/* DaVinci's Shadow customization - custom colourspace conversions - part 1 of 2 - Begin */
const double PI = 3.141592653589793;

typedef struct { float r, g, b; } RGB;
typedef struct { float y, u, v; } YUV;
typedef struct { float h, s, v; } HSV;
typedef struct { float h, s, l; } HSL;
typedef struct { float h, s, y; } HSY;
typedef struct { float c, m, y; } CMY;
typedef struct { float c, m, y, k; } KCMY;
RGB YUVToRGB(YUV yuv);
YUV RGBToYUV(RGB rgb);
RGB YNUVToRGB(YUV yuv);
YUV RGBToYNUV(RGB rgb);
RGB HSVToRGB(HSV hsv);
HSV RGBToHSV(RGB rgb);
RGB HSLToRGB(HSL hsl);
HSL RGBToHSL(RGB rgb);
RGB HSYToRGB(HSY hsy);
HSY RGBToHSY(RGB rgb);
RGB PDFHSYToRGB(HSY hsy);
HSY RGBToPDFHSY(RGB rgb);
RGB IHSLToRGB(HSY hsy);
HSY RGBToIHSL(RGB rgb);
/* DaVinci's Shadow customization - custom colourspace conversions - part 1 of 2 - End */

Var * fb_min(Var *obj, int axis, int direction);

Var *
ff_histogram(vfuncptr func, Var * arg)
{
	Var *obj = NULL, *compress = NULL, *normalize = NULL, *cumulative=NULL;
	int x,y,z, i, j, dsize;
	float *data;
	float v;

	float start = MAXFLOAT, size= MAXFLOAT;
	int steps = MAXINT;

	Alist alist[9];
	alist[0] = make_alist( "object",    ID_VAL,    NULL,    &obj);
	alist[1] = make_alist( "compress",  ID_VAL,    NULL,    &compress);
	alist[2] = make_alist( "normalize", ID_VAL,    NULL,    &normalize);
	alist[3] = make_alist( "cumulative",ID_VAL,    NULL,    &cumulative);
	alist[4] = make_alist( "start",     FLOAT,     NULL,    &start);
	alist[5] = make_alist( "size",      FLOAT,     NULL,    &size);
	alist[6] = make_alist( "steps",      INT,       NULL,    &steps);
	alist[7].name = NULL;

	if (parse_args(func, arg, alist) == 0) return(NULL);

	if (obj == NULL) {
		parse_error("%s: No object specified\n", func->name);
		return(NULL);
	}

	x = GetSamples(V_SIZE(obj), V_ORG(obj));
	y = GetLines(V_SIZE(obj), V_ORG(obj));
	z = GetBands(V_SIZE(obj), V_ORG(obj));
	dsize = V_DSIZE(obj);

	switch(V_FORMAT(obj)) {
		case BYTE: 
			if (start == MAXFLOAT) start = 0;
			if (size == MAXFLOAT) size = 1;
			if (steps == MAXINT) steps = 256;
			break;
		case SHORT: 
			if (start == MAXFLOAT) start = -32768;
			if (size == MAXFLOAT) size = 1;
			if (steps == MAXINT) steps = 65536;
			break;
		case INT:
			if (steps == MAXFLOAT) {
				parse_error("%s(...steps=...) required for INT format.", func->name);
				return(NULL);
			}
			break;
		case FLOAT:
			if (steps == MAXFLOAT) {
				parse_error("%s(...steps=...) required for FLOAT format.", func->name);
				return(NULL);
			}
			break;
		case DOUBLE:
			if (steps == MAXFLOAT) {
				parse_error("%s(...steps=...) required for DOUBLE format.", func->name);
				return(NULL);
			}
			break;
	}

	/*
	** Find minimum if necessary
	*/
	if (start == MAXFLOAT) {
		Var *vmin;
		vmin = fb_min(obj, 7, 0);
		if (vmin) start = V_FLOAT(vmin);
	}

	/*
	** find maximum if necessary
	*/
	if (size == MAXFLOAT) {
		Var *vmax;
		vmax = fb_min(obj, 7, 1);
		if (vmax) size = (V_FLOAT(vmax) - start) / steps;
	}

	if (start == MAXFLOAT || steps == MAXINT || steps == 0 || size == MAXFLOAT) {
		parse_error("Unable to determine start, steps or size");
		return(NULL);
	}

	data = (float *)calloc((steps+1) * 2, sizeof(float));

	/*
	** X axis
	*/
	for (i = 0 ; i < steps ; i++) {
		data[i*2] = start + i *size;
	}

	for (i = 0 ; i < dsize ; i++) {
		v = extract_float(obj, i);
		j = (v - start) / size;
		if (j < 0) {
			j = 0;
		} else if (j >= steps) {
			j = steps-1;
		}
		data[j*2+1]++;
	}

/*
** Exercise the options to compress, accumulate or normalize
*/

	j = steps;
	if (compress != NULL) {
		j = 0;
		for (i = j = 0 ; i < steps ; i++) {
			if (data[i*2+1] != 0) {
				data[j*2] = data[i*2];
				data[j*2+1] = data[i*2+1];
				j++;
			}
		}
	}

	if (normalize) {
		for (i = 0 ; i < j ; i++) {
			data[i*2+1] = (float)data[i*2+1]/(float)dsize;
		}
	}

	if (cumulative != NULL) {
		for (i = 1 ; i < j ; i++) {
			data[i*2+1] += data[(i-1)*2+1];
		}
	}

	return(newVal(BSQ, 2, j, 1, FLOAT, data));
}


Var *
ff_hstats(vfuncptr func, Var * arg)
{
	Var *obj = NULL;
	int x,y,z, i, j;
	Var *both, *avg, *stddev;
	double sum, sum2;
	double x1, y1;
	int n;

	Alist alist[2];
	alist[0] = make_alist( "object",    ID_VAL,    NULL,    &obj);
	alist[1].name = NULL;

	if (parse_args(func, arg, alist) == 0) return(NULL);

	if (obj == NULL) {
		parse_error("%s: No object specified\n", func->name);
		return(NULL);
	}

	x = GetSamples(V_SIZE(obj), V_ORG(obj));
	y = GetLines(V_SIZE(obj), V_ORG(obj));
	z = GetBands(V_SIZE(obj), V_ORG(obj));

	if (x != 2 || z != 1 || y < 2) {
		parse_error("Object does not look like a histogram.");
		return(NULL);
	}

	avg = newVal(BSQ, 1, 1, 1, DOUBLE, calloc(1, sizeof(double)));
	stddev = newVal(BSQ, 1, 1, 1, DOUBLE, calloc(1, sizeof(double)));

	/*
	** Use one pass method
	*/
	sum = 0;
	sum2 = 0;
	n = 0;
	for (i = 0 ; i < y ; i++) {
		j = cpos(0, i, 0, obj);

		x1 = extract_double(obj, j);
		y1 = extract_double(obj, j+1);

		sum += x1*y1;
		sum2 += (x1*x1)*y1;
		n += y1;
	}
	V_DOUBLE(stddev) = sqrt((sum2 - (sum*sum/n))/(n-1));
	V_DOUBLE(avg) = sum/n;

	both = new_struct(0);
	add_struct(both, "avg", avg);
	add_struct(both, "stddev", stddev);
	return(both);
}

Var *
ff_rgb2hsv(vfuncptr func, Var * arg)
{
	Var *obj = NULL, *maxval = NULL;
	float *data;
	double mval;
	int x,y,z,i,j,k1,k2,k3;

	RGB a;
	HSV b;

	Alist alist[3];
	alist[0] = make_alist( "object",    ID_VAL,    NULL,     &obj);
	alist[1] = make_alist( "maxval",  ID_VAL,    NULL,     &maxval);
	alist[2].name = NULL;

	if (parse_args(func, arg, alist) == 0) return(NULL);

	if (obj == NULL) {
		parse_error("%s: No object specified\n", func->name);
		return(NULL);
	}

	if (maxval == NULL) {
		switch (V_FORMAT(obj)) {
			case BYTE:		mval = (1 << 8)-1; break;
			case SHORT:		mval = MAXSHORT; break;
			case INT:		mval = MAXINT; break;
			case FLOAT:		mval = 1.0; break;
			case DOUBLE:	mval = 1.0; break;
		}
	} else {
		mval = extract_double(maxval, 0);
	}

	x = GetSamples(V_SIZE(obj), V_ORG(obj));
	y = GetLines(V_SIZE(obj), V_ORG(obj));
	z = GetBands(V_SIZE(obj), V_ORG(obj));

	if (z != 3) {
            parse_error("%s: Input must have 3 bands.\n", func->name);
	    return(NULL);
	}

	data = (float *)calloc(4, 3*x*y);
	for (i = 0 ; i < y ; i++) {
		for (j = 0 ; j < x ; j++) {
			k1 = cpos(j,i,0, obj);
			k2 = cpos(j,i,1, obj);
			k3 = cpos(j,i,2, obj);

			a.r = extract_double(obj, k1) / mval;
			a.g = extract_double(obj, k2) / mval;
			a.b = extract_double(obj, k3) / mval;

			b = RGBToHSV(a);

			data[k1] = b.h;
			data[k2] = b.s;
			data[k3] = b.v;
		}
	}
	return(newVal(V_ORG(obj), 
		V_SIZE(obj)[0], 
		V_SIZE(obj)[1], 
		V_SIZE(obj)[2], 
		FLOAT, data));
}

Var *
ff_hsv2rgb(vfuncptr func, Var * arg)
{
  Var *obj = NULL;
	float *data;
	double mval = 1.0;
	int x,y,z,i,j,k1,k2,k3;

	HSV a;
	RGB b;

	Alist alist[3];
	alist[0] = make_alist( "object",  ID_VAL,    NULL,     &obj);
	alist[1] = make_alist( "maxval",  DOUBLE,    NULL,     &mval);
	alist[2].name = NULL;

	if (parse_args(func, arg, alist) == 0) return(NULL);

	if (obj == NULL) {
		parse_error("%s: No object specified\n", func->name);
		return(NULL);
	}

	x = GetSamples(V_SIZE(obj), V_ORG(obj));
	y = GetLines(V_SIZE(obj), V_ORG(obj));
	z = GetBands(V_SIZE(obj), V_ORG(obj));

	if (z != 3) {
		parse_error("%s: Input must have 3 bands.\n", func->name);
		return(NULL);
	}

	data = (float *)calloc(4, 3*x*y);
	for (i = 0 ; i < y ; i++) {
		for (j = 0 ; j < x ; j++) {
			k1 = cpos(j,i,0, obj);
			k2 = cpos(j,i,1, obj);
			k3 = cpos(j,i,2, obj);

			a.h = extract_double(obj, k1);
			a.s = extract_double(obj, k2);
			a.v = extract_double(obj, k3);

			b = HSVToRGB(a);

			data[k1] = b.r*mval;
			data[k2] = b.g*mval;
			data[k3] = b.b*mval;
		}
	}
	return(newVal(V_ORG(obj), 
		V_SIZE(obj)[0], 
		V_SIZE(obj)[1], 
		V_SIZE(obj)[2], 
		FLOAT, data));
}


HSV
RGBToHSV(RGB rgb)
{
    HSV hsv;
    float   mn, mx;
    float   rc, gc, bc;
    
    mx = max(max(rgb.r, rgb.g), rgb.b);
    mn = min(min(rgb.r, rgb.g), rgb.b);
    hsv.v = mx;
    if (mx == 0.0)
        hsv.s = 0.0;
    else
        hsv.s = (mx - mn) / mx;
    if (hsv.s == 0.0)
        hsv.h = 0.0;
    else {
        rc = (mx - rgb.r) / (mx - mn);
        gc = (mx - rgb.g) / (mx - mn);
        bc = (mx - rgb.b) / (mx - mn);
        if (rgb.r == mx)
            hsv.h = bc - gc;
        else if (rgb.g == mx)
            hsv.h = 2.0 + rc - bc;
        else if (rgb.b == mx)
            hsv.h = 4.0 + gc - rc;
 
        if (hsv.h < 0.0)
            hsv.h += 6.0;
        hsv.h = hsv.h / 6.0;
    }
    return hsv;
}

RGB
HSVToRGB(HSV hsv)
{
    RGB rgb;
    float   p, q, t, f;
    int i;
    
    if (hsv.s == 0.0)
		rgb.r = rgb.b = rgb.g = hsv.v;
    else {
        if (hsv.s > 1.0) hsv.s = 1.0;
        if (hsv.s < 0.0) hsv.s = 0.0;
        if (hsv.v > 1.0) hsv.v = 1.0;
        if (hsv.v < 0.0) hsv.v = 0.0;
        while (hsv.h >= 1.0)
            hsv.h -= 1.0;

        hsv.h = 6.0 * hsv.h;
        i = (int) hsv.h;
        f = hsv.h - (float) i;
        p = hsv.v * (1.0 - hsv.s);
        q = hsv.v * (1.0 - (hsv.s * f));
        t = hsv.v * (1.0 - (hsv.s * (1.0 - f)));
 
        switch(i) {
			case 0: rgb.r = hsv.v; rgb.g = t; rgb.b = p; break;
			case 1: rgb.r = q; rgb.g = hsv.v; rgb.b = p; break;
			case 2: rgb.r = p; rgb.g = hsv.v; rgb.b = t; break;
			case 3: rgb.r = p; rgb.g = q; rgb.b = hsv.v; break;
			case 4: rgb.r = t; rgb.g = p; rgb.b = hsv.v; break;
			case 5: rgb.r = hsv.v; rgb.g = p; rgb.b = q; break;
		}
	}
	return rgb;
}


/* DaVinci's Shadow customization - custom colourspace conversions - part 2 of 2 - Begin */
Var *
ff_rgb2hsl(vfuncptr func, Var * arg)
{
	Var *obj = NULL, *maxval = NULL;
	float *data;
	double mval;
	int x,y,z,i,j,k1,k2,k3;

	RGB a;
	HSL b;

	Alist alist[3];
	alist[0] = make_alist( "object",    ID_VAL,    NULL,     &obj);
	alist[1] = make_alist( "maxval",  ID_VAL,    NULL,     &maxval);
	alist[2].name = NULL;

	if (parse_args(func, arg, alist) == 0) return(NULL);

	if (obj == NULL) {
		parse_error("%s: No object specified\n", func->name);
		return(NULL);
	}

	if (maxval == NULL) {
		switch (V_FORMAT(obj)) {
			case BYTE:		mval = (1 << 8)-1; break;
			case SHORT:		mval = MAXSHORT; break;
			case INT:		mval = MAXINT; break;
			case FLOAT:		mval = 1.0; break;
			case DOUBLE:		mval = 1.0; break;
		}
	} else {
		mval = extract_double(maxval, 0);
	}

	x = GetSamples(V_SIZE(obj), V_ORG(obj));
	y = GetLines(V_SIZE(obj), V_ORG(obj));
	z = GetBands(V_SIZE(obj), V_ORG(obj));

	if (z != 3) {
            parse_error("%s: Input must have 3 bands.\n", func->name);
	    return(NULL);
	}

	data = (float *)calloc(4, 3*x*y);
	for (i = 0 ; i < y ; i++) {
		for (j = 0 ; j < x ; j++) {
			k1 = cpos(j,i,0, obj);
			k2 = cpos(j,i,1, obj);
			k3 = cpos(j,i,2, obj);

			a.r = extract_double(obj, k1) / mval;
			a.g = extract_double(obj, k2) / mval;
			a.b = extract_double(obj, k3) / mval;

			b = RGBToHSL(a);

			data[k1] = b.h;
			data[k2] = b.s;
			data[k3] = b.l;
		}
	}
	return(newVal(V_ORG(obj), 
		V_SIZE(obj)[0], 
		V_SIZE(obj)[1], 
		V_SIZE(obj)[2], 
		FLOAT, data));
}

Var *
ff_hsl2rgb(vfuncptr func, Var * arg)
{
  Var *obj = NULL;
	float *data;
	double mval = 1.0;
	int x,y,z,i,j,k1,k2,k3;

	HSL a;
	RGB b;

	Alist alist[3];
	alist[0] = make_alist( "object",  ID_VAL,    NULL,     &obj);
	alist[1] = make_alist( "maxval",  DOUBLE,    NULL,     &mval);
	alist[2].name = NULL;

	if (parse_args(func, arg, alist) == 0) return(NULL);

	if (obj == NULL) {
		parse_error("%s: No object specified\n", func->name);
		return(NULL);
	}

	x = GetSamples(V_SIZE(obj), V_ORG(obj));
	y = GetLines(V_SIZE(obj), V_ORG(obj));
	z = GetBands(V_SIZE(obj), V_ORG(obj));

	if (z != 3) {
		parse_error("%s: Input must have 3 bands.\n", func->name);
		return(NULL);
	}

	data = (float *)calloc(4, 3*x*y);
	for (i = 0 ; i < y ; i++) {
		for (j = 0 ; j < x ; j++) {
			k1 = cpos(j,i,0, obj);
			k2 = cpos(j,i,1, obj);
			k3 = cpos(j,i,2, obj);

			a.h = extract_double(obj, k1);
			a.s = extract_double(obj, k2);
			a.l = extract_double(obj, k3);

			b = HSLToRGB(a);

			data[k1] = b.r*mval;
			data[k2] = b.g*mval;
			data[k3] = b.b*mval;
		}
	}
	return(newVal(V_ORG(obj), 
		V_SIZE(obj)[0], 
		V_SIZE(obj)[1], 
		V_SIZE(obj)[2], 
		FLOAT, data));
}

/* based on http://www.mjijackson.com/2008/02/rgb-to-hsl-and-rgb-to-hsv-color-model-conversion-algorithms-in-javascript */
HSL
RGBToHSL(RGB rgb)
{
	HSL hsl;
	float   mn, mx, av;
	float	delta, rhv;
	float   rc, gc, bc;
    
	mx = max(max(rgb.r, rgb.g), rgb.b);
	mn = min(min(rgb.r, rgb.g), rgb.b);
	av = (rgb.r + rgb.g + rgb.b) / 3.0;

	hsl.h = av;
	hsl.s = av;
	hsl.l = av;

	if (mx == mn)
	{
		hsl.h = 0.0;
		hsl.s = 0.0;
	}
	else
	{
		delta = mx - mn;
		if (hsl.l > 0.5)
		{
			hsl.s = delta / (2.0 - mx - mn);
		}
		else
		{
			hsl.s = delta / (mx + mn);
		}
		if (rgb.r == mx)
		{
			if (rgb.g < rgb.b)
			{
				rhv = 6.0;
			}
			else
			{
				rhv = 0.0;
			}
			hsl.h = (rgb.g - rgb.b) / delta + rhv;
		}
		else if (rgb.g == mx)
		{
			hsl.h = (rgb.b - rgb.r) / delta + 2.0;
		}
		else if (rgb.b == mx)
		{
			hsl.h = (rgb.r - rgb.g) / delta + 4.0;
		}
		hsl.h = hsl.h / 6.0;
	}

	return hsl;
}

/* based on http://www.mjijackson.com/2008/02/rgb-to-hsl-and-rgb-to-hsv-color-model-conversion-algorithms-in-javascript */
float
HueToChannel(float p, float q, float t)
{
	if (t < 0)
	{
		t += 1;
	}
	if (t > 1)
	{
		t -= 1;
	}
	if (t < (1.0/6.0))
	{
		return (float)(p + (q - p) * 6.0 * t);
	}
	if (t < 0.5)
	{
		return q;
	}
	if (t < (2.0/3.0))
	{
		return (float)(p + (q - p) * ((2.0/3.0) - t) * 6.0);
	}
	return p;
}

/* based on http://www.mjijackson.com/2008/02/rgb-to-hsl-and-rgb-to-hsv-color-model-conversion-algorithms-in-javascript */
RGB
HSLToRGB(HSL hsl)
{
	RGB rgb;
	float   p, q, t, f;
	int i;
    
	if (hsl.s == 0.0)
	{
		rgb.r = hsl.l;
		rgb.g = hsl.l;
		rgb.b = hsl.l;
	}
	else
	{
		if (hsl.l < 0.5)
		{
			q = hsl.l * (1.0 + hsl.s);
		}
		else
		{
			q = hsl.l + hsl.s - (hsl.l * hsl.s);
		}
		p = 2 * hsl.l - q;
		rgb.r = HueToChannel(p, q, hsl.h + (1.0/3.0));
		rgb.g = HueToChannel(p, q, hsl.h);
		rgb.b = HueToChannel(p, q, hsl.h - (1.0/3.0));
	}

	return rgb;
}


Var *
ff_rgb2hsy(vfuncptr func, Var * arg)
{
	Var *obj = NULL, *maxval = NULL;
	float *data;
	double mval;
	int x,y,z,i,j,k1,k2,k3;

	RGB a;
	HSY b;

	Alist alist[3];
	alist[0] = make_alist( "object",    ID_VAL,    NULL,     &obj);
	alist[1] = make_alist( "maxval",  ID_VAL,    NULL,     &maxval);
	alist[2].name = NULL;

	if (parse_args(func, arg, alist) == 0) return(NULL);

	if (obj == NULL) {
		parse_error("%s: No object specified\n", func->name);
		return(NULL);
	}

	if (maxval == NULL) {
		switch (V_FORMAT(obj)) {
			case BYTE:		mval = (1 << 8)-1; break;
			case SHORT:		mval = MAXSHORT; break;
			case INT:		mval = MAXINT; break;
			case FLOAT:		mval = 1.0; break;
			case DOUBLE:		mval = 1.0; break;
		}
	} else {
		mval = extract_double(maxval, 0);
	}

	x = GetSamples(V_SIZE(obj), V_ORG(obj));
	y = GetLines(V_SIZE(obj), V_ORG(obj));
	z = GetBands(V_SIZE(obj), V_ORG(obj));

	if (z != 3) {
            parse_error("%s: Input must have 3 bands.\n", func->name);
	    return(NULL);
	}

	data = (float *)calloc(4, 3*x*y);
	for (i = 0 ; i < y ; i++) {
		for (j = 0 ; j < x ; j++) {
			k1 = cpos(j,i,0, obj);
			k2 = cpos(j,i,1, obj);
			k3 = cpos(j,i,2, obj);

			a.r = extract_double(obj, k1) / mval;
			a.g = extract_double(obj, k2) / mval;
			a.b = extract_double(obj, k3) / mval;

			b = RGBToHSY(a);

			data[k1] = b.h;
			data[k2] = b.s;
			data[k3] = b.y;
		}
	}
	return(newVal(V_ORG(obj), 
		V_SIZE(obj)[0], 
		V_SIZE(obj)[1], 
		V_SIZE(obj)[2], 
		FLOAT, data));
}

Var *
ff_hsy2rgb(vfuncptr func, Var * arg)
{
  Var *obj = NULL;
	float *data;
	double mval = 1.0;
	int x,y,z,i,j,k1,k2,k3;

	HSY a;
	RGB b;

	Alist alist[3];
	alist[0] = make_alist( "object",  ID_VAL,    NULL,     &obj);
	alist[1] = make_alist( "maxval",  DOUBLE,    NULL,     &mval);
	alist[2].name = NULL;

	if (parse_args(func, arg, alist) == 0) return(NULL);

	if (obj == NULL) {
		parse_error("%s: No object specified\n", func->name);
		return(NULL);
	}

	x = GetSamples(V_SIZE(obj), V_ORG(obj));
	y = GetLines(V_SIZE(obj), V_ORG(obj));
	z = GetBands(V_SIZE(obj), V_ORG(obj));

	if (z != 3) {
		parse_error("%s: Input must have 3 bands.\n", func->name);
		return(NULL);
	}

	data = (float *)calloc(4, 3*x*y);
	for (i = 0 ; i < y ; i++) {
		for (j = 0 ; j < x ; j++) {
			k1 = cpos(j,i,0, obj);
			k2 = cpos(j,i,1, obj);
			k3 = cpos(j,i,2, obj);

			a.h = extract_double(obj, k1);
			a.s = extract_double(obj, k2);
			a.y = extract_double(obj, k3);

			b = HSYToRGB(a);

			data[k1] = b.r*mval;
			data[k2] = b.g*mval;
			data[k3] = b.b*mval;
		}
	}
	return(newVal(V_ORG(obj), 
		V_SIZE(obj)[0], 
		V_SIZE(obj)[1], 
		V_SIZE(obj)[2], 
		FLOAT, data));
}

Var *
ff_rgb2pdfhsy(vfuncptr func, Var * arg)
{
	Var *obj = NULL, *maxval = NULL;
	float *data;
	double mval;
	int x,y,z,i,j,k1,k2,k3;

	RGB a;
	HSY b;

	Alist alist[3];
	alist[0] = make_alist( "object",    ID_VAL,    NULL,     &obj);
	alist[1] = make_alist( "maxval",  ID_VAL,    NULL,     &maxval);
	alist[2].name = NULL;

	if (parse_args(func, arg, alist) == 0) return(NULL);

	if (obj == NULL) {
		parse_error("%s: No object specified\n", func->name);
		return(NULL);
	}

	if (maxval == NULL) {
		switch (V_FORMAT(obj)) {
			case BYTE:		mval = (1 << 8)-1; break;
			case SHORT:		mval = MAXSHORT; break;
			case INT:		mval = MAXINT; break;
			case FLOAT:		mval = 1.0; break;
			case DOUBLE:		mval = 1.0; break;
		}
	} else {
		mval = extract_double(maxval, 0);
	}

	x = GetSamples(V_SIZE(obj), V_ORG(obj));
	y = GetLines(V_SIZE(obj), V_ORG(obj));
	z = GetBands(V_SIZE(obj), V_ORG(obj));

	if (z != 3) {
            parse_error("%s: Input must have 3 bands.\n", func->name);
	    return(NULL);
	}

	data = (float *)calloc(4, 3*x*y);
	for (i = 0 ; i < y ; i++) {
		for (j = 0 ; j < x ; j++) {
			k1 = cpos(j,i,0, obj);
			k2 = cpos(j,i,1, obj);
			k3 = cpos(j,i,2, obj);

			a.r = extract_double(obj, k1) / mval;
			a.g = extract_double(obj, k2) / mval;
			a.b = extract_double(obj, k3) / mval;

			b = RGBToPDFHSY(a);

			data[k1] = b.h;
			data[k2] = b.s;
			data[k3] = b.y;
		}
	}
	return(newVal(V_ORG(obj), 
		V_SIZE(obj)[0], 
		V_SIZE(obj)[1], 
		V_SIZE(obj)[2], 
		FLOAT, data));
}

Var *
ff_pdfhsy2rgb(vfuncptr func, Var * arg)
{
  Var *obj = NULL;
	float *data;
	double mval = 1.0;
	int x,y,z,i,j,k1,k2,k3;

	HSY a;
	RGB b;

	Alist alist[3];
	alist[0] = make_alist( "object",  ID_VAL,    NULL,     &obj);
	alist[1] = make_alist( "maxval",  DOUBLE,    NULL,     &mval);
	alist[2].name = NULL;

	if (parse_args(func, arg, alist) == 0) return(NULL);

	if (obj == NULL) {
		parse_error("%s: No object specified\n", func->name);
		return(NULL);
	}

	x = GetSamples(V_SIZE(obj), V_ORG(obj));
	y = GetLines(V_SIZE(obj), V_ORG(obj));
	z = GetBands(V_SIZE(obj), V_ORG(obj));

	if (z != 3) {
		parse_error("%s: Input must have 3 bands.\n", func->name);
		return(NULL);
	}

	data = (float *)calloc(4, 3*x*y);
	for (i = 0 ; i < y ; i++) {
		for (j = 0 ; j < x ; j++) {
			k1 = cpos(j,i,0, obj);
			k2 = cpos(j,i,1, obj);
			k3 = cpos(j,i,2, obj);

			a.h = extract_double(obj, k1);
			a.s = extract_double(obj, k2);
			a.y = extract_double(obj, k3);

			b = PDFHSYToRGB(a);

			data[k1] = b.r*mval;
			data[k2] = b.g*mval;
			data[k3] = b.b*mval;
		}
	}
	return(newVal(V_ORG(obj), 
		V_SIZE(obj)[0], 
		V_SIZE(obj)[1], 
		V_SIZE(obj)[2], 
		FLOAT, data));
}

Var *
ff_rgb2ihsl(vfuncptr func, Var * arg)
{
	Var *obj = NULL, *maxval = NULL;
	float *data;
	double mval;
	int x,y,z,i,j,k1,k2,k3;

	RGB a;
	HSY b;

	Alist alist[3];
	alist[0] = make_alist( "object",    ID_VAL,    NULL,     &obj);
	alist[1] = make_alist( "maxval",  ID_VAL,    NULL,     &maxval);
	alist[2].name = NULL;

	if (parse_args(func, arg, alist) == 0) return(NULL);

	if (obj == NULL) {
		parse_error("%s: No object specified\n", func->name);
		return(NULL);
	}

	if (maxval == NULL) {
		switch (V_FORMAT(obj)) {
			case BYTE:		mval = (1 << 8)-1; break;
			case SHORT:		mval = MAXSHORT; break;
			case INT:		mval = MAXINT; break;
			case FLOAT:		mval = 1.0; break;
			case DOUBLE:		mval = 1.0; break;
		}
	} else {
		mval = extract_double(maxval, 0);
	}

	x = GetSamples(V_SIZE(obj), V_ORG(obj));
	y = GetLines(V_SIZE(obj), V_ORG(obj));
	z = GetBands(V_SIZE(obj), V_ORG(obj));

	if (z != 3) {
            parse_error("%s: Input must have 3 bands.\n", func->name);
	    return(NULL);
	}

	data = (float *)calloc(4, 3*x*y);
	for (i = 0 ; i < y ; i++) {
		for (j = 0 ; j < x ; j++) {
			k1 = cpos(j,i,0, obj);
			k2 = cpos(j,i,1, obj);
			k3 = cpos(j,i,2, obj);

			a.r = extract_double(obj, k1) / mval;
			a.g = extract_double(obj, k2) / mval;
			a.b = extract_double(obj, k3) / mval;

			b = RGBToIHSL(a);

			data[k1] = b.h;
			data[k2] = b.s;
			data[k3] = b.y;
		}
	}
	return(newVal(V_ORG(obj), 
		V_SIZE(obj)[0], 
		V_SIZE(obj)[1], 
		V_SIZE(obj)[2], 
		FLOAT, data));
}

Var *
ff_ihsl2rgb(vfuncptr func, Var * arg)
{
  Var *obj = NULL;
	float *data;
	double mval = 1.0;
	int x,y,z,i,j,k1,k2,k3;

	HSY a;
	RGB b;

	Alist alist[3];
	alist[0] = make_alist( "object",  ID_VAL,    NULL,     &obj);
	alist[1] = make_alist( "maxval",  DOUBLE,    NULL,     &mval);
	alist[2].name = NULL;

	if (parse_args(func, arg, alist) == 0) return(NULL);

	if (obj == NULL) {
		parse_error("%s: No object specified\n", func->name);
		return(NULL);
	}

	x = GetSamples(V_SIZE(obj), V_ORG(obj));
	y = GetLines(V_SIZE(obj), V_ORG(obj));
	z = GetBands(V_SIZE(obj), V_ORG(obj));

	if (z != 3) {
		parse_error("%s: Input must have 3 bands.\n", func->name);
		return(NULL);
	}

	data = (float *)calloc(4, 3*x*y);
	for (i = 0 ; i < y ; i++) {
		for (j = 0 ; j < x ; j++) {
			k1 = cpos(j,i,0, obj);
			k2 = cpos(j,i,1, obj);
			k3 = cpos(j,i,2, obj);

			a.h = extract_double(obj, k1);
			a.s = extract_double(obj, k2);
			a.y = extract_double(obj, k3);

			b = IHSLToRGB(a);

			data[k1] = b.r*mval;
			data[k2] = b.g*mval;
			data[k3] = b.b*mval;
		}
	}
	return(newVal(V_ORG(obj), 
		V_SIZE(obj)[0], 
		V_SIZE(obj)[1], 
		V_SIZE(obj)[2], 
		FLOAT, data));
}

/* RGBToHSY and HSYToRGB are based on the equations in "A 3D-polar Coordinate	*/
/* Colour Representation Suitable for Image Analysis" (2002), by Allan Hanbury	*/
/* and Jean Serra.								*/

HSY
_RGBToHSY(RGB rgb, float redco, float greenco, float blueco)
{
	HSY hsy;
	float   mn, mx;
	float	hueprime, hpnum, hpdenom;
    
	mx = max(max(rgb.r, rgb.g), rgb.b);
	mn = min(min(rgb.r, rgb.g), rgb.b);
	hsy.y = (redco * rgb.r) + (greenco * rgb.g) + (blueco * rgb.b);

	hsy.s = mx - mn;

	/* build the components of the hue prime equation (which is a bit long for one line)	*/
	hpnum = rgb.r - (0.5 * rgb.g) - (0.5 * rgb.b);
	hpdenom = pow((pow(rgb.r, 2.0) + pow(rgb.g, 2.0) + pow(rgb.b, 2.0) - (rgb.r * rgb.g) - (rgb.r * rgb.b) - (rgb.b * rgb.g)), 0.5);
	
	/* calculate hue prime and convert from radians to a value between 0 and 1		*/
	hueprime = (acos(hpnum / hpdenom) / (2.0 * PI));

	if (rgb.b > rgb.g)
	{
		hsy.h = 1 - hueprime;
	}
	else
	{
		hsy.h = hueprime;
	}
	/* if the hue is undefined, set it to 1.0 and saturation to zero as that seems to be	*/
	/* the least-offensive option.								*/
	if (hsy.h != hsy.h)
	{
		hsy.h = 1.0;
		hsy.s = 0.0;
		/*hsy.h = 0.852;
		hsy.s = 1.0;*/
	}

	/* clip any invalid data */
	hsy.h = min(1.0, hsy.h);
	hsy.h = max(0.0, hsy.h);
	hsy.s = min(1.0, hsy.s);
	hsy.s = max(0.0, hsy.s);
	hsy.y = min(1.0, hsy.y);
	hsy.y = max(0.0, hsy.y);

	return hsy;
}

HSY
RGBToHSY(RGB rgb)
{
	/* these are the standard weightings for calculating luminance		*/
	/* Adobe's PDF doc specifies them as rounded to 0.3, 0.59, and 0.11	*/
	/* Hanburry and Serra's paper specifies them as 0.2126, 0.7152, and	*/
	/* 0.0722.								*/
	/* Note that the same values must be used to calculate the inverted	*/
	/* matrix for the conversion back to RGB. I don't have time to figure	*/
	/* out how to invert a matrix in C++ code, so a similarly-flexible	*/
	/* HSYtoRGB function is left as an exercise to the reader.		*/
	return _RGBToHSY(rgb, 0.299, 0.587, 0.114);
}

HSY
RGBToPDFHSY(RGB rgb)
{
	/* these are the standard weightings for calculating luminance		*/
	/* Adobe's PDF doc specifies them as rounded to 0.3, 0.59, and 0.11	*/
	/* Hanburry and Serra's paper specifies them as 0.2126, 0.7152, and	*/
	/* 0.0722.								*/
	/* Note that the same values must be used to calculate the inverted	*/
	/* matrix for the conversion back to RGB. I don't have time to figure	*/
	/* out how to invert a matrix in C++ code, so a similarly-flexible	*/
	/* HSYtoRGB function is left as an exercise to the reader.		*/
	return _RGBToHSY(rgb, 0.3, 0.59, 0.11);
}

HSY
RGBToIHSL(RGB rgb)
{
	/* these are the standard weightings for calculating luminance		*/
	/* Adobe's PDF doc specifies them as rounded to 0.3, 0.59, and 0.11	*/
	/* Hanburry and Serra's paper specifies them as 0.2126, 0.7152, and	*/
	/* 0.0722.								*/
	/* Note that the same values must be used to calculate the inverted	*/
	/* matrix for the conversion back to RGB. I don't have time to figure	*/
	/* out how to invert a matrix in C++ code, so a similarly-flexible	*/
	/* HSYtoRGB function is left as an exercise to the reader.		*/
	return _RGBToHSY(rgb, 0.2126, 0.7152, 0.0722);
}

/* this is split out into a separate function so that future developers can implement the	*/
/* alternate RGB -> HSY algorithm (which also uses it) if so desired.				*/
float calcHStar(float hrad)
{
	int k;
	float hstar, sixtyDegrees;
	sixtyDegrees = (PI / 3.0);

	for (k = 0; k < 6; k++)
	{
		hstar = hrad - (k * sixtyDegrees);
		if ((hstar >= 0.0) && (hstar <= sixtyDegrees))
		{
			return hstar;
		}
	}

	return hstar;
}

/* this function is used instead of doing the calculation entirely in radians because	*/
/* the algorithm is *very* sensitive to rounding errors and will produce spurious 	*/
/* results (inverted hue, etc.) if this subsection is performed using radians.		*/
float calcHStarDegrees(float hrad)
{
	int k;
	float hstar, hDegrees;
	hDegrees = (hrad / (2.0 * PI)) * 360.0;

	while (hDegrees < 0.0)
	{
		hDegrees += 360.0;
	}
	while (hDegrees > 360.0)
	{
		hDegrees -= 360.0;
	}

	hstar = hDegrees;

	for (k = 0; k < 6; k++)
	{
		/* hstar = floor(hDegrees - (k * hDegrees)); */
		hstar = hDegrees - (k * hDegrees);
		if ((hstar >= 0.0) && (hstar <= 60.0))
		{
			break;
		}
	}

	/* convert back into radians */
	hstar = (hstar / 360) * (2.0 * PI);

	return hstar;
}

RGB
_HSYToRGB(HSY hsy, float ma1, float mb1, float mc1, float ma2, float mb2, float mc2, float ma3, float mb3, float mc3)
{
	RGB rgb;
	float hrad;
	float chroma, c1, c2;

	/* this isn't just an optimization. If the first check is not present, then areas of 0	*/
	/* saturation and 1 luminosity are black instead of white.				*/
	/* If the second check is not performed, some areas that should be grey are super-	*/
	/* saturated instead (due to a rounding error?).					*/
	/*if ((hsy.s == 0) || (hsy.s < 0.03))*/
	if (hsy.s == 0)
	{
		rgb.r = hsy.y;
		rgb.g = hsy.y;
		rgb.b = hsy.y;

		//for debugging conversion problems:
		//rgb.r = 1.0;
		//rgb.g = 0.0;
		//rgb.b = 0.0;
	}
	else
	{
		/* if hue is undefined, c1 and c2 are set to zero 	*/
		/* (per the paper).					*/
		if (hsy.h != hsy.h)
		{
			c1 = 0;
			c2 = 0;

			/* for debugging conversion problems: */
			/* c1 = 0.9;
			c2 = -0.9; */
		}
		else
		{
			/* convert h back to radians */
			hrad = hsy.h * (2.0 * PI);
		
			/* chroma = (sqrt(3.0) * hsy.s) / (2.0 * sin(((2.0 * PI) / 3.0) - calcHStarDegrees(hrad))); */
			chroma = (sqrt(3.0) * hsy.s) / (2.0 * sin(((2.0 * PI) / 3.0) - calcHStar(hrad)));
			c1 = chroma * cos(hrad);
			c2 = -chroma * sin(hrad);
		}
	
		rgb.r = (ma1 * hsy.y) + (mb1 * c1) + (mc1 * c2);
		rgb.g = (ma2 * hsy.y) + (mb2 * c1) + (mc2 * c2);
		rgb.b = (ma3 * hsy.y) + (mb3 * c1) + (mc3 * c2);
	}

	/* clip any invalid data */
	rgb.r = min(1.0, rgb.r);
	rgb.r = max(0.0, rgb.r);
	rgb.g = min(1.0, rgb.g);
	rgb.g = max(0.0, rgb.g);
	rgb.b = min(1.0, rgb.b);
	rgb.b = max(0.0, rgb.b);

	return rgb;
}

/*
float radiansToDegrees(float radians, bool mod360)
{
	float degrees;

	degrees = (radians / (2.0 * PI)) * 360.0;

	if (mod360)
	{
		while (degrees < 0.0)
		{
			degrees += 360.0;
		}
		while (degrees > 360.0)
		{
			degrees -= 360.0;
		}
	}

	return degrees;
}

float radiansToDegrees(float radians)
{
	return radiansToDegrees(radians, true);
}

float degreesToRadians(float degrees, bool mod1)
{
	float dg, radians;

	dg = degrees;

	if (mod1)
	{
		while (dg < 0.0)
		{
			dg += 360.0;
		}
		while (dg > 360.0)
		{
			dg -= 360.0;
		}
	}

	radians = (dg / 360.0) * (2.0 * PI);

	return radians;
}

float degreesToRadians(float degrees)
{
	return degreesToRadians(degrees, true);
}

*/

/* The following three functions use hardcoded values which are calculated by matrix	*/
/* inversion of the luminance channel weighting values in the corresponding RGB ->	*/
/* HSY functions, above.								*/

RGB
HSYToRGB(HSY hsy)
{
	/*	1	 0.701	 0.2731		*/
	/*	1	-0.2990	-0.3043		*/
	/*	1	-0.2990	 0.8504		*/

	return _HSYToRGB(hsy, 1.0, 0.701, 0.2731, 1.0, -0.2990, -0.3043, 1.0, -0.2990, 0.8504);
}

RGB
PDFHSYToRGB(HSY hsy)
{

	/*	1	 0.7	 0.2771		*/
	/*	1	-0.3	-0.3		*/
	/*	1	-0.3	 0.8545		*/

	return _HSYToRGB(hsy, 1, 0.7, 0.2771, 1, -0.3, -0.3, 1, -0.3, 0.8545);
}

RGB
IHSLToRGB(HSY hsy)
{

	/*	1	 0.7875	 0.3714		*/
	/*	1	-0.2125	-0.2059		*/
	/*	1	-0.2125	 0.9488		*/

	return _HSYToRGB(hsy, 1, 0.7875, 0.3714, 1, -0.2125, -0.2059, 1, -0.2125, 0.9488);
}





Var *
ff_rgb2yuv(vfuncptr func, Var * arg)
{
	Var *obj = NULL, *maxval = NULL;
	float *data;
	double mval;
	int x,y,z,i,j,k1,k2,k3;

	RGB a;
	YUV b;

	Alist alist[3];
	alist[0] = make_alist( "object",    ID_VAL,    NULL,     &obj);
	alist[1] = make_alist( "maxval",  ID_VAL,    NULL,     &maxval);
	alist[2].name = NULL;

	if (parse_args(func, arg, alist) == 0) return(NULL);

	if (obj == NULL) {
		parse_error("%s: No object specified\n", func->name);
		return(NULL);
	}

	if (maxval == NULL) {
		switch (V_FORMAT(obj)) {
			case BYTE:		mval = (1 << 8)-1; break;
			case SHORT:		mval = MAXSHORT; break;
			case INT:		mval = MAXINT; break;
			case FLOAT:		mval = 1.0; break;
			case DOUBLE:		mval = 1.0; break;
		}
	} else {
		mval = extract_double(maxval, 0);
	}

	x = GetSamples(V_SIZE(obj), V_ORG(obj));
	y = GetLines(V_SIZE(obj), V_ORG(obj));
	z = GetBands(V_SIZE(obj), V_ORG(obj));

	if (z != 3) {
            parse_error("%s: Input must have 3 bands.\n", func->name);
	    return(NULL);
	}

	data = (float *)calloc(4, 3*x*y);
	for (i = 0 ; i < y ; i++) {
		for (j = 0 ; j < x ; j++) {
			k1 = cpos(j,i,0, obj);
			k2 = cpos(j,i,1, obj);
			k3 = cpos(j,i,2, obj);

			a.r = extract_double(obj, k1) / mval;
			a.g = extract_double(obj, k2) / mval;
			a.b = extract_double(obj, k3) / mval;

			b = RGBToYUV(a);

			data[k1] = b.y;
			data[k2] = b.u;
			data[k3] = b.v;
		}
	}
	return(newVal(V_ORG(obj), 
		V_SIZE(obj)[0], 
		V_SIZE(obj)[1], 
		V_SIZE(obj)[2], 
		FLOAT, data));
}

Var *
ff_rgb2ynuv(vfuncptr func, Var * arg)
{
	Var *obj = NULL, *maxval = NULL;
	float *data;
	double mval;
	int x,y,z,i,j,k1,k2,k3;

	RGB a;
	YUV b;

	Alist alist[3];
	alist[0] = make_alist( "object",    ID_VAL,    NULL,     &obj);
	alist[1] = make_alist( "maxval",  ID_VAL,    NULL,     &maxval);
	alist[2].name = NULL;

	if (parse_args(func, arg, alist) == 0) return(NULL);

	if (obj == NULL) {
		parse_error("%s: No object specified\n", func->name);
		return(NULL);
	}

	if (maxval == NULL) {
		switch (V_FORMAT(obj)) {
			case BYTE:		mval = (1 << 8)-1; break;
			case SHORT:		mval = MAXSHORT; break;
			case INT:		mval = MAXINT; break;
			case FLOAT:		mval = 1.0; break;
			case DOUBLE:		mval = 1.0; break;
		}
	} else {
		mval = extract_double(maxval, 0);
	}

	x = GetSamples(V_SIZE(obj), V_ORG(obj));
	y = GetLines(V_SIZE(obj), V_ORG(obj));
	z = GetBands(V_SIZE(obj), V_ORG(obj));

	if (z != 3) {
            parse_error("%s: Input must have 3 bands.\n", func->name);
	    return(NULL);
	}

	data = (float *)calloc(4, 3*x*y);
	for (i = 0 ; i < y ; i++) {
		for (j = 0 ; j < x ; j++) {
			k1 = cpos(j,i,0, obj);
			k2 = cpos(j,i,1, obj);
			k3 = cpos(j,i,2, obj);

			a.r = extract_double(obj, k1) / mval;
			a.g = extract_double(obj, k2) / mval;
			a.b = extract_double(obj, k3) / mval;

			b = RGBToYNUV(a);

			data[k1] = b.y;
			data[k2] = b.u;
			data[k3] = b.v;
		}
	}
	return(newVal(V_ORG(obj), 
		V_SIZE(obj)[0], 
		V_SIZE(obj)[1], 
		V_SIZE(obj)[2], 
		FLOAT, data));
}

Var *
ff_yuv2rgb(vfuncptr func, Var * arg)
{
  Var *obj = NULL;
	float *data;
	double mval = 1.0;
	int x,y,z,i,j,k1,k2,k3;

	YUV a;
	RGB b;

	Alist alist[3];
	alist[0] = make_alist( "object",  ID_VAL,    NULL,     &obj);
	alist[1] = make_alist( "maxval",  DOUBLE,    NULL,     &mval);
	alist[2].name = NULL;

	if (parse_args(func, arg, alist) == 0) return(NULL);

	if (obj == NULL) {
		parse_error("%s: No object specified\n", func->name);
		return(NULL);
	}

	x = GetSamples(V_SIZE(obj), V_ORG(obj));
	y = GetLines(V_SIZE(obj), V_ORG(obj));
	z = GetBands(V_SIZE(obj), V_ORG(obj));

	if (z != 3) {
		parse_error("%s: Input must have 3 bands.\n", func->name);
		return(NULL);
	}

	data = (float *)calloc(4, 3*x*y);
	for (i = 0 ; i < y ; i++) {
		for (j = 0 ; j < x ; j++) {
			k1 = cpos(j,i,0, obj);
			k2 = cpos(j,i,1, obj);
			k3 = cpos(j,i,2, obj);

			a.y = extract_double(obj, k1);
			a.u = extract_double(obj, k2);
			a.v = extract_double(obj, k3);

			b = YUVToRGB(a);

			data[k1] = b.r*mval;
			data[k2] = b.g*mval;
			data[k3] = b.b*mval;
		}
	}
	return(newVal(V_ORG(obj), 
		V_SIZE(obj)[0], 
		V_SIZE(obj)[1], 
		V_SIZE(obj)[2], 
		FLOAT, data));
}

Var *
ff_ynuv2rgb(vfuncptr func, Var * arg)
{
  Var *obj = NULL;
	float *data;
	double mval = 1.0;
	int x,y,z,i,j,k1,k2,k3;

	YUV a;
	RGB b;

	Alist alist[3];
	alist[0] = make_alist( "object",  ID_VAL,    NULL,     &obj);
	alist[1] = make_alist( "maxval",  DOUBLE,    NULL,     &mval);
	alist[2].name = NULL;

	if (parse_args(func, arg, alist) == 0) return(NULL);

	if (obj == NULL) {
		parse_error("%s: No object specified\n", func->name);
		return(NULL);
	}

	x = GetSamples(V_SIZE(obj), V_ORG(obj));
	y = GetLines(V_SIZE(obj), V_ORG(obj));
	z = GetBands(V_SIZE(obj), V_ORG(obj));

	if (z != 3) {
		parse_error("%s: Input must have 3 bands.\n", func->name);
		return(NULL);
	}

	data = (float *)calloc(4, 3*x*y);
	for (i = 0 ; i < y ; i++) {
		for (j = 0 ; j < x ; j++) {
			k1 = cpos(j,i,0, obj);
			k2 = cpos(j,i,1, obj);
			k3 = cpos(j,i,2, obj);

			a.y = extract_double(obj, k1);
			a.u = extract_double(obj, k2);
			a.v = extract_double(obj, k3);

			b = YNUVToRGB(a);

			data[k1] = b.r*mval;
			data[k2] = b.g*mval;
			data[k3] = b.b*mval;
		}
	}
	return(newVal(V_ORG(obj), 
		V_SIZE(obj)[0], 
		V_SIZE(obj)[1], 
		V_SIZE(obj)[2], 
		FLOAT, data));
}


YUV
_RGBToYUV(RGB rgb, float weightRed, float weightBlue, int normalizeUV)
{
	YUV yuv;
	float weightGreen;
	float uMax;
	float vMax;

	weightGreen = 1.0 - weightRed - weightBlue;
	/* Not sure if these two are calculated or magic */
	uMax = 0.436;
	vMax = 0.615;

	yuv.y = (rgb.r * weightRed) + (rgb.g * weightGreen) + (rgb.b * weightBlue);
	yuv.u = uMax * ((rgb.b - yuv.y) / (1.0 - weightBlue));
	yuv.v = vMax * ((rgb.r - yuv.y) / (1.0 - weightRed));

	if (normalizeUV == 1)
	{
		yuv.u = yuv.u + uMax;
		yuv.u = yuv.u / (uMax * 2);
		yuv.v = yuv.v + vMax;
		yuv.v = yuv.v / (vMax * 2);
	}

	/* clip any invalid data */
	yuv.y = min(1.0, yuv.y);
	yuv.y = max(0.0, yuv.y);
	yuv.u = min(1.0, yuv.u);
	yuv.u = max(0.0, yuv.u);
	yuv.v = min(1.0, yuv.v);
	yuv.v = max(0.0, yuv.v);

	return yuv;
}

RGB
_YUVToRGB(YUV yuv, float weightRed, float weightBlue, int normalizedUV)
{
	RGB rgb;
	float weightGreen;
	float uMax;
	float vMax;
	float y;
	float u;
	float v;

	weightGreen = 1.0 - weightRed - weightBlue;
	/* Not sure if these two are calculated or magic */
	uMax = 0.436;
	vMax = 0.615;

	/* to avoid making changes to the source data */
	y = yuv.y;
	u = yuv.u;
	v = yuv.v;

	if (normalizedUV == 1)
	{
		u = u * (uMax * 2);
		u = u - uMax;
		v = v * (vMax * 2);
		v = v - vMax;
	}

	rgb.r = y + (v * ((1 - weightRed) / vMax));
	rgb.g = y - (u * ((weightBlue * (1 - weightBlue)) / (uMax * weightGreen))) - (v * ((weightRed * (1 - weightRed)) / (vMax * weightGreen)));
	rgb.b = y + (u * ((1 - weightBlue) / uMax));

	/* clip any invalid data */
	rgb.r = min(1.0, rgb.r);
	rgb.r = max(0.0, rgb.r);
	rgb.g = min(1.0, rgb.g);
	rgb.g = max(0.0, rgb.g);
	rgb.b = min(1.0, rgb.b);
	rgb.b = max(0.0, rgb.b);

	return rgb;
}

YUV
RGBToYUV(RGB rgb)
{
	/* these are the standard weightings for calculating luminance		*/
	/* Adobe's PDF doc specifies them as rounded to 0.3, 0.59, and 0.11	*/
	/* Hanburry and Serra's paper specifies them as 0.2126, 0.7152, and	*/
	/* 0.0722.								*/
	return _RGBToYUV(rgb, 0.299, 0.114, 0);
}

YUV
RGBToYNUV(RGB rgb)
{
	/* Normalizes the U and V channels for better use in DCS calculations.	*/
	return _RGBToYUV(rgb, 0.299, 0.114, 1);
}

RGB
YUVToRGB(YUV yuv)
{
	return _YUVToRGB(yuv, 0.299, 0.114, 0);
}

RGB
YNUVToRGB(YUV yuv)
{
	/* Assumes the the U and V channels were normalized.	*/
	return _YUVToRGB(yuv, 0.299, 0.114, 1);
}
/* DaVinci's Shadow customization - custom colourspace conversions - part 2 of 2 - End */

#define MAX_INTENSITY 255

RGB
CMYToRGB(cmy)
CMY cmy;

{
    RGB rgb;

    rgb.r = MAX_INTENSITY - cmy.c;
    rgb.g = MAX_INTENSITY - cmy.m;
    rgb.b = MAX_INTENSITY - cmy.y;
    return rgb;
}

/*
 * Convert an RGB to CMY.
 */

CMY
RGBToCMY(rgb)
RGB rgb;

{
    CMY cmy;

    cmy.c = MAX_INTENSITY - rgb.r;
    cmy.m = MAX_INTENSITY - rgb.g;
    cmy.y = MAX_INTENSITY - rgb.b;
    return cmy;
}

KCMY 
RGBtoKCMY(RGB rgb)
{
	KCMY kcmy;

    kcmy.c = MAX_INTENSITY - rgb.r;
    kcmy.m = MAX_INTENSITY - rgb.g;
    kcmy.y = MAX_INTENSITY - rgb.b;

	kcmy.k = min(min(kcmy.c, kcmy.m), kcmy.y);

	if (kcmy.k > 0) {
		kcmy.c = kcmy.c - kcmy.k;
		kcmy.m = kcmy.m - kcmy.k;
		kcmy.y = kcmy.y - kcmy.k;
	}
	return(kcmy);
}

RGB 
KCMYtoRGB(KCMY kcmy)
{
	RGB rgb;

	rgb.r = MAX_INTENSITY - (kcmy.c + kcmy.k);
	rgb.g = MAX_INTENSITY - (kcmy.m + kcmy.k);
	rgb.b = MAX_INTENSITY - (kcmy.y + kcmy.k);

	return(rgb);
}

/**
 ** Compute entropy of an image
 **/

Var *
ff_entropy(vfuncptr func, Var * arg)
{
	Var *obj = NULL;
	int i, dsize, count;
	void *a, *b, *data;
	int format, nbytes;
	float p, ent = 0;
	int (*cmp)(const void *, const void *);

	Alist alist[2];
	alist[0] = make_alist( "object",    ID_VAL,    NULL,     &obj);
	alist[1].name = NULL;

	if (parse_args(func, arg, alist) == 0) return(NULL);

	if (obj == NULL) {
		parse_error("%s: No object specified\n", func->name);
		return(NULL);
	}

	dsize = V_DSIZE(obj);
	format = V_FORMAT(obj);
	nbytes = NBYTES(V_FORMAT(obj));
	data = memdup(V_DATA(obj), dsize * NBYTES(V_FORMAT(obj)));

	switch(format) {
		case BYTE:          cmp = cmp_byte; break;
		case SHORT:         cmp = cmp_short; break;
		case INT:           cmp = cmp_int; break;
		case FLOAT:         cmp = cmp_float; break;
		case DOUBLE:        cmp = cmp_double; break;
	}
	qsort(data, V_DSIZE(obj), NBYTES(format), cmp);

	a = data;
	count = 0;
	for (i = 0 ; i < dsize ; i++) {
		b = ((char *)a) + nbytes;
		if (!cmp(a, b) && (i+1) < dsize) {
			count++;
		} else {
			p  = (float)(count+1)/(float)dsize;
			ent += p * log(p) / M_LN2;
			count = 0;
		}
		a = b;
	} 
	ent = -ent;
	free(data);
	return(newVal(BSQ, 1, 1, 1, FLOAT, memdup(&ent, sizeof(FLOAT))));
}
