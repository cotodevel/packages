///////////////////////////////////////////////////////////////////////
//	    C Implementation of Wu's Color Quantizer (v. 2)
//	    (see Graphics Gems vol. II, pp. 126-133)
//
// Author:	Xiaolin Wu
// Dept. of Computer Science
// Univ. of Western Ontario
// London, Ontario N6A 5B7
// wu@csd.uwo.ca
// 
// Algorithm: Greedy orthogonal bipartition of RGB space for variance
// 	   minimization aided by inclusion-exclusion tricks.
// 	   For speed no nearest neighbor search is done. Slightly
// 	   better performance can be expected by more sophisticated
// 	   but more expensive versions.
// 
// The author thanks Tom Lane at Tom_Lane@G.GP.CS.CMU.EDU for much of
// additional documentation and a cure to a previous bug.
// 
// Free to distribute, comments and suggestions are appreciated.
///////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////
// History
// -------
// July 2000:  C++ Implementation of Wu's Color Quantizer
//             and adaptation for the FreeImage 2 Library
//             Author: Hervé Drolon (drolon@infonie.fr)
// March 2004: Adaptation for the FreeImage 3 library (port to big endian processors)
//             Author: Hervé Drolon (drolon@infonie.fr)
///////////////////////////////////////////////////////////////////////

#include "Quantizers.h"
#include "FreeImage.h"
#include "Utilities.h"

///////////////////////////////////////////////////////////////////////

// Size of a 3D array : 33 x 33 x 33
#define SIZE_3D	35937

// 3D array indexation
#define INDEX(r, g, b)	((r << 10) + (r << 6) + r + (g << 5) + g + b)

#define MAXCOLOR	256

// Constructor / Destructor

WuQuantizer::WuQuantizer(FIBITMAP *dib) {
	width = FreeImage_GetWidth(dib);
	height = FreeImage_GetHeight(dib);
	pitch = FreeImage_GetPitch(dib);
	m_dib = dib;

	gm2 = NULL;
	wt = mr = mg = mb = NULL;
	Qadd = NULL;

	// Allocate 3D arrays
	gm2 = (float*)malloc(SIZE_3D * sizeof(float));
	wt = (LONG*)malloc(SIZE_3D * sizeof(LONG));
	mr = (LONG*)malloc(SIZE_3D * sizeof(LONG));
	mg = (LONG*)malloc(SIZE_3D * sizeof(LONG));
	mb = (LONG*)malloc(SIZE_3D * sizeof(LONG));

	// Allocate Qadd
	Qadd = (WORD *)malloc(sizeof(WORD) * width * height);

	if(!gm2 || !wt || !mr || !mg || !mb || !Qadd) {
		if(gm2)	free(gm2);
		if(wt)	free(wt);
		if(mr)	free(mr);
		if(mg)	free(mg);
		if(mb)	free(mb);
		if(Qadd)  free(Qadd);
		throw FI_MSG_ERROR_MEMORY;
	}
	memset(gm2, 0, SIZE_3D * sizeof(float));
	memset(wt, 0, SIZE_3D * sizeof(LONG));
	memset(mr, 0, SIZE_3D * sizeof(LONG));
	memset(mg, 0, SIZE_3D * sizeof(LONG));
	memset(mb, 0, SIZE_3D * sizeof(LONG));
	memset(Qadd, 0, sizeof(WORD) * width * height);
}


// Histogram is in elements 1..HISTSIZE along each axis,
// element 0 is for base or marginal value
// NB: these must start out 0!

// Build 3-D color histogram of counts, r/g/b, c^2
void 
WuQuantizer::Hist3D(LONG *vwt, LONG *vmr, LONG *vmg, LONG *vmb, float *m2, int ReserveSize, RGBQUAD *ReservePalette) {
	int ind = 0;
	int inr, ing, inb, table[256];
	int i;
	unsigned y, x;

	for(i = 0; i < 256; i++)
		table[i] = i * i;

	if (FreeImage_GetBPP(m_dib) == 24) {
		for(y = 0; y < height; y++) {
			BYTE *bits = FreeImage_GetScanLine(m_dib, y);

			for(x = 0; x < width; x++)	{
				inr = (bits[FI_RGBA_RED] >> 3) + 1;
				ing = (bits[FI_RGBA_GREEN] >> 3) + 1;
				inb = (bits[FI_RGBA_BLUE] >> 3) + 1;
				ind = INDEX(inr, ing, inb);
				Qadd[y*width + x] = (WORD)ind;
				// [inr][ing][inb]
				vwt[ind]++;
				vmr[ind] += bits[FI_RGBA_RED];
				vmg[ind] += bits[FI_RGBA_GREEN];
				vmb[ind] += bits[FI_RGBA_BLUE];
				m2[ind] += (float)(table[bits[FI_RGBA_RED]] + table[bits[FI_RGBA_GREEN]] + table[bits[FI_RGBA_BLUE]]);
				bits += 3;
			}
		}
	} else {
		for(y = 0; y < height; y++) {
			BYTE *bits = FreeImage_GetScanLine(m_dib, y);

			for(x = 0; x < width; x++)	{
				inr = (bits[FI_RGBA_RED] >> 3) + 1;
				ing = (bits[FI_RGBA_GREEN] >> 3) + 1;
				inb = (bits[FI_RGBA_BLUE] >> 3) + 1;
				ind = INDEX(inr, ing, inb);
				Qadd[y*width + x] = (WORD)ind;
				// [inr][ing][inb]
				vwt[ind]++;
				vmr[ind] += bits[FI_RGBA_RED];
				vmg[ind] += bits[FI_RGBA_GREEN];
				vmb[ind] += bits[FI_RGBA_BLUE];
				m2[ind] += (float)(table[bits[FI_RGBA_RED]] + table[bits[FI_RGBA_GREEN]] + table[bits[FI_RGBA_BLUE]]);
				bits += 4;
			}
		}
	}

	if( ReserveSize > 0 ) {
		int max = 0;
		for(i = 0; i < SIZE_3D; i++) {
			if( vwt[i] > max ) max = vwt[i];
		}
		max++;
		for(i = 0; i < ReserveSize; i++) {
			inr = (ReservePalette[i].rgbRed >> 3) + 1;
			ing = (ReservePalette[i].rgbGreen >> 3) + 1;
			inb = (ReservePalette[i].rgbBlue >> 3) + 1;
			ind = INDEX(inr, ing, inb);
			wt[ind] = max;
			mr[ind] = max * ReservePalette[i].rgbRed;
			mg[ind] = max * ReservePalette[i].rgbGreen;
			mb[ind] = max * ReservePalette[i].rgbBlue;
			gm2[ind] = (float)max * (float)(table[ReservePalette[i].rgbRed] + table[ReservePalette[i].rgbGreen] + table[ReservePalette[i].rgbBlue]);
		}
	}
}

// Wu Quantization algorithm
FIBITMAP *
WuQuantizer::Quantize(int PaletteSize, int ReserveSize, RGBQUAD *ReservePalette) {
	BYTE *tag = NULL;

	try {
		Box	cube[MAXCOLOR];
		int	next;
		LONG i, weight;
		int k;
		float vv[MAXCOLOR], temp;
		
		// Compute 3D histogram

		Hist3D(wt, mr, mg, mb, gm2, ReserveSize, ReservePalette);

		// Compute moments

		M3D(wt, mr, mg, mb, gm2);

		cube[0].r0 = cube[0].g0 = cube[0].b0 = 0;
		cube[0].r1 = cube[0].g1 = cube[0].b1 = 32;
		next = 0;

		for (i = 1; i < PaletteSize; i++) {
			if(Cut(&cube[next], &cube[i])) {
				// volume test ensures we won't try to cut one-cell box
				vv[next] = (cube[next].vol > 1) ? Var(&cube[next]) : 0;
				vv[i] = (cube[i].vol > 1) ? Var(&cube[i]) : 0;
			} else {
				  vv[next] = 0.0;   // don't try to split this box again
				  i--;              // didn't create box i
			}

			next = 0; temp = vv[0];

			for (k = 1; k <= i; k++) {
				if (vv[k] > temp) {
					temp = vv[k]; next = k;
				}
			}

			if (temp <= 0.0) {
				  PaletteSize = i + 1;

				  // Error: "Only got 'PaletteSize' boxes"

				  break;
			}
		}

		// Partition done

		// the space for array gm2 can be freed now

		free(gm2);

		gm2 = NULL;

		// Allocate a new dib

		FIBITMAP *new_dib = FreeImage_Allocate(width, height, 8);

		if (new_dib == NULL) {
			throw FI_MSG_ERROR_MEMORY;
		}

		// create an optimized palette

		RGBQUAD *new_pal = FreeImage_GetPalette(new_dib);

		tag = (BYTE*) malloc(SIZE_3D * sizeof(BYTE));
		if (tag == NULL) {
			throw FI_MSG_ERROR_MEMORY;
		}
		memset(tag, 0, SIZE_3D * sizeof(BYTE));

		for (k = 0; k < PaletteSize ; k++) {
			Mark(&cube[k], k, tag);
			weight = Vol(&cube[k], wt);

			if (weight) {
				new_pal[k].rgbRed	= (BYTE)(((float)Vol(&cube[k], mr) / (float)weight) + 0.5f);
				new_pal[k].rgbGreen = (BYTE)(((float)Vol(&cube[k], mg) / (float)weight) + 0.5f);
				new_pal[k].rgbBlue	= (BYTE)(((float)Vol(&cube[k], mb) / (float)weight) + 0.5f);
			} else {
				// Error: bogus box 'k'

				new_pal[k].rgbRed = new_pal[k].rgbGreen = new_pal[k].rgbBlue = 0;		
			}
		}

		int npitch = FreeImage_GetPitch(new_dib);

		for (unsigned y = 0; y < height; y++) {
			BYTE *new_bits = FreeImage_GetBits(new_dib) + (y * npitch);

			for (unsigned x = 0; x < width; x++) {
				new_bits[x] = tag[Qadd[y*width + x]];
			}
		}

		// output 'new_pal' as color look-up table contents,
		// 'new_bits' as the quantized image (array of table addresses).

		free(tag);

		return (FIBITMAP*) new_dib;
	} catch(...) {
		free(tag);
	}

	return NULL;
}
