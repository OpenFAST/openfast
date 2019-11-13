/* Copyright 2002-2010 Paul Hsieh
 * This file is part of Bstrlib.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *    1. Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *    2. Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *
 *    3. Neither the name of bstrlib nor the names of its contributors may be
 *       used to endorse or promote products derived from this software
 *       without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * Alternatively, the contents of this file may be used under the terms of
 * GNU General Public License Version 2 (the "GPL").
 */

/*
 * This file is the core module for implementing the bstring functions.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#if defined (_MSC_VER)
/* These warnings from MSVC++ are totally pointless. */
# define _CRT_SECURE_NO_WARNINGS
#endif

#include <stdio.h>
#include <stddef.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "bstrlib.h"

/* Just a length safe wrapper for memmove. */

#define bBlockCopy(D, S, L) \
do { \
	if ((L) > 0) { \
		memmove((D), (S), (L)); \
	} \
} while (0);

/**
 * Compute the snapped size for a given requested size.
 *
 * By snapping to powers of 2 like this, repeated reallocations are avoided.
 * */
static int
snapUpSize(int i) {
	if (i < 8) {
		i = 8;
	} else {
		unsigned int j;
		j = (unsigned int)i;
		j |= (j >>  1);
		j |= (j >>  2);
		j |= (j >>  4);
		j |= (j >>  8);	/* Ok, since int >= 16 bits */
#if (UINT_MAX != 0xffff)
		j |= (j >> 16);	/* For 32 bit int systems */
#if (UINT_MAX > 0xffffffffUL)
		j |= (j >> 32);	/* For 64 bit int systems */
#endif
#endif
		/* Least power of two greater than i */
		j++;
		if ((int)j >= i) {
			i = (int)j;
		}
	}
	return i;
}

int
balloc(bstring b, int olen)
{
	int len;
	if (b == NULL || b->data == NULL ||
	    b->slen < 0 || b->mlen <= 0 ||
	    b->mlen < b->slen || olen <= 0) {
		return BSTR_ERR;
	}
	if (olen >= b->mlen) {
		unsigned char *x;
		if ((len = snapUpSize(olen)) <= b->mlen) {
			return BSTR_OK;
		}
		/* Assume probability of a non-moving realloc is 0.125 */
		if (7 * b->mlen < 8 * b->slen) {
			/* If slen is close to mlen in size then use realloc
			 * to reduce the memory defragmentation
			 */
retry:
			x = realloc(b->data, len);
			if (x == NULL) {
				/* Since we failed, try mallocating the tighest
				 * possible mallocation
				 */
				len = olen;
				x = realloc(b->data, len);
				if (!x) {
					return BSTR_ERR;
				}
			}
		} else {
			/* If slen is not close to mlen then avoid the penalty
			 * of copying the extra bytes that are mallocated, but
			 * not considered part of the string
			 */
			x = malloc(len);
			if (!x) {
				/* Perhaps there is no available memory for the
				 * two mallocations to be in memory at once
				 */
				goto retry;
			} else {
				if (b->slen) {
					memcpy(x, b->data, b->slen);
				}
				free(b->data);
			}
		}
		b->data = x;
		b->mlen = len;
		b->data[b->slen] = (unsigned char)'\0';
	}
	return BSTR_OK;
}

int
ballocmin(bstring b, int len)
{
	unsigned char *s;
	if (b == NULL || b->data == NULL ||
	    (b->slen + 1) < 0 || b->mlen <= 0 ||
	    b->mlen < b->slen || len <= 0) {
		return BSTR_ERR;
	}
	if (len < b->slen + 1) {
		len = b->slen + 1;
	}
	if (len != b->mlen) {
		s = realloc(b->data, (size_t)len);
		if (NULL == s) {
			return BSTR_ERR;
		}
		s[b->slen] = (unsigned char)'\0';
		b->data = s;
		b->mlen = len;
	}

	return BSTR_OK;
}

bstring
bfromcstr(const char *str)
{
	bstring b;
	int i;
	size_t j;
	if (str == NULL) {
		return NULL;
	}
	j = strlen(str);
	i = snapUpSize((int)(j + (2 - (j != 0))));
	if (i <= (int)j) {
		return NULL;
	}
	b = malloc(sizeof(struct tagbstring));
	if (!b) {
		return NULL;
	}
	b->slen = (int)j;
	b->mlen = i;
	b->data = malloc(b->mlen);
	if (!b->data) {
		free (b);
		return NULL;
	}
	memcpy(b->data, str, j + 1);
	return b;
}

bstring
bfromcstralloc(int mlen, const char *str)
{
	bstring b;
	int i;
	size_t j;
	if (str == NULL) {
		return NULL;
	}
	j = strlen(str);
	i = snapUpSize((int)(j + (2 - (j != 0))));
	if (i <= (int) j) {
		return NULL;
	}
	b = malloc(sizeof(struct tagbstring));
	if (b == NULL) {
		return NULL;
	}
	b->slen = (int)j;
	if (i < mlen) {
		i = mlen;
	}
	b->mlen = i;
	b->data = malloc(b->mlen);
	if (!b->data) {
		free(b);
		return NULL;
	}
	memcpy(b->data, str, j + 1);
	return b;
}

bstring
blk2bstr(const void *blk, int len)
{
	bstring b;
	int i;
	if (blk == NULL || len < 0) {
		return NULL;
	}
	b = malloc(sizeof(struct tagbstring));
	if (b == NULL) {
		return NULL;
	}
	b->slen = len;
	i = len + (2 - (len != 0));
	i = snapUpSize(i);
	b->mlen = i;
	b->data = malloc(b->mlen);
	if (!b->data) {
		free(b);
		return NULL;
	}
	if (len > 0) {
		memcpy(b->data, blk, len);
	}
	b->data[len] = (unsigned char)'\0';
	return b;
}

char *
bstr2cstr(const bstring b, char z)
{
	int i, l;
	char *r;
	if (!b || b->slen < 0 || !b->data) {
		return NULL;
	}
	l = b->slen;
	r = malloc((size_t)(l + 1));
	if (r == NULL) {
		return r;
	}
	for (i = 0; i < l; i ++) {
		r[i] = (char)((b->data[i] == '\0') ? z : (char)(b->data[i]));
	}
	r[l] = (unsigned char)'\0';
	return r;
}

int
bcstrfree(char *s)
{
	free(s);
	return BSTR_OK;
}

int
bconcat(bstring b0, const bstring b1)
{
	int len, d;
	bstring aux = b1;
	if (!b0 || !b1 || !b0->data || !b1->data) {
		return BSTR_ERR;
	}
	d = b0->slen;
	len = b1->slen;
	if ((d | (b0->mlen - d) | len | (d + len)) < 0) {
		return BSTR_ERR;
	}
	if (b0->mlen <= d + len + 1) {
		ptrdiff_t pd = b1->data - b0->data;
		if (0 <= pd && pd < b0->mlen) {
			aux = bstrcpy(b1);
			if (!aux) {
				return BSTR_ERR;
			}
		}
		if (balloc(b0, d + len + 1) != BSTR_OK) {
			if (aux != b1) {
				bdestroy(aux);
			}
			return BSTR_ERR;
		}
	}
	bBlockCopy(&b0->data[d], &aux->data[0], len);
	b0->data[d + len] = (unsigned char)'\0';
	b0->slen = d + len;
	if (aux != b1) {
		bdestroy(aux);
	}
	return BSTR_OK;
}

int
bconchar(bstring b, char c)
{
	int d;
	if (!b) {
		return BSTR_ERR;
	}
	d = b->slen;
	if ((d | (b->mlen - d)) < 0 || balloc(b, d + 2) != BSTR_OK) {
		return BSTR_ERR;
	}
	b->data[d] = (unsigned char)c;
	b->data[d + 1] = (unsigned char)'\0';
	b->slen++;
	return BSTR_OK;
}

int
bcatcstr(bstring b, const char *s)
{
	char *d;
	int i, l;
	if (b == NULL || b->data == NULL ||
	    b->slen < 0 || b->mlen < b->slen
	    || b->mlen <= 0 || s == NULL) {
		return BSTR_ERR;
	}
	/* Optimistically concatenate directly */
	l = b->mlen - b->slen;
	d = (char *)&b->data[b->slen];
	for (i = 0; i < l; ++i) {
		if ((*d++ = *s++) == '\0') {
			b->slen += i;
			return BSTR_OK;
		}
	}
	b->slen += i;
	/* Need to explicitely resize and concatenate tail */
	return bcatblk(b, s, strlen(s));
}

int
bcatblk(bstring b, const void *s, int len)
{
	int nl;
	if (!b || !b->data ||
	    b->slen < 0 || b->mlen < b->slen ||
	    b->mlen <= 0 || !s || len < 0) {
		return BSTR_ERR;
	}
	if (0 > (nl = b->slen + len)) {
		/* Overflow? */
		return BSTR_ERR;
	}
	if (b->mlen <= nl && 0 > balloc(b, nl + 1)) {
		return BSTR_ERR;
	}
	bBlockCopy(&b->data[b->slen], s, len);
	b->slen = nl;
	b->data[nl] = (unsigned char)'\0';
	return BSTR_OK;
}

bstring
bstrcpy(const bstring b)
{
	bstring b0;
	int i, j;
	/* Attempted to copy an invalid string? */
	if (!b || b->slen < 0 || !b->data) {
		return NULL;
	}
	b0 = malloc(sizeof(struct tagbstring));
	if (!b0) {
		/* Unable to mallocate memory for string header */
		return NULL;
	}
	i = b->slen;
	j = snapUpSize(i + 1);
	b0->data = malloc(j);
	if (b0->data == NULL) {
		j = i + 1;
		b0->data = (unsigned char *)malloc(j);
		if (b0->data == NULL) {
			/* Unable to mallocate memory for string data */
			free(b0);
			return NULL;
		}
	}
	b0->mlen = j;
	b0->slen = i;
	if (i) {
		memcpy(b0->data, b->data, i);
	}
	b0->data[b0->slen] = (unsigned char)'\0';
	return b0;
}

int
bassign(bstring a, const bstring b)
{
	if (!b || !b->data || b->slen < 0) {
		return BSTR_ERR;
	}
	if (b->slen != 0) {
		if (balloc(a, b->slen) != BSTR_OK) {
			return BSTR_ERR;
		}
		memmove(a->data, b->data, b->slen);
	} else {
		if (!a || !a->data ||
		    a->mlen < a->slen ||
		    a->slen < 0 || a->mlen == 0) {
			return BSTR_ERR;
		}
	}
	a->data[b->slen] = (unsigned char)'\0';
	a->slen = b->slen;
	return BSTR_OK;
}

int
bassignmidstr(bstring a, const bstring b, int left, int len)
{
	if (!b || !b->data || b->slen < 0) {
		return BSTR_ERR;
	}
	if (left < 0) {
		len += left;
		left = 0;
	}
	if (len > b->slen - left) {
		len = b->slen - left;
	}
	if (!a || !a->data ||
	    a->mlen < a->slen ||
	    a->slen < 0 || a->mlen == 0) {
		return BSTR_ERR;
	}
	if (len > 0) {
		if (balloc(a, len) != BSTR_OK) {
			return BSTR_ERR;
		}
		memmove(a->data, b->data + left, len);
		a->slen = len;
	} else {
		a->slen = 0;
	}
	a->data[a->slen] = (unsigned char)'\0';
	return BSTR_OK;
}

int
bassigncstr(bstring a, const char *str)
{
	int i;
	size_t len;
	if (!a || !a->data ||
	    a->mlen < a->slen || a->slen < 0 ||
	    a->mlen == 0 || !str) {
		return BSTR_ERR;
	}
	for (i = 0; i < a->mlen; ++i) {
		if ('\0' == (a->data[i] = str[i])) {
			a->slen = i;
			return BSTR_OK;
		}
	}
	a->slen = i;
	len = strlen(str + i);
	if (len > INT_MAX || i + len + 1 > INT_MAX ||
	    0 > balloc(a, (int)(i + len + 1))) {
		return BSTR_ERR;
	}
	bBlockCopy(a->data + i, str + i, (size_t)len + 1);
	a->slen += (int)len;
	return BSTR_OK;
}

int
bassignblk(bstring a, const void *s, int len)
{
	if (!a || !a->data ||
	    a->mlen < a->slen || a->slen < 0 ||
	    a->mlen == 0 || !s ||
	    len + 1 < 1) {
		return BSTR_ERR;
	}
	if (len + 1 > a->mlen && 0 > balloc(a, len + 1)) {
		return BSTR_ERR;
	}
	bBlockCopy(a->data, s, len);
	a->data[len] = (unsigned char)'\0';
	a->slen = len;
	return BSTR_OK;
}

int
btrunc(bstring b, int n)
{
	if (n < 0 || !b ||
	    !b->data || b->mlen < b->slen ||
	    b->slen < 0 || b->mlen <= 0) {
		return BSTR_ERR;
	}
	if (b->slen > n) {
		b->slen = n;
		b->data[n] = (unsigned char)'\0';
	}
	return BSTR_OK;
}

#define upcase(c) \
	(toupper((unsigned char)c))

#define downcase(c) \
	(tolower((unsigned char)c))

#define wspace(c) \
	(isspace((unsigned char)c))

int
btoupper(bstring b)
{
	int i, len;
	if (b == NULL || b->data == NULL ||
	    b->mlen < b->slen || b->slen < 0 ||
	    b->mlen <= 0) {
		return BSTR_ERR;
	}
	for (i = 0, len = b->slen; i < len; i++) {
		b->data[i] = (unsigned char)upcase(b->data[i]);
	}
	return BSTR_OK;
}

int
btolower(bstring b)
{
	int i, len;
	if (b == NULL || b->data == NULL ||
	    b->mlen < b->slen || b->slen < 0 ||
	    b->mlen <= 0) {
		return BSTR_ERR;
	}
	for (i = 0, len = b->slen; i < len; i++) {
		b->data[i] = (unsigned char)downcase(b->data[i]);
	}
	return BSTR_OK;
}

int
bstricmp(const bstring b0, const bstring b1)
{
	int i, v, n;
	if (bdata (b0) == NULL || b0->slen < 0 ||
	    bdata (b1) == NULL || b1->slen < 0) {
		return SHRT_MIN;
	}
	if ((n = b0->slen) > b1->slen) {
		n = b1->slen;
	} else if (b0->slen == b1->slen && b0->data == b1->data) {
		return BSTR_OK;
	}
	for (i = 0; i < n; i ++) {
		v  = (char)downcase(b0->data[i]) - (char)downcase(b1->data[i]);
		if (0 != v) {
			return v;
		}
	}
	if (b0->slen > n) {
		v = (char)downcase(b0->data[n]);
		if (v) {
			return v;
		}
		return UCHAR_MAX + 1;
	}
	if (b1->slen > n) {
		v = - (char)downcase(b1->data[n]);
		if (v) {
			return v;
		}
		return -(int)(UCHAR_MAX + 1);
	}
	return BSTR_OK;
}

int
bstrnicmp(const bstring b0, const bstring b1, int n)
{
	int i, v, m;
	if (bdata (b0) == NULL || b0->slen < 0 ||
	    bdata (b1) == NULL || b1->slen < 0 ||
	    n < 0) {
		return SHRT_MIN;
	}
	m = n;
	if (m > b0->slen) {
		m = b0->slen;
	}
	if (m > b1->slen) {
		m = b1->slen;
	}
	if (b0->data != b1->data) {
		for (i = 0; i < m; i ++) {
			v  = (char)downcase(b0->data[i]);
			v -= (char)downcase(b1->data[i]);
			if (v != 0) {
				return b0->data[i] - b1->data[i];
			}
		}
	}
	if (n == m || b0->slen == b1->slen) {
		return BSTR_OK;
	}
	if (b0->slen > m) {
		v = (char)downcase(b0->data[m]);
		if (v) {
			return v;
		}
		return UCHAR_MAX + 1;
	}
	v = - (char)downcase(b1->data[m]);
	if (v) {
		return v;
	}
	return -(int)(UCHAR_MAX + 1);
}

int
biseqcaseless(const bstring b0, const bstring b1)
{
	int i, n;
	if (bdata (b0) == NULL || b0->slen < 0 ||
	    bdata (b1) == NULL || b1->slen < 0) {
		return BSTR_ERR;
	}
	if (b0->slen != b1->slen) {
		return BSTR_OK;
	}
	if (b0->data == b1->data || b0->slen == 0) {
		return 1;
	}
	for (i = 0, n = b0->slen; i < n; i++) {
		if (b0->data[i] != b1->data[i]) {
			unsigned char c = (unsigned char)downcase(b0->data[i]);
			if (c != (unsigned char)downcase(b1->data[i])) {
				return 0;
			}
		}
	}
	return 1;
}

int
bisstemeqcaselessblk(const bstring b0, const void *blk, int len)
{
	int i;
	if (bdata(b0) == NULL || b0->slen < 0 || NULL == blk || len < 0) {
		return BSTR_ERR;
	}
	if (b0->slen < len) {
		return BSTR_OK;
	}
	if (b0->data == (const unsigned char *)blk || len == 0) {
		return 1;
	}
	for (i = 0; i < len; i++) {
		if (b0->data[i] != ((const unsigned char *)blk)[i]) {
			if (downcase(b0->data[i]) !=
			    downcase(((const unsigned char *)blk)[i])) {
				return 0;
			}
		}
	}
	return 1;
}

int
bltrimws(bstring b)
{
	int i, len;
	if (!b || !b->data ||
	    b->mlen < b->slen || b->slen < 0 ||
	    b->mlen <= 0) {
		return BSTR_ERR;
	}
	for (len = b->slen, i = 0; i < len; i++) {
		if (!wspace(b->data[i])) {
			return bdelete(b, 0, i);
		}
	}
	b->data[0] = (unsigned char) '\0';
	b->slen = 0;
	return BSTR_OK;
}

int
brtrimws(bstring b)
{
	int i;
	if (b == NULL ||
	    b->data == NULL ||
	    b->mlen < b->slen ||
	    b->slen < 0 ||
	    b->mlen <= 0) {
		return BSTR_ERR;
	}
	for (i = b->slen - 1; i >= 0; i--) {
		if (!wspace(b->data[i])) {
			if (b->mlen > i) {
				b->data[i + 1] = (unsigned char)'\0';
			}
			b->slen = i + 1;
			return BSTR_OK;
		}
	}
	b->data[0] = (unsigned char)'\0';
	b->slen = 0;
	return BSTR_OK;
}

int
btrimws(bstring b)
{
	int i, j;
	if (b == NULL ||
	    b->data == NULL ||
	    b->mlen < b->slen ||
	    b->slen < 0 ||
	    b->mlen <= 0) {
		return BSTR_ERR;
	}
	for (i = b->slen - 1; i >= 0; i--) {
		if (!wspace(b->data[i])) {
			if (b->mlen > i) {
				b->data[i + 1] = (unsigned char)'\0';
			}
			b->slen = i + 1;
			for (j = 0; wspace (b->data[j]); j++)
				;
			return bdelete(b, 0, j);
		}
	}
	b->data[0] = (unsigned char)'\0';
	b->slen = 0;
	return BSTR_OK;
}

int
biseq(const bstring b0, const bstring b1)
{
	if (!b0 || !b1 ||
	    !b0->data || !b1->data ||
	    b0->slen < 0 || b1->slen < 0) {
		return BSTR_ERR;
	}
	if (b0->slen != b1->slen) {
		return BSTR_OK;
	}
	if (b0->data == b1->data || b0->slen == 0) {
		return 1;
	}
	return !memcmp(b0->data, b1->data, b0->slen);
}

int
bisstemeqblk(const bstring b0, const void *blk, int len)
{
	int i;
	if (!bdata(b0) || b0->slen < 0 || !blk || len < 0) {
		return BSTR_ERR;
	}
	if (b0->slen < len) {
		return BSTR_OK;
	}
	if (b0->data == (const unsigned char *)blk || len == 0) {
		return 1;
	}
	for (i = 0; i < len; i ++) {
		if (b0->data[i] != ((const unsigned char *)blk)[i]) {
			return BSTR_OK;
		}
	}
	return 1;
}

int
biseqcstr(const bstring b, const char *s)
{
	int i;
	if (!b || !s  || !b->data || b->slen < 0) {
		return BSTR_ERR;
	}
	for (i = 0; i < b->slen; i++) {
		if (s[i] == '\0' || b->data[i] != (unsigned char)s[i]) {
			return BSTR_OK;
		}
	}
	return s[i] == '\0';
}

int
biseqcstrcaseless(const bstring b, const char *s)
{
	int i;
	if (!b || !s || !b->data || b->slen < 0) {
		return BSTR_ERR;
	}
	for (i = 0; i < b->slen; i++) {
		if (s[i] == '\0' || (b->data[i] != (unsigned char)s[i] &&
		    downcase(b->data[i]) != (unsigned char)downcase(s[i]))) {
			return BSTR_OK;
		}
	}
	return s[i] == '\0';
}

int
bstrcmp(const bstring b0, const bstring b1)
{
	int i, v, n;
	if (!b0 || !b1 || !b0->data || !b1->data ||
	    b0->slen < 0 || b1->slen < 0) {
		return SHRT_MIN;
	}
	n = b0->slen;
	if (n > b1->slen) {
		n = b1->slen;
	}
	if (b0->slen == b1->slen && (b0->data == b1->data || b0->slen == 0)) {
		return BSTR_OK;
	}
	for (i = 0; i < n; i ++) {
		v = ((char)b0->data[i]) - ((char)b1->data[i]);
		if (v != 0) {
			return v;
		}
		if (b0->data[i] == (unsigned char)'\0') {
			return BSTR_OK;
		}
	}
	if (b0->slen > n) {
		return 1;
	}
	if (b1->slen > n) {
		return -1;
	}
	return BSTR_OK;
}

int
bstrncmp(const bstring b0, const bstring b1, int n)
{
	int i, v, m;
	if (!b0 || !b1 || !b0->data || !b1->data ||
	    b0->slen < 0 || b1->slen < 0) {
		return SHRT_MIN;
	}
	m = n;
	if (m > b0->slen) {
		m = b0->slen;
	}
	if (m > b1->slen) {
		m = b1->slen;
	}
	if (b0->data != b1->data) {
		for (i = 0; i < m; i++) {
			v = ((char)b0->data[i]) - ((char)b1->data[i]);
			if (v != 0) {
				return v;
			}
			if (b0->data[i] == (unsigned char)'\0') {
				return BSTR_OK;
			}
		}
	}
	if (n == m || b0->slen == b1->slen) {
		return BSTR_OK;
	}
	if (b0->slen > m) {
		return 1;
	}
	return -1;
}

bstring
bmidstr(const bstring b, int left, int len)
{
	if (b == NULL || b->slen < 0 || b->data == NULL) {
		return NULL;
	}
	if (left < 0) {
		len += left;
		left = 0;
	}
	if (len > b->slen - left) {
		len = b->slen - left;
	}
	if (len <= 0) {
		return bfromcstr("");
	}
	return blk2bstr(b->data + left, len);
}

int
bdelete(bstring b, int pos, int len)
{
	/* Clamp to left side of bstring */
	if (pos < 0) {
		len += pos;
		pos = 0;
	}
	if (len < 0 || b == NULL || b->data == NULL || b->slen < 0 ||
	    b->mlen < b->slen || b->mlen <= 0) {
		return BSTR_ERR;
	}
	if (len > 0 && pos < b->slen) {
		if (pos + len >= b->slen) {
			b->slen = pos;
		} else {
			bBlockCopy((char *)(b->data + pos),
				   (char *)(b->data + pos + len),
				   b->slen - (pos+len));
			b->slen -= len;
		}
		b->data[b->slen] = (unsigned char)'\0';
	}
	return BSTR_OK;
}

int
bdestroy(bstring b)
{
	if (b == NULL || b->slen < 0 ||
	    b->mlen <= 0 || b->mlen < b->slen ||
	    b->data == NULL) {
		return BSTR_ERR;
	}
	free(b->data);
	/* In case there is any stale usage, there is one more chance to
	 * notice this error.
	 */
	b->slen = -1;
	b->mlen = -__LINE__;
	b->data = NULL;
	free(b);
	return BSTR_OK;
}

int
binstr(const bstring b1, int pos, const bstring b2)
{
	int j, ii, ll, lf;
	unsigned char *d0;
	unsigned char c0;
	register unsigned char *d1;
	register unsigned char c1;
	register int i;
	if (b1 == NULL || b1->data == NULL || b1->slen < 0 ||
	    b2 == NULL || b2->data == NULL || b2->slen < 0) {
		return BSTR_ERR;
	}
	if (b1->slen == pos) {
		return (b2->slen == 0) ? pos : BSTR_ERR;
	}
	if (b1->slen < pos || pos < 0) {
		return BSTR_ERR;
	}
	if (b2->slen == 0) {
		return pos;
	}
	/* No space to find such a string? */
	if ((lf = b1->slen - b2->slen + 1) <= pos) {
		return BSTR_ERR;
	}
	/* An obvious alias case */
	if (b1->data == b2->data && pos == 0) {
		return 0;
	}
	i = pos;
	d0 = b2->data;
	d1 = b1->data;
	ll = b2->slen;
	/* Peel off the b2->slen == 1 case */
	c0 = d0[0];
	if (1 == ll) {
		for (; i < lf; i++) {
			if (c0 == d1[i]) {
				return i;
			}
		}
		return BSTR_ERR;
	}
	c1 = c0;
	j = 0;
	lf = b1->slen - 1;
	ii = -1;
	if (i < lf) {
		do {
			/* Unrolled current character test */
			if (c1 != d1[i]) {
				if (c1 != d1[1+i]) {
					i += 2;
					continue;
				}
				i++;
			}
			/* Take note if this is the start of a potential
			 * match
			 */
			if (0 == j) {
				ii = i;
			}
			/* Shift the test character down by one */
			j++;
			i++;
			/* If this isn't past the last character continue */
			if (j < ll) {
				c1 = d0[j];
				continue;
			}
N0:
			/* If no characters mismatched, then we matched */
			if (i == ii + j) {
				return ii;
			}
			/* Shift back to the beginning */
			i -= j;
			j  = 0;
			c1 = c0;
		} while (i < lf);
	}
	/* Deal with last case if unrolling caused a misalignment */
	if (i == lf && ll == j + 1 && c1 == d1[i]) {
		goto N0;
	}
	return BSTR_ERR;
}

int
binstrr(const bstring b1, int pos, const bstring b2)
{
	int j, i, l;
	unsigned char *d0, *d1;
	if (b1 == NULL || b1->data == NULL || b1->slen < 0 ||
	    b2 == NULL || b2->data == NULL || b2->slen < 0) {
		return BSTR_ERR;
	}
	if (b1->slen == pos && b2->slen == 0) {
		return pos;
	}
	if (b1->slen < pos || pos < 0) {
		return BSTR_ERR;
	}
	if (b2->slen == 0) {
		return pos;
	}
	/* Obvious alias case */
	if (b1->data == b2->data && pos == 0 && b2->slen <= b1->slen) {
		return 0;
	}
	i = pos;
	if ((l = b1->slen - b2->slen) < 0) {
		return BSTR_ERR;
	}
	/* If no space to find such a string then snap back */
	if (l + 1 <= i) {
		i = l;
	}
	j = 0;
	d0 = b2->data;
	d1 = b1->data;
	l  = b2->slen;
	while (1) {
		if (d0[j] == d1[i + j]) {
			j++;
			if (j >= l) {
				return i;
			}
		} else {
			i--;
			if (i < 0) {
				break;
			}
			j = 0;
		}
	}
	return BSTR_ERR;
}

int
binstrcaseless(const bstring b1, int pos, const bstring b2)
{
	int j, i, l, ll;
	unsigned char *d0, *d1;
	if (b1 == NULL || b1->data == NULL || b1->slen < 0 ||
	    b2 == NULL || b2->data == NULL || b2->slen < 0) {
		return BSTR_ERR;
	}
	if (b1->slen == pos) {
		return (b2->slen == 0) ? pos : BSTR_ERR;
	}
	if (b1->slen < pos || pos < 0) {
		return BSTR_ERR;
	}
	if (b2->slen == 0) {
		return pos;
	}
	l = b1->slen - b2->slen + 1;
	/* No space to find such a string? */
	if (l <= pos) {
		return BSTR_ERR;
	}
	/* An obvious alias case */
	if (b1->data == b2->data && pos == 0) {
		return BSTR_OK;
	}
	i = pos;
	j = 0;
	d0 = b2->data;
	d1 = b1->data;
	ll = b2->slen;
	while (1) {
		if (d0[j] == d1[i + j] ||
		    downcase(d0[j]) == downcase (d1[i + j])) {
			j++;
			if (j >= ll) {
				return i;
			}
		} else {
			i ++;
			if (i >= l) {
				break;
			}
			j = 0;
		}
	}
	return BSTR_ERR;
}

int
binstrrcaseless(const bstring b1, int pos, const bstring b2)
{
	int j, i, l;
	unsigned char *d0, *d1;
	if (b1 == NULL || b1->data == NULL || b1->slen < 0 ||
	    b2 == NULL || b2->data == NULL || b2->slen < 0) {
		return BSTR_ERR;
	}
	if (b1->slen == pos && b2->slen == 0) {
		return pos;
	}
	if (b1->slen < pos || pos < 0) {
		return BSTR_ERR;
	}
	if (b2->slen == 0) {
		return pos;
	}
	/* Obvious alias case */
	if (b1->data == b2->data && pos == 0 && b2->slen <= b1->slen) {
		return BSTR_OK;
	}
	i = pos;
	if ((l = b1->slen - b2->slen) < 0) {
		return BSTR_ERR;
	}
	/* If no space to find such a string then snap back */
	if (l + 1 <= i) {
		i = l;
	}
	j = 0;
	d0 = b2->data;
	d1 = b1->data;
	l = b2->slen;
	while (1) {
		if (d0[j] == d1[i + j] ||
		    downcase (d0[j]) == downcase(d1[i + j])){
			j++;
			if (j >= l) {
				return i;
			}
		} else {
			i--;
			if (i < 0) {
				break;
			}
			j = 0;
		}
	}
	return BSTR_ERR;
}

int
bstrchrp(const bstring b, int c, int pos)
{
	unsigned char *p;
	if (b == NULL || b->data == NULL || b->slen <= pos || pos < 0) {
		return BSTR_ERR;
	}
	p = (unsigned char *)memchr((b->data + pos), (unsigned char)c,
				    (b->slen - pos));
	if (p) {
		return (int)(p - b->data);
	}
	return BSTR_ERR;
}

int
bstrrchrp(const bstring b, int c, int pos)
{
	int i;
	if (b == NULL || b->data == NULL || b->slen <= pos || pos < 0) {
		return BSTR_ERR;
	}
	for (i = pos; i >= 0; i--) {
		if (b->data[i] == (unsigned char)c) {
			return i;
		}
	}
	return BSTR_ERR;
}

#ifndef BSTRLIB_AGGRESSIVE_MEMORY_FOR_SPEED_TRADEOFF
#define LONG_LOG_BITS_QTY (3)
#define LONG_BITS_QTY \
	(1 << LONG_LOG_BITS_QTY)
#define LONG_TYPE unsigned char
#define CFCLEN \
	((1 << CHAR_BIT) / LONG_BITS_QTY)
struct charField {
	LONG_TYPE content[CFCLEN];
};
#define testInCharField(cf, c) \
	((cf)->content[(c) >> LONG_LOG_BITS_QTY] & (((long)1) << ((c) & (LONG_BITS_QTY-1))))
#define setInCharField(cf, idx) \
do { \
	unsigned int c = (unsigned int)(idx); \
	(cf)->content[c >> LONG_LOG_BITS_QTY] |= (LONG_TYPE)(1ul << (c & (LONG_BITS_QTY-1))); \
} while (0)
#else /* BSTRLIB_AGGRESSIVE_MEMORY_FOR_SPEED_TRADEOFF */
#define CFCLEN \
	(1 << CHAR_BIT)
struct charField {
	unsigned char content[CFCLEN];
};
#define testInCharField(cf, c) \
	((cf)->content[(unsigned char) (c)])
#define setInCharField(cf, idx) \
	(cf)->content[(unsigned int) (idx)] = ~0
#endif /* BSTRLIB_AGGRESSIVE_MEMORY_FOR_SPEED_TRADEOFF */

/* Convert a bstring to charField */
static int
buildCharField(struct charField *cf, const bstring b)
{
	int i;
	if (b == NULL || b->data == NULL || b->slen <= 0) {
		return BSTR_ERR;
	}
	memset((void *)cf->content, 0, sizeof(struct charField));
	for (i = 0; i < b->slen; i++) {
		setInCharField(cf, b->data[i]);
	}
	return BSTR_OK;
}

static void
invertCharField (struct charField *cf)
{
	int i;
	for (i = 0; i < CFCLEN; i++) {
		cf->content[i] = ~cf->content[i];
	}
}

/* Inner engine for binchr */
static int
binchrCF(const unsigned char *data, int len, int pos,
	 const struct charField *cf)
{
	int i;
	for (i = pos; i < len; i++) {
		unsigned char c = (unsigned char)data[i];
		if (testInCharField(cf, c)) {
			return i;
		}
	}
	return BSTR_ERR;
}

int
binchr(const bstring b0, int pos, const bstring b1)
{
	struct charField chrs;
	if (pos < 0 || b0 == NULL ||
	    b0->data == NULL || b0->slen <= pos) {
		return BSTR_ERR;
	}
	if (1 == b1->slen) {
		return bstrchrp(b0, b1->data[0], pos);
	}
	if (0 > buildCharField (&chrs, b1)) {
		return BSTR_ERR;
	}
	return binchrCF(b0->data, b0->slen, pos, &chrs);
}

/* Inner engine for binchrr */
static int
binchrrCF(const unsigned char *data, int pos, const struct charField *cf)
{
	int i;
	for (i = pos; i >= 0; i--) {
		unsigned int c = (unsigned int)data[i];
		if (testInCharField(cf, c)) {
			return i;
		}
	}
	return BSTR_ERR;
}

int
binchrr(const bstring b0, int pos, const bstring b1)
{
	struct charField chrs;
	if (pos < 0 || b0 == NULL ||
	    b0->data == NULL || b1 == NULL ||
	    b0->slen < pos) {
		return BSTR_ERR;
	}
	if (pos == b0->slen) {
		pos--;
	}
	if (1 == b1->slen) {
		return bstrrchrp(b0, b1->data[0], pos);
	}
	if (0 > buildCharField(&chrs, b1)) {
		return BSTR_ERR;
	}
	return binchrrCF(b0->data, pos, &chrs);
}

int
bninchr(const bstring b0, int pos, const bstring b1)
{
	struct charField chrs;
	if (pos < 0 || b0 == NULL ||
	    b0->data == NULL || b0->slen <= pos) {
		return BSTR_ERR;
	}
	if (buildCharField(&chrs, b1) < 0) {
		return BSTR_ERR;
	}
	invertCharField(&chrs);
	return binchrCF(b0->data, b0->slen, pos, &chrs);
}

int
bninchrr(const bstring b0, int pos, const bstring b1)
{
	struct charField chrs;
	if (pos < 0 || b0 == NULL ||
	    b0->data == NULL || b0->slen < pos) {
		return BSTR_ERR;
	}
	if (pos == b0->slen) {
		pos--;
	}
	if (buildCharField(&chrs, b1) < 0) {
		return BSTR_ERR;
	}
	invertCharField(&chrs);
	return binchrrCF(b0->data, pos, &chrs);
}

int
bsetstr(bstring b0, int pos, const bstring b1, unsigned char fill)
{
	int d, newlen;
	ptrdiff_t pd;
	bstring aux = (bstring) b1;
	if (pos < 0 || b0 == NULL || b0->slen < 0 ||
	    NULL == b0->data || b0->mlen < b0->slen || b0->mlen <= 0) {
		return BSTR_ERR;
	}
	if (b1 != NULL && (b1->slen < 0 || b1->data == NULL)) {
		return BSTR_ERR;
	}
	d = pos;
	/* Aliasing case */
	if (NULL != aux) {
		if ((pd = (ptrdiff_t)(b1->data - b0->data)) >= 0 &&
		    pd < (ptrdiff_t) b0->mlen) {
			if (NULL == (aux = bstrcpy (b1))) return BSTR_ERR;
		}
		d += aux->slen;
	}
	/* Increase memory size if necessary */
	if (balloc(b0, d + 1) != BSTR_OK) {
		if (aux != b1) {
			bdestroy (aux);
		}
		return BSTR_ERR;
	}
	newlen = b0->slen;
	/* Fill in "fill" character as necessary */
	if (pos > newlen) {
		memset(b0->data + b0->slen, (int)fill,
		       (size_t)(pos - b0->slen));
		newlen = pos;
	}
	/* Copy b1 to position pos in b0. */
	if (aux != NULL) {
		bBlockCopy((char *)(b0->data + pos), (char *)aux->data,
			   aux->slen);
		if (aux != b1) {
			bdestroy(aux);
		}
	}
	/* Indicate the potentially increased size of b0 */
	if (d > newlen) {
		newlen = d;
	}
	b0->slen = newlen;
	b0->data[newlen] = (unsigned char)'\0';
	return BSTR_OK;
}

int
binsert(bstring b1, int pos, const bstring b2, unsigned char fill)
{
	int d, l;
	ptrdiff_t pd;
	bstring aux = (bstring) b2;
	if (pos < 0 || b1 == NULL || b2 == NULL || b1->slen < 0 ||
	    b2->slen < 0 || b1->mlen < b1->slen || b1->mlen <= 0) {
		return BSTR_ERR;
	}
	/* Aliasing case */
	if ((pd = (ptrdiff_t) (b2->data - b1->data)) >= 0 &&
	    pd < (ptrdiff_t) b1->mlen) {
		if (NULL == (aux = bstrcpy (b2))) {
			return BSTR_ERR;
		}
	}
	/* Compute the two possible end pointers */
	d = b1->slen + aux->slen;
	l = pos + aux->slen;
	if ((d|l) < 0) {
		return BSTR_ERR;
	}
	if (l > d) {
		/* Inserting past the end of the string */
		if (balloc(b1, l + 1) != BSTR_OK) {
			if (aux != b2) {
				bdestroy(aux);
			}
			return BSTR_ERR;
		}
		memset(b1->data + b1->slen, (int)fill,
		       (size_t)(pos - b1->slen));
		b1->slen = l;
	} else {
		/* Inserting in the middle of the string */
		if (balloc(b1, d + 1) != BSTR_OK) {
			if (aux != b2) {
				bdestroy(aux);
			}
			return BSTR_ERR;
		}
		bBlockCopy(b1->data + l, b1->data + pos, d - l);
		b1->slen = d;
	}
	bBlockCopy(b1->data + pos, aux->data, aux->slen);
	b1->data[b1->slen] = (unsigned char)'\0';
	if (aux != b2) {
		bdestroy(aux);
	}
	return BSTR_OK;
}

int
breplace(bstring b1, int pos, int len, const bstring b2, unsigned char fill)
{
	int pl, ret;
	ptrdiff_t pd;
	bstring aux = (bstring) b2;
	if (pos < 0 || len < 0 || (pl = pos + len) < 0 || b1 == NULL ||
	    b2 == NULL || b1->data == NULL || b2->data == NULL ||
	    b1->slen < 0 || b2->slen < 0 || b1->mlen < b1->slen ||
	    b1->mlen <= 0) {
		return BSTR_ERR;
	}
	/* Straddles the end? */
	if (pl >= b1->slen) {
		if ((ret = bsetstr (b1, pos, b2, fill)) < 0) {
			return ret;
		}
		if (pos + b2->slen < b1->slen) {
			b1->slen = pos + b2->slen;
			b1->data[b1->slen] = (unsigned char) '\0';
		}
		return ret;
	}
	/* Aliasing case */
	pd = (ptrdiff_t)(b2->data - b1->data);
	if (pd >= 0 && pd < (ptrdiff_t)b1->slen) {
		aux = bstrcpy(b2);
		if (!aux) {
			return BSTR_ERR;
		}
	}
	if (aux->slen > len) {
		if (balloc(b1, b1->slen + aux->slen - len) != BSTR_OK) {
			if (aux != b2) {
				bdestroy(aux);
			}
			return BSTR_ERR;
		}
	}
	if (aux->slen != len) {
		memmove(b1->data + pos + aux->slen, b1->data + pos + len,
			b1->slen - (pos + len));
	}
	memcpy(b1->data + pos, aux->data, aux->slen);
	b1->slen += aux->slen - len;
	b1->data[b1->slen] = (unsigned char)'\0';
	if (aux != b2) {
		bdestroy(aux);
	}
	return BSTR_OK;
}

typedef int (*instr_fnptr)(const bstring s1, int pos, const bstring s2);

#define INITIAL_STATIC_FIND_INDEX_COUNT 32

/*
 *  findreplaceengine is used to implement bfindreplace and
 *  bfindreplacecaseless. It works by breaking the three cases of
 *  expansion, reduction and replacement, and solving each of these
 *  in the most efficient way possible.
 */
static int
findreplaceengine(bstring b, const bstring find, const bstring repl,
		  int pos, instr_fnptr instr)
{
	int i, ret, slen, mlen, delta, acc;
	int *d;
	/* This +1 is unnecessary, but it shuts up LINT. */
	int static_d[INITIAL_STATIC_FIND_INDEX_COUNT + 1];
	ptrdiff_t pd;
	bstring auxf = (bstring) find;
	bstring auxr = (bstring) repl;
	if (!b || !b->data || !find ||
	    !find->data || !repl || !repl->data ||
	    pos < 0 || find->slen <= 0 || b->mlen < 0 ||
	    b->slen > b->mlen || b->mlen <= 0 || b->slen < 0 ||
	    repl->slen < 0) {
		return BSTR_ERR;
	}
	if (pos > b->slen - find->slen) {
		return BSTR_OK;
	}
	/* Alias with find string */
	pd = (ptrdiff_t)(find->data - b->data);
	if ((ptrdiff_t)(pos - find->slen) < pd && pd < (ptrdiff_t)b->slen) {
		auxf = bstrcpy(find);
		if (!auxf) {
			return BSTR_ERR;
		}
	}
	/* Alias with repl string */
	pd = (ptrdiff_t)(repl->data - b->data);
	if ((ptrdiff_t)(pos - repl->slen) < pd && pd < (ptrdiff_t)b->slen) {
		auxr = bstrcpy (repl);
		if (!auxr) {
			if (auxf != find) {
				bdestroy(auxf);
			}
			return BSTR_ERR;
		}
	}
	delta = auxf->slen - auxr->slen;
	/* in-place replacement since find and replace strings are of equal
	 * length
	 */
	if (delta == 0) {
		while ((pos = instr(b, pos, auxf)) >= 0) {
			memcpy(b->data + pos, auxr->data, auxr->slen);
			pos += auxf->slen;
		}
		if (auxf != find) {
			bdestroy (auxf);
		}
		if (auxr != repl) {
			bdestroy (auxr);
		}
		return BSTR_OK;
	}
	/* shrinking replacement since auxf->slen > auxr->slen */
	if (delta > 0) {
		acc = 0;
		while ((i = instr (b, pos, auxf)) >= 0) {
			if (acc && i > pos) {
				memmove(b->data + pos - acc, b->data + pos,
					i - pos);
			}
			if (auxr->slen) {
				memcpy(b->data + i - acc, auxr->data,
				       auxr->slen);
			}
			acc += delta;
			pos = i + auxf->slen;
		}

		if (acc) {
			i = b->slen;
			if (i > pos) {
				memmove(b->data + pos - acc, b->data + pos,
					i - pos);
			}
			b->slen -= acc;
			b->data[b->slen] = (unsigned char) '\0';
		}

		if (auxf != find) {
			bdestroy (auxf);
		}
		if (auxr != repl) {
			bdestroy (auxr);
		}
		return BSTR_OK;
	}
	/* expanding replacement since find->slen < repl->slen. Its a lot
	 * more complicated. This works by first finding all the matches and
	 * storing them to a growable array, then doing at most one resize of
	 * the destination bstring and then performing the direct memory
	 * transfers of the string segment pieces to form the final result. The
	 * growable array of matches uses a deferred doubling reallocing
	 * strategy. What this means is that it starts as a reasonably fixed
	 * sized auto array in the hopes that many if not most cases will never
	 * need to grow this array. But it switches as soon as the bounds of
	 * the array will be exceeded. An extra find result is always appended
	 * to this array that corresponds to the end of the destination string,
	 * so slen is checked against mlen - 1 rather than mlen before
	 * resizing.
	*/
	mlen = INITIAL_STATIC_FIND_INDEX_COUNT;
	d = (int *) static_d; /* Avoid malloc for trivial/initial cases */
	acc = slen = 0;
	while ((pos = instr(b, pos, auxf)) >= 0) {
		if (slen >= mlen - 1) {
			int sl, *t;
			mlen += mlen;
			sl = sizeof(int *) * mlen;
			if (static_d == d) {
				/* static_d cannot be realloced */
				d = NULL;
			}
			if (mlen <= 0 || sl < mlen ||
			    NULL == (t = (int *) realloc(d, sl))) {
				ret = BSTR_ERR;
				goto done;
			}
			if (NULL == d) {
				memcpy(t, static_d, sizeof (static_d));
			}
			d = t;
		}
		d[slen] = pos;
		slen++;
		acc -= delta;
		pos += auxf->slen;
		if (pos < 0 || acc < 0) {
			ret = BSTR_ERR;
			goto done;
		}
	}
	/* slen <= INITIAL_STATIC_INDEX_COUNT-1 or mlen-1 here. */
	d[slen] = b->slen;
	ret = balloc (b, b->slen + acc + 1);
	if (BSTR_OK == ret) {
		b->slen += acc;
		for (i = slen-1; i >= 0; i--) {
			int s, l;
			s = d[i] + auxf->slen;
			l = d[i+1] - s; /* d[slen] may be accessed here. */
			if (l) {
				memmove(b->data + s + acc, b->data + s, l);
			}
			if (auxr->slen) {
				memmove(b->data + s + acc - auxr->slen,
					auxr->data, auxr->slen);
			}
			acc += delta;
		}
		b->data[b->slen] = (unsigned char)'\0';
	}
done:
	if (static_d == d) {
		d = NULL;
	}
	free(d);
	if (auxf != find) {
		bdestroy(auxf);
	}
	if (auxr != repl) {
		bdestroy(auxr);
	}
	return ret;
}

int
bfindreplace(bstring b, const bstring find, const bstring repl, int pos)
{
	return findreplaceengine(b, find, repl, pos, binstr);
}

int
bfindreplacecaseless(bstring b, const bstring find, const bstring repl, int pos)
{
	return findreplaceengine(b, find, repl, pos, binstrcaseless);
}

int
binsertch(bstring b, int pos, int len, unsigned char fill)
{
	int d, l, i;
	if (pos < 0 || !b ||
	    b->slen < 0 || b->mlen < b->slen ||
	    b->mlen <= 0 || len < 0) {
		return BSTR_ERR;
	}
	/* Compute the two possible end pointers */
	d = b->slen + len;
	l = pos + len;
	if ((d|l) < 0) {
		return BSTR_ERR;
	}
	if (l > d) {
		/* Inserting past the end of the string */
		if (balloc(b, l + 1) != BSTR_OK) {
			return BSTR_ERR;
		}
		pos = b->slen;
		b->slen = l;
	} else {
		/* Inserting in the middle of the string */
		if (balloc(b, d + 1) != BSTR_OK) {
			return BSTR_ERR;
		}
		for (i = d - 1; i >= l; i--) {
			b->data[i] = b->data[i - len];
		}
		b->slen = d;
	}
	for (i = pos; i < l; i++) {
		b->data[i] = fill;
	}
	b->data[b->slen] = (unsigned char)'\0';
	return BSTR_OK;
}

int
bpattern(bstring b, int len)
{
	int i, d;
	d = blength(b);
	if (d <= 0 || len < 0 || balloc(b, len + 1) != BSTR_OK) {
		return BSTR_ERR;
	}
	if (len > 0) {
		if (d == 1) {
			return bsetstr(b, len, NULL, b->data[0]);
		}
		for (i = d; i < len; i++) {
			b->data[i] = b->data[i - d];
		}
	}
	b->data[len] = (unsigned char)'\0';
	b->slen = len;
	return BSTR_OK;
}

#define BS_BUFF_SZ (1024)

int
breada(bstring b, bNread readPtr, void *parm)
{
	int i, l, n;
	if (b == NULL || b->mlen <= 0 ||
	    b->slen < 0 || b->mlen < b->slen ||
	    b->mlen <= 0 || readPtr == NULL) {
		return BSTR_ERR;
	}
	i = b->slen;
	for (n = i + 16; ; n += ((n < BS_BUFF_SZ) ? n : BS_BUFF_SZ)) {
		if (BSTR_OK != balloc(b, n + 1)) {
			return BSTR_ERR;
		}
		l = (int)readPtr((void *)(b->data + i), 1, n - i, parm);
		i += l;
		b->slen = i;
		if (i < n) {
			break;
		}
	}
	b->data[i] = (unsigned char)'\0';
	return BSTR_OK;
}

bstring
bread(bNread readPtr, void *parm)
{
	bstring buff;
	if (0 > breada(buff = bfromcstr (""), readPtr, parm)) {
		bdestroy(buff);
		return NULL;
	}
	return buff;
}

int
bassigngets(bstring b, bNgetc getcPtr, void *parm, char terminator)
{
	int c, d, e;
	if (!b || b->mlen <= 0 ||
	    b->slen < 0 || b->mlen < b->slen ||
	    b->mlen <= 0 || getcPtr == NULL) {
		return BSTR_ERR;
	}
	d = 0;
	e = b->mlen - 2;
	while ((c = getcPtr(parm)) >= 0) {
		if (d > e) {
			b->slen = d;
			if (balloc (b, d + 2) != BSTR_OK) {
				return BSTR_ERR;
			}
			e = b->mlen - 2;
		}
		b->data[d] = (unsigned char)c;
		d++;
		if (c == terminator) {
			break;
		}
	}
	b->data[d] = (unsigned char)'\0';
	b->slen = d;
	return d == 0 && c < 0;
}

int
bgetsa(bstring b, bNgetc getcPtr, void *parm, char terminator)
{
	int c, d, e;
	if (!b || b->mlen <= 0 ||
	    b->slen < 0 || b->mlen < b->slen ||
	    b->mlen <= 0 || !getcPtr) {
		return BSTR_ERR;
	}
	d = b->slen;
	e = b->mlen - 2;
	while ((c = getcPtr(parm)) >= 0) {
		if (d > e) {
			b->slen = d;
			if (balloc(b, d + 2) != BSTR_OK) {
				return BSTR_ERR;
			}
			e = b->mlen - 2;
		}
		b->data[d] = (unsigned char) c;
		d++;
		if (c == terminator) {
			break;
		}
	}
	b->data[d] = (unsigned char)'\0';
	b->slen = d;
	return d == 0 && c < 0;
}

bstring
bgets(bNgetc getcPtr, void *parm, char terminator)
{
	bstring buff;
	if (0 > bgetsa(buff = bfromcstr (""), getcPtr, parm, terminator) ||
	    0 >= buff->slen) {
		bdestroy(buff);
		buff = NULL;
	}
	return buff;
}

struct bStream {
	bstring buff; /* Buffer for over-reads */
	void *parm; /* The stream handle for core stream */
	bNread readFnPtr; /* fread compatible fnptr for core stream */
	int isEOF; /* track file's EOF state */
	int maxBuffSz;
};

struct bStream *
bsopen (bNread readPtr, void *parm)
{
	struct bStream *s;
	if (readPtr == NULL) {
		return NULL;
	}
	s = malloc(sizeof (struct bStream));
	if (!s) {
		return NULL;
	}
	s->parm = parm;
	s->buff = bfromcstr ("");
	s->readFnPtr = readPtr;
	s->maxBuffSz = BS_BUFF_SZ;
	s->isEOF = 0;
	return s;
}

int bsbufflength(struct bStream *s, int sz)
{
	int oldSz;
	if (!s || sz < 0) {
		return BSTR_ERR;
	}
	oldSz = s->maxBuffSz;
	if (sz > 0) {
		s->maxBuffSz = sz;
	}
	return oldSz;
}

int
bseof(const struct bStream *s)
{
	if (!s || !s->readFnPtr) {
		return BSTR_ERR;
	}
	return s->isEOF && (s->buff->slen == 0);
}

void *
bsclose(struct bStream *s)
{
	void *parm;
	if (s == NULL) {
		return NULL;
	}
	s->readFnPtr = NULL;
	if (s->buff) {
		bdestroy(s->buff);
	}
	s->buff = NULL;
	parm = s->parm;
	s->parm = NULL;
	s->isEOF = 1;
	free(s);
	return parm;
}

int
bsreadlna(bstring r, struct bStream *s, char terminator)
{
	int i, l, ret, rlo;
	char *b;
	struct tagbstring x;
	if (!s || !s->buff ||
	    !r || r->mlen <= 0 ||
	    r->slen < 0 || r->mlen < r->slen) {
		return BSTR_ERR;
	}
	l = s->buff->slen;
	if (BSTR_OK != balloc(s->buff, s->maxBuffSz + 1)) {
		return BSTR_ERR;
	}
	b = (char *)s->buff->data;
	x.data = (unsigned char *)b;
	/* First check if the current buffer holds the terminator */
	b[l] = terminator; /* Set sentinel */
	for (i=0; b[i] != terminator; i++) ;
	if (i < l) {
		x.slen = i + 1;
		ret = bconcat(r, &x);
		s->buff->slen = l;
		if (BSTR_OK == ret) {
			bdelete(s->buff, 0, i + 1);
		}
		return BSTR_OK;
	}
	rlo = r->slen;
	/* If not then just concatenate the entire buffer to the output */
	x.slen = l;
	if (BSTR_OK != bconcat(r, &x)) {
		return BSTR_ERR;
	}
	/* Perform direct in-place reads into the destination to allow for
	 * the minimum of data-copies
	 */
	while (1) {
		if (BSTR_OK != balloc(r, r->slen + s->maxBuffSz + 1)) {
			return BSTR_ERR;
		}
		b = (char *) (r->data + r->slen);
		l = (int) s->readFnPtr(b, 1, s->maxBuffSz, s->parm);
		if (l <= 0) {
			r->data[r->slen] = (unsigned char)'\0';
			s->buff->slen = 0;
			s->isEOF = 1;
			/* If nothing was read return with an error message */
			return BSTR_ERR & -(r->slen == rlo);
		}
		b[l] = terminator; /* Set sentinel */
		for (i=0; b[i] != terminator; i++)
			;
		if (i < l) {
			break;
		}
		r->slen += l;
	}
	/* Terminator found, push over-read back to buffer */
	i++;
	r->slen += i;
	s->buff->slen = l - i;
	memcpy(s->buff->data, b + i, l - i);
	r->data[r->slen] = (unsigned char)'\0';
	return BSTR_OK;
}

int
bsreadlnsa(bstring r, struct bStream *s, const bstring term)
{
	int i, l, ret, rlo;
	unsigned char *b;
	struct tagbstring x;
	struct charField cf;
	if (!s || !s->buff || !r || !term ||
	    !term->data || r->mlen <= 0 || r->slen < 0 ||
	    r->mlen < r->slen) {
		return BSTR_ERR;
	}
	if (term->slen == 1) {
		return bsreadlna(r, s, term->data[0]);
	}
	if (term->slen < 1 || buildCharField(&cf, term)) {
		return BSTR_ERR;
	}
	l = s->buff->slen;
	if (BSTR_OK != balloc(s->buff, s->maxBuffSz + 1)) {
		return BSTR_ERR;
	}
	b = (unsigned char *)s->buff->data;
	x.data = b;
	/* First check if the current buffer holds the terminator */
	b[l] = term->data[0]; /* Set sentinel */
	for (i = 0; !testInCharField(&cf, b[i]); i++)
		;
	if (i < l) {
		x.slen = i + 1;
		ret = bconcat(r, &x);
		s->buff->slen = l;
		if (BSTR_OK == ret) {
			bdelete(s->buff, 0, i + 1);
		}
		return BSTR_OK;
	}
	rlo = r->slen;
	/* If not then just concatenate the entire buffer to the output */
	x.slen = l;
	if (BSTR_OK != bconcat(r, &x)) {
		return BSTR_ERR;
	}
	/* Perform direct in-place reads into the destination to allow for
	 * the minimum of data-copies
	 */
	while (1) {
		if (BSTR_OK != balloc(r, r->slen + s->maxBuffSz + 1)) {
			return BSTR_ERR;
		}
		b = (unsigned char *)(r->data + r->slen);
		l = (int) s->readFnPtr(b, 1, s->maxBuffSz, s->parm);
		if (l <= 0) {
			r->data[r->slen] = (unsigned char)'\0';
			s->buff->slen = 0;
			s->isEOF = 1;
			/* If nothing was read return with an error message */
			return BSTR_ERR & -(r->slen == rlo);
		}
		b[l] = term->data[0]; /* Set sentinel */
		for (i = 0; !testInCharField(&cf, b[i]); i++)
			;
		if (i < l) {
			break;
		}
		r->slen += l;
	}
	/* Terminator found, push over-read back to buffer */
	i++;
	r->slen += i;
	s->buff->slen = l - i;
	memcpy(s->buff->data, b + i, l - i);
	r->data[r->slen] = (unsigned char)'\0';
	return BSTR_OK;
}

int
bsreada(bstring r, struct bStream *s, int n)
{
	int l, ret, orslen;
	char *b;
	struct tagbstring x;
	if (!s || !s->buff || !r || r->mlen <= 0
	 || r->slen < 0 || r->mlen < r->slen || n <= 0) {
		return BSTR_ERR;
	}
	n += r->slen;
	if (n <= 0) {
		return BSTR_ERR;
	}
	l = s->buff->slen;
	orslen = r->slen;
	if (0 == l) {
		if (s->isEOF) {
			return BSTR_ERR;
		}
		if (r->mlen > n) {
			l = (int)s->readFnPtr(r->data + r->slen, 1,
					      n - r->slen, s->parm);
			if (0 >= l || l > n - r->slen) {
				s->isEOF = 1;
				return BSTR_ERR;
			}
			r->slen += l;
			r->data[r->slen] = (unsigned char)'\0';
			return 0;
		}
	}
	if (BSTR_OK != balloc(s->buff, s->maxBuffSz + 1)) {
		return BSTR_ERR;
	}
	b = (char *) s->buff->data;
	x.data = (unsigned char *)b;
	do {
		if (l + r->slen >= n) {
			x.slen = n - r->slen;
			ret = bconcat(r, &x);
			s->buff->slen = l;
			if (BSTR_OK == ret) {
				bdelete(s->buff, 0, x.slen);
			}
			return BSTR_ERR & -(r->slen == orslen);
		}
		x.slen = l;
		if (BSTR_OK != bconcat (r, &x)) {
			break;
		}
		l = n - r->slen;
		if (l > s->maxBuffSz) {
			l = s->maxBuffSz;
		}
		l = (int)s->readFnPtr(b, 1, l, s->parm);

	} while (l > 0);
	if (l < 0) {
		l = 0;
	}
	if (l == 0) {
		s->isEOF = 1;
	}
	s->buff->slen = l;
	return BSTR_ERR & -(r->slen == orslen);
}

int
bsreadln(bstring r, struct bStream *s, char terminator)
{
	if (!s || !s->buff || !r || r->mlen <= 0) {
		return BSTR_ERR;
	}
	if (BSTR_OK != balloc(s->buff, s->maxBuffSz + 1)) {
		return BSTR_ERR;
	}
	r->slen = 0;
	return bsreadlna(r, s, terminator);
}

int
bsreadlns(bstring r, struct bStream *s, const bstring term)
{
	if (!s || !s->buff || !r || !term || !term->data || r->mlen <= 0) {
		return BSTR_ERR;
	}
	if (term->slen == 1) {
		return bsreadln (r, s, term->data[0]);
	}
	if (term->slen < 1) {
		return BSTR_ERR;
	}
	if (BSTR_OK != balloc(s->buff, s->maxBuffSz + 1)) {
		return BSTR_ERR;
	}
	r->slen = 0;
	return bsreadlnsa(r, s, term);
}

int
bsread(bstring r, struct bStream *s, int n)
{
	if (!s || !s->buff || !r || r->mlen <= 0 || n <= 0) {
		return BSTR_ERR;
	}
	if (BSTR_OK != balloc(s->buff, s->maxBuffSz + 1)) {
		return BSTR_ERR;
	}
	r->slen = 0;
	return bsreada(r, s, n);
}

int
bsunread(struct bStream *s, const bstring b)
{
	if (!s || !s->buff) {
		return BSTR_ERR;
	}
	return binsert(s->buff, 0, b, (unsigned char)'?');
}

int
bspeek(bstring r, const struct bStream *s)
{
	if (!s || !s->buff) {
		return BSTR_ERR;
	}
	return bassign(r, s->buff);
}

bstring
bjoin(const struct bstrList *bl, const bstring sep)
{
	bstring b;
	int i, c, v;
	if (bl == NULL || bl->qty < 0) {
		return NULL;
	}
	if (sep != NULL && (sep->slen < 0 || sep->data == NULL)) {
		return NULL;
	}
	for (i = 0, c = 1; i < bl->qty; i++) {
		v = bl->entry[i]->slen;
		if (v < 0) {
			return NULL; /* Invalid input */
		}
		c += v;
		if (c < 0) {
			return NULL; /* Wrap around ?? */
		}
	}
	if (sep != NULL) {
		c += (bl->qty - 1) * sep->slen;
	}
	b = (bstring)malloc(sizeof(struct tagbstring));
	if (NULL == b) {
		return NULL; /* Out of memory */
	}
	b->data = (unsigned char *)malloc(c);
	if (b->data == NULL) {
		free (b);
		return NULL;
	}
	b->mlen = c;
	b->slen = c-1;
	for (i = 0, c = 0; i < bl->qty; i++) {
		if (i > 0 && sep != NULL) {
			memcpy(b->data + c, sep->data, sep->slen);
			c += sep->slen;
		}
		v = bl->entry[i]->slen;
		memcpy(b->data + c, bl->entry[i]->data, v);
		c += v;
	}
	b->data[c] = (unsigned char)'\0';
	return b;
}

#define BSSSC_BUFF_LEN (256)

int
bssplitscb(struct bStream *s, const bstring splitStr,
	   int (*cb)(void *parm, int ofs, const bstring entry),
	   void *parm)
{
	struct charField chrs;
	bstring buff;
	int i, p, ret;
	if (!cb || !s || !s->readFnPtr ||
	    !splitStr || splitStr->slen < 0) {
		return BSTR_ERR;
	}
	buff = bfromcstr ("");
	if (!buff) {
		return BSTR_ERR;
	}
	if (splitStr->slen == 0) {
		while (bsreada(buff, s, BSSSC_BUFF_LEN) >= 0)
			;
		if ((ret = cb(parm, 0, buff)) > 0) {
			ret = 0;
		}
	} else {
		buildCharField(&chrs, splitStr);
		ret = p = i = 0;
		while (1) {
			if (i >= buff->slen) {
				bsreada(buff, s, BSSSC_BUFF_LEN);
				if (i >= buff->slen) {
					if (0 < (ret = cb (parm, p, buff))) {
						ret = 0;
					}
					break;
				}
			}
			if (testInCharField(&chrs, buff->data[i])) {
				struct tagbstring t;
				unsigned char c;
				blk2tbstr(t, buff->data + i + 1,
					  buff->slen - (i + 1));
				if ((ret = bsunread(s, &t)) < 0) {
					break;
				}
				buff->slen = i;
				c = buff->data[i];
				buff->data[i] = (unsigned char)'\0';
				if ((ret = cb(parm, p, buff)) < 0) {
					break;
				}
				buff->data[i] = c;
				buff->slen = 0;
				p += i + 1;
				i = -1;
			}
			i++;
		}
	}
	bdestroy(buff);
	return ret;
}

int
bssplitstrcb(struct bStream *s, const bstring splitStr,
	     int (*cb)(void *parm, int ofs, const bstring entry),
	     void *parm)
{
	bstring buff;
	int i, p, ret;
	if (!cb || !s || !s->readFnPtr ||
	    !splitStr || splitStr->slen < 0) {
		return BSTR_ERR;
	}
	if (splitStr->slen == 1) {
		return bssplitscb(s, splitStr, cb, parm);
	}
	buff = bfromcstr ("");
	if (!buff) {
		return BSTR_ERR;
	}
	if (splitStr->slen == 0) {
		for (i = 0; bsreada(buff, s, BSSSC_BUFF_LEN) >= 0; i++) {
			if ((ret = cb (parm, 0, buff)) < 0) {
				bdestroy(buff);
				return ret;
			}
			buff->slen = 0;
		}
		return BSTR_OK;
	} else {
		ret = p = i = 0;
		for (i = p = 0; ;) {
			ret = binstr (buff, 0, splitStr);
			if (ret >= 0) {
				struct tagbstring t;
				blk2tbstr(t, buff->data, ret);
				i = ret + splitStr->slen;
				ret = cb (parm, p, &t);
				if (ret < 0) {
					break;
				}
				p += i;
				bdelete(buff, 0, i);
			} else {
				bsreada(buff, s, BSSSC_BUFF_LEN);
				if (bseof(s)) {
					ret = cb (parm, p, buff);
					if (ret > 0) {
						ret = 0;
					}
					break;
				}
			}
		}
	}
	bdestroy(buff);
	return ret;
}

struct bstrList *
bstrListCreate(void)
{
	struct bstrList *sl = malloc(sizeof(struct bstrList));
	if (sl) {
		sl->entry = (bstring *)malloc(1 * sizeof(bstring));
		if (!sl->entry) {
			free(sl);
			sl = NULL;
		} else {
			sl->qty = 0;
			sl->mlen = 1;
		}
	}
	return sl;
}

int
bstrListDestroy(struct bstrList *sl)
{
	int i;
	if (!sl || sl->qty < 0) {
		return BSTR_ERR;
	}
	for (i = 0; i < sl->qty; i++) {
		if (sl->entry[i]) {
			bdestroy(sl->entry[i]);
			sl->entry[i] = NULL;
		}
	}
	sl->qty  = -1;
	sl->mlen = -1;
	free(sl->entry);
	sl->entry = NULL;
	free(sl);
	return BSTR_OK;
}

int
bstrListAlloc(struct bstrList *sl, int msz)
{
	bstring *l;
	int smsz;
	size_t nsz;
	if (!sl || msz <= 0 ||
	    !sl->entry || sl->qty < 0 ||
	    sl->mlen <= 0 || sl->qty > sl->mlen) {
		return BSTR_ERR;
	}
	if (sl->mlen >= msz) {
		return BSTR_OK;
	}
	smsz = snapUpSize(msz);
	nsz = ((size_t)smsz) * sizeof(bstring);
	if (nsz < (size_t) smsz) {
		return BSTR_ERR;
	}
	l = realloc(sl->entry, nsz);
	if (!l) {
		smsz = msz;
		nsz = ((size_t)smsz) * sizeof(bstring);
		l = realloc(sl->entry, nsz);
		if (!l) {
			return BSTR_ERR;
		}
	}
	sl->mlen = smsz;
	sl->entry = l;
	return BSTR_OK;
}

int
bstrListAllocMin(struct bstrList *sl, int msz)
{
	bstring *l;
	size_t nsz;
	if (!sl || msz <= 0 ||
	    !sl->entry || sl->qty < 0 ||
	    sl->mlen <= 0 || sl->qty > sl->mlen) {
		return BSTR_ERR;
	}
	if (msz < sl->qty) {
		msz = sl->qty;
	}
	if (sl->mlen == msz) {
		return BSTR_OK;
	}
	nsz = ((size_t)msz) * sizeof(bstring);
	if (nsz < (size_t)msz) {
		return BSTR_ERR;
	}
	l = realloc(sl->entry, nsz);
	if (!l) {
		return BSTR_ERR;
	}
	sl->mlen = msz;
	sl->entry = l;
	return BSTR_OK;
}

int
bsplitcb(const bstring str, unsigned char splitChar, int pos,
	 int (*cb) (void *parm, int ofs, int len),
	 void *parm)
{
	int i, p, ret;
	if (!cb || !str || pos < 0 || pos > str->slen) {
		return BSTR_ERR;
	}
	p = pos;
	do {
		for (i = p; i < str->slen; i++) {
			if (str->data[i] == splitChar) {
				break;
			}
		}
		if ((ret = cb(parm, p, i - p)) < 0) {
			return ret;
		}
		p = i + 1;
	} while (p <= str->slen);
	return BSTR_OK;
}

int
bsplitscb(const bstring str, const bstring splitStr, int pos,
	  int (*cb)(void *parm, int ofs, int len),
	  void *parm)
{
	struct charField chrs;
	int i, p, ret;
	if (!cb || !str || pos < 0 || pos > str->slen ||
	    !splitStr || splitStr->slen < 0) {
		return BSTR_ERR;
	}
	if (splitStr->slen == 0) {
		if ((ret = cb (parm, 0, str->slen)) > 0) {
			ret = 0;
		}
		return ret;
	}
	if (splitStr->slen == 1) {
		return bsplitcb (str, splitStr->data[0], pos, cb, parm);
	}
	buildCharField(&chrs, splitStr);
	p = pos;
	do {
		for (i = p; i < str->slen; i++) {
			if (testInCharField(&chrs, str->data[i])) {
				break;
			}
		}
		if ((ret = cb(parm, p, i - p)) < 0) {
			return ret;
		}
		p = i + 1;
	} while (p <= str->slen);
	return BSTR_OK;
}

int
bsplitstrcb(const bstring str, const bstring splitStr, int pos,
	    int (*cb)(void *parm, int ofs, int len),
	    void *parm)
{
	int i, p, ret;
	if (!cb || !str || pos < 0 || pos > str->slen ||
	    !splitStr || splitStr->slen < 0) {
		return BSTR_ERR;
	}
	if (0 == splitStr->slen) {
		for (i = pos; i < str->slen; i++) {
			ret = cb (parm, i, 1);
			if (ret < 0) {
				return ret;
			}
		}
		return BSTR_OK;
	}
	if (splitStr->slen == 1) {
		return bsplitcb(str, splitStr->data[0], pos, cb, parm);
	}
	for (i = p = pos; i <= str->slen - splitStr->slen; i++) {
		ret = memcmp(splitStr->data, str->data + i, splitStr->slen);
		if (0 == ret) {
			ret = cb (parm, p, i - p);
			if (ret < 0) {
				return ret;
			}
			i += splitStr->slen;
			p = i;
		}
	}
	ret = cb (parm, p, str->slen - p);
	if (ret < 0) {
		return ret;
	}
	return BSTR_OK;
}

struct genBstrList {
	bstring b;
	struct bstrList *bl;
};

static int
bscb(void *parm, int ofs, int len)
{
	struct genBstrList *g = (struct genBstrList *)parm;
	if (g->bl->qty >= g->bl->mlen) {
		int mlen = g->bl->mlen * 2;
		bstring *tbl;
		while (g->bl->qty >= mlen) {
			if (mlen < g->bl->mlen) {
				return BSTR_ERR;
			}
			mlen += mlen;
		}
		tbl = (bstring *)realloc(g->bl->entry, sizeof(bstring) * mlen);
		if (tbl == NULL) {
			return BSTR_ERR;
		}
		g->bl->entry = tbl;
		g->bl->mlen = mlen;
	}
	g->bl->entry[g->bl->qty] = bmidstr(g->b, ofs, len);
	g->bl->qty++;
	return BSTR_OK;
}

struct bstrList *
bsplit(const bstring str, unsigned char splitChar)
{
	struct genBstrList g;
	if (!str || !str->data || str->slen < 0) {
		return NULL;
	}
	g.bl = malloc(sizeof(struct bstrList));
	if (!g.bl) {
		return NULL;
	}
	g.bl->mlen = 4;
	g.bl->entry = malloc(g.bl->mlen * sizeof(bstring));
	if (!g.bl->entry) {
		free(g.bl);
		return NULL;
	}

	g.b = (bstring)str;
	g.bl->qty = 0;
	if (bsplitcb(str, splitChar, 0, bscb, &g) < 0) {
		bstrListDestroy(g.bl);
		return NULL;
	}
	return g.bl;
}

struct bstrList *
bsplitstr(const bstring str, const bstring splitStr)
{
	struct genBstrList g;
	if (!str || !str->data || str->slen < 0) {
		return NULL;
	}
	g.bl = malloc(sizeof(struct bstrList));
	if (!g.bl) {
		return NULL;
	}
	g.bl->mlen = 4;
	g.bl->entry = malloc(g.bl->mlen * sizeof (bstring));
	if (!g.bl->entry) {
		free(g.bl);
		return NULL;
	}
	g.b = (bstring)str;
	g.bl->qty = 0;
	if (bsplitstrcb(str, splitStr, 0, bscb, &g) < 0) {
		bstrListDestroy(g.bl);
		return NULL;
	}
	return g.bl;
}

struct bstrList *
bsplits(const bstring str, const bstring splitStr)
{
	struct genBstrList g;
	if (!str || str->slen < 0 || !str->data ||
	    !splitStr || splitStr->slen < 0 || !splitStr->data) {
		return NULL;
	}
	g.bl = malloc(sizeof(struct bstrList));
	if (!g.bl) {
		return NULL;
	}
	g.bl->mlen = 4;
	g.bl->entry = malloc (g.bl->mlen * sizeof(bstring));
	if (!g.bl->entry) {
		free(g.bl);
		return NULL;
	}
	g.b = (bstring)str;
	g.bl->qty = 0;
	if (bsplitscb(str, splitStr, 0, bscb, &g) < 0) {
		bstrListDestroy(g.bl);
		return NULL;
	}
	return g.bl;
}

#define exvsnprintf(r, b, n, f, a) \
{ \
	r = vsnprintf(b, n, f, a); \
}

#define START_VSNBUFF (16)

/* On IRIX vsnprintf returns n-1 when the operation would overflow the target
 * buffer, WATCOM and MSVC both return -1, while C99 requires that the returned
 * value be exactly what the length would be if the buffer would be large
 * enough.  This leads to the idea that if the return value is larger than n,
 * then changing n to the return value will reduce the number of iterations
 * required.
 */

int
bformata(bstring b, const char *fmt, ...)
{
	va_list arglist;
	bstring buff;
	int n, r;
	if (!b || !fmt || !b->data || b->mlen <= 0 ||
	    b->slen < 0 || b->slen > b->mlen) {
		return BSTR_ERR;
	}
	/* Since the length is not determinable beforehand, a search is
	 * performed using the truncating "vsnprintf" call (to avoid buffer
	 * overflows) on increasing potential sizes for the output result.
	 */
	n = (int)(2 * strlen(fmt));
	if (n < START_VSNBUFF) {
		n = START_VSNBUFF;
	}
	buff = bfromcstralloc(n + 2, "");
	if (!buff) {
		n = 1;
		buff = bfromcstralloc(n + 2, "");
		if (!buff) {
			return BSTR_ERR;
		}
	}
	while (1) {
		va_start(arglist, fmt);
		exvsnprintf(r, (char *) buff->data, n + 1, fmt, arglist);
		va_end(arglist);
		buff->data[n] = (unsigned char) '\0';
		buff->slen = (int) (strlen) ((char *) buff->data);

		if (buff->slen < n) {
			break;
		}
		if (r > n) {
			n = r;
		} else {
			n += n;
		}
		if (BSTR_OK != balloc(buff, n + 2)) {
			bdestroy(buff);
			return BSTR_ERR;
		}
	}
	r = bconcat(b, buff);
	bdestroy(buff);
	return r;
}

int
bassignformat(bstring b, const char *fmt, ...)
{
	va_list arglist;
	bstring buff;
	int n, r;
	if (!b || !fmt || !b->data || b->mlen <= 0 ||
	    b->slen < 0 || b->slen > b->mlen) {
		return BSTR_ERR;
	}
	/* Since the length is not determinable beforehand, a search is
	 * performed using the truncating "vsnprintf" call (to avoid buffer
	 * overflows) on increasing potential sizes for the output result.
	 */
	n = (int)(2 * strlen(fmt));
	if (n < START_VSNBUFF) {
		n = START_VSNBUFF;
	}
	buff = bfromcstralloc (n + 2, "");
	if (!buff) {
		n = 1;
		buff = bfromcstralloc (n + 2, "");
		if (!buff) {
			return BSTR_ERR;
		}
	}
	while (1) {
		va_start(arglist, fmt);
		exvsnprintf(r, (char *)buff->data, n + 1, fmt, arglist);
		va_end(arglist);
		buff->data[n] = (unsigned char)'\0';
		buff->slen = (int)strlen((char *)buff->data);
		if (buff->slen < n) {
			break;
		}
		if (r > n) {
			n = r;
		} else {
			n += n;
		}
		if (BSTR_OK != balloc(buff, n + 2)) {
			bdestroy(buff);
			return BSTR_ERR;
		}
	}

	r = bassign(b, buff);
	bdestroy(buff);
	return r;
}

bstring
bformat(const char *fmt, ...)
{
	va_list arglist;
	bstring buff;
	int n, r;
	if (!fmt) {
		return NULL;
	}
	/* Since the length is not determinable beforehand, a search is
	 * performed using the truncating "vsnprintf" call (to avoid buffer
	 * overflows) on increasing potential sizes for the output result.
	 */
	n = (int)(2 * strlen(fmt));
	if (n < START_VSNBUFF) {
		n = START_VSNBUFF;
	}
	buff = bfromcstralloc(n + 2, "");
	if (!buff) {
		n = 1;
		buff = bfromcstralloc(n + 2, "");
		if (!buff) {
			return NULL;
		}
	}
	while (1) {
		va_start(arglist, fmt);
		exvsnprintf(r, (char *)buff->data, n + 1, fmt, arglist);
		va_end(arglist);
		buff->data[n] = (unsigned char)'\0';
		buff->slen = (int)strlen((char *)buff->data);
		if (buff->slen < n) {
			break;
		}
		if (r > n) {
			n = r;
		} else {
			n += n;
		}
		if (BSTR_OK != balloc(buff, n + 2)) {
			bdestroy(buff);
			return NULL;
		}
	}
	return buff;
}

int
bvcformata(bstring b, int count, const char *fmt, va_list arg)
{
	int n, r, l;
	if (!b || !fmt || count <= 0 || !b->data ||
	    b->mlen <= 0 || b->slen < 0 || b->slen > b->mlen) {
		return BSTR_ERR;
	}
	if (count > (n = b->slen + count) + 2) {
		return BSTR_ERR;
	}
	if (BSTR_OK != balloc(b, n + 2)) {
		return BSTR_ERR;
	}
	exvsnprintf(r, (char *)b->data + b->slen, count + 2, fmt, arg);
	/* Did the operation complete successfully within bounds? */
	for (l = b->slen; l <= n; l++) {
		if ('\0' == b->data[l]) {
			b->slen = l;
			return BSTR_OK;
		}
	}
	/* Abort, since the buffer was not large enough.  The return value
	 * tries to help set what the retry length should be.
	 */
	b->data[b->slen] = '\0';
	if (r > count + 1) {
		/* Does r specify a particular target length? */
		n = r;
	} else {
		/* If not, just double the size of count */
		n = count + count;
		if (count > n) {
			n = INT_MAX;
		}
	}
	n = -n;
	if (n > BSTR_ERR - 1) {
		n = BSTR_ERR - 1;
	}
	return n;
}
