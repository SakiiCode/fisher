#include <stdlib.h>
#include <stdio.h>
#include <setjmp.h>
#include <math.h>
#include <limits.h>

#ifdef __linux__
#include <sys/mman.h>
#endif

typedef signed int integer;
typedef double doublereal;
typedef signed int logical;
typedef float real;

static const integer c__3 = 3;
static const integer c__1 = 1;
static const integer c__5 = 5;
static const integer c__2 = 2;
static const integer c__4 = 4;
static const integer c__6 = 6;
static const integer c__30 = 30;
static const integer c__7 = 7;
static const integer c__40 = 40;
static const integer c__20 = 20;
static const integer c__501 = 501;
static const integer c__502 = 502;

jmp_buf err_buf;

int prterr_(const integer *icode, const char *mes)
{
    printf("FEXACT ERROR: %d %s\n", *icode, mes);
    longjmp(err_buf, *icode);
}

integer iwork_(integer *iwkmax, integer *iwkpt, integer *number, const integer *itype)
{

    integer ret_val;

    ret_val = *iwkpt;
    if (*itype == 2 || *itype == 3)
    {
        *iwkpt += *number;
    }
    else
    {
        if (ret_val % 2 != 0)
        {
            ++ret_val;
        }
        *iwkpt += *number << 1;
        ret_val /= 2;
    }
    if (*iwkpt > *iwkmax + 1)
    {
        prterr_(&c__40, "Out of workspace.");
    }
    return ret_val;
}

int isort_(integer *n, integer *ix)
{
    integer i, j, m, il[10], kl, it, iu[10], ku, ikey;

    --ix;

    m = 1;
    i = 1;
    j = *n;
L10:
    if (i >= j)
    {
        goto L40;
    }
    kl = i;
    ku = j;
    ikey = i;
    ++j;

L20:
    ++i;
    if (i < j)
    {
        if (ix[ikey] > ix[i])
        {
            goto L20;
        }
    }

L30:
    --j;
    if (ix[j] > ix[ikey])
    {
        goto L30;
    }

    if (i < j)
    {
        it = ix[i];
        ix[i] = ix[j];
        ix[j] = it;
        goto L20;
    }
    it = ix[ikey];
    ix[ikey] = ix[j];
    ix[j] = it;

    if (m < 11)
    {
        if (j - kl < ku - j)
        {
            il[m - 1] = j + 1;
            iu[m - 1] = ku;
            i = kl;
            --j;
        }
        else
        {
            il[m - 1] = kl;
            iu[m - 1] = j - 1;
            i = j + 1;
            j = ku;
        }
        ++m;
        goto L10;
    }
    else
    {
        prterr_(&c__20, "This should never occur.");
    }

L40:
    --m;
    if (m == 0)
    {
        goto L9000;
    }
    i = il[m - 1];
    j = iu[m - 1];
    goto L10;

L9000:
    return 0;
}

int f11act_(integer *irow, integer i1, integer i2, integer *new)
{

    integer i;

    --new;
    --irow;

    for (i = 1; i < i1; ++i)
    {
        new[i] = irow[i];
    }

    for (i = i1; i <= i2; ++i)
    {
        new[i] = irow[i + 1];
    }

    return 0;
}

int f8xact_(integer *irow, integer is, integer i1, integer izero, integer *new)
{

    integer i;

    --new;
    --irow;

    for (i = 1; i < i1; ++i)
    {
        new[i] = irow[i];
    }

    for (i = i1; i < izero; ++i)
    {
        if (is >= irow[i + 1])
        {
            goto L30;
        }
        new[i] = irow[i + 1];
    }

    i = izero;
L30:
    new[i] = is;
L40:
    ++i;
    if (i > izero)
    {
        return 0;
    }
    new[i] = irow[i];
    goto L40;
}

int f4xact_(integer nrow, integer *irow, integer ncol,
            integer *icol, doublereal *dsp, doublereal *fact, integer *icstk,
            integer *ncstk, integer *lstk, integer *mstk, integer *nstk, integer *nrstk, integer *irstk, doublereal *ystk, doublereal tol)
{

    integer icstk_dim1, icstk_offset, irstk_dim1, irstk_offset;

    integer i, j, k, l, m, n;
    doublereal y;
    integer mn, ic1, ir1, ict, nco;
    doublereal amx;
    integer irt, nro, istk;

    irstk_dim1 = nrow;
    irstk_offset = 1 + irstk_dim1;
    irstk -= irstk_offset;
    --irow;
    icstk_dim1 = ncol;
    icstk_offset = 1 + icstk_dim1;
    icstk -= icstk_offset;
    --icol;
    --ncstk;
    --lstk;
    --mstk;
    --nstk;
    --nrstk;
    --ystk;

    if (nrow == 1)
    {
        for (i = 1; i <= ncol; ++i)
        {
            *dsp -= fact[icol[i]];
        }
        goto L9000;
    }

    if (ncol == 1)
    {
        for (i = 1; i <= nrow; ++i)
        {
            *dsp -= fact[irow[i]];
        }
        goto L9000;
    }

    if (nrow * ncol == 4)
    {
        if (irow[2] <= icol[2])
        {
            *dsp = *dsp - fact[irow[2]] - fact[icol[1]] - fact[icol[2] - irow[2]];
        }
        else
        {
            *dsp = *dsp - fact[icol[2]] - fact[irow[1]] - fact[irow[2] - icol[2]];
        }
        goto L9000;
    }

    for (i = 1; i <= nrow; ++i)
    {
        irstk[i + irstk_dim1] = irow[nrow - i + 1];
    }

    for (j = 1; j <= ncol; ++j)
    {
        icstk[j + icstk_dim1] = icol[ncol - j + 1];
    }

    nro = nrow;
    nco = ncol;
    nrstk[1] = nro;
    ncstk[1] = nco;
    ystk[1] = 0.f;
    y = 0.f;
    istk = 1;
    l = 1;
    amx = 0.f;

L50:
    ir1 = irstk[istk * irstk_dim1 + 1];
    ic1 = icstk[istk * icstk_dim1 + 1];
    if (ir1 > ic1)
    {
        if (nro >= nco)
        {
            m = nco - 1;
            n = 2;
        }
        else
        {
            m = nro;
            n = 1;
        }
    }
    else if (ir1 < ic1)
    {
        if (nro <= nco)
        {
            m = nro - 1;
            n = 1;
        }
        else
        {
            m = nco;
            n = 2;
        }
    }
    else
    {
        if (nro <= nco)
        {
            m = nro - 1;
            n = 1;
        }
        else
        {
            m = nco - 1;
            n = 2;
        }
    }

L60:
    if (n == 1)
    {
        i = l;
        j = 1;
    }
    else
    {
        i = 1;
        j = l;
    }

    irt = irstk[i + istk * irstk_dim1];
    ict = icstk[j + istk * icstk_dim1];
    mn = irt;
    if (mn > ict)
    {
        mn = ict;
    }
    y += fact[mn];
    if (irt == ict)
    {
        --nro;
        --nco;
        f11act_(&irstk[istk * irstk_dim1 + 1], i, nro, &irstk[(istk + 1) * irstk_dim1 + 1]);
        f11act_(&icstk[istk * icstk_dim1 + 1], j, nco, &icstk[(istk + 1) * icstk_dim1 + 1]);
    }
    else if (irt > ict)
    {
        --nco;
        f11act_(&icstk[istk * icstk_dim1 + 1], j, nco, &icstk[(istk + 1) * icstk_dim1 + 1]);
        f8xact_(&irstk[istk * irstk_dim1 + 1], irt - ict, i, nro, &irstk[(istk + 1) * irstk_dim1 + 1]);
    }
    else
    {
        --nro;
        f11act_(&irstk[istk * irstk_dim1 + 1], i, nro, &irstk[(istk + 1) * irstk_dim1 + 1]);
        f8xact_(&icstk[istk * icstk_dim1 + 1], ict - irt, j, nco, &icstk[(istk + 1) * icstk_dim1 + 1]);
    }

    if (nro == 1)
    {
        for (k = 1; k <= nco; ++k)
        {
            y += fact[icstk[k + (istk + 1) * icstk_dim1]];
        }
        goto L90;
    }

    if (nco == 1)
    {
        for (k = 1; k <= nro; ++k)
        {
            y += fact[irstk[k + (istk + 1) * irstk_dim1]];
        }
        goto L90;
    }

    lstk[istk] = l;
    mstk[istk] = m;
    nstk[istk] = n;
    ++istk;
    nrstk[istk] = nro;
    ncstk[istk] = nco;
    ystk[istk] = y;
    l = 1;
    goto L50;

L90:
    if (y > amx)
    {
        amx = y;
        if (*dsp - amx <= tol)
        {
            *dsp = 0.f;
            goto L9000;
        }
    }

L100:
    --istk;
    if (istk == 0)
    {
        *dsp -= amx;
        if (*dsp - amx <= tol)
        {
            *dsp = 0.f;
        }
        goto L9000;
    }
    l = lstk[istk] + 1;

L110:
    if (l > mstk[istk])
    {
        goto L100;
    }
    n = nstk[istk];
    nro = nrstk[istk];
    nco = ncstk[istk];
    y = ystk[istk];
    if (n == 1)
    {
        if (irstk[l + istk * irstk_dim1] < irstk[l - 1 + istk * irstk_dim1])
        {
            goto L60;
        }
    }
    else if (n == 2)
    {
        if (icstk[l + istk * icstk_dim1] < icstk[l - 1 + istk * icstk_dim1])
        {
            goto L60;
        }
    }

    ++l;
    goto L110;
L9000:
    return 0;
}

doublereal alogam_(doublereal x, integer *ifault)
{

    doublereal a1 = .918938533204673;
    doublereal a2 = 5.95238095238e-4;
    doublereal a3 = 7.93650793651e-4;
    doublereal a4 = .002777777777778;
    doublereal a5 = .083333333333333;
    doublereal half = .5;
    doublereal zero = 0.;
    doublereal one = 1.;
    doublereal seven = 7.;

    doublereal ret_val;

    doublereal f, y, z;

    ret_val = zero;
    *ifault = 1;
    if (x < zero)
    {
        return ret_val;
    }
    *ifault = 0;
    y = x;
    f = zero;
    if (y >= seven)
    {
        goto L30;
    }
    f = y;
L10:
    y += one;
    if (y >= seven)
    {
        goto L20;
    }
    f *= y;
    goto L10;
L20:
    f = -log(f);
L30:
    z = one / (y * y);
    ret_val = f + (y - half) * log(y) - y + a1 + (((-a2 * z + a3) * z - a4) * z + a5) / y;
    return ret_val;
}

doublereal gammds_(doublereal y, doublereal p, integer *ifault)
{

    doublereal e = 1e-6;
    doublereal zero = 0.;
    doublereal one = 1.;

    doublereal ret_val, d__1, d__2;

    doublereal a, c, f;
    integer ifail;

    *ifault = 1;
    ret_val = zero;
    if (y <= zero || p <= zero)
    {
        return ret_val;
    }
    *ifault = 2;

    d__2 = p + one;
    d__1 = p * log(y) - alogam_(d__2, &ifail) - y;
    f = exp(d__1);
    if (f == zero)
    {
        return ret_val;
    }
    *ifault = 0;

    c = one;
    ret_val = one;
    a = p;
L10:
    a += one;
    c = c * y / a;
    ret_val += c;
    if (c / ret_val > e)
    {
        goto L10;
    }
    ret_val *= f;
    return ret_val;
}

int f5xact_(doublereal pastp, doublereal tol, integer *kval, integer *key, integer ldkey, integer *ipoin, doublereal *stp,
            integer ldstp, integer *ifrq, integer *npoin, integer *nr, integer *nl, integer *ifreq, integer *itop, logical ipsh)
{

    static integer itp;
    integer ird, ipn, itmp;
    doublereal test1, test2;

    --nl;
    --nr;
    --npoin;
    --ifrq;
    --stp;
    --ipoin;
    --key;

    if (ipsh)
    {

        ird = *kval % ldkey + 1;

        for (itp = ird; itp <= ldkey; ++itp)
        {
            if (key[itp] == *kval)
            {
                goto L40;
            }
            if (key[itp] < 0)
            {
                goto L30;
            }
        }

        for (itp = 1; itp < ird; ++itp)
        {
            if (key[itp] == *kval)
            {
                goto L40;
            }
            if (key[itp] < 0)
            {
                goto L30;
            }
        }

        prterr_(&c__6, "LDKEY is too small for this problem.  It is not poss"
                       "ible to estimate the value of LDKEY required, but twice the "
                       "current value may be sufficient.");

    L30:
        key[itp] = *kval;
        ++(*itop);
        ipoin[itp] = *itop;

        if (*itop > ldstp)
        {
            prterr_(&c__7, "LDSTP is too small for this problem.  It is not "
                           "possible to estimate the value of LDSTP required, but tw"
                           "ice the current value may be sufficient.");
        }

        npoin[*itop] = -1;
        nr[*itop] = -1;
        nl[*itop] = -1;
        stp[*itop] = pastp;
        ifrq[*itop] = *ifreq;
        goto L9000;
    }

L40:
    ipn = ipoin[itp];
    test1 = pastp - tol;
    test2 = pastp + tol;

L50:
    if (stp[ipn] < test1)
    {
        ipn = nl[ipn];
        if (ipn > 0)
        {
            goto L50;
        }
    }
    else if (stp[ipn] > test2)
    {
        ipn = nr[ipn];
        if (ipn > 0)
        {
            goto L50;
        }
    }
    else
    {
        ifrq[ipn] += *ifreq;
        goto L9000;
    }

    ++(*itop);
    if (*itop > ldstp)
    {
        prterr_(&c__7, "LDSTP is too small for this problem.  It is not poss"
                       "ible to estimate the value of LDSTP rerquired, but twice the"
                       " current value may be sufficient.");
        goto L9000;
    }

    ipn = ipoin[itp];
    itmp = ipn;
L60:
    if (stp[ipn] < test1)
    {
        itmp = ipn;
        ipn = nl[ipn];
        if (ipn > 0)
        {
            goto L60;
        }
        else
        {
            nl[itmp] = *itop;
        }
    }
    else if (stp[ipn] > test2)
    {
        itmp = ipn;
        ipn = nr[ipn];
        if (ipn > 0)
        {
            goto L60;
        }
        else
        {
            nr[itmp] = *itop;
        }
    }

    npoin[*itop] = npoin[itmp];
    npoin[itmp] = *itop;
    stp[*itop] = pastp;
    ifrq[*itop] = *ifreq;
    nl[*itop] = -1;
    nr[*itop] = -1;

L9000:
    return 0;
}

int f7xact_(integer nrow, integer *imax, integer *idif,
            integer *k, integer *ks, integer *iflag)
{

    integer i__1;

    integer i, m, k1, mm;

    --idif;
    --imax;

    *iflag = 0;

    if (*ks == 0)
    {
    L10:
        ++(*ks);
        if (idif[*ks] == imax[*ks])
        {
            goto L10;
        }
    }

    if (idif[*k] > 0 && *k > *ks)
    {
        --idif[*k];
    L30:
        --(*k);
        if (imax[*k] == 0)
        {
            goto L30;
        }
        m = *k;

    L40:
        if (idif[m] >= imax[m])
        {
            --m;
            goto L40;
        }
        ++idif[m];

        if (m == *ks)
        {
            if (idif[m] == imax[m])
            {
                *ks = *k;
            }
        }
    }
    else
    {

    L50:
        for (k1 = *k + 1; k1 <= nrow; ++k1)
        {
            if (idif[k1] > 0)
            {
                goto L70;
            }
        }
        *iflag = 1;
        goto L9000;

    L70:
        mm = 1;
        i__1 = *k;
        for (i = 1; i <= i__1; ++i)
        {
            mm += idif[i];
            idif[i] = 0;
        }
        *k = k1;
    L90:
        --(*k);

        i__1 = imax[*k];
        m = ((mm) <= (i__1) ? (mm) : (i__1));
        idif[*k] = m;
        mm -= m;
        if (mm > 0 && *k != 1)
        {
            goto L90;
        }

        if (mm > 0)
        {
            if (k1 != nrow)
            {
                *k = k1;
                goto L50;
            }
            *iflag = 1;
            goto L9000;
        }

        --idif[k1];
        *ks = 0;
    L100:
        ++(*ks);
        if (*ks > *k)
        {
            goto L9000;
        }
        if (idif[*ks] >= imax[*ks])
        {
            goto L100;
        }
    }

L9000:
    return 0;
}

int f6xact_(integer nrow, integer *irow, integer *iflag,
            integer *kyy, integer *key, integer ldkey, integer *last, integer *ipn)
{
    integer j, kval;

    --key;
    --kyy;
    --irow;

L10:
    ++(*last);
    if (*last <= ldkey)
    {
        if (key[*last] < 0)
        {
            goto L10;
        }

        kval = key[*last];
        key[*last] = -9999;
        for (j = nrow; j >= 2; --j)
        {
            irow[j] = kval / kyy[j];
            kval -= irow[j] * kyy[j];
        }
        irow[1] = kval;
        *ipn = *last;
    }
    else
    {
        *last = 0;
        *iflag = 3;
    }
    return 0;
}

doublereal f9xact_(integer n, integer mm, integer *ir, doublereal *fact)
{

    doublereal ret_val;

    integer k;

    ret_val = fact[mm];
    for (k = 0; k < n; ++k)
    {
        ret_val -= fact[ir[k]];
    }

    return ret_val;
}

int f10act_(integer nrow, integer *irow, integer ncol,
            integer *icol, doublereal *val, logical *xmin, doublereal *fact,
            integer *nd, integer *ne, integer *m)
{

    integer i, is, ix, nrw1;

    --m;
    --ne;
    --nd;
    --icol;
    --irow;

    for (i = 1; i < nrow; ++i)
    {
        nd[i] = 0;
    }

    is = icol[1] / nrow;
    ne[1] = is;
    ix = icol[1] - nrow * is;
    m[1] = ix;
    if (ix != 0)
    {
        ++nd[ix];
    }

    for (i = 2; i <= ncol; ++i)
    {
        ix = icol[i] / nrow;
        ne[i] = ix;
        is += ix;
        ix = icol[i] - nrow * ix;
        m[i] = ix;
        if (ix != 0)
        {
            ++nd[ix];
        }
    }

    for (i = nrow - 2; i >= 1; --i)
    {
        nd[i] += nd[i + 1];
    }

    ix = 0;
    nrw1 = nrow + 1;
    for (i = nrow; i >= 2; --i)
    {
        ix = ix + is + nd[nrw1 - i] - irow[i];
        if (ix < 0)
        {
            return 0;
        }
    }

    for (i = 1; i <= ncol; ++i)
    {
        ix = ne[i];
        is = m[i];
        *val = *val + is * fact[ix + 1] + (nrow - is) * fact[ix];
    }
    *xmin = 1;

    return 0;
}

int f3xact_(integer nrow, integer *irow, integer ncol,
            integer *icol, doublereal *dlp, integer mm, doublereal *fact,
            integer *ico, integer *iro, integer *it, integer *lb, integer *nr,
            integer *nt, integer *nu, integer *itc, integer *ist, doublereal *stv,
            doublereal *alen, doublereal tol)
{

    integer ldst = 200;
    integer nst = 0;
    integer nitc = 0;

    doublereal d__2;

    integer i, k;
    doublereal v;
    integer n11, n12, ii, nn, ks, ic1, ic2, nc1, nn1, nr1, nco;
    doublereal val;
    integer nct, ipn, irl, key = 0, lev, itp, nro;
    doublereal vmn;
    integer nrt, kyy, nc1s;
    logical xmin;

    --stv;
    --ist;
    --itc;
    --nu;
    --nt;
    --nr;
    --lb;
    --it;
    --iro;
    --ico;
    --icol;
    --irow;

    for (i = 0; i <= ncol; ++i)
    {
        alen[i] = 0.f;
    }
    for (i = 1; i <= 400; ++i)
    {
        ist[i] = -1;
    }

    if (nrow <= 1)
    {
        if (nrow > 0)
        {
            *dlp -= fact[icol[1]];
            for (i = 2; i <= ncol; ++i)
            {
                *dlp -= fact[icol[i]];
            }
        }
        goto L9000;
    }

    if (ncol <= 1)
    {
        if (ncol > 0)
        {
            *dlp = *dlp - fact[irow[1]] - fact[irow[2]];
            for (i = 3; i <= nrow; ++i)
            {
                *dlp -= fact[irow[i]];
            }
        }
        goto L9000;
    }

    if (nrow * ncol == 4)
    {
        n11 = (irow[1] + 1) * (icol[1] + 1) / (mm + 2);
        n12 = irow[1] - n11;
        *dlp = *dlp - fact[n11] - fact[n12] - fact[icol[1] - n11] - fact[icol[2] - n12];
        goto L9000;
    }

    val = 0.f;
    xmin = 0;
    if (irow[nrow] <= irow[1] + ncol)
    {
        f10act_(nrow, &irow[1], ncol, &icol[1], &val, &xmin, fact, &lb[1], &nu[1], &nr[1]);
    }
    if (!xmin)
    {
        if (icol[ncol] <= icol[1] + nrow)
        {
            f10act_(ncol, &icol[1], nrow, &irow[1], &val, &xmin, fact, &lb[1],
                    &nu[1], &nr[1]);
        }
    }

    if (xmin)
    {
        *dlp -= val;
        goto L9000;
    }

    nn = mm;

    if (nrow >= ncol)
    {
        nro = nrow;
        nco = ncol;

        for (i = 1; i <= nrow; ++i)
        {
            iro[i] = irow[i];
        }

        ico[1] = icol[1];
        nt[1] = nn - ico[1];
        for (i = 2; i <= ncol; ++i)
        {
            ico[i] = icol[i];
            nt[i] = nt[i - 1] - ico[i];
        }
    }
    else
    {
        nro = ncol;
        nco = nrow;

        ico[1] = irow[1];
        nt[1] = nn - ico[1];
        for (i = 2; i <= nrow; ++i)
        {
            ico[i] = irow[i];
            nt[i] = nt[i - 1] - ico[i];
        }

        for (i = 1; i <= ncol; ++i)
        {
            iro[i] = icol[i];
        }
    }

    vmn = 1e10;
    nc1s = nco - 1;
    irl = 1;
    ks = 0;
    k = ldst;
    kyy = ico[nco] + 1;
    goto L100;

L90:
    xmin = 0;
    if (iro[nro] <= iro[irl] + nco)
    {
        f10act_(nro, &iro[irl], nco, &ico[1], &val, &xmin, fact, &lb[1], &nu[1], &nr[1]);
    }
    if (!xmin)
    {
        if (ico[nco] <= ico[1] + nro)
        {
            f10act_(nco, &ico[1], nro, &iro[irl], &val, &xmin, fact, &lb[1],
                    &nu[1], &nr[1]);
        }
    }

    if (xmin)
    {
        if (val < vmn)
        {
            vmn = val;
        }
        goto L200;
    }

L100:
    lev = 1;
    nr1 = nro - 1;
    nrt = iro[irl];
    nct = ico[1];
    lb[1] = (integer)((doublereal)((nrt + 1) * (nct + 1)) / (doublereal)(nn + nr1 * nc1s + 1) - tol) - 1;
    nu[1] = (integer)((doublereal)((nrt + nc1s) * (nct + nr1)) / (doublereal)(nn + nr1 + nc1s)) - lb[1] + 1;
    nr[1] = nrt - lb[1];

L110:
    --nu[lev];
    if (nu[lev] == 0)
    {
        if (lev == 1)
        {
            goto L200;
        }
        --lev;
        goto L110;
    }
    ++lb[lev];
    --nr[lev];
L120:
    alen[lev] = alen[lev - 1] + fact[lb[lev]];
    if (lev < nc1s)
    {
        nn1 = nt[lev];
        nrt = nr[lev];
        ++lev;
        nc1 = nco - lev;
        nct = ico[lev];
        lb[lev] = (integer)((doublereal)((nrt + 1) * (nct + 1)) / (doublereal)(nn1 + nr1 * nc1 + 1) - tol);
        nu[lev] = (integer)((doublereal)((nrt + nc1) * (nct + nr1)) / (doublereal)(nn1 + nr1 + nc1) - lb[lev] + 1);
        nr[lev] = nrt - lb[lev];
        goto L120;
    }
    alen[nco] = alen[lev] + fact[nr[lev]];
    lb[nco] = nr[lev];

    v = val + alen[nco];
    if (nro == 2)
    {

        v = v + fact[ico[1] - lb[1]] + fact[ico[2] - lb[2]];
        for (i = 3; i <= nco; ++i)
        {
            v += fact[ico[i] - lb[i]];
        }
        if (v < vmn)
        {
            vmn = v;
        }
    }
    else if (nro == 3 && nco == 2)
    {

        nn1 = nn - iro[irl] + 2;
        ic1 = ico[1] - lb[1];
        ic2 = ico[2] - lb[2];
        n11 = (iro[irl + 1] + 1) * (ic1 + 1) / nn1;
        n12 = iro[irl + 1] - n11;
        v = v + fact[n11] + fact[n12] + fact[ic1 - n11] + fact[ic2 - n12];
        if (v < vmn)
        {
            vmn = v;
        }
    }
    else
    {

        for (i = 1; i <= nco; ++i)
        {
            it[i] = ico[i] - lb[i];
        }

        if (nco == 2)
        {
            if (it[1] > it[2])
            {
                ii = it[1];
                it[1] = it[2];
                it[2] = ii;
            }
        }
        else if (nco == 3)
        {
            ii = it[1];
            if (ii > it[3])
            {
                if (ii > it[2])
                {
                    if (it[2] > it[3])
                    {
                        it[1] = it[3];
                        it[3] = ii;
                    }
                    else
                    {
                        it[1] = it[2];
                        it[2] = it[3];
                        it[3] = ii;
                    }
                }
                else
                {
                    it[1] = it[3];
                    it[3] = it[2];
                    it[2] = ii;
                }
            }
            else if (ii > it[2])
            {
                it[1] = it[2];
                it[2] = ii;
            }
            else if (it[2] > it[3])
            {
                ii = it[2];
                it[2] = it[3];
                it[3] = ii;
            }
        }
        else
        {
            isort_(&nco, &it[1]);
        }

        double dkyy = (double)kyy;
        double dkey = it[1] * dkyy + it[2];
        for (i = 3; i <= nco; ++i)
        {
            dkey = it[i] + dkey * dkyy;
        }
        if (dkey > 0x7fffffff)
        {
            prterr_(&c__502, "The hash table key cannot be computed because the la"
                             "rgest key is larger than the largest representable integer. "
                             " The algorithm cannot proceed.");
        }
        else
        {
            key = (int)dkey;
        }

        ipn = key % ldst + 1;

        ii = ks + ipn;
        for (itp = ipn; itp <= ldst; ++itp)
        {
            if (ist[ii] < 0)
            {
                goto L180;
            }
            else if (ist[ii] == key)
            {
                goto L190;
            }
            ++ii;
        }

        ii = ks + 1;
        for (itp = 1; itp < ipn; ++itp)
        {
            if (ist[ii] < 0)
            {
                goto L180;
            }
            else if (ist[ii] == key)
            {
                goto L190;
            }
            ++ii;
        }

        prterr_(&c__30, "Stack length exceeded in f3xact.  This problem shou"
                        "ld not occur.");

    L180:
        ist[ii] = key;
        stv[ii] = v;
        ++nst;
        ii = nst + ks;
        itc[ii] = itp;
        goto L110;

    L190:

        d__2 = stv[ii];
        stv[ii] = ((v) <= (d__2) ? (v) : (d__2));
    }
    goto L110;

L200:
    if (nitc > 0)
    {

        itp = itc[nitc + k] + k;
        --nitc;
        val = stv[itp];
        key = ist[itp];
        ist[itp] = -1;

        for (i = nco; i >= 2; --i)
        {
            ico[i] = key % kyy;
            key /= kyy;
        }
        ico[1] = key;

        nt[1] = nn - ico[1];
        for (i = 2; i <= nco; ++i)
        {
            nt[i] = nt[i - 1] - ico[i];
        }
        goto L90;
    }
    else if (nro > 2 && nst > 0)
    {

        nitc = nst;
        nst = 0;
        k = ks;
        ks = ldst - ks;
        nn -= iro[irl];
        ++irl;
        --nro;
        goto L200;
    }

    *dlp -= vmn;
L9000:
    return 0;
}

int f2xact_(integer nrow, integer ncol, doublereal *table,
            integer *ldtabl, doublereal *expect, doublereal *percnt, doublereal *emin, doublereal *prt, doublereal *pre, doublereal *fact, integer *ico, integer *iro, integer *kyy, integer *idif, integer *irn, integer *key, integer ldkey, integer *ipoin, doublereal *stp, integer ldstp,
            integer *ifrq, doublereal *dlp, doublereal *dsp, doublereal *tm,
            integer *key2, integer *iwk, doublereal *rwk)
{

    integer imax = 2147483647;
    real amiss = -12345.f;
    doublereal tol = 3.45254e-7;
    real emx = 1e30f;

    integer table_dim1, table_offset, i__1;
    doublereal d__1, d__3, d__4;

    integer i, j, k, n, k1;
    doublereal dd, df;
    integer i31, i32, i33, i34, i35, i36, i37, i38, i39, i41, i42, i43,
        i44, i45, i46, i47, i48, ii, kb, kd, ks;
    doublereal pv;
    integer i310, i311;
    doublereal ddf;
    integer nco, nrb;
    doublereal emn, drn, dro, obs;
    integer ipn, ipo, itp = 0, nro;
    doublereal tmp = 0, obs2, obs3;
    integer nro2, kval, kmax, jkey, last;
    logical ipsh;
    integer itmp;
    doublereal dspt;
    integer itop, jstp, ntot, jstp2, jstp3, jstp4, iflag, ncell, ifreq;
    logical chisq = 0;
    integer ikkey;
    doublereal pastp;
    integer ikstp;
    integer ikstp2;
    integer ifault;

    table_dim1 = *ldtabl;
    table_offset = 1 + table_dim1;
    table -= table_offset;
    --ico;
    --iro;
    --kyy;
    --idif;
    --irn;
    --key;
    --ipoin;
    --stp;
    --ifrq;
    --dlp;
    --dsp;
    --tm;
    --key2;
    --iwk;
    --rwk;

    i__1 = 2 * ldkey;
    for (i = 1; i <= i__1; ++i)
    {
        key[i] = -9999;
        key2[i] = -9999;
    }

    *pre = 0.f;
    itop = 0;
    if (*expect > 0.)
    {
        emn = *emin;
    }
    else
    {
        emn = emx;
    }

    k = ((nrow) >= (ncol) ? (nrow) : (ncol));

    i31 = 1;
    i32 = i31 + k;
    i33 = i32 + k;
    i34 = i33 + k;
    i35 = i34 + k;
    i36 = i35 + k;
    i37 = i36 + k;
    i38 = i37 + k;
    i39 = i38 + 400;
    i310 = 1;
    i311 = 401;

    k = nrow + ncol + 1;
    i41 = 1;
    i42 = i41 + k;
    i43 = i42 + k;
    i44 = i43 + k;
    i45 = i44 + k;
    i46 = i45 + k;
    i47 = i46 + k * ((nrow) >= (ncol) ? (nrow) : (ncol));
    i48 = 1;

    if (nrow > *ldtabl)
    {
        prterr_(&c__1, "NROW must be less than or equal to LDTABL.");
    }
    if (ncol <= 1)
    {
        prterr_(&c__4, "NCOL must be greater than 1.0.");
    }

    ntot = 0;
    for (i = 1; i <= nrow; ++i)
    {
        iro[i] = 0;
        for (j = 1; j <= ncol; ++j)
        {
            if (table[i + j * table_dim1] < -1e-4)
            {
                prterr_(&c__2, "All elements of TABLE must be positive.");
            }
            iro[i] += (int)(table[i + j * table_dim1]);
            ntot += (int)(table[i + j * table_dim1]);
        }
    }

    if (ntot == 0)
    {
        prterr_(&c__3, "All elements of TABLE are zero.  PRT and PRE are set"
                       " to missing values (NaN, not a number).");
        *prt = amiss;
        *pre = amiss;
        goto L9000;
    }

    for (i = 1; i <= ncol; ++i)
    {
        ico[i] = 0;
        for (j = 1; j <= nrow; ++j)
        {
            ico[i] += (int)(table[j + i * table_dim1]);
        }
    }

    isort_(&nrow, &iro[1]);
    isort_(&ncol, &ico[1]);

    if (nrow > ncol)
    {
        nro = ncol;
        nco = nrow;

        for (i = 1; i <= nrow; ++i)
        {
            itmp = iro[i];
            if (i <= ncol)
            {
                iro[i] = ico[i];
            }
            ico[i] = itmp;
        }
    }
    else
    {
        nro = nrow;
        nco = ncol;
    }

    kyy[1] = 1;
    for (i = 2; i <= nro; ++i)
    {

        if (iro[i - 1] + 1 <= imax / kyy[i - 1])
        {
            kyy[i] = kyy[i - 1] * (iro[i - 1] + 1);
            j /= kyy[i - 1];
        }
        else
        {
            prterr_(&c__5, "The hash table key cannot be computed because th"
                           "e largest key is larger than the largest representable i"
                           "nteger.  The algorithm cannot proceed.");
        }
    }

    if (iro[nro - 1] + 1 <= imax / kyy[nro - 1])
    {
        kmax = (iro[nro] + 1) * kyy[nro - 1];
    }
    else
    {
        prterr_(&c__5, "The hash table key cannot be computed because the la"
                       "rgest key is larger than the largest representable integer. "
                       " The algorithm cannot proceed.");
        goto L9000;
    }

    if (iro[nro] + 1 > imax / kyy[nro])
    {
        prterr_(&c__501, "The hash table key cannot be computed because the la"
                         "rgest key is larger than the largest representable integer. "
                         " The algorithm cannot proceed.");
    }

    fact[0] = 0.;
    fact[1] = 0.;
    fact[2] = log(2.);
    for (i = 3; i <= ntot; i += 2)
    {
        fact[i] = fact[i - 1] + log((doublereal)i);
        j = i + 1;
        if (j <= ntot)
        {
            fact[j] = fact[i] + fact[2] + fact[j / 2] - fact[j / 2 - 1];
        }
    }

    obs = tol;
    ntot = 0;
    for (j = 1; j <= nco; ++j)
    {
        dd = 0.f;
        for (i = 1; i <= nro; ++i)
        {
            if (nrow <= ncol)
            {
                dd += fact[(int)(table[i + j * table_dim1])];
                ntot += (int)(table[i + j * table_dim1]);
            }
            else
            {
                dd += fact[(int)(table[j + i * table_dim1])];
                ntot += (int)(table[j + i * table_dim1]);
            }
        }
        obs = obs + fact[ico[j]] - dd;
    }

    dro = f9xact_(nro, ntot, &iro[1], fact);
    *prt = exp(obs - dro);

    k = nco;
    last = ldkey + 1;
    jkey = ldkey + 1;
    jstp = ldstp + 1;
    jstp2 = ldstp * 3 + 1;
    jstp3 = (ldstp << 2) + 1;
    jstp4 = ldstp * 5 + 1;
    ikkey = 0;
    ikstp = 0;
    ikstp2 = ldstp << 1;
    ipo = 1;
    ipoin[1] = 1;
    stp[1] = 0.f;
    ifrq[1] = 1;
    ifrq[ikstp2 + 1] = -1;

L110:
    kb = nco - k + 1;
    ks = 0;
    n = ico[kb];
    kd = nro + 1;
    kmax = nro;

    for (i = 1; i <= nro; ++i)
    {
        idif[i] = 0;
    }

L130:
    --kd;

    i__1 = iro[kd];
    ntot = ((n) <= (i__1) ? (n) : (i__1));
    idif[kd] = ntot;
    if (idif[kmax] == 0)
    {
        --kmax;
    }
    n -= ntot;
    if (n > 0 && kd != 1)
    {
        goto L130;
    }
    if (n != 0)
    {
        goto L310;
    }

    k1 = k - 1;
    n = ico[kb];
    ntot = 0;
    for (i = kb + 1; i <= nco; ++i)
    {
        ntot += ico[i];
    }

L150:
    for (i = 1; i <= nro; ++i)
    {
        irn[i] = iro[i] - idif[i];
    }

    if (k1 > 1)
    {
        if (nro == 2)
        {
            if (irn[1] > irn[2])
            {
                ii = irn[1];
                irn[1] = irn[2];
                irn[2] = ii;
            }
        }
        else if (nro == 3)
        {
            ii = irn[1];
            if (ii > irn[3])
            {
                if (ii > irn[2])
                {
                    if (irn[2] > irn[3])
                    {
                        irn[1] = irn[3];
                        irn[3] = ii;
                    }
                    else
                    {
                        irn[1] = irn[2];
                        irn[2] = irn[3];
                        irn[3] = ii;
                    }
                }
                else
                {
                    irn[1] = irn[3];
                    irn[3] = irn[2];
                    irn[2] = ii;
                }
            }
            else if (ii > irn[2])
            {
                irn[1] = irn[2];
                irn[2] = ii;
            }
            else if (irn[2] > irn[3])
            {
                ii = irn[2];
                irn[2] = irn[3];
                irn[3] = ii;
            }
        }
        else
        {
            for (j = 2; j <= nro; ++j)
            {
                i = j - 1;
                ii = irn[j];
            L170:
                if (ii < irn[i])
                {
                    irn[i + 1] = irn[i];
                    --i;
                    if (i > 0)
                    {
                        goto L170;
                    }
                }
                irn[i + 1] = ii;
            }
        }

        for (i = 1; i <= nro; ++i)
        {
            if (irn[i] != 0)
            {
                goto L200;
            }
        }
    L200:
        nrb = i;
        nro2 = nro - i + 1;
    }
    else
    {
        nrb = 1;
        nro2 = nro;
    }

    ddf = f9xact_(nro, n, &idif[1], fact);
    drn = f9xact_(nro2, ntot, &irn[nrb], fact) - dro + ddf;

    if (k1 > 1)
    {
        kval = irn[1] + irn[2] * kyy[2];
        for (i = 3; i <= nro; ++i)
        {
            kval += irn[i] * kyy[i];
        }

        i = kval % (ldkey << 1) + 1;

        i__1 = ldkey << 1;
        for (itp = i; itp <= i__1; ++itp)
        {
            ii = key2[itp];
            if (ii == kval)
            {
                goto L240;
            }
            else if (ii < 0)
            {
                key2[itp] = kval;
                dlp[itp] = 1.;
                dsp[itp] = 1.;
                goto L240;
            }
        }

        for (itp = 1; itp < i; ++itp)
        {
            ii = key2[itp];
            if (ii == kval)
            {
                goto L240;
            }
            else if (ii < 0)
            {
                key2[itp] = kval;
                dlp[itp] = 1.f;
                goto L240;
            }
        }

        prterr_(&c__6, "LDKEY is too small.  It is not possible to give thev"
                       "alue of LDKEY required, but you could try doubling LDKEY (an"
                       "d possibly LDSTP).");
    }

L240:
    ipsh = 1;

    ipn = ipoin[ipo + ikkey];
    pastp = stp[ipn + ikstp];
    ifreq = ifrq[ipn + ikstp];

    if (k1 > 1)
    {
        obs2 = obs - fact[ico[kb + 1]] - fact[ico[kb + 2]] - ddf;
        for (i = 3; i <= k1; ++i)
        {
            obs2 -= fact[ico[kb + i]];
        }

        if (dlp[itp] > 0.)
        {
            dspt = obs - obs2 - ddf;

            dlp[itp] = 0.;
            f3xact_(nro2, &irn[nrb], k1, &ico[kb + 1], &dlp[itp], ntot,
                    fact, &iwk[i31], &iwk[i32], &iwk[i33], &iwk[i34], &iwk[i35], &iwk[i36], &iwk[i37], &iwk[i38], &iwk[i39], &rwk[i310], &rwk[i311], tol);

            d__1 = dlp[itp];
            dlp[itp] = ((0.) <= (d__1) ? (0.) : (d__1));

            dsp[itp] = dspt;
            f4xact_(nro2, &irn[nrb], k1, &ico[kb + 1], &dsp[itp], fact, &iwk[i47], &iwk[i41], &iwk[i42], &iwk[i43], &iwk[i44], &iwk[i45], &iwk[i46], &rwk[i48], tol);

            d__1 = dsp[itp] - dspt;
            dsp[itp] = ((0.) <= (d__1) ? (0.) : (d__1));

            if ((doublereal)(irn[nrb] * ico[kb + 1]) / (doublereal)ntot >
                emn)
            {
                ncell = 0.f;
                for (i = 1; i <= nro2; ++i)
                {
                    for (j = 1; j <= k1; ++j)
                    {
                        if ((doublereal)(irn[nrb + i - 1] * ico[kb + j]) >=
                            ntot * *expect)
                        {
                            ++ncell;
                        }
                    }
                }
                if ((doublereal)(ncell * 100) >= k1 * nro2 * *percnt)
                {
                    tmp = 0.f;
                    for (i = 1; i <= nro2; ++i)
                    {
                        tmp = tmp + fact[irn[nrb + i - 1]] - fact[irn[nrb + i - 1] - 1];
                    }
                    tmp *= k1 - 1;
                    for (j = 1; j <= k1; ++j)
                    {
                        tmp += (nro2 - 1) * (fact[ico[kb + j]] - fact[ico[kb + j] - 1]);
                    }
                    df = (doublereal)((nro2 - 1) * (k1 - 1));
                    tmp += df * 1.83787706640934548356065947281;
                    tmp -= (nro2 * k1 - 1) * (fact[ntot] - fact[ntot - 1]);
                    tm[itp] = (obs - dro) * -2. - tmp;
                }
                else
                {

                    tm[itp] = -9876.;
                }
            }
            else
            {
                tm[itp] = -9876.;
            }
        }
        obs3 = obs2 - dlp[itp];
        obs2 -= dsp[itp];
        if (tm[itp] == -9876.)
        {
            chisq = 0;
        }
        else
        {
            chisq = 1;
            tmp = tm[itp];
        }
    }
    else
    {
        obs2 = obs - drn - dro;
        obs3 = obs2;
    }

L300:
    if (pastp <= obs3)
    {

        *pre += (doublereal)ifreq * exp(pastp + drn);
    }
    else if (pastp < obs2)
    {
        if (chisq)
        {
            df = (doublereal)((nro2 - 1) * (k1 - 1));

            d__3 = tmp + (pastp + drn) * 2.;
            d__1 = ((0.) >= (d__3) ? (0.) : (d__3)) / 2.;
            d__4 = df / 2.;
            pv = 1.f - gammds_(d__1, d__4, &ifault);
            *pre += (doublereal)ifreq * exp(pastp + drn) * pv;
        }
        else
        {

            f5xact_(pastp + ddf, tol, &kval, &key[jkey], ldkey, &ipoin[jkey], &stp[jstp], ldstp, &ifrq[jstp], &ifrq[jstp2], &ifrq[jstp3], &ifrq[jstp4], &ifreq, &itop, ipsh);
            ipsh = 0;
        }
    }

    ipn = ifrq[ipn + ikstp2];
    if (ipn > 0)
    {
        pastp = stp[ipn + ikstp];
        ifreq = ifrq[ipn + ikstp];
        goto L300;
    }

    f7xact_(kmax, &iro[1], &idif[1], &kd, &ks, &iflag);
    if (iflag != 1)
    {
        goto L150;
    }

L310:
    iflag = 1;
    f6xact_(nro, &iro[1], &iflag, &kyy[1], &key[ikkey + 1], ldkey, &last, &ipo);

    if (iflag == 3)
    {
        --k;
        itop = 0;
        ikkey = jkey - 1;
        ikstp = jstp - 1;
        ikstp2 = jstp2 - 1;
        jkey = ldkey - jkey + 2;
        jstp = ldstp - jstp + 2;
        jstp2 = (ldstp << 1) + jstp;
        i__1 = ldkey << 1;
        for (i = 1; i <= i__1; ++i)
        {
            key2[i] = -9999;
        }
        if (k >= 2)
        {
            goto L310;
        }
    }
    else
    {
        goto L110;
    }

L9000:
    return 0;
}

int fexact_(integer nrow, integer ncol, doublereal *table,
            integer ldtabl, doublereal *expect, doublereal *percnt, doublereal *emin, doublereal *prt, doublereal *pre, integer ws)
{

    integer table_dim1, table_offset, i__1, i__2, i__3;
    doublereal *equiv_1 = 0;
    size_t allocation = ws * sizeof(integer);
#ifdef __linux__
    posix_memalign((void **)&equiv_1, 1 << 21, allocation);
    if (equiv_1 == 0)
    {
        equiv_1 = malloc(allocation);
    }
    else
    {
        madvise(equiv_1, allocation,
                14);
    }
#else
    equiv_1 = malloc(allocation);
#endif
    integer i, j, k, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, kk,
        i3a, i3b, i3c, i9a, nco, nro, numb, iiwk;

    integer irwk;

    integer mult, ntot;

    integer ireal, ldkey;
    real amiss;
    integer ldstp;
    integer iwkpt;
    integer iwkmax;

    table_dim1 = ldtabl;
    table_offset = 1 + table_dim1;
    table -= table_offset;

    iwkmax = ws;

    mult = 30;

    ireal = 4;

    amiss = -12345.f;

    iwkpt = 1;

    int errorcode = setjmp(err_buf);
    if (errorcode != 0)
    {
        free(equiv_1);
        return -errorcode;
    }

    if (nrow > ldtabl)
    {
        prterr_(&c__1, "NROW must be less than or equal to LDTABL.");
    }
    ntot = 0;
    for (i = 1; i <= nrow; ++i)
    {
        for (j = 1; j <= ncol; ++j)
        {
            if (table[i + j * table_dim1] < 0.)
            {
                prterr_(&c__2, "All elements of TABLE must be positive.");
            }
            ntot = (integer)(ntot + table[i + j * table_dim1]);
        }
    }
    if (ntot == 0)
    {
        prterr_(&c__3, "All elements of TABLE are zero.  PRT and PRE are set"
                       " to missing values (NaN, not a number).");
        *prt = amiss;
        *pre = amiss;
        goto L9000;
    }

    nco = ((nrow) >= (ncol) ? (nrow) : (ncol));
    nro = nrow + ncol - nco;
    k = nrow + ncol + 1;
    kk = k * ((nrow) >= (ncol) ? (nrow) : (ncol));

    i__1 = ntot + 1;
    i1 = iwork_(&iwkmax, &iwkpt, &i__1, &ireal);
    i2 = iwork_(&iwkmax, &iwkpt, &nco, &c__2);
    i3 = iwork_(&iwkmax, &iwkpt, &nco, &c__2);
    i3a = iwork_(&iwkmax, &iwkpt, &nco, &c__2);
    i3b = iwork_(&iwkmax, &iwkpt, &nro, &c__2);
    i3c = iwork_(&iwkmax, &iwkpt, &nro, &c__2);

    i__2 = k * 5 + (kk << 1), i__3 = ((nrow) >= (ncol) ? (nrow) : (ncol)) * 7 + 800;
    i__1 = ((i__2) >= (i__3) ? (i__2) : (i__3));
    iiwk = iwork_(&iwkmax, &iwkpt, &i__1, &c__2);

    i__2 = ((nrow) >= (ncol) ? (nrow) : (ncol)) + 401;
    i__1 = ((i__2) >= (k) ? (i__2) : (k));
    irwk = iwork_(&iwkmax, &iwkpt, &i__1, &ireal);

    if (ireal == 4)
    {
        numb = mult * 10 + 18;
        ldkey = (iwkmax - iwkpt + 1) / numb;
    }
    else
    {

        numb = (mult << 3) + 12;
        ldkey = (iwkmax - iwkpt + 1) / numb;
    }

    ldstp = mult * ldkey;
    i__1 = ldkey << 1;
    i4 = iwork_(&iwkmax, &iwkpt, &i__1, &c__2);
    i__1 = ldkey << 1;
    i5 = iwork_(&iwkmax, &iwkpt, &i__1, &c__2);
    i__1 = ldstp << 1;
    i6 = iwork_(&iwkmax, &iwkpt, &i__1, &ireal);
    i__1 = ldstp * 6;
    i7 = iwork_(&iwkmax, &iwkpt, &i__1, &c__2);
    i__1 = ldkey << 1;
    i8 = iwork_(&iwkmax, &iwkpt, &i__1, &ireal);
    i__1 = ldkey << 1;
    i9 = iwork_(&iwkmax, &iwkpt, &i__1, &ireal);
    i__1 = ldkey << 1;
    i9a = iwork_(&iwkmax, &iwkpt, &i__1, &ireal);
    i__1 = ldkey << 1;
    i10 = iwork_(&iwkmax, &iwkpt, &i__1, &c__2);

    f2xact_(nrow, ncol, &table[table_offset], &ldtabl, expect, percnt, emin,
            prt, pre, &(equiv_1)[i1 - 1], &((integer *)equiv_1)[i2 - 1], &((integer *)equiv_1)[i3 - 1], &((integer *)equiv_1)[i3a - 1], &((integer *)equiv_1)[i3b - 1], &((integer *)equiv_1)[i3c - 1], &((integer *)equiv_1)[i4 - 1], ldkey, &((integer *)equiv_1)[i5 - 1], &(equiv_1)[i6 - 1], ldstp, &((integer *)equiv_1)[i7 - 1], &(equiv_1)[i8 - 1],
            &(equiv_1)[i9 - 1], &(equiv_1)[i9a - 1], &((integer *)equiv_1)[i10 - 1], &((integer *)equiv_1)[iiwk - 1], &(equiv_1)[irwk - 1]);

L9000:

    free(equiv_1);

    return 0;
}
