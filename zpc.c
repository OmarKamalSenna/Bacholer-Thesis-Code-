#include "f2c.h"
struct
{
	int nsevt, nevnt, nsbrun, ievt, isbrun;
} para3_;

#define para3_1 para3_

struct
{
	int mul;
} para1_;

#define para1_1 para1_

struct
{
	double xmp, xmu, alpha, rscut2, cutof2;
} para2_;

#define para2_1 para2_

struct
{
	int iftflg, ireflg, igeflg, ibstfg;
} para4_;

#define para4_1 para4_

struct
{
	int iconfg, iordsc;
} para5_;

#define para5_1 para5_

struct para6_1_
{
	double centy;
};

#define para6_1 (*(struct para6_1_ *)&para6_)

struct
{
	int ioscar;
} para7_;

#define para7_1 para7_

struct
{
	double gx0[400001], gy0[400001], gz0[400001], ft0[400001], px0[400001], py0[400001], pz0[400001], e0[400001], xmass0[400001];
	int ityp0[400001];
} prec1_;

#define prec1_1 prec1_

struct
{
	double gx[400001], gy[400001], gz[400001], ft[400001], px[400001], py[400001], pz[400001], e[400001], xmass[400001];
	int ityp[400001];
} prec2_;

#define prec2_1 prec2_

struct
{
	double gxs[400001], gys[400001], gzs[400001], fts[400001], pxs[400001], pys[400001], pzs[400001], es[400001], xmasss[400001];
	int ityps[400001];
} prec3_;

#define prec3_1 prec3_

struct
{
	double vx[400001], vy[400001], vz[400001];
} prec4_;

#define prec4_1 prec4_

struct
{
	double eta[400001], rap[400001], tau[400001];
} prec5_;

#define prec5_1 prec5_

struct
{
	double etas[400001], raps[400001], taus[400001];
} prec6_;

#define prec6_1 prec6_

struct
{
	int jxa, jya, jza;
} aurec1_;

#define aurec1_1 aurec1_

struct
{
	double dgxa[400001], dgya[400001], dgza[400001];
} aurec2_;

#define aurec2_1 aurec2_

struct
{
	int iscat, jscat, next[400001], last[400001], ictype, icsta[400001],
		nic[400001], icels[400001];
} ilist1_;

#define ilist1_1 ilist1_

struct
{
	int icell, icel[1000];
} ilist2_;

#define ilist2_1 ilist2_

struct
{
	double size1, size2, size3, v1, v2, v3, size;
} ilist3_;

#define ilist3_1 ilist3_

struct
{
	int ifmpt, ichkpt, indx[400001];
} ilist4_;

#define ilist4_1 ilist4_

struct
{
	double t;
	int iopern, icolln;
} ilist6_;

#define ilist6_1 ilist6_

struct
{
	int lstrg0[400001], lpart0[400001];
} ilist7_;

#define ilist7_1 ilist7_

struct
{
	int lstrg1[400001], lpart1[400001];
} ilist8_;

#define ilist8_1 ilist8_

struct rndm1_1_
{
	int number;
};

#define rndm1_1 (*(struct rndm1_1_ *)&rndm1_)

struct
{
	int iff;
} rndm2_;

#define rndm2_1 rndm2_

struct
{
	int iseedp;
} rndm3_;

#define rndm3_1 rndm3_

struct ana1_1_
{
	double ts[12];
};

#define ana1_1 (*(struct ana1_1_ *)&ana1_)

struct
{
	double det[12], dn[12], detdy[12], detdn[12], dndy[12], det1[12], dn1[12], detdy1[12], detdn1[12], dndy1[12], det2[12], dn2[12], detdy2[12], detdn2[12], dndy2[12];
} ana2_;

#define ana2_1 ana2_

struct
{
	double em[192];
} ana3_;

#define ana3_1 ana3_

struct
{
	double fdetdy[24], fdndy[24], fdndpt[12];
} ana4_;

#define ana4_1 ana4_

struct
{
	double enenew, pxnew, pynew, pznew;
} lor_;

#define lor_1 lor_

struct
{
	double ct[400001], ot[400001], tlarge;
} ilist5_;

#define ilist5_1 ilist5_

struct
{
	double formt;
} par1_;

#define par1_1 par1_

struct
{
	int nevent, isoft, isflag, izpc;
} anim_;

#define anim_1 anim_

struct
{
	double smearp, smearh;
} smearz_;

#define smearz_1 smearz_

struct
{
	double vxp0[400001], vyp0[400001], vzp0[400001];
} precpa_;

#define precpa_1 precpa_

struct
{
	double gxfrz[400001], gyfrz[400001], gzfrz[400001], ftfrz[400001],
		pxfrz[400001], pyfrz[400001], pzfrz[400001], efrz[400001], xmfrz[400001], tfrz[302];
	int ifrz[400001], idfrz[400001], itlast;
} frzprc_;

#define frzprc_1 frzprc_

union {
	struct
	{
		double xn1, xn2, xn3;
	} _1;
	struct
	{
		double vx3, vy3, vz3;
	} _2;
} cprod_;

#define cprod_1 (cprod_._1)
#define cprod_2 (cprod_._2)

struct
{
	int iaevt, iarun;
} arevt_;

#define arevt_1 arevt_

struct
{
	double e_1;
} para6_ = {0.};

struct
{
	int e_1;
} rndm1_ = {0};

struct
{
	double e_1[12];
} ana1_ = {.11, .12, .15, .2, .3, .4, .6, .8, 1., 2., 4., 6.};

static int c__9 = 9;
static int c__1 = 1;
static int c__3 = 3;
static int c__5 = 5;
static int c_b75 = 400001;
static int c__12 = 12;

int zpcmn_(void)
{

	int i__1, i__2;

	static int i__, j;
	extern int zpca1_(void), zpca2_(void), zpcou_(void),
		inievt_(void), inirun_(void), zpcrun_(void), zpstrg_(void);

	i__1 = para3_1.nevnt;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		para3_1.ievt = i__;

		inievt_();

		i__2 = para3_1.nsbrun;
		for (j = 1; j <= i__2; ++j)
		{
			para3_1.isbrun = j;

			inirun_();

		L3000:

			switch (zpcrun_())
			{
			case 1:
				goto L4000;
			}
			zpca1_();
			goto L3000;
		L4000:
			zpca2_();
		}
	}
	zpcou_();

	zpstrg_();
	return 0;
}

int zpcbdt_(void)
{
	return 0;
}

int inizpc_(void)
{
	extern int inian1_(void), readpa_(void), inipar_(void);

	readpa_();
	inipar_();
	inian1_();
	return 0;
}

int readpa_(void)
{

	int i__1;
	double d__1;
	olist o__1;
	cllist cl__1;

	int f_open(olist *), s_wsle(cilist *), do_lio(int *, int *, char *, ftnlen), e_wsle(void);
	int s_stop(char *, ftnlen);
	int f_clos(cllist *);

	static double a;
	static int i__;
	extern double ran1_(int *);
	static int iseed, iseed2, isedng, irused;

	static cilist io___4 = {0, 6, 0, 0, 0};
	static cilist io___5 = {0, 6, 0, 0, 0};

	iseed = rndm3_1.iseedp;

	o__1.oerr = 0;
	o__1.ounit = 25;
	o__1.ofnmlen = 11;
	o__1.ofnm = "ana/zpc.res";
	o__1.orl = 0;
	o__1.osta = "unknown";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);

	if (para7_1.ioscar == 1)
	{
		o__1.oerr = 0;
		o__1.ounit = 26;
		o__1.ofnmlen = 16;
		o__1.ofnm = "ana/parton.oscar";
		o__1.orl = 0;
		o__1.osta = "unknown";
		o__1.oacc = 0;
		o__1.ofm = 0;
		o__1.oblnk = 0;
		f_open(&o__1);
		o__1.oerr = 0;
		o__1.ounit = 19;
		o__1.ofnmlen = 16;
		o__1.ofnm = "ana/hadron.oscar";
		o__1.orl = 0;
		o__1.osta = "unknown";
		o__1.oacc = 0;
		o__1.ofm = 0;
		o__1.oblnk = 0;
		f_open(&o__1);
	}

	para2_1.xmp = 0.;

	d__1 = para2_1.alpha / para2_1.xmu;
	para2_1.cutof2 = d__1 * d__1 * 4.5;

	para2_1.rscut2 = .01;

	para3_1.nsevt = 1;

	para3_1.nevnt = 1;

	para3_1.nsbrun = 1;

	para4_1.iftflg = 0;

	para4_1.ireflg = 1;

	if (para4_1.ireflg == 0)
	{
		o__1.oerr = 0;
		o__1.ounit = 27;
		o__1.ofnmlen = 7;
		o__1.ofnm = "zpc.inp";
		o__1.orl = 0;
		o__1.osta = "UNKNOWN";
		o__1.oacc = 0;
		o__1.ofm = 0;
		o__1.oblnk = 0;
		f_open(&o__1);
	}

	para4_1.igeflg = 0;

	para4_1.ibstfg = 0;

	para5_1.iconfg = 1;

	para5_1.iordsc = 11;

	ilist3_1.v1 = .2;
	ilist3_1.v2 = .2;
	ilist3_1.v3 = .2;

	ilist3_1.size1 = 1.5;
	ilist3_1.size2 = 1.5;
	ilist3_1.size3 = .7;
	if (ilist3_1.size1 == 0. || ilist3_1.size2 == 0. || ilist3_1.size3 == 0.)
	{
		if (ilist3_1.size1 != 0. || ilist3_1.size2 != 0. || ilist3_1.size3 != 0. || ilist3_1.v1 != 0. || ilist3_1.v2 != 0. || ilist3_1.v3 != 0.)
		{
			s_wsle(&io___4);
			do_lio(&c__9, &c__1, "to get rid of space division:", (ftnlen)29);
			e_wsle();
			s_wsle(&io___5);
			do_lio(&c__9, &c__1, "set all sizes and vs to 0", (ftnlen)25);
			e_wsle();
			s_stop("chker", (ftnlen)5);
		}
	}

	d__1 = min(ilist3_1.size1, ilist3_1.size2);
	ilist3_1.size = min(d__1, ilist3_1.size3);

	rndm2_1.iff = -1;

	isedng = -iseed;

	a = ran1_(&isedng);

	irused = 2;
	i__1 = irused - 1;
	for (i__ = 1; i__ <= i__1; ++i__)
	{

		iseed2 = 2;
		a = ran1_(&iseed2);
	}

	if (para5_1.iconfg == 2 || para5_1.iconfg == 3)
	{
		ilist3_1.v1 = 0.;
		ilist3_1.v2 = 0.;
	}
	if (para5_1.iconfg == 4 || para5_1.iconfg == 5)
	{
		ilist3_1.v1 = 0.;
		ilist3_1.v2 = 0.;
		ilist3_1.v3 = 0.;
	}
	cl__1.cerr = 0;
	cl__1.cunit = 5;
	cl__1.csta = 0;
	f_clos(&cl__1);
	return 0;
}

int inipar_(void)
{

	if (para4_1.ibstfg != 0)
	{
		para6_1.centy = -6.;
	}
	return 0;
}

int inian1_(void)
{

	double cosh(double);

	static double a;
	static int i__;

	if (para4_1.ibstfg != 0)
	{
		a = cosh(6.);
		for (i__ = 1; i__ <= 12; ++i__)
		{
			ana1_1.ts[i__ - 1] *= a;
		}
	}
	return 0;
}

int inievt_(void)
{
	extern int readi_(void), genei_(void), boosti_(void);

	if (para4_1.ireflg == 0)
	{
		readi_();
	}
	if (para4_1.igeflg != 0)
	{
		genei_();
	}
	if (para4_1.ibstfg != 0)
	{
		boosti_();
	}
	return 0;
}

int readi_(void)
{

	int i__1;

	int s_rsle(cilist *), do_lio(int *, int *, char *, ftnlen),
		e_rsle(void);

	static int i__, neve, ntyp;
	static double field[9];

	static cilist io___16 = {0, 27, 1, 0, 0};

	for (i__ = 1; i__ <= 400001; ++i__)
	{
		if (para3_1.ievt != 1 && i__ == 1)
		{
			prec1_1.ityp0[i__ - 1] = ntyp;
			prec1_1.gx0[0] = field[0];
			prec1_1.gy0[0] = field[1];
			prec1_1.gz0[0] = field[2];
			prec1_1.ft0[0] = field[3];
			prec1_1.px0[0] = field[4];
			prec1_1.py0[0] = field[5];
			prec1_1.pz0[0] = field[6];
			prec1_1.e0[0] = field[7];
			prec1_1.xmass0[i__ - 1] = field[8];
			para1_1.mul = 1;
		}
		else
		{
		L900:
			i__1 = s_rsle(&io___16);
			if (i__1 != 0)
			{
				goto L1000;
			}
			i__1 = do_lio(&c__3, &c__1, (char *)&neve, (ftnlen)sizeof(int));
			if (i__1 != 0)
			{
				goto L1000;
			}
			i__1 = do_lio(&c__3, &c__1, (char *)&ntyp, (ftnlen)sizeof(int));
			if (i__1 != 0)
			{
				goto L1000;
			}
			i__1 = do_lio(&c__5, &c__9, (char *)&field[0], (ftnlen)sizeof(double));
			if (i__1 != 0)
			{
				goto L1000;
			}
			i__1 = e_rsle();
			if (i__1 != 0)
			{
				goto L1000;
			}
			if (neve < para3_1.nsevt)
			{
				goto L900;
			}
			if (neve > para3_1.nsevt + para3_1.ievt - 1)
			{
				goto L1000;
			}
			prec1_1.ityp0[i__ - 1] = ntyp;
			prec1_1.gx0[i__ - 1] = field[0];
			prec1_1.gy0[i__ - 1] = field[1];
			prec1_1.gz0[i__ - 1] = field[2];
			prec1_1.ft0[i__ - 1] = field[3];
			prec1_1.px0[i__ - 1] = field[4];
			prec1_1.py0[i__ - 1] = field[5];
			prec1_1.pz0[i__ - 1] = field[6];
			prec1_1.e0[i__ - 1] = field[7];
			prec1_1.xmass0[i__ - 1] = field[8];
			++para1_1.mul;
		}
	}
L1000:
	return 0;
}

int genei_(void)
{

	int i__1;
	double d__1, d__2;

	double sqrt(double), tanh(double), sinh(double), cosh(double);
	int s_wsle(cilist *), do_lio(int *, int *, char *, ftnlen),
		e_wsle(void);
	int s_stop(char *, ftnlen);

	static double e;
	static int i__;
	static double x, y, z__, r0, px, py, pz, bex, bey, bez;
	extern double ran1_(int *);
	static double tau0, deta, temp;
	static int iseed;
	extern int posit1_(double *, double *,
					   double *),
		posit2_(double *, double *), posit3_(double *, double *, double *);
	static double etamin, etamax;
	static int incmul;
	extern int energy_(double *, double *), momntm_(double *, double *, double *, double *), lorenz_(double *, double *, double *, double *, double *, double *, double *);

	static cilist io___37 = {0, 6, 0, 0, 0};

	iseed = rndm3_1.iseedp;
	incmul = 4000;
	temp = .5;
	etamin = -5.;
	etamax = 5.;
	r0 = 5.;
	tau0 = .1;
	deta = etamax - etamin;
	i__1 = para1_1.mul + incmul;
	for (i__ = para1_1.mul + 1; i__ <= i__1; ++i__)
	{
		prec1_1.ityp0[i__ - 1] = 21;
		prec1_1.xmass0[i__ - 1] = para2_1.xmp;
		energy_(&e, &temp);
		momntm_(&px, &py, &pz, &e);

		d__1 = e;

		d__2 = para2_1.xmp;
		e = sqrt(d__1 * d__1 + d__2 * d__2);
		if (para5_1.iconfg <= 3)
		{
			prec5_1.eta[i__ - 1] = etamin + deta * ran1_(&iseed);
			bex = 0.;
			bey = 0.;
			bez = -tanh(prec5_1.eta[i__ - 1]);
			lorenz_(&e, &px, &py, &pz, &bex, &bey, &bez);
			prec1_1.px0[i__ - 1] = lor_1.pxnew;
			prec1_1.py0[i__ - 1] = lor_1.pynew;
			prec1_1.pz0[i__ - 1] = lor_1.pznew;
			prec1_1.e0[i__ - 1] = lor_1.enenew;
		}
		else
		{
			prec1_1.px0[i__ - 1] = px;
			prec1_1.py0[i__ - 1] = py;
			prec1_1.pz0[i__ - 1] = pz;
			prec1_1.e0[i__ - 1] = e;
		}
	}
	i__1 = para1_1.mul + incmul;
	for (i__ = para1_1.mul + 1; i__ <= i__1; ++i__)
	{
		if (para5_1.iconfg <= 3)
		{
			prec1_1.gz0[i__ - 1] = tau0 * sinh(prec5_1.eta[i__ - 1]);
			prec1_1.ft0[i__ - 1] = tau0 * cosh(prec5_1.eta[i__ - 1]);
			if (para5_1.iconfg == 1)
			{
				posit1_(&x, &y, &r0);
				prec1_1.gx0[i__ - 1] = x + prec1_1.px0[i__ - 1] * prec1_1.ft0[i__ - 1] / prec1_1.e0[i__ - 1];
				prec1_1.gy0[i__ - 1] = y + prec1_1.py0[i__ - 1] * prec1_1.ft0[i__ - 1] / prec1_1.e0[i__ - 1];
			}
			else if (para5_1.iconfg == 2 || para5_1.iconfg == 3)
			{
				posit2_(&x, &y);
				prec1_1.gx0[i__ - 1] = x;
				prec1_1.gy0[i__ - 1] = y;
			}
		}
		else
		{
			prec1_1.ft0[i__ - 1] = 0.;
			posit3_(&x, &y, &z__);
			prec1_1.gx0[i__ - 1] = x;
			prec1_1.gy0[i__ - 1] = y;
			prec1_1.gz0[i__ - 1] = z__;
		}
	}
	para1_1.mul += incmul;

	if (para1_1.mul >= 400001 || para1_1.mul == 0)
	{
		s_wsle(&io___37);
		do_lio(&c__9, &c__1, "event", (ftnlen)5);
		do_lio(&c__3, &c__1, (char *)&para3_1.ievt, (ftnlen)sizeof(int));
		do_lio(&c__9, &c__1, "has", (ftnlen)3);
		do_lio(&c__3, &c__1, (char *)&para1_1.mul, (ftnlen)sizeof(int));
		do_lio(&c__9, &c__1, "number of gluon", (ftnlen)15);
		do_lio(&c__9, &c__1, "adjusting counting is necessary", (ftnlen)31);
		e_wsle();
		s_stop("adarr", (ftnlen)5);
	}
	return 0;
}

int posit1_(double *x, double *y, double *r0)
{

	double d__1, d__2;

	extern double ran1_(int *);
	static int iseed;

	iseed = rndm3_1.iseedp;
L10:
	*x = ran1_(&iseed) * 2. - 1.;
	*y = ran1_(&iseed) * 2. - 1.;

	d__1 = *x;

	d__2 = *y;
	if (d__1 * d__1 + d__2 * d__2 > 1.)
	{
		goto L10;
	}
	*x *= *r0;
	*y *= *r0;
	return 0;
}

int posit2_(double *x, double *y)
{
	extern double ran1_(int *);
	static int iseed;

	iseed = rndm3_1.iseedp;
	*x = ran1_(&iseed) * 2. - 1.;
	*y = ran1_(&iseed) * 2. - 1.;
	*x = *x * 5. * ilist3_1.size1;
	*y = *y * 5. * ilist3_1.size2;
	return 0;
}

int posit3_(double *x, double *y, double *z__)
{
	extern double ran1_(int *);
	static int iseed;

	iseed = rndm3_1.iseedp;
	*x = ran1_(&iseed) * 2. - 1.;
	*y = ran1_(&iseed) * 2. - 1.;
	*z__ = ran1_(&iseed) * 2. - 1.;
	*x = *x * 5. * ilist3_1.size1;
	*y = *y * 5. * ilist3_1.size2;
	*z__ = *z__ * 5. * ilist3_1.size3;
	return 0;
}

int energy_(double *e, double *temp)
{

	double d__1, d__2;

	double log(double), sqrt(double), exp(double);

	extern double ran1_(int *);
	static int iseed;

	iseed = rndm3_1.iseedp;
L1000:
	*e = ran1_(&iseed);
	*e *= ran1_(&iseed);
	*e *= ran1_(&iseed);
	if (*e <= 0.)
	{
		goto L1000;
	}
	*e = -(*temp) * log(*e);

	d__1 = *e;

	d__2 = para2_1.xmp;
	if (ran1_(&iseed) > exp((*e - sqrt(d__1 * d__1 + d__2 * d__2)) / *temp))
	{
		goto L1000;
	}
	return 0;
}

int momntm_(double *px, double *py, double *pz,
			double *e)
{

	double d__1;

	double sqrt(double), cos(double), sin(double);

	static double phi;
	extern double ran1_(int *);
	static double cost, sint;
	static int iseed;

	iseed = rndm3_1.iseedp;
	cost = ran1_(&iseed) * 2. - 1.;

	d__1 = cost;
	sint = sqrt(1. - d__1 * d__1);
	phi = ran1_(&iseed) * 6.28318530717958;
	*px = *e * sint * cos(phi);
	*py = *e * sint * sin(phi);
	*pz = *e * cost;
	return 0;
}

int boosti_(void)
{

	int i__1;

	double tanh(double);

	static int i__;
	static double e1, px1, py1, pz1, bex, bey, bez;
	extern int lorenz_(double *, double *,
					   double *, double *, double *, double *,
					   double *);

	bex = 0.;
	bey = 0.;
	bez = -tanh(para6_1.centy);

	i__1 = para1_1.mul;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		px1 = prec1_1.gx0[i__ - 1];
		py1 = prec1_1.gy0[i__ - 1];
		pz1 = prec1_1.gz0[i__ - 1];
		e1 = prec1_1.ft0[i__ - 1];
		lorenz_(&e1, &px1, &py1, &pz1, &bex, &bey, &bez);
		prec1_1.gx0[i__ - 1] = lor_1.pxnew;
		prec1_1.gy0[i__ - 1] = lor_1.pynew;
		prec1_1.gz0[i__ - 1] = lor_1.pznew;
		prec1_1.ft0[i__ - 1] = lor_1.enenew;
		px1 = prec1_1.px0[i__ - 1];
		py1 = prec1_1.py0[i__ - 1];
		pz1 = prec1_1.pz0[i__ - 1];
		e1 = prec1_1.e0[i__ - 1];
		lorenz_(&e1, &px1, &py1, &pz1, &bex, &bey, &bez);
		prec1_1.px0[i__ - 1] = lor_1.pxnew;
		prec1_1.py0[i__ - 1] = lor_1.pynew;
		prec1_1.pz0[i__ - 1] = lor_1.pznew;
		prec1_1.e0[i__ - 1] = lor_1.enenew;
	}
	return 0;
}

int inirun_(void)
{
	extern int ftime_(void), inian2_(void), inirec_(void),
		iilist_(void);

	ftime_();
	inirec_();
	iilist_();
	inian2_();
	return 0;
}

int ftime_(void)
{

	int i__1;
	double d__1, d__2, d__3;

	static int i__;
	static double xmt2;
	static int iseed;
	extern double ftime1_(int *);
	extern int index1_(int *, int *, double *,
					   int *);

	iseed = rndm3_1.iseedp;

	for (i__ = 1; i__ <= 400001; ++i__)
	{
		ilist5_1.ct[i__ - 1] = 0.;
		ilist5_1.ot[i__ - 1] = 0.;
	}
	ilist5_1.tlarge = 1e6;

	if (para4_1.iftflg == 0)
	{

		if (anim_1.isoft == 3 || anim_1.isoft == 4 || anim_1.isoft == 5)
		{
			i__1 = para1_1.mul;
			for (i__ = 1; i__ <= i__1; ++i__)
			{
				if (prec1_1.ft0[i__ - 1] > ilist5_1.tlarge)
				{
					prec1_1.ft0[i__ - 1] = ilist5_1.tlarge;
				}
			}
			goto L150;
		}
		else
		{

			for (i__ = 1; i__ <= 400001; ++i__)
			{
				prec1_1.ft0[i__ - 1] = ilist5_1.tlarge;
			}
			i__1 = para1_1.mul;
			for (i__ = 1; i__ <= i__1; ++i__)
			{

				d__1 = prec1_1.px0[i__ - 1];

				d__2 = prec1_1.py0[i__ - 1];

				d__3 = para2_1.xmp;
				xmt2 = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
				par1_1.formt = xmt2 / prec1_1.e0[i__ - 1];
				prec1_1.ft0[i__ - 1] = ftime1_(&iseed);
				if (prec1_1.ft0[i__ - 1] > ilist5_1.tlarge)
				{
					prec1_1.ft0[i__ - 1] = ilist5_1.tlarge;
				}
			}
		}
	}

L150:

	if (para1_1.mul > 1)
	{
		index1_(&c_b75, &para1_1.mul, prec1_1.ft0, ilist4_1.indx);
	}
	else
	{

		ilist4_1.indx[0] = 1;
	}

	return 0;
}

int inirec_(void)
{

	int i__1;
	double d__1;

	double log(double), cosh(double);

	static int i__;
	static double vxp[400001], vyp[400001], vzp[400001];
	static int iseed, indxi;
	static double formt, energy;
	extern int inifrz_(void);

	iseed = rndm3_1.iseedp;

	if (anim_1.isoft == 5)
	{
		frzprc_1.itlast = 0;
		inifrz_();
	}
	i__1 = para1_1.mul;
	for (i__ = 1; i__ <= i__1; ++i__)
	{

		indxi = ilist4_1.indx[i__ - 1];
		prec2_1.ityp[i__ - 1] = prec1_1.ityp0[indxi - 1];
		prec2_1.gx[i__ - 1] = prec1_1.gx0[indxi - 1];
		prec2_1.gy[i__ - 1] = prec1_1.gy0[indxi - 1];
		prec2_1.gz[i__ - 1] = prec1_1.gz0[indxi - 1];
		prec2_1.ft[i__ - 1] = prec1_1.ft0[indxi - 1];
		prec2_1.px[i__ - 1] = prec1_1.px0[indxi - 1];
		prec2_1.py[i__ - 1] = prec1_1.py0[indxi - 1];
		prec2_1.pz[i__ - 1] = prec1_1.pz0[indxi - 1];
		prec2_1.e[i__ - 1] = prec1_1.e0[indxi - 1];
		prec2_1.xmass[i__ - 1] = prec1_1.xmass0[indxi - 1];
		ilist8_1.lstrg1[i__ - 1] = ilist7_1.lstrg0[indxi - 1];
		ilist8_1.lpart1[i__ - 1] = ilist7_1.lpart0[indxi - 1];
		vxp[i__ - 1] = precpa_1.vxp0[indxi - 1];
		vyp[i__ - 1] = precpa_1.vyp0[indxi - 1];
		vzp[i__ - 1] = precpa_1.vzp0[indxi - 1];

		if (anim_1.isoft == 5)
		{
			frzprc_1.idfrz[i__ - 1] = prec2_1.ityp[i__ - 1];
			frzprc_1.gxfrz[i__ - 1] = prec2_1.gx[i__ - 1];
			frzprc_1.gyfrz[i__ - 1] = prec2_1.gy[i__ - 1];
			frzprc_1.gzfrz[i__ - 1] = prec2_1.gz[i__ - 1];
			frzprc_1.ftfrz[i__ - 1] = prec2_1.ft[i__ - 1];
			frzprc_1.pxfrz[i__ - 1] = prec2_1.px[i__ - 1];
			frzprc_1.pyfrz[i__ - 1] = prec2_1.py[i__ - 1];
			frzprc_1.pzfrz[i__ - 1] = prec2_1.pz[i__ - 1];
			frzprc_1.efrz[i__ - 1] = prec2_1.e[i__ - 1];
			frzprc_1.xmfrz[i__ - 1] = prec2_1.xmass[i__ - 1];
			frzprc_1.ifrz[i__ - 1] = 0;
		}
	}

	i__1 = para1_1.mul;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		prec3_1.ityps[i__ - 1] = prec2_1.ityp[i__ - 1];
		prec3_1.gxs[i__ - 1] = prec2_1.gx[i__ - 1];
		prec3_1.gys[i__ - 1] = prec2_1.gy[i__ - 1];
		prec3_1.gzs[i__ - 1] = prec2_1.gz[i__ - 1];
		prec3_1.fts[i__ - 1] = prec2_1.ft[i__ - 1];
		prec3_1.pxs[i__ - 1] = prec2_1.px[i__ - 1];
		prec3_1.pys[i__ - 1] = prec2_1.py[i__ - 1];
		prec3_1.pzs[i__ - 1] = prec2_1.pz[i__ - 1];
		prec3_1.es[i__ - 1] = prec2_1.e[i__ - 1];
		prec3_1.xmasss[i__ - 1] = prec2_1.xmass[i__ - 1];
	}
	i__1 = para1_1.mul;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		energy = prec2_1.e[i__ - 1];
		prec4_1.vx[i__ - 1] = prec2_1.px[i__ - 1] / energy;
		prec4_1.vy[i__ - 1] = prec2_1.py[i__ - 1] / energy;
		prec4_1.vz[i__ - 1] = prec2_1.pz[i__ - 1] / energy;
		if (para4_1.iftflg == 0)
		{
			formt = prec2_1.ft[i__ - 1];

			if (anim_1.isoft == 3 || anim_1.isoft == 4 || anim_1.isoft == 5)
			{
				prec2_1.gx[i__ - 1] += vxp[i__ - 1] * formt;
				prec2_1.gy[i__ - 1] += vyp[i__ - 1] * formt;
				prec2_1.gz[i__ - 1] += vzp[i__ - 1] * formt;
			}
			else
			{
				prec2_1.gx[i__ - 1] += prec4_1.vx[i__ - 1] * formt;
				prec2_1.gy[i__ - 1] += prec4_1.vy[i__ - 1] * formt;
				prec2_1.gz[i__ - 1] += prec4_1.vz[i__ - 1] * formt;
			}
		}
	}
	if (para5_1.iconfg <= 3)
	{
		i__1 = para1_1.mul;
		for (i__ = 1; i__ <= i__1; ++i__)
		{
			if (prec2_1.ft[i__ - 1] <= (d__1 = prec2_1.gz[i__ - 1], abs(d__1)))
			{
				prec5_1.eta[i__ - 1] = 1e6;
			}
			else
			{
				prec5_1.eta[i__ - 1] = log((prec2_1.ft[i__ - 1] + prec2_1.gz[i__ - 1]) / (prec2_1.ft[i__ - 1] - prec2_1.gz[i__ - 1])) * .5;
			}
			if (prec2_1.e[i__ - 1] <= (d__1 = prec2_1.pz[i__ - 1], abs(d__1)))
			{
				prec5_1.rap[i__ - 1] = 1e6;
			}
			else
			{
				prec5_1.rap[i__ - 1] = log((prec2_1.e[i__ - 1] + prec2_1.pz[i__ - 1]) / (prec2_1.e[i__ - 1] - prec2_1.pz[i__ - 1])) * .5;
			}
			prec5_1.tau[i__ - 1] = prec2_1.ft[i__ - 1] / cosh(prec5_1.eta[i__ - 1]);
		}
		i__1 = para1_1.mul;
		for (i__ = 1; i__ <= i__1; ++i__)
		{
			prec6_1.etas[i__ - 1] = prec5_1.eta[i__ - 1];
			prec6_1.raps[i__ - 1] = prec5_1.rap[i__ - 1];
			prec6_1.taus[i__ - 1] = prec5_1.tau[i__ - 1];
		}
	}
	return 0;
}

int iilist_(void)
{

	int i__1;

	static int i__, i1, i2, i3;

	ilist1_1.iscat = 400001;
	ilist1_1.jscat = 400001;
	i__1 = para1_1.mul;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		ilist1_1.next[i__ - 1] = 0;
		ilist1_1.last[i__ - 1] = 0;
		ilist1_1.icsta[i__ - 1] = 0;
		ilist1_1.nic[i__ - 1] = 0;
		ilist1_1.icels[i__ - 1] = 0;
	}
	ilist2_1.icell = 0;
	for (i1 = 1; i1 <= 10; ++i1)
	{
		for (i2 = 1; i2 <= 10; ++i2)
		{
			for (i3 = 1; i3 <= 10; ++i3)
			{
				ilist2_1.icel[i1 + (i2 + i3 * 10) * 10 - 111] = 0;
			}
		}
	}
	ilist4_1.ichkpt = 0;
	ilist4_1.ifmpt = 1;
	i__1 = para1_1.mul;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		ilist5_1.ct[i__ - 1] = ilist5_1.tlarge;
		ilist5_1.ot[i__ - 1] = ilist5_1.tlarge;
	}
	ilist6_1.iopern = 0;
	ilist6_1.icolln = 0;
	ilist6_1.t = 0.;
	return 0;
}

int inian2_(void)
{
	static int i__;

	if (para5_1.iconfg <= 3)
	{
		for (i__ = 1; i__ <= 12; ++i__)
		{
			ana2_1.det[i__ - 1] = 0.;
			ana2_1.dn[i__ - 1] = 0.;
			ana2_1.det1[i__ - 1] = 0.;
			ana2_1.dn1[i__ - 1] = 0.;
			ana2_1.det2[i__ - 1] = 0.;
			ana2_1.dn2[i__ - 1] = 0.;
		}
	}
	return 0;
}

int zpcrun_(void)
{

	int s_wsle(cilist *), do_lio(int *, int *, char *, ftnlen),
		e_wsle(void);
	int s_stop(char *, ftnlen);

	static double t1;
	extern int scat_(double *, int *, int *),
		local_(double *), ulist_(double *);
	static int iscat0, jscat0;
	extern int celasn_(void), cellre_(int *, double *), getict_(double *);
	static int niscat;
	extern int savrec_(int *);
	static int njscat;

	static cilist io___73 = {0, 6, 0, 0, 0};

	if (ilist1_1.ictype % 2 == 0)
	{
		savrec_(&ilist1_1.iscat);
		savrec_(&ilist1_1.jscat);
	}

	getict_(&t1);

	if (para5_1.iconfg == 1 && t1 > ilist5_1.tlarge / 2.)
	{
		return 1;
	}
	if (para5_1.iconfg == 2 || para5_1.iconfg == 3)
	{
		if (t1 > 250.)
		{
			return 1;
		}
	}
	if (para5_1.iconfg == 4 || para5_1.iconfg == 5)
	{
		if (t1 > 6.1)
		{
			return 1;
		}
	}

	if (anim_1.isoft == 5)
	{
		local_(&t1);
	}

	++ilist6_1.iopern;
	ilist6_1.t = t1;
	if (ilist1_1.ictype % 2 == 0)
	{
		++ilist6_1.icolln;
	}

	if (para5_1.iconfg == 1 || para5_1.iconfg == 2 || para5_1.iconfg == 4)
	{
		if (ilist1_1.ictype == 1 || ilist1_1.ictype == 2 || ilist1_1.ictype == 5 || ilist1_1.ictype == 6)
		{
			celasn_();
		}
	}

	if (ilist1_1.ictype != 1)
	{
		iscat0 = ilist1_1.iscat;
		jscat0 = ilist1_1.jscat;

		ilist1_1.iscat = max(iscat0, jscat0);
		ilist1_1.jscat = min(iscat0, jscat0);

		if (ilist1_1.jscat != 0)
		{
			if (ilist1_1.next[ilist1_1.jscat - 1] != ilist1_1.iscat)
			{
				s_wsle(&io___73);
				do_lio(&c__9, &c__1, "iscat=", (ftnlen)6);
				do_lio(&c__3, &c__1, (char *)&ilist1_1.iscat, (ftnlen)sizeof(int));
				do_lio(&c__9, &c__1, "jscat=", (ftnlen)6);
				do_lio(&c__3, &c__1, (char *)&ilist1_1.jscat, (ftnlen)sizeof(int));
				do_lio(&c__9, &c__1, "next(", (ftnlen)5);
				do_lio(&c__3, &c__1, (char *)&ilist1_1.jscat, (ftnlen)sizeof(int));
				do_lio(&c__9, &c__1, ")=", (ftnlen)2);
				do_lio(&c__3, &c__1, (char *)&ilist1_1.next[ilist1_1.jscat - 1], (ftnlen)sizeof(int));
				e_wsle();
				if (ilist5_1.ct[ilist1_1.iscat - 1] < ilist5_1.tlarge / 2.)
				{
					s_stop("tterr", (ftnlen)5);
				}
				if (ilist5_1.ct[ilist1_1.jscat - 1] < ilist5_1.tlarge / 2.)
				{
					s_stop("tterr", (ftnlen)5);
				}
			}
		}

		niscat = ilist1_1.iscat;
		njscat = ilist1_1.jscat;

		if (ilist1_1.icsta[ilist1_1.iscat - 1] != 0)
		{
			cellre_(&niscat, &ilist6_1.t);
		}
		if (ilist1_1.jscat != 0)
		{
			if (ilist1_1.icsta[ilist1_1.jscat - 1] != 0)
			{
				cellre_(&njscat, &ilist6_1.t);
			}
		}

		if (ilist1_1.ictype % 2 == 0)
		{
			scat_(&ilist6_1.t, &ilist1_1.iscat, &ilist1_1.jscat);
		}
	}

	ulist_(&ilist6_1.t);

	if (ilist4_1.ifmpt <= para1_1.mul)
	{
		if (ilist1_1.ictype != 0 && ilist1_1.ictype != 3 && ilist1_1.ictype != 4)
		{
			++ilist4_1.ichkpt;
			++ilist4_1.ifmpt;
		}
	}
	return 0;
}

int savrec_(int *i__)
{

	prec3_1.ityps[*i__ - 1] = prec2_1.ityp[*i__ - 1];
	prec3_1.gxs[*i__ - 1] = prec2_1.gx[*i__ - 1];
	prec3_1.gys[*i__ - 1] = prec2_1.gy[*i__ - 1];
	prec3_1.gzs[*i__ - 1] = prec2_1.gz[*i__ - 1];
	prec3_1.fts[*i__ - 1] = prec2_1.ft[*i__ - 1];
	prec3_1.pxs[*i__ - 1] = prec2_1.px[*i__ - 1];
	prec3_1.pys[*i__ - 1] = prec2_1.py[*i__ - 1];
	prec3_1.pzs[*i__ - 1] = prec2_1.pz[*i__ - 1];
	prec3_1.es[*i__ - 1] = prec2_1.e[*i__ - 1];
	prec3_1.xmasss[*i__ - 1] = prec2_1.xmass[*i__ - 1];
	prec6_1.etas[*i__ - 1] = prec5_1.eta[*i__ - 1];
	prec6_1.raps[*i__ - 1] = prec5_1.rap[*i__ - 1];
	prec6_1.taus[*i__ - 1] = prec5_1.tau[*i__ - 1];
	return 0;
}

int getict_(double *t1)
{

	int i__1;

	static int i__;

	*t1 = ilist5_1.tlarge;
	ilist1_1.iscat = 0;
	ilist1_1.jscat = 0;

	i__1 = ilist4_1.ichkpt;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		if (ilist5_1.ot[i__ - 1] < *t1)
		{
			*t1 = ilist5_1.ot[i__ - 1];
			ilist1_1.iscat = i__;
		}
	}
	if (ilist1_1.iscat != 0)
	{
		ilist1_1.jscat = ilist1_1.next[ilist1_1.iscat - 1];
	}

	if (ilist1_1.iscat != 0 && ilist1_1.jscat != 0)
	{
		if (ilist1_1.icsta[ilist1_1.iscat - 1] == 0 && ilist1_1.icsta[ilist1_1.jscat - 1] == 0)
		{
			ilist1_1.ictype = 0;
		}
		else
		{
			ilist1_1.ictype = 4;
		}
	}
	else if (ilist1_1.iscat != 0 || ilist1_1.jscat != 0)
	{
		ilist1_1.ictype = 3;
	}

	if (ilist4_1.ifmpt <= para1_1.mul)
	{
		if (prec2_1.ft[ilist4_1.ifmpt - 1] < *t1)
		{
			ilist1_1.ictype = 1;
			*t1 = prec2_1.ft[ilist4_1.ifmpt - 1];
		}
		else if (prec2_1.ft[ilist4_1.ifmpt - 1] == *t1)
		{
			if (ilist1_1.ictype == 0)
			{
				ilist1_1.ictype = 2;
			}
			if (ilist1_1.ictype == 3)
			{
				ilist1_1.ictype = 5;
			}
			if (ilist1_1.ictype == 4)
			{
				ilist1_1.ictype = 6;
			}
		}
	}
	return 0;
}

int celasn_(void)
{

	double d__1;

	static int i__, j, i1, i2, i3;
	static double td, tt;
	extern int integ_(double *);
	extern int newcre_(int *, int *);

	i__ = ilist4_1.ifmpt;
	tt = prec2_1.ft[i__ - 1];
	td = tt - ilist3_1.size;
	if (para5_1.iconfg == 1 && (ilist3_1.size1 == 0. || ilist3_1.size2 == 0. || ilist3_1.size3 == 0.))
	{
		i1 = 11;
		i2 = 11;
		i3 = 11;
	}
	else if (para5_1.iconfg == 4 || td <= 0.)
	{
		d__1 = prec2_1.gx[i__ - 1] / ilist3_1.size1;
		i1 = integ_(&d__1) + 6;
		d__1 = prec2_1.gy[i__ - 1] / ilist3_1.size2;
		i2 = integ_(&d__1) + 6;
		d__1 = prec2_1.gz[i__ - 1] / ilist3_1.size3;
		i3 = integ_(&d__1) + 6;
		d__1 = prec2_1.gx[i__ - 1] / ilist3_1.size1;
		if ((double)integ_(&d__1) == prec2_1.gx[i__ - 1] /
											 ilist3_1.size1 &&
			prec4_1.vx[i__ - 1] < 0.)
		{
			--i1;
		}
		d__1 = prec2_1.gy[i__ - 1] / ilist3_1.size2;
		if ((double)integ_(&d__1) == prec2_1.gy[i__ - 1] /
											 ilist3_1.size2 &&
			prec4_1.vy[i__ - 1] < 0.)
		{
			--i2;
		}
		d__1 = prec2_1.gz[i__ - 1] / ilist3_1.size3;
		if ((double)integ_(&d__1) == prec2_1.gz[i__ - 1] /
											 ilist3_1.size3 &&
			prec4_1.vz[i__ - 1] < 0.)
		{
			--i3;
		}
	}
	else
	{
		d__1 = prec2_1.gx[i__ - 1] / (ilist3_1.size1 + ilist3_1.v1 * td);
		i1 = integ_(&d__1) + 6;
		d__1 = prec2_1.gy[i__ - 1] / (ilist3_1.size2 + ilist3_1.v2 * td);
		i2 = integ_(&d__1) + 6;
		d__1 = prec2_1.gz[i__ - 1] / (ilist3_1.size3 + ilist3_1.v3 * td);
		i3 = integ_(&d__1) + 6;
		d__1 = prec2_1.gx[i__ - 1] / (ilist3_1.size1 + ilist3_1.v1 * td);
		if ((double)integ_(&d__1) == prec2_1.gx[i__ - 1] / (ilist3_1.size1 + ilist3_1.v1 * td) && prec4_1.vx[i__ - 1] < (i1 - 6) * ilist3_1.v1)
		{
			--i1;
		}
		d__1 = prec2_1.gy[i__ - 1] / (ilist3_1.size2 + ilist3_1.v2 * td);
		if ((double)integ_(&d__1) == prec2_1.gy[i__ - 1] / (ilist3_1.size2 + ilist3_1.v2 * td) && prec4_1.vy[i__ - 1] < (i2 - 6) * ilist3_1.v2)
		{
			--i2;
		}
		d__1 = prec2_1.gz[i__ - 1] / (ilist3_1.size3 + ilist3_1.v3 * td);
		if ((double)integ_(&d__1) == prec2_1.gz[i__ - 1] / (ilist3_1.size3 + ilist3_1.v3 * td) && prec4_1.vz[i__ - 1] < (i3 - 6) * ilist3_1.v3)
		{
			--i3;
		}
	}
	if (i1 <= 0 || i1 >= 11 || i2 <= 0 || i2 >= 11 || i3 <= 0 || i3 >= 11)
	{
		i1 = 11;
		i2 = 11;
		i3 = 11;
	}
	if (i1 == 11)
	{
		j = ilist2_1.icell;
		newcre_(&i__, &j);
		ilist2_1.icell = j;
		ilist1_1.icels[i__ - 1] = 111111;
	}
	else
	{
		j = ilist2_1.icel[i1 + (i2 + i3 * 10) * 10 - 111];
		newcre_(&i__, &j);
		ilist2_1.icel[i1 + (i2 + i3 * 10) * 10 - 111] = j;
		ilist1_1.icels[i__ - 1] = i1 * 10000 + i2 * 100 + i3;
	}
	return 0;
}

int integ_(double *x)
{

	int ret_val;

	if (*x < 0.)
	{
		ret_val = (int)(*x - 1.);
	}
	else
	{
		ret_val = (int)(*x);
	}
	return ret_val;
}

int cellre_(int *i__, double *t)
{

	int i__1;
	double d__1;

	static int j, k, i1, i2, i3;
	static double t0;
	static int ii;
	static double ddt, dtt;
	static logical good;
	static double ctmp, otmp, tmin1;
	extern int wallc_(int *, int *, int *,
					  int *, double *, double *);
	extern int integ_(double *);
	extern int dchin1_(int *, int *, int *,
					   int *, int *, double *),
		dchin2_(int *, int *,
				int *, int *, int *, double *),
		dchin3_(int *, int *, int *, int *, int *, double *);
	static int icels0;
	extern int oldcre_(int *), newcre_(int *, int *), dchout_(int *, int *, double *);

	t0 = *t;
L1000:
	if (para5_1.iconfg == 3 || para5_1.iconfg == 5)
	{
		k = ilist1_1.icsta[*i__ - 1] % 10;
		if (k == 1)
		{
			prec2_1.gx[*i__ - 1] -= ilist3_1.size1 * 10.;
			aurec2_1.dgxa[*i__ - 1] += ilist3_1.size1 * 10.;
			i__1 = ilist4_1.ichkpt;
			for (ii = 1; ii <= i__1; ++ii)
			{
				if (ilist1_1.next[ii - 1] == *i__)
				{
					aurec2_1.dgxa[ii - 1] -= ilist3_1.size1 * 10.;
				}
			}
		}
		if (k == 2)
		{
			prec2_1.gx[*i__ - 1] += ilist3_1.size1 * 10.;
			aurec2_1.dgxa[*i__ - 1] -= ilist3_1.size1 * 10.;
			i__1 = ilist4_1.ichkpt;
			for (ii = 1; ii <= i__1; ++ii)
			{
				if (ilist1_1.next[ii - 1] == *i__)
				{
					aurec2_1.dgxa[ii - 1] += ilist3_1.size1 * 10.;
				}
			}
		}
		if (k == 3)
		{
			prec2_1.gy[*i__ - 1] -= ilist3_1.size2 * 10.;
			aurec2_1.dgya[*i__ - 1] += ilist3_1.size2 * 10.;
			i__1 = ilist4_1.ichkpt;
			for (ii = 1; ii <= i__1; ++ii)
			{
				if (ilist1_1.next[ii - 1] == *i__)
				{
					aurec2_1.dgya[ii - 1] -= ilist3_1.size2 * 10.;
				}
			}
		}
		if (k == 4)
		{
			prec2_1.gy[*i__ - 1] += ilist3_1.size2 * 10.;
			aurec2_1.dgya[*i__ - 1] -= ilist3_1.size2 * 10.;
			i__1 = ilist4_1.ichkpt;
			for (ii = 1; ii <= i__1; ++ii)
			{
				if (ilist1_1.next[ii - 1] == *i__)
				{
					aurec2_1.dgya[ii - 1] += ilist3_1.size2 * 10.;
				}
			}
		}
		if (para5_1.iconfg == 5)
		{
			if (k == 5)
			{
				prec2_1.gz[*i__ - 1] -= ilist3_1.size3 * 10.;
				aurec2_1.dgza[*i__ - 1] += ilist3_1.size3 * 10.;
				i__1 = ilist4_1.ichkpt;
				for (ii = 1; ii <= i__1; ++ii)
				{
					if (ilist1_1.next[ii - 1] == *i__)
					{
						aurec2_1.dgza[ii - 1] -= ilist3_1.size3 * 10.;
					}
				}
			}
			if (k == 6)
			{
				prec2_1.gz[*i__ - 1] += ilist3_1.size3 * 10.;
				aurec2_1.dgza[*i__ - 1] -= ilist3_1.size3 * 10.;
				i__1 = ilist4_1.ichkpt;
				for (ii = 1; ii <= i__1; ++ii)
				{
					if (ilist1_1.next[ii - 1] == *i__)
					{
						aurec2_1.dgza[ii - 1] += ilist3_1.size3 * 10.;
					}
				}
			}
		}
	}
	else
	{
		icels0 = ilist1_1.icels[*i__ - 1];
		i1 = icels0 / 10000;
		i2 = (icels0 - i1 * 10000) / 100;
		i3 = icels0 - i1 * 10000 - i2 * 100;

		if (i1 >= 1 && i1 <= 10 && i2 >= 1 && i2 <= 10 && i3 >= 1 && i3 <= 10)
		{

			if (ilist2_1.icel[i1 + (i2 + i3 * 10) * 10 - 111] == *i__)
			{
				ilist2_1.icel[i1 + (i2 + i3 * 10) * 10 - 111] = ilist1_1.nic[*i__ - 1];
			}

			oldcre_(i__);

			k = ilist1_1.icsta[*i__ - 1] % 10;

			if (para5_1.iconfg == 1)
			{
				good = i1 == 1 && k == 2 || i1 == 10 && k == 1 || i2 == 1 && k == 4 || i2 == 10 && k == 3 || i3 == 1 && k == 6 ||
					   i3 == 10 && k == 5;
			}
			if (para5_1.iconfg == 2)
			{
				good = i3 == 1 && k == 6 || i3 == 10 && k == 5;
			}
			if (good)
			{

				newcre_(i__, &ilist2_1.icell);

				ilist1_1.icels[*i__ - 1] = 111111;
			}
			else
			{
				if (k == 1)
				{
					++i1;
				}
				if (k == 2)
				{
					--i1;
				}
				if (k == 3)
				{
					++i2;
				}
				if (k == 4)
				{
					--i2;
				}
				if (k == 5)
				{
					++i3;
				}
				if (k == 6)
				{
					--i3;
				}
				if (para5_1.iconfg == 2 || para5_1.iconfg == 4)
				{
					if (i1 == 0)
					{
						i1 = 10;
						prec2_1.gx[*i__ - 1] += ilist3_1.size1 * 10.;
					}
					if (i1 == 11)
					{
						i1 = 1;
						prec2_1.gx[*i__ - 1] -= ilist3_1.size1 * 10.;
					}
					if (i2 == 0)
					{
						i2 = 10;
						prec2_1.gy[*i__ - 1] += ilist3_1.size2 * 10.;
					}
					if (i2 == 11)
					{
						i2 = 1;
						prec2_1.gy[*i__ - 1] -= ilist3_1.size2 * 10.;
					}
					if (para5_1.iconfg == 4)
					{
						if (i3 == 0)
						{
							i3 = 10;
							prec2_1.gz[*i__ - 1] += ilist3_1.size3 * 10.;
						}
						if (i3 == 11)
						{
							i3 = 1;
							prec2_1.gz[*i__ - 1] -= ilist3_1.size3 * 10.;
						}
					}
				}
				j = ilist2_1.icel[i1 + (i2 + i3 * 10) * 10 - 111];
				newcre_(i__, &j);

				ilist2_1.icel[i1 + (i2 + i3 * 10) * 10 - 111] = j;
				ilist1_1.icels[*i__ - 1] = i1 * 10000 + i2 * 100 + i3;
			}
		}
		else
		{
			if (ilist2_1.icell == *i__)
			{
				ilist2_1.icell = ilist1_1.nic[*i__ - 1];
			}
			oldcre_(i__);
			k = ilist1_1.icsta[*i__ - 1] % 10;
			ddt = *t - prec2_1.ft[*i__ - 1];
			dtt = *t - ilist3_1.size;
			if (dtt <= 0.)
			{
				d__1 = (prec2_1.gx[*i__ - 1] + prec4_1.vx[*i__ - 1] * ddt) /
					   ilist3_1.size1;
				i1 = integ_(&d__1) + 6;
				d__1 = (prec2_1.gy[*i__ - 1] + prec4_1.vy[*i__ - 1] * ddt) /
					   ilist3_1.size2;
				i2 = integ_(&d__1) + 6;
				d__1 = (prec2_1.gz[*i__ - 1] + prec4_1.vz[*i__ - 1] * ddt) /
					   ilist3_1.size3;
				i3 = integ_(&d__1) + 6;
			}
			else
			{
				d__1 = (prec2_1.gx[*i__ - 1] + prec4_1.vx[*i__ - 1] * ddt) / (ilist3_1.size1 + ilist3_1.v1 * dtt);
				i1 = integ_(&d__1) + 6;
				d__1 = (prec2_1.gy[*i__ - 1] + prec4_1.vy[*i__ - 1] * ddt) / (ilist3_1.size2 + ilist3_1.v2 * dtt);
				i2 = integ_(&d__1) + 6;
				d__1 = (prec2_1.gz[*i__ - 1] + prec4_1.vz[*i__ - 1] * ddt) / (ilist3_1.size3 + ilist3_1.v3 * dtt);
				i3 = integ_(&d__1) + 6;
			}
			if (k == 1)
			{
				i1 = 1;
			}
			if (k == 2)
			{
				i1 = 10;
			}
			if (k == 3)
			{
				i2 = 1;
			}
			if (k == 4)
			{
				i2 = 10;
			}
			if (k == 5)
			{
				i3 = 1;
			}
			if (k == 6)
			{
				i3 = 10;
			}
			j = ilist2_1.icel[i1 + (i2 + i3 * 10) * 10 - 111];
			newcre_(i__, &j);
			ilist2_1.icel[i1 + (i2 + i3 * 10) * 10 - 111] = j;
			ilist1_1.icels[*i__ - 1] = i1 * 10000 + i2 * 100 + i3;
		}
	}
	if (ilist1_1.next[*i__ - 1] != 0)
	{
		otmp = ilist5_1.ot[ilist1_1.next[*i__ - 1] - 1];
		ctmp = ilist5_1.ct[ilist1_1.next[*i__ - 1] - 1];
	}
	if (i1 == 11 && i2 == 11 && i3 == 11)
	{
		dchout_(i__, &k, t);
	}
	else
	{
		if (para5_1.iconfg == 1)
		{
			dchin1_(i__, &k, &i1, &i2, &i3, t);
		}
		else if (para5_1.iconfg == 2)
		{
			dchin2_(i__, &k, &i1, &i2, &i3, t);
		}
		else if (para5_1.iconfg == 4)
		{
			dchin3_(i__, &k, &i1, &i2, &i3, t);
		}
	}
	if (ilist1_1.icsta[*i__ - 1] / 10 == 11)
	{
		ilist5_1.ot[ilist1_1.next[*i__ - 1] - 1] = otmp;
		ilist5_1.ct[ilist1_1.next[*i__ - 1] - 1] = ctmp;
		ilist1_1.next[ilist1_1.next[*i__ - 1] - 1] = *i__;
		wallc_(i__, &i1, &i2, &i3, &t0, &tmin1);
		if (tmin1 < ilist5_1.ct[*i__ - 1])
		{
			ilist1_1.icsta[*i__ - 1] += 10;
			t0 = tmin1;
			goto L1000;
		}
	}
	return 0;
}

int oldcre_(int *i__)
{
	static int j;

	if (ilist1_1.nic[*i__ - 1] == 0)
	{
		return 0;
	}
	j = ilist1_1.nic[*i__ - 1];
	if (ilist1_1.nic[j - 1] == *i__)
	{
		ilist1_1.nic[j - 1] = 0;
		return 0;
	}
	while (ilist1_1.nic[j - 1] != *i__)
	{
		j = ilist1_1.nic[j - 1];
	}
	ilist1_1.nic[j - 1] = ilist1_1.nic[*i__ - 1];
	return 0;
}

int newcre_(int *i__, int *k)
{
	static int j;

	if (*k == 0)
	{
		*k = *i__;
		ilist1_1.nic[*i__ - 1] = 0;
	}
	else if (ilist1_1.nic[*k - 1] == 0)
	{
		ilist1_1.nic[*k - 1] = *i__;
		ilist1_1.nic[*i__ - 1] = *k;
	}
	else
	{
		j = *k;
		while (ilist1_1.nic[j - 1] != *k)
		{
			j = ilist1_1.nic[j - 1];
		}
		ilist1_1.nic[j - 1] = *i__;
		ilist1_1.nic[*i__ - 1] = *k;
	}
	return 0;
}

int scat_(double *t, int *iscat, int *jscat)
{
	extern int newmom_(double *), newpos_(double *,
											  int *);

	newpos_(t, iscat);
	newpos_(t, jscat);
	newmom_(t);
	return 0;
}

int newpos_(double *t, int *i__)
{

	double d__1;

	double log(double), cosh(double);

	static double dt1;

	dt1 = ilist5_1.ct[*i__ - 1] - prec2_1.ft[*i__ - 1];
	prec2_1.gx[*i__ - 1] += prec4_1.vx[*i__ - 1] * dt1;
	prec2_1.gy[*i__ - 1] += prec4_1.vy[*i__ - 1] * dt1;
	prec2_1.gz[*i__ - 1] += prec4_1.vz[*i__ - 1] * dt1;
	prec2_1.ft[*i__ - 1] = ilist5_1.ct[*i__ - 1];
	if (para5_1.iconfg <= 3)
	{
		if (prec2_1.ft[*i__ - 1] <= (d__1 = prec2_1.gz[*i__ - 1], abs(d__1)))
		{
			prec5_1.eta[*i__ - 1] = 1e6;
		}
		else
		{
			prec5_1.eta[*i__ - 1] = log((prec2_1.ft[*i__ - 1] + prec2_1.gz[*i__ - 1]) / (prec2_1.ft[*i__ - 1] - prec2_1.gz[*i__ - 1])) * .5;
		}
		prec5_1.tau[*i__ - 1] = prec2_1.ft[*i__ - 1] / cosh(prec5_1.eta[*i__ - 1]);
	}
	return 0;
}

int newmom_(double *t)
{

	double d__1, d__2, d__3, d__4;

	double acos(double), sqrt(double), log(double);

	static double e1, e2;
	static int i1, j1, i2, j2, k1, k2;
	static double t1, t2, x1, y1, z1, x2, y2, z2, pp2, px1, py1, pz1, px2,
		py2, pz2, bex, bey, bez, rap1, rap2, rts2, that, theta;
	extern int getht_(int *, int *, double *,
					  double *);
	static int icels1, icels2;
	extern int cropro_(double *, double *,
					   double *, double *, double *, double *),
		lorenz_(
			double *, double *, double *, double *,
			double *, double *, double *),
		zprota_(double *,
				double *, double *, double *, double *,
				double *, double *),
		xnormv_(double *, double *,
				double *);

	if (anim_1.isoft == 5)
	{
		if (frzprc_1.ifrz[ilist1_1.iscat - 1] == 1 || frzprc_1.ifrz[ilist1_1.jscat - 1] == 1)
		{
			ilist1_1.last[ilist1_1.iscat - 1] = ilist1_1.jscat;
			ilist1_1.last[ilist1_1.jscat - 1] = ilist1_1.iscat;
			return 0;
		}
	}

	rndm2_1.iff = -rndm2_1.iff;
	if (para5_1.iconfg == 2 || para5_1.iconfg == 4)
	{
		icels1 = ilist1_1.icels[ilist1_1.iscat - 1];
		i1 = icels1 / 10000;
		j1 = (icels1 - i1 * 10000) / 100;
		icels2 = ilist1_1.icels[ilist1_1.jscat - 1];
		i2 = icels2 / 10000;
		j2 = (icels2 - i2 * 10000) / 100;
		if (para5_1.iconfg == 4)
		{
			k1 = icels1 - i1 * 10000 - j1 * 100;
			k2 = icels2 - i2 * 10000 - j2 * 100;
		}
	}
	px1 = prec2_1.px[ilist1_1.iscat - 1];
	py1 = prec2_1.py[ilist1_1.iscat - 1];
	pz1 = prec2_1.pz[ilist1_1.iscat - 1];
	e1 = prec2_1.e[ilist1_1.iscat - 1];
	x1 = prec2_1.gx[ilist1_1.iscat - 1];
	y1 = prec2_1.gy[ilist1_1.iscat - 1];
	z1 = prec2_1.gz[ilist1_1.iscat - 1];
	t1 = prec2_1.ft[ilist1_1.iscat - 1];
	px2 = prec2_1.px[ilist1_1.jscat - 1];
	py2 = prec2_1.py[ilist1_1.jscat - 1];
	pz2 = prec2_1.pz[ilist1_1.jscat - 1];
	e2 = prec2_1.e[ilist1_1.jscat - 1];
	if (para5_1.iconfg == 1)
	{
		x2 = prec2_1.gx[ilist1_1.jscat - 1];
		y2 = prec2_1.gy[ilist1_1.jscat - 1];
		z2 = prec2_1.gz[ilist1_1.jscat - 1];
	}
	else if (para5_1.iconfg == 2 || para5_1.iconfg == 4)
	{
		if (i1 - i2 > 5)
		{
			x2 = prec2_1.gx[ilist1_1.jscat - 1] + ilist3_1.size1 * 10.;
		}
		else if (i1 - i2 < -5)
		{
			x2 = prec2_1.gx[ilist1_1.jscat - 1] - ilist3_1.size1 * 10.;
		}
		else
		{
			x2 = prec2_1.gx[ilist1_1.jscat - 1];
		}
		if (j1 - j2 > 5)
		{
			y2 = prec2_1.gy[ilist1_1.jscat - 1] + ilist3_1.size2 * 10.;
		}
		else if (j1 - j2 < -5)
		{
			y2 = prec2_1.gy[ilist1_1.jscat - 1] - ilist3_1.size2 * 10.;
		}
		else
		{
			y2 = prec2_1.gy[ilist1_1.jscat - 1];
		}
		if (para5_1.iconfg == 4)
		{
			if (k1 - k2 > 5)
			{
				z2 = prec2_1.gz[ilist1_1.jscat - 1] + ilist3_1.size3 * 10.;
			}
			else if (k1 - k2 < -5)
			{
				z2 = prec2_1.gz[ilist1_1.jscat - 1] - ilist3_1.size3 * 10.;
			}
			else
			{
				z2 = prec2_1.gz[ilist1_1.jscat - 1];
			}
		}
		else
		{
			z2 = prec2_1.gz[ilist1_1.jscat - 1];
		}
	}
	else if (para5_1.iconfg == 3 || para5_1.iconfg == 5)
	{
		x2 = prec2_1.gx[ilist1_1.jscat - 1] + aurec2_1.dgxa[ilist1_1.jscat -
															1];
		y2 = prec2_1.gy[ilist1_1.jscat - 1] + aurec2_1.dgya[ilist1_1.jscat -
															1];
		if (para5_1.iconfg == 5)
		{
			z2 = prec2_1.gz[ilist1_1.jscat - 1] + aurec2_1.dgza[ilist1_1.jscat - 1];
		}
		else
		{
			z2 = prec2_1.gz[ilist1_1.jscat - 1];
		}
	}
	t2 = prec2_1.ft[ilist1_1.jscat - 1];

	d__1 = e1 + e2;

	d__2 = px1 + px2;

	d__3 = py1 + py2;

	d__4 = pz1 + pz2;
	rts2 = d__1 * d__1 - d__2 * d__2 - d__3 * d__3 - d__4 * d__4;

	bex = (px1 + px2) / (e1 + e2);
	bey = (py1 + py2) / (e1 + e2);
	bez = (pz1 + pz2) / (e1 + e2);
	lorenz_(&e1, &px1, &py1, &pz1, &bex, &bey, &bez);

	px1 = lor_1.pxnew;
	py1 = lor_1.pynew;
	pz1 = lor_1.pznew;
	e1 = lor_1.enenew;

	d__1 = lor_1.pxnew;

	d__2 = lor_1.pynew;

	d__3 = lor_1.pznew;
	pp2 = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
	getht_(&ilist1_1.iscat, &ilist1_1.jscat, &pp2, &that);
	theta = acos(that / (pp2 * 2.) + 1.);
	theta = (double)rndm2_1.iff * theta;

	lorenz_(&t1, &x1, &y1, &z1, &bex, &bey, &bez);
	x1 = lor_1.pxnew;
	y1 = lor_1.pynew;
	z1 = lor_1.pznew;
	lorenz_(&t2, &x2, &y2, &z2, &bex, &bey, &bez);
	x2 = lor_1.pxnew;
	y2 = lor_1.pynew;
	z2 = lor_1.pznew;

	d__1 = x1 - x2;
	d__2 = y1 - y2;
	d__3 = z1 - z2;
	cropro_(&d__1, &d__2, &d__3, &px1, &py1, &pz1);
	xnormv_(&cprod_1.xn1, &cprod_1.xn2, &cprod_1.xn3);

	zprota_(&cprod_1.xn1, &cprod_1.xn2, &cprod_1.xn3, &theta, &px1, &py1, &pz1);

	px2 = -px1;
	py2 = -py1;
	pz2 = -pz1;

	d__1 = px2;

	d__2 = py2;

	d__3 = pz2;

	d__4 = prec2_1.xmass[ilist1_1.jscat - 1];
	e2 = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3 + d__4 * d__4);

	d__1 = -bex;
	d__2 = -bey;
	d__3 = -bez;
	lorenz_(&e1, &px1, &py1, &pz1, &d__1, &d__2, &d__3);
	prec2_1.px[ilist1_1.iscat - 1] = lor_1.pxnew;
	prec2_1.py[ilist1_1.iscat - 1] = lor_1.pynew;
	prec2_1.pz[ilist1_1.iscat - 1] = lor_1.pznew;
	prec2_1.e[ilist1_1.iscat - 1] = lor_1.enenew;
	d__1 = -bex;
	d__2 = -bey;
	d__3 = -bez;
	lorenz_(&e2, &px2, &py2, &pz2, &d__1, &d__2, &d__3);
	prec2_1.px[ilist1_1.jscat - 1] = lor_1.pxnew;
	prec2_1.py[ilist1_1.jscat - 1] = lor_1.pynew;
	prec2_1.pz[ilist1_1.jscat - 1] = lor_1.pznew;
	prec2_1.e[ilist1_1.jscat - 1] = lor_1.enenew;
	prec4_1.vx[ilist1_1.iscat - 1] = prec2_1.px[ilist1_1.iscat - 1] /
									 prec2_1.e[ilist1_1.iscat - 1];
	prec4_1.vy[ilist1_1.iscat - 1] = prec2_1.py[ilist1_1.iscat - 1] /
									 prec2_1.e[ilist1_1.iscat - 1];
	prec4_1.vz[ilist1_1.iscat - 1] = prec2_1.pz[ilist1_1.iscat - 1] /
									 prec2_1.e[ilist1_1.iscat - 1];
	prec4_1.vx[ilist1_1.jscat - 1] = prec2_1.px[ilist1_1.jscat - 1] /
									 prec2_1.e[ilist1_1.jscat - 1];
	prec4_1.vy[ilist1_1.jscat - 1] = prec2_1.py[ilist1_1.jscat - 1] /
									 prec2_1.e[ilist1_1.jscat - 1];
	prec4_1.vz[ilist1_1.jscat - 1] = prec2_1.pz[ilist1_1.jscat - 1] /
									 prec2_1.e[ilist1_1.jscat - 1];
	ilist1_1.last[ilist1_1.iscat - 1] = ilist1_1.jscat;
	ilist1_1.last[ilist1_1.jscat - 1] = ilist1_1.iscat;
	if (para5_1.iconfg <= 3)
	{
		if (prec2_1.e[ilist1_1.iscat - 1] <= (d__1 = prec2_1.pz[ilist1_1.iscat - 1], abs(d__1)))
		{
			prec5_1.rap[ilist1_1.iscat - 1] = 1e6;
		}
		else
		{
			prec5_1.rap[ilist1_1.iscat - 1] = log((prec2_1.e[ilist1_1.iscat -
															 1] +
												   prec2_1.pz[ilist1_1.iscat - 1]) /
												  (prec2_1.e[ilist1_1.iscat - 1] - prec2_1.pz[ilist1_1.iscat - 1])) *
											  .5;
		}
		if (prec2_1.e[ilist1_1.jscat - 1] <= (d__1 = prec2_1.pz[ilist1_1.jscat - 1], abs(d__1)))
		{
			prec5_1.rap[ilist1_1.jscat - 1] = 1e6;
		}
		else
		{
			prec5_1.rap[ilist1_1.jscat - 1] = log((prec2_1.e[ilist1_1.jscat -
															 1] +
												   prec2_1.pz[ilist1_1.jscat - 1]) /
												  (prec2_1.e[ilist1_1.jscat - 1] - prec2_1.pz[ilist1_1.jscat - 1])) *
											  .5;
		}

		rap1 = prec5_1.rap[ilist1_1.iscat - 1];
		rap2 = prec5_1.rap[ilist1_1.jscat - 1];
		if (rap1 < para6_1.centy + .5 && rap1 > para6_1.centy - .5)
		{
		}
		if (rap2 < para6_1.centy + .5 && rap2 > para6_1.centy - .5)
		{
		}
	}
	return 0;
}

int getht_(int *iscat, int *jscat, double *pp2,
		   double *that)
{

	double d__1;

	static double rx, xm2;
	extern double ran1_(int *);
	static double xmp2, xmu2;
	static int iseed;

	iseed = rndm3_1.iseedp;

	d__1 = para2_1.xmu * .197327054;
	xmu2 = d__1 * d__1;

	d__1 = para2_1.xmp;
	xmp2 = d__1 * d__1;
	xm2 = xmu2 + xmp2;
	rx = ran1_(&iseed);
	*that = xm2 * (1. / ((1. - xm2 / (*pp2 * 4. + xm2)) * rx - 1.) + 1.);

	if (anim_1.izpc == 100)
	{
		*that = *pp2 * -4. * rx;
	}
	return 0;
}

int ulist_(double *t)
{
	static int l;
	extern int ulist1_(int *, double *);

	if (ilist1_1.ictype == 1 || ilist1_1.ictype == 2 || ilist1_1.ictype == 5 || ilist1_1.ictype == 6)
	{
		l = ilist4_1.ifmpt;
		ulist1_(&l, t);
	}
	if (ilist1_1.ictype != 1)
	{
		l = ilist1_1.iscat;
		ulist1_(&l, t);
		if (ilist1_1.jscat != 0)
		{
			l = ilist1_1.jscat;
			ulist1_(&l, t);
		}
	}
	return 0;
}

int ulist1_(int *l, double *t)
{
	static int k, i1, i2, i3, nc;
	static double tmin, tmin1;
	extern int wallc_(int *, int *, int *,
					  int *, double *, double *),
		chkin1_(int *,
				int *, int *, int *, double *, double *,
				int *),
		chkin2_(int *, int *, int *, int *,
				double *, double *, int *);
	static int icels0;
	extern int chkin3_(int *, int *, int *,
					   int *, double *, double *, int *),
		chkcel_(
			int *, int *, int *, int *, double *,
			double *, int *),
		chkout_(int *, double *,
				double *, int *),
		fixtim_(int *, double *,
				double *, double *, int *);

	icels0 = ilist1_1.icels[*l - 1];
	i1 = icels0 / 10000;
	i2 = (icels0 - i1 * 10000) / 100;
	i3 = icels0 - i1 * 10000 - i2 * 100;

	k = ilist1_1.icsta[*l - 1] % 10;
	wallc_(l, &i1, &i2, &i3, t, &tmin1);
	tmin = tmin1;
	nc = 0;
	if (i1 == 11 && i2 == 11 && i3 == 11)
	{
		chkout_(l, t, &tmin, &nc);
	}
	else
	{
		if (para5_1.iconfg == 1)
		{
			chkin1_(l, &i1, &i2, &i3, t, &tmin, &nc);
		}
		else if (para5_1.iconfg == 2)
		{
			chkin2_(l, &i1, &i2, &i3, t, &tmin, &nc);
		}
		else if (para5_1.iconfg == 4)
		{
			chkin3_(l, &i1, &i2, &i3, t, &tmin, &nc);
		}
		else if (para5_1.iconfg == 3 || para5_1.iconfg == 5)
		{
			chkcel_(l, &i1, &i2, &i3, t, &tmin, &nc);
		}
	}
	fixtim_(l, t, &tmin1, &tmin, &nc);
	return 0;
}

int wallc_(int *i__, int *i1, int *i2, int *i3, double *t, double *tmin)
{
	extern int wallc1_(int *, int *, int *,
					   int *, double *, double *),
		wallc2_(int *,
				int *, int *, int *, double *, double *),
		wallcb_(int *, double *, double *);

	*tmin = ilist5_1.tlarge;
	if (para5_1.iconfg <= 2 || para5_1.iconfg == 4)
	{

		if (*i1 >= 1 && *i1 <= 10 || *i2 >= 1 && *i2 <= 10 || *i3 >= 1 && *i3 <= 10)
		{
			wallc1_(i__, i1, i2, i3, t, tmin);
		}
		else
		{
			wallcb_(i__, t, tmin);
		}
	}
	else if (para5_1.iconfg == 3 || para5_1.iconfg == 5)
	{
		wallc2_(i__, i1, i2, i3, t, tmin);
	}
	return 0;
}

int wallc1_(int *i__, int *i1, int *i2, int *i3, double *t, double *tmin)
{

	double d__1;

	static double t1, t2, t3, tf, v1p, v2p, x1p, x2p, x3p, v3p;

	x1p = prec2_1.gx[*i__ - 1];
	x2p = prec2_1.gy[*i__ - 1];
	x3p = prec2_1.gz[*i__ - 1];
	tf = prec2_1.ft[*i__ - 1];
	v1p = prec4_1.vx[*i__ - 1];
	v2p = prec4_1.vy[*i__ - 1];
	v3p = prec4_1.vz[*i__ - 1];
	if (*t < ilist3_1.size && tf < ilist3_1.size)
	{
		if (v1p > 0.)
		{
			t1 = (((double)(*i1) - 5.) * ilist3_1.size1 - x1p) / v1p +
				 tf;
		}
		else if (v1p < 0.)
		{
			t1 = (((double)(*i1) - 6.) * ilist3_1.size1 - x1p) / v1p +
				 tf;
		}
		else
		{
			t1 = ilist5_1.tlarge;
		}
		if (v2p > 0.)
		{
			t2 = (((double)(*i2) - 5.) * ilist3_1.size2 - x2p) / v2p +
				 tf;
		}
		else if (v2p < 0.)
		{
			t2 = (((double)(*i2) - 6.) * ilist3_1.size2 - x2p) / v2p +
				 tf;
		}
		else
		{
			t2 = ilist5_1.tlarge;
		}
		if (v3p > 0.)
		{
			t3 = (((double)(*i3) - 5.) * ilist3_1.size3 - x3p) / v3p +
				 tf;
		}
		else if (v3p < 0.)
		{
			t3 = (((double)(*i3) - 6.) * ilist3_1.size3 - x3p) / v3p +
				 tf;
		}
		else
		{
			t3 = ilist5_1.tlarge;
		}

		d__1 = min(t1, t2);
		*tmin = min(d__1, t3);

		if (*tmin == t1)
		{
			if (v1p > 0.)
			{
				ilist1_1.icsta[*i__ - 1] = 101;
			}
			else
			{
				ilist1_1.icsta[*i__ - 1] = 102;
			}
		}
		if (*tmin == t2)
		{
			if (v2p > 0.)
			{
				ilist1_1.icsta[*i__ - 1] = 103;
			}
			else
			{
				ilist1_1.icsta[*i__ - 1] = 104;
			}
		}
		if (*tmin == t3)
		{
			if (v3p > 0.)
			{
				ilist1_1.icsta[*i__ - 1] = 105;
			}
			else
			{
				ilist1_1.icsta[*i__ - 1] = 106;
			}
		}
		if (*tmin <= ilist3_1.size)
		{
			return 0;
		}
	}
	if (v1p > (*i1 - 5) * ilist3_1.v1)
	{
		t1 = ((*i1 - 5) * (ilist3_1.size1 - ilist3_1.v1 * ilist3_1.size) +
			  v1p * tf - x1p) /
			 (v1p - (*i1 - 5) * ilist3_1.v1);
	}
	else if (v1p < (*i1 - 6) * ilist3_1.v1)
	{
		t1 = ((*i1 - 6) * (ilist3_1.size1 - ilist3_1.v1 * ilist3_1.size) +
			  v1p * tf - x1p) /
			 (v1p - (*i1 - 6) * ilist3_1.v1);
	}
	else
	{
		t1 = ilist5_1.tlarge;
	}
	if (v2p > (*i2 - 5) * ilist3_1.v2)
	{
		t2 = ((*i2 - 5) * (ilist3_1.size2 - ilist3_1.v2 * ilist3_1.size) +
			  v2p * tf - x2p) /
			 (v2p - (*i2 - 5) * ilist3_1.v2);
	}
	else if (v2p < (*i2 - 6) * ilist3_1.v2)
	{
		t2 = ((*i2 - 6) * (ilist3_1.size2 - ilist3_1.v2 * ilist3_1.size) +
			  v2p * tf - x2p) /
			 (v2p - (*i2 - 6) * ilist3_1.v2);
	}
	else
	{
		t2 = ilist5_1.tlarge;
	}
	if (v3p > (*i3 - 5) * ilist3_1.v3)
	{
		t3 = ((*i3 - 5) * (ilist3_1.size3 - ilist3_1.v3 * ilist3_1.size) +
			  v3p * tf - x3p) /
			 (v3p - (*i3 - 5) * ilist3_1.v3);
	}
	else if (v3p < (*i3 - 6) * ilist3_1.v3)
	{
		t3 = ((*i3 - 6) * (ilist3_1.size3 - ilist3_1.v3 * ilist3_1.size) +
			  v3p * tf - x3p) /
			 (v3p - (*i3 - 6) * ilist3_1.v3);
	}
	else
	{
		t3 = ilist5_1.tlarge;
	}

	d__1 = min(t1, t2);
	*tmin = min(d__1, t3);

	if (*tmin == t1)
	{
		if (v1p > (*i1 - 5) * ilist3_1.v1)
		{
			ilist1_1.icsta[*i__ - 1] = 101;
		}
		else
		{
			ilist1_1.icsta[*i__ - 1] = 102;
		}
	}
	if (*tmin == t2)
	{
		if (v2p > (*i2 - 5) * ilist3_1.v2)
		{
			ilist1_1.icsta[*i__ - 1] = 103;
		}
		else
		{
			ilist1_1.icsta[*i__ - 1] = 104;
		}
	}
	if (*tmin == t3)
	{
		if (v3p > (*i3 - 5) * ilist3_1.v3)
		{
			ilist1_1.icsta[*i__ - 1] = 105;
		}
		else
		{
			ilist1_1.icsta[*i__ - 1] = 106;
		}
	}
	return 0;
}

int wallc2_(int *i__, int *i1, int *i2, int *i3, double *t, double *tmin)
{

	double d__1;

	static double t1, t2, t3, tf, v1p, v2p, x1p, x2p, x3p, v3p;

	x1p = prec2_1.gx[*i__ - 1];
	x2p = prec2_1.gy[*i__ - 1];
	x3p = prec2_1.gz[*i__ - 1];
	tf = prec2_1.ft[*i__ - 1];
	v1p = prec4_1.vx[*i__ - 1];
	v2p = prec4_1.vy[*i__ - 1];
	v3p = prec4_1.vz[*i__ - 1];
	if (v1p > 0.)
	{
		t1 = (ilist3_1.size1 * 5. - x1p) / v1p + tf;
	}
	else if (v1p < 0.)
	{
		t1 = (ilist3_1.size1 * -5. - x1p) / v1p + tf;
	}
	else
	{
		t1 = ilist5_1.tlarge;
	}
	if (v2p > 0.)
	{
		t2 = (ilist3_1.size2 * 5. - x2p) / v2p + tf;
	}
	else if (v2p < 0.)
	{
		t2 = (ilist3_1.size2 * -5. - x2p) / v2p + tf;
	}
	else
	{
		t2 = ilist5_1.tlarge;
	}
	if (para5_1.iconfg == 5)
	{
		if (v3p > 0.)
		{
			t3 = (ilist3_1.size3 * 5. - x3p) / v3p + tf;
		}
		else if (v3p < 0.)
		{
			t3 = (ilist3_1.size3 * -5. - x3p) / v3p + tf;
		}
		else
		{
			t3 = ilist5_1.tlarge;
		}
	}
	else
	{
		t3 = ilist5_1.tlarge;
	}

	d__1 = min(t1, t2);
	*tmin = min(d__1, t3);

	if (*tmin == t1)
	{
		if (v1p > 0.)
		{
			ilist1_1.icsta[*i__ - 1] = 101;
		}
		else
		{
			ilist1_1.icsta[*i__ - 1] = 102;
		}
	}
	if (*tmin == t2)
	{
		if (v2p > 0.)
		{
			ilist1_1.icsta[*i__ - 1] = 103;
		}
		else
		{
			ilist1_1.icsta[*i__ - 1] = 104;
		}
	}
	if (*tmin == t3)
	{
		if (v3p > 0.)
		{
			ilist1_1.icsta[*i__ - 1] = 105;
		}
		else
		{
			ilist1_1.icsta[*i__ - 1] = 106;
		}
	}
	return 0;
}

int wallcb_(int *i__, double *t, double *tmin)
{

	double d__1;

	static double t1, t2, t3, tf, v1p, v2p, x1p, x2p, x3p, v3p, x1q, x2q,
		x3q, x1pp, x2pp, x3pp;
	static int icsta1, icsta2, icsta3;

	if (ilist3_1.size1 == 0. || ilist3_1.size2 == 0. || ilist3_1.size3 == 0.)
	{
		return 0;
	}
	x1p = prec2_1.gx[*i__ - 1];
	x2p = prec2_1.gy[*i__ - 1];
	x3p = prec2_1.gz[*i__ - 1];
	v1p = prec4_1.vx[*i__ - 1];
	v2p = prec4_1.vy[*i__ - 1];
	v3p = prec4_1.vz[*i__ - 1];
	tf = prec2_1.ft[*i__ - 1];
	if (*t < ilist3_1.size && tf < ilist3_1.size)
	{
		if (x1p < ilist3_1.size1 * -5. && v1p > 0.)
		{
			t1 = (ilist3_1.size1 * -5. - x1p) / v1p + tf;
		}
		else if (x1p > ilist3_1.size1 * 5. && v1p < 0.)
		{
			t1 = -(x1p - ilist3_1.size1 * 5.) / v1p + tf;
		}
		else
		{
			t1 = ilist5_1.tlarge;
		}
		if (t1 != ilist5_1.tlarge)
		{
			x2pp = x2p + v2p * (t1 - tf);
			x3pp = x3p + v3p * (t1 - tf);
			if (x2pp <= ilist3_1.size2 * -5. || x2pp >= ilist3_1.size2 * 5. ||
				x3pp <= ilist3_1.size3 * -5. || x3pp >= ilist3_1.size3 * 5.)
			{
				t1 = ilist5_1.tlarge;
			}
		}
		if (x2p < ilist3_1.size2 * -5. && v2p > 0.)
		{
			t2 = (ilist3_1.size2 * -5. - x2p) / v2p + tf;
		}
		else if (x2p > ilist3_1.size2 * 5. && v2p < 0.)
		{
			t2 = -(x2p - ilist3_1.size2 * 5.) / v2p + tf;
		}
		else
		{
			t2 = ilist5_1.tlarge;
		}
		if (t2 != ilist5_1.tlarge)
		{
			x1pp = x1p + v1p * (t2 - tf);
			x3pp = x3p + v3p * (t2 - tf);
			if (x1pp <= ilist3_1.size1 * -5. || x1pp >= ilist3_1.size1 * 5. ||
				x3pp <= ilist3_1.size3 * -5. || x3pp >= ilist3_1.size3 * 5.)
			{
				t2 = ilist5_1.tlarge;
			}
		}
		if (x3p < ilist3_1.size3 * -5. && v3p > 0.)
		{
			t3 = (ilist3_1.size3 * -5. - x3p) / v3p + tf;
		}
		else if (x3p > ilist3_1.size3 * 5. && v3p < 0.)
		{
			t3 = -(x3p - ilist3_1.size3 * 5.) / v3p + tf;
		}
		else
		{
			t3 = ilist5_1.tlarge;
		}
		if (t3 != ilist5_1.tlarge)
		{
			x1pp = x1p + v1p * (t3 - tf);
			x2pp = x2p + v2p * (t3 - tf);
			if (x1pp <= ilist3_1.size1 * -5. || x1pp >= ilist3_1.size1 * 5. ||
				x2pp <= ilist3_1.size2 * -5. || x2pp >= ilist3_1.size2 * 5.)
			{
				t3 = ilist5_1.tlarge;
			}
		}

		d__1 = min(t1, t2);
		*tmin = min(d__1, t3);

		if (*tmin == t1)
		{
			if (v1p > 0.)
			{
				ilist1_1.icsta[*i__ - 1] = 101;
			}
			else
			{
				ilist1_1.icsta[*i__ - 1] = 102;
			}
		}
		if (*tmin == t2)
		{
			if (v2p > 0.)
			{
				ilist1_1.icsta[*i__ - 1] = 103;
			}
			else
			{
				ilist1_1.icsta[*i__ - 1] = 104;
			}
		}
		if (*tmin == t3)
		{
			if (v3p > 0.)
			{
				ilist1_1.icsta[*i__ - 1] = 105;
			}
			else
			{
				ilist1_1.icsta[*i__ - 1] = 106;
			}
		}
		if (*tmin <= ilist3_1.size)
		{
			return 0;
		}
	}

	x1q = x1p + v1p * (*t - tf);
	x2q = x2p + v2p * (*t - tf);
	x3q = x3p + v3p * (*t - tf);
	if (x1q < (ilist3_1.size1 + ilist3_1.v1 * (*t - ilist3_1.size)) * -5. &&
		v1p > ilist3_1.v1 * -5.)
	{
		t1 = ((ilist3_1.size1 - ilist3_1.v1 * ilist3_1.size) * -5. + v1p * tf - x1p) / (v1p - ilist3_1.v1 * -5.);
		icsta1 = 101;
	}
	else if (x1q > (ilist3_1.size1 + ilist3_1.v1 * (*t - ilist3_1.size)) *
					   5. &&
			 v1p < ilist3_1.v1 * 5.)
	{
		t1 = ((ilist3_1.size1 - ilist3_1.v1 * ilist3_1.size) * 5. + v1p * tf - x1p) / (v1p - ilist3_1.v1 * 5.);
		icsta1 = 102;
	}
	else
	{
		t1 = ilist5_1.tlarge;
	}
	if (t1 != ilist5_1.tlarge)
	{
		x2pp = x2p + v2p * (t1 - tf);
		x3pp = x3p + v3p * (t1 - tf);
		if (x2pp <= (ilist3_1.size2 + ilist3_1.v2 * (t1 - ilist3_1.size)) *
						-5. ||
			x2pp >= (ilist3_1.size2 + ilist3_1.v2 * (t1 -
													 ilist3_1.size)) *
						5. ||
			x3pp <= (ilist3_1.size3 + ilist3_1.v3 * (t1 - ilist3_1.size)) * -5. || x3pp >= (ilist3_1.size3 + ilist3_1.v3 * (t1 - ilist3_1.size)) * 5.)
		{
			t1 = ilist5_1.tlarge;
		}
	}
	if (x2q < (ilist3_1.size2 + ilist3_1.v2 * (*t - ilist3_1.size)) * -5. &&
		v2p > ilist3_1.v2 * -5.)
	{
		t2 = ((ilist3_1.size2 - ilist3_1.v2 * ilist3_1.size) * -5. + v2p * tf - x2p) / (v2p - ilist3_1.v2 * -5.);
		icsta2 = 103;
	}
	else if (x2q > (ilist3_1.size2 + ilist3_1.v2 * (*t - ilist3_1.size)) *
					   5. &&
			 v2p < ilist3_1.v2 * 5.)
	{
		t2 = ((ilist3_1.size2 - ilist3_1.v2 * ilist3_1.size) * 5. + v2p * tf - x2p) / (v2p - ilist3_1.v2 * 5.);
		icsta2 = 104;
	}
	else
	{
		t2 = ilist5_1.tlarge;
	}
	if (t2 != ilist5_1.tlarge)
	{
		x1pp = x1p + v1p * (t2 - tf);
		x3pp = x3p + v3p * (t2 - tf);
		if (x1pp <= (ilist3_1.size1 + ilist3_1.v1 * (t2 - ilist3_1.size)) *
						-5. ||
			x1pp >= (ilist3_1.size1 + ilist3_1.v1 * (t2 -
													 ilist3_1.size)) *
						5. ||
			x3pp <= (ilist3_1.size3 + ilist3_1.v3 * (t2 - ilist3_1.size)) * -5. || x3pp >= (ilist3_1.size3 + ilist3_1.v3 * (t2 - ilist3_1.size)) * 5.)
		{
			t2 = ilist5_1.tlarge;
		}
	}
	if (x3q < (ilist3_1.size3 + ilist3_1.v3 * (*t - ilist3_1.size)) * -5. &&
		v3p > ilist3_1.v3 * -5.)
	{
		t3 = ((ilist3_1.size3 - ilist3_1.v3 * ilist3_1.size) * -5. + v3p * tf - x3p) / (v3p - ilist3_1.v3 * -5.);
		icsta3 = 105;
	}
	else if (x3q > (ilist3_1.size3 + ilist3_1.v3 * (*t - ilist3_1.size)) *
					   5. &&
			 v3p < ilist3_1.v3 * 5.)
	{
		t3 = ((ilist3_1.size3 - ilist3_1.v3 * ilist3_1.size) * 5. + v3p * tf - x3p) / (v3p - ilist3_1.v3 * 5.);
		icsta3 = 106;
	}
	else
	{
		t3 = ilist5_1.tlarge;
	}
	if (t3 != ilist5_1.tlarge)
	{
		x2pp = x2p + v2p * (t3 - tf);
		x1pp = x1p + v1p * (t3 - tf);
		if (x2pp <= (ilist3_1.size2 + ilist3_1.v2 * (t3 - ilist3_1.size)) *
						-5. ||
			x2pp >= (ilist3_1.size2 + ilist3_1.v2 * (t3 -
													 ilist3_1.size)) *
						5. ||
			x1pp <= (ilist3_1.size1 + ilist3_1.v1 * (t3 - ilist3_1.size)) * -5. || x1pp >= (ilist3_1.size1 + ilist3_1.v1 * (t3 - ilist3_1.size)) * 5.)
		{
			t3 = ilist5_1.tlarge;
		}
	}

	d__1 = min(t1, t2);
	*tmin = min(d__1, t3);

	if (*tmin == t1)
	{
		ilist1_1.icsta[*i__ - 1] = icsta1;
	}
	else if (*tmin == t2)
	{
		ilist1_1.icsta[*i__ - 1] = icsta2;
	}
	else if (*tmin == t3)
	{
		ilist1_1.icsta[*i__ - 1] = icsta3;
	}
	return 0;
}

int chkout_(int *l, double *t, double *tmin,
			int *nc)
{
	static int i__, j, k, m1, m2, m3;
	extern int chkcel_(int *, int *, int *,
					   int *, double *, double *, int *);

	m1 = 11;
	m2 = 11;
	m3 = 11;
	chkcel_(l, &m1, &m2, &m3, t, tmin, nc);
	for (i__ = 1; i__ <= 10; ++i__)
	{
		for (j = 1; j <= 10; ++j)
		{
			for (k = 1; k <= 10; ++k)
			{
				if (i__ == 1 || i__ == 10 || j == 1 || j == 10 || k == 1 || k == 10)
				{
					chkcel_(l, &i__, &j, &k, t, tmin, nc);
				}
			}
		}
	}
	return 0;
}

int chkin1_(int *l, int *i1, int *i2, int *i3, double *t, double *tmin, int *nc)
{

	int i__1, i__2, i__3;

	static int i__, j, k, m1, m2, m3, itest;
	extern int chkcel_(int *, int *, int *,
					   int *, double *, double *, int *);

	itest = 0;
	i__1 = *i1 + 1;
	for (i__ = *i1 - 1; i__ <= i__1; ++i__)
	{
		i__2 = *i2 + 1;
		for (j = *i2 - 1; j <= i__2; ++j)
		{
			i__3 = *i3 + 1;
			for (k = *i3 - 1; k <= i__3; ++k)
			{
				if (i__ >= 1 && i__ <= 10 && j >= 1 && j <= 10 && k >= 1 && k <= 10)
				{
					chkcel_(l, &i__, &j, &k, t, tmin, nc);
				}
				else if (itest == 0)
				{
					m1 = 11;
					m2 = 11;
					m3 = 11;
					chkcel_(l, &m1, &m2, &m3, t, tmin, nc);
					itest = 1;
				}
			}
		}
	}
	return 0;
}

int chkin2_(int *l, int *i1, int *i2, int *i3, double *t, double *tmin, int *nc)
{

	int i__1, i__2, i__3;

	static int i__, j, k, ia, ib, ic, itest;
	extern int chkcel_(int *, int *, int *,
					   int *, double *, double *, int *);

	itest = 0;
	i__1 = *i1 + 1;
	for (i__ = *i1 - 1; i__ <= i__1; ++i__)
	{
		i__2 = *i2 + 1;
		for (j = *i2 - 1; j <= i__2; ++j)
		{
			i__3 = *i3 + 1;
			for (k = *i3 - 1; k <= i__3; ++k)
			{
				ia = i__;
				ib = j;
				ic = k;
				if (k >= 1 && k <= 10)
				{
					if (i__ == 0)
					{
						ia = 10;
					}
					if (i__ == 11)
					{
						ia = 1;
					}
					if (j == 0)
					{
						ib = 10;
					}
					if (j == 11)
					{
						ib = 1;
					}
					chkcel_(l, &ia, &ib, &ic, t, tmin, nc);
				}
			}
		}
	}
	return 0;
}

int chkin3_(int *l, int *i1, int *i2, int *i3, double *t, double *tmin, int *nc)
{

	int i__1, i__2, i__3;

	static int i__, j, k, ia, ib, ic, itest;
	extern int chkcel_(int *, int *, int *,
					   int *, double *, double *, int *);

	itest = 0;
	i__1 = *i1 + 1;
	for (i__ = *i1 - 1; i__ <= i__1; ++i__)
	{
		i__2 = *i2 + 1;
		for (j = *i2 - 1; j <= i__2; ++j)
		{
			i__3 = *i3 + 1;
			for (k = *i3 - 1; k <= i__3; ++k)
			{
				if (i__ == 0)
				{
					ia = 10;
				}
				else if (i__ == 11)
				{
					ia = 1;
				}
				else
				{
					ia = i__;
				}
				if (j == 0)
				{
					ib = 10;
				}
				else if (j == 11)
				{
					ib = 1;
				}
				else
				{
					ib = j;
				}
				if (k == 0)
				{
					ic = 10;
				}
				else if (k == 11)
				{
					ic = 1;
				}
				else
				{
					ic = k;
				}
				chkcel_(l, &ia, &ib, &ic, t, tmin, nc);
			}
		}
	}
	return 0;
}

int chkcel_(int *il, int *i1, int *i2, int *i3, double *t, double *tmin, int *nc)
{

	int i__1;

	static int j, l;
	extern int ck_(int *, int *);
	static int jj;
	extern int ud2_(int *, int *, double *,
					double *, int *);
	static int ick, jud2;

	if (para5_1.iconfg == 3 || para5_1.iconfg == 5)
	{
		jj = ilist4_1.ichkpt;
		i__1 = jj;
		for (j = 1; j <= i__1; ++j)
		{
			ck_(&j, &ick);

			jud2 = j;

			if (ick == 1)
			{
				ud2_(&jud2, il, t, tmin, nc);
			}
		}
		return 0;
	}
	if (*i1 == 11 && *i2 == 11 && *i3 == 11)
	{
		l = ilist2_1.icell;
	}
	else
	{
		l = ilist2_1.icel[*i1 + (*i2 + *i3 * 10) * 10 - 111];
	}

	if (l == 0)
	{
		return 0;
	}
	j = ilist1_1.nic[l - 1];

	if (j == 0)
	{
		ck_(&l, &ick);
		if (ick == 1)
		{
			ud2_(&l, il, t, tmin, nc);
		}
	}
	else
	{

		ck_(&l, &ick);
		if (ick == 1)
		{
			ud2_(&l, il, t, tmin, nc);
		}
		while (j != l)
		{
			ck_(&j, &ick);
			if (ick == 1)
			{
				ud2_(&j, il, t, tmin, nc);
			}
			j = ilist1_1.nic[j - 1];
		}
	}
	return 0;
}

int ck_(int *l, int *ick)
{

	*ick = 1;
	if (ilist1_1.ictype == 1)
	{
		if (*l == ilist4_1.ifmpt)
		{
			*ick = 0;
		}
	}
	else if (ilist1_1.ictype == 0 || ilist1_1.ictype == 3 ||
			 ilist1_1.ictype == 4)
	{
		if (*l == ilist1_1.iscat || *l == ilist1_1.jscat)
		{
			*ick = 0;
		}
	}
	else
	{
		if (*l == ilist1_1.iscat || *l == ilist1_1.jscat || *l == ilist4_1.ifmpt)
		{
			*ick = 0;
		}
	}

	return 0;
}

int dchout_(int *l, int *ii, double *t)
{

	int i__1, i__2;
	double d__1;

	static int i__, j, k, i1, i2, i3;
	static double x1, x2, x3, td, tt;
	extern int integ_(double *);
	extern int dchcel_(int *, int *, int *,
					   int *, double *);

	tt = prec2_1.ft[*l - 1];
	td = *t - ilist3_1.size;
	x1 = prec2_1.gx[*l - 1] + prec4_1.vx[*l - 1] * (*t - tt);
	x2 = prec2_1.gy[*l - 1] + prec4_1.vy[*l - 1] * (*t - tt);
	x3 = prec2_1.gz[*l - 1] + prec4_1.vz[*l - 1] * (*t - tt);
	if (td <= 0.)
	{
		d__1 = x1 / ilist3_1.size1;
		i1 = integ_(&d__1) + 6;
		d__1 = x2 / ilist3_1.size2;
		i2 = integ_(&d__1) + 6;
		d__1 = x3 / ilist3_1.size3;
		i3 = integ_(&d__1) + 6;
		d__1 = x1 / ilist3_1.size1;
		if ((double)integ_(&d__1) == x1 / ilist3_1.size1 && prec4_1.vx[*l - 1] < 0.)
		{
			--i1;
		}
		d__1 = x2 / ilist3_1.size2;
		if ((double)integ_(&d__1) == x2 / ilist3_1.size2 && prec4_1.vy[*l - 1] < 0.)
		{
			--i2;
		}
		d__1 = x3 / ilist3_1.size3;
		if ((double)integ_(&d__1) == x3 / ilist3_1.size3 && prec4_1.vz[*l - 1] < 0.)
		{
			--i3;
		}
	}
	else
	{
		d__1 = x1 / (ilist3_1.size1 + ilist3_1.v1 * td);
		i1 = integ_(&d__1) + 6;
		d__1 = x2 / (ilist3_1.size2 + ilist3_1.v2 * td);
		i2 = integ_(&d__1) + 6;
		d__1 = x3 / (ilist3_1.size3 + ilist3_1.v3 * td);
		i3 = integ_(&d__1) + 6;

		d__1 = x1 / (ilist3_1.size1 + ilist3_1.v1 * td);
		if ((double)integ_(&d__1) == x1 / (ilist3_1.size1 + ilist3_1.v1 *
																	td) &&
			prec4_1.vx[*l - 1] < (i1 - 6) * ilist3_1.v1)
		{
			--i1;
		}

		d__1 = x2 / (ilist3_1.size2 + ilist3_1.v2 * td);
		if ((double)integ_(&d__1) == x2 / (ilist3_1.size2 + ilist3_1.v2 *
																	td) &&
			prec4_1.vy[*l - 1] < (i2 - 6) * ilist3_1.v2)
		{
			--i2;
		}

		d__1 = x3 / (ilist3_1.size3 + ilist3_1.v3 * td);
		if ((double)integ_(&d__1) == x3 / (ilist3_1.size3 + ilist3_1.v3 *
																	td) &&
			prec4_1.vz[*l - 1] < (i3 - 6) * ilist3_1.v3)
		{
			--i3;
		}
	}
	if (*ii == 1)
	{
		i__ = 9;
		i__1 = i2 + 1;
		for (j = i2 - 1; j <= i__1; ++j)
		{
			i__2 = i3 + 1;
			for (k = i3 - 1; k <= i__2; ++k)
			{
				if (j >= 1 && j <= 10 && k >= 1 && k <= 10)
				{
					dchcel_(l, &i__, &j, &k, t);
				}
			}
		}
	}
	if (*ii == 2)
	{
		i__ = 2;
		i__1 = i2 + 1;
		for (j = i2 - 1; j <= i__1; ++j)
		{
			i__2 = i3 + 1;
			for (k = i3 - 1; k <= i__2; ++k)
			{
				if (j >= 1 && j <= 10 && k >= 1 && k <= 10)
				{
					dchcel_(l, &i__, &j, &k, t);
				}
			}
		}
	}
	if (*ii == 3)
	{
		j = 9;
		i__1 = i1 + 1;
		for (i__ = i1 - 1; i__ <= i__1; ++i__)
		{
			i__2 = i3 + 1;
			for (k = i3 - 1; k <= i__2; ++k)
			{
				if (i__ >= 1 && i__ <= 10 && k >= 1 && k <= 10)
				{
					dchcel_(l, &i__, &j, &k, t);
				}
			}
		}
	}
	if (*ii == 4)
	{
		j = 2;
		i__1 = i1 + 1;
		for (i__ = i1 - 1; i__ <= i__1; ++i__)
		{
			i__2 = i3 + 1;
			for (k = i3 - 1; k <= i__2; ++k)
			{
				if (i__ >= 1 && i__ <= 10 && k >= 1 && k <= 10)
				{
					dchcel_(l, &i__, &j, &k, t);
				}
			}
		}
	}
	if (*ii == 5)
	{
		k = 9;
		i__1 = i1 + 1;
		for (i__ = i1 - 1; i__ <= i__1; ++i__)
		{
			i__2 = i2 + 1;
			for (j = i2 - 1; j <= i__2; ++j)
			{
				if (i__ >= 1 && i__ <= 10 && j >= 1 && j <= 10)
				{
					dchcel_(l, &i__, &j, &k, t);
				}
			}
		}
	}
	if (*ii == 6)
	{
		k = 2;
		i__1 = i1 + 1;
		for (i__ = i1 - 1; i__ <= i__1; ++i__)
		{
			i__2 = i2 + 1;
			for (j = i2 - 1; j <= i__2; ++j)
			{
				if (i__ >= 1 && i__ <= 10 && j >= 1 && j <= 10)
				{
					dchcel_(l, &i__, &j, &k, t);
				}
			}
		}
	}
	return 0;
}

int dchin1_(int *l, int *ii, int *i1, int *i2, int *i3, double *t)
{

	int i__1, i__2;

	static int i__, j, k, itest;
	extern int dchcel_(int *, int *, int *,
					   int *, double *);

	itest = 0;
	if (*ii == 1)
	{
		if (*i1 == 1)
		{
			goto L100;
		}
		if (*i1 == 2)
		{
			if (*i2 >= 2 && *i2 <= 9 && *i3 >= 2 && *i3 <= 9)
			{
				i__ = 11;
				j = 11;
				k = 11;
				dchcel_(l, &i__, &j, &k, t);
			}
			goto L100;
		}
		i__ = *i1 - 2;
		i__1 = *i2 + 1;
		for (j = *i2 - 1; j <= i__1; ++j)
		{
			i__2 = *i3 + 1;
			for (k = *i3 - 1; k <= i__2; ++k)
			{
				if (j >= 1 && j <= 10 && k >= 1 && k <= 10)
				{
					dchcel_(l, &i__, &j, &k, t);
				}
			}
		}
	}
	if (*ii == 2)
	{
		if (*i1 == 10)
		{
			goto L100;
		}
		if (*i1 == 9)
		{
			if (*i2 >= 2 && *i2 <= 9 && *i3 >= 2 && *i3 <= 9)
			{
				i__ = 11;
				j = 11;
				k = 11;
				dchcel_(l, &i__, &j, &k, t);
			}
			goto L100;
		}
		i__ = *i1 + 2;
		i__1 = *i2 + 1;
		for (j = *i2 - 1; j <= i__1; ++j)
		{
			i__2 = *i3 + 1;
			for (k = *i3 - 1; k <= i__2; ++k)
			{
				if (j >= 1 && j <= 10 && k >= 1 && k <= 10)
				{
					dchcel_(l, &i__, &j, &k, t);
				}
			}
		}
	}
	if (*ii == 3)
	{
		if (*i2 == 1)
		{
			goto L100;
		}
		if (*i2 == 2)
		{
			if (*i1 >= 2 && *i1 <= 9 && *i3 >= 2 && *i3 <= 9)
			{
				i__ = 11;
				j = 11;
				k = 11;
				dchcel_(l, &i__, &j, &k, t);
			}
			goto L100;
		}
		j = *i2 - 2;
		i__1 = *i1 + 1;
		for (i__ = *i1 - 1; i__ <= i__1; ++i__)
		{
			i__2 = *i3 + 1;
			for (k = *i3 - 1; k <= i__2; ++k)
			{
				if (i__ >= 1 && i__ <= 10 && k >= 1 && k <= 10)
				{
					dchcel_(l, &i__, &j, &k, t);
				}
			}
		}
	}
	if (*ii == 4)
	{
		if (*i2 == 10)
		{
			goto L100;
		}
		if (*i2 == 9)
		{
			if (*i1 >= 2 && *i1 <= 9 && *i3 >= 2 && *i3 <= 9)
			{
				i__ = 11;
				j = 11;
				k = 11;
				dchcel_(l, &i__, &j, &k, t);
			}
			goto L100;
		}
		j = *i2 + 2;
		i__1 = *i1 + 1;
		for (i__ = *i1 - 1; i__ <= i__1; ++i__)
		{
			i__2 = *i3 + 1;
			for (k = *i3 - 1; k <= i__2; ++k)
			{
				if (i__ >= 1 && i__ <= 10 && k >= 1 && k <= 10)
				{
					dchcel_(l, &i__, &j, &k, t);
				}
			}
		}
	}
	if (*ii == 5)
	{
		if (*i3 == 1)
		{
			goto L100;
		}
		if (*i3 == 2)
		{
			if (*i1 >= 2 && *i1 <= 9 && *i2 >= 2 && *i2 <= 9)
			{
				i__ = 11;
				j = 11;
				k = 11;
				dchcel_(l, &i__, &j, &k, t);
			}
			goto L100;
		}
		k = *i3 - 2;
		i__1 = *i1 + 1;
		for (i__ = *i1 - 1; i__ <= i__1; ++i__)
		{
			i__2 = *i2 + 1;
			for (j = *i2 - 1; j <= i__2; ++j)
			{
				if (i__ >= 1 && i__ <= 10 && j >= 1 && j <= 10)
				{
					dchcel_(l, &i__, &j, &k, t);
				}
			}
		}
	}
	if (*ii == 6)
	{
		if (*i3 == 10)
		{
			goto L100;
		}
		if (*i3 == 9)
		{
			if (*i1 >= 2 && *i1 <= 9 && *i2 >= 2 && *i2 <= 9)
			{
				i__ = 11;
				j = 11;
				k = 11;
				dchcel_(l, &i__, &j, &k, t);
			}
			goto L100;
		}
		k = *i3 + 2;
		i__1 = *i1 + 1;
		for (i__ = *i1 - 1; i__ <= i__1; ++i__)
		{
			i__2 = *i2 + 1;
			for (j = *i2 - 1; j <= i__2; ++j)
			{
				if (i__ >= 1 && i__ <= 10 && j >= 1 && j <= 10)
				{
					dchcel_(l, &i__, &j, &k, t);
				}
			}
		}
	}
L100:
	return 0;
}

int dchin2_(int *l, int *ii, int *i1, int *i2, int *i3, double *t)
{

	int i__1, i__2;

	static int i__, j, k, ia, ib, ic;
	extern int dchcel_(int *, int *, int *,
					   int *, double *);

	if (*ii == 1)
	{
		i__ = *i1 - 2;
		if (i__ <= 0)
		{
			i__ += 10;
		}
		ia = i__;
		i__1 = *i2 + 1;
		for (j = *i2 - 1; j <= i__1; ++j)
		{
			i__2 = *i3 + 1;
			for (k = *i3 - 1; k <= i__2; ++k)
			{
				ib = j;
				ic = k;
				if (j == 0)
				{
					ib = 10;
				}
				if (j == 11)
				{
					ib = 1;
				}
				if (k >= 1 && k <= 10)
				{
					dchcel_(l, &ia, &ib, &ic, t);
				}
			}
		}
	}
	if (*ii == 2)
	{
		i__ = *i1 + 2;
		if (i__ >= 11)
		{
			i__ += -10;
		}
		ia = i__;
		i__1 = *i2 + 1;
		for (j = *i2 - 1; j <= i__1; ++j)
		{
			i__2 = *i3 + 1;
			for (k = *i3 - 1; k <= i__2; ++k)
			{
				ib = j;
				ic = k;
				if (j == 0)
				{
					ib = 10;
				}
				if (j == 11)
				{
					ib = 1;
				}
				if (k >= 1 && k <= 10)
				{
					dchcel_(l, &ia, &ib, &ic, t);
				}
			}
		}
	}
	if (*ii == 3)
	{
		j = *i2 - 2;
		if (j <= 0)
		{
			j += 10;
		}
		ib = j;
		i__1 = *i1 + 1;
		for (i__ = *i1 - 1; i__ <= i__1; ++i__)
		{
			i__2 = *i3 + 1;
			for (k = *i3 - 1; k <= i__2; ++k)
			{
				ia = i__;
				ic = k;
				if (i__ == 0)
				{
					ia = 10;
				}
				if (i__ == 11)
				{
					ia = 1;
				}
				if (k >= 1 && k <= 10)
				{
					dchcel_(l, &ia, &ib, &ic, t);
				}
			}
		}
	}
	if (*ii == 4)
	{
		j = *i2 + 2;
		if (j >= 11)
		{
			j += -10;
		}
		ib = j;
		i__1 = *i1 + 1;
		for (i__ = *i1 - 1; i__ <= i__1; ++i__)
		{
			i__2 = *i3 + 1;
			for (k = *i3 - 1; k <= i__2; ++k)
			{
				ia = i__;
				ic = k;
				if (i__ == 0)
				{
					ia = 10;
				}
				if (i__ == 11)
				{
					ia = 1;
				}
				if (k >= 1 && k <= 10)
				{
					dchcel_(l, &ia, &ib, &ic, t);
				}
			}
		}
	}
	if (*ii == 5)
	{
		if (*i3 == 2)
		{
			goto L100;
		}
		k = *i3 - 2;
		ic = k;
		i__1 = *i1 + 1;
		for (i__ = *i1 - 1; i__ <= i__1; ++i__)
		{
			i__2 = *i2 + 1;
			for (j = *i2 - 1; j <= i__2; ++j)
			{
				ia = i__;
				ib = j;
				if (i__ == 0)
				{
					ia = 10;
				}
				if (i__ == 11)
				{
					ia = 1;
				}
				if (j == 0)
				{
					ib = 10;
				}
				if (j == 11)
				{
					ib = 1;
				}
				dchcel_(l, &ia, &ib, &ic, t);
			}
		}
	}
	if (*ii == 6)
	{
		if (*i3 == 9)
		{
			goto L100;
		}
		k = *i3 + 2;
		ic = k;
		i__1 = *i1 + 1;
		for (i__ = *i1 - 1; i__ <= i__1; ++i__)
		{
			i__2 = *i2 + 1;
			for (j = *i2 - 1; j <= i__2; ++j)
			{
				ia = i__;
				ib = j;
				if (i__ == 0)
				{
					ia = 10;
				}
				if (i__ == 11)
				{
					ia = 1;
				}
				if (j == 0)
				{
					ib = 10;
				}
				if (j == 11)
				{
					ib = 1;
				}
				dchcel_(l, &ia, &ib, &ic, t);
			}
		}
	}
L100:
	return 0;
}

int dchin3_(int *l, int *ii, int *i1, int *i2, int *i3, double *t)
{

	int i__1, i__2;

	static int i__, j, k, ia, ib, ic;
	extern int dchcel_(int *, int *, int *,
					   int *, double *);

	if (*ii == 1)
	{
		i__ = *i1 - 2;
		if (i__ <= 0)
		{
			i__ += 10;
		}
		ia = i__;
		i__1 = *i2 + 1;
		for (j = *i2 - 1; j <= i__1; ++j)
		{
			i__2 = *i3 + 1;
			for (k = *i3 - 1; k <= i__2; ++k)
			{
				ib = j;
				ic = k;
				if (j == 0)
				{
					ib = 10;
				}
				if (j == 11)
				{
					ib = 1;
				}
				if (k == 0)
				{
					ic = 10;
				}
				if (k == 11)
				{
					ic = 1;
				}
				dchcel_(l, &ia, &ib, &ic, t);
			}
		}
	}
	if (*ii == 2)
	{
		i__ = *i1 + 2;
		if (i__ >= 11)
		{
			i__ += -10;
		}
		ia = i__;
		i__1 = *i2 + 1;
		for (j = *i2 - 1; j <= i__1; ++j)
		{
			i__2 = *i3 + 1;
			for (k = *i3 - 1; k <= i__2; ++k)
			{
				ib = j;
				ic = k;
				if (j == 0)
				{
					ib = 10;
				}
				if (j == 11)
				{
					ib = 1;
				}
				if (k == 0)
				{
					ic = 10;
				}
				if (k == 11)
				{
					ic = 1;
				}
				dchcel_(l, &ia, &ib, &ic, t);
			}
		}
	}
	if (*ii == 3)
	{
		j = *i2 - 2;
		if (j <= 0)
		{
			j += 10;
		}
		ib = j;
		i__1 = *i1 + 1;
		for (i__ = *i1 - 1; i__ <= i__1; ++i__)
		{
			i__2 = *i3 + 1;
			for (k = *i3 - 1; k <= i__2; ++k)
			{
				ia = i__;
				ic = k;
				if (i__ == 0)
				{
					ia = 10;
				}
				if (i__ == 11)
				{
					ia = 1;
				}
				if (k == 0)
				{
					ic = 10;
				}
				if (k == 11)
				{
					ic = 1;
				}
				dchcel_(l, &ia, &ib, &ic, t);
			}
		}
	}
	if (*ii == 4)
	{
		j = *i2 + 2;
		if (j >= 11)
		{
			j += -10;
		}
		ib = j;
		i__1 = *i1 + 1;
		for (i__ = *i1 - 1; i__ <= i__1; ++i__)
		{
			i__2 = *i3 + 1;
			for (k = *i3 - 1; k <= i__2; ++k)
			{
				ia = i__;
				ic = k;
				if (i__ == 0)
				{
					ia = 10;
				}
				if (i__ == 11)
				{
					ia = 1;
				}
				if (k == 0)
				{
					ic = 10;
				}
				if (k == 11)
				{
					ic = 1;
				}
				dchcel_(l, &ia, &ib, &ic, t);
			}
		}
	}
	if (*ii == 5)
	{
		k = *i3 - 2;
		if (k <= 0)
		{
			k += 10;
		}
		ic = k;
		i__1 = *i1 + 1;
		for (i__ = *i1 - 1; i__ <= i__1; ++i__)
		{
			i__2 = *i2 + 1;
			for (j = *i2 - 1; j <= i__2; ++j)
			{
				ia = i__;
				ib = j;
				if (i__ == 0)
				{
					ia = 10;
				}
				if (i__ == 11)
				{
					ia = 1;
				}
				if (j == 0)
				{
					ib = 10;
				}
				if (j == 11)
				{
					ib = 1;
				}
				dchcel_(l, &ia, &ib, &ic, t);
			}
		}
	}
	if (*ii == 6)
	{
		k = *i3 + 2;
		if (k >= 11)
		{
			k += -10;
		}
		ic = k;
		i__1 = *i1 + 1;
		for (i__ = *i1 - 1; i__ <= i__1; ++i__)
		{
			i__2 = *i2 + 1;
			for (j = *i2 - 1; j <= i__2; ++j)
			{
				ia = i__;
				ib = j;
				if (i__ == 0)
				{
					ia = 10;
				}
				if (i__ == 11)
				{
					ia = 1;
				}
				if (j == 0)
				{
					ib = 10;
				}
				if (j == 11)
				{
					ib = 1;
				}
				dchcel_(l, &ia, &ib, &ic, t);
			}
		}
	}

	return 0;
}

int dchcel_(int *l, int *i__, int *j, int *k,
			double *t)
{

	int s_stop(char *, ftnlen);

	static int m, n;
	static double tm;
	extern int reor_(double *, double *, int *,
					 int *);
	static int last0;

	if (*i__ == 11 || *j == 11 || *k == 11)
	{
		if (!(*i__ == 11 && *j == 11 && *k == 11))
		{
			s_stop("cerr", (ftnlen)4);
		}
		m = ilist2_1.icell;
	}
	else
	{
		m = ilist2_1.icel[*i__ + (*j + *k * 10) * 10 - 111];
	}
	if (m == 0)
	{
		return 0;
	}
	if (ilist1_1.next[m - 1] == *l)
	{
		tm = ilist5_1.tlarge;
		last0 = 0;
		reor_(t, &tm, &m, &last0);
	}
	n = ilist1_1.nic[m - 1];
	if (n == 0)
	{
		return 0;
	}
	while (n != m)
	{
		if (ilist1_1.next[n - 1] == *l)
		{
			tm = ilist5_1.tlarge;
			last0 = 0;
			reor_(t, &tm, &n, &last0);
		}
		n = ilist1_1.nic[n - 1];
	}
	return 0;
}

int fixtim_(int *l, double *t, double *tmin1,
			double *tmin, int *nc)
{
	static int k;

	k = *nc;
	if (*tmin < *tmin1)
	{
		ilist5_1.ot[*l - 1] = *tmin;
		if (ilist5_1.ct[*l - 1] < *tmin1)
		{
			ilist1_1.icsta[*l - 1] = 0;
		}
		else
		{
			ilist1_1.icsta[*l - 1] += 10;
		}
		ilist1_1.next[*l - 1] = k;
	}
	else if (*tmin == *tmin1)
	{
		ilist5_1.ot[*l - 1] = *tmin;
		if (*nc == 0)
		{
			ilist1_1.next[*l - 1] = 0;
		}
		else
		{
			ilist1_1.icsta[*l - 1] += 10;
			ilist1_1.next[*l - 1] = k;
		}
	}
	else
	{
		ilist5_1.ot[*l - 1] = *tmin1;
		ilist1_1.next[*l - 1] = 0;
	}
	return 0;
}

int ud2_(int *i__, int *j, double *t, double *tmin, int *nc)
{
	static int i1, i2, i3;
	static double t1, t2, tm;
	extern int isco_(int *, int *, logical *,
					 double *, double *, double *),
		reor_(double *,
			  double *, int *, int *);
	static double tmin1;
	static logical allok;
	extern int wallc_(int *, int *, int *,
					  int *, double *, double *);
	static int icels0;
	extern int fixtim_(int *, double *, double *,
					   double *, int *);

	isco_(i__, j, &allok, &tm, &t1, &t2);
	if (allok)
	{

		if (tm < *tmin)
		{
			*tmin = tm;
			ilist5_1.ct[*j - 1] = t2;
			*nc = *i__;
			if (para5_1.iconfg == 3 || para5_1.iconfg == 5)
			{
				aurec2_1.dgxa[*j - 1] = aurec1_1.jxa * 10. * ilist3_1.size1;
				aurec2_1.dgya[*j - 1] = aurec1_1.jya * 10. * ilist3_1.size2;
				if (para5_1.iconfg == 5)
				{
					aurec2_1.dgza[*j - 1] = aurec1_1.jza * 10. *
											ilist3_1.size3;
				}
			}
		}
		if (tm <= ilist5_1.ot[*i__ - 1])
		{
			ilist5_1.ct[*i__ - 1] = t1;
			icels0 = ilist1_1.icels[*i__ - 1];
			i1 = icels0 / 10000;
			i2 = (icels0 - i1 * 10000) / 100;
			i3 = icels0 - i1 * 10000 - i2 * 100;
			wallc_(i__, &i1, &i2, &i3, t, &tmin1);
			fixtim_(i__, t, &tmin1, &tm, j);
			if (para5_1.iconfg == 3 || para5_1.iconfg == 5)
			{
				aurec2_1.dgxa[*i__ - 1] = -aurec1_1.jxa * 10. *
										  ilist3_1.size1;
				aurec2_1.dgya[*i__ - 1] = -aurec1_1.jya * 10. *
										  ilist3_1.size2;
				if (para5_1.iconfg == 5)
				{
					aurec2_1.dgza[*i__ - 1] = -aurec1_1.jza * 10. *
											  ilist3_1.size3;
				}
			}
		}
		if (tm > ilist5_1.ot[*i__ - 1] && ilist1_1.next[*i__ - 1] == *j)
		{
			ilist5_1.ct[*i__ - 1] = t1;
			reor_(t, &tm, i__, j);
		}
	}
	else if (ilist1_1.next[*i__ - 1] == *j)
	{
		tm = ilist5_1.tlarge;
		reor_(t, &tm, i__, j);
	}
	return 0;
}

int isco_(int *i__, int *j, logical *allok,
		  double *tm, double *t1, double *t2)
{
	extern int isco1_(int *, int *, logical *,
					  double *, double *, double *),
		isco2_(int *,
			   int *, logical *, double *, double *, double *),
		isco3_(int *, int *, logical *, double *, double *, double *), isco4_(int *, int *, logical *, double *, double *, double *), isco5_(int *, int *, logical *, double *, double *, double *),
		isco6_(int *, int *, logical *, double *, double *, double *), isco7_(int *, int *, logical *, double *, double *, double *), isco8_(int *, int *, logical *, double *, double *, double *),
		isco9_(int *, int *, logical *, double *, double *, double *), isco10_(int *, int *, logical *, double *, double *, double *), isco11_(int *, int *, logical *, double *, double *, double *),
		isco12_(int *, int *, logical *, double *, double *, double *);
	static int iorder;

	iorder = para5_1.iordsc / 10;
	if (para5_1.iconfg == 1)
	{
		if (iorder == 1)
		{
			isco1_(i__, j, allok, tm, t1, t2);
		}
		else if (iorder == 2)
		{
			isco2_(i__, j, allok, tm, t1, t2);
		}
		else if (iorder == 3)
		{
			isco3_(i__, j, allok, tm, t1, t2);
		}
	}
	else if (para5_1.iconfg == 2 || para5_1.iconfg == 4)
	{
		if (iorder == 1)
		{
			isco4_(i__, j, allok, tm, t1, t2);
		}
		else if (iorder == 2)
		{
			isco5_(i__, j, allok, tm, t1, t2);
		}
		else if (iorder == 3)
		{
			isco6_(i__, j, allok, tm, t1, t2);
		}
	}
	else if (para5_1.iconfg == 3)
	{
		if (iorder == 1)
		{
			isco7_(i__, j, allok, tm, t1, t2);
		}
		else if (iorder == 2)
		{
			isco8_(i__, j, allok, tm, t1, t2);
		}
		else if (iorder == 3)
		{
			isco9_(i__, j, allok, tm, t1, t2);
		}
	}
	else if (para5_1.iconfg == 5)
	{
		if (iorder == 1)
		{
			isco10_(i__, j, allok, tm, t1, t2);
		}
		else if (iorder == 2)
		{
			isco11_(i__, j, allok, tm, t1, t2);
		}
		else if (iorder == 3)
		{
			isco12_(i__, j, allok, tm, t1, t2);
		}
	}
	return 0;
}

int isco1_(int *i__, int *j, logical *allok,
		   double *tm, double *t1, double *t2)
{

	double d__1, d__2, d__3, d__4;

	static double a, b, c__, d__, f, g, h__;
	static int i1, i2;
	static double p1, p2, p3, p4, q4, q1, q2, q3, r4, r1, r2, r3, ee, vp,
		dm2, tc1, tc2, rts2;

	*allok = ilist1_1.last[*i__ - 1] != *j || ilist1_1.last[*j - 1] != *i__;

	i1 = *i__;
	i2 = *j;
	p4 = prec2_1.ft[i2 - 1] - prec2_1.ft[i1 - 1];
	p1 = prec2_1.gx[i2 - 1] - prec2_1.gx[i1 - 1];
	p2 = prec2_1.gy[i2 - 1] - prec2_1.gy[i1 - 1];
	p3 = prec2_1.gz[i2 - 1] - prec2_1.gz[i1 - 1];
	q4 = prec2_1.e[i1 - 1];
	q1 = prec2_1.px[i1 - 1];
	q2 = prec2_1.py[i1 - 1];
	q3 = prec2_1.pz[i1 - 1];
	r4 = prec2_1.e[i2 - 1];
	r1 = prec2_1.px[i2 - 1];
	r2 = prec2_1.py[i2 - 1];
	r3 = prec2_1.pz[i2 - 1];
	a = p4 * q4 - p1 * q1 - p2 * q2 - p3 * q3;
	b = p4 * r4 - p1 * r1 - p2 * r2 - p3 * r3;
	c__ = q4 * q4 - q1 * q1 - q2 * q2 - q3 * q3;
	d__ = r4 * r4 - r1 * r1 - r2 * r2 - r3 * r3;
	ee = q4 * r4 - q1 * r1 - q2 * r2 - q3 * r3;
	f = p4 * p4 - p1 * p1 - p2 * p2 - p3 * p3;

	h__ = a + b;
	if (h__ > 0.)
	{
		g = a;
		a = -b;
		b = -g;
		g = c__;
		c__ = d__;
		d__ = g;
		i1 = *j;
		i2 = *i__;
	}

	if (*allok)
	{
		vp = a * d__ - b * ee;
		*allok = *allok && vp < 0.;
	}

	if (*allok)
	{

		d__1 = a;

		d__2 = b;

		d__3 = ee;
		dm2 = -f - (d__1 * d__1 * d__ + d__2 * d__2 * c__ - a * 2. * b * ee) /
					   (d__3 * d__3 - c__ * d__);
		*allok = *allok && dm2 < para2_1.cutof2;
	}

	if (*allok)
	{

		d__1 = ee;
		tc1 = prec2_1.ft[i1 - 1] - prec2_1.e[i1 - 1] * (a * d__ - b * ee) / (d__1 * d__1 - c__ * d__);

		d__1 = ee;
		tc2 = prec2_1.ft[i2 - 1] + prec2_1.e[i2 - 1] * (b * c__ - a * ee) / (d__1 * d__1 - c__ * d__);
		*tm = (tc1 + tc2) * .5;
		*allok = *allok && *tm > prec2_1.ft[*i__ - 1] && *tm > prec2_1.ft[*j - 1];
	}

	if (*allok)
	{

		d__1 = q4 + r4;

		d__2 = q1 + r1;

		d__3 = q2 + r2;

		d__4 = q3 + r3;
		rts2 = d__1 * d__1 - d__2 * d__2 - d__3 * d__3 - d__4 * d__4;
		*allok = *allok && rts2 > para2_1.rscut2;
	}
	if (!(*allok))
	{
		*tm = ilist5_1.tlarge;
		*t1 = ilist5_1.tlarge;
		*t2 = ilist5_1.tlarge;
	}
	else if (h__ > 0.)
	{
		*t1 = *tm;
		*t2 = *tm;
	}
	else
	{
		*t1 = *tm;
		*t2 = *tm;
	}
	return 0;
}

int isco2_(int *i__, int *j, logical *allok,
		   double *tm, double *t1, double *t2)
{

	double d__1, d__2, d__3, d__4;

	static double a, b, c__, d__, f, g, h__;
	static int i1, i2;
	static double p1, p2, p3, p4, q4, q1, q2, q3, r4, r1, r2, r3, ee, vp,
		dm2, tc1, tc2, rts2;

	*allok = ilist1_1.last[*i__ - 1] != *j || ilist1_1.last[*j - 1] != *i__;

	i1 = *i__;
	i2 = *j;
	p4 = prec2_1.ft[i2 - 1] - prec2_1.ft[i1 - 1];
	p1 = prec2_1.gx[i2 - 1] - prec2_1.gx[i1 - 1];
	p2 = prec2_1.gy[i2 - 1] - prec2_1.gy[i1 - 1];
	p3 = prec2_1.gz[i2 - 1] - prec2_1.gz[i1 - 1];
	q4 = prec2_1.e[i1 - 1];
	q1 = prec2_1.px[i1 - 1];
	q2 = prec2_1.py[i1 - 1];
	q3 = prec2_1.pz[i1 - 1];
	r4 = prec2_1.e[i2 - 1];
	r1 = prec2_1.px[i2 - 1];
	r2 = prec2_1.py[i2 - 1];
	r3 = prec2_1.pz[i2 - 1];
	a = p4 * q4 - p1 * q1 - p2 * q2 - p3 * q3;
	b = p4 * r4 - p1 * r1 - p2 * r2 - p3 * r3;
	c__ = q4 * q4 - q1 * q1 - q2 * q2 - q3 * q3;
	d__ = r4 * r4 - r1 * r1 - r2 * r2 - r3 * r3;
	ee = q4 * r4 - q1 * r1 - q2 * r2 - q3 * r3;
	f = p4 * p4 - p1 * p1 - p2 * p2 - p3 * p3;

	h__ = a + b;
	if (h__ > 0.)
	{
		g = a;
		a = -b;
		b = -g;
		g = c__;
		c__ = d__;
		d__ = g;
		i1 = *j;
		i2 = *i__;
	}

	if (*allok)
	{
		vp = a * d__ - b * ee;
		*allok = *allok && vp < 0.;
	}

	if (*allok)
	{

		d__1 = a;

		d__2 = b;

		d__3 = ee;
		dm2 = -f - (d__1 * d__1 * d__ + d__2 * d__2 * c__ - a * 2. * b * ee) /
					   (d__3 * d__3 - c__ * d__);
		*allok = *allok && dm2 < para2_1.cutof2;
	}

	if (*allok)
	{

		d__1 = ee;
		tc1 = prec2_1.ft[i1 - 1] - prec2_1.e[i1 - 1] * (a * d__ - b * ee) / (d__1 * d__1 - c__ * d__);

		d__1 = ee;
		tc2 = prec2_1.ft[i2 - 1] + prec2_1.e[i2 - 1] * (b * c__ - a * ee) / (d__1 * d__1 - c__ * d__);
		if (para5_1.iordsc == 20)
		{
			*tm = min(tc1, tc2);
		}
		else if (para5_1.iordsc == 21)
		{
			*tm = (tc1 + tc2) * .5;
		}
		else
		{
			*tm = max(tc1, tc2);
		}
		*allok = *allok && *tm > prec2_1.ft[*i__ - 1] && *tm > prec2_1.ft[*j - 1];
	}

	if (*allok)
	{

		d__1 = q4 + r4;

		d__2 = q1 + r1;

		d__3 = q2 + r2;

		d__4 = q3 + r3;
		rts2 = d__1 * d__1 - d__2 * d__2 - d__3 * d__3 - d__4 * d__4;
		*allok = *allok && rts2 > para2_1.rscut2;
	}
	if (!(*allok))
	{
		*tm = ilist5_1.tlarge;
		*t1 = ilist5_1.tlarge;
		*t2 = ilist5_1.tlarge;
	}
	else if (h__ > 0.)
	{
		*t1 = tc2;
		*t2 = tc1;
	}
	else
	{
		*t1 = tc1;
		*t2 = tc2;
	}
	return 0;
}

int isco3_(int *i__, int *j, logical *allok,
		   double *tm, double *t1, double *t2)
{

	double d__1, d__2, d__3, d__4;

	static double e1, e2;
	static int i1, i2;
	static double v2, dt, dx, dy, dz, vp, dm2, px1, py1, pz1, px2, py2,
		pz2, vx1, vy1, vz1, dgx, dgy, dgz, dvx, dvy, dvz, rts2;

	*allok = ilist1_1.last[*i__ - 1] != *j || ilist1_1.last[*j - 1] != *i__;
	if (prec2_1.ft[*i__ - 1] >= prec2_1.ft[*j - 1])
	{
		i1 = *j;
		i2 = *i__;
	}
	else
	{
		i1 = *i__;
		i2 = *j;
	}
	if (*allok)
	{
		*t1 = prec2_1.ft[i1 - 1];
		vx1 = prec4_1.vx[i1 - 1];
		vy1 = prec4_1.vy[i1 - 1];
		vz1 = prec4_1.vz[i1 - 1];
		*t2 = prec2_1.ft[i2 - 1];
		dvx = prec4_1.vx[i2 - 1] - vx1;
		dvy = prec4_1.vy[i2 - 1] - vy1;
		dvz = prec4_1.vz[i2 - 1] - vz1;
		dt = *t2 - *t1;
		dx = prec2_1.gx[i2 - 1] - prec2_1.gx[i1 - 1] - vx1 * dt;
		dy = prec2_1.gy[i2 - 1] - prec2_1.gy[i1 - 1] - vy1 * dt;
		dz = prec2_1.gz[i2 - 1] - prec2_1.gz[i1 - 1] - vz1 * dt;
		vp = dvx * dx + dvy * dy + dvz * dz;
		*allok = *allok && vp < 0.;
	}
	if (*allok)
	{
		v2 = dvx * dvx + dvy * dvy + dvz * dvz;
		if (v2 == 0.)
		{
			*tm = ilist5_1.tlarge;
		}
		else
		{
			*tm = *t2 - vp / v2;
		}

		*allok = *allok && *tm > *t1 && *tm > *t2;
	}
	if (*allok)
	{
		dgx = dx - dvx * *t2;
		dgy = dy - dvy * *t2;
		dgz = dz - dvz * *t2;

		d__1 = *tm;
		dm2 = -v2 * (d__1 * d__1) + dgx * dgx + dgy * dgy + dgz * dgz;
		*allok = *allok && dm2 < para2_1.cutof2;
	}
	if (*allok)
	{
		e1 = prec2_1.e[i1 - 1];
		px1 = prec2_1.px[i1 - 1];
		py1 = prec2_1.py[i1 - 1];
		pz1 = prec2_1.pz[i1 - 1];
		e2 = prec2_1.e[i2 - 1];
		px2 = prec2_1.px[i2 - 1];
		py2 = prec2_1.py[i2 - 1];
		pz2 = prec2_1.pz[i2 - 1];

		d__1 = e1 + e2;

		d__2 = px1 + px2;

		d__3 = py1 + py2;

		d__4 = pz1 + pz2;
		rts2 = d__1 * d__1 - d__2 * d__2 - d__3 * d__3 - d__4 * d__4;
		*allok = *allok && rts2 > para2_1.rscut2;
	}
	if (!(*allok))
	{
		*tm = ilist5_1.tlarge;
		*t1 = ilist5_1.tlarge;
		*t2 = ilist5_1.tlarge;
	}
	else
	{
		*t1 = *tm;
		*t2 = *tm;
	}
	return 0;
}

int isco4_(int *i__, int *j, logical *allok,
		   double *tm, double *t1, double *t2)
{

	double d__1, d__2, d__3, d__4;

	static double a, b, c__, d__, f, g, h__;
	static int i1, i2;
	static double p1, p2, p3, p4, q4, q1, q2, q3, r4, r1, r2, r3, ee, vp;
	static int ii1, ii2, jj1, jj2, kk1, kk2;
	static double dm2, tc1, tc2, rts2;
	static int icels1, icels2;

	*allok = ilist1_1.last[*i__ - 1] != *j || ilist1_1.last[*j - 1] != *i__;

	icels1 = ilist1_1.icels[*i__ - 1];
	ii1 = icels1 / 10000;
	jj1 = (icels1 - ii1 * 10000) / 100;
	kk1 = icels1 - ii1 * 10000 - jj1 * 100;
	icels2 = ilist1_1.icels[*j - 1];
	ii2 = icels2 / 10000;
	jj2 = (icels2 - ii2 * 10000) / 100;
	kk2 = icels2 - ii2 * 10000 - jj2 * 100;
	i1 = *i__;
	i2 = *j;
	p4 = prec2_1.ft[i2 - 1] - prec2_1.ft[i1 - 1];
	p1 = prec2_1.gx[i2 - 1] - prec2_1.gx[i1 - 1];
	p2 = prec2_1.gy[i2 - 1] - prec2_1.gy[i1 - 1];
	p3 = prec2_1.gz[i2 - 1] - prec2_1.gz[i1 - 1];
	if (ii1 - ii2 > 5)
	{
		p1 += ilist3_1.size1 * 10.;
	}
	else if (ii1 - ii2 < -5)
	{
		p1 -= ilist3_1.size1 * 10.;
	}
	if (jj1 - jj2 > 5)
	{
		p2 += ilist3_1.size2 * 10.;
	}
	else if (jj1 - jj2 < -5)
	{
		p2 -= ilist3_1.size2 * 10.;
	}
	if (kk1 - kk2 > 5)
	{
		p3 += ilist3_1.size3 * 10.;
	}
	else if (kk1 - kk2 < -5)
	{
		p3 -= ilist3_1.size3 * 10.;
	}
	q4 = prec2_1.e[i1 - 1];
	q1 = prec2_1.px[i1 - 1];
	q2 = prec2_1.py[i1 - 1];
	q3 = prec2_1.pz[i1 - 1];
	r4 = prec2_1.e[i2 - 1];
	r1 = prec2_1.px[i2 - 1];
	r2 = prec2_1.py[i2 - 1];
	r3 = prec2_1.pz[i2 - 1];
	a = p4 * q4 - p1 * q1 - p2 * q2 - p3 * q3;
	b = p4 * r4 - p1 * r1 - p2 * r2 - p3 * r3;
	c__ = q4 * q4 - q1 * q1 - q2 * q2 - q3 * q3;
	d__ = r4 * r4 - r1 * r1 - r2 * r2 - r3 * r3;
	ee = q4 * r4 - q1 * r1 - q2 * r2 - q3 * r3;
	f = p4 * p4 - p1 * p1 - p2 * p2 - p3 * p3;

	h__ = a + b;
	if (h__ > 0.)
	{
		g = a;
		a = -b;
		b = -g;
		g = c__;
		c__ = d__;
		d__ = g;
		i1 = *j;
		i2 = *i__;
	}

	if (*allok)
	{
		vp = a * d__ - b * ee;
		*allok = *allok && vp < 0.;
	}

	if (*allok)
	{

		d__1 = a;

		d__2 = b;

		d__3 = ee;
		dm2 = -f - (d__1 * d__1 * d__ + d__2 * d__2 * c__ - a * 2. * b * ee) /
					   (d__3 * d__3 - c__ * d__);
		*allok = *allok && dm2 < para2_1.cutof2;
	}

	if (*allok)
	{

		d__1 = ee;
		tc1 = prec2_1.ft[i1 - 1] - prec2_1.e[i1 - 1] * (a * d__ - b * ee) / (d__1 * d__1 - c__ * d__);

		d__1 = ee;
		tc2 = prec2_1.ft[i2 - 1] + prec2_1.e[i2 - 1] * (b * c__ - a * ee) / (d__1 * d__1 - c__ * d__);
		*tm = (tc1 + tc2) * .5;
		*allok = *allok && *tm > prec2_1.ft[*i__ - 1] && *tm > prec2_1.ft[*j - 1];
	}

	if (*allok)
	{

		d__1 = q4 + r4;

		d__2 = q1 + r1;

		d__3 = q2 + r2;

		d__4 = q3 + r3;
		rts2 = d__1 * d__1 - d__2 * d__2 - d__3 * d__3 - d__4 * d__4;
		*allok = *allok && rts2 > para2_1.rscut2;
	}
	if (!(*allok))
	{
		*tm = ilist5_1.tlarge;
		*t1 = ilist5_1.tlarge;
		*t2 = ilist5_1.tlarge;
	}
	else if (h__ > 0.)
	{
		*t1 = *tm;
		*t2 = *tm;
	}
	else
	{
		*t1 = *tm;
		*t2 = *tm;
	}
	return 0;
}

int isco5_(int *i__, int *j, logical *allok,
		   double *tm, double *t1, double *t2)
{

	double d__1, d__2, d__3, d__4;

	static double a, b, c__, d__, f, g, h__;
	static int i1, i2;
	static double p1, p2, p3, p4, q4, q1, q2, q3, r4, r1, r2, r3, ee, vp;
	static int ii1, ii2, jj1, jj2, kk1, kk2;
	static double dm2, tc1, tc2, rts2;
	static int icels1, icels2;

	*allok = ilist1_1.last[*i__ - 1] != *j || ilist1_1.last[*j - 1] != *i__;

	icels1 = ilist1_1.icels[*i__ - 1];
	ii1 = icels1 / 10000;
	jj1 = (icels1 - ii1 * 10000) / 100;
	kk1 = icels1 - ii1 * 10000 - jj1 * 100;
	icels2 = ilist1_1.icels[*j - 1];
	ii2 = icels2 / 10000;
	jj2 = (icels2 - ii2 * 10000) / 100;
	kk2 = icels2 - ii2 * 10000 - jj2 * 100;
	i1 = *i__;
	i2 = *j;
	p4 = prec2_1.ft[i2 - 1] - prec2_1.ft[i1 - 1];
	p1 = prec2_1.gx[i2 - 1] - prec2_1.gx[i1 - 1];
	p2 = prec2_1.gy[i2 - 1] - prec2_1.gy[i1 - 1];
	p3 = prec2_1.gz[i2 - 1] - prec2_1.gz[i1 - 1];
	if (ii1 - ii2 > 5)
	{
		p1 += ilist3_1.size1 * 10.;
	}
	else if (ii1 - ii2 < -5)
	{
		p1 -= ilist3_1.size1 * 10.;
	}
	if (jj1 - jj2 > 5)
	{
		p2 += ilist3_1.size2 * 10.;
	}
	else if (jj1 - jj2 < -5)
	{
		p2 -= ilist3_1.size2 * 10.;
	}
	if (kk1 - kk2 > 5)
	{
		p3 += ilist3_1.size3 * 10.;
	}
	else if (kk1 - kk2 < -5)
	{
		p3 -= ilist3_1.size3 * 10.;
	}
	q4 = prec2_1.e[i1 - 1];
	q1 = prec2_1.px[i1 - 1];
	q2 = prec2_1.py[i1 - 1];
	q3 = prec2_1.pz[i1 - 1];
	r4 = prec2_1.e[i2 - 1];
	r1 = prec2_1.px[i2 - 1];
	r2 = prec2_1.py[i2 - 1];
	r3 = prec2_1.pz[i2 - 1];
	a = p4 * q4 - p1 * q1 - p2 * q2 - p3 * q3;
	b = p4 * r4 - p1 * r1 - p2 * r2 - p3 * r3;
	c__ = q4 * q4 - q1 * q1 - q2 * q2 - q3 * q3;
	d__ = r4 * r4 - r1 * r1 - r2 * r2 - r3 * r3;
	ee = q4 * r4 - q1 * r1 - q2 * r2 - q3 * r3;
	f = p4 * p4 - p1 * p1 - p2 * p2 - p3 * p3;

	h__ = a + b;
	if (h__ > 0.)
	{
		g = a;
		a = -b;
		b = -g;
		g = c__;
		c__ = d__;
		d__ = g;
		i1 = *j;
		i2 = *i__;
	}

	if (*allok)
	{
		vp = a * d__ - b * ee;
		*allok = *allok && vp < 0.;
	}

	if (*allok)
	{

		d__1 = a;

		d__2 = b;

		d__3 = ee;
		dm2 = -f - (d__1 * d__1 * d__ + d__2 * d__2 * c__ - a * 2. * b * ee) /
					   (d__3 * d__3 - c__ * d__);
		*allok = *allok && dm2 < para2_1.cutof2;
	}

	if (*allok)
	{

		d__1 = ee;
		tc1 = prec2_1.ft[i1 - 1] - prec2_1.e[i1 - 1] * (a * d__ - b * ee) / (d__1 * d__1 - c__ * d__);

		d__1 = ee;
		tc2 = prec2_1.ft[i2 - 1] + prec2_1.e[i2 - 1] * (b * c__ - a * ee) / (d__1 * d__1 - c__ * d__);
		if (para5_1.iordsc == 20)
		{
			*tm = min(tc1, tc2);
		}
		else if (para5_1.iordsc == 21)
		{
			*tm = (tc1 + tc2) * .5;
		}
		else
		{
			*tm = max(tc1, tc2);
		}
		*allok = *allok && *tm > prec2_1.ft[*i__ - 1] && *tm > prec2_1.ft[*j - 1];
	}

	if (*allok)
	{

		d__1 = q4 + r4;

		d__2 = q1 + r1;

		d__3 = q2 + r2;

		d__4 = q3 + r3;
		rts2 = d__1 * d__1 - d__2 * d__2 - d__3 * d__3 - d__4 * d__4;
		*allok = *allok && rts2 > para2_1.rscut2;
	}
	if (!(*allok))
	{
		*tm = ilist5_1.tlarge;
		*t1 = ilist5_1.tlarge;
		*t2 = ilist5_1.tlarge;
	}
	else if (h__ > 0.)
	{
		*t1 = tc2;
		*t2 = tc1;
	}
	else
	{
		*t1 = tc1;
		*t2 = tc2;
	}
	return 0;
}

int isco6_(int *i__, int *j, logical *allok,
		   double *tm, double *t1, double *t2)
{

	double d__1, d__2, d__3, d__4;

	static double e1, e2;
	static int i1, i2;
	static double dt, dx, dy, dz, vp;
	static int ii1, ii2, jj1, jj2, kk1, kk2;
	static double dm2, v2p, px1, py1, pz1, px2, py2, pz2, vx1, vy1, vz1,
		dgx, dgy, dgz, dvx, dvy, dvz, rts2;
	static int icels1, icels2;

	*allok = ilist1_1.last[*i__ - 1] != *j || ilist1_1.last[*j - 1] != *i__;
	if (prec2_1.ft[*i__ - 1] >= prec2_1.ft[*j - 1])
	{
		i1 = *j;
		i2 = *i__;
	}
	else
	{
		i1 = *i__;
		i2 = *j;
	}
	icels1 = ilist1_1.icels[i1 - 1];
	ii1 = icels1 / 10000;
	jj1 = (icels1 - ii1 * 10000) / 100;
	kk1 = icels1 - ii1 * 10000 - jj1 * 100;
	icels2 = ilist1_1.icels[i2 - 1];
	ii2 = icels2 / 10000;
	jj2 = (icels2 - ii2 * 10000) / 100;
	kk2 = icels2 - ii2 * 10000 - jj2 * 100;
	if (*allok)
	{
		*t1 = prec2_1.ft[i1 - 1];
		vx1 = prec4_1.vx[i1 - 1];
		vy1 = prec4_1.vy[i1 - 1];
		vz1 = prec4_1.vz[i1 - 1];
		*t2 = prec2_1.ft[i2 - 1];
		dvx = prec4_1.vx[i2 - 1] - vx1;
		dvy = prec4_1.vy[i2 - 1] - vy1;
		dvz = prec4_1.vz[i2 - 1] - vz1;
		dt = *t2 - *t1;
		dx = prec2_1.gx[i2 - 1] - prec2_1.gx[i1 - 1] - vx1 * dt;
		dy = prec2_1.gy[i2 - 1] - prec2_1.gy[i1 - 1] - vy1 * dt;
		dz = prec2_1.gz[i2 - 1] - prec2_1.gz[i1 - 1] - vz1 * dt;
		if (ii1 - ii2 > 5)
		{
			dx += ilist3_1.size1 * 10.;
		}
		else if (ii1 - ii2 < -5)
		{
			dx -= ilist3_1.size1 * 10.;
		}
		if (jj1 - jj2 > 5)
		{
			dy += ilist3_1.size2 * 10.;
		}
		else if (jj1 - jj2 < -5)
		{
			dy -= ilist3_1.size2 * 10.;
		}
		if (kk1 - kk2 > 5)
		{
			dz += ilist3_1.size3 * 10.;
		}
		else if (kk1 - kk2 < -5)
		{
			dz -= ilist3_1.size3 * 10.;
		}
		vp = dvx * dx + dvy * dy + dvz * dz;
		*allok = *allok && vp < 0.;
	}
	if (*allok)
	{
		v2p = dvx * dvx + dvy * dvy + dvz * dvz;
		if (v2p == 0.)
		{
			*tm = ilist5_1.tlarge;
		}
		else
		{
			*tm = *t2 - vp / v2p;
		}

		*allok = *allok && *tm > *t1 && *tm > *t2;
	}
	if (*allok)
	{
		dgx = dx - dvx * *t2;
		dgy = dy - dvy * *t2;
		dgz = dz - dvz * *t2;

		d__1 = *tm;
		dm2 = -v2p * (d__1 * d__1) + dgx * dgx + dgy * dgy + dgz * dgz;
		*allok = *allok && dm2 < para2_1.cutof2;
	}
	if (*allok)
	{
		e1 = prec2_1.e[i1 - 1];
		px1 = prec2_1.px[i1 - 1];
		py1 = prec2_1.py[i1 - 1];
		pz1 = prec2_1.pz[i1 - 1];
		e2 = prec2_1.e[i2 - 1];
		px2 = prec2_1.px[i2 - 1];
		py2 = prec2_1.py[i2 - 1];
		pz2 = prec2_1.pz[i2 - 1];

		d__1 = e1 + e2;

		d__2 = px1 + px2;

		d__3 = py1 + py2;

		d__4 = pz1 + pz2;
		rts2 = d__1 * d__1 - d__2 * d__2 - d__3 * d__3 - d__4 * d__4;
		*allok = *allok && rts2 > para2_1.rscut2;
	}
	if (!(*allok))
	{
		*tm = ilist5_1.tlarge;
		*t1 = ilist5_1.tlarge;
		*t2 = ilist5_1.tlarge;
	}
	else
	{
		*t1 = *tm;
		*t2 = *tm;
	}
	return 0;
}

int isco7_(int *i__, int *j, logical *allok,
		   double *tm, double *t1, double *t2)
{

	double d__1, d__2, d__3, d__4;

	static double a, b, c__, d__, f, g, h__;
	static int i1, i2;
	static double p1, p2, p3, p4, q4, q1, q2, q3, r4, r1, r2, r3, ee;
	static int ii, jj;
	static double vp, dm2, tc1, tc2, tmp, rts2;
	static logical allokp;

	*allok = ilist1_1.last[*i__ - 1] != *j || ilist1_1.last[*j - 1] != *i__;

	*tm = ilist5_1.tlarge;
	if (*allok)
	{
		for (ii = -1; ii <= 1; ++ii)
		{
			for (jj = -1; jj <= 1; ++jj)
			{
				allokp = TRUE_;
				i1 = *i__;
				i2 = *j;
				p4 = prec2_1.ft[*j - 1] - prec2_1.ft[*i__ - 1];
				p1 = prec2_1.gx[*j - 1] - prec2_1.gx[*i__ - 1];
				p2 = prec2_1.gy[*j - 1] - prec2_1.gy[*i__ - 1];
				p3 = prec2_1.gz[*j - 1] - prec2_1.gz[*i__ - 1];
				p1 += ii * 10. * ilist3_1.size1;
				p2 += jj * 10. * ilist3_1.size2;
				q4 = prec2_1.e[*i__ - 1];
				q1 = prec2_1.px[*i__ - 1];
				q2 = prec2_1.py[*i__ - 1];
				q3 = prec2_1.pz[*i__ - 1];
				r4 = prec2_1.e[*j - 1];
				r1 = prec2_1.px[*j - 1];
				r2 = prec2_1.py[*j - 1];
				r3 = prec2_1.pz[*j - 1];
				a = p4 * q4 - p1 * q1 - p2 * q2 - p3 * q3;
				b = p4 * r4 - p1 * r1 - p2 * r2 - p3 * r3;
				c__ = q4 * q4 - q1 * q1 - q2 * q2 - q3 * q3;
				d__ = r4 * r4 - r1 * r1 - r2 * r2 - r3 * r3;
				ee = q4 * r4 - q1 * r1 - q2 * r2 - q3 * r3;
				f = p4 * p4 - p1 * p1 - p2 * p2 - p3 * p3;

				h__ = a + b;
				if (h__ > 0.)
				{
					g = a;
					a = -b;
					b = -g;
					g = c__;
					c__ = d__;
					d__ = g;
					i1 = *j;
					i2 = *i__;
				}

				if (allokp)
				{
					vp = a * d__ - b * ee;
					allokp = allokp && vp < 0.;
				}

				if (allokp)
				{

					d__1 = a;

					d__2 = b;

					d__3 = ee;
					dm2 = -f - (d__1 * d__1 * d__ + d__2 * d__2 * c__ - a * 2. * b * ee) / (d__3 * d__3 - c__ * d__);
					allokp = allokp && dm2 < para2_1.cutof2;
				}

				if (allokp)
				{

					d__1 = ee;
					tc1 = prec2_1.ft[i1 - 1] - prec2_1.e[i1 - 1] * (a * d__ - b * ee) / (d__1 * d__1 - c__ * d__);

					d__1 = ee;
					tc2 = prec2_1.ft[i2 - 1] + prec2_1.e[i2 - 1] * (b * c__ - a * ee) / (d__1 * d__1 - c__ * d__);
					tmp = (tc1 + tc2) * .5;
					allokp = allokp && tmp > prec2_1.ft[*i__ - 1] && tmp > prec2_1.ft[*j - 1];
				}
				if (allokp && tmp < *tm)
				{
					*tm = tmp;
					aurec1_1.jxa = ii;
					aurec1_1.jya = jj;
				}
			}
		}
		if (*tm == ilist5_1.tlarge)
		{
			*allok = FALSE_;
		}
	}

	if (*allok)
	{
		q4 = prec2_1.e[i1 - 1];
		q1 = prec2_1.px[i1 - 1];
		q2 = prec2_1.py[i1 - 1];
		q3 = prec2_1.pz[i1 - 1];
		r4 = prec2_1.e[i2 - 1];
		r1 = prec2_1.px[i2 - 1];
		r2 = prec2_1.py[i2 - 1];
		r3 = prec2_1.pz[i2 - 1];

		d__1 = q4 + r4;

		d__2 = q1 + r1;

		d__3 = q2 + r2;

		d__4 = q3 + r3;
		rts2 = d__1 * d__1 - d__2 * d__2 - d__3 * d__3 - d__4 * d__4;
		*allok = *allok && rts2 > para2_1.rscut2;
	}
	if (!(*allok))
	{
		*tm = ilist5_1.tlarge;
		*t1 = ilist5_1.tlarge;
		*t2 = ilist5_1.tlarge;
	}
	else if (h__ > 0.)
	{
		*t1 = *tm;
		*t2 = *tm;
	}
	else
	{
		*t1 = *tm;
		*t2 = *tm;
	}
	return 0;
}

int isco8_(int *i__, int *j, logical *allok,
		   double *tm, double *t1, double *t2)
{

	double d__1, d__2, d__3, d__4;

	static double a, b, c__, d__, f, g, h__;
	static int i1, i2;
	static double p1, p2, p3, p4, q4, q1, q2, q3, r4, r1, r2, r3, ha, ee;
	static int ii, jj;
	static double vp, dm2, tc1, tc2, tmp, tc1a, tc2a, rts2;
	static logical allokp;

	*allok = ilist1_1.last[*i__ - 1] != *j || ilist1_1.last[*j - 1] != *i__;

	*tm = ilist5_1.tlarge;
	if (*allok)
	{
		for (ii = -1; ii <= 1; ++ii)
		{
			for (jj = -1; jj <= 1; ++jj)
			{
				allokp = TRUE_;
				i1 = *i__;
				i2 = *j;
				p4 = prec2_1.ft[*j - 1] - prec2_1.ft[*i__ - 1];
				p1 = prec2_1.gx[*j - 1] - prec2_1.gx[*i__ - 1];
				p2 = prec2_1.gy[*j - 1] - prec2_1.gy[*i__ - 1];
				p3 = prec2_1.gz[*j - 1] - prec2_1.gz[*i__ - 1];
				p1 += ii * 10. * ilist3_1.size1;
				p2 += jj * 10. * ilist3_1.size2;
				q4 = prec2_1.e[*i__ - 1];
				q1 = prec2_1.px[*i__ - 1];
				q2 = prec2_1.py[*i__ - 1];
				q3 = prec2_1.pz[*i__ - 1];
				r4 = prec2_1.e[*j - 1];
				r1 = prec2_1.px[*j - 1];
				r2 = prec2_1.py[*j - 1];
				r3 = prec2_1.pz[*j - 1];
				a = p4 * q4 - p1 * q1 - p2 * q2 - p3 * q3;
				b = p4 * r4 - p1 * r1 - p2 * r2 - p3 * r3;
				c__ = q4 * q4 - q1 * q1 - q2 * q2 - q3 * q3;
				d__ = r4 * r4 - r1 * r1 - r2 * r2 - r3 * r3;
				ee = q4 * r4 - q1 * r1 - q2 * r2 - q3 * r3;
				f = p4 * p4 - p1 * p1 - p2 * p2 - p3 * p3;

				h__ = a + b;
				if (h__ > 0.)
				{
					g = a;
					a = -b;
					b = -g;
					g = c__;
					c__ = d__;
					d__ = g;
					i1 = *j;
					i2 = *i__;
				}

				if (allokp)
				{
					vp = a * d__ - b * ee;
					allokp = allokp && vp < 0.;
				}

				if (allokp)
				{

					d__1 = a;

					d__2 = b;

					d__3 = ee;
					dm2 = -f - (d__1 * d__1 * d__ + d__2 * d__2 * c__ - a * 2. * b * ee) / (d__3 * d__3 - c__ * d__);
					allokp = allokp && dm2 < para2_1.cutof2;
				}

				if (allokp)
				{

					d__1 = ee;
					tc1 = prec2_1.ft[i1 - 1] - prec2_1.e[i1 - 1] * (a * d__ - b * ee) / (d__1 * d__1 - c__ * d__);

					d__1 = ee;
					tc2 = prec2_1.ft[i2 - 1] + prec2_1.e[i2 - 1] * (b * c__ - a * ee) / (d__1 * d__1 - c__ * d__);
					if (para5_1.iordsc == 20)
					{
						tmp = min(tc1, tc2);
					}
					else if (para5_1.iordsc == 21)
					{
						tmp = (tc1 + tc2) * .5;
					}
					else
					{
						tmp = max(tc1, tc2);
					}
					allokp = allokp && tmp > prec2_1.ft[*i__ - 1] && tmp > prec2_1.ft[*j - 1];
				}
				if (allokp && tmp < *tm)
				{
					*tm = tmp;
					aurec1_1.jxa = ii;
					aurec1_1.jya = jj;
					ha = h__;
					tc1a = tc1;
					tc2a = tc2;
				}
			}
		}
		if (*tm == ilist5_1.tlarge)
		{
			*allok = FALSE_;
		}
	}

	if (*allok)
	{
		q4 = prec2_1.e[i1 - 1];
		q1 = prec2_1.px[i1 - 1];
		q2 = prec2_1.py[i1 - 1];
		q3 = prec2_1.pz[i1 - 1];
		r4 = prec2_1.e[i2 - 1];
		r1 = prec2_1.px[i2 - 1];
		r2 = prec2_1.py[i2 - 1];
		r3 = prec2_1.pz[i2 - 1];

		d__1 = q4 + r4;

		d__2 = q1 + r1;

		d__3 = q2 + r2;

		d__4 = q3 + r3;
		rts2 = d__1 * d__1 - d__2 * d__2 - d__3 * d__3 - d__4 * d__4;
		*allok = *allok && rts2 > para2_1.rscut2;
	}
	if (!(*allok))
	{
		*tm = ilist5_1.tlarge;
		*t1 = ilist5_1.tlarge;
		*t2 = ilist5_1.tlarge;
	}
	else if (ha > 0.)
	{
		*t1 = tc2a;
		*t2 = tc1a;
	}
	else
	{
		*t1 = tc1a;
		*t2 = tc2a;
	}
	return 0;
}

int isco9_(int *i__, int *j, logical *allok,
		   double *tm, double *t1, double *t2)
{

	double d__1, d__2, d__3, d__4;

	static double e1, e2;
	static int i1, i2, ii, jj;
	static double dt, dx, dy, dz, vp, dm2, px1, py1, pz1, px2, py2, pz2,
		vx1, vy1, vz1, dgx, dgy, dgz, tmp, dvx, dvy, dvz, rts2;
	static int isign;
	static logical allokp;

	*allok = ilist1_1.last[*i__ - 1] != *j || ilist1_1.last[*j - 1] != *i__;
	if (prec2_1.ft[*i__ - 1] >= prec2_1.ft[*j - 1])
	{
		i1 = *j;
		i2 = *i__;
		isign = -1;
	}
	else
	{
		i1 = *i__;
		i2 = *j;
		isign = 1;
	}
	if (*allok)
	{
		*tm = ilist5_1.tlarge;
		*t1 = prec2_1.ft[i1 - 1];
		vx1 = prec4_1.vx[i1 - 1];
		vy1 = prec4_1.vy[i1 - 1];
		vz1 = prec4_1.vz[i1 - 1];
		*t2 = prec2_1.ft[i2 - 1];
		dvx = prec4_1.vx[i2 - 1] - vx1;
		dvy = prec4_1.vy[i2 - 1] - vy1;
		dvz = prec4_1.vz[i2 - 1] - vz1;
		dt = *t2 - *t1;
		for (ii = -1; ii <= 1; ++ii)
		{
			for (jj = -1; jj <= 1; ++jj)
			{
				allokp = TRUE_;
				dx = prec2_1.gx[i2 - 1] - prec2_1.gx[i1 - 1] - vx1 * dt;
				dy = prec2_1.gy[i2 - 1] - prec2_1.gy[i1 - 1] - vy1 * dt;
				dz = prec2_1.gz[i2 - 1] - prec2_1.gz[i1 - 1] - vz1 * dt;
				dx += ii * 10. * ilist3_1.size1;
				dy += jj * 10. * ilist3_1.size2;
				vp = dvx * dx + dvy * dy + dvz * dz;
				allokp = allokp && vp < 0.;
				if (allokp)
				{
					ilist3_1.v2 = dvx * dvx + dvy * dvy + dvz * dvz;
					if (ilist3_1.v2 == 0.)
					{
						tmp = ilist5_1.tlarge;
					}
					else
					{
						tmp = *t2 - vp / ilist3_1.v2;
					}

					allokp = allokp && tmp > *t1 && tmp > *t2;
				}
				if (allokp)
				{
					dgx = dx - dvx * *t2;
					dgy = dy - dvy * *t2;
					dgz = dz - dvz * *t2;

					d__1 = tmp;
					dm2 = -ilist3_1.v2 * (d__1 * d__1) + dgx * dgx + dgy * dgy + dgz * dgz;
					allokp = allokp && dm2 < para2_1.cutof2;
				}
				if (allokp && tmp < *tm)
				{
					*tm = tmp;
					aurec1_1.jxa = isign * ii;
					aurec1_1.jya = isign * jj;
				}
			}
		}
		if (*tm == ilist5_1.tlarge)
		{
			*allok = FALSE_;
		}
	}
	if (*allok)
	{
		e1 = prec2_1.e[i1 - 1];
		px1 = prec2_1.px[i1 - 1];
		py1 = prec2_1.py[i1 - 1];
		pz1 = prec2_1.pz[i1 - 1];
		e2 = prec2_1.e[i2 - 1];
		px2 = prec2_1.px[i2 - 1];
		py2 = prec2_1.py[i2 - 1];
		pz2 = prec2_1.pz[i2 - 1];

		d__1 = e1 + e2;

		d__2 = px1 + px2;

		d__3 = py1 + py2;

		d__4 = pz1 + pz2;
		rts2 = d__1 * d__1 - d__2 * d__2 - d__3 * d__3 - d__4 * d__4;
		*allok = *allok && rts2 > para2_1.rscut2;
	}
	if (!(*allok))
	{
		*tm = ilist5_1.tlarge;
		*t1 = ilist5_1.tlarge;
		*t2 = ilist5_1.tlarge;
	}
	else
	{
		*t1 = *tm;
		*t2 = *tm;
	}
	return 0;
}

int isco10_(int *i__, int *j, logical *allok,
			double *tm, double *t1, double *t2)
{

	double d__1, d__2, d__3, d__4;

	static double a, b, c__, d__, f, g, h__;
	static int i1, i2;
	static double p1, p2, p3, p4, q4, q1, q2, q3, r4, r1, r2, r3, ee;
	static int ii, jj, kk;
	static double vp, dm2, tc1, tc2, tmp, rts2;
	static logical allokp;

	*allok = ilist1_1.last[*i__ - 1] != *j || ilist1_1.last[*j - 1] != *i__;

	*tm = ilist5_1.tlarge;
	if (*allok)
	{
		for (ii = -1; ii <= 1; ++ii)
		{
			for (jj = -1; jj <= 1; ++jj)
			{
				for (kk = -1; kk <= 1; ++kk)
				{
					allokp = TRUE_;
					i1 = *i__;
					i2 = *j;
					p4 = prec2_1.ft[*j - 1] - prec2_1.ft[*i__ - 1];
					p1 = prec2_1.gx[*j - 1] - prec2_1.gx[*i__ - 1];
					p2 = prec2_1.gy[*j - 1] - prec2_1.gy[*i__ - 1];
					p3 = prec2_1.gz[*j - 1] - prec2_1.gz[*i__ - 1];
					p1 += ii * 10. * ilist3_1.size1;
					p2 += jj * 10. * ilist3_1.size2;
					p3 += kk * 10. * ilist3_1.size3;
					q4 = prec2_1.e[*i__ - 1];
					q1 = prec2_1.px[*i__ - 1];
					q2 = prec2_1.py[*i__ - 1];
					q3 = prec2_1.pz[*i__ - 1];
					r4 = prec2_1.e[*j - 1];
					r1 = prec2_1.px[*j - 1];
					r2 = prec2_1.py[*j - 1];
					r3 = prec2_1.pz[*j - 1];
					a = p4 * q4 - p1 * q1 - p2 * q2 - p3 * q3;
					b = p4 * r4 - p1 * r1 - p2 * r2 - p3 * r3;
					c__ = q4 * q4 - q1 * q1 - q2 * q2 - q3 * q3;
					d__ = r4 * r4 - r1 * r1 - r2 * r2 - r3 * r3;
					ee = q4 * r4 - q1 * r1 - q2 * r2 - q3 * r3;
					f = p4 * p4 - p1 * p1 - p2 * p2 - p3 * p3;

					h__ = a + b;
					if (h__ > 0.)
					{
						g = a;
						a = -b;
						b = -g;
						g = c__;
						c__ = d__;
						d__ = g;
						i1 = *j;
						i2 = *i__;
					}

					if (allokp)
					{
						vp = a * d__ - b * ee;
						allokp = allokp && vp < 0.;
					}

					if (allokp)
					{

						d__1 = a;

						d__2 = b;

						d__3 = ee;
						dm2 = -f - (d__1 * d__1 * d__ + d__2 * d__2 * c__ - a * 2. * b * ee) / (d__3 * d__3 - c__ * d__);
						allokp = allokp && dm2 < para2_1.cutof2;
					}

					if (allokp)
					{

						d__1 = ee;
						tc1 = prec2_1.ft[i1 - 1] - prec2_1.e[i1 - 1] * (a * d__ - b * ee) / (d__1 * d__1 - c__ * d__);

						d__1 = ee;
						tc2 = prec2_1.ft[i2 - 1] + prec2_1.e[i2 - 1] * (b * c__ - a * ee) / (d__1 * d__1 - c__ * d__);
						tmp = (tc1 + tc2) * .5;
						allokp = allokp && tmp > prec2_1.ft[*i__ - 1] && tmp > prec2_1.ft[*j - 1];
					}
					if (allokp && tmp < *tm)
					{
						*tm = tmp;
						aurec1_1.jxa = ii;
						aurec1_1.jya = jj;
						aurec1_1.jza = kk;
					}
				}
			}
		}
		if (*tm == ilist5_1.tlarge)
		{
			*allok = FALSE_;
		}
	}

	if (*allok)
	{
		q4 = prec2_1.e[i1 - 1];
		q1 = prec2_1.px[i1 - 1];
		q2 = prec2_1.py[i1 - 1];
		q3 = prec2_1.pz[i1 - 1];
		r4 = prec2_1.e[i2 - 1];
		r1 = prec2_1.px[i2 - 1];
		r2 = prec2_1.py[i2 - 1];
		r3 = prec2_1.pz[i2 - 1];

		d__1 = q4 + r4;

		d__2 = q1 + r1;

		d__3 = q2 + r2;

		d__4 = q3 + r3;
		rts2 = d__1 * d__1 - d__2 * d__2 - d__3 * d__3 - d__4 * d__4;
		*allok = *allok && rts2 > para2_1.rscut2;
	}
	if (!(*allok))
	{
		*tm = ilist5_1.tlarge;
		*t1 = ilist5_1.tlarge;
		*t2 = ilist5_1.tlarge;
	}
	else if (h__ > 0.)
	{
		*t1 = *tm;
		*t2 = *tm;
	}
	else
	{
		*t1 = *tm;
		*t2 = *tm;
	}
	return 0;
}

int isco11_(int *i__, int *j, logical *allok,
			double *tm, double *t1, double *t2)
{

	double d__1, d__2, d__3, d__4;

	static double a, b, c__, d__, f, g, h__;
	static int i1, i2;
	static double p1, p2, p3, p4, q4, q1, q2, q3, r4, r1, r2, r3, ha, ee;
	static int ii, jj, kk;
	static double vp, dm2, tc1, tc2, tmp, tc1a, tc2a, rts2;
	static logical allokp;

	*allok = ilist1_1.last[*i__ - 1] != *j || ilist1_1.last[*j - 1] != *i__;

	*tm = ilist5_1.tlarge;
	if (*allok)
	{
		for (ii = -1; ii <= 1; ++ii)
		{
			for (jj = -1; jj <= 1; ++jj)
			{
				for (kk = -1; kk <= 1; ++kk)
				{
					allokp = TRUE_;
					i1 = *i__;
					i2 = *j;
					p4 = prec2_1.ft[*j - 1] - prec2_1.ft[*i__ - 1];
					p1 = prec2_1.gx[*j - 1] - prec2_1.gx[*i__ - 1];
					p2 = prec2_1.gy[*j - 1] - prec2_1.gy[*i__ - 1];
					p3 = prec2_1.gz[*j - 1] - prec2_1.gz[*i__ - 1];
					p1 += ii * 10. * ilist3_1.size1;
					p2 += jj * 10. * ilist3_1.size2;
					p3 += kk * 10. * ilist3_1.size3;
					q4 = prec2_1.e[*i__ - 1];
					q1 = prec2_1.px[*i__ - 1];
					q2 = prec2_1.py[*i__ - 1];
					q3 = prec2_1.pz[*i__ - 1];
					r4 = prec2_1.e[*j - 1];
					r1 = prec2_1.px[*j - 1];
					r2 = prec2_1.py[*j - 1];
					r3 = prec2_1.pz[*j - 1];
					a = p4 * q4 - p1 * q1 - p2 * q2 - p3 * q3;
					b = p4 * r4 - p1 * r1 - p2 * r2 - p3 * r3;
					c__ = q4 * q4 - q1 * q1 - q2 * q2 - q3 * q3;
					d__ = r4 * r4 - r1 * r1 - r2 * r2 - r3 * r3;
					ee = q4 * r4 - q1 * r1 - q2 * r2 - q3 * r3;
					f = p4 * p4 - p1 * p1 - p2 * p2 - p3 * p3;

					h__ = a + b;
					if (h__ > 0.)
					{
						g = a;
						a = -b;
						b = -g;
						g = c__;
						c__ = d__;
						d__ = g;
						i1 = *j;
						i2 = *i__;
					}

					if (allokp)
					{
						vp = a * d__ - b * ee;
						allokp = allokp && vp < 0.;
					}

					if (allokp)
					{

						d__1 = a;

						d__2 = b;

						d__3 = ee;
						dm2 = -f - (d__1 * d__1 * d__ + d__2 * d__2 * c__ - a * 2. * b * ee) / (d__3 * d__3 - c__ * d__);
						allokp = allokp && dm2 < para2_1.cutof2;
					}

					if (allokp)
					{

						d__1 = ee;
						tc1 = prec2_1.ft[i1 - 1] - prec2_1.e[i1 - 1] * (a * d__ - b * ee) / (d__1 * d__1 - c__ * d__);

						d__1 = ee;
						tc2 = prec2_1.ft[i2 - 1] + prec2_1.e[i2 - 1] * (b * c__ - a * ee) / (d__1 * d__1 - c__ * d__);
						if (para5_1.iordsc == 20)
						{
							tmp = min(tc1, tc2);
						}
						else if (para5_1.iordsc == 21)
						{
							tmp = (tc1 + tc2) * .5;
						}
						else
						{
							tmp = max(tc1, tc2);
						}
						allokp = allokp && tmp > prec2_1.ft[*i__ - 1] && tmp > prec2_1.ft[*j - 1];
					}
					if (allokp && tmp < *tm)
					{
						*tm = tmp;
						aurec1_1.jxa = ii;
						aurec1_1.jya = jj;
						aurec1_1.jza = kk;
						ha = h__;
						tc1a = tc1;
						tc2a = tc2;
					}
				}
			}
		}
		if (*tm == ilist5_1.tlarge)
		{
			*allok = FALSE_;
		}
	}

	if (*allok)
	{
		q4 = prec2_1.e[i1 - 1];
		q1 = prec2_1.px[i1 - 1];
		q2 = prec2_1.py[i1 - 1];
		q3 = prec2_1.pz[i1 - 1];
		r4 = prec2_1.e[i2 - 1];
		r1 = prec2_1.px[i2 - 1];
		r2 = prec2_1.py[i2 - 1];
		r3 = prec2_1.pz[i2 - 1];

		d__1 = q4 + r4;

		d__2 = q1 + r1;

		d__3 = q2 + r2;

		d__4 = q3 + r3;
		rts2 = d__1 * d__1 - d__2 * d__2 - d__3 * d__3 - d__4 * d__4;
		*allok = *allok && rts2 > para2_1.rscut2;
	}
	if (!(*allok))
	{
		*tm = ilist5_1.tlarge;
		*t1 = ilist5_1.tlarge;
		*t2 = ilist5_1.tlarge;
	}
	else if (ha > 0.)
	{
		*t1 = tc2a;
		*t2 = tc1a;
	}
	else
	{
		*t1 = tc1a;
		*t2 = tc2a;
	}
	return 0;
}

int isco12_(int *i__, int *j, logical *allok,
			double *tm, double *t1, double *t2)
{

	double d__1, d__2, d__3, d__4;

	static double e1, e2;
	static int i1, i2, ii, jj, kk;
	static double dt, dx, dy, dz, vp, dm2, px1, py1, pz1, px2, py2, pz2,
		vx1, vy1, vz1, dgx, dgy, dgz, tmp, dvx, dvy, dvz, rts2;
	static int isign;
	static logical allokp;

	*allok = ilist1_1.last[*i__ - 1] != *j || ilist1_1.last[*j - 1] != *i__;
	if (prec2_1.ft[*i__ - 1] >= prec2_1.ft[*j - 1])
	{
		i1 = *j;
		i2 = *i__;
		isign = -1;
	}
	else
	{
		i1 = *i__;
		i2 = *j;
		isign = 1;
	}
	if (*allok)
	{
		*tm = ilist5_1.tlarge;
		*t1 = prec2_1.ft[i1 - 1];
		vx1 = prec4_1.vx[i1 - 1];
		vy1 = prec4_1.vy[i1 - 1];
		vz1 = prec4_1.vz[i1 - 1];
		*t2 = prec2_1.ft[i2 - 1];
		dvx = prec4_1.vx[i2 - 1] - vx1;
		dvy = prec4_1.vy[i2 - 1] - vy1;
		dvz = prec4_1.vz[i2 - 1] - vz1;
		dt = *t2 - *t1;
		for (ii = -1; ii <= 1; ++ii)
		{
			for (jj = -1; jj <= 1; ++jj)
			{
				for (kk = -1; kk <= 1; ++kk)
				{
					allokp = TRUE_;
					dx = prec2_1.gx[i2 - 1] - prec2_1.gx[i1 - 1] - vx1 * dt;
					dy = prec2_1.gy[i2 - 1] - prec2_1.gy[i1 - 1] - vy1 * dt;
					dz = prec2_1.gz[i2 - 1] - prec2_1.gz[i1 - 1] - vz1 * dt;
					dx += ii * 10. * ilist3_1.size1;
					dy += jj * 10. * ilist3_1.size2;
					dz += kk * 10. * ilist3_1.size3;
					vp = dvx * dx + dvy * dy + dvz * dz;
					allokp = allokp && vp < 0.;
					if (allokp)
					{
						ilist3_1.v2 = dvx * dvx + dvy * dvy + dvz * dvz;
						if (ilist3_1.v2 == 0.)
						{
							tmp = ilist5_1.tlarge;
						}
						else
						{
							tmp = *t2 - vp / ilist3_1.v2;
						}

						allokp = allokp && tmp > *t1 && tmp > *t2;
					}
					if (allokp)
					{
						dgx = dx - dvx * *t2;
						dgy = dy - dvy * *t2;
						dgz = dz - dvz * *t2;

						d__1 = tmp;
						dm2 = -ilist3_1.v2 * (d__1 * d__1) + dgx * dgx + dgy * dgy + dgz * dgz;
						allokp = allokp && dm2 < para2_1.cutof2;
					}
					if (allokp && tmp < *tm)
					{
						*tm = tmp;
						aurec1_1.jxa = isign * ii;
						aurec1_1.jya = isign * jj;
						aurec1_1.jza = isign * kk;
					}
				}
			}
		}
		if (*tm == ilist5_1.tlarge)
		{
			*allok = FALSE_;
		}
	}
	if (*allok)
	{
		e1 = prec2_1.e[i1 - 1];
		px1 = prec2_1.px[i1 - 1];
		py1 = prec2_1.py[i1 - 1];
		pz1 = prec2_1.pz[i1 - 1];
		e2 = prec2_1.e[i2 - 1];
		px2 = prec2_1.px[i2 - 1];
		py2 = prec2_1.py[i2 - 1];
		pz2 = prec2_1.pz[i2 - 1];

		d__1 = e1 + e2;

		d__2 = px1 + px2;

		d__3 = py1 + py2;

		d__4 = pz1 + pz2;
		rts2 = d__1 * d__1 - d__2 * d__2 - d__3 * d__3 - d__4 * d__4;
		*allok = *allok && rts2 > para2_1.rscut2;
	}
	if (!(*allok))
	{
		*tm = ilist5_1.tlarge;
		*t1 = ilist5_1.tlarge;
		*t2 = ilist5_1.tlarge;
	}
	else
	{
		*t1 = *tm;
		*t2 = *tm;
	}
	return 0;
}

int reor_(double *t, double *tmin, int *j,
		  int *last0)
{
	static int i1, i2, i3, nc;
	extern int chin1_(int *, int *, int *,
					  int *, int *, double *, double *, int *),
		chin2_(int *, int *, int *, int *, int *,
			   double *, double *, int *),
		chin3_(int *, int *, int *, int *, int *, double *, double *,
			   int *);
	static double tmin1;
	extern int wallc_(int *, int *, int *,
					  int *, double *, double *),
		chout_(int *, int *, double *, double *, int *);
	static int icels0;
	extern int chcell_(int *, int *, int *,
					   int *, int *, double *, double *, int *),
		fixtim_(int *, double *, double *, double *,
				int *);

	icels0 = ilist1_1.icels[*j - 1];
	i1 = icels0 / 10000;
	i2 = (icels0 - i1 * 10000) / 100;
	i3 = icels0 - i1 * 10000 - i2 * 100;
	wallc_(j, &i1, &i2, &i3, t, &tmin1);
	if (*tmin <= tmin1)
	{
		nc = *last0;
	}
	else
	{
		*tmin = tmin1;
		nc = 0;
	}
	if (para5_1.iconfg == 3 || para5_1.iconfg == 5)
	{
		chcell_(j, &i1, &i2, &i3, last0, t, tmin, &nc);
	}
	else
	{
		if (i1 == 11 && i2 == 11 && i3 == 11)
		{
			chout_(j, last0, t, tmin, &nc);
		}
		else
		{
			if (para5_1.iconfg == 1)
			{
				chin1_(j, &i1, &i2, &i3, last0, t, tmin, &nc);
			}
			else if (para5_1.iconfg == 2)
			{
				chin2_(j, &i1, &i2, &i3, last0, t, tmin, &nc);
			}
			else if (para5_1.iconfg == 4)
			{
				chin3_(j, &i1, &i2, &i3, last0, t, tmin, &nc);
			}
		}
	}
	fixtim_(j, t, &tmin1, tmin, &nc);
	return 0;
}

int chout_(int *l, int *last0, double *t,
		   double *tmin, int *nc)
{
	static int i__, j, k, m1, m2, m3;
	extern int chcell_(int *, int *, int *,
					   int *, int *, double *, double *, int *);

	m1 = 11;
	m2 = 11;
	m3 = 11;
	chcell_(l, &m1, &m2, &m3, last0, t, tmin, nc);
	for (i__ = 1; i__ <= 10; ++i__)
	{
		for (j = 1; j <= 10; ++j)
		{
			for (k = 1; k <= 10; ++k)
			{
				if (i__ == 1 || i__ == 10 || j == 1 || j == 10 || k == 1 || k == 10)
				{
					chcell_(l, &i__, &j, &k, last0, t, tmin, nc);
				}
			}
		}
	}
	return 0;
}

int chin1_(int *l, int *i1, int *i2, int *i3,
		   int *last0, double *t, double *tmin, int *nc)
{

	int i__1, i__2, i__3;

	static int i__, j, k, m1, m2, m3, itest;
	extern int chcell_(int *, int *, int *,
					   int *, int *, double *, double *, int *);

	itest = 0;
	i__1 = *i1 + 1;
	for (i__ = *i1 - 1; i__ <= i__1; ++i__)
	{
		i__2 = *i2 + 1;
		for (j = *i2 - 1; j <= i__2; ++j)
		{
			i__3 = *i3 + 1;
			for (k = *i3 - 1; k <= i__3; ++k)
			{
				if (i__ >= 1 && i__ <= 10 && j >= 1 && j <= 10 && k >= 1 && k <= 10)
				{
					chcell_(l, &i__, &j, &k, last0, t, tmin, nc);
				}
				else if (itest == 0)
				{
					m1 = 11;
					m2 = 11;
					m3 = 11;
					chcell_(l, &m1, &m2, &m3, last0, t, tmin, nc);
					itest = 1;
				}
			}
		}
	}
	return 0;
}

int chin2_(int *l, int *i1, int *i2, int *i3,
		   int *last0, double *t, double *tmin, int *nc)
{

	int i__1, i__2, i__3;

	static int i__, j, k, ia, ib, ic, itest;
	extern int chcell_(int *, int *, int *,
					   int *, int *, double *, double *, int *);

	itest = 0;
	i__1 = *i1 + 1;
	for (i__ = *i1 - 1; i__ <= i__1; ++i__)
	{
		i__2 = *i2 + 1;
		for (j = *i2 - 1; j <= i__2; ++j)
		{
			i__3 = *i3 + 1;
			for (k = *i3 - 1; k <= i__3; ++k)
			{
				ia = i__;
				ib = j;
				ic = k;
				if (k >= 1 && k <= 10)
				{
					if (i__ == 0)
					{
						ia = 10;
					}
					if (i__ == 11)
					{
						ia = 1;
					}
					if (j == 0)
					{
						ib = 10;
					}
					if (j == 11)
					{
						ib = 1;
					}
					chcell_(l, &ia, &ib, &ic, last0, t, tmin, nc);
				}
			}
		}
	}
	return 0;
}

int chin3_(int *l, int *i1, int *i2, int *i3,
		   int *last0, double *t, double *tmin, int *nc)
{

	int i__1, i__2, i__3;

	static int i__, j, k, ia, ib, ic, itest;
	extern int chcell_(int *, int *, int *,
					   int *, int *, double *, double *, int *);

	itest = 0;
	i__1 = *i1 + 1;
	for (i__ = *i1 - 1; i__ <= i__1; ++i__)
	{
		i__2 = *i2 + 1;
		for (j = *i2 - 1; j <= i__2; ++j)
		{
			i__3 = *i3 + 1;
			for (k = *i3 - 1; k <= i__3; ++k)
			{
				if (i__ == 0)
				{
					ia = 10;
				}
				else if (i__ == 11)
				{
					ia = 1;
				}
				else
				{
					ia = i__;
				}
				if (j == 0)
				{
					ib = 10;
				}
				else if (j == 11)
				{
					ib = 1;
				}
				else
				{
					ib = j;
				}
				if (k == 0)
				{
					ic = 10;
				}
				else if (k == 11)
				{
					ic = 1;
				}
				else
				{
					ic = k;
				}
				chcell_(l, &ia, &ib, &ic, last0, t, tmin, nc);
			}
		}
	}
	return 0;
}

int chcell_(int *il, int *i1, int *i2, int *i3, int *last0, double *t, double *tmin, int *nc)
{

	int i__1;

	static int j, l, jj;
	extern int mintm_(int *, int *, double *,
					  int *);
	static int jmintm;

	if (para5_1.iconfg == 3 || para5_1.iconfg == 5)
	{
		jj = ilist4_1.ichkpt;
		i__1 = jj;
		for (j = 1; j <= i__1; ++j)
		{

			jmintm = j;
			if (j != *il && j != *last0)
			{
				mintm_(il, &jmintm, tmin, nc);
			}
		}
		return 0;
	}

	if (*i1 == 11 && *i2 == 11 && *i3 == 11)
	{
		l = ilist2_1.icell;
	}
	else
	{
		l = ilist2_1.icel[*i1 + (*i2 + *i3 * 10) * 10 - 111];
	}
	if (l == 0)
	{
		return 0;
	}
	j = ilist1_1.nic[l - 1];

	if (j == 0)
	{

		if (l == *il || l == *last0)
		{
			return 0;
		}
		mintm_(il, &l, tmin, nc);
	}
	else
	{
		if (l != *il && l != *last0)
		{
			mintm_(il, &l, tmin, nc);
		}
		while (j != l)
		{
			if (j != *il && j != *last0)
			{
				mintm_(il, &j, tmin, nc);
			}
			j = ilist1_1.nic[j - 1];
		}
	}
	return 0;
}

int mintm_(int *i__, int *j, double *tmin,
		   int *nc)
{
	static double t1, t2, tm;
	extern int isco_(int *, int *, logical *,
					 double *, double *, double *);
	static logical allok;

	isco_(i__, j, &allok, &tm, &t1, &t2);
	if (allok && tm < *tmin)
	{
		*tmin = tm;
		ilist5_1.ct[*i__ - 1] = t1;
		*nc = *j;
		if (para5_1.iconfg == 3 || para5_1.iconfg == 5)
		{
			aurec2_1.dgxa[*i__ - 1] = -aurec1_1.jxa * 10. * ilist3_1.size1;
			aurec2_1.dgya[*i__ - 1] = -aurec1_1.jya * 10. * ilist3_1.size2;
			if (para5_1.iconfg == 5)
			{
				aurec2_1.dgza[*i__ - 1] = -aurec1_1.jza * 10. *
										  ilist3_1.size3;
			}
		}
	}
	return 0;
}

int zpca1_(void)
{
	extern int zpca1a_(int *);

	if (ilist1_1.ictype % 2 == 0)
	{
		zpca1a_(&ilist1_1.iscat);
		zpca1a_(&ilist1_1.jscat);
	}
	return 0;
}

int zpca1a_(int *i__)
{

	int i__1;
	double d__1, d__2, d__3;

	double sqrt(double);

	static double p0, p1, p2, p3, t1, t2, et;
	static int ian, ipic;
	static double rapi;
	extern int zpca1b_(double *, double *, int *), zpca1c_(double *, double *, double *, double *,
																	   int *);

	if (para5_1.iconfg == 1)
	{
		t1 = prec3_1.fts[*i__ - 1];
		t2 = prec2_1.ft[*i__ - 1];
		ipic = 11;
	}
	else if (para5_1.iconfg == 2 || para5_1.iconfg == 3)
	{

		t1 = prec6_1.taus[*i__ - 1];
		t2 = prec5_1.tau[*i__ - 1];
		ipic = 12;
	}
	else if (para5_1.iconfg == 4 || para5_1.iconfg == 5)
	{
		t1 = prec3_1.fts[*i__ - 1];
		t2 = prec2_1.ft[*i__ - 1];
		ipic = 12;
	}
	if (para5_1.iconfg <= 3)
	{
		i__1 = ipic;
		for (ian = 1; ian <= i__1; ++ian)
		{
			if (t1 <= ana1_1.ts[ian - 1] && t2 > ana1_1.ts[ian - 1])
			{
				rapi = prec6_1.raps[*i__ - 1];

				d__1 = prec3_1.pxs[*i__ - 1];

				d__2 = prec3_1.pys[*i__ - 1];

				d__3 = para2_1.xmp;
				et = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
				zpca1b_(&rapi, &et, &ian);
			}
		}
	}
	else
	{
		i__1 = ipic;
		for (ian = 1; ian <= i__1; ++ian)
		{
			if (t1 <= ana1_1.ts[ian - 1] && t2 > ana1_1.ts[ian - 1])
			{
				p0 = prec3_1.es[*i__ - 1];
				p1 = prec3_1.pxs[*i__ - 1];
				p2 = prec3_1.pys[*i__ - 1];
				p3 = prec3_1.pzs[*i__ - 1];
				zpca1c_(&p0, &p1, &p2, &p3, &ian);
			}
		}
	}
	return 0;
}

int zpca1b_(double *rapi, double *et, int *ian)
{

	if (*rapi > para6_1.centy - .5 && *rapi < para6_1.centy + .5)
	{
		ana2_1.det2[*ian - 1] += *et;
		ana2_1.dn2[*ian - 1] += 1.;

		if (*ian == 10)
		{
		}
		if (*ian == 11)
		{
		}
		if (*ian == 12)
		{
		}

		if (*rapi > para6_1.centy - .25 && *rapi < para6_1.centy + .25)
		{
			ana2_1.det1[*ian - 1] += *et;
			ana2_1.dn1[*ian - 1] += 1.;
			if (*rapi > para6_1.centy - .1 && *rapi < para6_1.centy + .1)
			{
				ana2_1.det[*ian - 1] += *et;
				ana2_1.dn[*ian - 1] += 1.;
			}
		}
	}
	return 0;
}

int zpca1c_(double *p0, double *p1, double *p2,
			double *p3, int *ian)
{
	static int i__, j;
	static double en[4];

	en[0] = *p0;
	en[1] = *p1;
	en[2] = *p2;
	en[3] = *p3;
	for (i__ = 1; i__ <= 4; ++i__)
	{
		for (j = 1; j <= 4; ++j)
		{
			ana3_1.em[i__ + (j + (*ian << 2) << 2) - 21] += en[i__ - 1] * en[j - 1] / *p0;
		}
	}
	return 0;
}

int zpca2_(void)
{

	int s_wsle(cilist *), do_lio(int *, int *, char *, ftnlen),
		e_wsle(void);

	extern int zpca2a_(void), zpca2b_(void), zpca2c_(void);

	static cilist io___696 = {0, 25, 0, 0, 0};
	static cilist io___697 = {0, 25, 0, 0, 0};
	static cilist io___698 = {0, 25, 0, 0, 0};
	static cilist io___699 = {0, 25, 0, 0, 0};
	static cilist io___700 = {0, 25, 0, 0, 0};
	static cilist io___701 = {0, 25, 0, 0, 0};

	if (para5_1.iconfg <= 3)
	{
		zpca2a_();
	}
	else
	{
		zpca2b_();
	}
	if (para7_1.ioscar == 1)
	{
		zpca2c_();
	}

	s_wsle(&io___696);
	do_lio(&c__9, &c__1, " Event ", (ftnlen)7);
	do_lio(&c__3, &c__1, (char *)&arevt_1.iaevt, (ftnlen)sizeof(int));
	do_lio(&c__9, &c__1, ", run ", (ftnlen)6);
	do_lio(&c__3, &c__1, (char *)&arevt_1.iarun, (ftnlen)sizeof(int));
	e_wsle();
	s_wsle(&io___697);
	do_lio(&c__9, &c__1, "    number of operations = ", (ftnlen)27);
	do_lio(&c__3, &c__1, (char *)&ilist6_1.iopern, (ftnlen)sizeof(int));
	e_wsle();
	s_wsle(&io___698);
	do_lio(&c__9, &c__1, "    number of collisions between particles = ", (ftnlen)45);
	do_lio(&c__3, &c__1, (char *)&ilist6_1.icolln, (ftnlen)sizeof(int));
	e_wsle();
	s_wsle(&io___699);
	do_lio(&c__9, &c__1, "    freezeout time=", (ftnlen)19);
	do_lio(&c__5, &c__1, (char *)&ilist6_1.t, (ftnlen)sizeof(double));
	e_wsle();
	s_wsle(&io___700);
	do_lio(&c__9, &c__1, "    ending at the ", (ftnlen)18);
	do_lio(&c__3, &c__1, (char *)&rndm1_1.number, (ftnlen)sizeof(int));
	do_lio(&c__9, &c__1, "th random number", (ftnlen)16);
	e_wsle();
	s_wsle(&io___701);
	do_lio(&c__9, &c__1, "    ending collision iff=", (ftnlen)25);
	do_lio(&c__3, &c__1, (char *)&rndm2_1.iff, (ftnlen)sizeof(int));
	e_wsle();
	return 0;
}

int zpca2a_(void)
{

	int i__1, i__2;
	double d__1, d__2, d__3;

	double sqrt(double);
	int s_wsle(cilist *), do_lio(int *, int *, char *, ftnlen),
		e_wsle(void);

	static int i__, j;
	static double t1, t2, et;
	static int ian, ipic;
	static double rapi;
	extern int zpca1b_(double *, double *, int *);

	static cilist io___710 = {0, 6, 0, 0, 0};
	static cilist io___711 = {0, 6, 0, 0, 0};

	i__1 = ilist4_1.ichkpt;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		rapi = prec5_1.rap[i__ - 1];

		d__1 = prec2_1.px[i__ - 1];

		d__2 = prec2_1.py[i__ - 1];

		d__3 = para2_1.xmp;
		et = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
		for (j = 1; j <= 24; ++j)
		{
			if (rapi > j + para6_1.centy - 13. && rapi < j + para6_1.centy -
															 12.)
			{
				ana4_1.fdetdy[j - 1] += et;
				ana4_1.fdndy[j - 1] += 1.;
			}
		}
		for (j = 1; j <= 12; ++j)
		{
			if (et > (j - 1) * .5 && et < j * .5)
			{
				ana4_1.fdndpt[j - 1] += 1.;
			}
		}
		if (para5_1.iconfg == 1)
		{
			t1 = prec2_1.ft[i__ - 1];
			t2 = ilist5_1.tlarge;
			ipic = 11;
		}
		else
		{
			t1 = prec5_1.tau[i__ - 1];
			t2 = ilist5_1.tlarge;
			ipic = 12;
		}
		i__2 = ipic;
		for (ian = 1; ian <= i__2; ++ian)
		{
			if (t1 <= ana1_1.ts[ian - 1] && t2 > ana1_1.ts[ian - 1])
			{
				zpca1b_(&rapi, &et, &ian);
			}
		}
		if (para5_1.iconfg == 1)
		{
			zpca1b_(&rapi, &et, &c__12);
		}
	}
	for (ian = 1; ian <= 12; ++ian)
	{
		if (ana2_1.dn[ian - 1] == 0. || ana2_1.dn1[ian - 1] == 0. ||
			ana2_1.dn2[ian - 1] == 0.)
		{
			s_wsle(&io___710);
			do_lio(&c__9, &c__1, "event=", (ftnlen)6);
			do_lio(&c__3, &c__1, (char *)&para3_1.ievt, (ftnlen)sizeof(int));
			e_wsle();
			s_wsle(&io___711);
			do_lio(&c__9, &c__1, "dn(", (ftnlen)3);
			do_lio(&c__3, &c__1, (char *)&ian, (ftnlen)sizeof(int));
			do_lio(&c__9, &c__1, ")=", (ftnlen)2);
			do_lio(&c__5, &c__1, (char *)&ana2_1.dn[ian - 1], (ftnlen)sizeof(double));
			do_lio(&c__9, &c__1, "dn1(", (ftnlen)4);
			do_lio(&c__3, &c__1, (char *)&ian, (ftnlen)sizeof(int));
			do_lio(&c__9, &c__1, ")=", (ftnlen)2);
			do_lio(&c__5, &c__1, (char *)&ana2_1.dn1[ian - 1], (ftnlen)sizeof(double));
			do_lio(&c__9, &c__1, "dn2(", (ftnlen)4);
			do_lio(&c__3, &c__1, (char *)&ian, (ftnlen)sizeof(int));
			do_lio(&c__9, &c__1, ")=", (ftnlen)2);
			do_lio(&c__5, &c__1, (char *)&ana2_1.dn2[ian - 1], (ftnlen)sizeof(double));
			e_wsle();
		}
		ana2_1.detdy[ian - 1] += ana2_1.det[ian - 1];
		if (ana2_1.dn[ian - 1] != 0.)
		{
			ana2_1.detdn[ian - 1] += ana2_1.det[ian - 1] / ana2_1.dn[ian - 1];
		}
		ana2_1.dndy[ian - 1] += ana2_1.dn[ian - 1];
		ana2_1.detdy1[ian - 1] += ana2_1.det1[ian - 1];
		if (ana2_1.dn1[ian - 1] != 0.)
		{
			ana2_1.detdn1[ian - 1] += ana2_1.det1[ian - 1] / ana2_1.dn1[ian -
																		1];
		}
		ana2_1.dndy1[ian - 1] += ana2_1.dn1[ian - 1];
		ana2_1.detdy2[ian - 1] += ana2_1.det2[ian - 1];
		if (ana2_1.dn2[ian - 1] != 0.)
		{
			ana2_1.detdn2[ian - 1] += ana2_1.det2[ian - 1] / ana2_1.dn2[ian -
																		1];
		}
		ana2_1.dndy2[ian - 1] += ana2_1.dn2[ian - 1];
	}
	return 0;
}

int zpca2b_(void)
{

	int i__1, i__2;

	static int i__;
	static double p0, p1, p2, p3, t1, t2;
	static int ian, ipic;
	extern int zpca1c_(double *, double *,
					   double *, double *, int *);

	i__1 = ilist4_1.ichkpt;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		t1 = prec2_1.ft[i__ - 1];
		t2 = ilist5_1.tlarge;
		ipic = 12;
		i__2 = ipic;
		for (ian = 1; ian <= i__2; ++ian)
		{
			if (t1 <= ana1_1.ts[ian - 1] && t2 > ana1_1.ts[ian - 1])
			{
				p0 = prec2_1.e[i__ - 1];
				p1 = prec2_1.px[i__ - 1];
				p2 = prec2_1.py[i__ - 1];
				p3 = prec2_1.pz[i__ - 1];
				zpca1c_(&p0, &p1, &p2, &p3, &ian);
			}
		}
	}
	return 0;
}

int zpca2c_(void)
{

	static int nff = 0;

	static char fmt_101[] = "(a12)";
	static char fmt_102[] = "(2(a8,2x),\002(\002,i3,\002,\002,i6,\002)+(\002"
							",i3,\002,\002,i6,\002)\002,2x,a4,2x,e10.4,2x,i8)";
	static char fmt_103[] = "(i10,2x,i10,2x,f8.3,2x,f8.3)";
	static char fmt_104[] = "(i10,2x,i10,2x,9(e12.6,2x))";

	int i__1;

	int s_wsfe(cilist *), do_fio(int *, char *, ftnlen), e_wsfe(void);
	int s_copy(char *, char *, ftnlen, ftnlen);

	static int i__;
	static double phi;
	static char code[8];
	static double bimp, ebeam;
	static int atarg, aproj, event, ztarg;
	static char versn[8];
	static int zproj;
	static char reffra[4];
	static int ntestp;

	static cilist io___722 = {0, 26, 0, fmt_101, 0};
	static cilist io___723 = {0, 26, 0, fmt_101, 0};
	static cilist io___733 = {0, 26, 0, fmt_102, 0};
	static cilist io___737 = {0, 26, 0, fmt_103, 0};
	static cilist io___739 = {0, 26, 0, fmt_104, 0};

	if (nff == 0)
	{
		s_wsfe(&io___722);
		do_fio(&c__1, "OSCAR1997A", (ftnlen)10);
		e_wsfe();
		s_wsfe(&io___723);
		do_fio(&c__1, "final_id_p_x", (ftnlen)12);
		e_wsfe();
		s_copy(code, "ZPC", (ftnlen)8, (ftnlen)3);
		s_copy(versn, "1.0.1", (ftnlen)8, (ftnlen)5);
		aproj = -1;
		zproj = -1;
		atarg = -1;
		ztarg = -1;
		s_copy(reffra, "cm", (ftnlen)4, (ftnlen)2);
		ebeam = 0.;
		ntestp = 1;
		s_wsfe(&io___733);
		do_fio(&c__1, code, (ftnlen)8);
		do_fio(&c__1, versn, (ftnlen)8);
		do_fio(&c__1, (char *)&aproj, (ftnlen)sizeof(int));
		do_fio(&c__1, (char *)&zproj, (ftnlen)sizeof(int));
		do_fio(&c__1, (char *)&atarg, (ftnlen)sizeof(int));
		do_fio(&c__1, (char *)&ztarg, (ftnlen)sizeof(int));
		do_fio(&c__1, reffra, (ftnlen)4);
		do_fio(&c__1, (char *)&ebeam, (ftnlen)sizeof(double));
		do_fio(&c__1, (char *)&ntestp, (ftnlen)sizeof(int));
		e_wsfe();
		nff = 1;
		event = 1;
		bimp = 0.;
		phi = 0.;
	}

	s_wsfe(&io___737);
	do_fio(&c__1, (char *)&event, (ftnlen)sizeof(int));
	do_fio(&c__1, (char *)&para1_1.mul, (ftnlen)sizeof(int));
	do_fio(&c__1, (char *)&bimp, (ftnlen)sizeof(double));
	do_fio(&c__1, (char *)&phi, (ftnlen)sizeof(double));
	e_wsfe();

	i__1 = para1_1.mul;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		s_wsfe(&io___739);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(int));
		do_fio(&c__1, (char *)&prec2_1.ityp[i__ - 1], (ftnlen)sizeof(int));
		do_fio(&c__1, (char *)&prec2_1.px[i__ - 1], (ftnlen)sizeof(double));
		do_fio(&c__1, (char *)&prec2_1.py[i__ - 1], (ftnlen)sizeof(double));
		do_fio(&c__1, (char *)&prec2_1.pz[i__ - 1], (ftnlen)sizeof(double));
		do_fio(&c__1, (char *)&prec2_1.e[i__ - 1], (ftnlen)sizeof(double));
		do_fio(&c__1, (char *)&prec2_1.xmass[i__ - 1], (ftnlen)sizeof(double));
		do_fio(&c__1, (char *)&prec2_1.gx[i__ - 1], (ftnlen)sizeof(double));
		do_fio(&c__1, (char *)&prec2_1.gy[i__ - 1], (ftnlen)sizeof(double));
		do_fio(&c__1, (char *)&prec2_1.gz[i__ - 1], (ftnlen)sizeof(double));
		do_fio(&c__1, (char *)&prec2_1.ft[i__ - 1], (ftnlen)sizeof(double));
		e_wsfe();
	}
	++event;
	return 0;
}

int zpcou_(void)
{
	extern int zpcou1_(void), zpcou2_(void);

	if (para5_1.iconfg <= 3)
	{
		zpcou1_();
	}
	else
	{
		zpcou2_();
	}
	return 0;
}

int zpcou1_(void)
{
	static double dy, dy1, dy2, dpt;
	static int ntotal;

	dpt = .5;
	dy2 = 1.;
	dy1 = .5;
	dy = .2;
	ntotal = para3_1.nevnt * para3_1.nsbrun;

	return 0;
}

int zpcou2_(void)
{

	double d__1, d__2, d__3, d__4;
	olist o__1;

	int f_open(olist *), s_wsle(cilist *), do_lio(int *, int *, char *, ftnlen), e_wsle(void);

	static int i__, ian;
	static double vol;
	static int ntotal;

	static cilist io___748 = {0, 28, 0, 0, 0};
	static cilist io___750 = {0, 28, 0, 0, 0};

	o__1.oerr = 0;
	o__1.ounit = 28;
	o__1.ofnmlen = 11;
	o__1.ofnm = "ana4/em.dat";
	o__1.orl = 0;
	o__1.osta = "unknown";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	vol = ilist3_1.size1 * 1e3 * ilist3_1.size2 * ilist3_1.size3;
	ntotal = para3_1.nevnt * para3_1.nsbrun;
	for (ian = 1; ian <= 12; ++ian)
	{
		s_wsle(&io___748);
		do_lio(&c__9, &c__1, "*** for time ", (ftnlen)13);
		do_lio(&c__5, &c__1, (char *)&ana1_1.ts[ian - 1], (ftnlen)sizeof(double));
		do_lio(&c__9, &c__1, "fm(s)", (ftnlen)5);
		e_wsle();
		for (i__ = 1; i__ <= 4; ++i__)
		{
			s_wsle(&io___750);
			d__1 = ana3_1.em[i__ + ((ian << 2) + 1 << 2) - 21] / vol / ntotal;
			do_lio(&c__5, &c__1, (char *)&d__1, (ftnlen)sizeof(double));
			d__2 = ana3_1.em[i__ + ((ian << 2) + 2 << 2) - 21] / vol / ntotal;
			do_lio(&c__5, &c__1, (char *)&d__2, (ftnlen)sizeof(double));
			d__3 = ana3_1.em[i__ + ((ian << 2) + 3 << 2) - 21] / vol / ntotal;
			do_lio(&c__5, &c__1, (char *)&d__3, (ftnlen)sizeof(double));
			d__4 = ana3_1.em[i__ + ((ian << 2) + 4 << 2) - 21] / vol / ntotal;
			do_lio(&c__5, &c__1, (char *)&d__4, (ftnlen)sizeof(double));
			e_wsle();
		}
	}
	return 0;
}

int lorenz_(double *energy, double *px, double *py, double *pz, double *bex, double *bey, double *bez)
{

	double d__1, d__2, d__3;

	int s_wsle(cilist *), do_lio(int *, int *, char *, ftnlen),
		e_wsle(void);
	double sqrt(double);

	static double gam, beta2;

	static cilist io___752 = {0, 6, 0, 0, 0};

	d__1 = *bex;

	d__2 = *bey;

	d__3 = *bez;
	beta2 = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
	if (beta2 == 0.)
	{
		lor_1.enenew = *energy;
		lor_1.pxnew = *px;
		lor_1.pynew = *py;
		lor_1.pznew = *pz;
	}
	else
	{
		if (beta2 > .999999999999999)
		{
			beta2 = .999999999999999;
			s_wsle(&io___752);
			do_lio(&c__9, &c__1, "beta2=0.999999999999999", (ftnlen)23);
			e_wsle();
		}

		gam = 1. / sqrt(1. - beta2);
		lor_1.enenew = gam * (*energy - *bex * *px - *bey * *py - *bez * *pz);

		d__1 = *bex;
		lor_1.pxnew = -gam * *bex * *energy + ((gam - 1.) * (d__1 * d__1) / beta2 + 1.) * *px + (gam - 1.) * *bex * *bey / beta2 * *py + (gam - 1.) * *bex * *bez / beta2 * *pz;

		d__1 = *bey;
		lor_1.pynew = -gam * *bey * *energy + (gam - 1.) * *bex * *bey / beta2 * *px + ((gam - 1.) * (d__1 * d__1) / beta2 + 1.) * *py + (gam - 1.) * *bey * *bez / beta2 * *pz;

		d__1 = *bez;
		lor_1.pznew = -gam * *bez * *energy + (gam - 1.) * *bex * *bez / beta2 * *px + (gam - 1.) * *bey * *bez / beta2 * *py + ((gam - 1.) * (d__1 * d__1) / beta2 + 1.) * *pz;
	}
	return 0;
}

int index1_(int *n, int *m, double *arrin,
			int *indx)
{

	int i__1;

	static int i__, j, l;
	static double q;
	static int ir, indxt;

	--indx;
	--arrin;

	i__1 = *m;
	for (j = 1; j <= i__1; ++j)
	{
		indx[j] = j;
	}
	l = *m / 2 + 1;
	ir = *m;
L10:
	if (l > 1)
	{
		--l;
		indxt = indx[l];
		q = arrin[indxt];
	}
	else
	{
		indxt = indx[ir];
		q = arrin[indxt];
		indx[ir] = indx[1];
		--ir;
		if (ir == 1)
		{
			indx[1] = indxt;
			return 0;
		}
	}
	i__ = l;
	j = l + l;
L20:
	if (j <= ir)
	{
		if (j < ir)
		{
			if (arrin[indx[j]] < arrin[indx[j + 1]])
			{
				++j;
			}
		}
		if (q < arrin[indx[j]])
		{
			indx[i__] = indx[j];
			i__ = j;
			j += j;
		}
		else
		{
			j = ir + 1;
		}
		goto L20;
	}
	indx[i__] = indxt;
	goto L10;
}

double ftime1_(int *iseed)
{

	double ret_val;

	double sqrt(double);

	static double aa;
	extern double ran1_(int *);

	aa = .197327054 / par1_1.formt;

	ret_val = aa * sqrt(1. / ran1_(iseed) - 1.);
	return ret_val;
}

int cropro_(double *vx1, double *vy1, double *vz1, double *vx2, double *vy2, double *vz2)
{

	cprod_2.vx3 = *vy1 * *vz2 - *vz1 * *vy2;
	cprod_2.vy3 = *vz1 * *vx2 - *vx1 * *vz2;
	cprod_2.vz3 = *vx1 * *vy2 - *vy1 * *vx2;
	return 0;
}

int xnormv_(double *vx, double *vy, double *vz)
{

	double d__1, d__2, d__3;

	double sqrt(double);

	static double vv;

	d__1 = *vx;

	d__2 = *vy;

	d__3 = *vz;
	vv = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
	*vx /= vv;
	*vy /= vv;
	*vz /= vv;
	return 0;
}

int zprota_(double *xn1, double *xn2, double *xn3, double *theta, double *v1, double *v2, double *v3)
{

	double d__1;

	double cos(double), sin(double);

	static double c__, s, a11, a12, a13, a21, a22, a23, a31, a32, a33, vx,
		vy, vz, omc;

	vx = *v1;
	vy = *v2;
	vz = *v3;
	c__ = cos(*theta);
	omc = 1. - c__;
	s = sin(*theta);

	d__1 = *xn1;
	a11 = d__1 * d__1 * omc + c__;
	a12 = *xn1 * *xn2 * omc - s * *xn3;
	a13 = *xn1 * *xn3 * omc + s * *xn2;
	a21 = *xn1 * *xn2 * omc + s * *xn3;

	d__1 = *xn2;
	a22 = d__1 * d__1 * omc + c__;
	a23 = *xn2 * *xn3 * omc - s * *xn1;
	a31 = *xn1 * *xn3 * omc - s * *xn2;
	a32 = *xn3 * *xn2 * omc + s * *xn1;

	d__1 = *xn3;
	a33 = d__1 * d__1 * omc + c__;
	*v1 = vx * a11 + vy * a12 + vz * a13;
	*v2 = vx * a21 + vy * a22 + vz * a23;
	*v3 = vx * a31 + vy * a32 + vz * a33;
	return 0;
}

double ran1_(int *idum)
{

	static int iff = 0;

	double ret_val;

	int s_wsle(cilist *), do_lio(int *, int *, char *, ftnlen),
		e_wsle(void);

	static int j;
	static double r__[97];
	static int ix1, ix2, ix3;

	static cilist io___783 = {0, 6, 0, 0, 0};

	if (*idum < 0 || iff == 0)
	{
		iff = 1;
		ix1 = (54773 - *idum) % 259200;
		ix1 = (ix1 * 7141 + 54773) % 259200;
		ix2 = ix1 % 134456;
		ix1 = (ix1 * 7141 + 54773) % 259200;
		ix3 = ix1 % 243000;
		for (j = 1; j <= 97; ++j)
		{
			ix1 = (ix1 * 7141 + 54773) % 259200;
			ix2 = (ix2 * 8121 + 28411) % 134456;
			r__[j - 1] = ((double)ix1 + (double)ix2 *
												7.4373772832748258e-6) *
						 3.8580246913580248e-6;
		}
		*idum = 1;
	}
	ix1 = (ix1 * 7141 + 54773) % 259200;
	ix2 = (ix2 * 8121 + 28411) % 134456;
	ix3 = (ix3 * 4561 + 51349) % 243000;

	j = ix3 * 97 / 243000 + 1;

	if (j > 97 || j < 1)
	{
		s_wsle(&io___783);
		do_lio(&c__9, &c__1, "In zpc ran1, j<1 or j>97", (ftnlen)24);
		do_lio(&c__3, &c__1, (char *)&j, (ftnlen)sizeof(int));
		e_wsle();
	}
	ret_val = r__[j - 1];
	r__[j - 1] = ((double)ix1 + (double)ix2 * 7.4373772832748258e-6) * 3.8580246913580248e-6;

	++rndm1_1.number;

	return ret_val;
}
