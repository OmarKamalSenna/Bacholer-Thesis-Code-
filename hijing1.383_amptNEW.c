#include "f2c.h"
struct hparnt_1_
{
	float hipr1[100];
	int ihpr2[50];
	float hint1[100];
	int ihnt2[50];
};

#define hparnt_1 (*(struct hparnt_1_ *)&hparnt_)

struct hjcrdn_1_
{
	float yp[900], yt[900];
};

#define hjcrdn_1 (*(struct hjcrdn_1_ *)&hjcrdn_)

struct
{
	int nelt, ninthj, nelp, ninp;
} hjglbr_;

#define hjglbr_1 hjglbr_

struct hmain1_1_
{
	int natt;
	float eatt;
	int jatt, nt, np, n0, n01, n10, n11;
};

#define hmain1_1 (*(struct hmain1_1_ *)&hmain1_)

struct hmain2_1_
{
	int katt[600004];
	float patt[600004];
};

#define hmain2_1 (*(struct hmain2_1_ *)&hmain2_)

struct hstrng_1_
{
	int nfp[4500];
	float pp[4500];
	int nft[4500];
	float pt[4500];
};

#define hstrng_1 (*(struct hstrng_1_ *)&hstrng_)

struct hjjet1_1_
{
	int npj[300], kfpj[150000];
	float pjpx[150000], pjpy[150000], pjpz[150000], pjpe[150000], pjpm[150000];
	int ntj[300], kftj[150000];
	float pjtx[150000], pjty[150000], pjtz[150000], pjte[150000], pjtm[150000];
};

#define hjjet1_1 (*(struct hjjet1_1_ *)&hjjet1_)

struct hjjet2_1_
{
	int nsg, njsg[150001], iasg[450003], k1sg[15000100], k2sg[15000100];
	float pxsg[15000100], pysg[15000100], pzsg[15000100], pesg[15000100], pmsg[15000100];
};

#define hjjet2_1 (*(struct hjjet2_1_ *)&hjjet2_)

struct
{
	int ndr, iadr[1800], kfdr[900];
	float pdr[4500];
} hjjet4_;

#define hjjet4_1 hjjet4_

struct
{
	float rtdr[300002];
} xydr_;

#define xydr_1 xydr_

struct
{
	int nseed;
} rndf77_;

#define rndf77_1 rndf77_

struct
{
	int n, k[45000];
	float p[45000], v[45000];
} lujets_;

#define lujets_1 lujets_

struct
{
	int mstu[200];
	float paru[200];
	int mstj[200];
	float parj[200];
} ludat1_;

#define ludat1_1 ludat1_

struct
{
	int itypar[150001];
	float gxar[150001], gyar[150001], gzar[150001], ftar[150001], pxar[150001],
		pyar[150001], pzar[150001], pear[150001], xmar[150001];
} arprc_;

#define arprc_1 arprc_

struct
{
	int mul;
} para1_;

#define para1_1 para1_

struct
{
	double gx0[400001], gy0[400001], gz0[400001], ft0[400001], px0[400001], py0[400001], pz0[400001], e0[400001], xmass0[400001];
	int ityp0[400001];
} prec1_;

#define prec1_1 prec1_

struct
{
	double gx5[400001], gy5[400001], gz5[400001], ft5[400001], px5[400001], py5[400001], pz5[400001], e5[400001], xmass5[400001];
	int ityp5[400001];
} prec2_;

#define prec2_1 prec2_

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

struct
{
	int nsp, nst, nsi;
} srec1_;

#define srec1_1 srec1_

struct
{
	double ataui[150001], zt1[150001], zt2[150001], zt3[150001];
} srec2_;

#define srec2_1 srec2_

struct
{
	double xnprod[30], etprod[30], xnfrz[30], etfrz[30], dnprod[30],
		detpro[30], dnfrz[30], detfrz[30];
} frzout_;

#define frzout_1 frzout_

struct
{
	int nevent, isoft, isflag, izpc;
} anim_;

#define anim_1 anim_

struct
{
	double pxsgs[450003], pysgs[450003], pzsgs[450003], pesgs[450003], pmsgs[450003], gxsgs[450003], gysgs[450003], gzsgs[450003], ftsgs[450003];
	int k1sgs[450003], k2sgs[450003], njsgs[150001];
} soft_;

#define soft_1 soft_

struct
{
	int nnozpc, itypn[4001];
	float gxn[4001], gyn[4001], gzn[4001], ftn[4001], pxn[4001], pyn[4001],
		pzn[4001], een[4001], xmn[4001];
} noprec_;

#define noprec_1 noprec_

struct
{
	int itimeh;
	float bimp;
} lastt_;

#define lastt_1 lastt_

struct hijdat_1_
{
	float hidat0[100], hidat[10];
};

#define hijdat_1 (*(struct hijdat_1_ *)&hijdat_)

struct
{
	int msel, msub[200], kfin[162];
	float ckin[200];
} pysubs_;

#define pysubs_1 pysubs_

struct
{
	int mstp[200];
	float parp[200];
	int msti[200];
	float pari[200];
} pypars_;

#define pypars_1 pypars_

struct
{
	int mint[400];
	float vint[400];
} pyint1_;

#define pyint1_1 pyint1_

struct
{
	int iset[200], kfpr[400];
	float coef[4000];
	int icol[320];
} pyint2_;

#define pyint2_1 pyint2_

struct
{
	int ngen[603];
	float xsec[603];
} pyint5_;

#define pyint5_1 pyint5_

struct hpint_1_
{
	int mint4, mint5;
	float atco[4000], atxs[201];
};

#define hpint_1 (*(struct hpint_1_ *)&hpint_)

struct
{
	int mdcy[1500], mdme[4000];
	float brat[2000];
	int kfdp[10000];
} ludat3_;

#define ludat3_1 ludat3_

struct
{
	int jjp, jjt;
	float amp, amt, apx0, atx0, ampn, amtn, amp0, amt0;
	int nfdp, nfdt;
	float wp, wm, sw, xremp, xremt, dpkc1, dpkc2, pp11, pp12, pt11, pt12, ptp2,
		ptt2;
} dpmcm1_;

#define dpmcm1_1 dpmcm1_

struct
{
	int ndpm, kdpm[40];
	float pdpm1[100], pdpm2[100];
} dpmcm2_;

#define dpmcm2_1 dpmcm2_

struct
{
	float r__, d__, fnorm, w;
} wood_;

#define wood_1 wood_

struct
{
	float rr[2010], xx[2010];
} hijhb_;

#define hijhb_1 hijhb_

struct
{
	int n, ipcrs;
} njet_;

#define njet_1 njet_

struct
{
	float x4;
} besel_;

#define besel_1 besel_

struct bveg1_1_
{
	double xl[10], xu[10], acc;
	int ndim, ncall, itmx, nprn;
};

#define bveg1_1 (*(struct bveg1_1_ *)&bveg1_)

struct
{
	double xi[500], si, si2, swgt, schi;
	int ndo, it;
} bveg2_;

#define bveg2_1 bveg2_

struct
{
	double f, ti, tsi;
} bveg3_;

#define bveg3_1 bveg3_

struct sedvax_1_
{
	int num1;
};

#define sedvax_1 (*(struct sedvax_1_ *)&sedvax_)

struct
{
	float e_1[100];
	int e_2[50];
	float e_3[100];
	int e_4[50];
} hparnt_ = {1.5f, .35f, .5f, .9f, 2.f, .1f, 1.5f, 2.f, -1.f, -2.25f,
			 2.f, .5f, 1.f, 2.f, .2f, 2.f, 2.5f, .3f, .1f, 1.4f, 1.6f, 1.f,
			 2.f, 0.f, 0.f, 0.f, 0.f, 0.f, .4f, 57.f, 28.5f, 3.9f, 0.f, 0.f,
			 0.f, 0.f, 0.f, 0.f, 0.f, 3.14159f, 0.f, .4f, .1f, 1.5f, .1f, .25f,
			 0.f, .5f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
			 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
			 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
			 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
			 0.f, 0.f, 1, 3, 0, 1, 1, 1, 1, 10, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0,
			 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
			 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
			 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
			 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
			 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
			 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
			 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
			 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
			 0.f, 0.f, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			 0, 0, 0, 0, 0, 0, 0, 0, 0};

struct
{
	double e_1[21];
	int fill_2[1];
	int e_3[3];
} bveg1_ = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., .01, {0}, 1000, 100, 0};

struct
{
	int e_1;
} sedvax_ = {30123984};

struct
{
	float e_1[1800];
} hjcrdn_ = {0.f};

struct
{
	int e_1;
	float e_2;
	int e_3[7];
} hmain1_ = {0, 0.f, 0, 0, 0, 0, 0, 0, 0};

struct
{
	int e_1[600004];
	float e_2[600004];
} hmain2_ = {0};

struct
{
	int e_1[4500];
	float e_2[4500];
	int e_3[4500];
	float e_4[4500];
} hstrng_ = {0};

struct
{
	float e_1[110];
} hijdat_ = {0.f, 2.f, 1.f, .35f, 23.8f, 0.f, 0.f, 0.f, 0.f, 5.f, 0.f,
			 3.f, .8f, .35f, 24.f, 0.f, 0.f, 0.f, 0.f, 20.f, 0.f, 5.f, .8f,
			 .3f, 26.f, 0.f, 0.f, 0.f, 0.f, 53.f, 0.f, 6.f, .7f, .3f, 26.2f,
			 0.f, 0.f, 0.f, 0.f, 62.f, 0.f, 7.f, .45f, .3f, 27.f, 0.f, 0.f,
			 0.f, 0.f, 100.f, 0.f, 8.f, .215f, .3f, 28.5f, 0.f, 0.f, 0.f, 0.f,
			 200.f, 2.25f, 8.f, .21f, .5f, 28.5f, 0.f, 0.f, 0.f, 0.f, 546.f,
			 2.5f, 10.f, .19f, .6f, 28.5f, 0.f, 0.f, 0.f, 0.f, 900.f, 4.f,
			 10.f, .19f, .6f, 28.5f, 0.f, 0.f, 0.f, 0.f, 1800.f, 4.1f, 10.f,
			 .19f, .6f, 28.5f, 0.f, 0.f, 0.f, 0.f, 4e3f, 0.f, 0.f, 0.f, 0.f,
			 0.f, 0.f, 0.f, 0.f, 0.f, 0.f};

struct
{
	int e_1[150300];
	float e_2[750000];
	int e_3[150300];
	float e_4[750000];
} hjjet1_ = {0};

struct
{
	int e_1[30600205];
	float e_2[75000500];
} hjjet2_ = {0};
struct
{
	int e_1[2];
	float e_2[4201];
} hpint_ = {0};

static int c__1 = 1;
static int c__2 = 2;
static int c__9 = 9;
static double c_b21 = .6666667;
static double c_b23 = .3333333;
static float c_b26 = 0.f;
static int c__4 = 4;
static int c__3 = 3;
static int c__0 = 0;
static int c__2112 = 2112;
static int c__2212 = 2212;
static float c_b296 = 36.f;
static int c__7 = 7;
static float c_b299 = 6.f;
static float c_b302 = 1.f;
static int c__5 = 5;
static int c__6 = 6;
static double c_b876 = .33333333333333331;
static double c_b877 = -.33333333333333331;
static float c_b879 = .001f;
static float c_b894 = 20.f;
static float c_b895 = .01f;
static double c_b923 = 2.5;
static double c_b940 = .16666666;
static double c_b941 = .33333333;
static double c_b944 = .33333;

int hijing_(char *frame, float *bmin0, float *bmax0, int frame_len)
{
	static char fmt_395[] = "(2i8,f10.4,4i5)";
	static char fmt_210[] = "(i6,2(1x,f8.3),1x,f10.3,1x,f6.3,4(1x,f8.2))";
	static char fmt_211[] = "(i6,2(1x,f8.3),1x,f10.3,1x,f6.3,4(1x,e8.2))";
	static char fmt_312[] = "(1x,3f10.3,i6,2f10.3,1x,i6,1x,i3)";
	static char fmt_512[] = "(1x,3f10.3,i6,f6.3,f10.3,1x,i6,1x,i3,1x,f10.3)";

	int i__1, i__2, i__3;
	float r__1, r__2, r__3, r__4;
	double d__1, d__2, d__3, d__4, d__5, d__6;

	double sqrt(double), cos(double), sin(double), log(double);
	int s_wsle(cilist *), do_lio(int *, int *, char *, int),e_wsle(void);
	int s_stop(char *, int);
	double pow_dd(double *, double *), exp(double);
	int s_wsfe(cilist *), do_fio(int *, char *, int), e_wsfe(void),s_cmp(char *, char *, int, int);
	static int i__, j;
	static float r__, x, b2, r2, y1, y2, y3, bb;
	static int i05, ii, jp, kp;
	static float cx, gs;
	static int kt, jt;
	static float tt, xr, sx;
	static int kp2, kt2;
	static float xr1, bbx, bby, phi;
	static int isg, jtp[3];
	static float tts;
	static int ntp, lsg;
	static float rnd1, rnd2, rnd3, rrb1, rrb2;
	static int jflg;
	static float dnbp, bmin, bmax, dnbt, scip[90000],sjip[90000], rnip[90000], dnbp1, dnbp2, dnbp3, dnbt1, dnbt2, dnbt3;
	static int miss, nlop;
	static float aphx1, aphx2;
	extern double romg_(float *);
	static int jout, njet, ijet, istr, npar, jjtp, nmom, nftp;
	extern int htop_(void), ptoh_(void);
	static int nhard;
	extern double hirnd_(int *);
	static int ipcol[90000];
	static float dengy;
	static int itcol[90000], nmini, ncolt;
	extern int zpcmn_(void);
	static int itest;
	static float ttrig, gstot;
	static int nstrg;
	extern int hjana1_(void);
	static int npart;
	extern int hjana2_(void);
	static int iityp, nsbst, idstr;
	static float gstot0;
	extern int hijcsc_(int *, int *);
	static int jphard;
	extern int hijini_(void), hijhrd_(int *, int *,int *, int *, int *);
	static int jthard;
	extern int hijfrg_(int *, int *, int *),jetini_(int *, int *, int *), quench_(int *, int *);
	extern double ranart_(int *);
	static int jpmini;
	extern int hijsft_(int *, int *, int *,int *);
	static int jtmini;
	extern int luedit_(int *), hboost_(void);
	static int ierror;
	static float rantot;
	static int nsbstr;
	extern int lulist_(int *);
	static cilist io___29 = {0, 6, 0, 0, 0};
	static cilist io___53 = {0, 6, 0, 0, 0};
	static cilist io___70 = {0, 6, 0, 0, 0};
	static cilist io___72 = {0, 6, 0, 0, 0};
	static cilist io___76 = {0, 14, 0, fmt_395, 0};
	static cilist io___77 = {0, 14, 0, fmt_210, 0};
	static cilist io___78 = {0, 14, 0, fmt_211, 0};
	static cilist io___85 = {0, 14, 0, fmt_395, 0};
	static cilist io___86 = {0, 14, 0, fmt_312, 0};
	static cilist io___93 = {0, 14, 0, fmt_395, 0};
	static cilist io___94 = {0, 14, 0, fmt_512, 0};
	static cilist io___95 = {0, 6, 0, 0, 0};
	static cilist io___96 = {0, 6, 0, 0, 0};
	static cilist io___98 = {0, 6, 0, 0, 0};
	static cilist io___99 = {0, 6, 0, 0, 0};
	static cilist io___101 = {0, 6, 0, 0, 0};
	static cilist io___102 = {0, 6, 0, 0, 0};
	static cilist io___103 = {0, 6, 0, 0, 0};
	static cilist io___104 = {0, 6, 0, 0, 0};

	r__1 = *bmax0, r__2 = hparnt_1.hipr1[33] + hparnt_1.hipr1[34];
	bmax = min(r__1, r__2);
	bmin = min(*bmin0, bmax);
	if (hparnt_1.ihnt2[0] <= 1 && hparnt_1.ihnt2[2] <= 1)
	{
		bmin = 0.f;
		bmax = sqrt(hparnt_1.hipr1[30] * .1f / hparnt_1.hipr1[39]) * 2.5f;
	}
	hjcrdn_1.yp[0] = 0.f;
	hjcrdn_1.yp[1] = 0.f;
	hjcrdn_1.yp[2] = 0.f;
	if (hparnt_1.ihnt2[0] <= 1)
	{
		goto L14;
	}
	i__1 = hparnt_1.ihnt2[0];
	for (kp = 1; kp <= i__1; ++kp)
	{
	L5:
		r__ = hirnd_(&c__1);
		x = ranart_(&rndf77_1.nseed);
		cx = x * 2.f - 1.f;
		sx = sqrt(1.f - cx * cx);
		phi = ranart_(&rndf77_1.nseed) * 2.f * hparnt_1.hipr1[39];
		hjcrdn_1.yp[kp * 3 - 3] = r__ * sx * cos(phi);
		hjcrdn_1.yp[kp * 3 - 2] = r__ * sx * sin(phi);
		hjcrdn_1.yp[kp * 3 - 1] = r__ * cx;
		i__2 = kp - 1;
		for (kp2 = 1; kp2 <= i__2; ++kp2)
		{
			r__1 = hjcrdn_1.yp[kp * 3 - 3] - hjcrdn_1.yp[kp2 * 3 - 3];
			dnbp1 = r__1 * r__1;
			r__1 = hjcrdn_1.yp[kp * 3 - 2] - hjcrdn_1.yp[kp2 * 3 - 2];
			dnbp2 = r__1 * r__1;
			r__1 = hjcrdn_1.yp[kp * 3 - 1] - hjcrdn_1.yp[kp2 * 3 - 1];
			dnbp3 = r__1 * r__1;
			dnbp = dnbp1 + dnbp2 + dnbp3;
			if (dnbp < hparnt_1.hipr1[28] * hparnt_1.hipr1[28])
			{
				goto L5;
			}
		}
	}
	if (hparnt_1.ihnt2[0] == 2)
	{
		r__1 = ranart_(&rndf77_1.nseed);
		rnd1 = max(r__1, 1e-20f);
		r__1 = ranart_(&rndf77_1.nseed);
		rnd2 = max(r__1, 1e-20f);
		r__1 = ranart_(&rndf77_1.nseed);
		rnd3 = max(r__1, 1e-20f);
		r__ = -(log(rnd1) * 4.38f / 2.f + log(rnd2) * .85f / 2.f + log(rnd3) * 3.7229999999999999f / 5.2299999999999995f);
		x = ranart_(&rndf77_1.nseed);
		cx = x * 2.f - 1.f;
		sx = sqrt(1.f - cx * cx);
		phi = ranart_(&rndf77_1.nseed) * 2.f * hparnt_1.hipr1[39];
		r__ /= 2.f;
		hjcrdn_1.yp[0] = r__ * sx * cos(phi);
		hjcrdn_1.yp[1] = r__ * sx * sin(phi);
		hjcrdn_1.yp[2] = r__ * cx;
		hjcrdn_1.yp[3] = -hjcrdn_1.yp[0];
		hjcrdn_1.yp[4] = -hjcrdn_1.yp[1];
		hjcrdn_1.yp[5] = -hjcrdn_1.yp[2];
	}
	i__1 = hparnt_1.ihnt2[0] - 1;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		i__2 = hparnt_1.ihnt2[0];
		for (j = i__ + 1; j <= i__2; ++j)
		{
			y1 = hjcrdn_1.yp[i__ * 3 - 3];
			y2 = hjcrdn_1.yp[i__ * 3 - 2];
			y3 = hjcrdn_1.yp[i__ * 3 - 1];
			hjcrdn_1.yp[i__ * 3 - 3] = hjcrdn_1.yp[j * 3 - 3];
			hjcrdn_1.yp[i__ * 3 - 2] = hjcrdn_1.yp[j * 3 - 2];
			hjcrdn_1.yp[i__ * 3 - 1] = hjcrdn_1.yp[j * 3 - 1];
			hjcrdn_1.yp[j * 3 - 3] = y1;
			hjcrdn_1.yp[j * 3 - 2] = y2;
			hjcrdn_1.yp[j * 3 - 1] = y3;
		}
	}

L14:
	hjcrdn_1.yt[0] = 0.f;
	hjcrdn_1.yt[1] = 0.f;
	hjcrdn_1.yt[2] = 0.f;
	if (hparnt_1.ihnt2[2] <= 1)
	{
		miss = -1;
	}
	i__2 = hparnt_1.ihnt2[2];
	for (kt = 1; kt <= i__2; ++kt)
	{
	L15:
		r__ = hirnd_(&c__2);
		x = ranart_(&rndf77_1.nseed);
		cx = x * 2.f - 1.f;
		sx = sqrt(1.f - cx * cx);
		phi = ranart_(&rndf77_1.nseed) * 2.f * hparnt_1.hipr1[39];
		hjcrdn_1.yt[kt * 3 - 3] = r__ * sx * cos(phi);
		hjcrdn_1.yt[kt * 3 - 2] = r__ * sx * sin(phi);
		hjcrdn_1.yt[kt * 3 - 1] = r__ * cx;
		i__1 = kt - 1;
		for (kt2 = 1; kt2 <= i__1; ++kt2)
		{
			r__1 = hjcrdn_1.yt[kt * 3 - 3] - hjcrdn_1.yt[kt2 * 3 - 3];
			dnbt1 = r__1 * r__1;
			r__1 = hjcrdn_1.yt[kt * 3 - 2] - hjcrdn_1.yt[kt2 * 3 - 2];
			dnbt2 = r__1 * r__1;
			r__1 = hjcrdn_1.yt[kt * 3 - 1] - hjcrdn_1.yt[kt2 * 3 - 1];
			dnbt3 = r__1 * r__1;
			dnbt = dnbt1 + dnbt2 + dnbt3;
			if (dnbt < hparnt_1.hipr1[28] * hparnt_1.hipr1[28])
			{
				goto L15;
			}
		}
	}

	if (hparnt_1.ihnt2[2] == 2)
	{
		r__1 = ranart_(&rndf77_1.nseed);
		rnd1 = max(r__1, 1e-20f);
		r__1 = ranart_(&rndf77_1.nseed);
		rnd2 = max(r__1, 1e-20f);
		r__1 = ranart_(&rndf77_1.nseed);
		rnd3 = max(r__1, 1e-20f);
		r__ = -(log(rnd1) * 4.38f / 2.f + log(rnd2) * .85f / 2.f + log(rnd3) * 3.7229999999999999f / 5.2299999999999995f);
		x = ranart_(&rndf77_1.nseed);
		cx = x * 2.f - 1.f;
		sx = sqrt(1.f - cx * cx);
		phi = ranart_(&rndf77_1.nseed) * 2.f * hparnt_1.hipr1[39];
		r__ /= 2.f;
		hjcrdn_1.yt[0] = r__ * sx * cos(phi);
		hjcrdn_1.yt[1] = r__ * sx * sin(phi);
		hjcrdn_1.yt[2] = r__ * cx;
		hjcrdn_1.yt[3] = -hjcrdn_1.yt[0];
		hjcrdn_1.yt[4] = -hjcrdn_1.yt[1];
		hjcrdn_1.yt[5] = -hjcrdn_1.yt[2];
	}

	i__2 = hparnt_1.ihnt2[2] - 1;
	for (i__ = 1; i__ <= i__2; ++i__)
	{
		i__1 = hparnt_1.ihnt2[2];
		for (j = i__ + 1; j <= i__1; ++j)
		{
			y1 = hjcrdn_1.yt[i__ * 3 - 3];
			y2 = hjcrdn_1.yt[i__ * 3 - 2];
			y3 = hjcrdn_1.yt[i__ * 3 - 1];
			hjcrdn_1.yt[i__ * 3 - 3] = hjcrdn_1.yt[j * 3 - 3];
			hjcrdn_1.yt[i__ * 3 - 2] = hjcrdn_1.yt[j * 3 - 2];
			hjcrdn_1.yt[i__ * 3 - 1] = hjcrdn_1.yt[j * 3 - 1];
			hjcrdn_1.yt[j * 3 - 3] = y1;
			hjcrdn_1.yt[j * 3 - 2] = y2;
			hjcrdn_1.yt[j * 3 - 1] = y3;
		}
	}	
L50:
	++miss;
	if (miss > 50)
	{
		s_wsle(&io___29);
		do_lio(&c__9, &c__1, "infinite loop happened in  HIJING", (int)33);
		e_wsle();
		s_stop("", (int)0);
	}
	itest = 0;
	hmain1_1.natt = 0;
	hmain1_1.jatt = 0;
	hmain1_1.eatt = 0.f;
	hijini_();
	nlop = 0;
L60:
	hmain1_1.nt = 0;
	hmain1_1.np = 0;
	hmain1_1.n0 = 0;
	hmain1_1.n01 = 0;
	hmain1_1.n10 = 0;
	hmain1_1.n11 = 0;
	hjglbr_1.nelt = 0;
	hjglbr_1.ninthj = 0;
	hjglbr_1.nelp = 0;
	hjglbr_1.ninp = 0;
	hjjet2_1.nsg = 0;
	ncolt = 0;
	r__1 = bmin;
	r__2 = bmax;
	r__3 = bmin;
	bb = sqrt(r__1 * r__1 + ranart_(&rndf77_1.nseed) * (r__2 * r__2 - r__3 *r__3));
	phi = 0.f;
	bbx = bb * cos(phi);
	bby = bb * sin(phi);
	hparnt_1.hint1[18] = bb;
	hparnt_1.hint1[19] = phi;

	i__1 = hparnt_1.ihnt2[0];
	for (jp = 1; jp <= i__1; ++jp)
	{
		i__2 = hparnt_1.ihnt2[2];
		for (jt = 1; jt <= i__2; ++jt)
		{
			scip[jp + jt * 300 - 301] = -1.f;
			r__1 = hjcrdn_1.yp[jp * 3 - 3] + bbx - hjcrdn_1.yt[jt * 3 - 3];
			r__2 = hjcrdn_1.yp[jp * 3 - 2] + bby - hjcrdn_1.yt[jt * 3 - 2];
			b2 = r__1 * r__1 + r__2 * r__2;
			r2 = b2 * hparnt_1.hipr1[39] / hparnt_1.hipr1[30] / .1f;
			r__2 = hjcrdn_1.yp[jp * 3 - 3];
			r__3 = hjcrdn_1.yp[jp * 3 - 2];
			d__1 = (double)((float)hparnt_1.ihnt2[0]);
			r__1 = (r__2 * r__2 + r__3 * r__3) / 1.4399999999999999f / pow_dd(&d__1, &c_b21);
			rrb1 = min(r__1, 1.f);
			r__2 = hjcrdn_1.yt[jt * 3 - 3];
			r__3 = hjcrdn_1.yt[jt * 3 - 2];
			d__1 = (double)((float)hparnt_1.ihnt2[2]);
			r__1 = (r__2 * r__2 + r__3 * r__3) / 1.4399999999999999f / pow_dd(&d__1, &c_b21);
			rrb2 = min(r__1, 1.f);
			d__1 = (double)hparnt_1.ihnt2[0];
			aphx1 = hparnt_1.hipr1[5] * 4.f / 3.f * (pow_dd(&d__1, &c_b23) - 1.f) * sqrt(1.f - rrb1);
			d__1 = (double)hparnt_1.ihnt2[2];
			aphx2 = hparnt_1.hipr1[5] * 4.f / 3.f * (pow_dd(&d__1, &c_b23) - 1.f) * sqrt(1.f - rrb2);
			hparnt_1.hint1[17] = hparnt_1.hint1[13] - aphx1 * hparnt_1.hint1[14] - aphx2 * hparnt_1.hint1[15] + aphx1 * aphx2 * hparnt_1.hint1[16];
			if (hparnt_1.ihpr2[13] == 0 || hparnt_1.ihnt2[0] == 1 &&hparnt_1.ihnt2[2] == 1)
			{
				gs = 1.f - exp(-(hparnt_1.hipr1[29] + hparnt_1.hint1[17]) *romg_(&r2) / hparnt_1.hipr1[30]);
				rantot = ranart_(&rndf77_1.nseed);
				scip[jp + jt * 300 - 301] = r2;
				rnip[jp + jt * 300 - 301] = rantot;
				sjip[jp + jt * 300 - 301] = hparnt_1.hint1[17];
				++ncolt;
				ipcol[ncolt - 1] = jp;
				itcol[ncolt - 1] = jt;
			}
			gstot0 = (1.f - exp(-(hparnt_1.hipr1[29] + hparnt_1.hint1[17]) /hparnt_1.hipr1[30] / 2.f * romg_(&c_b26))) *2.f;
			r2 /= gstot0;
			gs = 1.f - exp(-(hparnt_1.hipr1[29] + hparnt_1.hint1[17]) /hparnt_1.hipr1[30] * romg_(&r2));
			gstot = (1.f - sqrt(1.f - gs)) * 2.f;
			rantot = ranart_(&rndf77_1.nseed) * gstot0;
			if (rantot > gs)
			{
				hijcsc_(&jp, &jt);
			}		
		}
	}
	lastt_1.bimp = bb;
	s_wsle(&io___53);
	do_lio(&c__9, &c__1, "#impact parameter,nlop,ncolt=", (int)29);
	do_lio(&c__4, &c__1, (char *)&lastt_1.bimp, (int)sizeof(float));
	do_lio(&c__3, &c__1, (char *)&nlop, (int)sizeof(int));
	do_lio(&c__3, &c__1, (char *)&ncolt, (int)sizeof(int));
	e_wsle();
	if (ncolt == 0)
	{
		++nlop;
		if (nlop <= 20 || hparnt_1.ihnt2[0] == 1 && hparnt_1.ihnt2[2] == 1)
		{
			goto L60;
		}
		return 0;
	}
	if (hparnt_1.ihpr2[2] != 0)
	{
		nhard = (int)(ranart_(&rndf77_1.nseed) * (ncolt - 1) + .5f) + 1;
		nhard = min(nhard, ncolt);
		jphard = ipcol[nhard - 1];
		jthard = itcol[nhard - 1];
	}

	if (hparnt_1.ihpr2[8] == 1)
	{
		nmini = (int)(ranart_(&rndf77_1.nseed) * (ncolt - 1) + .5f) + 1;
		nmini = min(nmini, ncolt);
		jpmini = ipcol[nmini - 1];
		jtmini = itcol[nmini - 1];
	}
	i__2 = hparnt_1.ihnt2[0];
	for (jp = 1; jp <= i__2; ++jp)
	{
		i__1 = hparnt_1.ihnt2[2];
		for (jt = 1; jt <= i__1; ++jt)
		{
			++hstrng_1.nfp[jp + 2999];
			++hstrng_1.nft[jt + 2999];
			if (hstrng_1.nfp[jp + 1199] <= 1 && hstrng_1.nft[jt + 1199] > 1)
			{
				++hmain1_1.np;
				++hmain1_1.n01;
			}
			else if (hstrng_1.nfp[jp + 1199] > 1 && hstrng_1.nft[jt + 1199] <= 1)
			{
				++hmain1_1.nt;
				++hmain1_1.n10;
			}
			else if (hstrng_1.nfp[jp + 1199] <= 1 && hstrng_1.nft[jt + 1199] <= 1)
			{
				++hmain1_1.np;
				++hmain1_1.nt;
				++hmain1_1.n0;
			}
			else if (hstrng_1.nfp[jp + 1199] > 1 && hstrng_1.nft[jt + 1199] > 1)
			{
				++hmain1_1.n11;
			}
			jout = 0;
			hstrng_1.nfp[jp + 2699] = 0;
			hstrng_1.nft[jt + 2699] = 0;
			if (hparnt_1.ihpr2[7] == 0 && hparnt_1.ihpr2[2] == 0)
			{
				goto L160;
			}
			if (hstrng_1.nfp[jp + 1499] < 0 || hstrng_1.nft[jt + 1499] < 0)
			{
				goto L160;
			}
			r2 = scip[jp + jt * 300 - 301];
			hparnt_1.hint1[17] = sjip[jp + jt * 300 - 301];
			tt = romg_(&r2) * hparnt_1.hint1[17] / hparnt_1.hipr1[30];
			tts = hparnt_1.hipr1[29] * romg_(&r2) / hparnt_1.hipr1[30];
			njet = 0;
			if (hparnt_1.ihpr2[2] != 0 && jp == jphard && jt == jthard)
			{
				jetini_(&jp, &jt, &c__1);
				hijhrd_(&jp, &jt, &c__0, &jflg, &c__0);
				hparnt_1.hint1[25] = hparnt_1.hint1[46];
				hparnt_1.hint1[26] = hparnt_1.hint1[47];
				hparnt_1.hint1[27] = hparnt_1.hint1[48];
				hparnt_1.hint1[28] = hparnt_1.hint1[49];
				hparnt_1.hint1[35] = hparnt_1.hint1[66];
				hparnt_1.hint1[36] = hparnt_1.hint1[67];
				hparnt_1.hint1[37] = hparnt_1.hint1[68];
				hparnt_1.hint1[38] = hparnt_1.hint1[69];

				if (dabs(hparnt_1.hint1[45]) > hparnt_1.hipr1[10] && jflg ==2)
				{
					hstrng_1.nfp[jp + 1799] = 1;
				}
				if (dabs(hparnt_1.hint1[55]) > hparnt_1.hipr1[10] && jflg ==2)
				{
					hstrng_1.nft[jt + 1799] = 1;
				}
				r__1 = dabs(hparnt_1.hint1[45]), r__2 = dabs(hparnt_1.hint1[55]);
				if (max(r__1, r__2) > hparnt_1.hipr1[10] && jflg >= 3)
				{
					hjjet2_1.iasg[hjjet2_1.nsg + 300001] = 1;
				}
				hparnt_1.ihnt2[8] = hparnt_1.ihnt2[13];
				hparnt_1.ihnt2[9] = hparnt_1.ihnt2[14];
				for (i05 = 1; i05 <= 5; ++i05)
				{
					hparnt_1.hint1[i05 + 19] = hparnt_1.hint1[i05 + 39];
					hparnt_1.hint1[i05 + 29] = hparnt_1.hint1[i05 + 49];
				}
				jout = 1;
				if (hparnt_1.ihpr2[7] == 0)
				{
					goto L160;
				}
				r__2 = hjcrdn_1.yp[jp * 3 - 3];
				r__3 = hjcrdn_1.yp[jp * 3 - 2];
				d__1 = (double)((float)hparnt_1.ihnt2[0]);
				r__1 = (r__2 * r__2 + r__3 * r__3) / 1.4399999999999999f /
					   pow_dd(&d__1, &c_b21);
				rrb1 = min(r__1, 1.f);
				r__2 = hjcrdn_1.yt[jt * 3 - 3];
				r__3 = hjcrdn_1.yt[jt * 3 - 2];
				d__1 = (double)((float)hparnt_1.ihnt2[2]);
				r__1 = (r__2 * r__2 + r__3 * r__3) / 1.4399999999999999f /
					   pow_dd(&d__1, &c_b21);
				rrb2 = min(r__1, 1.f);
				d__1 = (double)hparnt_1.ihnt2[0];
				aphx1 = hparnt_1.hipr1[5] * 4.f / 3.f * (pow_dd(&d__1, &c_b23) - 1.f) * sqrt(1.f - rrb1);
				d__1 = (double)hparnt_1.ihnt2[2];
				aphx2 = hparnt_1.hipr1[5] * 4.f / 3.f * (pow_dd(&d__1, &c_b23) - 1.f) * sqrt(1.f - rrb2);
				hparnt_1.hint1[64] = hparnt_1.hint1[60] - aphx1 * hparnt_1.hint1[61] - aphx2 * hparnt_1.hint1[62] +aphx1 * aphx2 * hparnt_1.hint1[63];
				ttrig = romg_(&r2) * hparnt_1.hint1[64] / hparnt_1.hipr1[30];
				njet = -1;
				xr1 = -log(exp(-ttrig) + ranart_(&rndf77_1.nseed) * (1.f - exp(-ttrig)));
			L106:
				++njet;
				r__1 = ranart_(&rndf77_1.nseed);
				xr1 -= log((max(r__1, 1e-20f)));
				if (xr1 < ttrig)
				{
					goto L106;
				}
				xr = 0.f;
			L107:
				++njet;
				r__1 = ranart_(&rndf77_1.nseed);
				xr -= log((max(r__1, 1e-20f)));
				if (xr < tt - ttrig)
				{
					goto L107;
				}
				--njet;
				goto L112;
			}
			if (hparnt_1.ihpr2[8] == 1 && jp == jpmini && jt == jtmini)
			{
				goto L110;
			}
			if (hparnt_1.ihpr2[7] > 0 && rnip[jp + jt * 300 - 301] < exp(-tt) * (1.f - exp(-tts)))
			{
				goto L160;
			}
		L110:
			xr = -log(exp(-tt) + ranart_(&rndf77_1.nseed) * (1.f - exp(-tt)));
		L111:
			++njet;
			r__1 = ranart_(&rndf77_1.nseed);
			xr -= log((max(r__1, 1e-20f)));
			if (xr < tt)
			{
				goto L111;
			}
		L112:
			njet = min(njet, hparnt_1.ihpr2[7]);
			if (hparnt_1.ihpr2[7] < 0)
			{
				njet = abs(hparnt_1.ihpr2[7]);
			}
			i__3 = njet;
			for (ijet = 1; ijet <= i__3; ++ijet)
			{
				jetini_(&jp, &jt, &c__0);
				hijhrd_(&jp, &jt, &jout, &jflg, &c__1);
				if (jflg == 0)
				{
					goto L160;
				}
				if (jflg < 0)
				{
					if (hparnt_1.ihpr2[9] != 0)
					{
						s_wsle(&io___70);
						do_lio(&c__9, &c__1, "error occured in HIJHRD", (int)23);
						e_wsle();
					}
					goto L50;
				}
				++jout;
				if (dabs(hparnt_1.hint1[45]) > hparnt_1.hipr1[10] && jflg == 2)
				{
					hstrng_1.nfp[jp + 1799] = 1;
				}
				if (dabs(hparnt_1.hint1[55]) > hparnt_1.hipr1[10] && jflg ==2)
				{
					hstrng_1.nft[jt + 1799] = 1;
				}
				r__1 = dabs(hparnt_1.hint1[45]), r__2 = dabs(hparnt_1.hint1[55]);
				if (max(r__1, r__2) > hparnt_1.hipr1[10] && jflg >= 3)
				{
					hjjet2_1.iasg[hjjet2_1.nsg + 300001] = 1;
				}
			}
		L160:
			hijsft_(&jp, &jt, &jout, &ierror);
			if (ierror != 0)
			{
				if (hparnt_1.ihpr2[9] != 0)
				{
					s_wsle(&io___72);
					do_lio(&c__9, &c__1, "error occured in HIJSFT", (int)23);
					e_wsle();
				}
				goto L50;
			}

			hmain1_1.jatt += jout;
		}
	}

	i__1 = hparnt_1.ihnt2[0];
	for (jp = 1; jp <= i__1; ++jp)
	{
		if (hstrng_1.nfp[jp + 1199] > 2)
		{
			++hjglbr_1.ninp;
		}
		else if (hstrng_1.nfp[jp + 1199] == 2 || hstrng_1.nfp[jp + 1199] ==1)
		{
			++hjglbr_1.nelp;
		}
	}
	i__1 = hparnt_1.ihnt2[2];
	for (jt = 1; jt <= i__1; ++jt)
	{
		if (hstrng_1.nft[jt + 1199] > 2)
		{
			++hjglbr_1.ninthj;
		}
		else if (hstrng_1.nft[jt + 1199] == 2 || hstrng_1.nft[jt + 1199] ==
													 1)
		{
			++hjglbr_1.nelt;
		}
	}

	if ((hparnt_1.ihpr2[7] != 0 || hparnt_1.ihpr2[2] != 0) && hparnt_1.ihpr2[3] > 0 && hparnt_1.ihnt2[0] > 1 && hparnt_1.ihnt2[2] > 1)
	{
		i__1 = hparnt_1.ihnt2[0];
		for (i__ = 1; i__ <= i__1; ++i__)
		{
			if (hstrng_1.nfp[i__ + 1799] == 1)
			{
				quench_(&i__, &c__1);
			}
		}
		i__1 = hparnt_1.ihnt2[2];
		for (i__ = 1; i__ <= i__1; ++i__)
		{
			if (hstrng_1.nft[i__ + 1799] == 1)
			{
				quench_(&i__, &c__2);
			}
		}
		i__1 = hjjet2_1.nsg;
		for (isg = 1; isg <= i__1; ++isg)
		{
			if (hjjet2_1.iasg[isg + 300001] == 1)
			{
				quench_(&isg, &c__3);
			}
		}
	}
	if (anim_1.isoft == 1)
	{
		anim_1.isflag = 1;
		srec1_1.nsp = hparnt_1.ihnt2[0];
		srec1_1.nst = hparnt_1.ihnt2[2];
		srec1_1.nsi = hjjet2_1.nsg;
		istr = 0;
		npar = 0;
		i__1 = hparnt_1.ihnt2[0];
		for (i__ = 1; i__ <= i__1; ++i__)
		{
			++istr;
			i__2 = hjjet1_1.npj[i__ - 1];
			for (j = 1; j <= i__2; ++j)
			{
				if (hjjet1_1.kfpj[i__ + j * 300 - 301] == 21)
				{
					++npar;
					ilist7_1.lstrg0[npar - 1] = istr;
					ilist7_1.lpart0[npar - 1] = j;
					prec1_1.ityp0[npar - 1] = hjjet1_1.kfpj[i__ + j * 300 -301];
					prec1_1.gx0[npar - 1] = (double)(hjcrdn_1.yp[i__ * 3 - 3] + bb * .5f);
					prec1_1.gy0[npar - 1] = (double)hjcrdn_1.yp[i__ * 3 - 2];
					prec1_1.gz0[npar - 1] = 0.;
					prec1_1.ft0[npar - 1] = 0.;
					prec1_1.px0[npar - 1] = (double)hjjet1_1.pjpx[i__ +j * 300 - 301];
					prec1_1.py0[npar - 1] = (double)hjjet1_1.pjpy[i__ +j * 300 - 301];
					prec1_1.pz0[npar - 1] = (double)hjjet1_1.pjpz[i__ +j * 300 - 301];
					prec1_1.xmass0[npar - 1] = (double)hjjet1_1.pjpm[i__ + j * 300 - 301];
					d__1 = prec1_1.px0[npar - 1];
					d__2 = prec1_1.py0[npar - 1];
					d__3 = prec1_1.pz0[npar - 1];
					d__4 = prec1_1.xmass0[npar - 1];
					prec1_1.e0[npar - 1] = sqrt(d__1 * d__1 + d__2 * d__2 +d__3 * d__3 + d__4 * d__4);
				}
			}
		}
		i__1 = hparnt_1.ihnt2[2];
		for (i__ = 1; i__ <= i__1; ++i__)
		{
			++istr;
			i__2 = hjjet1_1.ntj[i__ - 1];
			for (j = 1; j <= i__2; ++j)
			{
				if (hjjet1_1.kftj[i__ + j * 300 - 301] == 21)
				{
					++npar;
					ilist7_1.lstrg0[npar - 1] = istr;
					ilist7_1.lpart0[npar - 1] = j;
					prec1_1.ityp0[npar - 1] = hjjet1_1.kftj[i__ + j * 300 -301];
					prec1_1.gx0[npar - 1] = (double)(hjcrdn_1.yt[i__ * 3 - 3] - bb * .5f);
					prec1_1.gy0[npar - 1] = (double)hjcrdn_1.yt[i__ * 3 - 2];
					prec1_1.gz0[npar - 1] = 0.;
					prec1_1.ft0[npar - 1] = 0.;
					prec1_1.px0[npar - 1] = (double)hjjet1_1.pjtx[i__ +j * 300 - 301];
					prec1_1.py0[npar - 1] = (double)hjjet1_1.pjty[i__ +j * 300 - 301];
					prec1_1.pz0[npar - 1] = (double)hjjet1_1.pjtz[i__ +j * 300 - 301];
					prec1_1.xmass0[npar - 1] = (double)hjjet1_1.pjtm[i__ + j * 300 - 301];
					d__1 = prec1_1.px0[npar - 1];
					d__2 = prec1_1.py0[npar - 1];
					d__3 = prec1_1.pz0[npar - 1];
					d__4 = prec1_1.xmass0[npar - 1];
					prec1_1.e0[npar - 1] = sqrt(d__1 * d__1 + d__2 * d__2 +d__3 * d__3 + d__4 * d__4);
				}
			}
		}
		i__1 = hjjet2_1.nsg;
		for (i__ = 1; i__ <= i__1; ++i__)
		{
			++istr;
			i__2 = hjjet2_1.njsg[i__ - 1];
			for (j = 1; j <= i__2; ++j)
			{
				if (hjjet2_1.k2sg[i__ + j * 150001 - 150002] == 21)
				{
					++npar;
					ilist7_1.lstrg0[npar - 1] = istr;
					ilist7_1.lpart0[npar - 1] = j;
					prec1_1.ityp0[npar - 1] = hjjet2_1.k2sg[i__ + j * 150001 - 150002];
					prec1_1.gx0[npar - 1] = (double)(hjcrdn_1.yp[hjjet2_1.iasg[i__ - 1] * 3 - 3] + hjcrdn_1.yt[hjjet2_1.iasg[i__ + 150000] * 3 - 3]) * .5;
					prec1_1.gy0[npar - 1] = (double)(hjcrdn_1.yp[hjjet2_1.iasg[i__ - 1] * 3 - 2] + hjcrdn_1.yt[hjjet2_1.iasg[i__ + 150000] * 3 - 2]) * .5;
					prec1_1.gz0[npar - 1] = 0.;
					prec1_1.ft0[npar - 1] = 0.;
					prec1_1.px0[npar - 1] = (double)hjjet2_1.pxsg[i__ +j * 150001 - 150002];
					prec1_1.py0[npar - 1] = (double)hjjet2_1.pysg[i__ +j * 150001 - 150002];
					prec1_1.pz0[npar - 1] = (double)hjjet2_1.pzsg[i__ +j * 150001 - 150002];
					prec1_1.xmass0[npar - 1] = (double)hjjet2_1.pmsg[i__ + j * 150001 - 150002];
					d__1 = prec1_1.px0[npar - 1];
					d__2 = prec1_1.py0[npar - 1];
					d__3 = prec1_1.pz0[npar - 1];
					d__4 = prec1_1.xmass0[npar - 1];
					prec1_1.e0[npar - 1] = sqrt(d__1 * d__1 + d__2 * d__2 +d__3 * d__3 + d__4 * d__4);
				}
			}
		}
		para1_1.mul = npar;
		hjana1_();
		zpcmn_();
		s_wsfe(&io___76);
		do_fio(&c__1, (char *)&itest, (int)sizeof(int));
		do_fio(&c__1, (char *)&para1_1.mul, (int)sizeof(int));
		do_fio(&c__1, (char *)&lastt_1.bimp, (int)sizeof(float));
		do_fio(&c__1, (char *)&hjglbr_1.nelp, (int)sizeof(int));
		do_fio(&c__1, (char *)&hjglbr_1.ninp, (int)sizeof(int));
		do_fio(&c__1, (char *)&hjglbr_1.nelt, (int)sizeof(int));
		do_fio(&c__1, (char *)&hjglbr_1.ninthj, (int)sizeof(int));
		e_wsfe();
		i__1 = para1_1.mul;
		for (i__ = 1; i__ <= i__1; ++i__)
		{
			d__5 = (d__1 = prec2_1.gx5[i__ - 1], abs(d__1)), d__6 = (d__2 = prec2_1.gy5[i__ - 1], abs(d__2)), d__5 = max(d__5, d__6),
			d__6 = (d__3 = prec2_1.gz5[i__ - 1], abs(d__3)), d__5 = max(d__5, d__6), d__6 = (d__4 = prec2_1.ft5[i__ - 1], abs(d__4));
			if (max(d__5, d__6) < 9999.)
			{
				s_wsfe(&io___77);
				do_fio(&c__1, (char *)&prec2_1.ityp5[i__ - 1], (int)sizeof(int));
				do_fio(&c__1, (char *)&prec2_1.px5[i__ - 1], (int)sizeof(double));
				do_fio(&c__1, (char *)&prec2_1.py5[i__ - 1], (int)sizeof(double));
				do_fio(&c__1, (char *)&prec2_1.pz5[i__ - 1], (int)sizeof(double));
				do_fio(&c__1, (char *)&prec2_1.xmass5[i__ - 1], (int)sizeof(double));
				do_fio(&c__1, (char *)&prec2_1.gx5[i__ - 1], (int)sizeof(double));
				do_fio(&c__1, (char *)&prec2_1.gy5[i__ - 1], (int)sizeof(double));
				do_fio(&c__1, (char *)&prec2_1.gz5[i__ - 1], (int)sizeof(double));
				do_fio(&c__1, (char *)&prec2_1.ft5[i__ - 1], (int)sizeof(double));
				e_wsfe();
			}
			else
			{
				s_wsfe(&io___78);
				do_fio(&c__1, (char *)&prec2_1.ityp5[i__ - 1], (int)sizeof(int));
				do_fio(&c__1, (char *)&prec2_1.px5[i__ - 1], (int)sizeof(double));
				do_fio(&c__1, (char *)&prec2_1.py5[i__ - 1], (int)sizeof(double));
				do_fio(&c__1, (char *)&prec2_1.pz5[i__ - 1], (int)sizeof(double));
				do_fio(&c__1, (char *)&prec2_1.xmass5[i__ - 1], (int)sizeof(double));
				do_fio(&c__1, (char *)&prec2_1.gx5[i__ - 1], (int)sizeof(double));
				do_fio(&c__1, (char *)&prec2_1.gy5[i__ - 1], (int)sizeof(double));
				do_fio(&c__1, (char *)&prec2_1.gz5[i__ - 1], (int)sizeof(double));
				do_fio(&c__1, (char *)&prec2_1.ft5[i__ - 1], (int)sizeof(double));
				e_wsfe();
			}
		}
		++itest;
		i__1 = para1_1.mul;
		for (i__ = 1; i__ <= i__1; ++i__)
		{
			if (ilist8_1.lstrg1[i__ - 1] <= srec1_1.nsp)
			{
				nstrg = ilist8_1.lstrg1[i__ - 1];
				npart = ilist8_1.lpart1[i__ - 1];
				hjjet1_1.kfpj[nstrg + npart * 300 - 301] = prec2_1.ityp5[i__ - 1];
				hjjet1_1.pjpx[nstrg + npart * 300 - 301] = (float)prec2_1.px5[i__ - 1];
				hjjet1_1.pjpy[nstrg + npart * 300 - 301] = (float)prec2_1.py5[i__ - 1];
				hjjet1_1.pjpz[nstrg + npart * 300 - 301] = (float)prec2_1.pz5[i__ - 1];
				hjjet1_1.pjpe[nstrg + npart * 300 - 301] = (float)prec2_1.e5[i__ - 1];
				hjjet1_1.pjpm[nstrg + npart * 300 - 301] = (float)
															   prec2_1.xmass5[i__ - 1];
			}
			else if (ilist8_1.lstrg1[i__ - 1] <= srec1_1.nsp + srec1_1.nst)
			{
				nstrg = ilist8_1.lstrg1[i__ - 1] - srec1_1.nsp;
				npart = ilist8_1.lpart1[i__ - 1];
				hjjet1_1.kftj[nstrg + npart * 300 - 301] = prec2_1.ityp5[i__ - 1];
				hjjet1_1.pjtx[nstrg + npart * 300 - 301] = (float)prec2_1.px5[i__ - 1];
				hjjet1_1.pjty[nstrg + npart * 300 - 301] = (float)prec2_1.py5[i__ - 1];
				hjjet1_1.pjtz[nstrg + npart * 300 - 301] = (float)prec2_1.pz5[i__ - 1];
				hjjet1_1.pjte[nstrg + npart * 300 - 301] = (float)prec2_1.e5[i__ - 1];
				hjjet1_1.pjtm[nstrg + npart * 300 - 301] = (float)
															   prec2_1.xmass5[i__ - 1];
			}
			else
			{
				nstrg = ilist8_1.lstrg1[i__ - 1] - srec1_1.nsp - srec1_1.nst;
				npart = ilist8_1.lpart1[i__ - 1];
				hjjet2_1.k2sg[nstrg + npart * 150001 - 150002] =
					prec2_1.ityp5[i__ - 1];
				hjjet2_1.pxsg[nstrg + npart * 150001 - 150002] = (float)
																	 prec2_1.px5[i__ - 1];
				hjjet2_1.pysg[nstrg + npart * 150001 - 150002] = (float)
																	 prec2_1.py5[i__ - 1];
				hjjet2_1.pzsg[nstrg + npart * 150001 - 150002] = (float)
																	 prec2_1.pz5[i__ - 1];
				hjjet2_1.pesg[nstrg + npart * 150001 - 150002] = (float)
																	 prec2_1.e5[i__ - 1];
				hjjet2_1.pmsg[nstrg + npart * 150001 - 150002] = (float)
																	 prec2_1.xmass5[i__ - 1];
			}
		}
		hjana2_();
	}
	else if (anim_1.isoft == 2)
	{
		srec1_1.nsp = hparnt_1.ihnt2[0];
		srec1_1.nst = hparnt_1.ihnt2[2];
		srec1_1.nsi = hjjet2_1.nsg;
		npar = 0;
		istr = 0;

		ludat1_1.mstj[0] = 0;
		hparnt_1.ihpr2[0] = 0;
		anim_1.isflag = 0;
		if (hparnt_1.ihpr2[19] != 0)
		{
			for (ntp = 1; ntp <= 2; ++ntp)
			{
				i__1 = hparnt_1.ihnt2[(ntp << 1) - 2];
				for (jjtp = 1; jjtp <= i__1; ++jjtp)
				{
					++istr;
					hijfrg_(&jjtp, &ntp, &ierror);
					if (ntp == 1)
					{
						i__2 = lujets_1.n - 2;
						hjjet1_1.npj[jjtp - 1] = max(i__2, 0);
					}
					else
					{
						i__2 = lujets_1.n - 2;
						hjjet1_1.ntj[jjtp - 1] = max(i__2, 0);
					}
					i__2 = lujets_1.n;
					for (ii = 1; ii <= i__2; ++ii)
					{
						++npar;
						ilist7_1.lstrg0[npar - 1] = istr;
						ilist7_1.lpart0[npar - 1] = ii;
						prec1_1.ityp0[npar - 1] = lujets_1.k[ii + 8999];
						prec1_1.gz0[npar - 1] = 0.;
						prec1_1.ft0[npar - 1] = 0.;
						prec1_1.px0[npar - 1] = (double)lujets_1.p[ii -
																   1];
						prec1_1.py0[npar - 1] = (double)lujets_1.p[ii +
																   8999];
						prec1_1.pz0[npar - 1] = (double)lujets_1.p[ii +
																   17999];
						prec1_1.xmass0[npar - 1] = (double)lujets_1.p[ii + 35999];
						d__1 = prec1_1.px0[npar - 1];
						d__2 = prec1_1.py0[npar - 1];
						d__3 = prec1_1.pz0[npar - 1];
						d__4 = prec1_1.xmass0[npar - 1];
						prec1_1.e0[npar - 1] = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3 + d__4 * d__4);
						if (ntp == 1)
						{
							prec1_1.gx0[npar - 1] = (double)(hjcrdn_1.yp[jjtp * 3 - 3] + bb * .5f);
							prec1_1.gy0[npar - 1] = (double)hjcrdn_1.yp[jjtp * 3 - 2];
							iityp = prec1_1.ityp0[npar - 1];
							nstrg = ilist7_1.lstrg0[npar - 1];
							if (iityp == 2112 || iityp == 2212)
							{
							}
							else if ((iityp == 1 || iityp == 2) && (ii == 1 || ii == lujets_1.n))
							{
								hstrng_1.pp[nstrg + 1499] = (float)
																prec1_1.px0[npar - 1];
								hstrng_1.pp[nstrg + 1799] = (float)
																prec1_1.py0[npar - 1];
								hstrng_1.pp[nstrg + 3899] = (float)
																prec1_1.xmass0[npar - 1];
							}
							else if ((iityp == 1103 || iityp == 2101 ||
									  iityp == 2103 || (float)iityp == 2203.f ||
									  iityp == 3101 || (float)iityp == 3103.f || iityp == 3201 || iityp == 3203 ||
									  iityp == 3303) &&
									 (ii == 1 || ii ==
													 lujets_1.n))
							{
								hstrng_1.pp[nstrg + 2099] = (float)
																prec1_1.px0[npar - 1];
								hstrng_1.pp[nstrg + 2399] = (float)
																prec1_1.py0[npar - 1];
								hstrng_1.pp[nstrg + 4199] = (float)
																prec1_1.xmass0[npar - 1];
							}
							else
							{
								npart = ilist7_1.lpart0[npar - 1] - 1;
								hjjet1_1.kfpj[nstrg + npart * 300 - 301] =
									prec1_1.ityp0[npar - 1];
								hjjet1_1.pjpx[nstrg + npart * 300 - 301] = (float)prec1_1.px0[npar - 1];
								hjjet1_1.pjpy[nstrg + npart * 300 - 301] = (float)prec1_1.py0[npar - 1];
								hjjet1_1.pjpz[nstrg + npart * 300 - 301] = (float)prec1_1.pz0[npar - 1];
								hjjet1_1.pjpe[nstrg + npart * 300 - 301] = (float)prec1_1.e0[npar - 1];
								hjjet1_1.pjpm[nstrg + npart * 300 - 301] = (float)prec1_1.xmass0[npar - 1];
							}
						}
						else
						{
							prec1_1.gx0[npar - 1] = (double)(hjcrdn_1.yt[jjtp * 3 - 3] - bb * .5f);
							prec1_1.gy0[npar - 1] = (double)hjcrdn_1.yt[jjtp * 3 - 2];
							iityp = prec1_1.ityp0[npar - 1];
							nstrg = ilist7_1.lstrg0[npar - 1] - srec1_1.nsp;
							if (iityp == 2112 || iityp == 2212)
							{
							}
							else if ((iityp == 1 || iityp == 2) && (ii == 1 || ii == lujets_1.n))
							{
								hstrng_1.pt[nstrg + 1499] = (float)
																prec1_1.px0[npar - 1];
								hstrng_1.pt[nstrg + 1799] = (float)
																prec1_1.py0[npar - 1];
								hstrng_1.pt[nstrg + 3899] = (float)
																prec1_1.xmass0[npar - 1];
							}
							else if ((iityp == 1103 || iityp == 2101 ||
									  iityp == 2103 || (float)iityp == 2203.f ||
									  iityp == 3101 || (float)iityp == 3103.f || iityp == 3201 || iityp == 3203 ||
									  iityp == 3303) &&
									 (ii == 1 || ii ==
													 lujets_1.n))
							{
								hstrng_1.pt[nstrg + 2099] = (float)
																prec1_1.px0[npar - 1];
								hstrng_1.pt[nstrg + 2399] = (float)
																prec1_1.py0[npar - 1];
								hstrng_1.pt[nstrg + 4199] = (float)
																prec1_1.xmass0[npar - 1];
							}
							else
							{
								npart = ilist7_1.lpart0[npar - 1] - 1;
								hjjet1_1.kftj[nstrg + npart * 300 - 301] =
									prec1_1.ityp0[npar - 1];
								hjjet1_1.pjtx[nstrg + npart * 300 - 301] = (float)prec1_1.px0[npar - 1];
								hjjet1_1.pjty[nstrg + npart * 300 - 301] = (float)prec1_1.py0[npar - 1];
								hjjet1_1.pjtz[nstrg + npart * 300 - 301] = (float)prec1_1.pz0[npar - 1];
								hjjet1_1.pjte[nstrg + npart * 300 - 301] = (float)prec1_1.e0[npar - 1];
								hjjet1_1.pjtm[nstrg + npart * 300 - 301] = (float)prec1_1.xmass0[npar - 1];
							}
						}
					}
				}
			}
			i__1 = hjjet2_1.nsg;
			for (isg = 1; isg <= i__1; ++isg)
			{
				++istr;
				hijfrg_(&isg, &c__3, &ierror);
				hjjet2_1.njsg[isg - 1] = lujets_1.n;

				i__2 = lujets_1.n;
				for (ii = 1; ii <= i__2; ++ii)
				{
					++npar;
					ilist7_1.lstrg0[npar - 1] = istr;
					ilist7_1.lpart0[npar - 1] = ii;
					prec1_1.ityp0[npar - 1] = lujets_1.k[ii + 8999];
					prec1_1.gx0[npar - 1] = (double)(hjcrdn_1.yp[hjjet2_1.iasg[isg - 1] * 3 - 3] + hjcrdn_1.yt[hjjet2_1.iasg[isg + 150000] * 3 - 3]) * .5;
					prec1_1.gy0[npar - 1] = (double)(hjcrdn_1.yp[hjjet2_1.iasg[isg - 1] * 3 - 2] + hjcrdn_1.yt[hjjet2_1.iasg[isg + 150000] * 3 - 2]) * .5;
					prec1_1.gz0[npar - 1] = 0.;
					prec1_1.ft0[npar - 1] = 0.;
					prec1_1.px0[npar - 1] = (double)lujets_1.p[ii - 1];
					prec1_1.py0[npar - 1] = (double)lujets_1.p[ii + 8999];
					prec1_1.pz0[npar - 1] = (double)lujets_1.p[ii +
															   17999];
					prec1_1.xmass0[npar - 1] = (double)lujets_1.p[ii +
																  35999];
					d__1 = prec1_1.px0[npar - 1];
					d__2 = prec1_1.py0[npar - 1];
					d__3 = prec1_1.pz0[npar - 1];
					d__4 = prec1_1.xmass0[npar - 1];
					prec1_1.e0[npar - 1] = sqrt(d__1 * d__1 + d__2 * d__2 +
												d__3 * d__3 + d__4 * d__4);
				}
			}
		}
		para1_1.mul = npar;
		hjana1_();
		zpcmn_();
		s_wsfe(&io___85);
		do_fio(&c__1, (char *)&itest, (int)sizeof(int));
		do_fio(&c__1, (char *)&para1_1.mul, (int)sizeof(int));
		do_fio(&c__1, (char *)&lastt_1.bimp, (int)sizeof(float));
		do_fio(&c__1, (char *)&hjglbr_1.nelp, (int)sizeof(int));
		do_fio(&c__1, (char *)&hjglbr_1.ninp, (int)sizeof(int));
		do_fio(&c__1, (char *)&hjglbr_1.nelt, (int)sizeof(int));
		do_fio(&c__1, (char *)&hjglbr_1.ninthj, (int)sizeof(int));
		e_wsfe();
		++itest;
		i__1 = para1_1.mul;
		for (i__ = 1; i__ <= i__1; ++i__)
		{
			s_wsfe(&io___86);
			do_fio(&c__1, (char *)&prec2_1.px5[i__ - 1], (int)sizeof(double));
			do_fio(&c__1, (char *)&prec2_1.py5[i__ - 1], (int)sizeof(double));
			do_fio(&c__1, (char *)&prec2_1.pz5[i__ - 1], (int)sizeof(double));
			do_fio(&c__1, (char *)&prec2_1.ityp5[i__ - 1], (int)sizeof(int));
			do_fio(&c__1, (char *)&prec2_1.xmass5[i__ - 1], (int)sizeof(double));
			do_fio(&c__1, (char *)&prec2_1.e5[i__ - 1], (int)sizeof(double));
			do_fio(&c__1, (char *)&ilist8_1.lstrg1[i__ - 1], (int)sizeof(int));
			do_fio(&c__1, (char *)&ilist8_1.lpart1[i__ - 1], (int)sizeof(int));
			e_wsfe();
		}
		for (nmom = 1; nmom <= 5; ++nmom)
		{
			i__1 = srec1_1.nsp;
			for (nstrg = 1; nstrg <= i__1; ++nstrg)
			{
				hstrng_1.pp[nstrg + nmom * 300 - 301] = 0.f;
			}
			i__1 = srec1_1.nst;
			for (nstrg = 1; nstrg <= i__1; ++nstrg)
			{
				hstrng_1.pt[nstrg + nmom * 300 - 301] = 0.f;
			}
		}
		i__1 = para1_1.mul;
		for (i__ = 1; i__ <= i__1; ++i__)
		{
			iityp = prec2_1.ityp5[i__ - 1];
			if (ilist8_1.lstrg1[i__ - 1] <= srec1_1.nsp)
			{
				nstrg = ilist8_1.lstrg1[i__ - 1];
				if (iityp == 2112 || iityp == 2212)
				{
					hstrng_1.pp[nstrg - 1] = (float)prec2_1.px5[i__ - 1];
					hstrng_1.pp[nstrg + 299] = (float)prec2_1.py5[i__ - 1];
					hstrng_1.pp[nstrg + 599] = (float)prec2_1.pz5[i__ - 1];
					hstrng_1.pp[nstrg + 899] = (float)prec2_1.e5[i__ - 1];
					hstrng_1.pp[nstrg + 1199] = (float)prec2_1.xmass5[i__ - 1];
				}
				else if ((iityp == 1 || iityp == 2) && (ilist8_1.lpart1[i__ - 1] == 1 || ilist8_1.lpart1[i__ - 1] == hjjet1_1.npj[nstrg - 1] + 2))
				{
					hstrng_1.pp[nstrg + 1499] = (float)prec2_1.px5[i__ - 1];
					hstrng_1.pp[nstrg + 1799] = (float)prec2_1.py5[i__ - 1];
					hstrng_1.pp[nstrg + 3899] = (float)prec2_1.xmass5[i__ - 1];
					hstrng_1.pp[nstrg - 1] += (float)prec2_1.px5[i__ - 1];
					hstrng_1.pp[nstrg + 299] += (float)prec2_1.py5[i__ - 1];
					hstrng_1.pp[nstrg + 599] += (float)prec2_1.pz5[i__ - 1];
					hstrng_1.pp[nstrg + 899] += (float)prec2_1.e5[i__ - 1];
					r__1 = hstrng_1.pp[nstrg + 899];
					r__2 = hstrng_1.pp[nstrg - 1];
					r__3 = hstrng_1.pp[nstrg + 299];
					r__4 = hstrng_1.pp[nstrg + 599];
					hstrng_1.pp[nstrg + 1199] = sqrt(r__1 * r__1 - r__2 * r__2 - r__3 * r__3 - r__4 * r__4);
				}
				else if ((iityp == 1103 || iityp == 2101 || iityp == 2103 ||
						  (float)iityp == 2203.f || iityp == 3101 || (float)iityp == 3103.f || iityp == 3201 || iityp == 3203 ||
						  iityp == 3303) &&
						 (ilist8_1.lpart1[i__ - 1] == 1 ||
						  ilist8_1.lpart1[i__ - 1] == hjjet1_1.npj[nstrg - 1] +
														  2))
				{
					hstrng_1.pp[nstrg + 2099] = (float)prec2_1.px5[i__ - 1];
					hstrng_1.pp[nstrg + 2399] = (float)prec2_1.py5[i__ - 1];
					hstrng_1.pp[nstrg + 4199] = (float)prec2_1.xmass5[i__ - 1];
					hstrng_1.pp[nstrg - 1] += (float)prec2_1.px5[i__ - 1];
					hstrng_1.pp[nstrg + 299] += (float)prec2_1.py5[i__ - 1];
					hstrng_1.pp[nstrg + 599] += (float)prec2_1.pz5[i__ - 1];
					hstrng_1.pp[nstrg + 899] += (float)prec2_1.e5[i__ - 1];
					r__1 = hstrng_1.pp[nstrg + 899];
					r__2 = hstrng_1.pp[nstrg - 1];
					r__3 = hstrng_1.pp[nstrg + 299];
					r__4 = hstrng_1.pp[nstrg + 599];
					hstrng_1.pp[nstrg + 1199] = sqrt(r__1 * r__1 - r__2 * r__2 - r__3 * r__3 - r__4 * r__4);
				}
				else
				{
					npart = ilist8_1.lpart1[i__ - 1] - 1;
					hjjet1_1.kfpj[nstrg + npart * 300 - 301] = prec2_1.ityp5[i__ - 1];
					hjjet1_1.pjpx[nstrg + npart * 300 - 301] = (float)
																   prec2_1.px5[i__ - 1];
					hjjet1_1.pjpy[nstrg + npart * 300 - 301] = (float)
																   prec2_1.py5[i__ - 1];
					hjjet1_1.pjpz[nstrg + npart * 300 - 301] = (float)
																   prec2_1.pz5[i__ - 1];
					hjjet1_1.pjpe[nstrg + npart * 300 - 301] = (float)
																   prec2_1.e5[i__ - 1];
					hjjet1_1.pjpm[nstrg + npart * 300 - 301] = (float)
																   prec2_1.xmass5[i__ - 1];
				}
			}
			else if (ilist8_1.lstrg1[i__ - 1] <= srec1_1.nsp + srec1_1.nst)
			{
				nstrg = ilist8_1.lstrg1[i__ - 1] - srec1_1.nsp;
				if (iityp == 2112 || iityp == 2212)
				{
					hstrng_1.pt[nstrg - 1] = (float)prec2_1.px5[i__ - 1];
					hstrng_1.pt[nstrg + 299] = (float)prec2_1.py5[i__ - 1];
					hstrng_1.pt[nstrg + 599] = (float)prec2_1.pz5[i__ - 1];
					hstrng_1.pt[nstrg + 899] = (float)prec2_1.e5[i__ - 1];
					hstrng_1.pt[nstrg + 1199] = (float)prec2_1.xmass5[i__ - 1];
				}
				else if ((iityp == 1 || iityp == 2) && (ilist8_1.lpart1[i__ - 1] == 1 || ilist8_1.lpart1[i__ - 1] == hjjet1_1.ntj[nstrg - 1] + 2))
				{
					hstrng_1.pt[nstrg + 1499] = (float)prec2_1.px5[i__ - 1];
					hstrng_1.pt[nstrg + 1799] = (float)prec2_1.py5[i__ - 1];
					hstrng_1.pt[nstrg + 3899] = (float)prec2_1.xmass5[i__ - 1];
					hstrng_1.pt[nstrg - 1] += (float)prec2_1.px5[i__ - 1];
					hstrng_1.pt[nstrg + 299] += (float)prec2_1.py5[i__ - 1];
					hstrng_1.pt[nstrg + 599] += (float)prec2_1.pz5[i__ - 1];
					hstrng_1.pt[nstrg + 899] += (float)prec2_1.e5[i__ - 1];
					r__1 = hstrng_1.pt[nstrg + 899];
					r__2 = hstrng_1.pt[nstrg - 1];
					r__3 = hstrng_1.pt[nstrg + 299];
					r__4 = hstrng_1.pt[nstrg + 599];
					hstrng_1.pt[nstrg + 1199] = sqrt(r__1 * r__1 - r__2 * r__2 - r__3 * r__3 - r__4 * r__4);
				}
				else if ((iityp == 1103 || iityp == 2101 || iityp == 2103 ||
						  (float)iityp == 2203.f || iityp == 3101 || (float)iityp == 3103.f || iityp == 3201 || iityp == 3203 ||
						  iityp == 3303) &&
						 (ilist8_1.lpart1[i__ - 1] == 1 ||
						  ilist8_1.lpart1[i__ - 1] == hjjet1_1.ntj[nstrg - 1] +
														  2))
				{
					hstrng_1.pt[nstrg + 2099] = (float)prec2_1.px5[i__ - 1];
					hstrng_1.pt[nstrg + 2399] = (float)prec2_1.py5[i__ - 1];
					hstrng_1.pt[nstrg + 4199] = (float)prec2_1.xmass5[i__ - 1];
					hstrng_1.pt[nstrg - 1] += (float)prec2_1.px5[i__ - 1];
					hstrng_1.pt[nstrg + 299] += (float)prec2_1.py5[i__ - 1];
					hstrng_1.pt[nstrg + 599] += (float)prec2_1.pz5[i__ - 1];
					hstrng_1.pt[nstrg + 899] += (float)prec2_1.e5[i__ - 1];
					r__1 = hstrng_1.pt[nstrg + 899];
					r__2 = hstrng_1.pt[nstrg - 1];
					r__3 = hstrng_1.pt[nstrg + 299];
					r__4 = hstrng_1.pt[nstrg + 599];
					hstrng_1.pt[nstrg + 1199] = sqrt(r__1 * r__1 - r__2 * r__2 - r__3 * r__3 - r__4 * r__4);
				}
				else
				{
					npart = ilist8_1.lpart1[i__ - 1] - 1;
					hjjet1_1.kftj[nstrg + npart * 300 - 301] = prec2_1.ityp5[i__ - 1];
					hjjet1_1.pjtx[nstrg + npart * 300 - 301] = (float)
																   prec2_1.px5[i__ - 1];
					hjjet1_1.pjty[nstrg + npart * 300 - 301] = (float)
																   prec2_1.py5[i__ - 1];
					hjjet1_1.pjtz[nstrg + npart * 300 - 301] = (float)
																   prec2_1.pz5[i__ - 1];
					hjjet1_1.pjte[nstrg + npart * 300 - 301] = (float)
																   prec2_1.e5[i__ - 1];
					hjjet1_1.pjtm[nstrg + npart * 300 - 301] = (float)
																   prec2_1.xmass5[i__ - 1];
				}
			}
			else
			{
				nstrg = ilist8_1.lstrg1[i__ - 1] - srec1_1.nsp - srec1_1.nst;
				npart = ilist8_1.lpart1[i__ - 1];
				hjjet2_1.k2sg[nstrg + npart * 150001 - 150002] =
					prec2_1.ityp5[i__ - 1];
				hjjet2_1.pxsg[nstrg + npart * 150001 - 150002] = (float)
																	 prec2_1.px5[i__ - 1];
				hjjet2_1.pysg[nstrg + npart * 150001 - 150002] = (float)
																	 prec2_1.py5[i__ - 1];
				hjjet2_1.pzsg[nstrg + npart * 150001 - 150002] = (float)
																	 prec2_1.pz5[i__ - 1];
				hjjet2_1.pesg[nstrg + npart * 150001 - 150002] = (float)
																	 prec2_1.e5[i__ - 1];
				hjjet2_1.pmsg[nstrg + npart * 150001 - 150002] = (float)
																	 prec2_1.xmass5[i__ - 1];
			}
		}
		ludat1_1.mstj[0] = 1;
		hparnt_1.ihpr2[0] = 1;
		anim_1.isflag = 1;
		hparnt_1.hipr1[0] = .94f;
		hjana2_();
	}
	else if (anim_1.isoft == 3 || anim_1.isoft == 4 || anim_1.isoft == 5)
	{
		anim_1.isflag = 0;
		if (hparnt_1.ihpr2[19] != 0)
		{
			i__1 = hjjet2_1.nsg;
			for (isg = 1; isg <= i__1; ++isg)
			{
				hijfrg_(&isg, &c__3, &ierror);

				nsbst = 1;
				idstr = 92;
				if (hparnt_1.ihpr2[20] == 0)
				{
					luedit_(&c__2);
				}
				else
				{
				L551:
					++nsbst;
					if (lujets_1.k[nsbst + 8999] < 91 || lujets_1.k[nsbst +
																	8999] > 93)
					{
						goto L551;
					}
					idstr = lujets_1.k[nsbst + 8999];
					++nsbst;
				}
				if (s_cmp(frame, "LAB", (int)8, (int)3) == 0)
				{
					hboost_();
				}
				nsbstr = 0;
				i__2 = lujets_1.n;
				for (i__ = nsbst; i__ <= i__2; ++i__)
				{
					if (lujets_1.k[i__ + 8999] == idstr)
					{
						++nsbstr;
					}
					lujets_1.k[i__ + 26999] = nsbstr;
					++hmain1_1.natt;
					hmain2_1.katt[hmain1_1.natt - 1] = lujets_1.k[i__ + 8999];
					hmain2_1.katt[hmain1_1.natt + 150000] = 20;
					hmain2_1.katt[hmain1_1.natt + 450002] = lujets_1.k[i__ -
																	   1];
					if (lujets_1.k[i__ + 17999] == 0)
					{
						hmain2_1.katt[hmain1_1.natt + 300001] = 0;
					}
					else if (lujets_1.k[i__ + 17999] != 0 && lujets_1.k[lujets_1.k[i__ + 17999] + 8999] == idstr)
					{
						hmain2_1.katt[hmain1_1.natt + 300001] = 0;
					}
					else
					{
						hmain2_1.katt[hmain1_1.natt + 300001] = hmain1_1.natt - i__ + lujets_1.k[i__ + 17999] + nsbstr -
																lujets_1.k[lujets_1.k[i__ + 17999] + 26999];
					}
					hmain2_1.patt[hmain1_1.natt - 1] = lujets_1.p[i__ - 1];
					hmain2_1.patt[hmain1_1.natt + 150000] = lujets_1.p[i__ +
																	   8999];
					hmain2_1.patt[hmain1_1.natt + 300001] = lujets_1.p[i__ +
																	   17999];
					hmain2_1.patt[hmain1_1.natt + 450002] = lujets_1.p[i__ +
																	   26999];
					hmain1_1.eatt += lujets_1.p[i__ + 26999];
					arprc_1.gxar[hmain1_1.natt - 1] = (hjcrdn_1.yp[hjjet2_1.iasg[isg - 1] * 3 - 3] + hjcrdn_1.yt[hjjet2_1.iasg[isg + 150000] * 3 - 3]) * .5f;
					arprc_1.gyar[hmain1_1.natt - 1] = (hjcrdn_1.yp[hjjet2_1.iasg[isg - 1] * 3 - 2] + hjcrdn_1.yt[hjjet2_1.iasg[isg + 150000] * 3 - 2]) * .5f;
					arprc_1.gzar[hmain1_1.natt - 1] = 0.f;
					arprc_1.ftar[hmain1_1.natt - 1] = 0.f;
					arprc_1.itypar[hmain1_1.natt - 1] = lujets_1.k[i__ + 8999];
					arprc_1.pxar[hmain1_1.natt - 1] = lujets_1.p[i__ - 1];
					arprc_1.pyar[hmain1_1.natt - 1] = lujets_1.p[i__ + 8999];
					arprc_1.pzar[hmain1_1.natt - 1] = lujets_1.p[i__ + 17999];
					arprc_1.pear[hmain1_1.natt - 1] = lujets_1.p[i__ + 26999];
					arprc_1.xmar[hmain1_1.natt - 1] = lujets_1.p[i__ + 35999];
				}
			}
			jtp[0] = hparnt_1.ihnt2[0];
			jtp[1] = hparnt_1.ihnt2[2];
			for (ntp = 1; ntp <= 2; ++ntp)
			{
				i__2 = jtp[ntp - 1];
				for (jjtp = 1; jjtp <= i__2; ++jjtp)
				{
					hijfrg_(&jjtp, &ntp, &ierror);

					nsbst = 1;
					idstr = 92;
					if (hparnt_1.ihpr2[20] == 0)
					{
						luedit_(&c__2);
					}
					else
					{
					L581:
						++nsbst;
						if (lujets_1.k[nsbst + 8999] < 91 || lujets_1.k[nsbst + 8999] > 93)
						{
							goto L581;
						}
						idstr = lujets_1.k[nsbst + 8999];
						++nsbst;
					}
					if (s_cmp(frame, "LAB", (int)8, (int)3) == 0)
					{
						hboost_();
					}
					nftp = hstrng_1.nfp[jjtp + 1199];
					if (ntp == 2)
					{
						nftp = hstrng_1.nft[jjtp + 1199] + 10;
					}
					nsbstr = 0;
					i__1 = lujets_1.n;
					for (i__ = nsbst; i__ <= i__1; ++i__)
					{
						if (lujets_1.k[i__ + 8999] == idstr)
						{
							++nsbstr;
						}
						lujets_1.k[i__ + 26999] = nsbstr;
						++hmain1_1.natt;
						hmain2_1.katt[hmain1_1.natt - 1] = lujets_1.k[i__ +
																	  8999];
						hmain2_1.katt[hmain1_1.natt + 150000] = nftp;
						hmain2_1.katt[hmain1_1.natt + 450002] = lujets_1.k[i__ - 1];
						if (lujets_1.k[i__ + 17999] == 0)
						{
							hmain2_1.katt[hmain1_1.natt + 300001] = 0;
						}
						else if (lujets_1.k[i__ + 17999] != 0 && lujets_1.k[lujets_1.k[i__ + 17999] + 8999] == idstr)
						{
							hmain2_1.katt[hmain1_1.natt + 300001] = 0;
						}
						else
						{
							hmain2_1.katt[hmain1_1.natt + 300001] =
								hmain1_1.natt - i__ + lujets_1.k[i__ + 17999] + nsbstr - lujets_1.k[lujets_1.k[i__ + 17999] + 26999];
						}
						hmain2_1.patt[hmain1_1.natt - 1] = lujets_1.p[i__ - 1];
						hmain2_1.patt[hmain1_1.natt + 150000] = lujets_1.p[i__ + 8999];
						hmain2_1.patt[hmain1_1.natt + 300001] = lujets_1.p[i__ + 17999];
						hmain2_1.patt[hmain1_1.natt + 450002] = lujets_1.p[i__ + 26999];
						hmain1_1.eatt += lujets_1.p[i__ + 26999];
						if (ntp == 1)
						{
							arprc_1.gxar[hmain1_1.natt - 1] = hjcrdn_1.yp[jjtp * 3 - 3] + bb * .5f;
							arprc_1.gyar[hmain1_1.natt - 1] = hjcrdn_1.yp[jjtp * 3 - 2];
						}
						else
						{
							arprc_1.gxar[hmain1_1.natt - 1] = hjcrdn_1.yt[jjtp * 3 - 3] - bb * .5f;
							arprc_1.gyar[hmain1_1.natt - 1] = hjcrdn_1.yt[jjtp * 3 - 2];
						}
						arprc_1.gzar[hmain1_1.natt - 1] = 0.f;
						arprc_1.ftar[hmain1_1.natt - 1] = 0.f;
						arprc_1.itypar[hmain1_1.natt - 1] = lujets_1.k[i__ +
																	   8999];
						arprc_1.pxar[hmain1_1.natt - 1] = lujets_1.p[i__ - 1];
						arprc_1.pyar[hmain1_1.natt - 1] = lujets_1.p[i__ +
																	 8999];
						arprc_1.pzar[hmain1_1.natt - 1] = lujets_1.p[i__ +
																	 17999];
						arprc_1.pear[hmain1_1.natt - 1] = lujets_1.p[i__ +
																	 26999];
						arprc_1.xmar[hmain1_1.natt - 1] = lujets_1.p[i__ +
																	 35999];
					}
				}
			}
		}
		if (hjjet4_1.ndr >= 1)
		{

			i__2 = hjjet4_1.ndr;
			for (i__ = 1; i__ <= i__2; ++i__)
			{
				++hmain1_1.natt;
				hmain2_1.katt[hmain1_1.natt - 1] = hjjet4_1.kfdr[i__ - 1];
				hmain2_1.katt[hmain1_1.natt + 150000] = 40;
				hmain2_1.katt[hmain1_1.natt + 300001] = 0;
				hmain2_1.patt[hmain1_1.natt - 1] = hjjet4_1.pdr[i__ - 1];
				hmain2_1.patt[hmain1_1.natt + 150000] = hjjet4_1.pdr[i__ +
																	 899];
				hmain2_1.patt[hmain1_1.natt + 300001] = hjjet4_1.pdr[i__ +
																	 1799];
				hmain2_1.patt[hmain1_1.natt + 450002] = hjjet4_1.pdr[i__ +
																	 2699];
				hmain1_1.eatt += hjjet4_1.pdr[i__ + 2699];
				arprc_1.gxar[hmain1_1.natt - 1] = xydr_1.rtdr[i__ - 1];
				arprc_1.gyar[hmain1_1.natt - 1] = xydr_1.rtdr[i__ + 150000];
				arprc_1.gzar[hmain1_1.natt - 1] = 0.f;
				arprc_1.ftar[hmain1_1.natt - 1] = 0.f;
				arprc_1.itypar[hmain1_1.natt - 1] = hmain2_1.katt[hmain1_1.natt - 1];
				arprc_1.pxar[hmain1_1.natt - 1] = hmain2_1.patt[hmain1_1.natt - 1];
				arprc_1.pyar[hmain1_1.natt - 1] = hmain2_1.patt[hmain1_1.natt + 150000];
				arprc_1.pzar[hmain1_1.natt - 1] = hmain2_1.patt[hmain1_1.natt + 300001];
				arprc_1.pear[hmain1_1.natt - 1] = hmain2_1.patt[hmain1_1.natt + 450002];
				arprc_1.xmar[hmain1_1.natt - 1] = hjjet4_1.pdr[i__ + 3599];
			}
		}

		hjana1_();
		htop_();
		srec1_1.nsp = 0;
		srec1_1.nst = 0;
		hjjet2_1.nsg = hmain1_1.natt;
		srec1_1.nsi = hjjet2_1.nsg;
		zpcmn_();
		s_wsfe(&io___93);
		do_fio(&c__1, (char *)&itest, (int)sizeof(int));
		do_fio(&c__1, (char *)&para1_1.mul, (int)sizeof(int));
		do_fio(&c__1, (char *)&lastt_1.bimp, (int)sizeof(float));
		do_fio(&c__1, (char *)&hjglbr_1.nelp, (int)sizeof(int));
		do_fio(&c__1, (char *)&hjglbr_1.ninp, (int)sizeof(int));
		do_fio(&c__1, (char *)&hjglbr_1.nelt, (int)sizeof(int));
		do_fio(&c__1, (char *)&hjglbr_1.ninthj, (int)sizeof(int));
		e_wsfe();
		++itest;
		i__2 = para1_1.mul;
		for (i__ = 1; i__ <= i__2; ++i__)
		{
			s_wsfe(&io___94);
			do_fio(&c__1, (char *)&prec2_1.px5[i__ - 1], (int)sizeof(double));
			do_fio(&c__1, (char *)&prec2_1.py5[i__ - 1], (int)sizeof(double));
			do_fio(&c__1, (char *)&prec2_1.pz5[i__ - 1], (int)sizeof(double));
			do_fio(&c__1, (char *)&prec2_1.ityp5[i__ - 1], (int)sizeof(int));
			do_fio(&c__1, (char *)&prec2_1.xmass5[i__ - 1], (int)sizeof(double));
			do_fio(&c__1, (char *)&prec2_1.e5[i__ - 1], (int)sizeof(double));
			do_fio(&c__1, (char *)&ilist8_1.lstrg1[i__ - 1], (int)sizeof(int));
			do_fio(&c__1, (char *)&ilist8_1.lpart1[i__ - 1], (int)sizeof(int));
			do_fio(&c__1, (char *)&prec2_1.ft5[i__ - 1], (int)sizeof(double));
			e_wsfe();
		}
		for (i__ = 1; i__ <= 150001; ++i__)
		{
			for (j = 1; j <= 3; ++j)
			{
				soft_1.k1sgs[i__ + j * 150001 - 150002] = 0;
				soft_1.k2sgs[i__ + j * 150001 - 150002] = 0;
				soft_1.pxsgs[i__ + j * 150001 - 150002] = 0.;
				soft_1.pysgs[i__ + j * 150001 - 150002] = 0.;
				soft_1.pzsgs[i__ + j * 150001 - 150002] = 0.;
				soft_1.pesgs[i__ + j * 150001 - 150002] = 0.;
				soft_1.pmsgs[i__ + j * 150001 - 150002] = 0.;
				soft_1.gxsgs[i__ + j * 150001 - 150002] = 0.;
				soft_1.gysgs[i__ + j * 150001 - 150002] = 0.;
				soft_1.gzsgs[i__ + j * 150001 - 150002] = 0.;
				soft_1.ftsgs[i__ + j * 150001 - 150002] = 0.;
			}
		}
		i__2 = para1_1.mul;
		for (i__ = 1; i__ <= i__2; ++i__)
		{
			iityp = prec2_1.ityp5[i__ - 1];
			nstrg = ilist8_1.lstrg1[i__ - 1];
			npart = ilist8_1.lpart1[i__ - 1];
			soft_1.k2sgs[nstrg + npart * 150001 - 150002] = prec2_1.ityp5[i__ - 1];
			soft_1.pxsgs[nstrg + npart * 150001 - 150002] = prec2_1.px5[i__ -
																		1];
			soft_1.pysgs[nstrg + npart * 150001 - 150002] = prec2_1.py5[i__ -
																		1];
			soft_1.pzsgs[nstrg + npart * 150001 - 150002] = prec2_1.pz5[i__ -
																		1];
			soft_1.pmsgs[nstrg + npart * 150001 - 150002] = prec2_1.xmass5[i__ - 1];
			d__1 = prec2_1.px5[i__ - 1];
			d__2 = prec2_1.py5[i__ - 1];
			d__3 = prec2_1.pz5[i__ - 1];
			d__4 = prec2_1.xmass5[i__ - 1];
			prec2_1.e5[i__ - 1] = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3 + d__4 * d__4);
			soft_1.pesgs[nstrg + npart * 150001 - 150002] = prec2_1.e5[i__ -
																	   1];
			soft_1.gxsgs[nstrg + npart * 150001 - 150002] = prec2_1.gx5[i__ -
																		1];
			soft_1.gysgs[nstrg + npart * 150001 - 150002] = prec2_1.gy5[i__ -
																		1];
			soft_1.gzsgs[nstrg + npart * 150001 - 150002] = prec2_1.gz5[i__ -
																		1];
			soft_1.ftsgs[nstrg + npart * 150001 - 150002] = prec2_1.ft5[i__ -
																		1];
		}
		hjana2_();
	}
	if (anim_1.isoft == 3 || anim_1.isoft == 4 || anim_1.isoft == 5)
	{
		hmain1_1.natt = 0;
		hmain1_1.eatt = 0.f;
		ptoh_();
		i__2 = noprec_1.nnozpc;
		for (i__ = 1; i__ <= i__2; ++i__)
		{
			++hmain1_1.natt;
			hmain2_1.katt[hmain1_1.natt - 1] = noprec_1.itypn[i__ - 1];
			hmain2_1.patt[hmain1_1.natt - 1] = noprec_1.pxn[i__ - 1];
			hmain2_1.patt[hmain1_1.natt + 150000] = noprec_1.pyn[i__ - 1];
			hmain2_1.patt[hmain1_1.natt + 300001] = noprec_1.pzn[i__ - 1];
			hmain2_1.patt[hmain1_1.natt + 450002] = noprec_1.een[i__ - 1];
			hmain1_1.eatt += noprec_1.een[i__ - 1];
			arprc_1.gxar[hmain1_1.natt - 1] = noprec_1.gxn[i__ - 1];
			arprc_1.gyar[hmain1_1.natt - 1] = noprec_1.gyn[i__ - 1];
			arprc_1.gzar[hmain1_1.natt - 1] = noprec_1.gzn[i__ - 1];
			arprc_1.ftar[hmain1_1.natt - 1] = noprec_1.ftn[i__ - 1];
			arprc_1.itypar[hmain1_1.natt - 1] = noprec_1.itypn[i__ - 1];
			arprc_1.pxar[hmain1_1.natt - 1] = noprec_1.pxn[i__ - 1];
			arprc_1.pyar[hmain1_1.natt - 1] = noprec_1.pyn[i__ - 1];
			arprc_1.pzar[hmain1_1.natt - 1] = noprec_1.pzn[i__ - 1];
			arprc_1.pear[hmain1_1.natt - 1] = noprec_1.een[i__ - 1];
			arprc_1.xmar[hmain1_1.natt - 1] = noprec_1.xmn[i__ - 1];
		}
		goto L565;
	}
	if (hparnt_1.ihpr2[19] != 0)
	{
		i__2 = hjjet2_1.nsg;
		for (isg = 1; isg <= i__2; ++isg)
		{
			hijfrg_(&isg, &c__3, &ierror);
			if (ludat1_1.mstu[23] != 0 || ierror > 0)
			{
				ludat1_1.mstu[23] = 0;
				ludat1_1.mstu[27] = 0;
				if (hparnt_1.ihpr2[9] != 0)
				{
					lulist_(&c__2);
					s_wsle(&io___95);
					do_lio(&c__9, &c__1, "error occured ISG, repeat the event", (int)35);
					e_wsle();
					s_wsle(&io___96);
					do_lio(&c__3, &c__1, (char *)&isg, (int)sizeof(int));
					e_wsle();
				}
				goto L50;
			}
			nsbst = 1;
			idstr = 92;
			if (hparnt_1.ihpr2[20] == 0)
			{
				luedit_(&c__2);
			}
			else
			{
			L351:
				++nsbst;
				if (lujets_1.k[nsbst + 8999] < 91 || lujets_1.k[nsbst + 8999] > 93)
				{
					goto L351;
				}
				idstr = lujets_1.k[nsbst + 8999];
				++nsbst;
			}

			if (s_cmp(frame, "LAB", (int)8, (int)3) == 0)
			{
				hboost_();
			}
			nsbstr = 0;
			i__1 = lujets_1.n;
			for (i__ = nsbst; i__ <= i__1; ++i__)
			{
				if (lujets_1.k[i__ + 8999] == idstr)
				{
					++nsbstr;
				}
				lujets_1.k[i__ + 26999] = nsbstr;
				++hmain1_1.natt;
				hmain2_1.katt[hmain1_1.natt - 1] = lujets_1.k[i__ + 8999];
				hmain2_1.katt[hmain1_1.natt + 150000] = 20;
				hmain2_1.katt[hmain1_1.natt + 450002] = lujets_1.k[i__ - 1];
				if (lujets_1.k[i__ + 17999] == 0)
				{
					hmain2_1.katt[hmain1_1.natt + 300001] = 0;
				}
				else if (lujets_1.k[i__ + 17999] != 0 && lujets_1.k[lujets_1.k[i__ + 17999] + 8999] == idstr)
				{
					hmain2_1.katt[hmain1_1.natt + 300001] = 0;
				}
				else
				{
					hmain2_1.katt[hmain1_1.natt + 300001] = hmain1_1.natt -
															i__ + lujets_1.k[i__ + 17999] + nsbstr -
															lujets_1.k[lujets_1.k[i__ + 17999] + 26999];
				}
				hmain2_1.patt[hmain1_1.natt - 1] = lujets_1.p[i__ - 1];
				hmain2_1.patt[hmain1_1.natt + 150000] = lujets_1.p[i__ + 8999];
				hmain2_1.patt[hmain1_1.natt + 300001] = lujets_1.p[i__ +
																   17999];
				hmain2_1.patt[hmain1_1.natt + 450002] = lujets_1.p[i__ +
																   26999];
				hmain1_1.eatt += lujets_1.p[i__ + 26999];
				lsg = srec1_1.nsp + srec1_1.nst + isg;
				arprc_1.gxar[hmain1_1.natt - 1] = (float)srec2_1.zt1[lsg - 1];
				arprc_1.gyar[hmain1_1.natt - 1] = (float)srec2_1.zt2[lsg - 1];
				arprc_1.gzar[hmain1_1.natt - 1] = (float)srec2_1.zt3[lsg - 1];
				arprc_1.ftar[hmain1_1.natt - 1] = (float)srec2_1.ataui[lsg -
																	   1];
				arprc_1.itypar[hmain1_1.natt - 1] = lujets_1.k[i__ + 8999];
				arprc_1.pxar[hmain1_1.natt - 1] = lujets_1.p[i__ - 1];
				arprc_1.pyar[hmain1_1.natt - 1] = lujets_1.p[i__ + 8999];
				arprc_1.pzar[hmain1_1.natt - 1] = lujets_1.p[i__ + 17999];
				arprc_1.pear[hmain1_1.natt - 1] = lujets_1.p[i__ + 26999];
				arprc_1.xmar[hmain1_1.natt - 1] = lujets_1.p[i__ + 35999];
			}
		}
		jtp[0] = hparnt_1.ihnt2[0];
		jtp[1] = hparnt_1.ihnt2[2];
		for (ntp = 1; ntp <= 2; ++ntp)
		{
			i__1 = jtp[ntp - 1];
			for (jjtp = 1; jjtp <= i__1; ++jjtp)
			{
				hijfrg_(&jjtp, &ntp, &ierror);
				if (ludat1_1.mstu[23] != 0 || ierror > 0)
				{
					ludat1_1.mstu[23] = 0;
					ludat1_1.mstu[27] = 0;
					if (hparnt_1.ihpr2[9] != 0)
					{
						lulist_(&c__2);
						s_wsle(&io___98);
						do_lio(&c__9, &c__1, "error occured P&T, repeat the "
											 "event",
							   (int)35);
						e_wsle();
						s_wsle(&io___99);
						do_lio(&c__3, &c__1, (char *)&ntp, (int)sizeof(int));
						do_lio(&c__3, &c__1, (char *)&jjtp, (int)sizeof(int));
						e_wsle();
					}
					goto L50;
				}
				nsbst = 1;
				idstr = 92;
				if (hparnt_1.ihpr2[20] == 0)
				{
					luedit_(&c__2);
				}
				else
				{
				L381:
					++nsbst;
					if (lujets_1.k[nsbst + 8999] < 91 || lujets_1.k[nsbst +
																	8999] > 93)
					{
						goto L381;
					}
					idstr = lujets_1.k[nsbst + 8999];
					++nsbst;
				}
				if (s_cmp(frame, "LAB", (int)8, (int)3) == 0)
				{
					hboost_();
				}
				nftp = hstrng_1.nfp[jjtp + 1199];
				if (ntp == 2)
				{
					nftp = hstrng_1.nft[jjtp + 1199] + 10;
				}
				nsbstr = 0;
				i__2 = lujets_1.n;
				for (i__ = nsbst; i__ <= i__2; ++i__)
				{
					if (lujets_1.k[i__ + 8999] == idstr)
					{
						++nsbstr;
					}
					lujets_1.k[i__ + 26999] = nsbstr;
					++hmain1_1.natt;
					hmain2_1.katt[hmain1_1.natt - 1] = lujets_1.k[i__ + 8999];
					hmain2_1.katt[hmain1_1.natt + 150000] = nftp;
					hmain2_1.katt[hmain1_1.natt + 450002] = lujets_1.k[i__ -
																	   1];
					if (lujets_1.k[i__ + 17999] == 0)
					{
						hmain2_1.katt[hmain1_1.natt + 300001] = 0;
					}
					else if (lujets_1.k[i__ + 17999] != 0 && lujets_1.k[lujets_1.k[i__ + 17999] + 8999] == idstr)
					{
						hmain2_1.katt[hmain1_1.natt + 300001] = 0;
					}
					else
					{
						hmain2_1.katt[hmain1_1.natt + 300001] = hmain1_1.natt - i__ + lujets_1.k[i__ + 17999] + nsbstr -
																lujets_1.k[lujets_1.k[i__ + 17999] + 26999];
					}
					hmain2_1.patt[hmain1_1.natt - 1] = lujets_1.p[i__ - 1];
					hmain2_1.patt[hmain1_1.natt + 150000] = lujets_1.p[i__ +
																	   8999];
					hmain2_1.patt[hmain1_1.natt + 300001] = lujets_1.p[i__ +
																	   17999];
					hmain2_1.patt[hmain1_1.natt + 450002] = lujets_1.p[i__ +
																	   26999];
					hmain1_1.eatt += lujets_1.p[i__ + 26999];
					if (ntp == 1)
					{
						lsg = jjtp;
					}
					else
					{
						lsg = jjtp + srec1_1.nsp;
					}
					arprc_1.gxar[hmain1_1.natt - 1] = (float)srec2_1.zt1[lsg - 1];
					arprc_1.gyar[hmain1_1.natt - 1] = (float)srec2_1.zt2[lsg - 1];
					arprc_1.gzar[hmain1_1.natt - 1] = (float)srec2_1.zt3[lsg - 1];
					arprc_1.ftar[hmain1_1.natt - 1] = (float)srec2_1.ataui[lsg - 1];
					arprc_1.itypar[hmain1_1.natt - 1] = lujets_1.k[i__ + 8999];
					arprc_1.pxar[hmain1_1.natt - 1] = lujets_1.p[i__ - 1];
					arprc_1.pyar[hmain1_1.natt - 1] = lujets_1.p[i__ + 8999];
					arprc_1.pzar[hmain1_1.natt - 1] = lujets_1.p[i__ + 17999];
					arprc_1.pear[hmain1_1.natt - 1] = lujets_1.p[i__ + 26999];
					arprc_1.xmar[hmain1_1.natt - 1] = lujets_1.p[i__ + 35999];
				}
			}
		}
	}
	i__1 = hjjet4_1.ndr;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		++hmain1_1.natt;
		hmain2_1.katt[hmain1_1.natt - 1] = hjjet4_1.kfdr[i__ - 1];
		hmain2_1.katt[hmain1_1.natt + 150000] = 40;
		hmain2_1.katt[hmain1_1.natt + 300001] = 0;
		hmain2_1.patt[hmain1_1.natt - 1] = hjjet4_1.pdr[i__ - 1];
		hmain2_1.patt[hmain1_1.natt + 150000] = hjjet4_1.pdr[i__ + 899];
		hmain2_1.patt[hmain1_1.natt + 300001] = hjjet4_1.pdr[i__ + 1799];
		hmain2_1.patt[hmain1_1.natt + 450002] = hjjet4_1.pdr[i__ + 2699];
		hmain1_1.eatt += hjjet4_1.pdr[i__ + 2699];
		arprc_1.gxar[hmain1_1.natt - 1] = xydr_1.rtdr[i__ - 1];
		arprc_1.gyar[hmain1_1.natt - 1] = xydr_1.rtdr[i__ + 150000];
		arprc_1.gzar[hmain1_1.natt - 1] = 0.f;
		arprc_1.ftar[hmain1_1.natt - 1] = 0.f;
		arprc_1.itypar[hmain1_1.natt - 1] = hmain2_1.katt[hmain1_1.natt - 1];
		arprc_1.pxar[hmain1_1.natt - 1] = hmain2_1.patt[hmain1_1.natt - 1];
		arprc_1.pyar[hmain1_1.natt - 1] = hmain2_1.patt[hmain1_1.natt +
														150000];
		arprc_1.pzar[hmain1_1.natt - 1] = hmain2_1.patt[hmain1_1.natt +
														300001];
		arprc_1.pear[hmain1_1.natt - 1] = hmain2_1.patt[hmain1_1.natt +
														450002];
		arprc_1.xmar[hmain1_1.natt - 1] = hjjet4_1.pdr[i__ + 3599];
	}
L565:
	dengy = hmain1_1.eatt / (hparnt_1.ihnt2[0] * hparnt_1.hint1[5] +
							 hparnt_1.ihnt2[2] * hparnt_1.hint1[6]) -
			1.f;
	if (dabs(dengy) > hparnt_1.hipr1[42] && hparnt_1.ihpr2[19] != 0 &&
		hparnt_1.ihpr2[20] == 0)
	{
		if (hparnt_1.ihpr2[9] != 0)
		{
			s_wsle(&io___101);
			do_lio(&c__9, &c__1, "Energy not conserved, repeat the event", (int)38);
			e_wsle();
		}
		s_wsle(&io___102);
		do_lio(&c__9, &c__1, "violated:EATT,NATT,B=", (int)21);
		do_lio(&c__4, &c__1, (char *)&hmain1_1.eatt, (int)sizeof(float));
		do_lio(&c__3, &c__1, (char *)&hmain1_1.natt, (int)sizeof(int));
		do_lio(&c__4, &c__1, (char *)&lastt_1.bimp, (int)sizeof(float));
		e_wsle();
		goto L50;
	}
	s_wsle(&io___103);
	do_lio(&c__9, &c__1, "satisfied:EATT,NATT,B=", (int)22);
	do_lio(&c__4, &c__1, (char *)&hmain1_1.eatt, (int)sizeof(float));
	do_lio(&c__3, &c__1, (char *)&hmain1_1.natt, (int)sizeof(int));
	do_lio(&c__4, &c__1, (char *)&lastt_1.bimp, (int)sizeof(float));
	e_wsle();
	s_wsle(&io___104);
	do_lio(&c__9, &c__1, " ", (int)1);
	e_wsle();
	return 0;
}

int hijset_(float *efrm, char *frame, char *proj, char *targ,
			int *iap, int *izp, int *iat, int *izt, int frame_len, int proj_len, int targ_len)
{
	static char fmt_100[] = "(//10x,\002*********************************"
							"*\002//\002****************\002/10x,\002*\002,48x,\002*\002/10x"
							",\002*         HIJING has been initialized at         *\002/10x"
							",\002*\002,13x,a4,\002= \002,f10.2,\002 GeV/n\002,13x,\002*\002/"
							"10x,\002*\002,48x,\002*\002/10x,\002*\002,8x,\002for \002,a4,"
							"\002(\002,i3,\002,\002,i3,\002)\002,\002 + \002,a4,\002(\002,i3"
							",\002,\002,i3,\002)\002,7x,\002*\002/10x,\002*******************"
							"*******************************\002)";

	float r__1, r__2;
	double d__1, d__2;

	int s_cmp(char *, char *, int, int), s_wsle(cilist *), do_lio(int *, int *, char *, int), e_wsle(void);
	int s_stop(char *, int);
	double pow_dd(double *, double *), sqrt(double), log(double);
	int s_copy(char *, char *, int, int);
	int s_wsfe(cilist *), do_fio(int *, char *, int), e_wsfe(void);

	static int i__, j;
	static double dd1, dd2, dd3, dd4;
	static float rkp, rmax;
	extern int fnkc2_();
	extern int hifun_(int *, float *, float *, int *),
		title_(void);
	static char eframe[4];
	extern int fnkick_();
	extern int hijcrs_(void), hijwds_(int *, int *, float *), lugive_(char *, int);
	extern double ulmass_(int *);
	extern int fnstrm_(), fnstrs_(), fnstru_();

	static cilist io___105 = {0, 6, 0, 0, 0};
	static cilist io___106 = {0, 6, 0, 0, 0};
	static cilist io___116 = {0, 6, 0, fmt_100, 0};

	title_();
	hparnt_1.ihnt2[0] = *iap;
	hparnt_1.ihnt2[1] = *izp;
	hparnt_1.ihnt2[2] = *iat;
	hparnt_1.ihnt2[3] = *izt;
	hparnt_1.ihnt2[4] = 0;
	hparnt_1.ihnt2[5] = 0;

	r__1 = ulmass_(&c__2112), r__2 = ulmass_(&c__2212);
	hparnt_1.hint1[7] = max(r__1, r__2);
	hparnt_1.hint1[8] = hparnt_1.hint1[7];

	if (s_cmp(proj, "A", (int)4, (int)1) != 0)
	{
		if (s_cmp(proj, "P", (int)4, (int)1) == 0)
		{
			hparnt_1.ihnt2[4] = 2212;
		}
		else if (s_cmp(proj, "PBAR", (int)4, (int)4) == 0)
		{
			hparnt_1.ihnt2[4] = -2212;
		}
		else if (s_cmp(proj, "PI+", (int)4, (int)3) == 0)
		{
			hparnt_1.ihnt2[4] = 211;
		}
		else if (s_cmp(proj, "PI-", (int)4, (int)3) == 0)
		{
			hparnt_1.ihnt2[4] = -211;
		}
		else if (s_cmp(proj, "K+", (int)4, (int)2) == 0)
		{
			hparnt_1.ihnt2[4] = 321;
		}
		else if (s_cmp(proj, "K-", (int)4, (int)2) == 0)
		{
			hparnt_1.ihnt2[4] = -321;
		}
		else if (s_cmp(proj, "N", (int)4, (int)1) == 0)
		{
			hparnt_1.ihnt2[4] = 2112;
		}
		else if (s_cmp(proj, "NBAR", (int)4, (int)4) == 0)
		{
			hparnt_1.ihnt2[4] = -2112;
		}
		else
		{
			s_wsle(&io___105);
			do_lio(&c__9, &c__1, proj, (int)4);
			do_lio(&c__9, &c__1, "wrong or unavailable proj name", (int)30);
			e_wsle();
			s_stop("", (int)0);
		}
		hparnt_1.hint1[7] = ulmass_(&hparnt_1.ihnt2[4]);
	}
	if (s_cmp(targ, "A", (int)4, (int)1) != 0)
	{
		if (s_cmp(targ, "P", (int)4, (int)1) == 0)
		{
			hparnt_1.ihnt2[5] = 2212;
		}
		else if (s_cmp(targ, "PBAR", (int)4, (int)4) == 0)
		{
			hparnt_1.ihnt2[5] = -2212;
		}
		else if (s_cmp(targ, "PI+", (int)4, (int)3) == 0)
		{
			hparnt_1.ihnt2[5] = 211;
		}
		else if (s_cmp(targ, "PI-", (int)4, (int)3) == 0)
		{
			hparnt_1.ihnt2[5] = -211;
		}
		else if (s_cmp(targ, "K+", (int)4, (int)2) == 0)
		{
			hparnt_1.ihnt2[5] = 321;
		}
		else if (s_cmp(targ, "K-", (int)4, (int)2) == 0)
		{
			hparnt_1.ihnt2[5] = -321;
		}
		else if (s_cmp(targ, "N", (int)4, (int)1) == 0)
		{
			hparnt_1.ihnt2[5] = 2112;
		}
		else if (s_cmp(targ, "NBAR", (int)4, (int)4) == 0)
		{
			hparnt_1.ihnt2[5] = -2112;
		}
		else
		{
			s_wsle(&io___106);
			do_lio(&c__9, &c__1, targ, (int)4);
			do_lio(&c__9, &c__1, "wrong or unavailable targ name", (int)30);
			e_wsle();
			s_stop("", (int)0);
		}
		hparnt_1.hint1[8] = ulmass_(&hparnt_1.ihnt2[5]);
	}
	if (hparnt_1.ihpr2[11] > 0)
	{
		lugive_("MDCY(C221,1)=0", (int)14);
		lugive_("MDCY(C313,1)=0", (int)14);
		lugive_("MDCY(C-313,1)=0", (int)15);
		lugive_("MDCY(C323,1)=0", (int)14);
		lugive_("MDCY(C-323,1)=0", (int)15);
		lugive_("MDCY(C311,1)=0", (int)14);
		lugive_("MDCY(C-311,1)=0", (int)15);
		lugive_("MDCY(C1114,1)=0", (int)15);
		lugive_("MDCY(C2114,1)=0", (int)15);
		lugive_("MDCY(C2214,1)=0", (int)15);
		lugive_("MDCY(C2224,1)=0", (int)15);
		lugive_("MDCY(C-1114,1)=0", (int)16);
		lugive_("MDCY(C-2114,1)=0", (int)16);
		lugive_("MDCY(C-2214,1)=0", (int)16);
		lugive_("MDCY(C-2224,1)=0", (int)16);
		lugive_("MDCY(C213,1)=0", (int)14);
		lugive_("MDCY(C-213,1)=0", (int)15);
		lugive_("MDCY(C113,1)=0", (int)14);
		lugive_("MDCY(C223,1)=0", (int)14);
		lugive_("MDCY(C333,1)=0", (int)14);
		lugive_("MDCY(C111,1)=0", (int)14);
		lugive_("MDCY(C310,1)=0", (int)14);
		lugive_("MDCY(C411,1)=0;MDCY(C-411,1)=0", (int)30);
		lugive_("MDCY(C421,1)=0;MDCY(C-421,1)=0", (int)30);
		lugive_("MDCY(C431,1)=0;MDCY(C-431,1)=0", (int)30);
		lugive_("MDCY(C511,1)=0;MDCY(C-511,1)=0", (int)30);
		lugive_("MDCY(C521,1)=0;MDCY(C-521,1)=0", (int)30);
		lugive_("MDCY(C531,1)=0;MDCY(C-531,1)=0", (int)30);
		lugive_("MDCY(C3122,1)=0;MDCY(C-3122,1)=0", (int)32);
		lugive_("MDCY(C3112,1)=0;MDCY(C-3112,1)=0", (int)32);
		lugive_("MDCY(C3212,1)=0;MDCY(C-3212,1)=0", (int)32);
		lugive_("MDCY(C3222,1)=0;MDCY(C-3222,1)=0", (int)32);
		lugive_("MDCY(C3312,1)=0;MDCY(C-3312,1)=0", (int)32);
		lugive_("MDCY(C3322,1)=0;MDCY(C-3322,1)=0", (int)32);
		lugive_("MDCY(C3334,1)=0;MDCY(C-3334,1)=0", (int)32);
	}
	ludat1_1.mstu[11] = 0;
	ludat1_1.mstu[20] = 1;
	if (hparnt_1.ihpr2[9] == 0)
	{
		ludat1_1.mstu[21] = 0;
		ludat1_1.mstu[24] = 0;
		ludat1_1.mstu[25] = 0;
	}
	ludat1_1.mstj[11] = hparnt_1.ihpr2[10];
	ludat1_1.parj[20] = hparnt_1.hipr1[1];
	rkp = hparnt_1.hipr1[3] * (hparnt_1.hipr1[2] + 2) / ludat1_1.parj[41] / (ludat1_1.parj[40] + 2);
	d__1 = (double)ludat1_1.parj[1];
	d__2 = (double)(1.f / rkp);
	ludat1_1.parj[1] = pow_dd(&d__1, &d__2);
	ludat1_1.parj[20] *= sqrt(rkp);
	hparnt_1.hipr1[1] = ludat1_1.parj[20];
	if (s_cmp(frame, "LAB", (int)4, (int)3) == 0)
	{
		dd1 = (double)(*efrm);
		dd2 = (double)hparnt_1.hint1[7];
		dd3 = (double)hparnt_1.hint1[8];
		r__1 = hparnt_1.hint1[7];
		r__2 = hparnt_1.hint1[8];
		hparnt_1.hint1[0] = sqrt(r__1 * r__1 + hparnt_1.hint1[8] * 2.f * *efrm + r__2 * r__2);
		d__1 = dd1;
		d__2 = dd2;
		dd4 = sqrt(d__1 * d__1 - d__2 * d__2) / (dd1 + dd3);
		hparnt_1.hint1[1] = (float)dd4;
		hparnt_1.hint1[2] = (float)log((dd4 + 1.) / (1. - dd4)) * .5f;
		d__1 = dd1;
		d__2 = dd2;
		dd4 = sqrt(d__1 * d__1 - d__2 * d__2) / dd1;
		hparnt_1.hint1[3] = (float)log((dd4 + 1.) / (1. - dd4)) * .5f;
		hparnt_1.hint1[4] = 0.f;
		hparnt_1.hint1[5] = *efrm;
		hparnt_1.hint1[6] = hparnt_1.hint1[8];
	}
	else if (s_cmp(frame, "CMS", (int)4, (int)3) == 0)
	{
		hparnt_1.hint1[0] = *efrm;
		hparnt_1.hint1[1] = 0.f;
		hparnt_1.hint1[2] = 0.f;
		dd1 = (double)hparnt_1.hint1[0];
		dd2 = (double)hparnt_1.hint1[7];
		dd3 = (double)hparnt_1.hint1[8];
		d__1 = dd2;
		d__2 = dd1;
		dd4 = sqrt(1. - d__1 * d__1 * 4. / (d__2 * d__2));
		hparnt_1.hint1[3] = (float)log((dd4 + 1.) / (1. - dd4)) * .5f;
		d__1 = dd3;
		d__2 = dd1;
		dd4 = sqrt(1. - d__1 * d__1 * 4. / (d__2 * d__2));
		hparnt_1.hint1[4] = (float)log((dd4 + 1.) / (1. - dd4)) * -.5f;
		hparnt_1.hint1[5] = hparnt_1.hint1[0] / 2.f;
		hparnt_1.hint1[6] = hparnt_1.hint1[0] / 2.f;
	}
	if (hparnt_1.ihnt2[0] > 1)
	{
		hijwds_(hparnt_1.ihnt2, &c__1, &rmax);
		hparnt_1.hipr1[33] = rmax;
	}
	if (hparnt_1.ihnt2[2] > 1)
	{
		hijwds_(&hparnt_1.ihnt2[2], &c__2, &rmax);
		hparnt_1.hipr1[34] = rmax;
	}

	i__ = 0;
L20:
	++i__;
	if (i__ == 10)
	{
		goto L30;
	}
	if (hijdat_1.hidat0[i__ * 10 - 1] <= hparnt_1.hint1[0])
	{
		goto L20;
	}
L30:
	if (i__ == 1)
	{
		i__ = 2;
	}
	for (j = 1; j <= 9; ++j)
	{
		hijdat_1.hidat[j - 1] = hijdat_1.hidat0[j + (i__ - 1) * 10 - 11] + (hijdat_1.hidat0[j + i__ * 10 - 11] - hijdat_1.hidat0[j + (i__ - 1) * 10 - 11]) * (hparnt_1.hint1[0] - hijdat_1.hidat0[(i__ - 1) * 10 - 1]) / (hijdat_1.hidat0[i__ * 10 - 1] - hijdat_1.hidat0[(i__ - 1) * 10 - 1]);
	}
	hparnt_1.hipr1[30] = hijdat_1.hidat[4];
	hparnt_1.hipr1[29] = hijdat_1.hidat[4] * 2.f;

	hijcrs_();

	if (hparnt_1.ihpr2[4] != 0)
	{
		hifun_(&c__3, &c_b26, &c_b296, (int *)fnkick_);
	}
	hifun_(&c__7, &c_b26, &c_b299, (int *)fnkc2_);
	hifun_(&c__4, &c_b26, &c_b302, (int *)fnstru_);
	hifun_(&c__5, &c_b26, &c_b302, (int *)fnstrm_);
	hifun_(&c__6, &c_b26, &c_b302, (int *)fnstrs_);
	s_copy(eframe, "Ecm", (int)4, (int)3);
	if (s_cmp(frame, "LAB", (int)4, (int)3) == 0)
	{
		s_copy(eframe, "Elab", (int)4, (int)4);
	}
	s_wsfe(&io___116);
	do_fio(&c__1, eframe, (int)4);
	do_fio(&c__1, (char *)&(*efrm), (int)sizeof(float));
	do_fio(&c__1, proj, (int)4);
	do_fio(&c__1, (char *)&hparnt_1.ihnt2[0], (int)sizeof(int));
	do_fio(&c__1, (char *)&hparnt_1.ihnt2[1], (int)sizeof(int));
	do_fio(&c__1, targ, (int)4);
	do_fio(&c__1, (char *)&hparnt_1.ihnt2[2], (int)sizeof(int));
	do_fio(&c__1, (char *)&hparnt_1.ihnt2[3], (int)sizeof(int));
	e_wsfe();
	return 0;
}

int fnkick_(float *x)
{
	float ret_val, r__1, r__2;
	double sqrt(double), exp(double);
	r__1 = hparnt_1.hipr1[18];
	r__2 = hparnt_1.hipr1[19];
	ret_val = 1.f / (*x + r__1 * r__1) / (*x + r__2 * r__2) / (exp((sqrt(*x) - hparnt_1.hipr1[19]) / .4f) + 1);
	return ret_val;
}

int fnkc2_(float *x)
{
	float ret_val;
	double exp(double);
	ret_val = *x * exp(*x * -2.f / hparnt_1.hipr1[41]);
	return ret_val;
}

int fnstru_(float *x)
{
	float ret_val, r__1, r__2, r__3;
	double d__1, d__2, d__3, d__4;
	double pow_dd(double *, double *);
	d__1 = (double)(1.f - *x);
	d__2 = (double)hparnt_1.hipr1[43];
	r__1 = *x;
	r__2 = hparnt_1.hipr1[44];
	r__3 = hparnt_1.hint1[0];
	d__3 = (double)(r__1 * r__1 + r__2 * r__2 / (r__3 * r__3));
	d__4 = (double)hparnt_1.hipr1[45];
	ret_val = pow_dd(&d__1, &d__2) / pow_dd(&d__3, &d__4);
	return ret_val;
}

int fnstrm_(float *x)
{
	float ret_val, r__1, r__2, r__3, r__4, r__5, r__6;
	double d__1, d__2, d__3, d__4;
	double pow_dd(double *, double *);
	r__1 = 1.f - *x;
	r__2 = hparnt_1.hipr1[44];
	r__3 = hparnt_1.hint1[0];
	d__1 = (double)(r__1 * r__1 + r__2 * r__2 / (r__3 * r__3));
	d__2 = (double)hparnt_1.hipr1[45];
	r__4 = *x;
	r__5 = hparnt_1.hipr1[44];
	r__6 = hparnt_1.hint1[0];
	d__3 = (double)(r__4 * r__4 + r__5 * r__5 / (r__6 * r__6));
	d__4 = (double)hparnt_1.hipr1[45];
	ret_val = 1.f / pow_dd(&d__1, &d__2) / pow_dd(&d__3, &d__4);
	return ret_val;
}

int fnstrs_(float *x)
{
	float ret_val, r__1, r__2, r__3;
	double d__1, d__2, d__3, d__4;
	double pow_dd(double *, double *);
	d__1 = (double)(1.f - *x);
	d__2 = (double)hparnt_1.hipr1[46];
	r__1 = *x;
	r__2 = hparnt_1.hipr1[44];
	r__3 = hparnt_1.hint1[0];
	d__3 = (double)(r__1 * r__1 + r__2 * r__2 / (r__3 * r__3));
	d__4 = (double)hparnt_1.hipr1[47];
	ret_val = pow_dd(&d__1, &d__2) / pow_dd(&d__3, &d__4);
	return ret_val;
}

int hboost_(void)
{
	int i__1;
	float r__1, r__2, r__3;
	double d__1;

	int s_wsle(cilist *), do_lio(int *, int *, char *, int),e_wsle(void);
	double sqrt(double), log(double), sinh(double), cosh(double);

	static int i__;
	static float y;
	static double db, dp3, dp4, dga;
	static float amt;
	static double dbeta;

	static cilist io___120 = {0, 6, 0, 0, 0};

	i__1 = lujets_1.n;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		dbeta = (double)(lujets_1.p[i__ + 17999] / lujets_1.p[i__ +26999]);
		if (abs(dbeta) >= 1.)
		{
			db = (double)hparnt_1.hint1[1];
			if (db > .99999999)
			{
				s_wsle(&io___120);
				do_lio(&c__9, &c__1, "(HIBOOT:) boost vector too large", (int)32);
				e_wsle();
				db = .99999999;
			}
			d__1 = db;
			dga = 1. / sqrt(1. - d__1 * d__1);
			dp3 = (double)lujets_1.p[i__ + 17999];
			dp4 = (double)lujets_1.p[i__ + 26999];
			lujets_1.p[i__ + 17999] = (float)((dp3 + db * dp4) * dga);
			lujets_1.p[i__ + 26999] = (float)((dp4 + db * dp3) * dga);
		}
		y = (float)log((dbeta + 1.) / (1. - dbeta)) * .5f;
		r__1 = lujets_1.p[i__ - 1];
		r__2 = lujets_1.p[i__ + 8999];
		r__3 = lujets_1.p[i__ + 35999];
		amt = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
		lujets_1.p[i__ + 17999] = amt * sinh(y + hparnt_1.hint1[2]);
		lujets_1.p[i__ + 26999] = amt * cosh(y + hparnt_1.hint1[2]);
	}
	return 0;
}

int quench_(int *jpjt, int *ntp)
{
	int i__1, i__2, i__3;
	float r__1, r__2, r__3, r__4;

	double cos(double), sin(double), sqrt(double), exp(double);

	static int i__, i2, j2;
	static float r0, v1, v2, v3, bb, de, dp, rd;
	static int jp, kp;
	static float dx, dy;
	static int lq, kt, mp, mt, nq;
	static float rn;
	static int jt;
	static float xj, yj, dp1, rd0, dp2, dp3, bbx, bby, phi;
	static int isg;
	static float rdp[300], rdt[300];
	static int lqp[300];
	static float drr;
	static int lqt[300];
	static float dphi, phip, phit, phiq, ptot, amshu, ershu, prshu, ptjet0;
	extern double ulangl_(float *, float *), ranart_(int *);

	bb = hparnt_1.hint1[18];
	phi = hparnt_1.hint1[19];
	bbx = bb * cos(phi);
	bby = bb * sin(phi);

	if (*ntp == 2)
	{
		goto L400;
	}
	if (*ntp == 3)
	{
		goto L2000;
	}
	if (hstrng_1.nfp[*jpjt + 1799] != 1)
	{
		return 0;
	}
	jp = *jpjt;
	i__1 = hjjet1_1.npj[jp - 1];
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		r__1 = hjjet1_1.pjpx[jp + i__ * 300 - 301];
		r__2 = hjjet1_1.pjpy[jp + i__ * 300 - 301];
		ptjet0 = sqrt(r__1 * r__1 + r__2 * r__2);
		if (ptjet0 <= hparnt_1.hipr1[10])
		{
			goto L290;
		}
		r__1 = hjjet1_1.pjpz[jp + i__ * 300 - 301];
		ptot = sqrt(ptjet0 * ptjet0 + r__1 * r__1);
		if (ptot < hparnt_1.hipr1[7])
		{
			goto L290;
		}
		phip = ulangl_(&hjjet1_1.pjpx[jp + i__ * 300 - 301], &hjjet1_1.pjpy[jp + i__ * 300 - 301]);
		kp = 0;
		i__2 = hparnt_1.ihnt2[0];
		for (i2 = 1; i2 <= i__2; ++i2)
		{
			dx = hjcrdn_1.yp[i2 * 3 - 3] - hjcrdn_1.yp[jp * 3 - 3];
			dy = hjcrdn_1.yp[i2 * 3 - 2] - hjcrdn_1.yp[jp * 3 - 2];
			phi = ulangl_(&dx, &dy);
			dphi = (r__1 = phi - phip, dabs(r__1));
			if (dphi >= hparnt_1.hipr1[39])
			{
				dphi = hparnt_1.hipr1[39] * 2.f - dphi;
			}
			rd0 = sqrt(dx * dx + dy * dy);
			++kp;
			lqp[kp - 1] = i2;
			rdp[kp - 1] = cos(dphi) * rd0;
		}
		i__2 = kp - 1;
		for (i2 = 1; i2 <= i__2; ++i2)
		{
			i__3 = kp;
			for (j2 = i2 + 1; j2 <= i__3; ++j2)
			{
				rd = rdp[i2 - 1];
				lq = lqp[i2 - 1];
				rdp[i2 - 1] = rdp[j2 - 1];
				lqp[i2 - 1] = lqp[j2 - 1];
				rdp[j2 - 1] = rd;
				lqp[j2 - 1] = lq;
			}
		}
		kt = 0;
		i__3 = hparnt_1.ihnt2[2];
		for (i2 = 1; i2 <= i__3; ++i2)
		{
			dx = hjcrdn_1.yt[i2 * 3 - 3] - hjcrdn_1.yp[jp * 3 - 3] - bbx;
			dy = hjcrdn_1.yt[i2 * 3 - 2] - hjcrdn_1.yp[jp * 3 - 2] - bby;
			phi = ulangl_(&dx, &dy);
			dphi = (r__1 = phi - phip, dabs(r__1));
			if (dphi >= hparnt_1.hipr1[39])
			{
				dphi = hparnt_1.hipr1[39] * 2.f - dphi;
			}
			rd0 = sqrt(dx * dx + dy * dy);
			++kt;
			lqt[kt - 1] = i2;
			rdt[kt - 1] = cos(dphi) * rd0;
		}
		i__3 = kt - 1;
		for (i2 = 1; i2 <= i__3; ++i2)
		{
			i__2 = kt;
			for (j2 = i2 + 1; j2 <= i__2; ++j2)
			{
				rd = rdt[i2 - 1];
				lq = lqt[i2 - 1];
				rdt[i2 - 1] = rdt[j2 - 1];
				lqt[i2 - 1] = lqt[j2 - 1];
				rdt[j2 - 1] = rd;
				lqt[j2 - 1] = lq;
			}
		}
		mp = 0;
		mt = 0;
		r0 = 0.f;
		nq = 0;
		dp = 0.f;
		r__1 = hjjet1_1.pjpx[jp + i__ * 300 - 301];
		r__2 = hjjet1_1.pjpy[jp + i__ * 300 - 301];
		r__3 = hjjet1_1.pjpz[jp + i__ * 300 - 301];
		ptot = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
		v1 = hjjet1_1.pjpx[jp + i__ * 300 - 301] / ptot;
		v2 = hjjet1_1.pjpy[jp + i__ * 300 - 301] / ptot;
		v3 = hjjet1_1.pjpz[jp + i__ * 300 - 301] / ptot;
	L200:
		rn = ranart_(&rndf77_1.nseed);
	L210:
		if (mt >= kt)
		{
			goto L220;
		}
		if (mp >= kp)
		{
			goto L240;
		}
		if (rdp[mp] > rdt[mt])
		{
			goto L240;
		}
	L220:
		++mp;
		drr = rdp[mp - 1] - r0;
		if (rn >= 1.f - exp(-drr / hparnt_1.hipr1[12]))
		{
			goto L210;
		}
		dp = drr * hparnt_1.hipr1[13];
		if (hjjet1_1.kfpj[jp + i__ * 300 - 301] != 21)
		{
			dp *= .5f;
		}
		if (dp <= .2f)
		{
			goto L210;
		}
		if (ptot <= .4f)
		{
			goto L290;
		}
		if (ptot <= dp)
		{
			dp = ptot - .2f;
		}
		de = dp;
		if (hjjet1_1.kfpj[jp + i__ * 300 - 301] != 21)
		{
			r__1 = hstrng_1.pp[lqp[mp - 1] - 1];
			r__2 = hstrng_1.pp[lqp[mp - 1] + 299];
			r__3 = hstrng_1.pp[lqp[mp - 1] + 599];
			prshu = r__1 * r__1 + r__2 * r__2 + r__3 * r__3;
			r__1 = hjjet1_1.pjpm[jp + i__ * 300 - 301];
			r__2 = ptot;
			r__3 = hjjet1_1.pjpm[jp + i__ * 300 - 301];
			r__4 = ptot - dp;
			de = sqrt(r__1 * r__1 + r__2 * r__2) - sqrt(r__3 * r__3 + r__4 *
																		  r__4);
			r__1 = hstrng_1.pp[lqp[mp - 1] + 899] + de - dp;
			ershu = r__1 * r__1;
			amshu = ershu - prshu;
			if (amshu < hparnt_1.hipr1[0] * hparnt_1.hipr1[0])
			{
				goto L210;
			}
			hstrng_1.pp[lqp[mp - 1] + 899] = sqrt(ershu);
			hstrng_1.pp[lqp[mp - 1] + 1199] = sqrt(amshu);
		}
		r0 = rdp[mp - 1];
		dp1 = dp * v1;
		dp2 = dp * v2;
		dp3 = dp * v3;
		++hjjet1_1.npj[lqp[mp - 1] - 1];
		hjjet1_1.kfpj[lqp[mp - 1] + hjjet1_1.npj[lqp[mp - 1] - 1] * 300 - 301] = 21;
		hjjet1_1.pjpx[lqp[mp - 1] + hjjet1_1.npj[lqp[mp - 1] - 1] * 300 - 301] = dp1;
		hjjet1_1.pjpy[lqp[mp - 1] + hjjet1_1.npj[lqp[mp - 1] - 1] * 300 - 301] = dp2;
		hjjet1_1.pjpz[lqp[mp - 1] + hjjet1_1.npj[lqp[mp - 1] - 1] * 300 - 301] = dp3;
		hjjet1_1.pjpe[lqp[mp - 1] + hjjet1_1.npj[lqp[mp - 1] - 1] * 300 - 301] = dp;
		hjjet1_1.pjpm[lqp[mp - 1] + hjjet1_1.npj[lqp[mp - 1] - 1] * 300 - 301] = 0.f;
		goto L260;
	L240:
		++mt;
		drr = rdt[mt - 1] - r0;
		if (rn >= 1.f - exp(-drr / hparnt_1.hipr1[12]))
		{
			goto L210;
		}
		dp = drr * hparnt_1.hipr1[13];
		if (dp <= .2f)
		{
			goto L210;
		}
		if (ptot <= dp)
		{
			dp = ptot - .2f;
		}
		de = dp;
		if (hjjet1_1.kfpj[jp + i__ * 300 - 301] != 21)
		{
			r__1 = hstrng_1.pt[lqt[mt - 1] - 1];
			r__2 = hstrng_1.pt[lqt[mt - 1] + 299];
			r__3 = hstrng_1.pt[lqt[mt - 1] + 599];
			prshu = r__1 * r__1 + r__2 * r__2 + r__3 * r__3;
			r__1 = hjjet1_1.pjpm[jp + i__ * 300 - 301];
			r__2 = ptot;
			r__3 = hjjet1_1.pjpm[jp + i__ * 300 - 301];
			r__4 = ptot - dp;
			de = sqrt(r__1 * r__1 + r__2 * r__2) - sqrt(r__3 * r__3 + r__4 *
																		  r__4);
			r__1 = hstrng_1.pt[lqt[mt - 1] + 899] + de - dp;
			ershu = r__1 * r__1;
			amshu = ershu - prshu;
			if (amshu < hparnt_1.hipr1[0] * hparnt_1.hipr1[0])
			{
				goto L210;
			}
			hstrng_1.pt[lqt[mt - 1] + 899] = sqrt(ershu);
			hstrng_1.pt[lqt[mt - 1] + 1199] = sqrt(amshu);
		}
		r0 = rdt[mt - 1];
		dp1 = dp * v1;
		dp2 = dp * v2;
		dp3 = dp * v3;
		++hjjet1_1.ntj[lqt[mt - 1] - 1];
		hjjet1_1.kftj[lqt[mt - 1] + hjjet1_1.ntj[lqt[mt - 1] - 1] * 300 - 301] = 21;
		hjjet1_1.pjtx[lqt[mt - 1] + hjjet1_1.ntj[lqt[mt - 1] - 1] * 300 - 301] = dp1;
		hjjet1_1.pjty[lqt[mt - 1] + hjjet1_1.ntj[lqt[mt - 1] - 1] * 300 - 301] = dp2;
		hjjet1_1.pjtz[lqt[mt - 1] + hjjet1_1.ntj[lqt[mt - 1] - 1] * 300 - 301] = dp3;
		hjjet1_1.pjte[lqt[mt - 1] + hjjet1_1.ntj[lqt[mt - 1] - 1] * 300 - 301] = dp;
		hjjet1_1.pjtm[lqt[mt - 1] + hjjet1_1.ntj[lqt[mt - 1] - 1] * 300 - 301] = 0.f;
	L260:
		hjjet1_1.pjpx[jp + i__ * 300 - 301] = (ptot - dp) * v1;
		hjjet1_1.pjpy[jp + i__ * 300 - 301] = (ptot - dp) * v2;
		hjjet1_1.pjpz[jp + i__ * 300 - 301] = (ptot - dp) * v3;
		hjjet1_1.pjpe[jp + i__ * 300 - 301] -= de;
		ptot -= dp;
		++nq;
		goto L200;
	}
	return 0;
L400:
	if (hstrng_1.nft[*jpjt + 1799] != 1)
	{
		return 0;
	}
	jt = *jpjt;
	i__1 = hjjet1_1.ntj[jt - 1];
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		r__1 = hjjet1_1.pjtx[jt + i__ * 300 - 301];
		r__2 = hjjet1_1.pjty[jt + i__ * 300 - 301];
		ptjet0 = sqrt(r__1 * r__1 + r__2 * r__2);
		if (ptjet0 <= hparnt_1.hipr1[10])
		{
			goto L690;
		}
		r__1 = hjjet1_1.pjtz[jt + i__ * 300 - 301];
		ptot = sqrt(ptjet0 * ptjet0 + r__1 * r__1);
		if (ptot < hparnt_1.hipr1[7])
		{
			goto L690;
		}
		phit = ulangl_(&hjjet1_1.pjtx[jt + i__ * 300 - 301], &hjjet1_1.pjty[jt + i__ * 300 - 301]);
		kp = 0;
		i__2 = hparnt_1.ihnt2[0];
		for (i2 = 1; i2 <= i__2; ++i2)
		{
			dx = hjcrdn_1.yp[i2 * 3 - 3] + bbx - hjcrdn_1.yt[jt * 3 - 3];
			dy = hjcrdn_1.yp[i2 * 3 - 2] + bby - hjcrdn_1.yt[jt * 3 - 2];
			phi = ulangl_(&dx, &dy);
			dphi = (r__1 = phi - phit, dabs(r__1));
			if (dphi >= hparnt_1.hipr1[39])
			{
				dphi = hparnt_1.hipr1[39] * 2.f - dphi;
			}
			rd0 = sqrt(dx * dx + dy * dy);
			++kp;
			lqp[kp - 1] = i2;
			rdp[kp - 1] = cos(dphi) * rd0;
		}
		i__2 = kp - 1;
		for (i2 = 1; i2 <= i__2; ++i2)
		{
			i__3 = kp;
			for (j2 = i2 + 1; j2 <= i__3; ++j2)
			{
				if (rdp[i2 - 1] < rdp[j2 - 1])
				{
					goto L510;
				}
				rd = rdp[i2 - 1];
				lq = lqp[i2 - 1];
				rdp[i2 - 1] = rdp[j2 - 1];
				lqp[i2 - 1] = lqp[j2 - 1];
				rdp[j2 - 1] = rd;
				lqp[j2 - 1] = lq;
			L510:;
			}
		}
		kt = 0;
		i__3 = hparnt_1.ihnt2[2];
		for (i2 = 1; i2 <= i__3; ++i2)
		{
			dx = hjcrdn_1.yt[i2 * 3 - 3] - hjcrdn_1.yt[jt * 3 - 3];
			dy = hjcrdn_1.yt[i2 * 3 - 2] - hjcrdn_1.yt[jt * 3 - 2];
			phi = ulangl_(&dx, &dy);
			dphi = (r__1 = phi - phit, dabs(r__1));
			if (dphi >= hparnt_1.hipr1[39])
			{
				dphi = hparnt_1.hipr1[39] * 2.f - dphi;
			}
			rd0 = sqrt(dx * dx + dy * dy);
			++kt;
			lqt[kt - 1] = i2;
			rdt[kt - 1] = cos(dphi) * rd0;
		}
		i__3 = kt - 1;
		for (i2 = 1; i2 <= i__3; ++i2)
		{
			i__2 = kt;
			for (j2 = i2 + 1; j2 <= i__2; ++j2)
			{
				rd = rdt[i2 - 1];
				lq = lqt[i2 - 1];
				rdt[i2 - 1] = rdt[j2 - 1];
				lqt[i2 - 1] = lqt[j2 - 1];
				rdt[j2 - 1] = rd;
				lqt[j2 - 1] = lq;
			}
		}
		mp = 0;
		mt = 0;
		nq = 0;
		dp = 0.f;
		r0 = 0.f;
		r__1 = hjjet1_1.pjtx[jt + i__ * 300 - 301];
		r__2 = hjjet1_1.pjty[jt + i__ * 300 - 301];
		r__3 = hjjet1_1.pjtz[jt + i__ * 300 - 301];
		ptot = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
		v1 = hjjet1_1.pjtx[jt + i__ * 300 - 301] / ptot;
		v2 = hjjet1_1.pjty[jt + i__ * 300 - 301] / ptot;
		v3 = hjjet1_1.pjtz[jt + i__ * 300 - 301] / ptot;
	L600:
		rn = ranart_(&rndf77_1.nseed);
	L610:
		if (mt >= kt)
		{
			goto L620;
		}
		if (mp >= kp)
		{
			goto L640;
		}
		if (rdp[mp] > rdt[mt])
		{
			goto L640;
		}
	L620:
		++mp;
		drr = rdp[mp - 1] - r0;
		if (rn >= 1.f - exp(-drr / hparnt_1.hipr1[12]))
		{
			goto L610;
		}
		dp = drr * hparnt_1.hipr1[13];
		if (hjjet1_1.kftj[jt + i__ * 300 - 301] != 21)
		{
			dp *= .5f;
		}
		if (dp <= .2f)
		{
			goto L610;
		}
		if (ptot <= dp)
		{
			dp = ptot - .2f;
		}
		de = dp;

		if (hjjet1_1.kftj[jt + i__ * 300 - 301] != 21)
		{
			r__1 = hstrng_1.pp[lqp[mp - 1] - 1];
			r__2 = hstrng_1.pp[lqp[mp - 1] + 299];
			r__3 = hstrng_1.pp[lqp[mp - 1] + 599];
			prshu = r__1 * r__1 + r__2 * r__2 + r__3 * r__3;
			r__1 = hjjet1_1.pjtm[jt + i__ * 300 - 301];
			r__2 = ptot;
			r__3 = hjjet1_1.pjtm[jt + i__ * 300 - 301];
			r__4 = ptot - dp;
			de = sqrt(r__1 * r__1 + r__2 * r__2) - sqrt(r__3 * r__3 + r__4 *
																		  r__4);
			r__1 = hstrng_1.pp[lqp[mp - 1] + 899] + de - dp;
			ershu = r__1 * r__1;
			amshu = ershu - prshu;
			if (amshu < hparnt_1.hipr1[0] * hparnt_1.hipr1[0])
			{
				goto L610;
			}
			hstrng_1.pp[lqp[mp - 1] + 899] = sqrt(ershu);
			hstrng_1.pp[lqp[mp - 1] + 1199] = sqrt(amshu);
		}
		r0 = rdp[mp - 1];
		dp1 = dp * v1;
		dp2 = dp * v2;
		dp3 = dp * v3;
		++hjjet1_1.npj[lqp[mp - 1] - 1];
		hjjet1_1.kfpj[lqp[mp - 1] + hjjet1_1.npj[lqp[mp - 1] - 1] * 300 - 301] = 21;
		hjjet1_1.pjpx[lqp[mp - 1] + hjjet1_1.npj[lqp[mp - 1] - 1] * 300 - 301] = dp1;
		hjjet1_1.pjpy[lqp[mp - 1] + hjjet1_1.npj[lqp[mp - 1] - 1] * 300 - 301] = dp2;
		hjjet1_1.pjpz[lqp[mp - 1] + hjjet1_1.npj[lqp[mp - 1] - 1] * 300 - 301] = dp3;
		hjjet1_1.pjpe[lqp[mp - 1] + hjjet1_1.npj[lqp[mp - 1] - 1] * 300 - 301] = dp;
		hjjet1_1.pjpm[lqp[mp - 1] + hjjet1_1.npj[lqp[mp - 1] - 1] * 300 - 301] = 0.f;
		goto L660;
	L640:
		++mt;
		drr = rdt[mt - 1] - r0;
		if (rn >= 1.f - exp(-drr / hparnt_1.hipr1[12]))
		{
			goto L610;
		}
		dp = drr * hparnt_1.hipr1[13];
		if (dp <= .2f)
		{
			goto L610;
		}
		if (ptot <= dp)
		{
			dp = ptot - .2f;
		}
		de = dp;
		if (hjjet1_1.kftj[jt + i__ * 300 - 301] != 21)
		{
			r__1 = hstrng_1.pt[lqt[mt - 1] - 1];
			r__2 = hstrng_1.pt[lqt[mt - 1] + 299];
			r__3 = hstrng_1.pt[lqt[mt - 1] + 599];
			prshu = r__1 * r__1 + r__2 * r__2 + r__3 * r__3;
			r__1 = hjjet1_1.pjtm[jt + i__ * 300 - 301];
			r__2 = ptot;
			r__3 = hjjet1_1.pjtm[jt + i__ * 300 - 301];
			r__4 = ptot - dp;
			de = sqrt(r__1 * r__1 + r__2 * r__2) - sqrt(r__3 * r__3 + r__4 *
																		  r__4);
			r__1 = hstrng_1.pt[lqt[mt - 1] + 899] + de - dp;
			ershu = r__1 * r__1;
			amshu = ershu - prshu;
			if (amshu < hparnt_1.hipr1[0] * hparnt_1.hipr1[0])
			{
				goto L610;
			}
			hstrng_1.pt[lqt[mt - 1] + 899] = sqrt(ershu);
			hstrng_1.pt[lqt[mt - 1] + 1199] = sqrt(amshu);
		}
		r0 = rdt[mt - 1];
		dp1 = dp * v1;
		dp2 = dp * v2;
		dp3 = dp * v3;
		++hjjet1_1.ntj[lqt[mt - 1] - 1];
		hjjet1_1.kftj[lqt[mt - 1] + hjjet1_1.ntj[lqt[mt - 1] - 1] * 300 - 301] = 21;
		hjjet1_1.pjtx[lqt[mt - 1] + hjjet1_1.ntj[lqt[mt - 1] - 1] * 300 - 301] = dp1;
		hjjet1_1.pjty[lqt[mt - 1] + hjjet1_1.ntj[lqt[mt - 1] - 1] * 300 - 301] = dp2;
		hjjet1_1.pjtz[lqt[mt - 1] + hjjet1_1.ntj[lqt[mt - 1] - 1] * 300 - 301] = dp3;
		hjjet1_1.pjte[lqt[mt - 1] + hjjet1_1.ntj[lqt[mt - 1] - 1] * 300 - 301] = dp;
		hjjet1_1.pjtm[lqt[mt - 1] + hjjet1_1.ntj[lqt[mt - 1] - 1] * 300 - 301] = 0.f;
	L660:
		hjjet1_1.pjtx[jt + i__ * 300 - 301] = (ptot - dp) * v1;
		hjjet1_1.pjty[jt + i__ * 300 - 301] = (ptot - dp) * v2;
		hjjet1_1.pjtz[jt + i__ * 300 - 301] = (ptot - dp) * v3;
		hjjet1_1.pjte[jt + i__ * 300 - 301] -= de;
		ptot -= dp;
		++nq;
		goto L600;
	}
	return 0;
L2000:
	isg = *jpjt;
	if (hjjet2_1.iasg[isg + 300001] != 1)
	{
		return 0;
	}

	jp = hjjet2_1.iasg[isg - 1];
	jt = hjjet2_1.iasg[isg + 150000];
	xj = (hjcrdn_1.yp[jp * 3 - 3] + bbx + hjcrdn_1.yt[jt * 3 - 3]) / 2.f;
	yj = (hjcrdn_1.yp[jp * 3 - 2] + bby + hjcrdn_1.yt[jt * 3 - 2]) / 2.f;
	i__1 = hjjet2_1.njsg[isg - 1];
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		r__1 = hjjet2_1.pxsg[isg + i__ * 150001 - 150002];
		r__2 = hjjet2_1.pysg[isg + i__ * 150001 - 150002];
		ptjet0 = sqrt(r__1 * r__1 + r__2 * r__2);
		if (ptjet0 <= hparnt_1.hipr1[10] || hjjet2_1.pesg[isg + i__ * 150001 - 150002] < hparnt_1.hipr1[0])
		{
			goto L2690;
		}
		r__1 = hjjet2_1.pzsg[isg + i__ * 150001 - 150002];
		ptot = sqrt(ptjet0 * ptjet0 + r__1 * r__1);
		if (ptot < max(hparnt_1.hipr1[0], hparnt_1.hipr1[7]))
		{
			goto L2690;
		}
		phiq = ulangl_(&hjjet2_1.pxsg[isg + i__ * 150001 - 150002], &hjjet2_1.pysg[isg + i__ * 150001 - 150002]);
		kp = 0;
		i__2 = hparnt_1.ihnt2[0];
		for (i2 = 1; i2 <= i__2; ++i2)
		{
			dx = hjcrdn_1.yp[i2 * 3 - 3] + bbx - xj;
			dy = hjcrdn_1.yp[i2 * 3 - 2] + bby - yj;
			phi = ulangl_(&dx, &dy);
			dphi = (r__1 = phi - phiq, dabs(r__1));
			if (dphi >= hparnt_1.hipr1[39])
			{
				dphi = hparnt_1.hipr1[39] * 2.f - dphi;
			}
			rd0 = sqrt(dx * dx + dy * dy);
			++kp;
			lqp[kp - 1] = i2;
			rdp[kp - 1] = cos(dphi) * rd0;
		}
		i__2 = kp - 1;
		for (i2 = 1; i2 <= i__2; ++i2)
		{
			i__3 = kp;
			for (j2 = i2 + 1; j2 <= i__3; ++j2)
			{
				rd = rdp[i2 - 1];
				lq = lqp[i2 - 1];
				rdp[i2 - 1] = rdp[j2 - 1];
				lqp[i2 - 1] = lqp[j2 - 1];
				rdp[j2 - 1] = rd;
				lqp[j2 - 1] = lq;
			}
		}
		kt = 0;
		i__3 = hparnt_1.ihnt2[2];
		for (i2 = 1; i2 <= i__3; ++i2)
		{
			if (hstrng_1.nft[i2 + 1199] != 3 || i2 == jt)
			{
				goto L2520;
			}
			dx = hjcrdn_1.yt[i2 * 3 - 3] - xj;
			dy = hjcrdn_1.yt[i2 * 3 - 2] - yj;
			phi = ulangl_(&dx, &dy);
			dphi = (r__1 = phi - phiq, dabs(r__1));
			if (dphi >= hparnt_1.hipr1[39])
			{
				dphi = hparnt_1.hipr1[39] * 2.f - dphi;
			}
			if (dphi > hparnt_1.hipr1[39] / 2.f)
			{
				goto L2520;
			}
			rd0 = sqrt(dx * dx + dy * dy);
			if (rd0 * sin(dphi) > hparnt_1.hipr1[11])
			{
				goto L2520;
			}
			++kt;
			lqt[kt - 1] = i2;
			rdt[kt - 1] = cos(dphi) * rd0;
		L2520:;
		}
		i__3 = kt - 1;
		for (i2 = 1; i2 <= i__3; ++i2)
		{
			i__2 = kt;
			for (j2 = i2 + 1; j2 <= i__2; ++j2)
			{
				if (rdt[i2 - 1] < rdt[j2 - 1])
				{
					goto L2530;
				}
				rd = rdt[i2 - 1];
				lq = lqt[i2 - 1];
				rdt[i2 - 1] = rdt[j2 - 1];
				lqt[i2 - 1] = lqt[j2 - 1];
				rdt[j2 - 1] = rd;
				lqt[j2 - 1] = lq;
			L2530:;
			}
		}
		mp = 0;
		mt = 0;
		nq = 0;
		dp = 0.f;
		r0 = 0.f;
		r__1 = hjjet2_1.pxsg[isg + i__ * 150001 - 150002];
		r__2 = hjjet2_1.pysg[isg + i__ * 150001 - 150002];
		r__3 = hjjet2_1.pzsg[isg + i__ * 150001 - 150002];
		ptot = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
		v1 = hjjet2_1.pxsg[isg + i__ * 150001 - 150002] / ptot;
		v2 = hjjet2_1.pysg[isg + i__ * 150001 - 150002] / ptot;
		v3 = hjjet2_1.pzsg[isg + i__ * 150001 - 150002] / ptot;
	L2600:
		rn = ranart_(&rndf77_1.nseed);
	L2610:
		if (mt >= kt)
		{
			goto L2620;
		}
		if (mp >= kp)
		{
			goto L2640;
		}
		if (rdp[mp] > rdt[mt])
		{
			goto L2640;
		}
	L2620:
		++mp;
		drr = rdp[mp - 1] - r0;
		if (rn >= 1.f - exp(-drr / hparnt_1.hipr1[12]))
		{
			goto L2610;
		}
		dp = drr * hparnt_1.hipr1[13] / 2.f;
		if (dp <= .2f)
		{
			goto L2610;
		}
		if (ptot <= dp)
		{
			dp = ptot - .2f;
		}
		de = dp;

		if (hjjet2_1.k2sg[isg + i__ * 150001 - 150002] != 21)
		{
			r__1 = hstrng_1.pp[lqp[mp - 1] - 1];
			r__2 = hstrng_1.pp[lqp[mp - 1] + 299];
			r__3 = hstrng_1.pp[lqp[mp - 1] + 599];
			prshu = r__1 * r__1 + r__2 * r__2 + r__3 * r__3;
			r__1 = hjjet2_1.pmsg[isg + i__ * 150001 - 150002];
			r__2 = ptot;
			r__3 = hjjet2_1.pmsg[isg + i__ * 150001 - 150002];
			r__4 = ptot - dp;
			de = sqrt(r__1 * r__1 + r__2 * r__2) - sqrt(r__3 * r__3 + r__4 * r__4);
			r__1 = hstrng_1.pp[lqp[mp - 1] + 899] + de - dp;
			ershu = r__1 * r__1;
			amshu = ershu - prshu;
			if (amshu < hparnt_1.hipr1[0] * hparnt_1.hipr1[0])
			{
				goto L2610;
			}
			hstrng_1.pp[lqp[mp - 1] + 899] = sqrt(ershu);
			hstrng_1.pp[lqp[mp - 1] + 1199] = sqrt(amshu);
		}
		r0 = rdp[mp - 1];
		dp1 = dp * v1;
		dp2 = dp * v2;
		dp3 = dp * v3;
		++hjjet1_1.npj[lqp[mp - 1] - 1];
		hjjet1_1.kfpj[lqp[mp - 1] + hjjet1_1.npj[lqp[mp - 1] - 1] * 300 - 301] = 21;
		hjjet1_1.pjpx[lqp[mp - 1] + hjjet1_1.npj[lqp[mp - 1] - 1] * 300 - 301] = dp1;
		hjjet1_1.pjpy[lqp[mp - 1] + hjjet1_1.npj[lqp[mp - 1] - 1] * 300 - 301] = dp2;
		hjjet1_1.pjpz[lqp[mp - 1] + hjjet1_1.npj[lqp[mp - 1] - 1] * 300 - 301] = dp3;
		hjjet1_1.pjpe[lqp[mp - 1] + hjjet1_1.npj[lqp[mp - 1] - 1] * 300 - 301] = dp;
		hjjet1_1.pjpm[lqp[mp - 1] + hjjet1_1.npj[lqp[mp - 1] - 1] * 300 - 301] = 0.f;
		goto L2660;
	L2640:
		++mt;
		drr = rdt[mt - 1] - r0;
		if (rn >= 1.f - exp(-drr / hparnt_1.hipr1[12]))
		{
			goto L2610;
		}
		dp = drr * hparnt_1.hipr1[13];
		if (dp <= .2f)
		{
			goto L2610;
		}
		if (ptot <= dp)
		{
			dp = ptot - .2f;
		}
		de = dp;
		if (hjjet2_1.k2sg[isg + i__ * 150001 - 150002] != 21)
		{
			r__1 = hstrng_1.pt[lqt[mt - 1] - 1];
			r__2 = hstrng_1.pt[lqt[mt - 1] + 299];
			r__3 = hstrng_1.pt[lqt[mt - 1] + 599];
			prshu = r__1 * r__1 + r__2 * r__2 + r__3 * r__3;
			r__1 = hjjet2_1.pmsg[isg + i__ * 150001 - 150002];
			r__2 = ptot;
			r__3 = hjjet2_1.pmsg[isg + i__ * 150001 - 150002];
			r__4 = ptot - dp;
			de = sqrt(r__1 * r__1 + r__2 * r__2) - sqrt(r__3 * r__3 + r__4 *
																		  r__4);
			r__1 = hstrng_1.pt[lqt[mt - 1] + 899] + de - dp;
			ershu = r__1 * r__1;
			amshu = ershu - prshu;
			if (amshu < hparnt_1.hipr1[0] * hparnt_1.hipr1[0])
			{
				goto L2610;
			}
			hstrng_1.pt[lqt[mt - 1] + 899] = sqrt(ershu);
			hstrng_1.pt[lqt[mt - 1] + 1199] = sqrt(amshu);
		}
		r0 = rdt[mt - 1];
		dp1 = dp * v1;
		dp2 = dp * v2;
		dp3 = dp * v3;
		++hjjet1_1.ntj[lqt[mt - 1] - 1];
		hjjet1_1.kftj[lqt[mt - 1] + hjjet1_1.ntj[lqt[mt - 1] - 1] * 300 - 301] = 21;
		hjjet1_1.pjtx[lqt[mt - 1] + hjjet1_1.ntj[lqt[mt - 1] - 1] * 300 - 301] = dp1;
		hjjet1_1.pjty[lqt[mt - 1] + hjjet1_1.ntj[lqt[mt - 1] - 1] * 300 - 301] = dp2;
		hjjet1_1.pjtz[lqt[mt - 1] + hjjet1_1.ntj[lqt[mt - 1] - 1] * 300 - 301] = dp3;
		hjjet1_1.pjte[lqt[mt - 1] + hjjet1_1.ntj[lqt[mt - 1] - 1] * 300 - 301] = dp;
		hjjet1_1.pjtm[lqt[mt - 1] + hjjet1_1.ntj[lqt[mt - 1] - 1] * 300 - 301] = 0.f;
	L2660:
		hjjet2_1.pxsg[isg + i__ * 150001 - 150002] = (ptot - dp) * v1;
		hjjet2_1.pysg[isg + i__ * 150001 - 150002] = (ptot - dp) * v2;
		hjjet2_1.pzsg[isg + i__ * 150001 - 150002] = (ptot - dp) * v3;
		hjjet2_1.pesg[isg + i__ * 150001 - 150002] -= de;
		ptot -= dp;
		++nq;
		goto L2600;
	}
	return 0;
}

int hijfrg_(int *jtp, int *ntp, int *ierror)
{
	int i__1;
	float r__1, r__2, r__3, r__4;

	int s_wsle(cilist *), do_lio(int *, int *, char *, int),
		e_wsle(void);
	double sqrt(double);

	static int i__, j, i0, ii, jj;
	static float am1, am2;
	static int kf1, kf2;
	static float pb1, pb2, pb3;
	static int kk1;
	static float amt;
	static int isg;
	static float pq21, pq22, pq11, pq12;
	static int iex;
	static float btz, amt1, amt2, pecm, pzcm, hpr150, pmax1, pmax2, pmax3,
		hdat20;
	static int jetot;
	extern int attrad_(int *), hirobo_(float *, float *, float *, float *, float *), luexec_(void), luedit_(int *);
	extern double ranart_(int *);

	static cilist io___187 = {0, 6, 0, 0, 0};
	static cilist io___188 = {0, 6, 0, 0, 0};
	static cilist io___189 = {0, 6, 0, 0, 0};
	static cilist io___190 = {0, 6, 0, 0, 0};

	*ierror = 0;
	luedit_(&c__0);
	lujets_1.n = 0;
	if (*ntp == 3)
	{
		isg = *jtp;
		lujets_1.n = hjjet2_1.njsg[isg - 1];
		i__1 = hjjet2_1.njsg[isg - 1];
		for (i__ = 1; i__ <= i__1; ++i__)
		{
			lujets_1.k[i__ - 1] = hjjet2_1.k1sg[isg + i__ * 150001 - 150002];
			lujets_1.k[i__ + 8999] = hjjet2_1.k2sg[isg + i__ * 150001 -
												   150002];
			lujets_1.p[i__ - 1] = hjjet2_1.pxsg[isg + i__ * 150001 - 150002];
			lujets_1.p[i__ + 8999] = hjjet2_1.pysg[isg + i__ * 150001 -
												   150002];
			lujets_1.p[i__ + 17999] = hjjet2_1.pzsg[isg + i__ * 150001 -
													150002];
			lujets_1.p[i__ + 26999] = hjjet2_1.pesg[isg + i__ * 150001 -
													150002];
			lujets_1.p[i__ + 35999] = hjjet2_1.pmsg[isg + i__ * 150001 -
													150002];
		}
		if (anim_1.isoft != 2 || anim_1.isflag != 0)
		{
			luexec_();
		}
		return 0;
	}

	if (*ntp == 2)
	{
		goto L200;
	}
	if (*jtp > hparnt_1.ihnt2[0])
	{
		return 0;
	}
	if (hstrng_1.nfp[*jtp + 1199] != 3 && hstrng_1.nfp[*jtp + 599] != 0 &&
		hjjet1_1.npj[*jtp - 1] == 0 && hstrng_1.nfp[*jtp + 2699] == 0)
	{
		goto L1000;
	}
	if (hstrng_1.nfp[*jtp + 4199] == -1)
	{
		kf1 = hstrng_1.nfp[*jtp + 299];
		kf2 = hstrng_1.nfp[*jtp - 1];
		pq21 = hstrng_1.pp[*jtp + 1499];
		pq22 = hstrng_1.pp[*jtp + 1799];
		pq11 = hstrng_1.pp[*jtp + 2099];
		pq12 = hstrng_1.pp[*jtp + 2399];
		am1 = hstrng_1.pp[*jtp + 4199];
		am2 = hstrng_1.pp[*jtp + 3899];
	}
	else
	{
		kf1 = hstrng_1.nfp[*jtp - 1];
		kf2 = hstrng_1.nfp[*jtp + 299];
		pq21 = hstrng_1.pp[*jtp + 2099];
		pq22 = hstrng_1.pp[*jtp + 2399];
		pq11 = hstrng_1.pp[*jtp + 1499];
		pq12 = hstrng_1.pp[*jtp + 1799];
		am1 = hstrng_1.pp[*jtp + 3899];
		am2 = hstrng_1.pp[*jtp + 4199];
	}
	pb1 = pq11 + pq21;
	pb2 = pq12 + pq22;
	pb3 = hstrng_1.pp[*jtp + 599];
	pecm = hstrng_1.pp[*jtp + 1199];
	btz = pb3 / hstrng_1.pp[*jtp + 899];
	if (((r__1 = pb1 - hstrng_1.pp[*jtp - 1], dabs(r__1)) > .01f || (r__2 =
																		 pb2 - hstrng_1.pp[*jtp + 299],
																	 dabs(r__2)) > .01f) &&
		hparnt_1.ihpr2[9] != 0)
	{
		s_wsle(&io___187);
		do_lio(&c__9, &c__1, "  Pt of Q and QQ do not sum to the total", (int)40);
		do_lio(&c__3, &c__1, (char *)&(*jtp), (int)sizeof(int));
		do_lio(&c__3, &c__1, (char *)&(*ntp), (int)sizeof(int));
		do_lio(&c__4, &c__1, (char *)&pq11, (int)sizeof(float));
		do_lio(&c__4, &c__1, (char *)&pq21, (int)sizeof(float));
		do_lio(&c__4, &c__1, (char *)&pb1, (int)sizeof(float));
		do_lio(&c__9, &c__1, "*", (int)1);
		do_lio(&c__4, &c__1, (char *)&pq12, (int)sizeof(float));
		do_lio(&c__4, &c__1, (char *)&pq22, (int)sizeof(float));
		do_lio(&c__4, &c__1, (char *)&pb2, (int)sizeof(float));
		do_lio(&c__9, &c__1, "*", (int)1);
		do_lio(&c__4, &c__1, (char *)&hstrng_1.pp[*jtp - 1], (int)sizeof(float));
		do_lio(&c__4, &c__1, (char *)&hstrng_1.pp[*jtp + 299], (int)sizeof(float));
		e_wsle();
	}
	goto L300;
L200:
	if (*jtp > hparnt_1.ihnt2[2])
	{
		return 0;
	}
	if (hstrng_1.nft[*jtp + 1199] != 3 && hstrng_1.nft[*jtp + 599] != 0 &&
		hjjet1_1.ntj[*jtp - 1] == 0 && hstrng_1.nft[*jtp + 2699] == 0)
	{
		goto L1200;
	}
	if (hstrng_1.nft[*jtp + 4199] == 1)
	{
		kf1 = hstrng_1.nft[*jtp - 1];
		kf2 = hstrng_1.nft[*jtp + 299];
		pq11 = hstrng_1.pt[*jtp + 1499];
		pq12 = hstrng_1.pt[*jtp + 1799];
		pq21 = hstrng_1.pt[*jtp + 2099];
		pq22 = hstrng_1.pt[*jtp + 2399];
		am1 = hstrng_1.pt[*jtp + 3899];
		am2 = hstrng_1.pt[*jtp + 4199];
	}
	else
	{
		kf1 = hstrng_1.nft[*jtp + 299];
		kf2 = hstrng_1.nft[*jtp - 1];
		pq11 = hstrng_1.pt[*jtp + 2099];
		pq12 = hstrng_1.pt[*jtp + 2399];
		pq21 = hstrng_1.pt[*jtp + 1499];
		pq22 = hstrng_1.pt[*jtp + 1799];
		am1 = hstrng_1.pt[*jtp + 4199];
		am2 = hstrng_1.pt[*jtp + 3899];
	}
	pb1 = pq11 + pq21;
	pb2 = pq12 + pq22;
	pb3 = hstrng_1.pt[*jtp + 599];
	pecm = hstrng_1.pt[*jtp + 1199];
	btz = pb3 / hstrng_1.pt[*jtp + 899];
	if (((r__1 = pb1 - hstrng_1.pt[*jtp - 1], dabs(r__1)) > .01f || (r__2 =
																		 pb2 - hstrng_1.pt[*jtp + 299],
																	 dabs(r__2)) > .01f) &&
		hparnt_1.ihpr2[9] != 0)
	{
		s_wsle(&io___188);
		do_lio(&c__9, &c__1, "  Pt of Q and QQ do not sum to the total", (int)40);
		do_lio(&c__3, &c__1, (char *)&(*jtp), (int)sizeof(int));
		do_lio(&c__3, &c__1, (char *)&(*ntp), (int)sizeof(int));
		do_lio(&c__4, &c__1, (char *)&pq11, (int)sizeof(float));
		do_lio(&c__4, &c__1, (char *)&pq21, (int)sizeof(float));
		do_lio(&c__4, &c__1, (char *)&pb1, (int)sizeof(float));
		do_lio(&c__9, &c__1, "*", (int)1);
		do_lio(&c__4, &c__1, (char *)&pq12, (int)sizeof(float));
		do_lio(&c__4, &c__1, (char *)&pq22, (int)sizeof(float));
		do_lio(&c__4, &c__1, (char *)&pb2, (int)sizeof(float));
		do_lio(&c__9, &c__1, "*", (int)1);
		do_lio(&c__4, &c__1, (char *)&hstrng_1.pt[*jtp - 1], (int)sizeof(float));
		do_lio(&c__4, &c__1, (char *)&hstrng_1.pt[*jtp + 299], (int)sizeof(float));
		e_wsle();
	}
L300:
	if (pecm < hparnt_1.hipr1[0])
	{
		*ierror = 1;
		if (hparnt_1.ihpr2[9] == 0)
		{
			return 0;
		}
		s_wsle(&io___189);
		do_lio(&c__9, &c__1, " ECM=", (int)5);
		do_lio(&c__4, &c__1, (char *)&pecm, (int)sizeof(float));
		do_lio(&c__9, &c__1, " energy of the string is too small", (int)34);
		e_wsle();
		s_wsle(&io___190);
		do_lio(&c__9, &c__1, "JTP,NTP,pq=", (int)11);
		do_lio(&c__3, &c__1, (char *)&(*jtp), (int)sizeof(int));
		do_lio(&c__3, &c__1, (char *)&(*ntp), (int)sizeof(int));
		do_lio(&c__4, &c__1, (char *)&pq11, (int)sizeof(float));
		do_lio(&c__4, &c__1, (char *)&pq12, (int)sizeof(float));
		do_lio(&c__4, &c__1, (char *)&pq21, (int)sizeof(float));
		do_lio(&c__4, &c__1, (char *)&pq22, (int)sizeof(float));
		e_wsle();
		return 0;
	}
	r__1 = pecm;
	r__2 = pb1;
	r__3 = pb2;
	amt = r__1 * r__1 + r__2 * r__2 + r__3 * r__3;
	r__1 = am1;
	r__2 = pq11;
	r__3 = pq12;
	amt1 = r__1 * r__1 + r__2 * r__2 + r__3 * r__3;
	r__1 = am2;
	r__2 = pq21;
	r__3 = pq22;
	amt2 = r__1 * r__1 + r__2 * r__2 + r__3 * r__3;
	r__2 = amt;
	r__3 = amt1;
	r__4 = amt2;
	pzcm = sqrt((r__1 = r__2 * r__2 + r__3 * r__3 + r__4 * r__4 - amt * 2.f * amt1 - amt * 2.f * amt2 - amt1 * 2.f * amt2, dabs(r__1))) / 2.f /
		   sqrt(amt);
	lujets_1.k[0] = 2;
	lujets_1.k[9000] = kf1;
	lujets_1.p[0] = pq11;
	lujets_1.p[9000] = pq12;
	lujets_1.p[18000] = pzcm;
	r__1 = pzcm;
	lujets_1.p[27000] = sqrt(amt1 + r__1 * r__1);
	lujets_1.p[36000] = am1;
	lujets_1.k[1] = 1;
	lujets_1.k[9001] = kf2;
	lujets_1.p[1] = pq21;
	lujets_1.p[9001] = pq22;
	lujets_1.p[18001] = -pzcm;
	r__1 = pzcm;
	lujets_1.p[27001] = sqrt(amt2 + r__1 * r__1);
	lujets_1.p[36001] = am2;
	lujets_1.n = 2;
	hirobo_(&c_b26, &c_b26, &c_b26, &c_b26, &btz);
	jetot = 0;
	r__1 = pq21;
	r__2 = pq22;
	r__3 = pq11;
	r__4 = pq12;
	if (r__1 * r__1 + r__2 * r__2 > r__3 * r__3 + r__4 * r__4)
	{
		pmax1 = lujets_1.p[1];
		pmax2 = lujets_1.p[9001];
		pmax3 = lujets_1.p[18001];
	}
	else
	{
		pmax1 = lujets_1.p[0];
		pmax2 = lujets_1.p[9000];
		pmax3 = lujets_1.p[18000];
	}
	if (*ntp == 1)
	{
		hstrng_1.pp[*jtp + 2699] = pmax1;
		hstrng_1.pp[*jtp + 2999] = pmax2;
		hstrng_1.pp[*jtp + 3299] = pmax3;
	}
	else if (*ntp == 2)
	{
		hstrng_1.pt[*jtp + 2699] = pmax1;
		hstrng_1.pt[*jtp + 2999] = pmax2;
		hstrng_1.pt[*jtp + 3299] = pmax3;
	}
	if (*ntp == 1 && hjjet1_1.npj[*jtp - 1] != 0)
	{
		jetot = hjjet1_1.npj[*jtp - 1];
		iex = 0;
		if (abs(kf1) > 1000 && kf1 < 0 || abs(kf1) < 1000 && kf1 > 0)
		{
			iex = 1;
		}
		for (i__ = lujets_1.n; i__ >= 2; --i__)
		{
			for (j = 1; j <= 5; ++j)
			{
				ii = hjjet1_1.npj[*jtp - 1] + i__;
				lujets_1.k[ii + j * 9000 - 9001] = lujets_1.k[i__ + j * 9000 - 9001];
				lujets_1.p[ii + j * 9000 - 9001] = lujets_1.p[i__ + j * 9000 - 9001];
				lujets_1.v[ii + j * 9000 - 9001] = lujets_1.v[i__ + j * 9000 - 9001];
			}
		}
		i__1 = hjjet1_1.npj[*jtp - 1];
		for (i__ = 1; i__ <= i__1; ++i__)
		{
			for (j = 1; j <= 5; ++j)
			{
				lujets_1.k[i__ + 1 + j * 9000 - 9001] = 0;
				lujets_1.v[i__ + 1 + j * 9000 - 9001] = 0.f;
			}
			i0 = i__;
			if (iex == 1 && (anim_1.isoft != 2 || anim_1.isflag != 0))
			{
				i0 = hjjet1_1.npj[*jtp - 1] - i__ + 1;
			}
			kk1 = hjjet1_1.kfpj[*jtp + i0 * 300 - 301];
			lujets_1.k[i__] = 2;
			lujets_1.k[i__ + 9000] = kk1;
			if (kk1 != 21 && kk1 != 0)
			{
				lujets_1.k[i__] = (abs(kk1) + ((iex << 1) - 1) * kk1) / 2 /
									  abs(kk1) +
								  1;
			}
			lujets_1.p[i__] = hjjet1_1.pjpx[*jtp + i0 * 300 - 301];
			lujets_1.p[i__ + 9000] = hjjet1_1.pjpy[*jtp + i0 * 300 - 301];
			lujets_1.p[i__ + 18000] = hjjet1_1.pjpz[*jtp + i0 * 300 - 301];
			lujets_1.p[i__ + 27000] = hjjet1_1.pjpe[*jtp + i0 * 300 - 301];
			lujets_1.p[i__ + 36000] = hjjet1_1.pjpm[*jtp + i0 * 300 - 301];
		}
		lujets_1.n += hjjet1_1.npj[*jtp - 1];
	}
	else if (*ntp == 2 && hjjet1_1.ntj[*jtp - 1] != 0)
	{
		jetot = hjjet1_1.ntj[*jtp - 1];
		iex = 1;
		if (abs(kf2) > 1000 && kf2 < 0 || abs(kf2) < 1000 && kf2 > 0)
		{
			iex = 0;
		}
		for (i__ = lujets_1.n; i__ >= 2; --i__)
		{
			for (j = 1; j <= 5; ++j)
			{
				ii = hjjet1_1.ntj[*jtp - 1] + i__;
				lujets_1.k[ii + j * 9000 - 9001] = lujets_1.k[i__ + j * 9000 - 9001];
				lujets_1.p[ii + j * 9000 - 9001] = lujets_1.p[i__ + j * 9000 - 9001];
				lujets_1.v[ii + j * 9000 - 9001] = lujets_1.v[i__ + j * 9000 - 9001];
			}
		}
		i__1 = hjjet1_1.ntj[*jtp - 1];
		for (i__ = 1; i__ <= i__1; ++i__)
		{
			for (j = 1; j <= 5; ++j)
			{
				lujets_1.k[i__ + 1 + j * 9000 - 9001] = 0;
				lujets_1.v[i__ + 1 + j * 9000 - 9001] = 0.f;
			}
			i0 = i__;
			if (iex == 1 && (anim_1.isoft != 2 || anim_1.isflag != 0))
			{
				i0 = hjjet1_1.ntj[*jtp - 1] - i__ + 1;
			}
			kk1 = hjjet1_1.kftj[*jtp + i0 * 300 - 301];
			lujets_1.k[i__] = 2;
			lujets_1.k[i__ + 9000] = kk1;
			if (kk1 != 21 && kk1 != 0)
			{
				lujets_1.k[i__] = (abs(kk1) + ((iex << 1) - 1) * kk1) / 2 /
									  abs(kk1) +
								  1;
			}
			lujets_1.p[i__] = hjjet1_1.pjtx[*jtp + i0 * 300 - 301];
			lujets_1.p[i__ + 9000] = hjjet1_1.pjty[*jtp + i0 * 300 - 301];
			lujets_1.p[i__ + 18000] = hjjet1_1.pjtz[*jtp + i0 * 300 - 301];
			lujets_1.p[i__ + 27000] = hjjet1_1.pjte[*jtp + i0 * 300 - 301];
			lujets_1.p[i__ + 36000] = hjjet1_1.pjtm[*jtp + i0 * 300 - 301];
		}
		lujets_1.n += hjjet1_1.ntj[*jtp - 1];
	}
	if (hparnt_1.ihpr2[0] > 0 && ranart_(&rndf77_1.nseed) <= hijdat_1.hidat[2])
	{
		hdat20 = hijdat_1.hidat[1];
		hpr150 = hparnt_1.hipr1[4];
		if (hparnt_1.ihpr2[7] == 0 && hparnt_1.ihpr2[2] == 0 &&
			hparnt_1.ihpr2[8] == 0)
		{
			hijdat_1.hidat[1] = 2.f;
		}
		if (hparnt_1.hint1[0] >= 1e3f && jetot == 0)
		{
			hijdat_1.hidat[1] = 3.f;
			hparnt_1.hipr1[4] = 5.f;
		}
		attrad_(ierror);
		hijdat_1.hidat[1] = hdat20;
		hparnt_1.hipr1[4] = hpr150;
	}
	else if (jetot == 0 && hparnt_1.ihpr2[0] > 0 && hparnt_1.hint1[0] >= 1e3f && ranart_(&rndf77_1.nseed) <= .8f)
	{
		hdat20 = hijdat_1.hidat[1];
		hpr150 = hparnt_1.hipr1[4];
		hijdat_1.hidat[1] = 3.f;
		hparnt_1.hipr1[4] = 5.f;
		if (hparnt_1.ihpr2[7] == 0 && hparnt_1.ihpr2[2] == 0 &&
			hparnt_1.ihpr2[8] == 0)
		{
			hijdat_1.hidat[1] = 2.f;
		}
		attrad_(ierror);
		hijdat_1.hidat[1] = hdat20;
		hparnt_1.hipr1[4] = hpr150;
	}
	if (*ierror != 0)
	{
		return 0;
	}

	if (anim_1.isoft != 2 || anim_1.isflag != 0)
	{
		luexec_();
	}
	return 0;
L1000:
	lujets_1.n = 1;
	lujets_1.k[0] = 1;
	lujets_1.k[9000] = hstrng_1.nfp[*jtp + 599];
	for (jj = 1; jj <= 5; ++jj)
	{
		lujets_1.p[jj * 9000 - 9000] = hstrng_1.pp[*jtp + jj * 300 - 301];
	}
	if (anim_1.isoft != 2 || anim_1.isflag != 0)
	{
		luexec_();
	}
	return 0;

L1200:
	lujets_1.n = 1;
	lujets_1.k[0] = 1;
	lujets_1.k[9000] = hstrng_1.nft[*jtp + 599];
	for (jj = 1; jj <= 5; ++jj)
	{
		lujets_1.p[jj * 9000 - 9000] = hstrng_1.pt[*jtp + jj * 300 - 301];
	}
	if (anim_1.isoft != 2 || anim_1.isflag != 0)
	{
		luexec_();
	}
	return 0;
}

int attrad_(int *ierror)
{
	int i__1;
	float r__1, r__2, r__3, r__4;

	double sqrt(double), acos(double), exp(double);

	static int i__, j, m;
	static float s, p1, p2, p3, p4, x1, x3;
	static int jl;
	static float sm, wp, bex, bey, bez, cth, phi, btt, ptg, pbt1, pbt2, pbt3,
		pbt4, ptg1, ptg2, ptg3;
	static int imin, imax;
	static float theta;
	extern int ar3jet_(float *, float *, float *, int *);
	static float fmfact;
	extern int arorie_(float *, float *, float *, int *);
	extern double ulangl_(float *, float *);
	extern int atrobo_(float *, float *, float *, float *, float *, int *, int *, int *);
	extern double ranart_(int *);

	*ierror = 0;
L40:
	sm = 0.f;
	jl = 1;
	i__1 = lujets_1.n - 1;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		r__1 = lujets_1.p[i__ + 35999];
		r__2 = lujets_1.p[i__ + 36000];
		s = (lujets_1.p[i__ + 26999] * lujets_1.p[i__ + 27000] - lujets_1.p[i__ - 1] * lujets_1.p[i__] - lujets_1.p[i__ + 8999] * lujets_1.p[i__ + 9000] - lujets_1.p[i__ + 17999] * lujets_1.p[i__ + 18000]) * 2.f + r__1 * r__1 + r__2 * r__2;
		if (s < 0.f)
		{
			s = 0.f;
		}
		wp = sqrt(s) - (lujets_1.p[i__ + 35999] + lujets_1.p[i__ + 36000]) *
						   1.5f;
		if (wp > sm)
		{
			pbt1 = lujets_1.p[i__ - 1] + lujets_1.p[i__];
			pbt2 = lujets_1.p[i__ + 8999] + lujets_1.p[i__ + 9000];
			pbt3 = lujets_1.p[i__ + 17999] + lujets_1.p[i__ + 18000];
			pbt4 = lujets_1.p[i__ + 26999] + lujets_1.p[i__ + 27000];
			r__1 = pbt1;
			r__2 = pbt2;
			r__3 = pbt3;
			r__4 = pbt4;
			btt = (r__1 * r__1 + r__2 * r__2 + r__3 * r__3) / (r__4 * r__4);
			jl = i__;
			sm = wp;
		}
	}
	r__1 = sm + (lujets_1.p[jl + 35999] + lujets_1.p[jl + 36000]) * 1.5f;
	s = r__1 * r__1;
	if (sm < hparnt_1.hipr1[4])
	{
		goto L2;
	}
	if (jl + 1 == lujets_1.n)
	{
		goto L190;
	}
	i__1 = jl + 2;
	for (j = lujets_1.n; j >= i__1; --j)
	{
		lujets_1.k[j] = lujets_1.k[j - 1];
		lujets_1.k[j + 9000] = lujets_1.k[j + 8999];
		for (m = 1; m <= 5; ++m)
		{
			lujets_1.p[j + 1 + m * 9000 - 9001] = lujets_1.p[j + m * 9000 -9001];
		}
	}
L190:
	++lujets_1.n;
	p1 = lujets_1.p[jl - 1] + lujets_1.p[jl];
	p2 = lujets_1.p[jl + 8999] + lujets_1.p[jl + 9000];
	p3 = lujets_1.p[jl + 17999] + lujets_1.p[jl + 18000];
	p4 = lujets_1.p[jl + 26999] + lujets_1.p[jl + 27000];
	bex = -p1 / p4;
	bey = -p2 / p4;
	bez = -p3 / p4;
	imin = jl;
	imax = jl + 1;
	atrobo_(&c_b26, &c_b26, &bex, &bey, &bez, &imin, &imax, ierror);
	if (*ierror != 0)
	{
		return 0;
	}
	r__1 = lujets_1.p[jl + 26999];
	r__2 = lujets_1.p[jl + 35999];
	cth = lujets_1.p[jl + 17999] / sqrt(r__1 * r__1 - r__2 * r__2);
	if (dabs(cth) > 1.f)
	{
		r__1 = -1.f, r__2 = min(1.f, cth);
		cth = max(r__1, r__2);
	}
	theta = acos(cth);
	phi = ulangl_(&lujets_1.p[jl - 1], &lujets_1.p[jl + 8999]);
	r__1 = -phi;
	atrobo_(&c_b26, &r__1, &c_b26, &c_b26, &c_b26, &imin, &imax, ierror);
	r__1 = -theta;
	atrobo_(&r__1, &c_b26, &c_b26, &c_b26, &c_b26, &imin, &imax, ierror);
L1:
	ar3jet_(&s, &x1, &x3, &jl);
	arorie_(&s, &x1, &x3, &jl);
	if (hijdat_1.hidat[1] > 0.f)
	{
		r__1 = lujets_1.p[jl - 1];
		r__2 = lujets_1.p[jl + 8999];
		ptg1 = sqrt(r__1 * r__1 + r__2 * r__2);
		r__1 = lujets_1.p[jl];
		r__2 = lujets_1.p[jl + 9000];
		ptg2 = sqrt(r__1 * r__1 + r__2 * r__2);
		r__1 = lujets_1.p[jl + 1];
		r__2 = lujets_1.p[jl + 9001];
		ptg3 = sqrt(r__1 * r__1 + r__2 * r__2);
		r__1 = max(ptg1, ptg2);
		ptg = max(r__1, ptg3);
		if (ptg > hijdat_1.hidat[1])
		{
			r__1 = ptg;
			r__2 = hijdat_1.hidat[1];
			r__3 = hparnt_1.hipr1[1];
			fmfact = exp(-(r__1 * r__1 - r__2 * r__2) / (r__3 * r__3));
			if (ranart_(&rndf77_1.nseed) > fmfact)
			{
				goto L1;
			}
		}
	}
	imin = jl;
	imax = jl + 2;
	r__1 = -bex;
	r__2 = -bey;
	r__3 = -bez;
	atrobo_(&theta, &phi, &r__1, &r__2, &r__3, &imin, &imax, ierror);
	if (*ierror != 0)
	{
		return 0;
	}
	lujets_1.k[jl + 1] = lujets_1.k[jl];
	lujets_1.k[jl + 9001] = lujets_1.k[jl + 9000];
	lujets_1.k[jl + 18001] = lujets_1.k[jl + 18000];
	lujets_1.k[jl + 27001] = lujets_1.k[jl + 27000];
	lujets_1.k[jl + 36001] = lujets_1.k[jl + 36000];
	lujets_1.p[jl + 36001] = lujets_1.p[jl + 36000];
	lujets_1.k[jl] = 2;
	lujets_1.k[jl + 9000] = 21;
	lujets_1.k[jl + 18000] = 0;
	lujets_1.k[jl + 27000] = 0;
	lujets_1.k[jl + 36000] = 0;
	lujets_1.p[jl + 36000] = 0.f;
	if (sm >= hparnt_1.hipr1[4])
	{
		goto L40;
	}
L2:
	lujets_1.k[0] = 2;
	lujets_1.k[18000] = 0;
	lujets_1.k[27000] = 0;
	lujets_1.k[36000] = 0;
	lujets_1.k[lujets_1.n - 1] = 1;
	lujets_1.k[lujets_1.n + 17999] = 0;
	lujets_1.k[lujets_1.n + 26999] = 0;
	lujets_1.k[lujets_1.n + 35999] = 0;
	return 0;
}

int ar3jet_(float *s, float *x1, float *x3, int *jl)
{
	float r__1;
	double d__1, d__2, d__3, d__4;

	double sqrt(double), log(double), pow_dd(double *, double *), exp(double);

	static float a, c__, d__, y, x2, fg, sm1, sm3, xt2;
	static int neg;
	static float yma, exp1, exp3, xt2m, ymax;
	static int ntry;
	extern double ranart_(int *);

	c__ = .33333333333333331f;
	if (lujets_1.k[*jl + 8999] != 21 && lujets_1.k[*jl + 9000] != 21)
	{
		c__ = .29629629629629628f;
	}
	exp1 = 3.f;
	exp3 = 3.f;
	if (lujets_1.k[*jl + 8999] != 21)
	{
		exp1 = 2.f;
	}
	if (lujets_1.k[*jl + 9000] != 21)
	{
		exp3 = 2.f;
	}
	a = .057599999999999998f / *s;
	yma = log(.5f / sqrt(a) + sqrt(.25f / a - 1));
	d__ = c__ * 4.f * yma;
	r__1 = lujets_1.p[*jl + 35999];
	sm1 = r__1 * r__1 / *s;
	r__1 = lujets_1.p[*jl + 36000];
	sm3 = r__1 * r__1 / *s;
	xt2m = (1.f - sqrt(sm1) * 2.f + sm1 - sm3) * (1.f - sqrt(sm3) * 2.f - sm1 + sm3);
	xt2m = min(.25f, xt2m);
	ntry = 0;
L1:
	if (ntry == 5000)
	{
		*x1 = (sqrt(sm1) * 2.f + 1.f + sm1 - sm3) * .5f;
		*x3 = (sqrt(sm3) * 2.f + 1.f - sm1 + sm3) * .5f;
		return 0;
	}
	++ntry;
	d__1 = (double)(xt2m / a);
	d__3 = (double)ranart_(&rndf77_1.nseed);
	d__4 = (double)(1.f / d__);
	d__2 = (double)pow_dd(&d__3, &d__4);
	xt2 = a * pow_dd(&d__1, &d__2);
	ymax = log(.5f / sqrt(xt2) + sqrt(.25f / xt2 - 1.f));
	y = (ranart_(&rndf77_1.nseed) * 2.f - 1.f) * ymax;
	*x1 = 1.f - sqrt(xt2) * exp(y);
	*x3 = 1.f - sqrt(xt2) * exp(-y);
	x2 = 2.f - *x1 - *x3;
	neg = 0;
	if (lujets_1.k[*jl + 8999] != 21 || lujets_1.k[*jl + 9000] != 21)
	{
		if ((1.f - *x1) * (1.f - x2) * (1.f - *x3) - x2 * sm1 * (1.f - *x1) -
					x2 * sm3 * (1.f - *x3) <=
				0.f ||
			*x1 <= sqrt(sm1) * 2.f - sm1 + sm3 || *x3 <= sqrt(sm3) * 2.f - sm3 + sm1)
		{
			neg = 1;
		}
		*x1 = *x1 + sm1 - sm3;
		*x3 = *x3 - sm1 + sm3;
	}
	if (neg == 1)
	{
		goto L1;
	}
	d__1 = (double)(*x1);
	d__2 = (double)exp1;
	d__3 = (double)(*x3);
	d__4 = (double)exp3;
	fg = ymax * 2.f * c__ * (pow_dd(&d__1, &d__2) + pow_dd(&d__3, &d__4)) /
		 d__;
	xt2m = xt2;
	if (fg < ranart_(&rndf77_1.nseed))
	{
		goto L1;
	}
	return 0;
}

int arorie_(float *s, float *x1, float *x3, int *jl)
{
	float r__1, r__2, r__3, r__4, r__5;

	double sqrt(double), acos(double), cos(double), sin(double);

	static float w, e1, e3, p1, p3, x2, pt1, pt3, pz1, pz3, del, bet, psi,
		cbet;
	extern double ulangl_(float *, float *), ranart_(int *);

	w = sqrt(*s);
	x2 = 2.f - *x1 - *x3;
	e1 = *x1 * .5f * w;
	e3 = *x3 * .5f * w;
	r__1 = e1;
	r__2 = lujets_1.p[*jl + 35999];
	p1 = sqrt(r__1 * r__1 - r__2 * r__2);
	r__1 = e3;
	r__2 = lujets_1.p[*jl + 36000];
	p3 = sqrt(r__1 * r__1 - r__2 * r__2);
	cbet = 1.f;
	if (p1 > 0.f && p3 > 0.f)
	{
		r__1 = lujets_1.p[*jl + 35999];
		r__2 = lujets_1.p[*jl + 36000];
		cbet = (r__1 * r__1 + r__2 * r__2 + e1 * 2.f * e3 - *s * (1.f - x2)) /
			   (p1 * 2.f * p3);
	}
	if (dabs(cbet) > 1.f)
	{
		r__1 = -1.f, r__2 = min(1.f, cbet);
		cbet = max(r__1, r__2);
	}
	bet = acos(cbet);
	if (p1 >= p3)
	{
		r__2 = p1;
		r__3 = p3;
		r__1 = r__2 * r__2 + r__3 * r__3 * cos(bet * 2.f);
		r__5 = p3;
		r__4 = -(r__5 * r__5) * sin(bet * 2.f);
		psi = ulangl_(&r__1, &r__4) * .5f;
		pt1 = p1 * sin(psi);
		pz1 = p1 * cos(psi);
		pt3 = p3 * sin(psi + bet);
		pz3 = p3 * cos(psi + bet);
	}
	else if (p3 > p1)
	{
		r__2 = p3;
		r__3 = p1;
		r__1 = r__2 * r__2 + r__3 * r__3 * cos(bet * 2.f);
		r__5 = p1;
		r__4 = -(r__5 * r__5) * sin(bet * 2.f);
		psi = ulangl_(&r__1, &r__4) * .5f;
		pt1 = p1 * sin(bet + psi);
		pz1 = -p1 * cos(bet + psi);
		pt3 = p3 * sin(psi);
		pz3 = -p3 * cos(psi);
	}
	del = hparnt_1.hipr1[39] * 2.f * ranart_(&rndf77_1.nseed);
	lujets_1.p[*jl + 26999] = e1;
	lujets_1.p[*jl - 1] = pt1 * sin(del);
	lujets_1.p[*jl + 8999] = -pt1 * cos(del);
	lujets_1.p[*jl + 17999] = pz1;
	lujets_1.p[*jl + 27001] = e3;
	lujets_1.p[*jl + 1] = pt3 * sin(del);
	lujets_1.p[*jl + 9001] = -pt3 * cos(del);
	lujets_1.p[*jl + 18001] = pz3;
	lujets_1.p[*jl + 27000] = w - e1 - e3;
	lujets_1.p[*jl] = -lujets_1.p[*jl - 1] - lujets_1.p[*jl + 1];
	lujets_1.p[*jl + 9000] = -lujets_1.p[*jl + 8999] - lujets_1.p[*jl + 9001];
	lujets_1.p[*jl + 18000] = -lujets_1.p[*jl + 17999] - lujets_1.p[*jl +
																	18001];
	return 0;
}

int atrobo_(float *the, float *phi, float *bex, float *bey, float *bez, int *imin, int *imax, int *ierror)
{
	int i__1;
	float r__1, r__2, r__3;
	double d__1, d__2, d__3;

	double cos(double), sin(double), sqrt(double);

	static int i__, j;
	static double dp[4];
	static float pv[3];
	static double dga;
	static float rot[9];
	static double dga2, dbep, dbex, dbey, dbez, dgabep;

	*ierror = 0;
	if (*imin <= 0 || *imax > lujets_1.n || *imin > *imax)
	{
		return 0;
	}
	r__1 = *the;
	r__2 = *phi;
	if (r__1 * r__1 + r__2 * r__2 > 1e-20f)
	{
		rot[0] = cos(*the) * cos(*phi);
		rot[3] = -sin(*phi);
		rot[6] = sin(*the) * cos(*phi);
		rot[1] = cos(*the) * sin(*phi);
		rot[4] = cos(*phi);
		rot[7] = sin(*the) * sin(*phi);
		rot[2] = -sin(*the);
		rot[5] = 0.f;
		rot[8] = cos(*the);
		i__1 = *imax;
		for (i__ = *imin; i__ <= i__1; ++i__)
		{
			for (j = 1; j <= 3; ++j)
			{
				pv[j - 1] = lujets_1.p[i__ + j * 9000 - 9001];
			}
			for (j = 1; j <= 3; ++j)
			{
				lujets_1.p[i__ + j * 9000 - 9001] = rot[j - 1] * pv[0] + rot[j + 2] * pv[1] + rot[j + 5] * pv[2];
			}
		}
	}
	r__1 = *bex;
	r__2 = *bey;
	r__3 = *bez;
	if (r__1 * r__1 + r__2 * r__2 + r__3 * r__3 > 1e-20f)
	{
		dbex = (double)(*bex);
		dbey = (double)(*bey);
		dbez = (double)(*bez);
		d__1 = dbex;
		d__2 = dbey;
		d__3 = dbez;
		dga2 = 1. - d__1 * d__1 - d__2 * d__2 - d__3 * d__3;
		if (dga2 <= 0.)
		{
			*ierror = 1;
			return 0;
		}
		dga = 1. / sqrt(dga2);
		i__1 = *imax;
		for (i__ = *imin; i__ <= i__1; ++i__)
		{
			for (j = 1; j <= 4; ++j)
			{
				dp[j - 1] = (double)lujets_1.p[i__ + j * 9000 - 9001];
			}
			dbep = dbex * dp[0] + dbey * dp[1] + dbez * dp[2];
			dgabep = dga * (dga * dbep / (dga + 1.) + dp[3]);
			lujets_1.p[i__ - 1] = (float)(dp[0] + dgabep * dbex);
			lujets_1.p[i__ + 8999] = (float)(dp[1] + dgabep * dbey);
			lujets_1.p[i__ + 17999] = (float)(dp[2] + dgabep * dbez);
			lujets_1.p[i__ + 26999] = (float)(dga * (dp[3] + dbep));
		}
	}
	return 0;
}

int hijhrd_(int *jp, int *jt, int *jout, int *jflg, int *iopjt0)
{
	int i__1, i__2;
	float r__1, r__2, r__3, r__4;

	double sqrt(double), exp(double);
	int s_wsle(cilist *), do_lio(int *, int *, char *, int),
		e_wsle(void);

	static int i__, j, l0, jj, ip[200], lp, it[200];
	static float qm;
	static int lt;
	static float wm, wp, sw;
	static int ip1, ip2, it1, it2, is7, is8, ipb[50], lpb, itb[50];
	static float epm;
	static int ltb;
	static float epp;
	static int ipq[50];
	static float etp;
	static int itq[50];
	static float etm, pep, pet;
	static int jpp, lpq, jtt, ltq;
	static float pxp, pyp, pzp, pxt, pyt, pzt, sxx;
	static int mxjt, mxsg, mxsj, miss, misp, mist;
	static float ampx, amtx, ecut1, ecut2;
	static int isub11, isub12, isub28;
	static float qmass2;
	static int iinird;
	static float pinird;
	extern double ranart_(int *);
	static int iopjet;
	extern int pythia_(void);
	extern double ulmass_(int *);

	static cilist io___330 = {0, 6, 0, 0, 0};
	static cilist io___338 = {0, 6, 0, 0, 0};
	static cilist io___339 = {0, 6, 0, 0, 0};
	static cilist io___340 = {0, 6, 0, 0, 0};
	static cilist io___342 = {0, 6, 0, 0, 0};
	static cilist io___343 = {0, 6, 0, 0, 0};
	static cilist io___344 = {0, 6, 0, 0, 0};
	static cilist io___345 = {0, 6, 0, 0, 0};
	static cilist io___346 = {0, 6, 0, 0, 0};
	static cilist io___347 = {0, 6, 0, 0, 0};
	static cilist io___348 = {0, 6, 0, 0, 0};

	mxjt = 500;
	mxsg = 900;
	mxsj = 100;
	*jflg = 0;
	hparnt_1.ihnt2[10] = *jp;
	hparnt_1.ihnt2[11] = *jt;

	iopjet = *iopjt0;
	if (iopjet == 1 && (hstrng_1.nfp[*jp + 1499] != 0 || hstrng_1.nft[*jt +
																	  1499] != 0))
	{
		iopjet = 0;
	}
	if (*jp > hparnt_1.ihnt2[0] || *jt > hparnt_1.ihnt2[2])
	{
		return 0;
	}
	if (hstrng_1.nfp[*jp + 1499] < 0 || hstrng_1.nft[*jt + 1499] < 0)
	{
		return 0;
	}
	if (*jout == 0)
	{
		epp = hstrng_1.pp[*jp + 899] + hstrng_1.pp[*jp + 599];
		epm = hstrng_1.pp[*jp + 899] - hstrng_1.pp[*jp + 599];
		etp = hstrng_1.pt[*jt + 899] + hstrng_1.pt[*jt + 599];
		etm = hstrng_1.pt[*jt + 899] - hstrng_1.pt[*jt + 599];
		if (epp < 0.f)
		{
			goto L1000;
		}
		if (epm < 0.f)
		{
			goto L1000;
		}
		if (etp < 0.f)
		{
			goto L1000;
		}
		if (etm < 0.f)
		{
			goto L1000;
		}
		if (epp / (epm + .01f) <= etp / (etm + .01f))
		{
			return 0;
		}
	}
	ecut1 = hparnt_1.hipr1[0] + hparnt_1.hipr1[7] + hstrng_1.pp[*jp + 3899] +
			hstrng_1.pp[*jp + 4199];
	ecut2 = hparnt_1.hipr1[0] + hparnt_1.hipr1[7] + hstrng_1.pt[*jt + 3899] +
			hstrng_1.pt[*jt + 4199];
	if (hstrng_1.pp[*jp + 899] <= ecut1)
	{
		hstrng_1.nfp[*jp + 1499] = -(i__1 = hstrng_1.nfp[*jp + 1499], abs(i__1));
		return 0;
	}
	if (hstrng_1.pt[*jt + 899] <= ecut2)
	{
		hstrng_1.nft[*jt + 1499] = -(i__1 = hstrng_1.nft[*jt + 1499], abs(i__1));
		return 0;
	}
	miss = 0;
	misp = 0;
	mist = 0;

	if (hstrng_1.nfp[*jp + 2699] == 0 && hstrng_1.nft[*jt + 2699] == 0)
	{
		pyint1_1.mint[43] = hpint_1.mint4;
		pyint1_1.mint[44] = hpint_1.mint5;
		pyint5_1.xsec[0] = hpint_1.atxs[0];
		pyint5_1.xsec[11] = hpint_1.atxs[11];
		pyint5_1.xsec[12] = hpint_1.atxs[12];
		pyint5_1.xsec[28] = hpint_1.atxs[28];
		for (i__ = 1; i__ <= 20; ++i__)
		{
			pyint2_1.coef[i__ * 200 - 190] = hpint_1.atco[i__ * 200 - 190];
			pyint2_1.coef[i__ * 200 - 189] = hpint_1.atco[i__ * 200 - 189];
			pyint2_1.coef[i__ * 200 - 173] = hpint_1.atco[i__ * 200 - 173];
		}
	}
	else
	{
		isub11 = 0;
		isub12 = 0;
		isub28 = 0;
		if (pyint5_1.xsec[11] != 0.f)
		{
			isub11 = 1;
		}
		if (pyint5_1.xsec[12] != 0.f)
		{
			isub12 = 1;
		}
		if (pyint5_1.xsec[28] != 0.f)
		{
			isub28 = 1;
		}
		pyint1_1.mint[43] = hpint_1.mint4 - isub11 - isub12 - isub28;
		pyint1_1.mint[44] = hpint_1.mint5 - isub11 - isub12 - isub28;
		pyint5_1.xsec[0] = hpint_1.atxs[0] - hpint_1.atxs[11] - hpint_1.atxs[12] - hpint_1.atxs[28];
		pyint5_1.xsec[11] = 0.f;
		pyint5_1.xsec[12] = 0.f;
		pyint5_1.xsec[28] = 0.f;
		for (i__ = 1; i__ <= 20; ++i__)
		{
			pyint2_1.coef[i__ * 200 - 190] = 0.f;
			pyint2_1.coef[i__ * 200 - 189] = 0.f;
			pyint2_1.coef[i__ * 200 - 173] = 0.f;
		}
	}
L155:
	pythia_();
	jj = pyint1_1.mint[30];
	if (jj != 1)
	{
		goto L155;
	}
	if (lujets_1.k[9006] == -lujets_1.k[9007])
	{
		r__1 = lujets_1.p[27006] + lujets_1.p[27007];
		r__2 = lujets_1.p[6] + lujets_1.p[7];
		r__3 = lujets_1.p[9006] + lujets_1.p[9007];
		r__4 = lujets_1.p[18006] + lujets_1.p[18007];
		qmass2 = r__1 * r__1 - r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
		qm = ulmass_(&lujets_1.k[9006]);
		r__1 = qm * 2.f + hparnt_1.hipr1[0];
		if (qmass2 < r__1 * r__1)
		{
			goto L155;
		}
	}
	pxp = hstrng_1.pp[*jp - 1] - lujets_1.p[2];
	pyp = hstrng_1.pp[*jp + 299] - lujets_1.p[9002];
	pzp = hstrng_1.pp[*jp + 599] - lujets_1.p[18002];
	pep = hstrng_1.pp[*jp + 899] - lujets_1.p[27002];
	pxt = hstrng_1.pt[*jt - 1] - lujets_1.p[3];
	pyt = hstrng_1.pt[*jt + 299] - lujets_1.p[9003];
	pzt = hstrng_1.pt[*jt + 599] - lujets_1.p[18003];
	pet = hstrng_1.pt[*jt + 899] - lujets_1.p[27003];
	if (pep <= ecut1)
	{
		++misp;
		if (misp < 50)
		{
			goto L155;
		}
		hstrng_1.nfp[*jp + 1499] = -(i__1 = hstrng_1.nfp[*jp + 1499], abs(i__1));
		return 0;
	}
	if (pet <= ecut2)
	{
		++mist;
		if (mist < 50)
		{
			goto L155;
		}
		hstrng_1.nft[*jt + 1499] = -(i__1 = hstrng_1.nft[*jt + 1499], abs(i__1));
		return 0;
	}
	wp = pep + pzp + pet + pzt;
	wm = pep - pzp + pet - pzt;
	if (wp < 0.f || wm < 0.f)
	{
		++miss;
		if (miss < 50)
		{
			goto L155;
		}
		return 0;
	}
	sw = wp * wm;
	r__1 = ecut1 - hparnt_1.hipr1[7];
	r__2 = pxp;
	r__3 = pyp;
	ampx = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + .01f);
	r__1 = ecut2 - hparnt_1.hipr1[7];
	r__2 = pxt;
	r__3 = pyt;
	amtx = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + .01f);
	r__1 = ampx + amtx;
	sxx = r__1 * r__1;
	if (sw < sxx || pyint1_1.vint[42] < hparnt_1.hipr1[0])
	{
		++miss;
		if (miss < 50)
		{
			goto L155;
		}
		return 0;
	}
	hparnt_1.hint1[40] = lujets_1.p[6];
	hparnt_1.hint1[41] = lujets_1.p[9006];
	hparnt_1.hint1[42] = lujets_1.p[18006];
	hparnt_1.hint1[43] = lujets_1.p[27006];
	hparnt_1.hint1[44] = lujets_1.p[36006];
	r__1 = lujets_1.p[6];
	r__2 = lujets_1.p[9006];
	hparnt_1.hint1[45] = sqrt(r__1 * r__1 + r__2 * r__2);
	hparnt_1.hint1[50] = lujets_1.p[7];
	hparnt_1.hint1[51] = lujets_1.p[9007];
	hparnt_1.hint1[52] = lujets_1.p[18007];
	hparnt_1.hint1[53] = lujets_1.p[27007];
	hparnt_1.hint1[54] = lujets_1.p[36007];
	r__1 = lujets_1.p[7];
	r__2 = lujets_1.p[9007];
	hparnt_1.hint1[55] = sqrt(r__1 * r__1 + r__2 * r__2);
	hparnt_1.ihnt2[13] = lujets_1.k[9006];
	hparnt_1.ihnt2[14] = lujets_1.k[9007];

	pinird = (1.f - exp((pyint1_1.vint[46] - hijdat_1.hidat[0]) * -2.f)) / (exp((pyint1_1.vint[46] - hijdat_1.hidat[0]) * -2.f) + 1.f);
	iinird = 0;
	if (ranart_(&rndf77_1.nseed) <= pinird)
	{
		iinird = 1;
	}
	if (lujets_1.k[9006] == -lujets_1.k[9007])
	{
		goto L190;
	}
	if (lujets_1.k[9006] == 21 && lujets_1.k[9007] == 21 && iopjet == 1)
	{
		goto L190;
	}
	*jflg = 2;
	jpp = 0;
	lpq = 0;
	lpb = 0;
	jtt = 0;
	ltq = 0;
	ltb = 0;
	is7 = 0;
	is8 = 0;
	hparnt_1.hint1[46] = 0.f;
	hparnt_1.hint1[47] = 0.f;
	hparnt_1.hint1[48] = 0.f;
	hparnt_1.hint1[49] = 0.f;
	hparnt_1.hint1[66] = 0.f;
	hparnt_1.hint1[67] = 0.f;
	hparnt_1.hint1[68] = 0.f;
	hparnt_1.hint1[69] = 0.f;
	i__1 = lujets_1.n;
	for (i__ = 9; i__ <= i__1; ++i__)
	{
		if (lujets_1.k[i__ + 17999] == 7)
		{
			hparnt_1.hint1[46] += lujets_1.p[i__ - 1];
			hparnt_1.hint1[47] += lujets_1.p[i__ + 8999];
			hparnt_1.hint1[48] += lujets_1.p[i__ + 17999];
			hparnt_1.hint1[49] += lujets_1.p[i__ + 26999];
		}
		if (lujets_1.k[i__ + 17999] == 8)
		{
			hparnt_1.hint1[66] += lujets_1.p[i__ - 1];
			hparnt_1.hint1[67] += lujets_1.p[i__ + 8999];
			hparnt_1.hint1[68] += lujets_1.p[i__ + 17999];
			hparnt_1.hint1[69] += lujets_1.p[i__ + 26999];
		}
		if (lujets_1.k[i__ + 8999] > 21 && lujets_1.k[i__ + 8999] <= 30)
		{
			++hjjet4_1.ndr;
			hjjet4_1.iadr[hjjet4_1.ndr - 1] = *jp;
			hjjet4_1.iadr[hjjet4_1.ndr + 899] = *jt;
			hjjet4_1.kfdr[hjjet4_1.ndr - 1] = lujets_1.k[i__ + 8999];
			hjjet4_1.pdr[hjjet4_1.ndr - 1] = lujets_1.p[i__ - 1];
			hjjet4_1.pdr[hjjet4_1.ndr + 899] = lujets_1.p[i__ + 8999];
			hjjet4_1.pdr[hjjet4_1.ndr + 1799] = lujets_1.p[i__ + 17999];
			hjjet4_1.pdr[hjjet4_1.ndr + 2699] = lujets_1.p[i__ + 26999];
			hjjet4_1.pdr[hjjet4_1.ndr + 3599] = lujets_1.p[i__ + 35999];
			xydr_1.rtdr[hjjet4_1.ndr - 1] = (hjcrdn_1.yp[*jp * 3 - 3] +
											 hjcrdn_1.yt[*jt * 3 - 3]) *
											.5f;
			xydr_1.rtdr[hjjet4_1.ndr + 150000] = (hjcrdn_1.yp[*jp * 3 - 2] +
												  hjcrdn_1.yt[*jt * 3 - 2]) *
												 .5f;
		}
		if (lujets_1.k[i__ + 17999] == 7 || lujets_1.k[i__ + 17999] == 3)
		{
			if (lujets_1.k[i__ + 17999] == 7 && lujets_1.k[i__ + 8999] != 21 && lujets_1.k[i__ + 8999] == lujets_1.k[9006] && is7 == 0)
			{
				hstrng_1.pp[*jp + 2699] = lujets_1.p[i__ - 1];
				hstrng_1.pp[*jp + 2999] = lujets_1.p[i__ + 8999];
				hstrng_1.pp[*jp + 3299] = lujets_1.p[i__ + 17999];
				pzp += lujets_1.p[i__ + 17999];
				pep += lujets_1.p[i__ + 26999];
				hstrng_1.nfp[*jp + 2699] = 1;
				is7 = 1;
			}
			if (lujets_1.k[i__ + 17999] == 3 && (lujets_1.k[i__ + 8999] != 21 || iinird == 0))
			{
				pxp += lujets_1.p[i__ - 1];
				pyp += lujets_1.p[i__ + 8999];
				pzp += lujets_1.p[i__ + 17999];
				pep += lujets_1.p[i__ + 26999];
			}
			++jpp;
			ip[jpp - 1] = i__;
			ip[jpp + 99] = 0;
			if (lujets_1.k[i__ + 8999] != 21)
			{
				if (lujets_1.k[i__ + 8999] > 0)
				{
					++lpq;
					ipq[lpq - 1] = jpp;
					ip[jpp + 99] = lpq;
				}
				else if (lujets_1.k[i__ + 8999] < 0)
				{
					++lpb;
					ipb[lpb - 1] = jpp;
					ip[jpp + 99] = -lpb;
				}
			}
		}
		else if (lujets_1.k[i__ + 17999] == 8 || lujets_1.k[i__ + 17999] ==
													 4)
		{
			if (lujets_1.k[i__ + 17999] == 8 && lujets_1.k[i__ + 8999] != 21 && lujets_1.k[i__ + 8999] == lujets_1.k[9007] && is8 == 0)
			{
				hstrng_1.pt[*jt + 2699] = lujets_1.p[i__ - 1];
				hstrng_1.pt[*jt + 2999] = lujets_1.p[i__ + 8999];
				hstrng_1.pt[*jt + 3299] = lujets_1.p[i__ + 17999];
				pzt += lujets_1.p[i__ + 17999];
				pet += lujets_1.p[i__ + 26999];
				hstrng_1.nft[*jt + 2699] = 1;
				is8 = 1;
			}
			if (lujets_1.k[i__ + 17999] == 4 && (lujets_1.k[i__ + 8999] != 21 || iinird == 0))
			{
				pxt += lujets_1.p[i__ - 1];
				pyt += lujets_1.p[i__ + 8999];
				pzt += lujets_1.p[i__ + 17999];
				pet += lujets_1.p[i__ + 26999];
			}
			++jtt;
			it[jtt - 1] = i__;
			it[jtt + 99] = 0;
			if (lujets_1.k[i__ + 8999] != 21)
			{
				if (lujets_1.k[i__ + 8999] > 0)
				{
					++ltq;
					itq[ltq - 1] = jtt;
					it[jtt + 99] = ltq;
				}
				else if (lujets_1.k[i__ + 8999] < 0)
				{
					++ltb;
					itb[ltb - 1] = jtt;
					it[jtt + 99] = -ltb;
				}
			}
		}
	}

	if (lpq != lpb || ltq != ltb)
	{
		++miss;
		if (miss <= 50)
		{
			goto L155;
		}
		s_wsle(&io___330);
		do_lio(&c__9, &c__1, " Q -QBAR NOT MATCHED IN HIJHRD", (int)30);
		e_wsle();
		*jflg = 0;
		return 0;
	}
	j = 0;
L181:
	++j;
	if (j > jpp)
	{
		goto L182;
	}
	if (ip[j + 99] == 0)
	{
		goto L181;
	}
	else if (ip[j + 99] != 0)
	{
		lp = (i__1 = ip[j + 99], abs(i__1));
		ip1 = ip[j - 1];
		ip2 = ip[j + 99];
		ip[j - 1] = ip[ipq[lp - 1] - 1];
		ip[j + 99] = ip[ipq[lp - 1] + 99];
		ip[ipq[lp - 1] - 1] = ip1;
		ip[ipq[lp - 1] + 99] = ip2;
		if (ip2 > 0)
		{
			ipq[ip2 - 1] = ipq[lp - 1];
		}
		else if (ip2 < 0)
		{
			ipb[-ip2 - 1] = ipq[lp - 1];
		}
		ip1 = ip[j];
		ip2 = ip[j + 100];
		ip[j] = ip[ipb[lp - 1] - 1];
		ip[j + 100] = ip[ipb[lp - 1] + 99];
		ip[ipb[lp - 1] - 1] = ip1;
		ip[ipb[lp - 1] + 99] = ip2;
		if (ip2 > 0)
		{
			ipq[ip2 - 1] = ipb[lp - 1];
		}
		else if (ip2 < 0)
		{
			ipb[-ip2 - 1] = ipb[lp - 1];
		}
		++j;
		goto L181;
	}
L182:
	j = 0;
L183:
	++j;
	if (j > jtt)
	{
		goto L184;
	}
	if (it[j + 99] == 0)
	{
		goto L183;
	}
	else if (it[j + 99] != 0)
	{
		lt = (i__1 = it[j + 99], abs(i__1));
		it1 = it[j - 1];
		it2 = it[j + 99];
		it[j - 1] = it[itq[lt - 1] - 1];
		it[j + 99] = it[itq[lt - 1] + 99];
		it[itq[lt - 1] - 1] = it1;
		it[itq[lt - 1] + 99] = it2;
		if (it2 > 0)
		{
			itq[it2 - 1] = itq[lt - 1];
		}
		else if (it2 < 0)
		{
			itb[-it2 - 1] = itq[lt - 1];
		}
		it1 = it[j];
		it2 = it[j + 100];
		it[j] = it[itb[lt - 1] - 1];
		it[j + 100] = it[itb[lt - 1] + 99];
		it[itb[lt - 1] - 1] = it1;
		it[itb[lt - 1] + 99] = it2;
		if (it2 > 0)
		{
			itq[it2 - 1] = itb[lt - 1];
		}
		else if (it2 < 0)
		{
			itb[-it2 - 1] = itb[lt - 1];
		}
		++j;
		goto L183;
	}
L184:
	if (hjjet1_1.npj[*jp - 1] + jpp > mxjt || hjjet1_1.ntj[*jt - 1] + jtt >
												  mxjt)
	{
		*jflg = 0;
		s_wsle(&io___338);
		do_lio(&c__9, &c__1, "number of partons per string exceeds", (int)36);
		e_wsle();
		s_wsle(&io___339);
		do_lio(&c__9, &c__1, "the common block size", (int)21);
		e_wsle();
		return 0;
	}
	i__1 = jpp;
	for (j = 1; j <= i__1; ++j)
	{
		hjjet1_1.kfpj[*jp + (hjjet1_1.npj[*jp - 1] + j) * 300 - 301] =
			lujets_1.k[ip[j - 1] + 8999];
		hjjet1_1.pjpx[*jp + (hjjet1_1.npj[*jp - 1] + j) * 300 - 301] =
			lujets_1.p[ip[j - 1] - 1];
		hjjet1_1.pjpy[*jp + (hjjet1_1.npj[*jp - 1] + j) * 300 - 301] =
			lujets_1.p[ip[j - 1] + 8999];
		hjjet1_1.pjpz[*jp + (hjjet1_1.npj[*jp - 1] + j) * 300 - 301] =
			lujets_1.p[ip[j - 1] + 17999];
		hjjet1_1.pjpe[*jp + (hjjet1_1.npj[*jp - 1] + j) * 300 - 301] =
			lujets_1.p[ip[j - 1] + 26999];
		hjjet1_1.pjpm[*jp + (hjjet1_1.npj[*jp - 1] + j) * 300 - 301] =
			lujets_1.p[ip[j - 1] + 35999];
	}
	hjjet1_1.npj[*jp - 1] += jpp;
	i__1 = jtt;
	for (j = 1; j <= i__1; ++j)
	{
		hjjet1_1.kftj[*jt + (hjjet1_1.ntj[*jt - 1] + j) * 300 - 301] =
			lujets_1.k[it[j - 1] + 8999];
		hjjet1_1.pjtx[*jt + (hjjet1_1.ntj[*jt - 1] + j) * 300 - 301] =
			lujets_1.p[it[j - 1] - 1];
		hjjet1_1.pjty[*jt + (hjjet1_1.ntj[*jt - 1] + j) * 300 - 301] =
			lujets_1.p[it[j - 1] + 8999];
		hjjet1_1.pjtz[*jt + (hjjet1_1.ntj[*jt - 1] + j) * 300 - 301] =
			lujets_1.p[it[j - 1] + 17999];
		hjjet1_1.pjte[*jt + (hjjet1_1.ntj[*jt - 1] + j) * 300 - 301] =
			lujets_1.p[it[j - 1] + 26999];
		hjjet1_1.pjtm[*jt + (hjjet1_1.ntj[*jt - 1] + j) * 300 - 301] =
			lujets_1.p[it[j - 1] + 35999];
	}
	hjjet1_1.ntj[*jt - 1] += jtt;
	goto L900;
L190:
	*jflg = 3;
	if (lujets_1.k[9006] != 21 && lujets_1.k[9007] != 21 && lujets_1.k[9006] * lujets_1.k[9007] > 0)
	{
		goto L155;
	}
	jpp = 0;
	lpq = 0;
	lpb = 0;
	i__1 = lujets_1.n;
	for (i__ = 9; i__ <= i__1; ++i__)
	{
		if (lujets_1.k[i__ + 8999] > 21 && lujets_1.k[i__ + 8999] <= 30)
		{
			++hjjet4_1.ndr;
			hjjet4_1.iadr[hjjet4_1.ndr - 1] = *jp;
			hjjet4_1.iadr[hjjet4_1.ndr + 899] = *jt;
			hjjet4_1.kfdr[hjjet4_1.ndr - 1] = lujets_1.k[i__ + 8999];
			hjjet4_1.pdr[hjjet4_1.ndr - 1] = lujets_1.p[i__ - 1];
			hjjet4_1.pdr[hjjet4_1.ndr + 899] = lujets_1.p[i__ + 8999];
			hjjet4_1.pdr[hjjet4_1.ndr + 1799] = lujets_1.p[i__ + 17999];
			hjjet4_1.pdr[hjjet4_1.ndr + 2699] = lujets_1.p[i__ + 26999];
			hjjet4_1.pdr[hjjet4_1.ndr + 3599] = lujets_1.p[i__ + 35999];
			xydr_1.rtdr[hjjet4_1.ndr - 1] = (hjcrdn_1.yp[*jp * 3 - 3] +
											 hjcrdn_1.yt[*jt * 3 - 3]) *
											.5f;
			xydr_1.rtdr[hjjet4_1.ndr + 150000] = (hjcrdn_1.yp[*jp * 3 - 2] +
												  hjcrdn_1.yt[*jt * 3 - 2]) *
												 .5f;
		}
		if (lujets_1.k[i__ + 17999] == 3 && (lujets_1.k[i__ + 8999] != 21 ||
											 iinird == 0))
		{
			pxp += lujets_1.p[i__ - 1];
			pyp += lujets_1.p[i__ + 8999];
			pzp += lujets_1.p[i__ + 17999];
			pep += lujets_1.p[i__ + 26999];
		}
		if (lujets_1.k[i__ + 17999] == 4 && (lujets_1.k[i__ + 8999] != 21 ||
											 iinird == 0))
		{
			pxt += lujets_1.p[i__ - 1];
			pyt += lujets_1.p[i__ + 8999];
			pzt += lujets_1.p[i__ + 17999];
			pet += lujets_1.p[i__ + 26999];
		}
		++jpp;
		ip[jpp - 1] = i__;
		ip[jpp + 99] = 0;
		if (lujets_1.k[i__ + 8999] != 21)
		{
			if (lujets_1.k[i__ + 8999] > 0)
			{
				++lpq;
				ipq[lpq - 1] = jpp;
				ip[jpp + 99] = lpq;
			}
			else if (lujets_1.k[i__ + 8999] < 0)
			{
				++lpb;
				ipb[lpb - 1] = jpp;
				ip[jpp + 99] = -lpb;
			}
		}
	}
	if (lpq != lpb)
	{
		++miss;
		if (miss <= 50)
		{
			goto L155;
		}
		s_wsle(&io___340);
		do_lio(&c__3, &c__1, (char *)&lpq, (int)sizeof(int));
		do_lio(&c__3, &c__1, (char *)&lpb, (int)sizeof(int));
		do_lio(&c__9, &c__1, "Q-QBAR NOT CONSERVED OR NOT MATCHED", (int)35);
		e_wsle();
		*jflg = 0;
		return 0;
	}
	j = 0;
L220:
	++j;
	if (j > jpp)
	{
		goto L222;
	}
	if (ip[j + 99] == 0)
	{
		goto L220;
	}
	lp = (i__1 = ip[j + 99], abs(i__1));
	ip1 = ip[j - 1];
	ip2 = ip[j + 99];
	ip[j - 1] = ip[ipq[lp - 1] - 1];
	ip[j + 99] = ip[ipq[lp - 1] + 99];
	ip[ipq[lp - 1] - 1] = ip1;
	ip[ipq[lp - 1] + 99] = ip2;
	if (ip2 > 0)
	{
		ipq[ip2 - 1] = ipq[lp - 1];
	}
	else if (ip2 < 0)
	{
		ipb[-ip2 - 1] = ipq[lp - 1];
	}
	ipq[lp - 1] = j;
	ip1 = ip[j];
	ip2 = ip[j + 100];
	ip[j] = ip[ipb[lp - 1] - 1];
	ip[j + 100] = ip[ipb[lp - 1] + 99];
	ip[ipb[lp - 1] - 1] = ip1;
	ip[ipb[lp - 1] + 99] = ip2;
	if (ip2 > 0)
	{
		ipq[ip2 - 1] = ipb[lp - 1];
	}
	else if (ip2 < 0)
	{
		ipb[-ip2 - 1] = ipb[lp - 1];
	}
	ipb[lp - 1] = j + 1;
	++j;
	goto L220;
L222:
	if (lpq >= 1)
	{
		i__1 = lpq;
		for (l0 = 2; l0 <= i__1; ++l0)
		{
			ip1 = ip[(l0 << 1) - 4];
			ip2 = ip[(l0 << 1) + 96];
			ip[(l0 << 1) - 4] = ip[ipq[l0 - 1] - 1];
			ip[(l0 << 1) + 96] = ip[ipq[l0 - 1] + 99];
			ip[ipq[l0 - 1] - 1] = ip1;
			ip[ipq[l0 - 1] + 99] = ip2;
			if (ip2 > 0)
			{
				ipq[ip2 - 1] = ipq[l0 - 1];
			}
			else if (ip2 < 0)
			{
				ipb[-ip2 - 1] = ipq[l0 - 1];
			}
			ipq[l0 - 1] = (l0 << 1) - 3;

			ip1 = ip[(l0 << 1) - 3];
			ip2 = ip[(l0 << 1) + 97];
			ip[(l0 << 1) - 3] = ip[ipb[l0 - 1] - 1];
			ip[(l0 << 1) + 97] = ip[ipb[l0 - 1] + 99];
			ip[ipb[l0 - 1] - 1] = ip1;
			ip[ipb[l0 - 1] + 99] = ip2;
			if (ip2 > 0)
			{
				ipq[ip2 - 1] = ipb[l0 - 1];
			}
			else if (ip2 < 0)
			{
				ipb[-ip2 - 1] = ipb[l0 - 1];
			}
			ipb[l0 - 1] = (l0 << 1) - 2;
		}
		ip1 = ip[(lpq << 1) - 2];
		ip2 = ip[(lpq << 1) + 98];
		ip[(lpq << 1) - 2] = ip[ipq[0] - 1];
		ip[(lpq << 1) + 98] = ip[ipq[0] + 99];
		ip[ipq[0] - 1] = ip1;
		ip[ipq[0] + 99] = ip2;
		if (ip2 > 0)
		{
			ipq[ip2 - 1] = ipq[0];
		}
		else if (ip2 < 0)
		{
			ipb[-ip2 - 1] = ipq[0];
		}
		ipq[0] = (lpq << 1) - 1;
		ip1 = ip[jpp - 1];
		ip2 = ip[jpp + 99];
		ip[jpp - 1] = ip[ipb[0] - 1];
		ip[jpp + 99] = ip[ipb[0] + 99];
		ip[ipb[0] - 1] = ip1;
		ip[ipb[0] + 99] = ip2;
		if (ip2 > 0)
		{
			ipq[ip2 - 1] = ipb[0];
		}
		else if (ip2 < 0)
		{
			ipb[-ip2 - 1] = ipb[0];
		}
		ipb[0] = jpp;
	}
	if (hjjet2_1.nsg >= mxsg)
	{
		*jflg = 0;
		s_wsle(&io___342);
		do_lio(&c__9, &c__1, "number of jets forming single strings exceeds",
			   (int)45);
		e_wsle();
		s_wsle(&io___343);
		do_lio(&c__9, &c__1, "the common block size", (int)21);
		e_wsle();
		return 0;
	}
	if (jpp > mxsj)
	{
		*jflg = 0;
		s_wsle(&io___344);
		do_lio(&c__9, &c__1, "number of partons per single jet system", (int)39);
		e_wsle();
		s_wsle(&io___345);
		do_lio(&c__9, &c__1, "exceeds the common block size", (int)29);
		e_wsle();
		return 0;
	}
	++hjjet2_1.nsg;
	hjjet2_1.njsg[hjjet2_1.nsg - 1] = jpp;
	hjjet2_1.iasg[hjjet2_1.nsg - 1] = *jp;
	hjjet2_1.iasg[hjjet2_1.nsg + 150000] = *jt;
	hjjet2_1.iasg[hjjet2_1.nsg + 300001] = 0;
	i__1 = jpp;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		hjjet2_1.k1sg[hjjet2_1.nsg + i__ * 150001 - 150002] = 2;
		hjjet2_1.k2sg[hjjet2_1.nsg + i__ * 150001 - 150002] = lujets_1.k[ip[i__ - 1] + 8999];
		if (hjjet2_1.k2sg[hjjet2_1.nsg + i__ * 150001 - 150002] < 0)
		{
			hjjet2_1.k1sg[hjjet2_1.nsg + i__ * 150001 - 150002] = 1;
		}
		hjjet2_1.pxsg[hjjet2_1.nsg + i__ * 150001 - 150002] = lujets_1.p[ip[i__ - 1] - 1];
		hjjet2_1.pysg[hjjet2_1.nsg + i__ * 150001 - 150002] = lujets_1.p[ip[i__ - 1] + 8999];
		hjjet2_1.pzsg[hjjet2_1.nsg + i__ * 150001 - 150002] = lujets_1.p[ip[i__ - 1] + 17999];
		hjjet2_1.pesg[hjjet2_1.nsg + i__ * 150001 - 150002] = lujets_1.p[ip[i__ - 1] + 26999];
		hjjet2_1.pmsg[hjjet2_1.nsg + i__ * 150001 - 150002] = lujets_1.p[ip[i__ - 1] + 35999];
	}
	hjjet2_1.k1sg[hjjet2_1.nsg - 1] = 2;
	hjjet2_1.k1sg[hjjet2_1.nsg + jpp * 150001 - 150002] = 1;
L900:
	hstrng_1.pp[*jp - 1] = pxp;
	hstrng_1.pp[*jp + 299] = pyp;
	hstrng_1.pp[*jp + 599] = pzp;
	hstrng_1.pp[*jp + 899] = pep;
	hstrng_1.pp[*jp + 1199] = 0.f;
	hstrng_1.pt[*jt - 1] = pxt;
	hstrng_1.pt[*jt + 299] = pyt;
	hstrng_1.pt[*jt + 599] = pzt;
	hstrng_1.pt[*jt + 899] = pet;
	hstrng_1.pt[*jt + 1199] = 0.f;
	++hstrng_1.nfp[*jp + 1499];
	++hstrng_1.nft[*jt + 1499];
	return 0;

L1000:
	*jflg = -1;
	if (hparnt_1.ihpr2[9] == 0)
	{
		return 0;
	}
	s_wsle(&io___346);
	do_lio(&c__9, &c__1, "Fatal HIJHRD error", (int)18);
	e_wsle();
	s_wsle(&io___347);
	do_lio(&c__3, &c__1, (char *)&(*jp), (int)sizeof(int));
	do_lio(&c__9, &c__1, " proj E+,E-", (int)11);
	do_lio(&c__4, &c__1, (char *)&epp, (int)sizeof(float));
	do_lio(&c__4, &c__1, (char *)&epm, (int)sizeof(float));
	do_lio(&c__9, &c__1, " status", (int)7);
	do_lio(&c__3, &c__1, (char *)&hstrng_1.nfp[*jp + 1199], (int)sizeof(int));
	e_wsle();
	s_wsle(&io___348);
	do_lio(&c__3, &c__1, (char *)&(*jt), (int)sizeof(int));
	do_lio(&c__9, &c__1, " targ E+,E_", (int)11);
	do_lio(&c__4, &c__1, (char *)&etp, (int)sizeof(float));
	do_lio(&c__4, &c__1, (char *)&etm, (int)sizeof(float));
	do_lio(&c__9, &c__1, " status", (int)7);
	do_lio(&c__3, &c__1, (char *)&hstrng_1.nft[*jt + 1199], (int)sizeof(int));
	e_wsle();
	return 0;
}

int jetini_(int *jp, int *jt, int *itrig)
{
	static int ini[8] = {0, 0, 0, 0, 0, 0, 0, 0};
	static int ilast = -1;

	int i__1;

	int s_copy(char *, char *, int, int);
	int s_wsle(cilist *), do_lio(int *, int *, char *, int),
		e_wsle(void);

	static int i__, j;
	static char beam[16], targ[16];
	static int isel, isub;
	static float coef0[32000], xsec0[1608];
	static int mint44[8], mint45[8], itype;
	extern int pyinit_(char *, char *, char *, float *,
					   int, int, int);

	static cilist io___356 = {0, 6, 0, 0, 0};
	static cilist io___358 = {0, 6, 0, 0, 0};

	hparnt_1.ihnt2[10] = *jp;
	hparnt_1.ihnt2[11] = *jt;
	if (hparnt_1.ihnt2[4] != 0 && hparnt_1.ihnt2[5] != 0)
	{
		itype = 1;
	}
	else if (hparnt_1.ihnt2[4] != 0 && hparnt_1.ihnt2[5] == 0)
	{
		itype = 1;
		if (hstrng_1.nft[*jt + 899] == 2112)
		{
			itype = 2;
		}
	}
	else if (hparnt_1.ihnt2[4] == 0 && hparnt_1.ihnt2[5] != 0)
	{
		itype = 1;
		if (hstrng_1.nfp[*jp + 899] == 2112)
		{
			itype = 2;
		}
	}
	else
	{
		if (hstrng_1.nfp[*jp + 899] == 2212 && hstrng_1.nft[*jt + 899] ==
												   2212)
		{
			itype = 1;
		}
		else if (hstrng_1.nfp[*jp + 899] == 2212 && hstrng_1.nft[*jt + 899] == 2112)
		{
			itype = 2;
		}
		else if (hstrng_1.nfp[*jp + 899] == 2112 && hstrng_1.nft[*jt + 899] == 2212)
		{
			itype = 3;
		}
		else
		{
			itype = 4;
		}
	}

	if (*itrig != 0)
	{
		goto L160;
	}
	if (*itrig == ilast)
	{
		goto L150;
	}
	pypars_1.mstp[1] = 2;
	pypars_1.mstp[32] = 1;
	pypars_1.parp[30] = hparnt_1.hipr1[16];
	pypars_1.mstp[50] = 3;
	pypars_1.mstp[60] = 1;
	pypars_1.mstp[70] = 1;
	if (hparnt_1.ihpr2[1] == 0 || hparnt_1.ihpr2[1] == 2)
	{
		pypars_1.mstp[60] = 0;
	}
	if (hparnt_1.ihpr2[1] == 0 || hparnt_1.ihpr2[1] == 1)
	{
		pypars_1.mstp[70] = 0;
	}

	pypars_1.mstp[80] = 0;
	pypars_1.mstp[81] = 1;
	pypars_1.mstp[110] = 0;
	if (hparnt_1.ihpr2[9] == 0)
	{
		pypars_1.mstp[121] = 0;
	}
	pypars_1.parp[80] = hparnt_1.hipr1[7];
	pysubs_1.ckin[4] = hparnt_1.hipr1[7];
	pysubs_1.ckin[2] = hparnt_1.hipr1[7];
	pysubs_1.ckin[3] = hparnt_1.hipr1[8];
	if (hparnt_1.hipr1[8] <= hparnt_1.hipr1[7])
	{
		pysubs_1.ckin[3] = -1.f;
	}
	pysubs_1.ckin[8] = -10.f;
	pysubs_1.ckin[9] = 10.f;
	pysubs_1.msel = 0;
	for (isub = 1; isub <= 200; ++isub)
	{
		pysubs_1.msub[isub - 1] = 0;
	}
	pysubs_1.msub[10] = 1;
	pysubs_1.msub[11] = 1;
	pysubs_1.msub[12] = 1;
	pysubs_1.msub[27] = 1;
	pysubs_1.msub[52] = 1;
	pysubs_1.msub[67] = 1;
	pysubs_1.msub[80] = 1;
	pysubs_1.msub[81] = 1;
	i__1 = min(8, ludat3_1.mdcy[1020]);
	for (j = 1; j <= i__1; ++j)
	{
		ludat3_1.mdme[ludat3_1.mdcy[520] + j - 2] = 0;
	}
	isel = 4;
	if (hparnt_1.hint1[0] >= 20.f && hparnt_1.ihpr2[17] == 1)
	{
		isel = 5;
	}
	ludat3_1.mdme[ludat3_1.mdcy[520] + isel - 2] = 1;
	pysubs_1.msub[13] = 1;
	pysubs_1.msub[17] = 1;
	pysubs_1.msub[28] = 1;
L150:
	if (ini[itype - 1] != 0)
	{
		goto L800;
	}
	goto L400;

L160:
	itype += 4;
	if (*itrig == ilast)
	{
		goto L260;
	}
	pypars_1.parp[80] = dabs(hparnt_1.hipr1[9]) - .25f;
	pysubs_1.ckin[4] = dabs(hparnt_1.hipr1[9]) - .25f;
	pysubs_1.ckin[2] = dabs(hparnt_1.hipr1[9]) - .25f;
	pysubs_1.ckin[3] = dabs(hparnt_1.hipr1[9]) + .25f;
	if (hparnt_1.hipr1[9] < hparnt_1.hipr1[7])
	{
		pysubs_1.ckin[3] = -1.f;
	}

	pysubs_1.msel = 0;
	for (isub = 1; isub <= 200; ++isub)
	{
		pysubs_1.msub[isub - 1] = 0;
	}
	if (hparnt_1.ihpr2[2] == 1)
	{
		pysubs_1.msub[10] = 1;
		pysubs_1.msub[11] = 1;
		pysubs_1.msub[12] = 1;
		pysubs_1.msub[27] = 1;
		pysubs_1.msub[52] = 1;
		pysubs_1.msub[67] = 1;
		pysubs_1.msub[80] = 1;
		pysubs_1.msub[81] = 1;
		pysubs_1.msub[13] = 1;
		pysubs_1.msub[17] = 1;
		pysubs_1.msub[28] = 1;
		i__1 = min(8, ludat3_1.mdcy[1020]);
		for (j = 1; j <= i__1; ++j)
		{
			ludat3_1.mdme[ludat3_1.mdcy[520] + j - 2] = 0;
		}
		isel = 4;
		if (hparnt_1.hint1[0] >= 20.f && hparnt_1.ihpr2[17] == 1)
		{
			isel = 5;
		}
		ludat3_1.mdme[ludat3_1.mdcy[520] + isel - 2] = 1;
	}
	else if (hparnt_1.ihpr2[2] == 2)
	{
		pysubs_1.msub[13] = 1;
		pysubs_1.msub[17] = 1;
		pysubs_1.msub[28] = 1;
	}
	else if (hparnt_1.ihpr2[2] == 3)
	{
		pysubs_1.ckin[2] = max(0.f, hparnt_1.hipr1[9]);
		pysubs_1.ckin[4] = hparnt_1.hipr1[7];
		pypars_1.parp[80] = hparnt_1.hipr1[7];
		pysubs_1.msub[80] = 1;
		pysubs_1.msub[81] = 1;
		i__1 = min(8, ludat3_1.mdcy[1020]);
		for (j = 1; j <= i__1; ++j)
		{
			ludat3_1.mdme[ludat3_1.mdcy[520] + j - 2] = 0;
		}
		isel = 4;
		if (hparnt_1.hint1[0] >= 20.f && hparnt_1.ihpr2[17] == 1)
		{
			isel = 5;
		}
		ludat3_1.mdme[ludat3_1.mdcy[520] + isel - 2] = 1;
	}
L260:
	if (ini[itype - 1] != 0)
	{
		goto L800;
	}

L400:
	ini[itype - 1] = 1;
	if (hparnt_1.ihpr2[9] == 0)
	{
		pypars_1.mstp[121] = 0;
	}
	if (hstrng_1.nfp[*jp + 899] == 2212)
	{
		s_copy(beam, "P", (int)16, (int)1);
	}
	else if (hstrng_1.nfp[*jp + 899] == -2212)
	{
		s_copy(beam, "P~", (int)16, (int)2);
	}
	else if (hstrng_1.nfp[*jp + 899] == 2112)
	{
		s_copy(beam, "N", (int)16, (int)1);
	}
	else if (hstrng_1.nfp[*jp + 899] == -2112)
	{
		s_copy(beam, "N~", (int)16, (int)2);
	}
	else if (hstrng_1.nfp[*jp + 899] == 211)
	{
		s_copy(beam, "PI+", (int)16, (int)3);
	}
	else if (hstrng_1.nfp[*jp + 899] == -211)
	{
		s_copy(beam, "PI-", (int)16, (int)3);
	}
	else if (hstrng_1.nfp[*jp + 899] == 321)
	{
		s_copy(beam, "PI+", (int)16, (int)3);
	}
	else if (hstrng_1.nfp[*jp + 899] == -321)
	{
		s_copy(beam, "PI-", (int)16, (int)3);
	}
	else
	{
		s_wsle(&io___356);
		do_lio(&c__9, &c__1, "unavailable beam type", (int)21);
		do_lio(&c__3, &c__1, (char *)&hstrng_1.nfp[*jp + 899], (int)sizeof(int));
		e_wsle();
	}
	if (hstrng_1.nft[*jt + 899] == 2212)
	{
		s_copy(targ, "P", (int)16, (int)1);
	}
	else if (hstrng_1.nft[*jt + 899] == -2212)
	{
		s_copy(targ, "P~", (int)16, (int)2);
	}
	else if (hstrng_1.nft[*jt + 899] == 2112)
	{
		s_copy(targ, "N", (int)16, (int)1);
	}
	else if (hstrng_1.nft[*jt + 899] == -2112)
	{
		s_copy(targ, "N~", (int)16, (int)2);
	}
	else if (hstrng_1.nft[*jt + 899] == 211)
	{
		s_copy(targ, "PI+", (int)16, (int)3);
	}
	else if (hstrng_1.nft[*jt + 899] == -211)
	{
		s_copy(targ, "PI-", (int)16, (int)3);
	}
	else if (hstrng_1.nft[*jt + 899] == 321)
	{
		s_copy(targ, "PI+", (int)16, (int)3);
	}
	else if (hstrng_1.nft[*jt + 899] == -321)
	{
		s_copy(targ, "PI-", (int)16, (int)3);
	}
	else
	{
		s_wsle(&io___358);
		do_lio(&c__9, &c__1, "unavailable target type", (int)23);
		do_lio(&c__3, &c__1, (char *)&hstrng_1.nft[*jt + 899], (int)sizeof(int));
		e_wsle();
	}

	hparnt_1.ihnt2[15] = 1;
	pyinit_("CMS", beam, targ, hparnt_1.hint1, (int)3, (int)16, (int)16);
	hpint_1.mint4 = pyint1_1.mint[43];
	hpint_1.mint5 = pyint1_1.mint[44];
	mint44[itype - 1] = pyint1_1.mint[43];
	mint45[itype - 1] = pyint1_1.mint[44];
	hpint_1.atxs[0] = pyint5_1.xsec[0];
	xsec0[itype - 1] = pyint5_1.xsec[0];
	for (i__ = 1; i__ <= 200; ++i__)
	{
		hpint_1.atxs[i__] = pyint5_1.xsec[i__];
		xsec0[itype + (i__ << 3) - 1] = pyint5_1.xsec[i__];
		for (j = 1; j <= 20; ++j)
		{
			hpint_1.atco[i__ + j * 200 - 201] = pyint2_1.coef[i__ + j * 200 -
															  201];
			coef0[itype + (i__ + j * 200 << 3) - 1609] = pyint2_1.coef[i__ +
																	   j * 200 - 201];
		}
	}

	hparnt_1.ihnt2[15] = 0;

	return 0;

L800:
	pyint1_1.mint[43] = mint44[itype - 1];
	pyint1_1.mint[44] = mint45[itype - 1];
	hpint_1.mint4 = pyint1_1.mint[43];
	hpint_1.mint5 = pyint1_1.mint[44];
	pyint5_1.xsec[0] = xsec0[itype - 1];
	hpint_1.atxs[0] = pyint5_1.xsec[0];
	for (i__ = 1; i__ <= 200; ++i__)
	{
		pyint5_1.xsec[i__] = xsec0[itype + (i__ << 3) - 1];
		hpint_1.atxs[i__] = pyint5_1.xsec[i__];
		for (j = 1; j <= 20; ++j)
		{
			pyint2_1.coef[i__ + j * 200 - 201] = coef0[itype + (i__ + j * 200 << 3) - 1609];
			hpint_1.atco[i__ + j * 200 - 201] = pyint2_1.coef[i__ + j * 200 -
															  201];
		}
	}
	ilast = *itrig;
	pyint1_1.mint[10] = hstrng_1.nfp[*jp + 899];
	pyint1_1.mint[11] = hstrng_1.nft[*jt + 899];
	return 0;
}

int hijini_(void)
{
	int i__1, i__2;
	float r__1, r__2;

	double sqrt(double);

	static int i__, idq, ipp, ipt, idqq;
	extern double ranart_(int *);
	extern int attflv_(int *, int *, int *);
	extern double ulmass_(int *);

	hjjet2_1.nsg = 0;
	hjjet4_1.ndr = 0;
	ipp = 2212;
	ipt = 2212;
	if (hparnt_1.ihnt2[4] != 0)
	{
		ipp = hparnt_1.ihnt2[4];
	}
	if (hparnt_1.ihnt2[5] != 0)
	{
		ipt = hparnt_1.ihnt2[5];
	}
	i__1 = hparnt_1.ihnt2[0];
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		hstrng_1.pp[i__ - 1] = 0.f;
		hstrng_1.pp[i__ + 299] = 0.f;
		r__1 = hparnt_1.hint1[0];
		r__2 = hparnt_1.hint1[7];
		hstrng_1.pp[i__ + 599] = sqrt(r__1 * r__1 / 4.f - r__2 * r__2);
		hstrng_1.pp[i__ + 899] = hparnt_1.hint1[0] / 2;
		hstrng_1.pp[i__ + 1199] = hparnt_1.hint1[7];
		hstrng_1.pp[i__ + 1499] = 0.f;
		hstrng_1.pp[i__ + 1799] = 0.f;
		hstrng_1.pp[i__ + 2099] = 0.f;
		hstrng_1.pp[i__ + 2399] = 0.f;
		hstrng_1.pp[i__ + 2699] = 0.f;
		hstrng_1.pp[i__ + 2999] = 0.f;
		hstrng_1.pp[i__ + 3299] = 0.f;
		hstrng_1.nfp[i__ + 599] = ipp;
		hstrng_1.nfp[i__ + 899] = ipp;
		hstrng_1.nfp[i__ + 1199] = 0;
		hstrng_1.nfp[i__ + 1499] = 0;
		hstrng_1.nfp[i__ + 1799] = 0;
		hstrng_1.nfp[i__ + 2099] = 0;
		hstrng_1.nfp[i__ + 2399] = 0;
		hstrng_1.nfp[i__ + 2699] = 0;
		hstrng_1.nfp[i__ + 2999] = 0;
		hjjet1_1.npj[i__ - 1] = 0;
		if (i__ > abs(hparnt_1.ihnt2[1]))
		{
			hstrng_1.nfp[i__ + 599] = 2112;
		}
		attflv_(&hstrng_1.nfp[i__ + 599], &idq, &idqq);
		hstrng_1.nfp[i__ - 1] = idq;
		hstrng_1.nfp[i__ + 299] = idqq;
		hstrng_1.nfp[i__ + 4199] = -1;
		if (abs(idq) > 1000 || (i__2 = idq * idqq, abs(i__2)) < 100 &&
								   ranart_(&rndf77_1.nseed) < .5f)
		{
			hstrng_1.nfp[i__ + 4199] = 1;
		}
		hstrng_1.pp[i__ + 3899] = ulmass_(&idq);
		hstrng_1.pp[i__ + 4199] = ulmass_(&idqq);
	}

	i__1 = hparnt_1.ihnt2[2];
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		hstrng_1.pt[i__ - 1] = 0.f;
		hstrng_1.pt[i__ + 299] = 0.f;
		r__1 = hparnt_1.hint1[0];
		r__2 = hparnt_1.hint1[8];
		hstrng_1.pt[i__ + 599] = -sqrt(r__1 * r__1 / 4.f - r__2 * r__2);
		hstrng_1.pt[i__ + 899] = hparnt_1.hint1[0] / 2.f;
		hstrng_1.pt[i__ + 1199] = hparnt_1.hint1[8];
		hstrng_1.pt[i__ + 1499] = 0.f;
		hstrng_1.pt[i__ + 1799] = 0.f;
		hstrng_1.pt[i__ + 2099] = 0.f;
		hstrng_1.pt[i__ + 2399] = 0.f;
		hstrng_1.pt[i__ + 2699] = 0.f;
		hstrng_1.pt[i__ + 2999] = 0.f;
		hstrng_1.pt[i__ + 3299] = 0.f;
		hstrng_1.nft[i__ + 599] = ipt;
		hstrng_1.nft[i__ + 899] = ipt;
		hstrng_1.nft[i__ + 1199] = 0;
		hstrng_1.nft[i__ + 1499] = 0;
		hstrng_1.nft[i__ + 1799] = 0;
		hstrng_1.nft[i__ + 2099] = 0;
		hstrng_1.nft[i__ + 2399] = 0;
		hstrng_1.nft[i__ + 2699] = 0;
		hstrng_1.nft[i__ + 2999] = 0;
		hjjet1_1.ntj[i__ - 1] = 0;
		if (i__ > abs(hparnt_1.ihnt2[3]))
		{
			hstrng_1.nft[i__ + 599] = 2112;
		}
		attflv_(&hstrng_1.nft[i__ + 599], &idq, &idqq);
		hstrng_1.nft[i__ - 1] = idq;
		hstrng_1.nft[i__ + 299] = idqq;
		hstrng_1.nft[i__ + 4199] = 1;
		if (abs(idq) > 1000 || (i__2 = idq * idqq, abs(i__2)) < 100 &&
								   ranart_(&rndf77_1.nseed) < .5f)
		{
			hstrng_1.nft[i__ + 4199] = -1;
		}
		hstrng_1.pt[i__ + 3899] = ulmass_(&idq);
		hstrng_1.pt[i__ + 4199] = ulmass_(&idqq);
	}
	return 0;
}

int attflv_(int *id, int *idq, int *idqq)
{
	static float x;
	static int id0, id00, nsign;
	extern double ranart_(int *);

	if (abs(*id) < 100)
	{
		nsign = 1;
		*idq = *id / 100;
		*idqq = -(*id) / 10 + *idq * 10;
		if (abs(*idq) == 3)
		{
			nsign = -1;
		}
		*idq = nsign * *idq;
		*idqq = nsign * *idqq;
		if (*idq < 0)
		{
			id0 = *idq;
			*idq = *idqq;
			*idqq = id0;
		}
		return 0;
	}
	*idq = 2;
	if (abs(*id) == 2112)
	{
		*idq = 1;
	}
	*idqq = 2101;
	x = ranart_(&rndf77_1.nseed);
	if (x <= .5f)
	{
		goto L30;
	}
	if (x > .666667f)
	{
		goto L10;
	}
	*idqq = 2103;
	goto L30;
L10:
	*idq = 1;
	*idqq = 2203;
	if (abs(*id) == 2112)
	{
		*idq = 2;
		*idqq = 1103;
	}
L30:
	if (*id < 0)
	{
		id00 = *idqq;
		*idqq = -(*idq);
		*idq = -id00;
	}
	return 0;
}

int hijcsc_(int *jp, int *jt)
{
	int i__1, i__2, i__3;
	float r__1, r__2, r__3, r__4;

	double sqrt(double), exp(double);

	static int i__, k;
	static float r2, bb, bx, by, bz, dx, dy, dz, gs, gs0, dis, dpp1, dpp2,
		psc1[5], psc2[5], dpt1, dpt2, pabs;
	extern double romg_(float *);
	extern int hijels_(float *, float *);
	extern double ranart_(int *);

	if (*jp == 0 || *jt == 0)
	{
		goto L25;
	}
	for (i__ = 1; i__ <= 5; ++i__)
	{
		psc1[i__ - 1] = hstrng_1.pp[*jp + i__ * 300 - 301];
		psc2[i__ - 1] = hstrng_1.pt[*jt + i__ * 300 - 301];
	}
	hijels_(psc1, psc2);
	dpp1 = psc1[0] - hstrng_1.pp[*jp - 1];
	dpp2 = psc1[1] - hstrng_1.pp[*jp + 299];
	dpt1 = psc2[0] - hstrng_1.pt[*jt - 1];
	dpt2 = psc2[1] - hstrng_1.pt[*jt + 299];
	hstrng_1.pp[*jp + 1499] += dpp1 / 2.f;
	hstrng_1.pp[*jp + 1799] += dpp2 / 2.f;
	hstrng_1.pp[*jp + 2099] += dpp1 / 2.f;
	hstrng_1.pp[*jp + 2399] += dpp2 / 2.f;
	hstrng_1.pt[*jt + 1499] += dpt1 / 2.f;
	hstrng_1.pt[*jt + 1799] += dpt2 / 2.f;
	hstrng_1.pt[*jt + 2099] += dpt1 / 2.f;
	hstrng_1.pt[*jt + 2399] += dpt2 / 2.f;
	for (i__ = 1; i__ <= 4; ++i__)
	{
		hstrng_1.pp[*jp + i__ * 300 - 301] = psc1[i__ - 1];
		hstrng_1.pt[*jt + i__ * 300 - 301] = psc2[i__ - 1];
	}
	i__1 = 1, i__2 = hstrng_1.nfp[*jp + 1199];
	hstrng_1.nfp[*jp + 1199] = max(i__1, i__2);
	i__1 = 1, i__2 = hstrng_1.nft[*jt + 1199];
	hstrng_1.nft[*jt + 1199] = max(i__1, i__2);
	return 0;
L25:
	if (*jp == 0)
	{
		goto L45;
	}
	r__1 = hstrng_1.pp[*jp - 1];
	r__2 = hstrng_1.pp[*jp + 299];
	r__3 = hstrng_1.pp[*jp + 599];
	pabs = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
	bx = hstrng_1.pp[*jp - 1] / pabs;
	by = hstrng_1.pp[*jp + 299] / pabs;
	bz = hstrng_1.pp[*jp + 599] / pabs;
	i__1 = hparnt_1.ihnt2[0];
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		dx = hjcrdn_1.yp[i__ * 3 - 3] - hjcrdn_1.yp[*jp * 3 - 3];
		dy = hjcrdn_1.yp[i__ * 3 - 2] - hjcrdn_1.yp[*jp * 3 - 2];
		dz = hjcrdn_1.yp[i__ * 3 - 1] - hjcrdn_1.yp[*jp * 3 - 1];
		dis = dx * bx + dy * by + dz * bz;
		r__1 = dx;
		r__2 = dy;
		r__3 = dz;
		r__4 = dis;
		bb = r__1 * r__1 + r__2 * r__2 + r__3 * r__3 - r__4 * r__4;
		r2 = bb * hparnt_1.hipr1[39] / hparnt_1.hipr1[30] / .1f;
		r__1 = exp(-(hparnt_1.hipr1[29] + hparnt_1.hint1[10]) /
				   hparnt_1.hipr1[30] / 2.f * romg_(&r2));
		gs = 1.f - r__1 * r__1;
		r__1 = exp(-(hparnt_1.hipr1[29] + hparnt_1.hint1[10]) /
				   hparnt_1.hipr1[30] / 2.f * romg_(&c_b26));
		gs0 = 1.f - r__1 * r__1;
		for (k = 1; k <= 5; ++k)
		{
			psc1[k - 1] = hstrng_1.pp[*jp + k * 300 - 301];
			psc2[k - 1] = hstrng_1.pp[i__ + k * 300 - 301];
		}
		hijels_(psc1, psc2);
		dpp1 = psc1[0] - hstrng_1.pp[*jp - 1];
		dpp2 = psc1[1] - hstrng_1.pp[*jp + 299];
		dpt1 = psc2[0] - hstrng_1.pp[i__ - 1];
		dpt2 = psc2[1] - hstrng_1.pp[i__ + 299];
		hstrng_1.pp[*jp + 1499] += dpp1 / 2.f;
		hstrng_1.pp[*jp + 1799] += dpp2 / 2.f;
		hstrng_1.pp[*jp + 2099] += dpp1 / 2.f;
		hstrng_1.pp[*jp + 2399] += dpp2 / 2.f;
		hstrng_1.pp[i__ + 1499] += dpt1 / 2.f;
		hstrng_1.pp[i__ + 1799] += dpt2 / 2.f;
		hstrng_1.pp[i__ + 2099] += dpt1 / 2.f;
		hstrng_1.pp[i__ + 2399] += dpt2 / 2.f;
		for (k = 1; k <= 5; ++k)
		{
			hstrng_1.pp[*jp + k * 300 - 301] = psc1[k - 1];
			hstrng_1.pp[i__ + k * 300 - 301] = psc2[k - 1];
		}
		i__2 = 1, i__3 = hstrng_1.nfp[i__ + 1199];
		hstrng_1.nfp[i__ + 1199] = max(i__2, i__3);
		goto L45;
	}
L45:
	if (*jt == 0)
	{
		return 0;
	}
	r__1 = hstrng_1.pt[*jt - 1];
	r__2 = hstrng_1.pt[*jt + 299];
	r__3 = hstrng_1.pt[*jt + 599];
	pabs = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
	bx = hstrng_1.pt[*jt - 1] / pabs;
	by = hstrng_1.pt[*jt + 299] / pabs;
	bz = hstrng_1.pt[*jt + 599] / pabs;
	i__1 = hparnt_1.ihnt2[2];
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		dx = hjcrdn_1.yt[i__ * 3 - 3] - hjcrdn_1.yt[*jt * 3 - 3];
		dy = hjcrdn_1.yt[i__ * 3 - 2] - hjcrdn_1.yt[*jt * 3 - 2];
		dz = hjcrdn_1.yt[i__ * 3 - 1] - hjcrdn_1.yt[*jt * 3 - 1];
		dis = dx * bx + dy * by + dz * bz;
		r__1 = dx;
		r__2 = dy;
		r__3 = dz;
		r__4 = dis;
		bb = r__1 * r__1 + r__2 * r__2 + r__3 * r__3 - r__4 * r__4;
		r2 = bb * hparnt_1.hipr1[39] / hparnt_1.hipr1[30] / .1f;
		r__1 = 1.f - exp(-(hparnt_1.hipr1[29] + hparnt_1.hint1[10]) /
						 hparnt_1.hipr1[30] / 2.f * romg_(&r2));
		gs = r__1 * r__1;
		r__1 = 1.f - exp(-(hparnt_1.hipr1[29] + hparnt_1.hint1[10]) /
						 hparnt_1.hipr1[30] / 2.f * romg_(&c_b26));
		gs0 = r__1 * r__1;
		for (k = 1; k <= 5; ++k)
		{
			psc1[k - 1] = hstrng_1.pt[*jt + k * 300 - 301];
			psc2[k - 1] = hstrng_1.pt[i__ + k * 300 - 301];
		}
		hijels_(psc1, psc2);
		dpp1 = psc1[0] - hstrng_1.pt[*jt - 1];
		dpp2 = psc1[1] - hstrng_1.pt[*jt + 299];
		dpt1 = psc2[0] - hstrng_1.pt[i__ - 1];
		dpt2 = psc2[1] - hstrng_1.pt[i__ + 299];
		hstrng_1.pt[*jt + 1499] += dpp1 / 2.f;
		hstrng_1.pt[*jt + 1799] += dpp2 / 2.f;
		hstrng_1.pt[*jt + 2099] += dpp1 / 2.f;
		hstrng_1.pt[*jt + 2399] += dpp2 / 2.f;
		hstrng_1.pt[i__ + 1499] += dpt1 / 2.f;
		hstrng_1.pt[i__ + 1799] += dpt2 / 2.f;
		hstrng_1.pt[i__ + 2099] += dpt1 / 2.f;
		hstrng_1.pt[i__ + 2399] += dpt2 / 2.f;
		for (k = 1; k <= 5; ++k)
		{
			hstrng_1.pt[*jt + k * 300 - 301] = psc1[k - 1];
			hstrng_1.pt[i__ + k * 300 - 301] = psc2[k - 1];
		}
		i__2 = 1, i__3 = hstrng_1.nft[i__ + 1199];
		hstrng_1.nft[i__ + 1199] = max(i__2, i__3);
		return 0;
	}
}

int hijels_(float *psc1, float *psc2)
{
	float r__1, r__2, r__3, r__4;
	double d__1, d__2, d__3;

	double sqrt(double), exp(double);
	int s_wsle(cilist *), do_lio(int *, int *, char *, int),
		e_wsle(void);
	double sin(double), cos(double);

	static float bb, cc;
	static double db;
	static float ep, rr, tt, am1, am2;
	static double dp1, dp2, dp3, dp4, dga;
	static float ecm;
	static double dbp;
	static float amm;
	static double dbx, dby, dbz;
	static float phi, els, pcm1, pcm2, pcm3, els0, pmax;
	static double dgabp;
	extern double ranart_(int *);

	static cilist io___413 = {0, 6, 0, 0, 0};

	--psc2;
	--psc1;

	cc = 1.f - hparnt_1.hint1[11] / hparnt_1.hint1[12];
	rr = (1.f - cc) * hparnt_1.hint1[12] / hparnt_1.hint1[11] / (1.f - hparnt_1.hipr1[32]) - 1.f;
	r__1 = rr;
	bb = (rr + 3.f + sqrt(rr * 10.f + 9.f + r__1 * r__1)) * .5f;
	r__1 = psc1[1] - psc2[1];
	r__2 = psc1[2] - psc2[2];
	r__3 = psc1[3] - psc2[3];
	ep = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
	if (ep <= .1f)
	{
		return 0;
	}
	r__1 = rr + 1.f;
	els0 = 98.f / ep + r__1 * r__1 * 52.f;
	pcm1 = psc1[1] + psc2[1];
	pcm2 = psc1[2] + psc2[2];
	pcm3 = psc1[3] + psc2[3];
	ecm = psc1[4] + psc2[4];
	r__1 = psc1[5];
	am1 = r__1 * r__1;
	r__1 = psc2[5];
	am2 = r__1 * r__1;
	r__1 = ecm;
	r__2 = pcm1;
	r__3 = pcm2;
	r__4 = pcm3;
	amm = r__1 * r__1 - r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
	if (amm <= psc1[5] + psc2[5])
	{
		return 0;
	}
	r__1 = amm;
	r__2 = am1;
	r__3 = am2;
	pmax = (r__1 * r__1 + r__2 * r__2 + r__3 * r__3 - amm * 2.f * am1 - amm * 2.f * am2 - am1 * 2.f * am2) / 4.f / amm;
	pmax = dabs(pmax);
L20:
	tt = ranart_(&rndf77_1.nseed) * min(pmax, 1.5f);
	r__1 = rr * exp((bb - 1.f) * -4.6f * tt) + 1.f;
	els = exp(tt * -2.8f) * 98.f / ep + exp(tt * -9.2f) * 52.f * (r__1 * r__1);
	if (ranart_(&rndf77_1.nseed) > els / els0)
	{
		goto L20;
	}
	phi = hparnt_1.hipr1[39] * 2.f * ranart_(&rndf77_1.nseed);

	dbx = (double)(pcm1 / ecm);
	dby = (double)(pcm2 / ecm);
	dbz = (double)(pcm3 / ecm);
	d__1 = dbx;
	d__2 = dby;
	d__3 = dbz;
	db = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
	if (db > .99999999)
	{
		dbx *= .99999999 / db;
		dby *= .99999999 / db;
		dbz *= .99999999 / db;
		db = .99999999;
		s_wsle(&io___413);
		do_lio(&c__9, &c__1, " (HIJELS) boost vector too large", (int)32);
		e_wsle();
	}
	d__1 = db;
	dga = 1. / sqrt(1. - d__1 * d__1);

	dp1 = (double)(sqrt(tt) * sin(phi));
	dp2 = (double)(sqrt(tt) * cos(phi));
	dp3 = (double)sqrt(pmax - tt);
	dp4 = (double)sqrt(pmax + am1);
	dbp = dbx * dp1 + dby * dp2 + dbz * dp3;
	dgabp = dga * (dga * dbp / (dga + 1.) + dp4);
	psc1[1] = (float)(dp1 + dgabp * dbx);
	psc1[2] = (float)(dp2 + dgabp * dby);
	psc1[3] = (float)(dp3 + dgabp * dbz);
	psc1[4] = (float)(dga * (dp4 + dbp));

	dp1 = -((double)(sqrt(tt) * sin(phi)));
	dp2 = -((double)(sqrt(tt) * cos(phi)));
	dp3 = -((double)sqrt(pmax - tt));
	dp4 = (double)sqrt(pmax + am2);
	dbp = dbx * dp1 + dby * dp2 + dbz * dp3;
	dgabp = dga * (dga * dbp / (dga + 1.) + dp4);
	psc2[1] = (float)(dp1 + dgabp * dbx);
	psc2[2] = (float)(dp2 + dgabp * dby);
	psc2[3] = (float)(dp3 + dgabp * dbz);
	psc2[4] = (float)(dga * (dp4 + dbp));
	return 0;
}

int hijsft_(int *jp, int *jt, int *jout, int *ierror)
{
	int i__1, i__2;
	float r__1, r__2, r__3, r__4, r__5, r__6;
	double d__1;

	double sqrt(double), cos(double), sin(double), pow_dd(double *, double *), exp(double), log(double);
	int s_wsle(cilist *), do_lio(int *, int *, char *, int),
		e_wsle(void);

	static float d1, d2, r1, r2, x1, x2, bb, cc, dd, bx, by, bb1, bb2, dd1,
		dd2, dd3, dd4, dp1, dp2, xp0, yp0, xt0, yt0, xx1, xx2, dpd, sdd,
		dtd, amq, pkc, phi, epp, epm, etp, etm, amx, snn, dpn, dtn, dpx,
		dtx;
	static int jsb;
	static float psb;
	static int nsb;
	static float xxp, sxx, xxt, dpe1, dpe2, pkc1, pkc2, dpm0, phi1, phi2, dtm0,
		phi0;
	static int nfp3, nfp5, nft3;
	static float dpx1, dpy1, dpz1, dpx2, dpy2, dpz2;
	static int nft5;
	static float pkc11, pkc12, pkc21, pkc22, ampd, amtd;
	static int isng;
	static float ptp02, ampx, ptt02;
	static int miss;
	static float amtx, xmin, xmax, xmin1, xmin2, xmax1;
	static int miss4;
	static float xmax2, dpkc11, dpkc12, dpkc21, dpkc22;
	static int kcdip, kcdit;
	static float cthep, cthet, ppjet, pkcmx, ptjet, spntd, spdtn, swptd, spdtx,
		spxtd, swptn, spntx, spxtn, swptx;
	extern double hirnd2_(int *, float *, float *), ulangl_(float *, float *), ranart_(int *);
	static float xminhi, xmaxhi, epmprm, eppprm, etmprm;
	extern double ulmass_(int *);
	static float etpprm;

	static cilist io___490 = {0, 6, 0, 0, 0};
	static cilist io___513 = {0, 6, 0, 0, 0};
	static cilist io___534 = {0, 6, 0, 0, 0};
	static cilist io___535 = {0, 6, 0, 0, 0};
	static cilist io___536 = {0, 6, 0, 0, 0};
	static cilist io___537 = {0, 6, 0, 0, 0};
	static cilist io___538 = {0, 6, 0, 0, 0};
	static cilist io___539 = {0, 6, 0, 0, 0};
	static cilist io___540 = {0, 6, 0, 0, 0};
	static cilist io___541 = {0, 6, 0, 0, 0};
	static cilist io___542 = {0, 6, 0, 0, 0};
	static cilist io___543 = {0, 6, 0, 0, 0};
	static cilist io___544 = {0, 6, 0, 0, 0};
	static cilist io___545 = {0, 6, 0, 0, 0};
	static cilist io___546 = {0, 6, 0, 0, 0};
	static cilist io___547 = {0, 6, 0, 0, 0};
	static cilist io___548 = {0, 6, 0, 0, 0};
	static cilist io___549 = {0, 6, 0, 0, 0};
	static cilist io___550 = {0, 6, 0, 0, 0};
	static cilist io___551 = {0, 6, 0, 0, 0};
	static cilist io___552 = {0, 6, 0, 0, 0};
	static cilist io___553 = {0, 6, 0, 0, 0};
	static cilist io___554 = {0, 6, 0, 0, 0};
	static cilist io___555 = {0, 6, 0, 0, 0};
	static cilist io___556 = {0, 6, 0, 0, 0};
	static cilist io___557 = {0, 6, 0, 0, 0};
	static cilist io___558 = {0, 6, 0, 0, 0};
	static cilist io___559 = {0, 6, 0, 0, 0};
	static cilist io___560 = {0, 6, 0, 0, 0};

	*ierror = 0;
	dpmcm1_1.jjp = *jp;
	dpmcm1_1.jjt = *jt;
	dpmcm2_1.ndpm = 0;
	if (*jp > hparnt_1.ihnt2[0] || *jt > hparnt_1.ihnt2[2])
	{
		return 0;
	}
	epp = hstrng_1.pp[*jp + 899] + hstrng_1.pp[*jp + 599];
	epm = hstrng_1.pp[*jp + 899] - hstrng_1.pp[*jp + 599];
	etp = hstrng_1.pt[*jt + 899] + hstrng_1.pt[*jt + 599];
	etm = hstrng_1.pt[*jt + 899] - hstrng_1.pt[*jt + 599];
	dpmcm1_1.wp = epp + etp;
	dpmcm1_1.wm = epm + etm;
	dpmcm1_1.sw = dpmcm1_1.wp * dpmcm1_1.wm;
	if (dpmcm1_1.wp < 0.f || dpmcm1_1.wm < 0.f)
	{
		goto L1000;
	}
	if (*jout == 0)
	{
		if (epp < 0.f)
		{
			goto L1000;
		}
		if (epm < 0.f)
		{
			goto L1000;
		}
		if (etp < 0.f)
		{
			goto L1000;
		}
		if (etm < 0.f)
		{
			goto L1000;
		}
		if (epp / (epm + .01f) <= etp / (etm + .01f))
		{
			return 0;
		}
	}
	hparnt_1.ihnt2[10] = *jp;
	hparnt_1.ihnt2[11] = *jt;

	miss = 0;
	pkc1 = 0.f;
	pkc2 = 0.f;
	pkc11 = 0.f;
	pkc12 = 0.f;
	pkc21 = 0.f;
	pkc22 = 0.f;
	dpkc11 = 0.f;
	dpkc12 = 0.f;
	dpkc21 = 0.f;
	dpkc22 = 0.f;
	if (hstrng_1.nfp[*jp + 2699] == 1 || hstrng_1.nft[*jt + 2699] == 1)
	{
		if (hstrng_1.nfp[*jp + 2699] == 1)
		{
			phi1 = ulangl_(&hstrng_1.pp[*jp + 2699], &hstrng_1.pp[*jp + 2999]);
			r__1 = hstrng_1.pp[*jp + 2699];
			r__2 = hstrng_1.pp[*jp + 2999];
			ppjet = sqrt(r__1 * r__1 + r__2 * r__2);
			pkc1 = ppjet;
			pkc11 = hstrng_1.pp[*jp + 2699];
			pkc12 = hstrng_1.pp[*jp + 2999];
		}
		if (hstrng_1.nft[*jt + 2699] == 1)
		{
			phi2 = ulangl_(&hstrng_1.pt[*jt + 2699], &hstrng_1.pt[*jt + 2999]);
			r__1 = hstrng_1.pt[*jt + 2699];
			r__2 = hstrng_1.pt[*jt + 2999];
			ptjet = sqrt(r__1 * r__1 + r__2 * r__2);
			pkc2 = ptjet;
			pkc21 = hstrng_1.pt[*jt + 2699];
			pkc22 = hstrng_1.pt[*jt + 2999];
		}
		if (hparnt_1.ihpr2[3] > 0 && hparnt_1.ihnt2[0] > 1 && hparnt_1.ihnt2[2] > 1)
		{
			if (hstrng_1.nfp[*jp + 2699] == 0)
			{
				phi = -phi2;
			}
			else if (hstrng_1.nft[*jt + 2699] == 0)
			{
				phi = phi1;
			}
			else
			{
				phi = (phi1 + phi2 - hparnt_1.hipr1[39]) / 2.f;
			}
			bx = hparnt_1.hint1[18] * cos(hparnt_1.hint1[19]);
			by = hparnt_1.hint1[18] * sin(hparnt_1.hint1[19]);
			xp0 = hjcrdn_1.yp[*jp * 3 - 3];
			yp0 = hjcrdn_1.yp[*jp * 3 - 2];
			xt0 = hjcrdn_1.yt[*jt * 3 - 3] + bx;
			yt0 = hjcrdn_1.yt[*jt * 3 - 2] + by;
			d__1 = (double)hparnt_1.ihnt2[0];
			r__3 = xp0;
			r__4 = yp0;
			r__1 = pow_dd(&d__1, &c_b23) * 1.2f, r__2 = sqrt(r__3 * r__3 +
															 r__4 * r__4);
			r1 = max(r__1, r__2);
			d__1 = (double)hparnt_1.ihnt2[2];
			r__3 = xt0 - bx;
			r__4 = yt0 - by;
			r__1 = pow_dd(&d__1, &c_b23) * 1.2f, r__2 = sqrt(r__3 * r__3 +
															 r__4 * r__4);
			r2 = max(r__1, r__2);
			if ((r__1 = cos(phi), dabs(r__1)) < 1e-5f)
			{
				dd1 = r1;
				dd2 = r1;
				r__2 = r2;
				r__3 = xp0 - bx;
				dd3 = (r__1 = by + sqrt(r__2 * r__2 - r__3 * r__3) - yp0,
					   dabs(r__1));
				r__2 = r2;
				r__3 = xp0 - bx;
				dd4 = (r__1 = by - sqrt(r__2 * r__2 - r__3 * r__3) - yp0,
					   dabs(r__1));
				goto L5;
			}
			bb = sin(phi) * 2.f * (cos(phi) * yp0 - sin(phi) * xp0);
			r__1 = yp0;
			r__2 = r1;
			r__3 = cos(phi);
			cc = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3) + xp0 * sin(phi) * (xp0 * sin(phi) - yp0 * 2.f * cos(phi));
			r__1 = bb;
			dd = r__1 * r__1 - cc * 4.f;
			if (dd < 0.f)
			{
				goto L10;
			}
			xx1 = (-bb + sqrt(dd)) / 2.f;
			xx2 = (-bb - sqrt(dd)) / 2.f;
			dd1 = (r__1 = (xx1 - xp0) / cos(phi), dabs(r__1));
			dd2 = (r__1 = (xx2 - xp0) / cos(phi), dabs(r__1));

			bb = sin(phi) * 2.f * (cos(phi) * (yt0 - by) - sin(phi) * xt0) -
				 bx * 2.f;
			r__1 = bx;
			r__2 = yt0 - by;
			r__3 = r2;
			r__4 = cos(phi);
			cc = (r__1 * r__1 + r__2 * r__2 - r__3 * r__3) * (r__4 * r__4) +
				 xt0 * sin(phi) * (xt0 * sin(phi) - cos(phi) * 2.f * (yt0 - by)) - bx * 2.f * sin(phi) * (cos(phi) * (yt0 - by) - sin(phi) * xt0);
			r__1 = bb;
			dd = r__1 * r__1 - cc * 4.f;
			if (dd < 0.f)
			{
				goto L10;
			}
			xx1 = (-bb + sqrt(dd)) / 2.f;
			xx2 = (-bb - sqrt(dd)) / 2.f;
			dd3 = (r__1 = (xx1 - xt0) / cos(phi), dabs(r__1));
			dd4 = (r__1 = (xx2 - xt0) / cos(phi), dabs(r__1));

		L5:
			dd1 = min(dd1, dd3);
			dd2 = min(dd2, dd4);
			if (dd1 < hparnt_1.hipr1[12])
			{
				dd1 = 0.f;
			}
			if (dd2 < hparnt_1.hipr1[12])
			{
				dd2 = 0.f;
			}
			if (hstrng_1.nfp[*jp + 2699] == 1 && ppjet > hparnt_1.hipr1[10])
			{
				dp1 = dd1 * hparnt_1.hipr1[13] / 2.f;
				r__1 = dp1, r__2 = ppjet - hparnt_1.hipr1[10];
				dp1 = min(r__1, r__2);
				pkc1 = ppjet - dp1;
				dpx1 = cos(phi1) * dp1;
				dpy1 = sin(phi1) * dp1;
				pkc11 = hstrng_1.pp[*jp + 2699] - dpx1;
				pkc12 = hstrng_1.pp[*jp + 2999] - dpy1;
				if (dp1 > 0.f)
				{
					r__1 = hstrng_1.pp[*jp + 3299];
					r__2 = ppjet;
					cthep = hstrng_1.pp[*jp + 3299] / sqrt(r__1 * r__1 + r__2 * r__2);
					r__1 = cthep;
					dpz1 = dp1 * cthep / sqrt(1.f - r__1 * r__1);
					r__1 = dpx1;
					r__2 = dpy1;
					r__3 = dpz1;
					dpe1 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
					eppprm = hstrng_1.pp[*jp + 899] + hstrng_1.pp[*jp + 599] - dpe1 - dpz1;
					epmprm = hstrng_1.pp[*jp + 899] - hstrng_1.pp[*jp + 599] - dpe1 + dpz1;
					if (eppprm <= 0.f || epmprm <= 0.f)
					{
						goto L15;
					}
					epp = eppprm;
					epm = epmprm;
					hstrng_1.pp[*jp + 2699] = pkc11;
					hstrng_1.pp[*jp + 2999] = pkc12;
					++hjjet1_1.npj[*jp - 1];
					hjjet1_1.kfpj[*jp + hjjet1_1.npj[*jp - 1] * 300 - 301] =
						21;
					hjjet1_1.pjpx[*jp + hjjet1_1.npj[*jp - 1] * 300 - 301] =
						dpx1;
					hjjet1_1.pjpy[*jp + hjjet1_1.npj[*jp - 1] * 300 - 301] =
						dpy1;
					hjjet1_1.pjpz[*jp + hjjet1_1.npj[*jp - 1] * 300 - 301] =
						dpz1;
					hjjet1_1.pjpe[*jp + hjjet1_1.npj[*jp - 1] * 300 - 301] =
						dpe1;
					hjjet1_1.pjpm[*jp + hjjet1_1.npj[*jp - 1] * 300 - 301] =
						0.f;
					hstrng_1.pp[*jp + 599] -= dpz1;
					hstrng_1.pp[*jp + 899] -= dpe1;
				}
			}
		L15:
			if (hstrng_1.nft[*jt + 2699] == 1 && ptjet > hparnt_1.hipr1[10])
			{
				dp2 = dd2 * hparnt_1.hipr1[13] / 2.f;
				r__1 = dp2, r__2 = ptjet - hparnt_1.hipr1[10];
				dp2 = min(r__1, r__2);
				pkc2 = ptjet - dp2;
				dpx2 = cos(phi2) * dp2;
				dpy2 = sin(phi2) * dp2;
				pkc21 = hstrng_1.pt[*jt + 2699] - dpx2;
				pkc22 = hstrng_1.pt[*jt + 2999] - dpy2;
				if (dp2 > 0.f)
				{
					r__1 = hstrng_1.pt[*jt + 3299];
					r__2 = ptjet;
					cthet = hstrng_1.pt[*jt + 3299] / sqrt(r__1 * r__1 + r__2 * r__2);
					r__1 = cthet;
					dpz2 = dp2 * cthet / sqrt(1.f - r__1 * r__1);
					r__1 = dpx2;
					r__2 = dpy2;
					r__3 = dpz2;
					dpe2 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
					etpprm = hstrng_1.pt[*jt + 899] + hstrng_1.pt[*jt + 599] - dpe2 - dpz2;
					etmprm = hstrng_1.pt[*jt + 899] - hstrng_1.pt[*jt + 599] - dpe2 + dpz2;
					if (etpprm <= 0.f || etmprm <= 0.f)
					{
						goto L16;
					}
					etp = etpprm;
					etm = etmprm;
					hstrng_1.pt[*jt + 2699] = pkc21;
					hstrng_1.pt[*jt + 2999] = pkc22;
					++hjjet1_1.ntj[*jt - 1];
					hjjet1_1.kftj[*jt + hjjet1_1.ntj[*jt - 1] * 300 - 301] =
						21;
					hjjet1_1.pjtx[*jt + hjjet1_1.ntj[*jt - 1] * 300 - 301] =
						dpx2;
					hjjet1_1.pjty[*jt + hjjet1_1.ntj[*jt - 1] * 300 - 301] =
						dpy2;
					hjjet1_1.pjtz[*jt + hjjet1_1.ntj[*jt - 1] * 300 - 301] =
						dpz2;
					hjjet1_1.pjte[*jt + hjjet1_1.ntj[*jt - 1] * 300 - 301] =
						dpe2;
					hjjet1_1.pjtm[*jt + hjjet1_1.ntj[*jt - 1] * 300 - 301] =
						0.f;
					hstrng_1.pt[*jt + 599] -= dpz2;
					hstrng_1.pt[*jt + 899] -= dpe2;
				}
			}
		L16:
			dpkc11 = -(hstrng_1.pp[*jp + 2699] - pkc11) / 2.f;
			dpkc12 = -(hstrng_1.pp[*jp + 2999] - pkc12) / 2.f;
			dpkc21 = -(hstrng_1.pt[*jt + 2699] - pkc21) / 2.f;
			dpkc22 = -(hstrng_1.pt[*jt + 2999] - pkc22) / 2.f;
			dpmcm1_1.wp = epp + etp;
			dpmcm1_1.wm = epm + etm;
			dpmcm1_1.sw = dpmcm1_1.wp * dpmcm1_1.wm;
		}
	}

L10:
	r__1 = hstrng_1.pp[*jp - 1];
	r__2 = hstrng_1.pp[*jp + 299];
	ptp02 = r__1 * r__1 + r__2 * r__2;
	r__1 = hstrng_1.pt[*jt - 1];
	r__2 = hstrng_1.pt[*jt + 299];
	ptt02 = r__1 * r__1 + r__2 * r__2;

	r__1 = hstrng_1.pp[*jp + 3899] + hstrng_1.pp[*jp + 4199], r__2 =
																  hstrng_1.pt[*jt + 3899] + hstrng_1.pt[*jt + 4199];
	amq = max(r__1, r__2);
	amx = hparnt_1.hipr1[0] + amq;
	dpmcm1_1.amp0 = amx;
	dpm0 = amx;
	dpmcm1_1.nfdp = 0;
	if (hstrng_1.nfp[*jp + 1199] <= 2 && hstrng_1.nfp[*jp + 599] != 0)
	{
		dpmcm1_1.amp0 = ulmass_(&hstrng_1.nfp[*jp + 599]);
		dpmcm1_1.nfdp = hstrng_1.nfp[*jp + 599] + (hstrng_1.nfp[*jp + 599] << 1) / (i__1 = hstrng_1.nfp[*jp + 599], abs(i__1));
		dpm0 = ulmass_(&dpmcm1_1.nfdp);
		if (dpm0 <= 0.f)
		{
			dpmcm1_1.nfdp -= (dpmcm1_1.nfdp << 1) / abs(dpmcm1_1.nfdp);
			dpm0 = ulmass_(&dpmcm1_1.nfdp);
		}
	}
	dpmcm1_1.amt0 = amx;
	dtm0 = amx;
	dpmcm1_1.nfdt = 0;
	if (hstrng_1.nft[*jt + 1199] <= 2 && hstrng_1.nft[*jt + 599] != 0)
	{
		dpmcm1_1.amt0 = ulmass_(&hstrng_1.nft[*jt + 599]);
		dpmcm1_1.nfdt = hstrng_1.nft[*jt + 599] + (hstrng_1.nft[*jt + 599] << 1) / (i__1 = hstrng_1.nft[*jt + 599], abs(i__1));
		dtm0 = ulmass_(&dpmcm1_1.nfdt);
		if (dtm0 <= 0.f)
		{
			dpmcm1_1.nfdt -= (dpmcm1_1.nfdt << 1) / abs(dpmcm1_1.nfdt);
			dtm0 = ulmass_(&dpmcm1_1.nfdt);
		}
	}

	r__1 = dpmcm1_1.amp0;
	dpmcm1_1.ampn = sqrt(r__1 * r__1 + ptp02);
	r__1 = dpmcm1_1.amt0;
	dpmcm1_1.amtn = sqrt(r__1 * r__1 + ptt02);
	r__1 = dpmcm1_1.ampn + dpmcm1_1.amtn;
	snn = r__1 * r__1 + .001f;

	if (dpmcm1_1.sw < snn + .001f)
	{
		goto L4000;
	}
	r__1 = max(dpmcm1_1.amp0, dpmcm1_1.amt0);
	swptn = (r__1 * r__1 + max(ptp02, ptt02)) * 4.f;
	r__1 = max(dpm0, dtm0);
	swptd = (r__1 * r__1 + max(ptp02, ptt02)) * 4.f;
	r__1 = amx;
	swptx = (r__1 * r__1 + max(ptp02, ptt02)) * 4.f;
	if (dpmcm1_1.sw <= swptn)
	{
		pkcmx = 0.f;
	}
	else if (dpmcm1_1.sw > swptn && dpmcm1_1.sw <= swptd && hjjet1_1.npj[*jp - 1] == 0 && hjjet1_1.ntj[*jt - 1] == 0)
	{
		r__1 = max(dpmcm1_1.amp0, dpmcm1_1.amt0);
		pkcmx = sqrt(dpmcm1_1.sw / 4.f - r__1 * r__1) - sqrt((max(ptp02,
																  ptt02)));
	}
	else if (dpmcm1_1.sw > swptd && dpmcm1_1.sw <= swptx && hjjet1_1.npj[*jp - 1] == 0 && hjjet1_1.ntj[*jt - 1] == 0)
	{
		r__1 = max(dpm0, dtm0);
		pkcmx = sqrt(dpmcm1_1.sw / 4.f - r__1 * r__1) - sqrt((max(ptp02,
																  ptt02)));
	}
	else if (dpmcm1_1.sw > swptx)
	{
		r__1 = amx;
		pkcmx = sqrt(dpmcm1_1.sw / 4.f - r__1 * r__1) - sqrt((max(ptp02,
																  ptt02)));
	}
	if (hstrng_1.nfp[*jp + 2699] == 1 || hstrng_1.nft[*jt + 2699] == 1)
	{
		if (pkc1 > pkcmx)
		{
			pkc1 = pkcmx;
			pkc11 = pkc1 * cos(phi1);
			pkc12 = pkc1 * sin(phi1);
			dpkc11 = -(hstrng_1.pp[*jp + 2699] - pkc11) / 2.f;
			dpkc12 = -(hstrng_1.pp[*jp + 2999] - pkc12) / 2.f;
		}
		if (pkc2 > pkcmx)
		{
			pkc2 = pkcmx;
			pkc21 = pkc2 * cos(phi2);
			pkc22 = pkc2 * sin(phi2);
			dpkc21 = -(hstrng_1.pt[*jt + 2699] - pkc21) / 2.f;
			dpkc22 = -(hstrng_1.pt[*jt + 2999] - pkc22) / 2.f;
		}
		dpmcm1_1.dpkc1 = dpkc11 + dpkc21;
		dpmcm1_1.dpkc2 = dpkc12 + dpkc22;
		hstrng_1.nfp[*jp + 2699] = -hstrng_1.nfp[*jp + 2699];
		hstrng_1.nft[*jt + 2699] = -hstrng_1.nft[*jt + 2699];
		goto L40;
	}
	isng = 0;
	if (hparnt_1.ihpr2[12] != 0 && ranart_(&rndf77_1.nseed) <= hijdat_1.hidat[3])
	{
		isng = 1;
	}
	if (hstrng_1.nfp[*jp + 1199] == 3 || hstrng_1.nft[*jt + 1199] == 3 || (hjjet1_1.npj[*jp - 1] != 0 || hstrng_1.nfp[*jp + 2699] != 0) || (hjjet1_1.ntj[*jt - 1] != 0 || hstrng_1.nft[*jt + 2699] != 0))
	{
		isng = 0;
	}

	if (hparnt_1.ihpr2[4] == 0)
	{
		r__1 = pkcmx;
		r__2 = hparnt_1.hipr1[1];
		pkc = hparnt_1.hipr1[1] * sqrt(-log(1.f - ranart_(&rndf77_1.nseed) * (1.f - exp(-(r__1 * r__1) / (r__2 * r__2)))));
		goto L30;
	}
	xminhi = 0.f;
	r__1 = pkcmx;
	xmaxhi = r__1 * r__1;
	pkc = hirnd2_(&c__3, &xminhi, &xmaxhi);
	pkc = sqrt(pkc);
	if (pkc > hparnt_1.hipr1[19])
	{
		r__1 = hparnt_1.hipr1[19];
		r__2 = hparnt_1.hipr1[1];
		r__3 = hparnt_1.hipr1[19];
		r__4 = hparnt_1.hipr1[1];
		r__5 = pkcmx;
		r__6 = hparnt_1.hipr1[1];
		pkc = hparnt_1.hipr1[1] * sqrt(-log(exp(-(r__1 * r__1) / (r__2 * r__2)) - ranart_(&rndf77_1.nseed) * (exp(-(r__3 * r__3) / (r__4 *
																																	r__4)) -
																											  exp(-(r__5 * r__5) / (r__6 * r__6)))));
	}

	if (isng == 1)
	{
		r__1 = pkcmx;
		pkc = sqrt(-log(1.f - ranart_(&rndf77_1.nseed) * (1.f - exp(-(r__1 *
																	  r__1) /
																	.42250000000000004f)))) *
			  .65f;
	}
L30:
	phi0 = hparnt_1.hipr1[39] * 2.f * ranart_(&rndf77_1.nseed);
	pkc11 = pkc * sin(phi0);
	pkc12 = pkc * cos(phi0);
	pkc21 = -pkc11;
	pkc22 = -pkc12;
	dpmcm1_1.dpkc1 = 0.f;
	dpmcm1_1.dpkc2 = 0.f;
L40:
	dpmcm1_1.pp11 = hstrng_1.pp[*jp - 1] + pkc11 - dpmcm1_1.dpkc1;
	dpmcm1_1.pp12 = hstrng_1.pp[*jp + 299] + pkc12 - dpmcm1_1.dpkc2;
	dpmcm1_1.pt11 = hstrng_1.pt[*jt - 1] + pkc21 - dpmcm1_1.dpkc1;
	dpmcm1_1.pt12 = hstrng_1.pt[*jt + 299] + pkc22 - dpmcm1_1.dpkc2;
	r__1 = dpmcm1_1.pp11;
	r__2 = dpmcm1_1.pp12;
	dpmcm1_1.ptp2 = r__1 * r__1 + r__2 * r__2;
	r__1 = dpmcm1_1.pt11;
	r__2 = dpmcm1_1.pt12;
	dpmcm1_1.ptt2 = r__1 * r__1 + r__2 * r__2;

	r__1 = dpmcm1_1.amp0;
	dpmcm1_1.ampn = sqrt(r__1 * r__1 + dpmcm1_1.ptp2);
	r__1 = dpmcm1_1.amt0;
	dpmcm1_1.amtn = sqrt(r__1 * r__1 + dpmcm1_1.ptt2);
	r__1 = dpmcm1_1.ampn + dpmcm1_1.amtn;
	snn = r__1 * r__1 + .001f;
	dpmcm1_1.wp = epp + etp;
	dpmcm1_1.wm = epm + etm;
	dpmcm1_1.sw = dpmcm1_1.wp * dpmcm1_1.wm;
	if (dpmcm1_1.sw < snn)
	{
		++miss;
		if (miss <= 100)
		{
			pkc = 0.f;
			goto L30;
		}
		if (hparnt_1.ihpr2[9] != 0)
		{
			s_wsle(&io___490);
			do_lio(&c__9, &c__1, "Error occured in Pt kick section of HIJSFT",
				   (int)42);
			e_wsle();
		}
		goto L4000;
	}
	r__1 = dpm0;
	ampd = sqrt(r__1 * r__1 + dpmcm1_1.ptp2);
	r__1 = dtm0;
	amtd = sqrt(r__1 * r__1 + dpmcm1_1.ptt2);
	r__1 = amx;
	ampx = sqrt(r__1 * r__1 + dpmcm1_1.ptp2);
	r__1 = amx;
	amtx = sqrt(r__1 * r__1 + dpmcm1_1.ptt2);
	r__1 = dpmcm1_1.ampn;
	dpn = r__1 * r__1 / dpmcm1_1.sw;
	r__1 = dpmcm1_1.amtn;
	dtn = r__1 * r__1 / dpmcm1_1.sw;
	r__1 = ampd;
	dpd = r__1 * r__1 / dpmcm1_1.sw;
	r__1 = amtd;
	dtd = r__1 * r__1 / dpmcm1_1.sw;
	r__1 = ampx;
	dpx = r__1 * r__1 / dpmcm1_1.sw;
	r__1 = amtx;
	dtx = r__1 * r__1 / dpmcm1_1.sw;

	r__1 = dpmcm1_1.ampn + amtd;
	spntd = r__1 * r__1;
	r__1 = dpmcm1_1.ampn + amtx;
	spntx = r__1 * r__1;
	r__1 = ampd + dpmcm1_1.amtn;
	spdtn = r__1 * r__1;
	r__1 = ampx + dpmcm1_1.amtn;
	spxtn = r__1 * r__1;
	r__1 = ampd + amtx;
	spdtx = r__1 * r__1;
	r__1 = ampx + amtd;
	spxtd = r__1 * r__1;
	r__1 = ampd + amtd;
	sdd = r__1 * r__1;
	r__1 = ampx + amtx;
	sxx = r__1 * r__1;

	if (dpmcm1_1.sw > sxx + .001f)
	{
		if (isng == 0)
		{
			d1 = dpx;
			d2 = dtx;
			nfp3 = 0;
			nft3 = 0;
			goto L400;
		}
		else
		{
			if (hstrng_1.nfp[*jp + 1199] == 3 && hstrng_1.nft[*jt + 1199] ==
													 3 ||
				(hjjet1_1.npj[*jp - 1] != 0 || hstrng_1.nfp[*jp +
															2699] != 0) ||
				(hjjet1_1.ntj[*jt - 1] != 0 ||
				 hstrng_1.nft[*jt + 2699] != 0))
			{
				d1 = dpx;
				d2 = dtx;
				nfp3 = 0;
				nft3 = 0;
				goto L400;
			}
			if (ranart_(&rndf77_1.nseed) > .5f || (hstrng_1.nft[*jt + 1199] >
													   2 ||
												   hjjet1_1.ntj[*jt - 1] != 0 || hstrng_1.nft[*jt + 2699] != 0))
			{
				d1 = dpn;
				d2 = dtx;
				nfp3 = hstrng_1.nfp[*jp + 599];
				nft3 = 0;
				goto L220;
			}
			else
			{
				d1 = dpx;
				d2 = dtn;
				nfp3 = 0;
				nft3 = hstrng_1.nft[*jt + 599];
				goto L240;
			}
		}
	}
	else if (dpmcm1_1.sw > max(spdtx, spxtd) + .001f && dpmcm1_1.sw <= sxx + .001f)
	{
		if ((hjjet1_1.npj[*jp - 1] == 0 && hjjet1_1.ntj[*jt - 1] == 0 &&
				 ranart_(&rndf77_1.nseed) > .5f ||
			 hjjet1_1.npj[*jp - 1] == 0 && hjjet1_1.ntj[*jt - 1] != 0) &&
			hstrng_1.nfp[*jp + 1199] <=
				2)
		{
			d1 = dpd;
			d2 = dtx;
			nfp3 = dpmcm1_1.nfdp;
			nft3 = 0;
			goto L220;
		}
		else if (hjjet1_1.ntj[*jt - 1] == 0 && hstrng_1.nft[*jt + 1199] <=
												   2)
		{
			d1 = dpx;
			d2 = dtd;
			nfp3 = 0;
			nft3 = dpmcm1_1.nfdt;
			goto L240;
		}
		goto L4000;
	}
	else if (dpmcm1_1.sw > min(spdtx, spxtd) + .001f && dpmcm1_1.sw <= max(
																		   spdtx, spxtd) +
																		   .001f)
	{
		if (spdtx <= spxtd && hjjet1_1.npj[*jp - 1] == 0 && hstrng_1.nfp[*jp + 1199] <= 2)
		{
			d1 = dpd;
			d2 = dtx;
			nfp3 = dpmcm1_1.nfdp;
			nft3 = 0;
			goto L220;
		}
		else if (spdtx > spxtd && hjjet1_1.ntj[*jt - 1] == 0 &&
				 hstrng_1.nft[*jt + 1199] <= 2)
		{
			d1 = dpx;
			d2 = dtd;
			nfp3 = 0;
			nft3 = dpmcm1_1.nfdt;
			goto L240;
		}
		if ((hjjet1_1.npj[*jp - 1] == 0 && hjjet1_1.ntj[*jt - 1] == 0 &&
				 ranart_(&rndf77_1.nseed) > .5f ||
			 hjjet1_1.npj[*jp - 1] == 0 && hjjet1_1.ntj[*jt - 1] != 0) &&
			hstrng_1.nfp[*jp + 1199] <=
				2)
		{
			d1 = dpn;
			d2 = dtx;
			nfp3 = hstrng_1.nfp[*jp + 599];
			nft3 = 0;
			goto L220;
		}
		else if (hjjet1_1.ntj[*jt - 1] == 0 && hstrng_1.nft[*jt + 1199] <=
												   2)
		{
			d1 = dpx;
			d2 = dtn;
			nfp3 = 0;
			nft3 = hstrng_1.nft[*jt + 599];
			goto L240;
		}
		goto L4000;
	}
	else if (dpmcm1_1.sw > max(spntx, spxtn) + .001f && dpmcm1_1.sw <= min(
																		   spdtx, spxtd) +
																		   .001f)
	{
		if ((hjjet1_1.npj[*jp - 1] == 0 && hjjet1_1.ntj[*jt - 1] == 0 &&
				 ranart_(&rndf77_1.nseed) > .5f ||
			 hjjet1_1.npj[*jp - 1] == 0 && hjjet1_1.ntj[*jt - 1] != 0) &&
			hstrng_1.nfp[*jp + 1199] <=
				2)
		{
			d1 = dpn;
			d2 = dtx;
			nfp3 = hstrng_1.nfp[*jp + 599];
			nft3 = 0;
			goto L220;
		}
		else if (hjjet1_1.ntj[*jt - 1] == 0 && hstrng_1.nft[*jt + 1199] <=
												   2)
		{
			d1 = dpx;
			d2 = dtn;
			nfp3 = 0;
			nft3 = hstrng_1.nft[*jt + 599];
			goto L240;
		}
		goto L4000;
	}
	else if (dpmcm1_1.sw > min(spntx, spxtn) + .001f && dpmcm1_1.sw <= max(
																		   spntx, spxtn) +
																		   .001f)
	{
		if (spntx <= spxtn && hjjet1_1.npj[*jp - 1] == 0 && hstrng_1.nfp[*jp + 1199] <= 2)
		{
			d1 = dpn;
			d2 = dtx;
			nfp3 = hstrng_1.nfp[*jp + 599];
			nft3 = 0;
			goto L220;
		}
		else if (spntx > spxtn && hjjet1_1.ntj[*jt - 1] == 0 &&
				 hstrng_1.nft[*jt + 1199] <= 2)
		{
			d1 = dpx;
			d2 = dtn;
			nfp3 = 0;
			nft3 = hstrng_1.nft[*jt + 599];
			goto L240;
		}
		goto L4000;
	}
	else if (dpmcm1_1.sw <= min(spntx, spxtn) + .001f && (hjjet1_1.npj[*jp - 1] != 0 || hjjet1_1.ntj[*jt - 1] != 0))
	{
		goto L4000;
	}
	else if (dpmcm1_1.sw <= min(spntx, spxtn) + .001f && hstrng_1.nfp[*jp + 1199] > 2 && hstrng_1.nft[*jt + 1199] > 2)
	{
		goto L4000;
	}
	else if (dpmcm1_1.sw > sdd + .001f && dpmcm1_1.sw <= min(spntx, spxtn) + .001f)
	{
		d1 = dpd;
		d2 = dtd;
		nfp3 = dpmcm1_1.nfdp;
		nft3 = dpmcm1_1.nfdt;
		goto L100;
	}
	else if (dpmcm1_1.sw > max(spntd, spdtn) + .001f && dpmcm1_1.sw <= sdd + .001f)
	{
		if (ranart_(&rndf77_1.nseed) > .5f)
		{
			d1 = dpd;
			d2 = dtn;
			nfp3 = dpmcm1_1.nfdp;
			nft3 = hstrng_1.nft[*jt + 599];
			goto L100;
		}
		else
		{
			d1 = dpn;
			d2 = dtd;
			nfp3 = hstrng_1.nfp[*jp + 599];
			nft3 = dpmcm1_1.nfdt;
			goto L100;
		}
	}
	else if (dpmcm1_1.sw > min(spntd, spdtn) + .001f && dpmcm1_1.sw <= max(
																		   spntd, spdtn) +
																		   .001f)
	{
		if (spntd > spdtn)
		{
			d1 = dpd;
			d2 = dtn;
			nfp3 = dpmcm1_1.nfdp;
			nft3 = hstrng_1.nft[*jt + 599];
			goto L100;
		}
		else
		{
			d1 = dpn;
			d2 = dtd;
			nfp3 = hstrng_1.nfp[*jp + 599];
			nft3 = dpmcm1_1.nfdt;
			goto L100;
		}
	}
	else if (dpmcm1_1.sw <= min(spntd, spdtn) + .001f)
	{
		d1 = dpn;
		d2 = dtn;
		nfp3 = hstrng_1.nfp[*jp + 599];
		nft3 = hstrng_1.nft[*jt + 599];
		goto L100;
	}
	s_wsle(&io___513);
	do_lio(&c__9, &c__1, " Error in HIJSFT: There is no path to here", (int)42);
	e_wsle();
	return 0;

L100:
	i__1 = 2, i__2 = hstrng_1.nfp[*jp + 1199];
	nfp5 = max(i__1, i__2);
	i__1 = 2, i__2 = hstrng_1.nft[*jt + 1199];
	nft5 = max(i__1, i__2);
	bb1 = d1 + 1.f - d2;
	bb2 = d2 + 1.f - d1;
	r__1 = bb1;
	r__2 = bb2;
	if (r__1 * r__1 < d1 * 4.f || r__2 * r__2 < d2 * 4.f)
	{
		++miss;
		if (miss > 100 || pkc == 0.f)
		{
			goto L3000;
		}
		pkc *= .5f;
		goto L30;
	}
	if (ranart_(&rndf77_1.nseed) < .5f)
	{
		r__1 = bb1;
		x1 = (bb1 - sqrt(r__1 * r__1 - d1 * 4.f)) / 2.f;
		r__1 = bb2;
		x2 = (bb2 - sqrt(r__1 * r__1 - d2 * 4.f)) / 2.f;
	}
	else
	{
		r__1 = bb1;
		x1 = (bb1 + sqrt(r__1 * r__1 - d1 * 4.f)) / 2.f;
		r__1 = bb2;
		x2 = (bb2 + sqrt(r__1 * r__1 - d2 * 4.f)) / 2.f;
	}
	hparnt_1.ihnt2[12] = 2;
	goto L600;

L220:
	i__1 = 2, i__2 = hstrng_1.nfp[*jp + 1199];
	nfp5 = max(i__1, i__2);
	nft5 = 3;
	if (nfp3 == 0)
	{
		nfp5 = 3;
	}
	bb2 = d2 + 1.f - d1;
	r__1 = bb2;
	if (r__1 * r__1 < d2 * 4.f)
	{
		++miss;
		if (miss > 100 || pkc == 0.f)
		{
			goto L3000;
		}
		pkc *= .5f;
		goto L30;
	}
	r__1 = bb2;
	xmin = (bb2 - sqrt(r__1 * r__1 - d2 * 4.f)) / 2.f;
	r__1 = bb2;
	xmax = (bb2 + sqrt(r__1 * r__1 - d2 * 4.f)) / 2.f;
	miss4 = 0;
L222:
	x2 = hirnd2_(&c__6, &xmin, &xmax);
	x1 = d1 / (1.f - x2);
	if (x2 * (1.f - x1) < d2 + 1e-4f / dpmcm1_1.sw)
	{
		++miss4;
		if (miss4 <= 1000)
		{
			goto L222;
		}
		goto L5000;
	}
	hparnt_1.ihnt2[12] = 2;
	goto L600;
L240:
	nfp5 = 3;
	i__1 = 2, i__2 = hstrng_1.nft[*jt + 1199];
	nft5 = max(i__1, i__2);
	if (nft3 == 0)
	{
		nft5 = 3;
	}
	bb1 = d1 + 1.f - d2;
	r__1 = bb1;
	if (r__1 * r__1 < d1 * 4.f)
	{
		++miss;
		if (miss > 100 || pkc == 0.f)
		{
			goto L3000;
		}
		pkc *= .5f;
		goto L30;
	}
	r__1 = bb1;
	xmin = (bb1 - sqrt(r__1 * r__1 - d1 * 4.f)) / 2.f;
	r__1 = bb1;
	xmax = (bb1 + sqrt(r__1 * r__1 - d1 * 4.f)) / 2.f;
	miss4 = 0;
L242:
	x1 = hirnd2_(&c__6, &xmin, &xmax);
	x2 = d2 / (1.f - x1);
	if (x1 * (1.f - x2) < d1 + 1e-4f / dpmcm1_1.sw)
	{
		++miss4;
		if (miss4 <= 1000)
		{
			goto L242;
		}
		goto L5000;
	}
	hparnt_1.ihnt2[12] = 2;
	goto L600;
L400:
	nfp5 = 3;
	nft5 = 3;
	bb1 = d1 + 1.f - d2;
	bb2 = d2 + 1.f - d1;
	r__1 = bb1;
	r__2 = bb2;
	if (r__1 * r__1 < d1 * 4.f || r__2 * r__2 < d2 * 4.f)
	{
		++miss;
		if (miss > 100 || pkc == 0.f)
		{
			goto L3000;
		}
		pkc *= .5f;
		goto L30;
	}
	r__1 = bb1;
	xmin1 = (bb1 - sqrt(r__1 * r__1 - d1 * 4.f)) / 2.f;
	r__1 = bb1;
	xmax1 = (bb1 + sqrt(r__1 * r__1 - d1 * 4.f)) / 2.f;
	r__1 = bb2;
	xmin2 = (bb2 - sqrt(r__1 * r__1 - d2 * 4.f)) / 2.f;
	r__1 = bb2;
	xmax2 = (bb2 + sqrt(r__1 * r__1 - d2 * 4.f)) / 2.f;
	miss4 = 0;
L410:
	x1 = hirnd2_(&c__4, &xmin1, &xmax1);
	x2 = hirnd2_(&c__4, &xmin2, &xmax2);
	if (hstrng_1.nfp[*jp + 1199] == 3 || hstrng_1.nft[*jt + 1199] == 3)
	{
		x1 = hirnd2_(&c__6, &xmin1, &xmax1);
		x2 = hirnd2_(&c__6, &xmin2, &xmax2);
	}
	if ((i__1 = hstrng_1.nfp[*jp - 1] * hstrng_1.nfp[*jp + 299], abs(i__1)) >
			1000000 ||
		(i__2 = hstrng_1.nfp[*jp - 1] * hstrng_1.nfp[*jp + 299], abs(i__2)) < 100)
	{
		x1 = hirnd2_(&c__5, &xmin1, &xmax1);
	}
	if ((i__1 = hstrng_1.nft[*jt - 1] * hstrng_1.nft[*jt + 299], abs(i__1)) >
			1000000 ||
		(i__2 = hstrng_1.nft[*jt - 1] * hstrng_1.nft[*jt + 299], abs(i__2)) < 100)
	{
		x2 = hirnd2_(&c__5, &xmin2, &xmax2);
	}
	if ((i__1 = hstrng_1.nfp[*jp - 1] * hstrng_1.nfp[*jp + 299], abs(i__1)) >
		1000000)
	{
		x1 = 1.f - x1;
	}
	xxp = x1 * (1.f - x2);
	xxt = x2 * (1.f - x1);
	if (xxp < d1 + 1e-4f / dpmcm1_1.sw || xxt < d2 + 1e-4f / dpmcm1_1.sw)
	{
		++miss4;
		if (miss4 <= 1000)
		{
			goto L410;
		}
		goto L5000;
	}
	hparnt_1.ihnt2[12] = 3;
L600:
	r__1 = dpmcm1_1.ampn;
	r__2 = dpmcm1_1.amtn;
	if (x1 * (1.f - x2) < (r__1 * r__1 - 1e-4f) / dpmcm1_1.sw || x2 * (1.f -
																	   x1) <
																	 (r__2 * r__2 - 1e-4f) / dpmcm1_1.sw)
	{
		++miss;
		if (miss > 100 || pkc == 0.f)
		{
			goto L2000;
		}
		pkc = 0.f;
		goto L30;
	}

	epp = (1.f - x2) * dpmcm1_1.wp;
	epm = x1 * dpmcm1_1.wm;
	etp = x2 * dpmcm1_1.wp;
	etm = (1.f - x1) * dpmcm1_1.wm;
	hstrng_1.pp[*jp + 599] = (epp - epm) / 2.f;
	hstrng_1.pp[*jp + 899] = (epp + epm) / 2.f;
	if (epp * epm - dpmcm1_1.ptp2 < 0.f)
	{
		goto L6000;
	}
	hstrng_1.pp[*jp + 1199] = sqrt(epp * epm - dpmcm1_1.ptp2);
	hstrng_1.nfp[*jp + 599] = nfp3;
	hstrng_1.nfp[*jp + 1199] = nfp5;
	hstrng_1.pt[*jt + 599] = (etp - etm) / 2.f;
	hstrng_1.pt[*jt + 899] = (etp + etm) / 2.f;
	if (etp * etm - dpmcm1_1.ptt2 < 0.f)
	{
		goto L6000;
	}
	hstrng_1.pt[*jt + 1199] = sqrt(etp * etm - dpmcm1_1.ptt2);
	hstrng_1.nft[*jt + 599] = nft3;
	hstrng_1.nft[*jt + 1199] = nft5;
	hstrng_1.pp[*jp - 1] = dpmcm1_1.pp11 - pkc11;
	hstrng_1.pp[*jp + 299] = dpmcm1_1.pp12 - pkc12;
	kcdip = 1;
	kcdit = 1;
	if ((i__1 = hstrng_1.nfp[*jp - 1] * hstrng_1.nfp[*jp + 299], abs(i__1)) >
			1000000 ||
		(i__2 = hstrng_1.nfp[*jp - 1] * hstrng_1.nfp[*jp + 299], abs(i__2)) < 100)
	{
		kcdip = 0;
	}
	if ((i__1 = hstrng_1.nft[*jt - 1] * hstrng_1.nft[*jt + 299], abs(i__1)) >
			1000000 ||
		(i__2 = hstrng_1.nft[*jt - 1] * hstrng_1.nft[*jt + 299], abs(i__2)) < 100)
	{
		kcdit = 0;
	}
	r__1 = pkc11;
	r__2 = pkc12;
	r__3 = hparnt_1.hipr1[21];
	if (kcdip == 0 && ranart_(&rndf77_1.nseed) < .5f || kcdip != 0 && ranart_(
																		  &rndf77_1.nseed) < .5f / ((r__1 * r__1 + r__2 * r__2) / (r__3 *
																																   r__3) +
																									1.f))
	{
		hstrng_1.pp[*jp + 1499] = (hstrng_1.pp[*jp - 1] - hstrng_1.pp[*jp + 1499] - hstrng_1.pp[*jp + 2099] - dpmcm1_1.dpkc1) / 2.f +
								  hstrng_1.pp[*jp + 1499];
		hstrng_1.pp[*jp + 1799] = (hstrng_1.pp[*jp + 299] - hstrng_1.pp[*jp + 1799] - hstrng_1.pp[*jp + 2399] - dpmcm1_1.dpkc2) / 2.f +
								  hstrng_1.pp[*jp + 1799];
		hstrng_1.pp[*jp + 2099] = (hstrng_1.pp[*jp - 1] - hstrng_1.pp[*jp + 1499] - hstrng_1.pp[*jp + 2099] - dpmcm1_1.dpkc1) / 2.f +
								  hstrng_1.pp[*jp + 2099] + pkc11;
		hstrng_1.pp[*jp + 2399] = (hstrng_1.pp[*jp + 299] - hstrng_1.pp[*jp + 1799] - hstrng_1.pp[*jp + 2399] - dpmcm1_1.dpkc2) / 2.f +
								  hstrng_1.pp[*jp + 2399] + pkc12;
	}
	else
	{
		hstrng_1.pp[*jp + 2099] = (hstrng_1.pp[*jp - 1] - hstrng_1.pp[*jp + 1499] - hstrng_1.pp[*jp + 2099] - dpmcm1_1.dpkc1) / 2.f +
								  hstrng_1.pp[*jp + 2099];
		hstrng_1.pp[*jp + 2399] = (hstrng_1.pp[*jp + 299] - hstrng_1.pp[*jp + 1799] - hstrng_1.pp[*jp + 2399] - dpmcm1_1.dpkc2) / 2.f +
								  hstrng_1.pp[*jp + 2399];
		hstrng_1.pp[*jp + 1499] = (hstrng_1.pp[*jp - 1] - hstrng_1.pp[*jp + 1499] - hstrng_1.pp[*jp + 2099] - dpmcm1_1.dpkc1) / 2.f +
								  hstrng_1.pp[*jp + 1499] + pkc11;
		hstrng_1.pp[*jp + 1799] = (hstrng_1.pp[*jp + 299] - hstrng_1.pp[*jp + 1799] - hstrng_1.pp[*jp + 2399] - dpmcm1_1.dpkc2) / 2.f +
								  hstrng_1.pp[*jp + 1799] + pkc12;
	}
	hstrng_1.pp[*jp - 1] = hstrng_1.pp[*jp + 1499] + hstrng_1.pp[*jp + 2099];
	hstrng_1.pp[*jp + 299] = hstrng_1.pp[*jp + 1799] + hstrng_1.pp[*jp + 2399];
	hstrng_1.pt[*jt - 1] = dpmcm1_1.pt11 - pkc21;
	hstrng_1.pt[*jt + 299] = dpmcm1_1.pt12 - pkc22;
	r__1 = pkc21;
	r__2 = pkc22;
	r__3 = hparnt_1.hipr1[21];
	if (kcdit == 0 && ranart_(&rndf77_1.nseed) < .5f || kcdit != 0 && ranart_(
																		  &rndf77_1.nseed) < .5f / ((r__1 * r__1 + r__2 * r__2) / (r__3 *
																																   r__3) +
																									1.f))
	{
		hstrng_1.pt[*jt + 1499] = (hstrng_1.pt[*jt - 1] - hstrng_1.pt[*jt + 1499] - hstrng_1.pt[*jt + 2099] - dpmcm1_1.dpkc1) / 2.f +
								  hstrng_1.pt[*jt + 1499];
		hstrng_1.pt[*jt + 1799] = (hstrng_1.pt[*jt + 299] - hstrng_1.pt[*jt + 1799] - hstrng_1.pt[*jt + 2399] - dpmcm1_1.dpkc2) / 2.f +
								  hstrng_1.pt[*jt + 1799];
		hstrng_1.pt[*jt + 2099] = (hstrng_1.pt[*jt - 1] - hstrng_1.pt[*jt + 1499] - hstrng_1.pt[*jt + 2099] - dpmcm1_1.dpkc1) / 2.f +
								  hstrng_1.pt[*jt + 2099] + pkc21;
		hstrng_1.pt[*jt + 2399] = (hstrng_1.pt[*jt + 299] - hstrng_1.pt[*jt + 1799] - hstrng_1.pt[*jt + 2399] - dpmcm1_1.dpkc2) / 2.f +
								  hstrng_1.pt[*jt + 2399] + pkc22;
	}
	else
	{
		hstrng_1.pt[*jt + 2099] = (hstrng_1.pt[*jt - 1] - hstrng_1.pt[*jt + 1499] - hstrng_1.pt[*jt + 2099] - dpmcm1_1.dpkc1) / 2.f +
								  hstrng_1.pt[*jt + 2099];
		hstrng_1.pt[*jt + 2399] = (hstrng_1.pt[*jt + 299] - hstrng_1.pt[*jt + 1799] - hstrng_1.pt[*jt + 2399] - dpmcm1_1.dpkc2) / 2.f +
								  hstrng_1.pt[*jt + 2399];
		hstrng_1.pt[*jt + 1499] = (hstrng_1.pt[*jt - 1] - hstrng_1.pt[*jt + 1499] - hstrng_1.pt[*jt + 2099] - dpmcm1_1.dpkc1) / 2.f +
								  hstrng_1.pt[*jt + 1499] + pkc21;
		hstrng_1.pt[*jt + 1799] = (hstrng_1.pt[*jt + 299] - hstrng_1.pt[*jt + 1799] - hstrng_1.pt[*jt + 2399] - dpmcm1_1.dpkc2) / 2.f +
								  hstrng_1.pt[*jt + 1799] + pkc22;
	}
	hstrng_1.pt[*jt - 1] = hstrng_1.pt[*jt + 1499] + hstrng_1.pt[*jt + 2099];
	hstrng_1.pt[*jt + 299] = hstrng_1.pt[*jt + 1799] + hstrng_1.pt[*jt + 2399];
	if (hjjet1_1.npj[*jp - 1] != 0)
	{
		hstrng_1.nfp[*jp + 1199] = 3;
	}
	if (hjjet1_1.ntj[*jt - 1] != 0)
	{
		hstrng_1.nft[*jt + 1199] = 3;
	}
	if (epp / (epm + 1e-4f) < etp / (etm + 1e-4f) && (i__1 = hstrng_1.nfp[*jp - 1] * hstrng_1.nfp[*jp + 299], abs(i__1)) < 1000000)
	{
		for (jsb = 1; jsb <= 15; ++jsb)
		{
			psb = hstrng_1.pp[*jp + jsb * 300 - 301];
			hstrng_1.pp[*jp + jsb * 300 - 301] = hstrng_1.pt[*jt + jsb * 300 - 301];
			hstrng_1.pt[*jt + jsb * 300 - 301] = psb;
			nsb = hstrng_1.nfp[*jp + jsb * 300 - 301];
			hstrng_1.nfp[*jp + jsb * 300 - 301] = hstrng_1.nft[*jt + jsb * 300 - 301];
			hstrng_1.nft[*jt + jsb * 300 - 301] = nsb;
		}
	}

	return 0;
L1000:
	*ierror = 1;
	if (hparnt_1.ihpr2[9] == 0)
	{
		return 0;
	}
	s_wsle(&io___534);
	do_lio(&c__9, &c__1, "     Fatal HIJSFT start error,abandon this event", (int)48);
	e_wsle();
	s_wsle(&io___535);
	do_lio(&c__9, &c__1, "     PROJ E+,E-,W+", (int)18);
	do_lio(&c__4, &c__1, (char *)&epp, (int)sizeof(float));
	do_lio(&c__4, &c__1, (char *)&epm, (int)sizeof(float));
	do_lio(&c__4, &c__1, (char *)&dpmcm1_1.wp, (int)sizeof(float));
	e_wsle();
	s_wsle(&io___536);
	do_lio(&c__9, &c__1, "     TARG E+,E-,W-", (int)18);
	do_lio(&c__4, &c__1, (char *)&etp, (int)sizeof(float));
	do_lio(&c__4, &c__1, (char *)&etm, (int)sizeof(float));
	do_lio(&c__4, &c__1, (char *)&dpmcm1_1.wm, (int)sizeof(float));
	e_wsle();
	s_wsle(&io___537);
	do_lio(&c__9, &c__1, "     W+*W-, (APN+ATN)^2", (int)23);
	do_lio(&c__4, &c__1, (char *)&dpmcm1_1.sw, (int)sizeof(float));
	do_lio(&c__4, &c__1, (char *)&snn, (int)sizeof(float));
	e_wsle();
	return 0;
L2000:
	*ierror = 0;
	if (hparnt_1.ihpr2[9] == 0)
	{
		return 0;
	}
	s_wsle(&io___538);
	do_lio(&c__9, &c__1, "     (2)energy partition fail,", (int)30);
	e_wsle();
	s_wsle(&io___539);
	do_lio(&c__9, &c__1, "     HIJSFT not performed, but continue", (int)39);
	e_wsle();
	s_wsle(&io___540);
	do_lio(&c__9, &c__1, "     MP1,MPN", (int)12);
	r__1 = x1 * (1.f - x2) * dpmcm1_1.sw;
	do_lio(&c__4, &c__1, (char *)&r__1, (int)sizeof(float));
	r__3 = dpmcm1_1.ampn;
	r__2 = r__3 * r__3;
	do_lio(&c__4, &c__1, (char *)&r__2, (int)sizeof(float));
	e_wsle();
	s_wsle(&io___541);
	do_lio(&c__9, &c__1, "     MT2,MTN", (int)12);
	r__1 = x2 * (1.f - x1) * dpmcm1_1.sw;
	do_lio(&c__4, &c__1, (char *)&r__1, (int)sizeof(float));
	r__3 = dpmcm1_1.amtn;
	r__2 = r__3 * r__3;
	do_lio(&c__4, &c__1, (char *)&r__2, (int)sizeof(float));
	e_wsle();
	return 0;
L3000:
	*ierror = 0;
	if (hparnt_1.ihpr2[9] == 0)
	{
		return 0;
	}
	s_wsle(&io___542);
	do_lio(&c__9, &c__1, "     (3)something is wrong with the pt kick, ", (int)45);
	e_wsle();
	s_wsle(&io___543);
	do_lio(&c__9, &c__1, "     HIJSFT not performed, but continue", (int)39);
	e_wsle();
	s_wsle(&io___544);
	do_lio(&c__9, &c__1, "     D1=", (int)8);
	do_lio(&c__4, &c__1, (char *)&d1, (int)sizeof(float));
	do_lio(&c__9, &c__1, " D2=", (int)4);
	do_lio(&c__4, &c__1, (char *)&d2, (int)sizeof(float));
	do_lio(&c__9, &c__1, " SW=", (int)4);
	do_lio(&c__4, &c__1, (char *)&dpmcm1_1.sw, (int)sizeof(float));
	e_wsle();
	s_wsle(&io___545);
	do_lio(&c__9, &c__1, "     HISTORY NFP5=", (int)18);
	do_lio(&c__3, &c__1, (char *)&hstrng_1.nfp[*jp + 1199], (int)sizeof(int));
	do_lio(&c__9, &c__1, " NFT5=", (int)6);
	do_lio(&c__3, &c__1, (char *)&hstrng_1.nft[*jt + 1199], (int)sizeof(int));
	e_wsle();
	s_wsle(&io___546);
	do_lio(&c__9, &c__1, "     THIS COLLISON NFP5=", (int)24);
	do_lio(&c__3, &c__1, (char *)&nfp5, (int)sizeof(int));
	do_lio(&c__9, &c__1, " NFT5=", (int)6);
	do_lio(&c__3, &c__1, (char *)&nft5, (int)sizeof(int));
	e_wsle();
	s_wsle(&io___547);
	do_lio(&c__9, &c__1, "     # OF JET IN PROJ", (int)21);
	do_lio(&c__3, &c__1, (char *)&hjjet1_1.npj[*jp - 1], (int)sizeof(int));
	do_lio(&c__9, &c__1, " IN TARG", (int)8);
	do_lio(&c__3, &c__1, (char *)&hjjet1_1.ntj[*jt - 1], (int)sizeof(int));
	e_wsle();
	return 0;
L4000:
	*ierror = 0;
	if (hparnt_1.ihpr2[9] == 0)
	{
		return 0;
	}
	s_wsle(&io___548);
	do_lio(&c__9, &c__1, "     (4)unable to choose process, but not harmful",
		   (int)49);
	e_wsle();
	s_wsle(&io___549);
	do_lio(&c__9, &c__1, "     HIJSFT not performed, but continue", (int)39);
	e_wsle();
	s_wsle(&io___550);
	do_lio(&c__9, &c__1, "     PTP=", (int)9);
	r__1 = sqrt(dpmcm1_1.ptp2);
	do_lio(&c__4, &c__1, (char *)&r__1, (int)sizeof(float));
	do_lio(&c__9, &c__1, " PTT=", (int)5);
	r__2 = sqrt(dpmcm1_1.ptt2);
	do_lio(&c__4, &c__1, (char *)&r__2, (int)sizeof(float));
	do_lio(&c__9, &c__1, " SW=", (int)4);
	do_lio(&c__4, &c__1, (char *)&dpmcm1_1.sw, (int)sizeof(float));
	e_wsle();
	s_wsle(&io___551);
	do_lio(&c__9, &c__1, "     AMCUT=", (int)11);
	do_lio(&c__4, &c__1, (char *)&amx, (int)sizeof(float));
	do_lio(&c__9, &c__1, " JP=", (int)4);
	do_lio(&c__3, &c__1, (char *)&(*jp), (int)sizeof(int));
	do_lio(&c__9, &c__1, " JT=", (int)4);
	do_lio(&c__3, &c__1, (char *)&(*jt), (int)sizeof(int));
	e_wsle();
	s_wsle(&io___552);
	do_lio(&c__9, &c__1, "     HISTORY NFP5=", (int)18);
	do_lio(&c__3, &c__1, (char *)&hstrng_1.nfp[*jp + 1199], (int)sizeof(int));
	do_lio(&c__9, &c__1, " NFT5=", (int)6);
	do_lio(&c__3, &c__1, (char *)&hstrng_1.nft[*jt + 1199], (int)sizeof(int));
	e_wsle();
	return 0;
L5000:
	*ierror = 0;
	if (hparnt_1.ihpr2[9] == 0)
	{
		return 0;
	}
	s_wsle(&io___553);
	do_lio(&c__9, &c__1, "     energy partition failed(5),for limited try", (int)47);
	e_wsle();
	s_wsle(&io___554);
	do_lio(&c__9, &c__1, "     HIJSFT not performed, but continue", (int)39);
	e_wsle();
	s_wsle(&io___555);
	do_lio(&c__9, &c__1, "     NFP5=", (int)10);
	do_lio(&c__3, &c__1, (char *)&nfp5, (int)sizeof(int));
	do_lio(&c__9, &c__1, " NFT5=", (int)6);
	do_lio(&c__3, &c__1, (char *)&nft5, (int)sizeof(int));
	e_wsle();
	s_wsle(&io___556);
	do_lio(&c__9, &c__1, "     D1", (int)7);
	do_lio(&c__4, &c__1, (char *)&d1, (int)sizeof(float));
	do_lio(&c__9, &c__1, " X1(1-X2)", (int)9);
	r__1 = x1 * (1.f - x2);
	do_lio(&c__4, &c__1, (char *)&r__1, (int)sizeof(float));
	e_wsle();
	s_wsle(&io___557);
	do_lio(&c__9, &c__1, "     D2", (int)7);
	do_lio(&c__4, &c__1, (char *)&d2, (int)sizeof(float));
	do_lio(&c__9, &c__1, " X2(1-X1)", (int)9);
	r__1 = x2 * (1.f - x1);
	do_lio(&c__4, &c__1, (char *)&r__1, (int)sizeof(float));
	e_wsle();
	return 0;
L6000:
	pkc = 0.f;
	++miss;
	if (miss < 100)
	{
		goto L30;
	}
	*ierror = 1;
	if (hparnt_1.ihpr2[9] == 0)
	{
		return 0;
	}
	s_wsle(&io___558);
	do_lio(&c__9, &c__1, " ERROR OCCURED, HIJSFT NOT PERFORMED", (int)36);
	e_wsle();
	s_wsle(&io___559);
	do_lio(&c__9, &c__1, " Abort this event", (int)17);
	e_wsle();
	s_wsle(&io___560);
	do_lio(&c__9, &c__1, "MTP,PTP2", (int)8);
	r__1 = epp * epm;
	do_lio(&c__4, &c__1, (char *)&r__1, (int)sizeof(float));
	do_lio(&c__4, &c__1, (char *)&dpmcm1_1.ptp2, (int)sizeof(float));
	do_lio(&c__9, &c__1, "  MTT,PTT2", (int)10);
	r__2 = etp * etm;
	do_lio(&c__4, &c__1, (char *)&r__2, (int)sizeof(float));
	do_lio(&c__4, &c__1, (char *)&dpmcm1_1.ptt2, (int)sizeof(float));
	e_wsle();
	return 0;
}

int hijwds_(int *ia, int *idh, float *xhigh)
{
	static int iaa[20] = {2, 4, 12, 16, 27, 32, 40, 56, 63, 93, 184, 197, 208, 0, 0, 0,
						  0, 0, 0, 0};
	static float rr[20] = {.01f, .964f, 2.355f, 2.608f, 2.84f, 3.458f, 3.766f,
						   3.971f, 4.214f, 4.87f, 6.51f, 6.38f, 6.624f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
						   0.f};
	static float dd[20] = {.5882f, .322f, .522f, .513f, .569f, .61f, .586f, .5935f,
						   .586f, .573f, .535f, .535f, .549f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
	static float ww[20] = {0.f, .517f, -.149f, -.051f, 0.f, -.208f, -.161f, 0.f, 0.f,
						   0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f};

	double d__1, d__2;

	double pow_dd(double *, double *), sqrt(double);

	static float a;
	static int i__;
	static float xlow, fgaus;
	extern int hifun_(int *, float *, float *, int *);
	double gauss1_(int *, float *, float *, float *);
	extern int rwdsax_();

	a = (float)(*ia);

	wood_1.d__ = .54f;
	d__1 = (double)a;
	d__2 = (double)a;
	wood_1.r__ = pow_dd(&d__1, &c_b876) * 1.19f - pow_dd(&d__2, &c_b877) *
													  1.61f;
	wood_1.w = 0.f;
	for (i__ = 1; i__ <= 13; ++i__)
	{
		if (*ia == iaa[i__ - 1])
		{
			wood_1.r__ = rr[i__ - 1];
			wood_1.d__ = dd[i__ - 1];
			wood_1.w = ww[i__ - 1];
		}
	}
	wood_1.fnorm = 1.f;
	xlow = 0.f;
	*xhigh = wood_1.r__ + wood_1.d__ * 12.f;
	if (wood_1.w < -.01f)
	{
		if (*xhigh > wood_1.r__ / sqrt((dabs(wood_1.w))))
		{
			*xhigh = wood_1.r__ / sqrt((dabs(wood_1.w)));
		}
	}
	fgaus = gauss1_((int *)rwdsax_, &xlow, xhigh, &c_b879);
	wood_1.fnorm = 1.f / fgaus;

	if (*idh == 1)
	{
		hparnt_1.hint1[71] = wood_1.r__;
		hparnt_1.hint1[72] = wood_1.d__;
		hparnt_1.hint1[73] = wood_1.w;
		hparnt_1.hint1[74] = wood_1.fnorm / 4.f / hparnt_1.hipr1[39];
	}
	else if (*idh == 2)
	{
		hparnt_1.hint1[75] = wood_1.r__;
		hparnt_1.hint1[76] = wood_1.d__;
		hparnt_1.hint1[77] = wood_1.w;
		hparnt_1.hint1[78] = wood_1.fnorm / 4.f / hparnt_1.hipr1[39];
	}

	hifun_(idh, &xlow, xhigh, (int *)rwdsax_);
	return 0;
}

double wdsax_(float *x)
{
	float ret_val, r__1;

	double exp(double), sqrt(double);

	r__1 = *x / wood_1.r__;
	ret_val = wood_1.fnorm * (wood_1.w * (r__1 * r__1) + 1.f) / (exp((*x - wood_1.r__) / wood_1.d__) + 1);
	if (wood_1.w < 0.f)
	{
		if (*x >= wood_1.r__ / sqrt((dabs(wood_1.w))))
		{
			ret_val = 0.f;
		}
	}
	return ret_val;
}

int rwdsax_(float *x)
{
	float ret_val;

	extern double wdsax_(float *);

	ret_val = *x * *x * wdsax_(x);
	return ret_val;
}

int hifun_(int *i__, float *xmin, float *xmax, int *fhb)
{
	static int j;
	static float xdd, fnorm;
	double gauss1_(int *, float *, float *, float *);

	fnorm = gauss1_((int *)fhb, xmin, xmax, &c_b879);
	for (j = 1; j <= 201; ++j)
	{
		hijhb_1.xx[*i__ + j * 10 - 11] = *xmin + (*xmax - *xmin) * (j - 1) /
													 200.f;
		xdd = hijhb_1.xx[*i__ + j * 10 - 11];
		hijhb_1.rr[*i__ + j * 10 - 11] = gauss1_((int *)fhb, xmin, &xdd, &c_b879) / fnorm;
	}
	return 0;
}

double hirnd_(int *i__)
{
	float ret_val;

	static int j, jl, jm, ju;
	static float rx;
	extern double ranart_(int *);

	rx = ranart_(&rndf77_1.nseed);
	jl = 0;
	ju = 202;
L10:
	if (ju - jl > 1)
	{
		jm = (ju + jl) / 2;
		if (hijhb_1.rr[*i__ + 1999] > hijhb_1.rr[*i__ - 1] == rx > hijhb_1.rr[*i__ + jm * 10 - 11])
		{
			jl = jm;
		}
		else
		{
			ju = jm;
		}
		goto L10;
	}
	j = jl;
	if (j < 1)
	{
		j = 1;
	}
	if (j >= 201)
	{
		j = 200;
	}
	ret_val = (hijhb_1.xx[*i__ + j * 10 - 11] + hijhb_1.xx[*i__ + (j + 1) * 10 - 11]) / 2.f;
	return ret_val;
}

double hirnd2_(int *i__, float *xmin, float *xmax)
{
	float ret_val;

	static int j, jl, jm, ju;
	static float rx;
	static int jmin, jmax;
	extern double ranart_(int *);

	if (*xmin < hijhb_1.xx[*i__ - 1])
	{
		*xmin = hijhb_1.xx[*i__ - 1];
	}
	if (*xmax > hijhb_1.xx[*i__ + 1999])
	{
		*xmax = hijhb_1.xx[*i__ + 1999];
	}
	jmin = (int)((*xmin - hijhb_1.xx[*i__ - 1]) * 200 / (hijhb_1.xx[*i__ + 1999] - hijhb_1.xx[*i__ - 1])) + 1;
	jmax = (int)((*xmax - hijhb_1.xx[*i__ - 1]) * 200 / (hijhb_1.xx[*i__ + 1999] - hijhb_1.xx[*i__ - 1])) + 1;
	rx = hijhb_1.rr[*i__ + jmin * 10 - 11] + (hijhb_1.rr[*i__ + jmax * 10 -
														 11] -
											  hijhb_1.rr[*i__ + jmin * 10 - 11]) *
												 ranart_(&rndf77_1.nseed);
	jl = 0;
	ju = 202;
L10:
	if (ju - jl > 1)
	{
		jm = (ju + jl) / 2;
		if (hijhb_1.rr[*i__ + 1999] > hijhb_1.rr[*i__ - 1] == rx > hijhb_1.rr[*i__ + jm * 10 - 11])
		{
			jl = jm;
		}
		else
		{
			ju = jm;
		}
		goto L10;
	}
	j = jl;
	if (j < 1)
	{
		j = 1;
	}
	if (j >= 201)
	{
		j = 200;
	}
	ret_val = (hijhb_1.xx[*i__ + j * 10 - 11] + hijhb_1.xx[*i__ + (j + 1) * 10 - 11]) / 2.f;
	return ret_val;
}

int hijcrs_(void)
{
	float r__1, r__2;
	double d__1;

	double pow_dd(double *, double *), log(double);

	static int i__;
	extern int fhin_(), ftot_();
	static float aphx1, aphx2;
	extern int fnjet_();
	double gauss1_(int *, float *, float *, float *);
	extern int crsjet_(void);
	extern int ftotrg_(), ftotjt_();

	if (hparnt_1.hint1[0] >= 10.f)
	{
		crsjet_();
	}
	d__1 = (double)hparnt_1.ihnt2[0];
	aphx1 = hparnt_1.hipr1[5] * (pow_dd(&d__1, &c_b23) - 1.f);
	d__1 = (double)hparnt_1.ihnt2[2];
	aphx2 = hparnt_1.hipr1[5] * (pow_dd(&d__1, &c_b23) - 1.f);
	hparnt_1.hint1[10] = hparnt_1.hint1[13] - aphx1 * hparnt_1.hint1[14] -
						 aphx2 * hparnt_1.hint1[15] + aphx1 * aphx2 * hparnt_1.hint1[16];
	hparnt_1.hint1[9] = gauss1_((int *)ftotjt_, &c_b26, &c_b894, &c_b895);
	hparnt_1.hint1[11] = gauss1_((int *)fhin_, &c_b26, &c_b894, &c_b895);
	hparnt_1.hint1[12] = gauss1_((int *)ftot_, &c_b26, &c_b894, &c_b895);
	hparnt_1.hint1[59] = hparnt_1.hint1[60] - aphx1 * hparnt_1.hint1[61] -
						 aphx2 * hparnt_1.hint1[62] + aphx1 * aphx2 * hparnt_1.hint1[63];
	hparnt_1.hint1[58] = gauss1_((int *)ftotrg_, &c_b26, &c_b894, &c_b895);
	if (hparnt_1.hint1[58] == 0.f)
	{
		hparnt_1.hint1[58] = hparnt_1.hint1[59];
	}
	if (hparnt_1.hint1[0] >= 10.f)
	{
		for (i__ = 0; i__ <= 20; ++i__)
		{
			njet_1.n = i__;
			hparnt_1.hint1[i__ + 79] = gauss1_((int *)fnjet_, &c_b26, &c_b894,
											   &c_b895) /
									   hparnt_1.hint1[11];
		}
	}
	hparnt_1.hint1[9] *= hparnt_1.hipr1[30];
	hparnt_1.hint1[11] *= hparnt_1.hipr1[30];
	hparnt_1.hint1[12] *= hparnt_1.hipr1[30];
	hparnt_1.hint1[58] *= hparnt_1.hipr1[30];
	if (hparnt_1.ihpr2[12] != 0)
	{
		r__1 = hparnt_1.hint1[0];
		r__2 = hparnt_1.hint1[0];
		hparnt_1.hipr1[32] = (36.f / (r__1 * r__1) + 1.f) * 1.36f * log(r__2 * r__2 * .1f + .6f);
		hparnt_1.hipr1[32] /= hparnt_1.hint1[11];
	}
	return 0;
}

int ftot_(float *x)
{
	float ret_val;

	double exp(double);

	static float omg;
	extern double omg0_(float *);

	omg = omg0_(x) * (hparnt_1.hipr1[29] + hparnt_1.hint1[10]) /
		  hparnt_1.hipr1[30] / 2.f;
	ret_val = (1.f - exp(-omg)) * 2.f;
	return ret_val;
}

int fhin_(float *x)
{
	float ret_val;

	double exp(double);

	static float omg;
	extern double omg0_(float *);

	omg = omg0_(x) * (hparnt_1.hipr1[29] + hparnt_1.hint1[10]) /
		  hparnt_1.hipr1[30] / 2.f;
	ret_val = 1.f - exp(omg * -2.f);
	return ret_val;
}

int ftotjt_(float *x)
{
	float ret_val;

	double exp(double);

	static float omg;
	extern double omg0_(float *);

	omg = omg0_(x) * hparnt_1.hint1[10] / hparnt_1.hipr1[30] / 2.f;
	ret_val = 1.f - exp(omg * -2.f);
	return ret_val;
}

int ftotrg_(float *x)
{
	float ret_val;

	double exp(double);

	static float omg;
	extern double omg0_(float *);

	omg = omg0_(x) * hparnt_1.hint1[59] / hparnt_1.hipr1[30] / 2.f;
	ret_val = 1.f - exp(omg * -2.f);
	return ret_val;
}

int fnjet_(float *x)
{
	int i__1;
	float ret_val;

	double log(double), exp(double);

	static float c0;
	extern double omg0_(float *);
	static float omg1;
	extern double sgmin_(int *);

	omg1 = omg0_(x) * hparnt_1.hint1[10] / hparnt_1.hipr1[30];
	i__1 = njet_1.n + 1;
	c0 = exp(njet_1.n * log(omg1) - sgmin_(&i__1));
	if (njet_1.n == 0)
	{
		c0 = 1.f - exp(omg0_(x) * -2.f * hparnt_1.hipr1[29] / hparnt_1.hipr1[30] / 2.f);
	}
	ret_val = c0 * exp(-omg1);
	return ret_val;
}

double sgmin_(int *n)
{
	int i__1;
	float ret_val;

	double log(double);

	static int i__;
	static float z__, ga;

	ga = 0.f;
	if (*n <= 2)
	{
		goto L20;
	}
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		z__ = (float)i__;
		ga += log(z__);
	}
L20:
	ret_val = ga;
	return ret_val;
}

double omg0_(float *x)
{
	float ret_val, r__1, r__2;

	double sqrt(double);

	extern int bk_();
	extern double gauss2_(int *, float *, float *, float *);

	besel_1.x4 = hparnt_1.hipr1[31] * sqrt(*x);
	r__1 = hparnt_1.hipr1[31];
	r__2 = besel_1.x4 + 20.f;
	ret_val = r__1 * r__1 * gauss2_((int *)bk_, &besel_1.x4, &r__2, &c_b895) /
			  96.f;
	return ret_val;
}

double romg_(float *x)
{
	static int i0 = 0;

	float ret_val;

	static int i__;
	static float fr[1001];
	static int ix;
	static float xr;
	extern double omg0_(float *);

	if (i0 != 0)
	{
		goto L100;
	}
	for (i__ = 1; i__ <= 1001; ++i__)
	{
		xr = (i__ - 1) * .01f;
		fr[i__ - 1] = omg0_(&xr);
	}
L100:
	i0 = 1;
	if (*x >= 10.f)
	{
		ret_val = 0.f;
		return ret_val;
	}
	ix = (int)(*x * 100);
	ret_val = (fr[ix] * ((ix + 1) * .01f - *x) + fr[ix + 1] * (*x - ix * .01f)) / .01f;
	return ret_val;
}

int bk_(float *x)
{
	float ret_val, r__1, r__2;
	double d__1;

	double exp(double), pow_dd(double *, double *);

	r__1 = *x;
	r__2 = besel_1.x4;
	d__1 = (double)(r__1 * r__1 - r__2 * r__2);
	ret_val = exp(-(*x)) * pow_dd(&d__1, &c_b923) / 15.f;
	return ret_val;
}

int crsjet_(void)
{
	static double sd;
	extern int fjet_();
	static double avgi, chi2a;
	extern int vegas_(int *, double *, double *,
					  double *);
	extern int fjetrg_();

	bveg1_1.ndim = 3;
	njet_1.ipcrs = 0;
	vegas_((int *)fjet_, &avgi, &sd, &chi2a);
	hparnt_1.hint1[13] = (float)avgi / 2.5682f;
	if (hparnt_1.ihpr2[5] == 1 && hparnt_1.ihnt2[0] > 1)
	{
		njet_1.ipcrs = 1;
		vegas_((int *)fjet_, &avgi, &sd, &chi2a);
		hparnt_1.hint1[14] = (float)avgi / 2.5682f;
	}
	if (hparnt_1.ihpr2[5] == 1 && hparnt_1.ihnt2[2] > 1)
	{
		njet_1.ipcrs = 2;
		vegas_((int *)fjet_, &avgi, &sd, &chi2a);
		hparnt_1.hint1[15] = (float)avgi / 2.5682f;
	}
	if (hparnt_1.ihpr2[5] == 1 && hparnt_1.ihnt2[0] > 1 && hparnt_1.ihnt2[2] > 1)
	{
		njet_1.ipcrs = 3;
		vegas_((int *)fjet_, &avgi, &sd, &chi2a);
		hparnt_1.hint1[16] = (float)avgi / 2.5682f;
	}
	if (hparnt_1.ihpr2[2] != 0)
	{
		njet_1.ipcrs = 0;
		vegas_((int *)fjetrg_, &avgi, &sd, &chi2a);
		hparnt_1.hint1[60] = (float)avgi / 2.5682f;
		if (hparnt_1.ihpr2[5] == 1 && hparnt_1.ihnt2[0] > 1)
		{
			njet_1.ipcrs = 1;
			vegas_((int *)fjetrg_, &avgi, &sd, &chi2a);
			hparnt_1.hint1[61] = (float)avgi / 2.5682f;
		}
		if (hparnt_1.ihpr2[5] == 1 && hparnt_1.ihnt2[2] > 1)
		{
			njet_1.ipcrs = 2;
			vegas_((int *)fjetrg_, &avgi, &sd, &chi2a);
			hparnt_1.hint1[62] = (float)avgi / 2.5682f;
		}
		if (hparnt_1.ihpr2[5] == 1 && hparnt_1.ihnt2[0] > 1 && hparnt_1.ihnt2[2] > 1)
		{
			njet_1.ipcrs = 3;
			vegas_((int *)fjetrg_, &avgi, &sd, &chi2a);
			hparnt_1.hint1[63] = (float)avgi / 2.5682f;
		}
	}
	return 0;
}

int fjet_(double *x, double *wgt)
{
	float r__1, r__2;
	double ret_val, d__1;

	double sqrt(double), log(double), exp(double);

	extern double g_(double *, double *, double *);
	static double y1, y2, xt, pt2, ymn2, ymx1, ymx2;

	--x;

	r__1 = hparnt_1.hint1[0];
	r__2 = hparnt_1.hipr1[7];
	d__1 = (double)hparnt_1.hipr1[7];
	pt2 = (double)(r__1 * r__1 / 4.f - r__2 * r__2) * x[1] + d__1 * d__1;
	xt = sqrt(pt2) * 2. / (double)hparnt_1.hint1[0];
	d__1 = xt;
	ymx1 = log(1. / xt + sqrt(1. / (d__1 * d__1) - 1.));
	y1 = ymx1 * 2. * x[2] - ymx1;
	ymx2 = log(2. / xt - exp(y1));
	ymn2 = log(2. / xt - exp(-y1));
	y2 = (ymx2 + ymn2) * x[3] - ymn2;
	r__1 = hparnt_1.hint1[0];
	r__2 = hparnt_1.hipr1[7];
	ret_val = ymx1 * 2. * (ymx2 + ymn2) * (double)(r__1 * r__1 / 4.f - r__2 * r__2) * g_(&y1, &y2, &pt2) / 2.;
	return ret_val;
}

int fjetrg_(double *x, double *wgt)
{
	float r__1, r__2;
	double ret_val, d__1;

	double sqrt(double), log(double), exp(double);

	extern double g_(double *, double *, double *);
	static double y1, y2, xt, am2, pt2, amt2, ymn2, ymx1, ymx2;
	extern double ghvq_(double *, double *, double *);
	static double gtrig;
	static float ptmin, ptmax;
	extern double gphotn_(double *, double *, double *);

	--x;

	ptmin = dabs(hparnt_1.hipr1[9]) - .25f;
	ptmin = max(ptmin, hparnt_1.hipr1[7]);
	am2 = 0.;
	if (hparnt_1.ihpr2[2] == 3)
	{
		r__1 = hparnt_1.hipr1[6];
		am2 = (double)(r__1 * r__1);
		ptmin = max(0.f, hparnt_1.hipr1[9]);
	}
	ptmax = dabs(hparnt_1.hipr1[9]) + .25f;
	if (hparnt_1.hipr1[9] <= 0.f)
	{
		ptmax = hparnt_1.hint1[0] / 2.f - (float)am2;
	}
	if (ptmax <= ptmin)
	{
		ptmax = ptmin + .25f;
	}
	r__1 = ptmax;
	r__2 = ptmin;
	d__1 = (double)ptmin;
	pt2 = (double)(r__1 * r__1 - r__2 * r__2) * x[1] + d__1 * d__1;
	amt2 = pt2 + am2;
	xt = sqrt(amt2) * 2. / (double)hparnt_1.hint1[0];
	d__1 = xt;
	ymx1 = log(1. / xt + sqrt(1. / (d__1 * d__1) - 1.));
	y1 = ymx1 * 2. * x[2] - ymx1;
	ymx2 = log(2. / xt - exp(y1));
	ymn2 = log(2. / xt - exp(-y1));
	y2 = (ymx2 + ymn2) * x[3] - ymn2;
	if (hparnt_1.ihpr2[2] == 3)
	{
		gtrig = ghvq_(&y1, &y2, &amt2) * 2.;
	}
	else if (hparnt_1.ihpr2[2] == 2)
	{
		gtrig = gphotn_(&y1, &y2, &pt2) * 2.;
	}
	else
	{
		gtrig = g_(&y1, &y2, &pt2);
	}
	r__1 = ptmax;
	r__2 = ptmin;
	ret_val = ymx1 * 2. * (ymx2 + ymn2) * (double)(r__1 * r__1 - r__2 * r__2) * gtrig / 2.;
	return ret_val;
}

double ghvq_(double *y1, double *y2, double *amt2)
{
	double ret_val, d__1, d__2, d__3;

	double sqrt(double), exp(double), log(double), cosh(double);

	static double f[14], x1, x2, af, ss, xt, ggg,
		aph, gqq, dlam;
	int parton_(double *, double *,
				double *, double *);

	xt = sqrt(*amt2) * 2. / (double)hparnt_1.hint1[0];
	x1 = xt * .5 * (exp(*y1) + exp(*y2));
	x2 = xt * .5 * (exp(-(*y1)) + exp(-(*y2)));
	d__1 = (double)hparnt_1.hint1[0];
	ss = x1 * x2 * (d__1 * d__1);
	af = 4.;
	if (hparnt_1.ihpr2[17] != 0)
	{
		af = 5.;
	}
	dlam = (double)hparnt_1.hipr1[14];
	d__1 = dlam;
	aph = 37.699111200000004 / (33. - af * 2.) / log(*amt2 / (d__1 * d__1));

	parton_(f, &x1, &x2, amt2);

	d__1 = (double)hparnt_1.hipr1[6];
	gqq = (cosh(*y1 - *y2) + d__1 * d__1 / *amt2) * 4. / (cosh(*y1 - *y2) + 1.) / 9. * (f[0] * f[3] + f[2] * f[1] + f[4] * f[7] + f[6] * f[5] + f[8] * f[11] + f[10] * f[9]);
	d__1 = (double)hparnt_1.hipr1[6];
	d__2 = (double)hparnt_1.hipr1[6], d__2 *= d__2;
	d__3 = *amt2;
	ggg = (cosh(*y1 - *y2) * 8. - 1.) * (cosh(*y1 - *y2) + d__1 * d__1 * 2. / *amt2 - d__2 * d__2 * 2. / (d__3 * d__3)) / (cosh(*y1 - *y2) + 1.) / 24. * f[12] * f[13];

	d__1 = aph;
	d__2 = ss;
	ret_val = (gqq + ggg) * (double)hparnt_1.hipr1[22] * 3.14159 * (d__1 * d__1) / (d__2 * d__2);
	return ret_val;
}

double gphotn_(double *y1, double *y2, double *pt2)
{
	double ret_val, d__1, d__2;

	double sqrt(double), exp(double), log(double);

	static double f[14], t, u, z__, g2, x1, x2, af,
		g11, g12, ss, xt, aph, dlam, aphem;
	int parton_(double *, double *,
				double *, double *);

	xt = sqrt(*pt2) * 2. / (double)hparnt_1.hint1[0];
	x1 = xt * .5 * (exp(*y1) + exp(*y2));
	x2 = xt * .5 * (exp(-(*y1)) + exp(-(*y2)));
	d__1 = xt;
	z__ = sqrt(1. - d__1 * d__1 / x1 / x2);
	d__1 = (double)hparnt_1.hint1[0];
	ss = x1 * x2 * (d__1 * d__1);
	t = -(1. - z__) / 2.;
	u = -(z__ + 1.) / 2.;
	af = 3.;
	dlam = (double)hparnt_1.hipr1[14];
	d__1 = dlam;
	aph = 37.699111200000004 / (33. - af * 2.) / log(*pt2 / (d__1 * d__1));
	aphem = .0072992700729927005;

	parton_(f, &x1, &x2, pt2);

	d__1 = u;
	g11 = -(d__1 * d__1 + 1.) / u / 3. * f[12] * (f[1] * 4. + f[3] * 4. + f[5] + f[7] + f[9] + f[11]) / 9.;
	d__1 = t;
	g12 = -(d__1 * d__1 + 1.) / t / 3. * f[13] * (f[0] * 4. + f[2] * 4. + f[4] + f[6] + f[8] + f[10]) / 9.;
	d__1 = u;
	d__2 = t;
	g2 = (d__1 * d__1 + d__2 * d__2) * 8. / u / t / 9. * (f[0] * 4. * f[3] + f[2] * 4. * f[1] + f[4] * f[7] + f[6] * f[5] + f[8] * f[11] + f[10] * f[9]) / 9.;

	d__1 = ss;
	ret_val = (g11 + g12 + g2) * (double)hparnt_1.hipr1[22] * 3.14159 *
			  aph * aphem / (d__1 * d__1);
	return ret_val;
}

double g_(double *y1, double *y2, double *pt2)
{
	double ret_val, d__1, d__2;

	double sqrt(double), exp(double), log(double);

	static double f[14], t, u, z__, g2, g4, g5, g7,
		x1, x2, af, g11, g12, g13, g31, g32, g61, g62, ss, xt, aph, dlam;
	extern double subcr1_(double *, double *), subcr2_(double *, double *), subcr3_(double *, double *), subcr4_(double *, double *), subcr5_(double *, double *),
		subcr6_(double *, double *), subcr7_(double *, double *);
	int parton_(double *, double *,
				double *, double *);

	xt = sqrt(*pt2) * 2. / (double)hparnt_1.hint1[0];
	x1 = xt * .5 * (exp(*y1) + exp(*y2));
	x2 = xt * .5 * (exp(-(*y1)) + exp(-(*y2)));
	d__1 = xt;
	z__ = sqrt(1. - d__1 * d__1 / x1 / x2);
	d__1 = (double)hparnt_1.hint1[0];
	ss = x1 * x2 * (d__1 * d__1);
	t = -(1. - z__) / 2.;
	u = -(z__ + 1.) / 2.;
	af = 3.;
	dlam = (double)hparnt_1.hipr1[14];
	d__1 = dlam;
	aph = 37.699111200000004 / (33. - af * 2.) / log(*pt2 / (d__1 * d__1));

	parton_(f, &x1, &x2, pt2);

	g11 = ((f[0] + f[2]) * (f[5] + f[7] + f[9] + f[11]) + (f[4] + f[6]) * (f[9] + f[11])) * subcr1_(&t, &u);

	g12 = ((f[1] + f[3]) * (f[4] + f[6] + f[8] + f[10]) + (f[5] + f[7]) * (f[8] + f[10])) * subcr1_(&u, &t);

	g13 = (f[0] * f[1] + f[2] * f[3] + f[4] * f[5] + f[6] * f[7] + f[8] * f[9] + f[10] * f[11]) * (subcr1_(&u, &t) + subcr1_(&t, &u) - 8. / t / u / 27.);

	g2 = (af - 1) * (f[0] * f[3] + f[1] * f[2] + f[4] * f[7] + f[5] * f[6] + f[8] * f[11] + f[9] * f[10]) * subcr2_(&t, &u);

	g31 = (f[0] * f[3] + f[4] * f[7] + f[8] * f[11]) * subcr3_(&t, &u);
	g32 = (f[1] * f[2] + f[5] * f[6] + f[9] * f[10]) * subcr3_(&u, &t);

	g4 = (f[0] * f[3] + f[1] * f[2] + f[4] * f[7] + f[5] * f[6] + f[8] * f[11] + f[9] * f[10]) * subcr4_(&t, &u);

	g5 = af * f[12] * f[13] * subcr5_(&t, &u);

	g61 = f[12] * (f[1] + f[3] + f[5] + f[7] + f[9] + f[11]) * subcr6_(&t, &u);
	g62 = f[13] * (f[0] + f[2] + f[4] + f[6] + f[8] + f[10]) * subcr6_(&u, &t);

	g7 = f[12] * f[13] * subcr7_(&t, &u);

	d__1 = aph;
	d__2 = ss;
	ret_val = (g11 + g12 + g13 + g2 + g31 + g32 + g4 + g5 + g61 + g62 + g7) *
			  (double)hparnt_1.hipr1[16] * 3.14159 * (d__1 * d__1) / (d__2 * d__2);
	return ret_val;
}

double subcr1_(double *t, double *u)
{
	double ret_val, d__1, d__2;

	d__1 = *u;
	d__2 = *t;
	ret_val = (d__1 * d__1 + 1.) * .44444444444444442 / (d__2 * d__2);
	return ret_val;
}

double subcr2_(double *t, double *u)
{
	double ret_val, d__1, d__2;

	d__1 = *t;
	d__2 = *u;
	ret_val = (d__1 * d__1 + d__2 * d__2) * .44444444444444442;
	return ret_val;
}

double subcr3_(double *t, double *u)
{
	double ret_val, d__1, d__2, d__3, d__4, d__5;

	d__1 = *t;
	d__2 = *u;
	d__3 = *u;
	d__4 = *t;
	d__5 = *u;
	ret_val = (d__1 * d__1 + d__2 * d__2 + (d__3 * d__3 + 1.) / (d__4 * d__4) - d__5 * d__5 * 2. / 3. / *t) * .44444444444444442;
	return ret_val;
}

double subcr4_(double *t, double *u)
{
	double ret_val, d__1, d__2;

	d__1 = *t;
	d__2 = *u;
	ret_val = (d__1 * d__1 + d__2 * d__2) * 2.6666666666666665 * (.44444444444444442 / *t / *u - 1.);
	return ret_val;
}

double subcr5_(double *t, double *u)
{
	double ret_val, d__1, d__2;

	d__1 = *t;
	d__2 = *u;
	ret_val = (d__1 * d__1 + d__2 * d__2) * .375 * (.44444444444444442 / *t / *u - 1.);
	return ret_val;
}

double subcr6_(double *t, double *u)
{
	double ret_val, d__1, d__2;

	d__1 = *u;
	d__2 = *t;
	ret_val = (d__1 * d__1 + 1.) * (1. / (d__2 * d__2) - 4. / *u / 9.);
	return ret_val;
}

double subcr7_(double *t, double *u)
{
	double ret_val, d__1, d__2;

	d__1 = *t;
	d__2 = *u;
	ret_val = (3. - *t * *u - *u / (d__1 * d__1) - *t / (d__2 * d__2)) * 4.5;
	return ret_val;
}

int parton_(double *f, double *x1, double *x2,
			double *qq)
{
	double d__1, d__2, d__3, d__4;

	double log(double), exp(double), pow_dd(double *, double *), sqrt(double);

	static int i__;
	static double s, q0, b12, ag, bg, b34, as, bs, at1, at2, at3, at4, fs1, fs2, cag, cnd, cas, gmd, aax, gmg, gms, rrx, fud1, fud2, dlam, btag, aphg, btas;
	extern double gmre_(double *);
	static double aphs, gmud, cnud;

	f -= 3;

	dlam = (double)hparnt_1.hipr1[14];
	q0 = (double)hparnt_1.hipr1[15];
	d__1 = dlam;
	d__2 = q0;
	d__3 = dlam;
	s = log(log(*qq / (d__1 * d__1)) / log(d__2 * d__2 / (d__3 * d__3)));
	if (hparnt_1.ihpr2[6] == 2)
	{
		goto L200;
	}
	d__1 = s;
	at1 = s * .004 + .419 - d__1 * d__1 * .007;
	d__1 = s;
	at2 = s * .724 + 3.46 - d__1 * d__1 * .066;
	d__1 = s;
	gmud = 4.4 - s * 4.86 + d__1 * d__1 * 1.33;
	d__1 = s;
	at3 = .763 - s * .237 + d__1 * d__1 * .026;
	d__1 = s;
	at4 = s * .627 + 4. - d__1 * d__1 * .019;
	d__1 = s;
	gmd = s * -.421 + d__1 * d__1 * .033;
	d__1 = s;
	cas = 1.265 - s * 1.132 + d__1 * d__1 * .293;
	d__1 = s;
	as = s * -.372 - d__1 * d__1 * .029;
	d__1 = s;
	bs = s * 1.59 + 8.05 - d__1 * d__1 * .153;
	d__1 = s;
	aphs = s * 6.31 - d__1 * d__1 * .273;
	d__1 = s;
	btas = s * -10.5 - d__1 * d__1 * 3.17;
	d__1 = s;
	gms = s * 14.7 + d__1 * d__1 * 9.8;
	d__1 = s;
	cag = 1.56 - s * 1.71 + d__1 * d__1 * .638;
	d__1 = s;
	ag = s * -.949 + d__1 * d__1 * .325;
	d__1 = s;
	bg = s * 1.44 + 6. - d__1 * d__1 * 1.05;
	d__1 = s;
	aphg = 9. - s * 7.19 + d__1 * d__1 * .255;
	d__1 = s;
	btag = s * -16.5 + d__1 * d__1 * 10.9;
	d__1 = s;
	gmg = s * 15.3 - d__1 * d__1 * 10.1;
	goto L300;
L200:
	at1 = s * .014 + .374;
	d__1 = s;
	at2 = s * .753 + 3.33 - d__1 * d__1 * .076;
	d__1 = s;
	gmud = 6.03 - s * 6.22 + d__1 * d__1 * 1.56;
	d__1 = s;
	at3 = .761 - s * .232 + d__1 * d__1 * .023;
	d__1 = s;
	at4 = s * .627 + 3.83 - d__1 * d__1 * .019;
	d__1 = s;
	gmd = s * -.418 + d__1 * d__1 * .036;
	d__1 = s;
	cas = 1.67 - s * 1.92 + d__1 * d__1 * .582;
	d__1 = s;
	as = s * -.273 - d__1 * d__1 * .164;
	d__1 = s;
	bs = s * .53 + 9.15 - d__1 * d__1 * .763;
	d__1 = s;
	aphs = s * 15.7 - d__1 * d__1 * 2.83;
	d__1 = s;
	btas = s * -101. + d__1 * d__1 * 44.7;
	d__1 = s;
	gms = s * 223. - d__1 * d__1 * 117.;
	d__1 = s;
	cag = .879 - s * .971 + d__1 * d__1 * .434;
	d__1 = s;
	ag = s * -1.16 + d__1 * d__1 * .476;
	d__1 = s;
	bg = s * 1.23 + 4. - d__1 * d__1 * .254;
	d__1 = s;
	aphg = 9. - s * 5.64 - d__1 * d__1 * .817;
	d__1 = s;
	btag = s * -7.54 + d__1 * d__1 * 5.5;
	d__1 = s;
	gmg = s * -.596 + d__1 * d__1 * 1.26;
L300:
	d__1 = at2 + 1.;
	d__2 = at1 + at2 + 1.;
	b12 = exp(gmre_(&at1) + gmre_(&d__1) - gmre_(&d__2));
	d__1 = at4 + 1.;
	d__2 = at3 + at4 + 1.;
	b34 = exp(gmre_(&at3) + gmre_(&d__1) - gmre_(&d__2));
	cnud = 3. / b12 / (gmud * at1 / (at1 + at2 + 1.) + 1.);
	cnd = 1. / b34 / (gmd * at3 / (at3 + at4 + 1.) + 1.);
	d__1 = 1. - *x1;
	fud1 = cnud * pow_dd(x1, &at1) * pow_dd(&d__1, &at2) * (gmud * *x1 + 1.);
	d__1 = 1. - *x1;
	d__2 = *x1;
	d__3 = *x1;
	fs1 = cas * pow_dd(x1, &as) * pow_dd(&d__1, &bs) * (aphs * *x1 + 1. + btas * (d__2 * d__2) + gms * (d__3 * (d__3 * d__3)));
	d__1 = 1. - *x1;
	f[7] = cnd * pow_dd(x1, &at3) * pow_dd(&d__1, &at4) * (gmd * *x1 + 1.) +
		   fs1 / 6.;
	f[3] = fud1 - f[7] + fs1 / 3.;
	f[5] = fs1 / 6.;
	f[9] = fs1 / 6.;
	f[11] = fs1 / 6.;
	f[13] = fs1 / 6.;
	d__1 = 1. - *x1;
	d__2 = *x1;
	d__3 = *x1;
	f[15] = cag * pow_dd(x1, &ag) * pow_dd(&d__1, &bg) * (aphg * *x1 + 1. + btag * (d__2 * d__2) + gmg * (d__3 * (d__3 * d__3)));

	d__1 = 1. - *x2;
	fud2 = cnud * pow_dd(x2, &at1) * pow_dd(&d__1, &at2) * (gmud * *x2 + 1.);
	d__1 = 1. - *x2;
	d__2 = *x2;
	d__3 = *x2;
	fs2 = cas * pow_dd(x2, &as) * pow_dd(&d__1, &bs) * (aphs * *x2 + 1. + btas * (d__2 * d__2) + gms * (d__3 * (d__3 * d__3)));
	d__1 = 1. - *x2;
	f[8] = cnd * pow_dd(x2, &at3) * pow_dd(&d__1, &at4) * (gmd * *x2 + 1.) +
		   fs2 / 6.;
	f[4] = fud2 - f[8] + fs2 / 3.;
	f[6] = fs2 / 6.;
	f[10] = fs2 / 6.;
	f[12] = fs2 / 6.;
	f[14] = fs2 / 6.;
	d__1 = 1. - *x2;
	d__2 = *x2;
	d__3 = *x2;
	f[16] = cag * pow_dd(x2, &ag) * pow_dd(&d__1, &bg) * (aphg * *x2 + 1. + btag * (d__2 * d__2) + gmg * (d__3 * (d__3 * d__3)));
	if (hparnt_1.ihpr2[5] == 1 && hparnt_1.ihnt2[0] > 1)
	{
		d__1 = (double)log((float)hparnt_1.ihnt2[0]);
		aax = (double)pow_dd(&d__1, &c_b940) * 1.193;
		d__1 = *x1;
		d__2 = *x1;
		d__3 = (double)((float)hparnt_1.ihnt2[0]);
		d__4 = *x1;
		rrx = aax * (d__1 * (d__1 * d__1) - d__2 * d__2 * 1.2 + *x1 * .21) +
			  1. + (double)((pow_dd(&d__3, &c_b941) - 1.f) * 1.079f) / (double)log((float)hparnt_1.ihnt2[0] + 1.f) * sqrt(*x1) * exp(-(d__4 * d__4) / .01);
		if (njet_1.ipcrs == 1 || njet_1.ipcrs == 3)
		{
			d__1 = *x1;
			rrx = exp(-(d__1 * d__1) / .01);
		}
		for (i__ = 1; i__ <= 7; ++i__)
		{
			f[(i__ << 1) + 1] = rrx * f[(i__ << 1) + 1];
		}
	}
	if (hparnt_1.ihpr2[5] == 1 && hparnt_1.ihnt2[2] > 1)
	{
		d__1 = (double)log((float)hparnt_1.ihnt2[2]);
		aax = (double)pow_dd(&d__1, &c_b940) * 1.193;
		d__1 = *x2;
		d__2 = *x2;
		d__3 = (double)((float)hparnt_1.ihnt2[2]);
		d__4 = *x2;
		rrx = aax * (d__1 * (d__1 * d__1) - d__2 * d__2 * 1.2 + *x2 * .21) +
			  1. + (double)((pow_dd(&d__3, &c_b944) - 1.f) * 1.079f) / (double)log((float)hparnt_1.ihnt2[2] + 1.f) * sqrt(*x2) * exp(-(d__4 * d__4) / .01);
		if (njet_1.ipcrs == 2 || njet_1.ipcrs == 3)
		{
			d__1 = *x2;
			rrx = exp(-(d__1 * d__1) / .01);
		}
		for (i__ = 1; i__ <= 7; ++i__)
		{
			f[(i__ << 1) + 2] = rrx * f[(i__ << 1) + 2];
		}
	}

	return 0;
}

double gmre_(double *x)
{
	double ret_val, d__1, d__2, d__3;

	double log(double);

	static double z__;

	z__ = *x;
	if (*x > 3.)
	{
		goto L10;
	}
	z__ = *x + 3.;
L10:
	d__1 = z__;
	d__2 = z__;
	d__3 = z__, d__3 *= d__3;
	ret_val = log(6.2831853000000004 / z__) * .5 + z__ * log(z__) - z__ + log(.083333333333333329 / z__ + 1. + .003472222222222222 / (d__1 * d__1) - .0026813271604938273 / (d__2 * (d__2 * d__2)) - 2.2947209362139917e-4 / (d__3 * d__3));
	if (z__ == *x)
	{
		return ret_val;
	}
	ret_val = ret_val - log(z__ - 1.) - log(z__ - 2.) - log(z__ - 3.);
}

int hidata_(void)
{
	return 0;
}

int vegas_0_(int n__, int *fxn, double *avgi, double *sd, double *chi2a)
{
	static int ndmx = 50;
	static double alph = 1.5;
	static double one = 1.;
	static int mds = -1;
	//extern int fxn;

	int i__1, i__2;
	double d__1, d__2;

	double pow_dd(double *, double *);
	int pow_ii(int *, int *);
	double pow_di(double *, int *), sqrt(double);
	int s_wsle(cilist *), do_lio(int *, int *, char *, int),
		e_wsle(void);
	double log(double);

	static double d__[500];
	static int i__, j, k;
	static double r__[50], x[10], f2, fb;
	static int ia[10];
	static double di[500];
	static int kg[10], nd, ng;
	static double rc, dt[10], dr, dx[10], xn, xo, f2b, ti2;
	static int ndm;
	static double dxg;
	static int npg;
	static double xnd, xin[50], wgt, dv2g, xjac;
	static float qran[10];
	extern int aran9_(float *, int *);
	static double calls;

	static cilist io___738 = {0, 6, 0, 0, 0};

	switch (n__)
	{
	case 1:
		goto L_vegas1;
	case 2:
		goto L_vegas2;
	case 3:
		goto L_vegas3;
	}

	bveg2_1.ndo = 1;
	i__1 = bveg1_1.ndim;
	for (j = 1; j <= i__1; ++j)
	{
		bveg2_1.xi[j * 50 - 50] = one;
	}

L_vegas1:
	bveg2_1.it = 0;
	bveg2_1.si = 0.;
	bveg2_1.si2 = bveg2_1.si;
	bveg2_1.swgt = bveg2_1.si;
	bveg2_1.schi = bveg2_1.si;

L_vegas2:
	nd = ndmx;
	ng = 1;
	if (mds == 0)
	{
		goto L2;
	}
	d__1 = (double)((float)bveg1_1.ncall / 2.f);
	d__2 = (double)(1.f / (float)bveg1_1.ndim);
	ng = (int)pow_dd(&d__1, &d__2);
	mds = 1;
	if ((ng << 1) - ndmx < 0)
	{
		goto L2;
	}
	mds = -1;
	npg = ng / ndmx + 1;
	nd = ng / npg;
	ng = npg * nd;
L2:
	k = pow_ii(&ng, &bveg1_1.ndim);
	npg = bveg1_1.ncall / k;
	if (npg < 2)
	{
		npg = 2;
	}
	calls = (double)(npg * k);
	dxg = one / ng;
	d__1 = calls * pow_di(&dxg, &bveg1_1.ndim);
	dv2g = d__1 * d__1 / npg / npg / (npg - one);
	xnd = (double)nd;
	ndm = nd - 1;
	dxg *= xnd;
	xjac = one / calls;
	i__1 = bveg1_1.ndim;
	for (j = 1; j <= i__1; ++j)
	{
		dx[j - 1] = bveg1_1.xu[j - 1] - bveg1_1.xl[j - 1];
		xjac *= dx[j - 1];
	}

	if (nd == bveg2_1.ndo)
	{
		goto L8;
	}
	rc = bveg2_1.ndo / xnd;
	i__1 = bveg1_1.ndim;
	for (j = 1; j <= i__1; ++j)
	{
		k = 0;
		xn = 0.;
		dr = xn;
		i__ = k;
	L4:
		++k;
		dr += one;
		xo = xn;
		xn = bveg2_1.xi[k + j * 50 - 51];
	L5:
		if (rc > dr)
		{
			goto L4;
		}
		++i__;
		dr -= rc;
		xin[i__ - 1] = xn - (xn - xo) * dr;
		if (i__ < ndm)
		{
			goto L5;
		}
		i__2 = ndm;
		for (i__ = 1; i__ <= i__2; ++i__)
		{
			bveg2_1.xi[i__ + j * 50 - 51] = xin[i__ - 1];
		}
		bveg2_1.xi[nd + j * 50 - 51] = one;
	}
	bveg2_1.ndo = nd;

L8:

L_vegas3:
L9:
	++bveg2_1.it;
	bveg3_1.ti = 0.;
	bveg3_1.tsi = bveg3_1.ti;
	i__1 = bveg1_1.ndim;
	for (j = 1; j <= i__1; ++j)
	{
		kg[j - 1] = 1;
		i__2 = nd;
		for (i__ = 1; i__ <= i__2; ++i__)
		{
			d__[i__ + j * 50 - 51] = bveg3_1.ti;
			di[i__ + j * 50 - 51] = bveg3_1.ti;
		}
	}

L11:
	fb = 0.;
	f2b = fb;
	k = 0;
L12:
	++k;
	aran9_(qran, &bveg1_1.ndim);
	wgt = xjac;
	i__2 = bveg1_1.ndim;
	for (j = 1; j <= i__2; ++j)
	{
		xn = (double)((float)kg[j - 1] - qran[j - 1]) * dxg + one;
		ia[j - 1] = (int)xn;
		if (ia[j - 1] > 1)
		{
			goto L13;
		}
		xo = bveg2_1.xi[ia[j - 1] + j * 50 - 51];
		rc = (xn - ia[j - 1]) * xo;
		goto L14;
	L13:
		xo = bveg2_1.xi[ia[j - 1] + j * 50 - 51] - bveg2_1.xi[ia[j - 1] - 1 +
															  j * 50 - 51];
		rc = bveg2_1.xi[ia[j - 1] - 1 + j * 50 - 51] + (xn - ia[j - 1]) * xo;
	L14:
		x[j - 1] = bveg1_1.xl[j - 1] + rc * dx[j - 1];
		wgt = wgt * xo * xnd;
	}

	bveg3_1.f = wgt;
	//bveg3_1.f *= (*fxn)(x, &wgt);
	f2 = bveg3_1.f * bveg3_1.f;
	fb += bveg3_1.f;
	f2b += f2;
	i__2 = bveg1_1.ndim;
	for (j = 1; j <= i__2; ++j)
	{
		di[ia[j - 1] + j * 50 - 51] += bveg3_1.f;
		if (mds >= 0)
		{
			d__[ia[j - 1] + j * 50 - 51] += f2;
		}
	}
	if (k < npg)
	{
		goto L12;
	}

	f2b = sqrt(f2b * npg);
	f2b = (f2b - fb) * (f2b + fb);
	bveg3_1.ti += fb;
	bveg3_1.tsi += f2b;
	if (mds >= 0)
	{
		goto L18;
	}
	i__2 = bveg1_1.ndim;
	for (j = 1; j <= i__2; ++j)
	{
		d__[ia[j - 1] + j * 50 - 51] += f2b;
	}
L18:
	k = bveg1_1.ndim;
L19:
	kg[k - 1] = kg[k - 1] % ng + 1;
	if (kg[k - 1] != 1)
	{
		goto L11;
	}
	--k;
	if (k > 0)
	{
		goto L19;
	}

	bveg3_1.tsi *= dv2g;
	ti2 = bveg3_1.ti * bveg3_1.ti;
	wgt = ti2 / (bveg3_1.tsi + 1e-37);
	bveg2_1.si += bveg3_1.ti * wgt;
	bveg2_1.si2 += ti2;
	bveg2_1.swgt += wgt;
	bveg2_1.swgt += 1e-37;
	bveg2_1.si2 += 1e-37;
	bveg2_1.schi += ti2 * wgt;
	*avgi = bveg2_1.si / bveg2_1.swgt;
	*sd = bveg2_1.swgt * bveg2_1.it / bveg2_1.si2;
	*chi2a = *sd * (bveg2_1.schi / bveg2_1.swgt - *avgi * *avgi) / (double)((float)bveg2_1.it - .999f);
	*sd = sqrt(one / *sd);
	if (bveg1_1.nprn == 0)
	{
		goto L21;
	}
	bveg3_1.tsi = sqrt(bveg3_1.tsi);
L21:
	i__2 = bveg1_1.ndim;
	for (j = 1; j <= i__2; ++j)
	{
		xo = d__[j * 50 - 50];
		xn = d__[j * 50 - 49];
		d__[j * 50 - 50] = (xo + xn) / 2.;
		dt[j - 1] = d__[j * 50 - 50];
		i__1 = ndm;
		for (i__ = 2; i__ <= i__1; ++i__)
		{
			d__[i__ + j * 50 - 51] = xo + xn;
			xo = xn;
			xn = d__[i__ + 1 + j * 50 - 51];
			d__[i__ + j * 50 - 51] = (d__[i__ + j * 50 - 51] + xn) / 3.;
			dt[j - 1] += d__[i__ + j * 50 - 51];
		}
		d__[nd + j * 50 - 51] = (xn + xo) / 2.;
		dt[j - 1] += d__[nd + j * 50 - 51];
	}

	i__2 = bveg1_1.ndim;
	for (j = 1; j <= i__2; ++j)
	{
		rc = 0.;
		i__1 = nd;
		for (i__ = 1; i__ <= i__1; ++i__)
		{
			r__[i__ - 1] = 0.;
			if (dt[j - 1] >= 1e18)
			{
				s_wsle(&io___738);
				do_lio(&c__9, &c__1, "************** A SINGULARITY >1.0D18", (int)36);
				e_wsle();
			}
			if (d__[i__ + j * 50 - 51] <= 1e-18)
			{
				goto L24;
			}
			xo = dt[j - 1] / d__[i__ + j * 50 - 51];
			d__1 = (xo - one) / xo / log(xo);
			r__[i__ - 1] = pow_dd(&d__1, &alph);
		L24:
			rc += r__[i__ - 1];
		}
		rc /= xnd;
		k = 0;
		xn = 0.;
		dr = xn;
		i__ = k;
	L25:
		++k;
		dr += r__[k - 1];
		xo = xn;
		xn = bveg2_1.xi[k + j * 50 - 51];
	L26:
		if (rc > dr)
		{
			goto L25;
		}
		++i__;
		dr -= rc;
		xin[i__ - 1] = xn - (xn - xo) * dr / (r__[k - 1] + 1e-30);
		if (i__ < ndm)
		{
			goto L26;
		}
		i__1 = ndm;
		for (i__ = 1; i__ <= i__1; ++i__)
		{
			bveg2_1.xi[i__ + j * 50 - 51] = xin[i__ - 1];
		}
		bveg2_1.xi[nd + j * 50 - 51] = one;
	}

	if (bveg2_1.it < bveg1_1.itmx && bveg1_1.acc * abs(*avgi) < *sd)
	{
		goto L9;
	}
	return 0;
}

int vegas_(int *fxn, double *avgi, double *sd,
		   double *chi2a)
{
	return vegas_0_(0, fxn, avgi, sd, chi2a);
}

int vegas1_(int *fxn, double *avgi, double *sd,
			double *chi2a)
{
	return vegas_0_(1, fxn, avgi, sd, chi2a);
}

int vegas2_(int *fxn, double *avgi, double *sd,
			double *chi2a)
{
	return vegas_0_(2, fxn, avgi, sd, chi2a);
}

int vegas3_(int *fxn, double *avgi, double *sd,
			double *chi2a)
{
	return vegas_0_(3, fxn, avgi, sd, chi2a);
}

int aran9_(float *qran, int *ndim)
{
	int i__1;

	static int i__;
	extern double ranart_(int *);

	--qran;

	i__1 = *ndim;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		qran[i__] = ranart_(&sedvax_1.num1);
	}
	return 0;
}

double gauss1_(int *f, float *a, float *b, float *eps)
{
	static float const__ = 1e-12f;
	static float w[12] = {.1012285f, .222381f, .3137067f, .3623838f, .0271525f,
						  .0622535f, .0951585f, .124629f, .149596f, .1691565f, .1826034f,
						  .1894506f};
	static float x[12] = {.9602899f, .7966665f, .5255324f, .1834346f, .9894009f,
						  .944575f, .8656312f, .7554044f, .6178762f, .4580168f, .2816036f,
						  .0950125f};

	static char fmt_7[] = "(1x,\002GAUSS1....TOO HIGH ACURACY REQUIRED\002)";

	float ret_val, r__1, r__2;

	int s_wsfe(cilist *), e_wsfe(void);

	static int i__;
	static float u, y, c1, c2, s8, aa, bb, s16, delta;
	//static int s8_1,s8_2,s8_3,s16_1,s16_2,s16_3;

	static cilist io___753 = {0, 6, 0, fmt_7, 0};

	delta = const__ * (r__1 = *a - *b, dabs(r__1));
	ret_val = 0.f;
	aa = *a;
L5:
	y = *b - aa;
	if (dabs(y) <= delta)
	{
		return ret_val;
	}
L2:
	bb = aa + y;
	c1 = (aa + bb) * .5f;
	c2 = c1 - aa;
	s8 = 0.f;
	s16 = 0.f;
	for (i__ = 1; i__ <= 4; ++i__)
	{
		u = x[i__ - 1] * c2;
		r__1 = c1 + u;
		r__2 = c1 - u;
		//s8_1=(int*)(&r__1);
		//s8_2=(int*)(&r__2);
		//s8_3= s8_1+s8_2;
		s8 += w[i__ - 1] * ((r__1)*(r__2));
	}
	for (i__ = 5; i__ <= 12; ++i__)
	{
		u = x[i__ - 1] * c2;
		r__1 = c1 + u;
		r__2 = c1 - u;
		//s16_1=(int*)(&r__1);
		//s16_2=(int*)(&r__2);
		//s16_3= s8_1+s8_2;
		s16 += w[i__ - 1] * ((r__1)*(r__2));
	}
	s8 *= c2;
	s16 *= c2;
	if ((r__1 = s16 - s8, dabs(r__1)) > *eps * (dabs(s16) + 1.f))
	{
		goto L4;
	}
	ret_val += s16;
	aa = bb;
	goto L5;
L4:
	y *= .5f;
	if (dabs(y) > delta)
	{
		goto L2;
	}
	s_wsfe(&io___753);
	e_wsfe();
	ret_val = 0.f;
	return ret_val;
}

double gauss2_(int *f, float *a, float *b, float *eps)
{
	static float const__ = 1e-12f;
	static float w[12] = {.1012285f, .222381f, .3137067f, .3623838f, .0271525f,
						  .0622535f, .0951585f, .124629f, .149596f, .1691565f, .1826034f,
						  .1894506f};
	static float x[12] = {.9602899f, .7966665f, .5255324f, .1834346f, .9894009f,
						  .944575f, .8656312f, .7554044f, .6178762f, .4580168f, .2816036f,
						  .0950125f};

	static char fmt_7[] = "(1x,\002GAUSS2....TOO HIGH ACURACY REQUIRED\002)";

	float ret_val, r__1, r__2;

	int s_wsfe(cilist *), e_wsfe(void);

	static int i__;
	static float u, y, c1, c2, s8, aa, bb, s16, delta;

	static cilist io___767 = {0, 6, 0, fmt_7, 0};

	delta = const__ * (r__1 = *a - *b, dabs(r__1));
	ret_val = 0.f;
	aa = *a;
L5:
	y = *b - aa;
	if (dabs(y) <= delta)
	{
		return ret_val;
	}
L2:
	bb = aa + y;
	c1 = (aa + bb) * .5f;
	c2 = c1 - aa;
	s8 = 0.f;
	s16 = 0.f;
	for (i__ = 1; i__ <= 4; ++i__)
	{
		u = x[i__ - 1] * c2;
		r__1 = c1 + u;
		r__2 = c1 - u;
		s8 += w[i__ - 1] * ((r__1) +(r__2));
	}
	for (i__ = 5; i__ <= 12; ++i__)
	{
		u = x[i__ - 1] * c2;
		r__1 = c1 + u;
		r__2 = c1 - u;
		s16 += w[i__ - 1] * ((r__1) + (r__2));
	}
	s8 *= c2;
	s16 *= c2;
	if ((r__1 = s16 - s8, dabs(r__1)) > *eps * (dabs(s16) + 1.f))
	{
		goto L4;
	}
	ret_val += s16;
	aa = bb;
	goto L5;
L4:
	y *= .5f;
	if (dabs(y) > delta)
	{
		goto L2;
	}
	s_wsfe(&io___767);
	e_wsfe();
	ret_val = 0.f;
	return ret_val;
}

int title_(void)
{
	static char fmt_200[] = "(//10x,\002************************************"
							"**************\002/10x,\002*     |      |       _______      /  "
							"------/     *\002/10x,\002*   ----- ------     |_____|     /_/  "
							"   /       *\002/10x,\002*    |||    /        |_____|      /    "
							"/ |       *\002/10x,\002*    /| |  /_/       /_______    /_  /  "
							"  |      *\002/10x,\002*   / |     / /     /  /  / |        ----"
							"---     *\002/10x,\002*     |    / /|       /  /  |     /     | "
							"       *\002/10x,\002*     |   / /  |     /  /  _|    /   ------"
							"-     *\002/10x,\002*                                           "
							"     *\002/10x,\002*********************************************"
							"*****\002/10x,\002                      HIJING                  "
							"    \002/10x,\002       Heavy Ion Jet INteraction Generator     "
							"   \002/10x,\002                        by                      "
							"  \002/10x,\002            X. N. Wang  and  M. Gyulassy         "
							"  \002/10x,\002             Lawrence Berkeley Laboratory        "
							"   \002//)";

	int s_wsfe(cilist *), e_wsfe(void);

	static cilist io___768 = {0, 6, 0, fmt_200, 0};

	s_wsfe(&io___768);
	e_wsfe();
	return 0;
}
