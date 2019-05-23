#include "f2c.h"
#include <math.h>
struct
{
	float r__[450003];
} aa_;

#define aa_1 aa_

struct
{
	float p[450003];
} bb_;

#define bb_1 bb_

struct
{
	float e[150001];
} cc_;

#define cc_1 cc_

struct
{
	float rho[82369], rhop[82369],rhon[82369];
} dd_;

#define dd_1 dd_

struct
{
	int id[150001], lb[150001];
} ee_;

#define ee_1 ee_

struct
{
	float proper[150001];
} hh_;

#define hh_1 hh_

struct
{
	float f[2342277];
} ff_;

#define ff_1 ff_

struct
{
	float dx, dy, dz, dpx, dpy, dpz;
} gg_;

#define gg_1 gg_

struct
{
	int nstar, ndirct;
	float dir;
} input_;

#define input_1 input_

struct
{
	float prho[2009];
} pp_;

#define pp_1 pp_

struct
{
	float phrho[2401];
} qq_;

#define qq_1 qq_

struct
{
	int massr[2];
} rr_;

#define rr_1 rr_

struct
{
	int inout[20];
} ss_;

#define ss_1 ss_

struct
{
	int zta, zpr;
} zz_;

#define zz_1 zz_

struct
{
	int num;
} run_;

#define run_1 run_

struct
{
	float tkaon[7], ekaon[14007];
} kkk_;

#define kkk_1 kkk_

struct
{
	float ak[5400], speck[12600];
	int mf;
} kaon_;

#define kaon_1 kaon_

struct
{
	float xarray[1001], earray[1001];
} table_;

#define table_1 table_

struct
{
	int masspr, massta, iseed, iavoid;
	float dt;
} input1_;

#define input1_1 input1_

struct
{
	float pirho[82369];
} ddpi_;

#define ddpi_1 ddpi_

struct
{
	float pel[82369], rxy[82369];
} tt_;

#define tt_1 tt_

struct
{
	int ilab, manyb, ntmax, icoll, insys, ipot, mode, imomen, nfreq,icflow, icrho, icou, kpoten, kmul;
} input2_;

#define input2_1 input2_

struct
{
	float plab, elab, zeropt, b0, bi, bm, dencut, cycbox;
} input3_;

#define input3_1 input3_

struct
{
	float arpar1[100];
	int iapar2[50];
	float arint1[100];
	int iaint2[50];
} arprnt_;

#define arprnt_1 arprnt_

struct
{
	float pro1[150001];
} arercp_;

#define arercp_1 arercp_

struct
{
	int multi1[1];
} arerc1_;

#define arerc1_1 arerc1_

struct
{
	int ityp1[150001];
	float gx1[150001], gy1[150001],gz1[150001], ft1[150001],px1[150001], py1[150001], pz1[150001],ee1[150001], xm1[150001];
} arprc1_;

#define arprc1_1 arprc1_

struct
{
	int itimeh;
	float bimp;
} lastt_;

#define lastt_1 lastt_

struct
{
	float efrm;
	int npart1, npart2;
} snn_;

#define snn_1 snn_

struct
{
	int lblast[150001];
	float xlast[600004], plast[600004];
	int nlast;
} hbt_;

#define hbt_1 hbt_

struct
{
	int nsav, iksdcy;
} resdcy_;

#define resdcy_1 resdcy_

struct
{
	int nseed;
} rndf77_;

#define rndf77_1 rndf77_

struct
{
	float ftsv[150001], ftsvt[150001];
} ftmax_;

#define ftmax_1 ftmax_

struct
{
	float dpertt[150001], dpertp[150001], dplast[150001],dpdcy[150001], dpdpi[150001], dpt[150001], dpp1[150001],dppion[150001];
} dpert_;

#define dpert_1 dpert_

struct
{
	float hipr1[100];
	int ihpr2[50];
	float hint1[100];
	int ihnt2[50];
} hparnt_;

#define hparnt_1 hparnt_

struct
{
	int nnn;
} nn_;

#define nn_1 nn_

struct
{
	float betax, betay, betaz, gamma;
} bg_;

#define bg_1 bg_

struct
{
	float rpion[450003];
} pa_;

#define pa_1 pa_

struct
{
	float ppion[450003];
} pb_;

#define pb_1 pb_

struct
{
	float epion[150001];
} pc_;

#define pc_1 pc_

struct
{
	int lpion[150001];
} pd_;

#define pd_1 pd_

struct
{
	float propi[150001];
} pe_;

#define pe_1 pe_

struct
{
	int lb1;
	float px1, py1, pz1, em1, e1, xfnl, yfnl, zfnl, tfnl, px1n, py1n, pz1n,dp1n;
} leadng_;

#define leadng_1 leadng_

struct
{
	float tfdcy[150001], tfdpi[150001], tft[150001];
} tdecay_;

#define tdecay_1 tdecay_

struct ppbmas_1_
{
	int niso[15], nstate;
	float ppbm[30], thresh[15], weight[15];
};

#define ppbmas_1 (*(struct ppbmas_1_ *)&ppbmas_)

struct ppb1_1_
{
	float ene, factr2[6], fsum, ppinnb, s, wtot;
};

#define ppb1_1 (*(struct ppb1_1_ *)&ppb1_)

struct
{
	float pprr, ppee, pppe, rpre, xopoe, rree;
} ppmm_;

#define ppmm_1 ppmm_

struct
{
	float em2;
	int lb2;
} dpi_;

#define dpi_1 dpi_

struct
{
	int iphidcy;
} phidcy_;

#define phidcy_1 phidcy_

struct
{
	int idpert, npertd, idxsec;
} para8_;

#define para8_1 para8_

struct
{
	float bxx[82369], byy[82369],bzz[82369];
} bbb_;

#define bbb_1 bbb_

struct
{
	int iaevt, iarun;
} arevt_;

#define arevt_1 arevt_

struct
{
	int lbnn1, lbnn2, lbnd1, lbnd2, lbns1, lbns2, lbnp1, lbnp2, lbdd1,lbdd2, lbds1, lbds2, lbdp1, lbdp2, lbss1, lbss2, lbsp1, lbsp2,lbpp1, lbpp2;
} dpifsl_;

#define dpifsl_1 dpifsl_

struct
{
	float xmnn1, xmnn2, xmnd1, xmnd2, xmns1, xmns2, xmnp1, xmnp2, xmdd1, xmdd2,xmds1, xmds2, xmdp1, xmdp2, xmss1, xmss2, xmsp1, xmsp2, xmpp1,xmpp2;
} dpifsm_;

#define dpifsm_1 dpifsm_

struct
{
	float sdmel, sdmnn, sdmnd, sdmns, sdmnp, sdmdd, sdmds, sdmdp, sdmss, sdmsp,sdmpp;
} dpisig_;

#define dpisig_1 dpisig_

struct
{
	int e_1[15];
	int fill_2[1];
	float e_3[45];
	int fill_4[15];
} ppbmas_ = {1, 2, 1, 16, 16, 4, 4, 64, 4, 4, 32, 32, 4, 8, 4, {0}, .93828f, .93828f, .939457f, .93828f, .939457f, .93828f, .939457f, 1.232f, .93828f, .939457f, 1.232f, 1.232f, 1.44f, 1.44f, 1.535f, .93828f, .939457f, .939457f, 1.232f, 1.232f, 1.44f, 1.44f, 1.232f, 1.535f, 1.535f, 1.44f, 1.535f, 1.44f, 1.535f, 1.535f, 1.87656f, 1.877737f, 1.878914f, 2.17028f, 2.171457f, 2.37828f, 2.379457f, 2.464f, 2.47328f, 2.474457f, 2.672f, 2.767f, 2.88f, 2.975f, 3.07f};

struct
{
	int fill_1[1];
	float e_2[6];
	int fill_3[4];
} ppb1_ = {{0}, 0.f, 1.f, .117f, .00327f, 3.58e-5f, 1.93e-7f};

static double c_b5 = .33333333333333331;
static int c__1 = 1;
static int c__0 = 0;
static int c__2 = 2;
static float c_b102 = 5.f;
static int c__9 = 9;
static int c__3 = 3;
static float c_b118 = .6f;
static float c_b119 = .3f;
static float c_b172 = 0.f;
static float c_b173 = 1.f;
static float c_b185 = 1.232f;
static float c_b187 = 1.08f;
static float c_b191 = 1.44f;
static float c_b234 = 2.f;
static float c_b235 = -1.f;
static double c_b342 = .16666666666666666;
static double c_b343 = -.333;
static double c_b352 = .67;
static double c_b358 = .333;
static float c_b503 = .77f;
static float c_b505 = .28f;
static int c__4 = 4;
static double c_b529 = .141;
static double c_b530 = 1.9;
static float c_b561 = 2.01f;
static double c_b646 = .7;
static double c_b651 = .048383999999999983;
static double c_b652 = 1.5;
static double c_b654 = -.519;
static double c_b655 = -.534;
static double c_b656 = -.304;
static double c_b657 = -.786;
static double c_b658 = -.277;
static double c_b736 = 3.3;

int artmn_(void)
{
	static float zet[91] = {1.f, 0.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
						   0.f, 0.f, 0.f, -1.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
						   -1.f, 0.f, 1.f, 0.f, -1.f, 0.f, -1.f, 0.f, -2.f, -1.f, 0.f, 1.f, 0.f, 0.f, 0.f,
						   0.f, -1.f, 0.f, 1.f, 0.f, -1.f, 0.f, 1.f, -1.f, 0.f, 1.f, 2.f, 0.f, 1.f, 0.f,
						   1.f, 0.f, -1.f, 0.f, 1.f, 0.f, 0.f, 0.f, -1.f, 0.f, 1.f, 0.f, -1.f, 0.f, 1.f,
						   0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, -1.f, 0.f, 0.f, 0.f,
						   0.f, -1.f};

	int i__1, i__2, i__3, i__4, i__5, i__6;
	float r__1, r__2, r__3, r__4;
	double d__1;

	double pow_dd(double *, double *);
	int s_wsfe(cilist *), e_wsfe(void);
	int s_stop(char *, ftnlen);
	double sqrt(double), log(double), exp(double);
	int i_nint(float *);

	static float b;
	static int i__, j;
	static float s, b2;
	static int i0, l0, ia, ld, le, ic, ie, il, jl;
	static float ct, et[150001], bx, by;
	static int ir, lt[150001], ix, iy, iz;
	static float pt[450003], rt[450003];
	static int is, np, nt;
	static float sp, sk;
	static int lp, ln;
	static float tz;
	static int ib, ii, lb1, ld1, ld2, ld3, ld4, ln1, lp1, lp2, lp3;
	static float pr0, pr1;
	static int ln2, ln5, np1;
	static float xx1, xx2, add, ald, ale, akg;
	static int ldd;
	static float apd, bkg, den, eta, rdd, arh, aom, akn;
	static int lpd, kkk;
	static float app;
	static int lkn;
	static float rpd;
	static int npi[1];
	static float epr, rkn, asy, rpp, rpn;
	static int lpp, lpn;
	static float apn;
	static int iso;
	static float alp, aln, ddx, ddy, ddz, grp, spt, spz, drr, udt;
	static int iss, nsh;
	static float ska0, ska1, aln5, addk, facl, cden;
	static int lddk;
	static float andk;
	static int mean, ncen;
	static float rddk, bmax;
	extern int dens_(int *, int *, int *,int *);
	static int lndk;
	static float ares, adou, rndk, appk, annk;
	static int lnnk, mass;
	extern int init_(int *, int *, int *, float *,float *, float *, float *, float *, int *, int *, int *);
	static float temp[450003];
	static int irun;
	static float pzta, prot[150001], rppk, rnnk,rres;
	static int lppk;
	static float pzpr;
	static int lrho, lres, ldou, mrun;
	static float rnsg, ecor;
	static int idir;
	static float ekin, rads, zras, vols, ecms0, engs;
	extern int flow_(int *);
	static int lrun;
	static float betac, acnnd, gammc, acndn, akaon, radta;
	static int lbloc, lcnne, ikaon;
	static float esbin, radpr;
	extern int gradu_(int *, int *, int *,int *, float *, float *, float *);
	static float skaon[7];
	static int imany;
	static float gradx, grady, gradz, rcnne, rcnnd, zdist, rcndn, rcoll, rbloc,rdirt;
	static int lkaon, numnt, lcoll, lcnnd, lcndn, ldirt, lnnom;
	static float adirt, acoll, annom, rdiff;
	static int nlost;
	static float denst, radut;
	static int ipart;
	static float adecay, betata;
	static int ldecay;
	static float addrho;
	extern int tablem_(void);
	static int lomega;
	static float gammta, rdecay, alkaon, gammas;
	static int lddrho;
	static float betapr;
	extern int graduk_(int *, int *, int *, float *, float *, float *);
	static float sekaon[14007];
	extern int iarflv_(int *);
	static float gammpr, psqare;
	static int ntotal;
	extern double ranart_(int *);
	extern int invflv_(int *);
	extern int coulin_(int *, int *, int *);
	static int lkaons, outpar;
	extern int relcol_(int *, int *, int *,int *, int *, int *, int *, int *, int *,int *, int *, int *, int *, int *, int *,int *, int *, int *, int *, int *, int *,int *, int *, int *, int *, float *, float *, float *);
	static int lnnrho;
	static float annrho, sumene, etotal, atotal;
	static int nquark, nbaryn;
	static float edenst, gradxk, gradyk, gradzk, gradxn, gradyn, gradzn,gradxp, gradyp, gradzp;
	extern int gradup_(int *, int *, int *, float *, float *, float *), gradun_(int *, int *, int *, float *, float *, float *);
	static float ctlong;
	static cilist io___14 = {0, 12, 0, "(//10X,'**** FATAL ERROR: TOO MANY TEST PART. **** ')", 0};
	hbt_1.nlast = 0;
	for (i__ = 1; i__ <= 150001; ++i__)
	{
		ftmax_1.ftsv[i__ - 1] = 0.f;
		for (irun = 1; irun <= 1; ++irun)
		{
			ftmax_1.ftsvt[i__ + irun * 150001 - 150002] = 0.f;
		}
		hbt_1.lblast[i__ - 1] = 999;
		for (j = 1; j <= 4; ++j)
		{
			hbt_1.xlast[j + (i__ << 2) - 5] = 0.f;
			hbt_1.plast[j + (i__ << 2) - 5] = 0.f;
		}
	}
	tablem_();
	ikaon = 1;
	input_1.nstar = 1;
	input_1.ndirct = 0;
	input_1.dir = .02f;
	asy = .032f;
	esbin = .04f;
	kaon_1.mf = 36;
	d__1 = (double)((float)input1_1.massta);
	radta = pow_dd(&d__1, &c_b5) * 1.124f;
	d__1 = (double)((float)input1_1.masspr);
	radpr = pow_dd(&d__1, &c_b5) * 1.124f;
	zdist = radta + radpr;
	bmax = radta + radpr;
	mass = input1_1.massta + input1_1.masspr;
	ntotal = run_1.num * mass;

	if (ntotal > 150001)
	{
		s_wsfe(&io___14);
		e_wsfe();
		s_stop("", (ftnlen)0);
	}

	eta = (float)input1_1.massta * .9383f;
	pzta = 0.f;
	betata = 0.f;
	gammta = 1.f;

	epr = (float)input1_1.masspr * (input3_1.elab * .001f + .9383f);
	r__1 = epr;
	r__2 = (float)input1_1.masspr * .9383f;
	pzpr = sqrt(r__1 * r__1 - r__2 * r__2);
	betapr = pzpr / epr;
	r__1 = betapr;
	gammpr = 1.f / sqrt(1.f - r__1 * r__1);

	betac = (pzpr + pzta) / (epr + eta);
	r__1 = betac;
	gammc = 1.f / sqrt(1.f - r__1 * r__1);

	if (input2_1.insys != 0)
	{

		r__1 = epr + eta;
		r__2 = pzpr;
		s = r__1 * r__1 - r__2 * r__2;
		xx1 = log((float)input1_1.massta) * 4.f;
		xx2 = log((float)input1_1.masspr) * 4.f;
		xx1 = exp(xx1);
		xx2 = exp(xx2);
		r__1 = s;
		i__1 = input1_1.massta;
		i__2 = input1_1.masspr;
		i__3 = input1_1.massta;
		i__4 = input1_1.masspr;
		psqare = (r__1 * r__1 + (xx1 + xx2) * .77511629195947218f - s * 2.f * .88040689000000005f * (float)(i__1 * i__1 + i__2 * i__2) - (float)(i__3 * i__3 * (i__4 * i__4)) * 2.f * .77511629195947218f) / (s * 4.f);

		r__1 = (float)input1_1.massta * .9383f;
		eta = sqrt(psqare + r__1 * r__1);
		pzta = -sqrt(psqare);
		betata = pzta / eta;
		r__1 = betata;
		gammta = 1.f / sqrt(1.f - r__1 * r__1);

		r__1 = (float)input1_1.masspr * .9383f;
		epr = sqrt(psqare + r__1 * r__1);
		pzpr = sqrt(psqare);
		betapr = pzpr / epr;
		r__1 = betapr;
		gammpr = 1.f / sqrt(1.f - r__1 * r__1);
	}
	else
	{
	}
	pzta /= (float)input1_1.massta;
	pzpr /= (float)input1_1.masspr;
	ecms0 = eta + epr;
	i__1 = input2_1.manyb;
	for (imany = 1; imany <= i__1; ++imany)
	{
		if (input2_1.manyb > 1)
		{
		L111:
			bx = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
			by = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
			b2 = bx * bx + by * by;
			if (b2 > 1.f)
			{
				goto L111;
			}
			b = sqrt(b2) * (input3_1.bm - input3_1.bi) + input3_1.bi;
		}
		else
		{
			b = input3_1.b0;
		}
		coulin_(&input1_1.masspr, &input1_1.massta, &run_1.num);
		r__1 = b / 2.f;
		r__2 = input3_1.zeropt + zdist / 2.f;
		init_(&c__1, &input1_1.massta, &run_1.num, &radta, &r__1, &r__2, &pzta, &gammta, &input1_1.iseed, &mass, &input2_1.imomen);
		i__2 = input1_1.massta + 1;
		r__1 = -b / 2.f;
		r__2 = input3_1.zeropt - zdist / 2.f;
		init_(&i__2, &mass, &run_1.num, &radpr, &r__1, &r__2, &pzpr, &gammpr,
			  &input1_1.iseed, &mass, &input2_1.imomen);
		outpar = 0;
		rr_1.massr[0] = 0;
		i__2 = run_1.num;
		for (ir = 1; ir <= i__2; ++ir)
		{
			rr_1.massr[ir] = mass;
		}
		dens_(&input2_1.ipot, &mass, &run_1.num, &outpar);

		if (input2_1.icoll != -1)
		{
			i__2 = ntotal;
			for (i__ = 1; i__ <= i__2; ++i__)
			{
				ix = i_nint(&aa_1.r__[i__ * 3 - 3]);
				iy = i_nint(&aa_1.r__[i__ * 3 - 2]);
				iz = i_nint(&aa_1.r__[i__ * 3 - 1]);
				if (ix >= 20 || iy >= 20 || iz >= 24 || ix <= -20 || iy <= -20 || iz <= -24)
				{
					goto L700;
				}
				gradu_(&input2_1.ipot, &ix, &iy, &iz, &gradx, &grady, &gradz);
				bb_1.p[i__ * 3 - 3] -= input1_1.dt * .5f * gradx;
				bb_1.p[i__ * 3 - 2] -= input1_1.dt * .5f * grady;
				bb_1.p[i__ * 3 - 1] -= input1_1.dt * .5f * gradz;
			L700:;
			}
		}
		rcnne = 0.f;
		rdd = 0.f;
		rpp = 0.f;
		rppk = 0.f;
		rpn = 0.f;
		rpd = 0.f;
		rkn = 0.f;
		rnnk = 0.f;
		rddk = 0.f;
		rndk = 0.f;
		rcnnd = 0.f;
		rcndn = 0.f;
		rcoll = 0.f;
		rbloc = 0.f;
		rdirt = 0.f;
		rdecay = 0.f;
		rres = 0.f;
		for (kkk = 1; kkk <= 5; ++kkk)
		{
			skaon[kkk - 1] = 0.f;
			for (is = 1; is <= 2000; ++is)
			{
				sekaon[kkk + is * 7 - 1] = 0.f;
			}
		}
		pr0 = 0.f;
		pr1 = 0.f;
		ska0 = 0.f;
		ska1 = 0.f;
		if (arprnt_1.iapar2[0] != 1)
		{
			for (i__ = 1; i__ <= 150001; ++i__)
			{
				for (j = 1; j <= 3; ++j)
				{
					aa_1.r__[j + i__ * 3 - 4] = 0.f;
					bb_1.p[j + i__ * 3 - 4] = 0.f;
				}
				cc_1.e[i__ - 1] = 0.f;
				ee_1.lb[i__ - 1] = 0;
				ee_1.id[i__ - 1] = 0;
				hh_1.proper[i__ - 1] = 1.f;
			}
			mass = 0;
			np = 0;
			i__2 = run_1.num;
			for (j = 1; j <= i__2; ++j)
			{
				rr_1.massr[j] = 0;
				npi[j - 1] = 1;
			}
			for (i__ = 1; i__ <= 1; ++i__)
			{
				for (j = 1; j <= 150001; ++j)
				{
					rt[(j + i__ * 150001) * 3 - 450006] = 0.f;
					rt[(j + i__ * 150001) * 3 - 450005] = 0.f;
					rt[(j + i__ * 150001) * 3 - 450004] = 0.f;
					pt[(j + i__ * 150001) * 3 - 450006] = 0.f;
					pt[(j + i__ * 150001) * 3 - 450005] = 0.f;
					pt[(j + i__ * 150001) * 3 - 450004] = 0.f;
					et[j + i__ * 150001 - 150002] = 0.f;
					lt[j + i__ * 150001 - 150002] = 0;
					prot[j + i__ * 150001 - 150002] = 1.f;
				}
			}
		}
		i__2 = input2_1.ntmax;
		for (nt = 1; nt <= i__2; ++nt)
		{
			lp1 = 0;
			lp2 = 0;
			lp3 = 0;
			ld1 = 0;
			ld2 = 0;
			ld3 = 0;
			ld4 = 0;
			ln1 = 0;
			ln2 = 0;
			ln5 = 0;
			le = 0;
			lkaon = 0;
			lkaons = 0;
			if (input2_1.icoll != 1)
			{
				numnt = nt;
				relcol_(&lcoll, &lbloc, &lcnne, &ldd, &lpp, &lppk, &lpn, &lpd,&lrho, &lomega, &lkn, &lnnk, &lddk, &lndk, &lcnnd, &lcndn, &ldirt, &ldecay, &lres, &ldou, &lddrho, &lnnrho, &lnnom, &numnt, &input2_1.ntmax, &sp, &akaon,&sk);
				rcoll += (float)lcoll / run_1.num;
				rbloc += (float)lbloc / run_1.num;
				rcnne += (float)lcnne / run_1.num;
				rdd += (float)ldd / run_1.num;
				rpp += (float)lpp / run_1.num;
				rppk += (float)lppk / run_1.num;
				rpn += (float)lpn / run_1.num;
				rpd += (float)lpd / run_1.num;
				rkn += (float)lkn / run_1.num;
				rnnk += (float)lnnk / run_1.num;
				rddk += (float)lddk / run_1.num;
				rndk += (float)lndk / run_1.num;
				rcnnd += (float)lcnnd / run_1.num;
				rcndn += (float)lcndn / run_1.num;
				rdirt += (float)ldirt / run_1.num;
				rdecay += (float)ldecay / run_1.num;
				rres += (float)lres / run_1.num;
				adirt = ldirt / input1_1.dt / run_1.num;
				acoll = (lcoll - lbloc) / input1_1.dt / run_1.num;
				acnnd = lcnnd / input1_1.dt / run_1.num;
				acndn = lcndn / input1_1.dt / run_1.num;
				adecay = ldecay / input1_1.dt / run_1.num;
				ares = lres / input1_1.dt / run_1.num;
				adou = ldou / input1_1.dt / run_1.num;
				addrho = lddrho / input1_1.dt / run_1.num;
				annrho = lnnrho / input1_1.dt / run_1.num;
				annom = lnnom / input1_1.dt / run_1.num;
				add = ldd / input1_1.dt / run_1.num;
				app = lpp / input1_1.dt / run_1.num;
				appk = lppk / input1_1.dt / run_1.num;
				apn = lpn / input1_1.dt / run_1.num;
				apd = lpd / input1_1.dt / run_1.num;
				arh = lrho / input1_1.dt / run_1.num;
				aom = lomega / input1_1.dt / run_1.num;
				akn = lkn / input1_1.dt / run_1.num;
				annk = lnnk / input1_1.dt / run_1.num;
				addk = lddk / input1_1.dt / run_1.num;
				andk = lndk / input1_1.dt / run_1.num;
			}

			dens_(&input2_1.ipot, &mass, &run_1.num, &outpar);

			sumene = 0.f;
			iso = 0;
			i__3 = run_1.num;
			for (mrun = 1; mrun <= i__3; ++mrun)
			{
				iso += rr_1.massr[mrun - 1];
				i__4 = rr_1.massr[mrun];
				for (i0 = 1; i0 <= i__4; ++i0)
				{
					i__ = i0 + iso;
					r__1 = cc_1.e[i__ - 1];
					r__2 = bb_1.p[i__ * 3 - 3];
					r__3 = bb_1.p[i__ * 3 - 2];
					r__4 = bb_1.p[i__ * 3 - 1];
					etotal = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 +
								  r__4 * r__4);
					sumene += etotal;
					if (input2_1.kpoten != 0 && ee_1.lb[i__ - 1] == 23)
					{
						den = 0.f;
						ix = i_nint(&aa_1.r__[i__ * 3 - 3]);
						iy = i_nint(&aa_1.r__[i__ * 3 - 2]);
						iz = i_nint(&aa_1.r__[i__ * 3 - 1]);
						if (ix < 20 && iy < 20 && iz < 24 && ix > -20 && iy > -20 && iz > -24)
						{
							den = dd_1.rho[ix + (iy + iz * 41) * 41 + 41184];
						}
						akg = .1727f;
						bkg = .333f;
						rnsg = den;
						r__1 = bkg * den;
						ecor = -akg * rnsg + r__1 * r__1;
						r__1 = etotal;
						etotal = sqrt(r__1 * r__1 + ecor);
					}

					if (input2_1.kpoten != 0 && ee_1.lb[i__ - 1] == 21)
					{
						den = 0.f;
						ix = i_nint(&aa_1.r__[i__ * 3 - 3]);
						iy = i_nint(&aa_1.r__[i__ * 3 - 2]);
						iz = i_nint(&aa_1.r__[i__ * 3 - 1]);
						if (ix < 20 && iy < 20 && iz < 24 && ix > -20 && iy > -20 && iz > -24)
						{
							den = dd_1.rho[ix + (iy + iz * 41) * 41 + 41184];
						}
						akg = .1727f;
						bkg = .333f;
						rnsg = den;
						r__1 = bkg * den;
						ecor = -akg * rnsg + r__1 * r__1;
						r__1 = etotal;
						etotal = sqrt(r__1 * r__1 + ecor);
					}

					aa_1.r__[i__ * 3 - 3] += input1_1.dt * bb_1.p[i__ * 3 - 3] / etotal;
					aa_1.r__[i__ * 3 - 2] += input1_1.dt * bb_1.p[i__ * 3 - 2] / etotal;
					aa_1.r__[i__ * 3 - 1] += input1_1.dt * bb_1.p[i__ * 3 - 1] / etotal;
					if (input3_1.cycbox != 0.f)
					{
						if (aa_1.r__[i__ * 3 - 3] > input3_1.cycbox / 2)
						{
							aa_1.r__[i__ * 3 - 3] -= input3_1.cycbox;
						}
						if (aa_1.r__[i__ * 3 - 3] <= -input3_1.cycbox / 2)
						{
							aa_1.r__[i__ * 3 - 3] += input3_1.cycbox;
						}
						if (aa_1.r__[i__ * 3 - 2] > input3_1.cycbox / 2)
						{
							aa_1.r__[i__ * 3 - 2] -= input3_1.cycbox;
						}
						if (aa_1.r__[i__ * 3 - 2] <= -input3_1.cycbox / 2)
						{
							aa_1.r__[i__ * 3 - 2] += input3_1.cycbox;
						}
						if (aa_1.r__[i__ * 3 - 1] > input3_1.cycbox / 2)
						{
							aa_1.r__[i__ * 3 - 1] -= input3_1.cycbox;
						}
						if (aa_1.r__[i__ * 3 - 1] <= -input3_1.cycbox / 2)
						{
							aa_1.r__[i__ * 3 - 1] += input3_1.cycbox;
						}
					}
					lb1 = ee_1.lb[i__ - 1];
					if (lb1 == 9)
					{
						++ld1;
					}
					if (lb1 == 8)
					{
						++ld2;
					}
					if (lb1 == 7)
					{
						++ld3;
					}
					if (lb1 == 6)
					{
						++ld4;
					}
					if (lb1 == 11)
					{
						++ln1;
					}
					if (lb1 == 10)
					{
						++ln2;
					}
					if (lb1 == 13 || lb1 == 12)
					{
						++ln5;
					}
					if (lb1 == 0)
					{
						++le;
					}
					if (lb1 == 23)
					{
						++lkaon;
					}
					if (lb1 == 30)
					{
						++lkaons;
					}
					if (lb1 == 5)
					{
						++lp1;
					}
					if (lb1 == 4)
					{
						++lp2;
					}
					if (lb1 == 3)
					{
						++lp3;
					}
				}
			}
			lp = lp1 + lp2 + lp3;
			ld = ld1 + ld2 + ld3 + ld4;
			ln = ln1 + ln2;
			alp = (float)lp / (float)run_1.num;
			ald = (float)ld / (float)run_1.num;
			aln = (float)ln / (float)run_1.num;
			aln5 = (float)ln5 / (float)run_1.num;
			atotal = alp + ald + aln + aln5 * .5f;
			ale = (float)le / (float)run_1.num;
			alkaon = (float)lkaon / (float)run_1.num;
			if (input2_1.icou == 1)
			{
				iso = 0;
				i__4 = run_1.num;
				for (irun = 1; irun <= i__4; ++irun)
				{
					iso += rr_1.massr[irun - 1];
					i__3 = rr_1.massr[irun];
					for (il = 1; il <= i__3; ++il)
					{
						temp[il * 3 - 3] = 0.f;
						temp[il * 3 - 2] = 0.f;
						temp[il * 3 - 1] = 0.f;
					}
					i__3 = rr_1.massr[irun];
					for (il = 1; il <= i__3; ++il)
					{
						i__ = iso + il;
						if (zet[ee_1.lb[i__ - 1] + 45] != 0.f)
						{
							i__5 = il - 1;
							for (jl = 1; jl <= i__5; ++jl)
							{
								j = iso + jl;
								if (zet[ee_1.lb[j - 1] + 45] != 0.f)
								{
									ddx = aa_1.r__[i__ * 3 - 3] - aa_1.r__[j *3 -3];
									ddy = aa_1.r__[i__ * 3 - 2] - aa_1.r__[j *3 -2];
									ddz = aa_1.r__[i__ * 3 - 1] - aa_1.r__[j *3 -1];
									r__1 = ddx;
									r__2 = ddy;
									r__3 = ddz;
									rdiff = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
									if (rdiff <= 1.f)
									{
										rdiff = 1.f;
									}
									r__1 = rdiff;
									grp = zet[ee_1.lb[i__ - 1] + 45] * zet[ee_1.lb[j - 1] + 45] / (r__1 * (r__1 * r__1));
									ddx *= grp;
									ddy *= grp;
									ddz *= grp;
									temp[il * 3 - 3] += ddx;
									temp[il * 3 - 2] += ddy;
									temp[il * 3 - 1] += ddz;
									temp[jl * 3 - 3] -= ddx;
									temp[jl * 3 - 2] -= ddy;
									temp[jl * 3 - 1] -= ddz;
								}
							}
						}
					}
					i__3 = rr_1.massr[irun];
					for (il = 1; il <= i__3; ++il)
					{
						i__ = iso + il;
						if (zet[ee_1.lb[i__ - 1] + 45] != 0.f)
						{
							for (idir = 1; idir <= 3; ++idir)
							{
								bb_1.p[idir + i__ * 3 - 4] += temp[idir + il * 3 - 4] * input1_1.dt * .00144f;
							}
						}
					}
				}
			}
			spt = 0.f;
			spz = 0.f;
			ncen = 0;
			ekin = 0.f;
			nlost = 0;
			mean = 0;
			nquark = 0;
			nbaryn = 0;
			rads = 2.f;
			zras = .1f;
			denst = 0.f;
			edenst = 0.f;
			i__4 = run_1.num;
			for (irun = 1; irun <= i__4; ++irun)
			{
				mean += rr_1.massr[irun - 1];
				i__3 = rr_1.massr[irun];
				for (j = 1; j <= i__3; ++j)
				{
					i__ = j + mean;

					r__1 = aa_1.r__[i__ * 3 - 3];
					r__2 = aa_1.r__[i__ * 3 - 2];
					radut = sqrt(r__1 * r__1 + r__2 * r__2);
					if (radut <= rads)
					{
						if ((r__1 = aa_1.r__[i__ * 3 - 1], dabs(r__1)) <=
							zras * nt * input1_1.dt)
						{
							r__1 = rads;
							vols = r__1 * r__1 * 3.14159f * zras;
							r__1 = bb_1.p[i__ * 3 - 3];
							r__2 = bb_1.p[i__ * 3 - 2];
							r__3 = bb_1.p[i__ * 3 - 1];
							r__4 = cc_1.e[i__ - 1];
							engs = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
							gammas = 1.f;
							if (cc_1.e[i__ - 1] != 0.f)
							{
								gammas = engs / cc_1.e[i__ - 1];
							}
							denst += 1.f / gammas / vols;
							edenst += engs / gammas / gammas / vols;
						}
					}
					r__1 = aa_1.r__[i__ * 3 - 3];
					r__2 = aa_1.r__[i__ * 3 - 2];
					r__3 = aa_1.r__[i__ * 3 - 1];
					drr = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
					if (drr <= 2.f)
					{
						r__1 = bb_1.p[i__ * 3 - 3];
						r__2 = bb_1.p[i__ * 3 - 2];
						spt = spt + r__1 * r__1 + r__2 * r__2;
						r__1 = bb_1.p[i__ * 3 - 1];
						spz += r__1 * r__1;
						++ncen;
						r__1 = bb_1.p[i__ * 3 - 3];
						r__2 = bb_1.p[i__ * 3 - 2];
						r__3 = bb_1.p[i__ * 3 - 1];
						r__4 = cc_1.e[i__ - 1];
						ekin = ekin + sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4) - cc_1.e[i__ - 1];
					}
					ix = i_nint(&aa_1.r__[i__ * 3 - 3]);
					iy = i_nint(&aa_1.r__[i__ * 3 - 2]);
					iz = i_nint(&aa_1.r__[i__ * 3 - 1]);
					if (ix < 20 && iy < 20 && iz < 24 && ix > -20 && iy > -20 && iz > -24)
					{
						if (dd_1.rho[ix + (iy + iz * 41) * 41 + 41184] /
								.168f >
							input3_1.dencut)
						{
							goto L5800;
						}
						if (dd_1.rho[ix + (iy + iz * 41) * 41 + 41184] /
									.168f >
								5.f &&
							cc_1.e[i__ - 1] > .9f)
						{
							++nbaryn;
						}
						if (tt_1.pel[ix + (iy + iz * 41) * 41 + 41184] > 2.f)
						{
							++nquark;
						}
					}
					if (input2_1.kpoten != 0 && ee_1.lb[i__ - 1] == 23)
					{
						den = 0.f;
						if (ix < 20 && iy < 20 && iz < 24 && ix > -20 && iy > -20 && iz > -24)
						{
							den = dd_1.rho[ix + (iy + iz * 41) * 41 + 41184];
							akg = .1727f;
							bkg = .333f;
							rnsg = den;
							r__1 = bkg * den;
							ecor = -akg * rnsg + r__1 * r__1;
							r__1 = bb_1.p[i__ * 3 - 3];
							r__2 = bb_1.p[i__ * 3 - 2];
							r__3 = bb_1.p[i__ * 3 - 1];
							r__4 = cc_1.e[i__ - 1];
							etotal = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4 + ecor);
							r__1 = bkg;
							ecor = -akg + r__1 * r__1 * 2.f * den + bkg * 2.f * etotal;
							graduk_(&ix, &iy, &iz, &gradxk, &gradyk, &gradzk);
							bb_1.p[i__ * 3 - 3] -= input1_1.dt * gradxk *ecor / (etotal * 2.f);
							bb_1.p[i__ * 3 - 2] -= input1_1.dt * gradyk *ecor / (etotal * 2.f);
							bb_1.p[i__ * 3 - 1] -= input1_1.dt * gradzk *ecor / (etotal * 2.f);
						}
					}

					if (input2_1.kpoten != 0 && ee_1.lb[i__ - 1] == 21)
					{
						den = 0.f;
						if (ix < 20 && iy < 20 && iz < 24 && ix > -20 && iy > -20 && iz > -24)
						{
							den = dd_1.rho[ix + (iy + iz * 41) * 41 + 41184];
							graduk_(&ix, &iy, &iz, &gradxk, &gradyk, &gradzk);
							akg = .1727f;
							bkg = .333f;
							rnsg = den;
							r__1 = bkg * den;
							ecor = -akg * rnsg + r__1 * r__1;
							r__1 = bb_1.p[i__ * 3 - 3];
							r__2 = bb_1.p[i__ * 3 - 2];
							r__3 = bb_1.p[i__ * 3 - 1];
							r__4 = cc_1.e[i__ - 1];
							etotal = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4 + ecor);
							r__1 = bkg;
							ecor = -akg + r__1 * r__1 * 2.f * den - bkg * 2.f * etotal;
							bb_1.p[i__ * 3 - 3] -= input1_1.dt * gradxk *ecor / (etotal * 2.f);
							bb_1.p[i__ * 3 - 2] -= input1_1.dt * gradyk *ecor / (etotal * 2.f);
							bb_1.p[i__ * 3 - 1] -= input1_1.dt * gradzk *ecor / (etotal * 2.f);
						}
					}

					if (j > mass)
					{
						goto L5800;
					}
					if (input2_1.icoll != -1)
					{
						if (ix < 20 && iy < 20 && iz < 24 && ix > -20 && iy > -20 && iz > -24)
						{
							gradu_(&input2_1.ipot, &ix, &iy, &iz, &gradx, &grady, &gradz);
							tz = 0.f;
							gradxn = 0.f;
							gradyn = 0.f;
							gradzn = 0.f;
							gradxp = 0.f;
							gradyp = 0.f;
							gradzp = 0.f;
							if (input2_1.icou == 1)
							{
								gradup_(&ix, &iy, &iz, &gradxp, &gradyp, &gradzp);
								gradun_(&ix, &iy, &iz, &gradxn, &gradyn, &gradzn);
								if (zet[ee_1.lb[i__ - 1] + 45] != 0.f)
								{
									tz = -1.f;
								}
								if (zet[ee_1.lb[i__ - 1] + 45] == 0.f)
								{
									tz = 1.f;
								}
							}
							if ((i__5 = ee_1.lb[i__ - 1], abs(i__5)) >= 14 &&
								(i__6 = ee_1.lb[i__ - 1], abs(i__6)) <=
									17)
							{
								facl = .66666666666666663f;
							}
							else if ((i__5 = ee_1.lb[i__ - 1], abs(i__5)) ==
										 40 ||
									 (i__6 = ee_1.lb[i__ - 1], abs(i__6)) == 41)
							{
								facl = .33333333333333331f;
							}
							else
							{
								facl = 1.f;
							}
							bb_1.p[i__ * 3 - 3] -= facl * input1_1.dt * (gradx + asy * (gradxn - gradxp) * tz);
							bb_1.p[i__ * 3 - 2] -= facl * input1_1.dt * (grady + asy * (gradyn - gradyp) * tz);
							bb_1.p[i__ * 3 - 1] -= facl * input1_1.dt * (gradz + asy * (gradzn - gradzp) * tz);
						}
					}
				L5800:;
				}
			}
			cden = dd_1.rho[41184] / .168f;
			if (nt / input2_1.nfreq * input2_1.nfreq == nt)
			{
				if (input2_1.icflow == 1)
				{
					flow_(&nt);
				}
			}
			if (arprnt_1.iapar2[0] != 1)
			{
				ct = nt * input1_1.dt;
				ia = 0;
				i__4 = run_1.num;
				for (irun = 1; irun <= i__4; ++irun)
				{
					i__3 = rr_1.massr[irun];
					for (ic = 1; ic <= i__3; ++ic)
					{
						ie = ia + ic;
						rt[(ic + irun * 150001) * 3 - 450006] = aa_1.r__[ie * 3 -3];
						rt[(ic + irun * 150001) * 3 - 450005] = aa_1.r__[ie * 3 -2];
						rt[(ic + irun * 150001) * 3 - 450004] = aa_1.r__[ie * 3 -1];
						pt[(ic + irun * 150001) * 3 - 450006] = bb_1.p[ie * 3 - 3];
						pt[(ic + irun * 150001) * 3 - 450005] = bb_1.p[ie * 3 - 2];
						pt[(ic + irun * 150001) * 3 - 450004] = bb_1.p[ie * 3 - 1];
						et[ic + irun * 150001 - 150002] = cc_1.e[ie - 1];
						lt[ic + irun * 150001 - 150002] = ee_1.lb[ie - 1];
						prot[ic + irun * 150001 - 150002] = hh_1.proper[ie -1];
						dpert_1.dpertt[ic + irun * 150001 - 150002] = dpert_1.dpertp[ie - 1];
					}
					np = rr_1.massr[irun];
					np1 = npi[irun - 1];
					ctlong = ct;
					if (nt == input2_1.ntmax - 1)
					{
						ctlong = 1e30f;
					}
					else if (nt == input2_1.ntmax)
					{
						goto L1111;
					}
					while (np1 <= arerc1_1.multi1[irun - 1] && arprc1_1.ft1[np1 + irun * 150001 - 150002] > ct - input1_1.dt && arprc1_1.ft1[np1 + irun * 150001 - 150002] <= ctlong)
					{
						++np;
						udt = (ct - arprc1_1.ft1[np1 + irun * 150001 - 150002]) / arprc1_1.ee1[np1 + irun * 150001 - 150002];
						if (nt == input2_1.ntmax - 1)
						{
							ftmax_1.ftsvt[np + irun * 150001 - 150002] =
								arprc1_1.ft1[np1 + irun * 150001 - 150002];
							if (arprc1_1.ft1[np1 + irun * 150001 - 150002] >
								ct)
							{
								udt = 0.f;
							}
						}
						rt[(np + irun * 150001) * 3 - 450006] = arprc1_1.gx1[np1 + irun * 150001 - 150002] + arprc1_1.px1[np1 + irun * 150001 - 150002] * udt;
						rt[(np + irun * 150001) * 3 - 450005] = arprc1_1.gy1[np1 + irun * 150001 - 150002] + arprc1_1.py1[np1 + irun * 150001 - 150002] * udt;
						rt[(np + irun * 150001) * 3 - 450004] = arprc1_1.gz1[np1 + irun * 150001 - 150002] + arprc1_1.pz1[np1 + irun * 150001 - 150002] * udt;
						pt[(np + irun * 150001) * 3 - 450006] = arprc1_1.px1[np1 + irun * 150001 - 150002];
						pt[(np + irun * 150001) * 3 - 450005] = arprc1_1.py1[np1 + irun * 150001 - 150002];
						pt[(np + irun * 150001) * 3 - 450004] = arprc1_1.pz1[np1 + irun * 150001 - 150002];
						et[np + irun * 150001 - 150002] = arprc1_1.xm1[np1 +irun * 150001 - 150002];
						lt[np + irun * 150001 - 150002] = iarflv_(&arprc1_1.ityp1[np1 + irun * 150001 - 150002]);
						dpert_1.dpertt[np + irun * 150001 - 150002] = dpert_1.dpp1[np1 + irun * 150001 - 150002];
						++np1;
						prot[np + irun * 150001 - 150002] = 1.f;
					}

				L1111:
					npi[irun - 1] = np1;
					ia += rr_1.massr[irun];
					rr_1.massr[irun] = np;
				}
				ia = 0;
				i__4 = run_1.num;
				for (irun = 1; irun <= i__4; ++irun)
				{
					ia += rr_1.massr[irun - 1];
					i__3 = rr_1.massr[irun];
					for (ic = 1; ic <= i__3; ++ic)
					{
						ie = ia + ic;
						aa_1.r__[ie * 3 - 3] = rt[(ic + irun * 150001) * 3 - 450006];
						aa_1.r__[ie * 3 - 2] = rt[(ic + irun * 150001) * 3 - 450005];
						aa_1.r__[ie * 3 - 1] = rt[(ic + irun * 150001) * 3 - 450004];
						bb_1.p[ie * 3 - 3] = pt[(ic + irun * 150001) * 3 - 450006];
						bb_1.p[ie * 3 - 2] = pt[(ic + irun * 150001) * 3 - 450005];
						bb_1.p[ie * 3 - 1] = pt[(ic + irun * 150001) * 3 - 450004];
						cc_1.e[ie - 1] = et[ic + irun * 150001 - 150002];
						ee_1.lb[ie - 1] = lt[ic + irun * 150001 - 150002];
						hh_1.proper[ie - 1] = prot[ic + irun * 150001 - 150002];
						if (nt == input2_1.ntmax - 1)
						{
							ftmax_1.ftsv[ie - 1] = ftmax_1.ftsvt[ic + irun * 150001 - 150002];
						}
						dpert_1.dpertp[ie - 1] = dpert_1.dpertt[ic + irun * 150001 - 150002];
					}
				}
			}
		}
		iss = 0;
		i__2 = run_1.num;
		for (lrun = 1; lrun <= i__2; ++lrun)
		{
			iss += rr_1.massr[lrun - 1];
			i__4 = rr_1.massr[lrun];
			for (l0 = 1; l0 <= i__4; ++l0)
			{
				ipart = iss + l0;
			}
		}
		if (arprnt_1.iapar2[0] != 1)
		{
			ia = 0;
			i__2 = run_1.num;
			for (irun = 1; irun <= i__2; ++irun)
			{
				ia += rr_1.massr[irun - 1];
				np1 = npi[irun - 1];
				nsh = rr_1.massr[irun] - np1 + 1;
				arerc1_1.multi1[irun - 1] += nsh;
				if (nsh > 0)
				{
					ib = arerc1_1.multi1[irun - 1];
					ie = rr_1.massr[irun] + 1;
					ii = -1;
				}
				else if (nsh < 0)
				{
					ib = rr_1.massr[irun] + 1;
					ie = arerc1_1.multi1[irun - 1];
					ii = 1;
				}
				if (nsh != 0)
				{
					i__4 = ie;
					i__3 = ii;
					for (i__ = ib; i__3 < 0 ? i__ >= i__4 : i__ <= i__4; i__ += i__3)
					{
						j = i__ - nsh;
						arprc1_1.ityp1[i__ + irun * 150001 - 150002] =arprc1_1.ityp1[j + irun * 150001 - 150002];
						arprc1_1.gx1[i__ + irun * 150001 - 150002] =arprc1_1.gx1[j + irun * 150001 - 150002];
						arprc1_1.gy1[i__ + irun * 150001 - 150002] =arprc1_1.gy1[j + irun * 150001 - 150002];
						arprc1_1.gz1[i__ + irun * 150001 - 150002] =arprc1_1.gz1[j + irun * 150001 - 150002];
						arprc1_1.ft1[i__ + irun * 150001 - 150002] =arprc1_1.ft1[j + irun * 150001 - 150002];
						arprc1_1.px1[i__ + irun * 150001 - 150002] =arprc1_1.px1[j + irun * 150001 - 150002];
						arprc1_1.py1[i__ + irun * 150001 - 150002] =arprc1_1.py1[j + irun * 150001 - 150002];
						arprc1_1.pz1[i__ + irun * 150001 - 150002] =arprc1_1.pz1[j + irun * 150001 - 150002];
						arprc1_1.ee1[i__ + irun * 150001 - 150002] =arprc1_1.ee1[j + irun * 150001 - 150002];
						arprc1_1.xm1[i__ + irun * 150001 - 150002] =arprc1_1.xm1[j + irun * 150001 - 150002];
						arercp_1.pro1[i__ + irun * 150001 - 150002] =arercp_1.pro1[j + irun * 150001 - 150002];
						dpert_1.dpp1[i__ + irun * 150001 - 150002] =dpert_1.dpp1[j + irun * 150001 - 150002];
					}
				}
				i__3 = rr_1.massr[irun];
				for (i__ = 1; i__ <= i__3; ++i__)
				{
					ib = ia + i__;
					arprc1_1.ityp1[i__ + irun * 150001 - 150002] = invflv_(&ee_1.lb[ib - 1]);
					arprc1_1.gx1[i__ + irun * 150001 - 150002] = aa_1.r__[ib *
																			  3 -
																		  3];
					arprc1_1.gy1[i__ + irun * 150001 - 150002] = aa_1.r__[ib *
																			  3 -
																		  2];
					arprc1_1.gz1[i__ + irun * 150001 - 150002] = aa_1.r__[ib *
																			  3 -
																		  1];
					if (arprc1_1.ft1[i__ + irun * 150001 - 150002] < ct)
					{
						arprc1_1.ft1[i__ + irun * 150001 - 150002] = ct;
					}
					arprc1_1.px1[i__ + irun * 150001 - 150002] = bb_1.p[ib *
																			3 -
																		3];
					arprc1_1.py1[i__ + irun * 150001 - 150002] = bb_1.p[ib *
																			3 -
																		2];
					arprc1_1.pz1[i__ + irun * 150001 - 150002] = bb_1.p[ib *
																			3 -
																		1];
					arprc1_1.xm1[i__ + irun * 150001 - 150002] = cc_1.e[ib -
																		1];
					r__1 = arprc1_1.px1[i__ + irun * 150001 - 150002];
					r__2 = arprc1_1.py1[i__ + irun * 150001 - 150002];
					r__3 = arprc1_1.pz1[i__ + irun * 150001 - 150002];
					r__4 = arprc1_1.xm1[i__ + irun * 150001 - 150002];
					arprc1_1.ee1[i__ + irun * 150001 - 150002] = sqrt(r__1 *
																		  r__1 +
																	  r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
					arercp_1.pro1[i__ + irun * 150001 - 150002] = hh_1.proper[ib - 1];
				}
			}
		}
	}

	return 0;
}

int coulin_(int *masspr, int *massta, int *num)
{
	int i__1, i__2;

	static int i__, mass, irun;

	mass = *massta + *masspr;
	i__1 = *num;
	for (irun = 1; irun <= i__1; ++irun)
	{
		i__2 = zz_1.zta + (irun - 1) * mass;
		for (i__ = (irun - 1) * mass + 1; i__ <= i__2; ++i__)
		{
			ee_1.lb[i__ - 1] = 1;
		}
		i__2 = *massta + (irun - 1) * mass;
		for (i__ = zz_1.zta + 1 + (irun - 1) * mass; i__ <= i__2; ++i__)
		{
			ee_1.lb[i__ - 1] = 2;
		}
		i__2 = *massta + zz_1.zpr + (irun - 1) * mass;
		for (i__ = *massta + 1 + (irun - 1) * mass; i__ <= i__2; ++i__)
		{
			ee_1.lb[i__ - 1] = 1;
		}
		i__2 = *massta + *masspr + (irun - 1) * mass;
		for (i__ = *massta + zz_1.zpr + 1 + (irun - 1) * mass; i__ <= i__2;
			 ++i__)
		{
			ee_1.lb[i__ - 1] = 2;
		}
	}
	return 0;
}

int relcol_(int *lcoll, int *lbloc, int *lcnne,
			int *ldd, int *lpp, int *lppk, int *lpn, int *lpd,
			int *lrho, int *lomega, int *lkn, int *lnnk, int *lddk, int *lndk, int *lcnnd, int *lcndn, int *ldirt,
			int *ldecay, int *lres, int *ldou, int *lddrho,
			int *lnnrho, int *lnnom, int *nt, int *ntmax, float *sp, float *akaon, float *sk)
{
	static float zet[91] = {1.f, 0.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
						   0.f, 0.f, 0.f, -1.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
						   -1.f, 0.f, 1.f, 0.f, -1.f, 0.f, -1.f, 0.f, -2.f, -1.f, 0.f, 1.f, 0.f, 0.f, 0.f,
						   0.f, -1.f, 0.f, 1.f, 0.f, -1.f, 0.f, 1.f, -1.f, 0.f, 1.f, 2.f, 0.f, 1.f, 0.f,
						   1.f, 0.f, -1.f, 0.f, 1.f, 0.f, 0.f, 0.f, -1.f, 0.f, 1.f, 0.f, -1.f, 0.f, 1.f,
						   0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, -1.f, 0.f, 0.f, 0.f,
						   0.f, -1.f};

	int i__1, i__2, i__3, i__4, i__5, i__6, i__7;
	float r__1, r__2, r__3, r__4, r__5, r__6, r__7, r__8;

	double sqrt(double), exp(double);
	int i_nint(float *), s_wsle(cilist *), do_lio(int *, int *, char *, ftnlen), e_wsle(void);
	double log(double);

	static int i__, j, n;
	static float e2;
	static int i0, i1, j1, j2, i2, n0;
	static float t0, x1, y1, z1, x2, y2, z2, ec;
	static int ia, j10, ic, ib, ie, ig, il, im;
	static float ds, et[150001];
	static int kp, lt[150001], is;
	static float pt[450003], rt[450003];
	static int in;
	static float xx, yy, ec0;
	extern int sdbelastic_(float *, float *);
	static int id1;
	static float am1;
	static int id2;
	static float am2, er1, er2, pk0;
	extern double pp1_(float *);
	static int ix1, iy1, iz1, ix2, iy2, iz2;
	extern double pp2_(float *);
	static float px2, py2, pz2, xx0, ece, sdb;
	static int lbm;
	static float dse;
	extern double w1440_(float *);
	extern int cms_(int *, int *, float *, float *,float *, float *);
	static float wid;
	extern double w1535_(float *);
	static float sig, pcx, pcy, pcz, srt;
	static int iss;
	extern int xnd_(float *, float *, float *, float *, int *, int *, float *, float *, float *, float *, float *, float *, float *);
	static float xmm;
	extern double ppt_(float *), xpp_(float *), xnp_(float *);
	static float dsr, sdm;
	static int ipp;
	static float ert;
	static int ikk;
	static float e1cm;
	static int ilb1, ilb2, lb1i, lb2i;
	static float em1i, em2i, e2cm;
	static int lbp1;
	static float emm1, emm2;
	static int lbp2, ipx1, ipy1, ipz1;
	static float px1i, py1i, pz1i, px2i, py2i, pz2i;
	static int ipx2, ipy2, ipz2;
	static float xsk1, xsk2, xsk3, xsk4, xsk5, xsk6, xsk7, xky1, xky2;
	extern double reab_(int *, int *, float *, int *);
	static float xky3, xky4;
	extern int crdd_(int *, float *, float *, float *, float *, int *, int *, int *, int *, float *, float *,
					 int *, int *);
	static float xky5;
	extern int crhb_(float *, float *, float *, float *, int *, int *, int *);
	static float xky6, xky7, eini, brel;
	extern int crnd_(int *, float *, float *, float *, float *, int *, int *, int *, float *, float *, float *, float *, float *, float *, float *, float *, int *, int *), cren_(
																																															float *, float *, float *, float *, int *, int *, int *);
	static int ntag;
	extern int crrd_(float *, float *, float *, float *, int *, int *, int *, float *, float *, float *, float *), crpd_(
																														   float *, float *, float *, float *, int *, int *, int *,
																														   float *, float *, float *, float *);
	static float sigk;
	static int lpdr;
	static float prot[150001];
	static int msum, irun, mass;
	static float xeta;
	extern double xn1535_(int *, int *, int *);
	static float xphi, xmax;
	extern int wida1_(float *, float *, float *, int *);
	extern double xnpi_(int *, int *, int *, float *);
	static float xres;
	extern int crpn_(float *, float *, float *, float *, int *, int *, int *, float *, float *, float *, float *), crnn_(
																														   int *, float *, float *, float *, float *, int *, int *,
																														   int *, int *, float *, float *, int *, int *);
	static float ppel;
	extern int ppxs_(int *, int *, float *, float *,
					 float *, int *);
	static float ppin, dspp;
	extern int crpp_(float *, float *, float *, float *, int *, int *, int *, float *, float *, float *, int *);
	static float dskk, dskk0, sigp, pzrt;
	static int ikkg;
	static float sigr0;
	static int ikkl;
	static float dshn, xky8, px1cm, py1cm, pz1cm;
	extern double pipp1_(float *);
	static float xky9, xky10, xky11, xky12, xky13, xky14, xky15, xky16, xky17;
	static int ikmp;
	static float dskn;
	extern int crkn_(float *, float *, float *, float *, int *, int *, int *);
	static float pt1i1, pt2i1, pt3i1, pt1i2, pt2i2, pt3i2;
	static int icase;
	extern int decay_(int *, int *, int *,
					  int *, float *, int *),
		sbbdm_(float *, float *, int *,
			   int *, float *, float *),
		sdmbb_(float *, float *, int *);
	extern double akpel_(float *), width_(float *);
	static float drmax;
	static int lrhor;
	static float rhomp, pxini, pyini, pzini, spipi, bmass, pkaon;
	extern int decay2_(int *, int *, int *,
					   int *, float *, int *);
	extern double aknel_(float *);
	static float brsgm, brsig;
	static int nchrg;
	extern int newka_(int *, int *, int *, float *, int *, int *, int *, int *, float *, float *,
					  float *, float *, int *);
	static int ictrl;
	extern double pnlka_(float *), pnska_(float *);
	static float sigma0, xkaon, xphin, xmaxn, xnpin;
	extern int distc0_(float *, float *, float *, int *,
					   float *, float *, float *, float *, float *, float *, float *, float *,
					   float *, float *, float *, float *, float *, float *, float *, float *,
					   float *);
	extern double dirct1_(float *), dirct2_(float *), dirct3_(float *);
	static float deltr0, dr0max, xnpid, xreab;
	extern double dpion_(float *, float *, int *, int *, float *);
	extern int crdir_(float *, float *, float *, float *,
					  int *, int *, int *);
	static float xkaon0, spika;
	extern double erhon_(float *, float *, int *, int *, float *);
	static float signn, xinel;
	static int ianti, ipert1;
	static float signn0;
	extern int xddin_(float *, float *, float *, float *,
					  int *, int *, float *, float *, float *, float *, float *,
					  float *, float *);
	static float ppsig, ppink;
	extern double pipik_(float *);
	static float xmaxn1, xnpin1;
	extern int spprr_(int *, int *, float *), sppee_(int *, int *, float *), spppe_(int *, int *, float *), srpre_(int *, int *, float *), sopoe_(int *, int *, float *), srree_(int *, int *, float *);
	static float dsppr, dsppb;
	static int icheck;
	extern int crlaba_(float *, float *, float *, float *, float *, float *, int *, int *, int *, int *, int *,
					   int *),
		xphib_(int *, int *, float *, float *, float *,
			   float *, float *, float *, float *, float *, float *),
		crdmbb_(float *,
				float *, float *, float *, int *, int *, int *, int *, float *, int *, int *);
	static float dshnr;
	extern int crdbel_(float *, float *, float *, float *,
					   int *, int *, int *, int *, float *, int *,
					   int *);
	static float sigkp;
	extern int crlan_(float *, float *, float *, float *,
					  int *, int *, int *);
	static int idecay;
	extern int lambar_(int *, int *, float *, float *);
	static int ichann;
	static float siglab;
	static int iblock;
	static float pdecay, gfactr;
	extern int resdec_(int *, int *, int *, float *, int *);
	static float sigela;
	extern double akplam_(float *), aknlam_(float *);
	static float xdecay;
	static int inewka;
	extern int inidcy_(void);
	static float dptemp[150001], fttemp[150001];
	static int massrn[2];
	static float resona, ftpisv[150001];
	static int nodelt, lomgar;
	extern double ranart_(int *);
	static float pnstar;
	static int nnnini;
	static float rppmax, rsqare;
	static int ifirst;
	static float sigsgm;
	extern double akpsgm_(float *), aknsgm_(float *);
	static float deltre;
	extern int distce_(int *, int *, float *, float *,
					   float *, float *, float *, int *, float *, float *, float *),
		pibphi_(float *, int *, int *, float *, float *, float *,
				float *);
	extern double pionpp_(float *);
	static float xdirct, sumsrt, xnelas, deltar;
	extern int dreson_(int *, int *), crkpla_(float *, float *, float *, float *, float *, float *, float *, float *, int *, int *, int *, int *, int *, float *), ksreso_(int *, int *);
	static float xelstc, cutoff, sdprod, pfinal, dspert, spprho;
	extern int getnst_(float *);
	extern double ppbbar_(float *), prbbar_(float *), rrbbar_(float *),
		pobbar_(float *), robbar_(float *), oobbar_(float *), xppbar_(float *);
	static float dsppbr;
	extern int crppba_(float *, float *, float *, float *,
					   int *, int *, int *),
		pertur_(float *, float *, float *,
				float *, int *, int *, int *, int *, int *,
				int *);
	static int icontp;
	extern int crphib_(float *, float *, float *, float *,
					   int *, int *, float *, float *, float *, float *, float *,
					   float *, int *),
		phimes_(int *, int *, float *, float *,
				float *, float *, float *, float *, float *, float *, float *);
	static float sigphi;
	extern int crphim_(float *, float *, float *, float *,
					   int *, int *, float *, float *, float *, float *, float *,
					   float *, float *, int *, int *, int *),
		xkhype_(int
					*,
				int *, float *, float *, float *, float *, float *, float *,
				float *, float *, float *, float *, float *, float *, float *, float *,
				float *, float *, float *, float *, float *),
		crkhyp_(float *, float *,
				float *, float *, int *, int *, float *, float *, float *,
				float *, float *, float *, float *, float *, float *, float *, float *,
				float *, float *, float *, float *, float *, float *, float *, int *,
				int *),
		crkphi_(float *, float *, float *, float *, float *,
				int *, float *, float *, int *, int *, int *,
				int *, int *, int *, float *, float *),
		crksph_(float *,
				float *, float *, float *, float *, float *, float *, int *,
				int *, int *, int *, int *, int *, int *,
				int *, float *);
	static float dsknr, p1beta, transf;
	extern int hbtout_(int *, int *, int *);
	static int ipdflag;
	static float dsrpert;

	static cilist io___379 = {0, 6, 0, 0, 0};

	inidcy_();
	resona = 5.f;
	nodelt = 0;
	sumsrt = 0.f;
	*lcoll = 0;
	*lbloc = 0;
	*lcnne = 0;
	*ldd = 0;
	*lpp = 0;
	*lpd = 0;
	lpdr = 0;
	*lrho = 0;
	lrhor = 0;
	*lomega = 0;
	lomgar = 0;
	*lpn = 0;
	*lkn = 0;
	*lnnk = 0;
	*lddk = 0;
	*lndk = 0;
	*lppk = 0;
	*lcnnd = 0;
	*lcndn = 0;
	*ldirt = 0;
	*ldecay = 0;
	*lres = 0;
	*ldou = 0;
	*lddrho = 0;
	*lnnrho = 0;
	*lnnom = 0;
	msum = 0;
	massrn[0] = 0;
	for (il = 1; il <= 5; ++il)
	{
		kkk_1.tkaon[il - 1] = 0.f;
		for (is = 1; is <= 2000; ++is)
		{
			kkk_1.ekaon[il + is * 7 - 1] = 0.f;
		}
	}
	i__1 = run_1.num;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		for (j = 1; j <= 150001; ++j)
		{
			pe_1.propi[j + i__ * 150001 - 150002] = 1.f;
		}
	}
	for (i__ = 1; i__ <= 150001; ++i__)
	{
		fttemp[i__ - 1] = 0.f;
		for (irun = 1; irun <= 1; ++irun)
		{
			ftpisv[i__ + irun * 150001 - 150002] = 0.f;
		}
	}
	*sp = 0.f;
	*akaon = 0.f;
	*sk = 0.f;
	mass = 0;
	i__1 = run_1.num;
	for (irun = 1; irun <= i__1; ++irun)
	{
		nn_1.nnn = 0;
		msum += rr_1.massr[irun - 1];
		j10 = 2;
		if (*nt == *ntmax)
		{
			j10 = 1;
		}

		i__2 = rr_1.massr[irun];
		for (j1 = j10; j1 <= i__2; ++j1)
		{
			i1 = j1 + msum;
			if (cc_1.e[i1 - 1] == 0.f)
			{
				goto L800;
			}
			if (ee_1.lb[i1 - 1] < -45 || ee_1.lb[i1 - 1] > 45)
			{
				goto L800;
			}
			x1 = aa_1.r__[i1 * 3 - 3];
			y1 = aa_1.r__[i1 * 3 - 2];
			z1 = aa_1.r__[i1 * 3 - 1];
			leadng_1.px1 = bb_1.p[i1 * 3 - 3];
			leadng_1.py1 = bb_1.p[i1 * 3 - 2];
			leadng_1.pz1 = bb_1.p[i1 * 3 - 1];
			leadng_1.em1 = cc_1.e[i1 - 1];
			am1 = leadng_1.em1;
			r__1 = leadng_1.em1;
			r__2 = leadng_1.px1;
			r__3 = leadng_1.py1;
			r__4 = leadng_1.pz1;
			leadng_1.e1 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
			id1 = ee_1.id[i1 - 1];
			leadng_1.lb1 = ee_1.lb[i1 - 1];
			if (*nt == *ntmax && (leadng_1.lb1 == 21 || leadng_1.lb1 == 23))
			{
				pk0 = ranart_(&rndf77_1.nseed);
				if (pk0 < .25f)
				{
					ee_1.lb[i1 - 1] = 22;
				}
				else if (pk0 < .5f)
				{
					ee_1.lb[i1 - 1] = 24;
				}
				leadng_1.lb1 = ee_1.lb[i1 - 1];
			}
			if (leadng_1.lb1 == 0 || leadng_1.lb1 == 25 || leadng_1.lb1 == 26 || leadng_1.lb1 == 27 || leadng_1.lb1 == 28 ||
				leadng_1.lb1 == 29 || abs(leadng_1.lb1) == 30 || abs(leadng_1.lb1) >= 6 && abs(leadng_1.lb1) <= 13 ||
				resdcy_1.iksdcy == 1 && leadng_1.lb1 == 24 || abs(leadng_1.lb1) == 16)
			{
			}
			else
			{
				goto L1;
			}
			if (leadng_1.lb1 >= 25 && leadng_1.lb1 <= 27)
			{
				wid = .151f;
			}
			else if (leadng_1.lb1 == 28)
			{
				wid = .00841f;
			}
			else if (leadng_1.lb1 == 29)
			{
				wid = .00443f;
			}
			else if (abs(leadng_1.lb1) == 30)
			{
				wid = .051f;
			}
			else if (leadng_1.lb1 == 0)
			{
				wid = 1.18e-6f;
			}
			else if (resdcy_1.iksdcy == 1 && leadng_1.lb1 == 24)
			{
				wid = 7.36e-15f;
			}
			else if (abs(leadng_1.lb1) == 16)
			{
				wid = 8.87e-6f;
			}
			else if (leadng_1.lb1 == 32)
			{
				wida1_(&leadng_1.em1, &rhomp, &wid, &input1_1.iseed);
			}
			else if (abs(leadng_1.lb1) >= 6 && abs(leadng_1.lb1) <= 9)
			{
				wid = width_(&leadng_1.em1);
			}
			else if (abs(leadng_1.lb1) == 10 || abs(leadng_1.lb1) == 11)
			{
				wid = w1440_(&leadng_1.em1);
			}
			else if (abs(leadng_1.lb1) == 12 || abs(leadng_1.lb1) == 13)
			{
				wid = w1535_(&leadng_1.em1);
			}
			if (*nt == *ntmax)
			{
				pdecay = 1.1f;
				if (phidcy_1.iphidcy == 0 && abs(leadng_1.lb1) == 29)
				{
					pdecay = 0.f;
				}
			}
			else
			{
				t0 = .19733f / wid;
				gfactr = leadng_1.e1 / leadng_1.em1;
				t0 *= gfactr;
				if (t0 > 0.f)
				{
					pdecay = 1.f - exp(-input1_1.dt / t0);
				}
				else
				{
					pdecay = 0.f;
				}
			}
			xdecay = ranart_(&rndf77_1.nseed);
			if (xdecay < pdecay)
			{
				idecay = irun;
				leadng_1.tfnl = *nt * input1_1.dt;
				if (*nt == *ntmax && ftmax_1.ftsv[i1 - 1] > (*ntmax - 1) *
																input1_1.dt)
				{
					leadng_1.tfnl = ftmax_1.ftsv[i1 - 1];
				}
				leadng_1.xfnl = x1;
				leadng_1.yfnl = y1;
				leadng_1.zfnl = z1;
				if (leadng_1.lb1 == 0 || leadng_1.lb1 == 25 || leadng_1.lb1 == 26 || leadng_1.lb1 == 27 || leadng_1.lb1 == 28 ||
					leadng_1.lb1 == 29 || abs(leadng_1.lb1) == 30 || abs(leadng_1.lb1) >= 6 && abs(leadng_1.lb1) <= 9 ||
					resdcy_1.iksdcy == 1 && leadng_1.lb1 == 24 || abs(leadng_1.lb1) == 16)
				{
					resdec_(&i1, nt, &nn_1.nnn, &wid, &idecay);
					bb_1.p[i1 * 3 - 3] = leadng_1.px1n;
					bb_1.p[i1 * 3 - 2] = leadng_1.py1n;
					bb_1.p[i1 * 3 - 1] = leadng_1.pz1n;
					dpert_1.dpertp[i1 - 1] = leadng_1.dp1n;
					if (*nt == *ntmax)
					{
						aa_1.r__[i1 * 3 - 3] = leadng_1.xfnl;
						aa_1.r__[i1 * 3 - 2] = leadng_1.yfnl;
						aa_1.r__[i1 * 3 - 1] = leadng_1.zfnl;
						tdecay_1.tfdcy[i1 - 1] = leadng_1.tfnl;
					}

					if (abs(leadng_1.lb1) >= 6 && abs(leadng_1.lb1) <= 9)
					{
						++(*ldecay);
					}
				}
				else if (abs(leadng_1.lb1) == 10 || abs(leadng_1.lb1) == 11)
				{
					++nn_1.nnn;
					++(*ldecay);
					pnstar = 1.f;
					if (cc_1.e[i1 - 1] > 1.22f)
					{
						pnstar = .6f;
					}
					if (ranart_(&rndf77_1.nseed) <= pnstar)
					{
						decay_(&idecay, &i1, &nn_1.nnn, &input1_1.iseed, &wid,
							   nt);
					}
					else
					{
						decay2_(&idecay, &i1, &nn_1.nnn, &input1_1.iseed, &wid, nt);
						++nn_1.nnn;
					}
				}
				else if (abs(leadng_1.lb1) == 12 || abs(leadng_1.lb1) == 13)
				{
					++nn_1.nnn;
					decay_(&idecay, &i1, &nn_1.nnn, &input1_1.iseed, &wid, nt);
					++(*ldecay);
				}

				if (*nt == *ntmax)
				{
					if (ee_1.lb[i1 - 1] == 25 || ee_1.lb[i1 - 1] == 26 ||
						ee_1.lb[i1 - 1] == 27)
					{
						wid = .151f;
					}
					else if (ee_1.lb[i1 - 1] == 0)
					{
						wid = 1.18e-6f;
					}
					else if (ee_1.lb[i1 - 1] == 24 && resdcy_1.iksdcy == 1)
					{
						wid = 7.36e-17f;
					}
					else
					{
						goto L9000;
					}
					leadng_1.lb1 = ee_1.lb[i1 - 1];
					leadng_1.px1 = bb_1.p[i1 * 3 - 3];
					leadng_1.py1 = bb_1.p[i1 * 3 - 2];
					leadng_1.pz1 = bb_1.p[i1 * 3 - 1];
					leadng_1.em1 = cc_1.e[i1 - 1];
					r__1 = leadng_1.em1;
					r__2 = leadng_1.px1;
					r__3 = leadng_1.py1;
					r__4 = leadng_1.pz1;
					leadng_1.e1 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
					resdec_(&i1, nt, &nn_1.nnn, &wid, &idecay);
					bb_1.p[i1 * 3 - 3] = leadng_1.px1n;
					bb_1.p[i1 * 3 - 2] = leadng_1.py1n;
					bb_1.p[i1 * 3 - 1] = leadng_1.pz1n;
					aa_1.r__[i1 * 3 - 3] = leadng_1.xfnl;
					aa_1.r__[i1 * 3 - 2] = leadng_1.yfnl;
					aa_1.r__[i1 * 3 - 1] = leadng_1.zfnl;
					tdecay_1.tfdcy[i1 - 1] = leadng_1.tfnl;
					dpert_1.dpertp[i1 - 1] = leadng_1.dp1n;
				}
			L9000:
				goto L800;
			}
		L1:
			if (*nt == *ntmax)
			{
				goto L800;
			}
			x1 = aa_1.r__[i1 * 3 - 3];
			y1 = aa_1.r__[i1 * 3 - 2];
			z1 = aa_1.r__[i1 * 3 - 1];

			i__3 = j1 - 1;
			for (j2 = 1; j2 <= i__3; ++j2)
			{
				i2 = j2 + msum;
				if (cc_1.e[i2 - 1] == 0.f)
				{
					goto L600;
				}
				if (cc_1.e[i1 - 1] == 0.f)
				{
					goto L800;
				}
				if (ee_1.lb[i2 - 1] < -45 || ee_1.lb[i2 - 1] > 45)
				{
					goto L600;
				}
				x2 = aa_1.r__[i2 * 3 - 3];
				y2 = aa_1.r__[i2 * 3 - 2];
				z2 = aa_1.r__[i2 * 3 - 1];
				dr0max = 5.f;
				ilb1 = (i__4 = ee_1.lb[i1 - 1], abs(i__4));
				ilb2 = (i__4 = ee_1.lb[i2 - 1], abs(i__4));
				if (ilb1 == 42 || ilb2 == 42)
				{
					if (ilb1 >= 1 && ilb1 <= 2 || ilb1 >= 6 && ilb1 <= 13 ||
						ilb2 >= 1 && ilb2 <= 2 || ilb2 >= 6 && ilb2 <= 13)
					{
						if (ee_1.lb[i1 - 1] * ee_1.lb[i2 - 1] > 0)
						{
							dr0max = 10.f;
						}
					}
				}

				r__1 = x1 - x2;
				r__2 = y1 - y2;
				r__3 = z1 - z2;
				r__4 = dr0max;
				if (r__1 * r__1 + r__2 * r__2 + r__3 * r__3 > r__4 * r__4)
				{
					goto L600;
				}
				if (ee_1.id[i1 - 1] * ee_1.id[i2 - 1] == input1_1.iavoid)
				{
					goto L400;
				}
				id1 = ee_1.id[i1 - 1];
				id2 = ee_1.id[i2 - 1];

				r__1 = x1 / gg_1.dx;
				ix1 = i_nint(&r__1);
				r__1 = y1 / gg_1.dy;
				iy1 = i_nint(&r__1);
				r__1 = z1 / gg_1.dz;
				iz1 = i_nint(&r__1);
				leadng_1.px1 = bb_1.p[i1 * 3 - 3];
				leadng_1.py1 = bb_1.p[i1 * 3 - 2];
				leadng_1.pz1 = bb_1.p[i1 * 3 - 1];
				leadng_1.em1 = cc_1.e[i1 - 1];
				am1 = leadng_1.em1;
				leadng_1.lb1 = ee_1.lb[i1 - 1];
				r__1 = leadng_1.em1;
				r__2 = leadng_1.px1;
				r__3 = leadng_1.py1;
				r__4 = leadng_1.pz1;
				leadng_1.e1 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 +
								   r__4 * r__4);
				r__1 = leadng_1.px1 / gg_1.dpx;
				ipx1 = i_nint(&r__1);
				r__1 = leadng_1.py1 / gg_1.dpy;
				ipy1 = i_nint(&r__1);
				r__1 = leadng_1.pz1 / gg_1.dpz;
				ipz1 = i_nint(&r__1);
				dpi_1.lb2 = ee_1.lb[i2 - 1];
				px2 = bb_1.p[i2 * 3 - 3];
				py2 = bb_1.p[i2 * 3 - 2];
				pz2 = bb_1.p[i2 * 3 - 1];
				dpi_1.em2 = cc_1.e[i2 - 1];
				am2 = dpi_1.em2;
				lb1i = ee_1.lb[i1 - 1];
				lb2i = ee_1.lb[i2 - 1];
				px1i = bb_1.p[i1 * 3 - 3];
				py1i = bb_1.p[i1 * 3 - 2];
				pz1i = bb_1.p[i1 * 3 - 1];
				em1i = cc_1.e[i1 - 1];
				px2i = bb_1.p[i2 * 3 - 3];
				py2i = bb_1.p[i2 * 3 - 2];
				pz2i = bb_1.p[i2 * 3 - 1];
				em2i = cc_1.e[i2 - 1];
				r__1 = cc_1.e[i1 - 1];
				r__2 = bb_1.p[i1 * 3 - 3];
				r__3 = bb_1.p[i1 * 3 - 2];
				r__4 = bb_1.p[i1 * 3 - 1];
				r__5 = cc_1.e[i2 - 1];
				r__6 = bb_1.p[i2 * 3 - 3];
				r__7 = bb_1.p[i2 * 3 - 2];
				r__8 = bb_1.p[i2 * 3 - 1];
				eini = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4) + sqrt(r__5 * r__5 + r__6 * r__6 + r__7 * r__7 + r__8 * r__8);
				pxini = bb_1.p[i1 * 3 - 3] + bb_1.p[i2 * 3 - 3];
				pyini = bb_1.p[i1 * 3 - 2] + bb_1.p[i2 * 3 - 2];
				pzini = bb_1.p[i1 * 3 - 1] + bb_1.p[i2 * 3 - 1];
				nnnini = nn_1.nnn;

				iblock = 0;

				deltr0 = 3.f;
				if (abs(leadng_1.lb1) >= 14 && abs(leadng_1.lb1) <= 17 || abs(
																			  leadng_1.lb1) >= 30 &&
																			  abs(leadng_1.lb1) <= 45)
				{
					deltr0 = 5.f;
				}
				if (abs(dpi_1.lb2) >= 14 && abs(dpi_1.lb2) <= 17 || abs(
																		dpi_1.lb2) >= 30 &&
																		abs(dpi_1.lb2) <= 45)
				{
					deltr0 = 5.f;
				}
				if (leadng_1.lb1 == 28 && dpi_1.lb2 == 28)
				{
					deltr0 = 4.84f;
				}
				if (leadng_1.lb1 >= 3 && leadng_1.lb1 <= 5 && (dpi_1.lb2 >= 3 && dpi_1.lb2 <= 5))
				{
					r__1 = dpi_1.em2;
					r__2 = px2;
					r__3 = py2;
					r__4 = pz2;
					e2 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
					r__1 = leadng_1.e1 + e2;
					r__2 = leadng_1.px1 + px2;
					r__3 = leadng_1.py1 + py2;
					r__4 = leadng_1.pz1 + pz2;
					spipi = r__1 * r__1 - r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
					if (spipi >= 2.3715999999999999f)
					{
						deltr0 = 3.5f;
					}
				}
				if (leadng_1.lb1 == 23 && (dpi_1.lb2 >= 14 && dpi_1.lb2 <= 17))
				{
					goto L3699;
				}
				if (dpi_1.lb2 == 23 && (leadng_1.lb1 >= 14 && leadng_1.lb1 <=
																  17))
				{
					goto L3699;
				}
				if (leadng_1.lb1 == 21 && dpi_1.lb2 == 23)
				{
					goto L3699;
				}
				if (dpi_1.lb2 == 21 && leadng_1.lb1 == 23)
				{
					goto L3699;
				}
				if (leadng_1.lb1 == 30 && dpi_1.lb2 == 21)
				{
					goto L3699;
				}
				if (dpi_1.lb2 == 30 && leadng_1.lb1 == 21)
				{
					goto L3699;
				}
				if (leadng_1.lb1 == -30 && dpi_1.lb2 == 23)
				{
					goto L3699;
				}
				if (dpi_1.lb2 == -30 && leadng_1.lb1 == 23)
				{
					goto L3699;
				}
				if (leadng_1.lb1 == -30 && dpi_1.lb2 == 30)
				{
					goto L3699;
				}
				if (dpi_1.lb2 == -30 && leadng_1.lb1 == 30)
				{
					goto L3699;
				}

				if (leadng_1.lb1 == 21 || leadng_1.lb1 == 23)
				{
					if (dpi_1.lb2 == 0 || dpi_1.lb2 >= 25 && dpi_1.lb2 <= 28)
					{
						goto L3699;
					}
				}
				else if (dpi_1.lb2 == 21 || dpi_1.lb2 == 23)
				{
					if (leadng_1.lb1 == 0 || leadng_1.lb1 >= 25 &&
												 leadng_1.lb1 <= 28)
					{
						goto L3699;
					}
				}
				if (abs(leadng_1.lb1) == 30 && (dpi_1.lb2 == 0 || dpi_1.lb2 >= 25 && dpi_1.lb2 <= 28 || dpi_1.lb2 >= 3 && dpi_1.lb2 <= 5))
				{
					goto L3699;
				}
				else if (abs(dpi_1.lb2) == 30 && (leadng_1.lb1 == 0 ||
												  leadng_1.lb1 >= 25 && leadng_1.lb1 <= 28 ||
												  leadng_1.lb1 >= 3 && leadng_1.lb1 <= 5))
				{
					goto L3699;
				}
				else if (abs(leadng_1.lb1) == 30 && (abs(dpi_1.lb2) == 1 ||
													 abs(dpi_1.lb2) == 2 || abs(dpi_1.lb2) >= 6 && abs(dpi_1.lb2) <= 13))
				{
					goto L3699;
				}
				if (abs(dpi_1.lb2) == 30 && (abs(leadng_1.lb1) == 1 || abs(leadng_1.lb1) == 2 || abs(leadng_1.lb1) >= 6 && abs(leadng_1.lb1) <= 13))
				{
					goto L3699;
				}
				if ((leadng_1.lb1 == 23 || leadng_1.lb1 == 21) && (abs(
																	   dpi_1.lb2) == 1 ||
																   abs(dpi_1.lb2) == 2 || abs(dpi_1.lb2) >= 6 && abs(dpi_1.lb2) <= 13))
				{
					goto L3699;
				}
				else if ((dpi_1.lb2 == 23 || dpi_1.lb2 == 21) && (abs(
																	  leadng_1.lb1) == 1 ||
																  abs(leadng_1.lb1) == 2 || abs(leadng_1.lb1) >= 6 && abs(leadng_1.lb1) <= 13))
				{
					goto L3699;
				}

				rppmax = 3.57f;
				if ((leadng_1.lb1 == -1 || leadng_1.lb1 == -2 || leadng_1.lb1 >= -13 && leadng_1.lb1 <= -6) && (dpi_1.lb2 == 1 ||
																												dpi_1.lb2 == 2 || dpi_1.lb2 >= 6 && dpi_1.lb2 <= 13))
				{
					deltr0 = rppmax;
					goto L2699;
				}
				else if ((dpi_1.lb2 == -1 || dpi_1.lb2 == -2 || dpi_1.lb2 >= -13 && dpi_1.lb2 <= -6) && (leadng_1.lb1 == 1 ||
																										 leadng_1.lb1 == 2 || leadng_1.lb1 >= 6 && leadng_1.lb1 <= 13))
				{
					deltr0 = rppmax;
					goto L2699;
				}
				if (abs(leadng_1.lb1) >= 14 && abs(leadng_1.lb1) <= 17 || abs(
																			  dpi_1.lb2) >= 14 &&
																			  abs(dpi_1.lb2) <= 17)
				{
					goto L3699;
				}

				if (abs(leadng_1.lb1) == 42 || abs(dpi_1.lb2) == 42)
				{
					ilb1 = abs(leadng_1.lb1);
					ilb2 = abs(dpi_1.lb2);
					if (ilb1 >= 1 && ilb1 <= 2 || ilb1 >= 6 && ilb1 <= 13 ||
						ilb2 >= 1 && ilb2 <= 2 || ilb2 >= 6 && ilb2 <= 13)
					{
						if (leadng_1.lb1 * dpi_1.lb2 > 0)
						{
							deltr0 = 9.5f;
						}
					}
				}

				if (abs(leadng_1.lb1) >= 40 && abs(leadng_1.lb1) <= 45 || abs(
																			  dpi_1.lb2) >= 40 &&
																			  abs(dpi_1.lb2) <= 45)
				{
					goto L3699;
				}

				if (leadng_1.lb1 == 29 && (dpi_1.lb2 >= 1 && dpi_1.lb2 <= 13 || dpi_1.lb2 >= 21 && dpi_1.lb2 <= 28 || abs(dpi_1.lb2) == 30) || dpi_1.lb2 == 29 && (leadng_1.lb1 >= 1 && leadng_1.lb1 <= 13 || leadng_1.lb1 >= 21 && leadng_1.lb1 <= 28 || abs(leadng_1.lb1) == 30))
				{
					deltr0 = 3.f;
					goto L3699;
				}

				if (abs(leadng_1.lb1) == 30 || abs(dpi_1.lb2) == 30)
				{
					goto L400;
				}
				if (leadng_1.lb1 == 23 && (dpi_1.lb2 < 1 || dpi_1.lb2 > 17))
				{
					goto L400;
				}
				if (dpi_1.lb2 == 23 && (leadng_1.lb1 < 1 || leadng_1.lb1 > 17))
				{
					goto L400;
				}

				if (leadng_1.lb1 <= -1 && leadng_1.lb1 >= -13 && (dpi_1.lb2 == 0 || dpi_1.lb2 >= 3 && dpi_1.lb2 <= 5 || dpi_1.lb2 >= 25 && dpi_1.lb2 <= 28) || dpi_1.lb2 <= -1 &&
																																								   dpi_1.lb2 >= -13 && (leadng_1.lb1 == 0 || leadng_1.lb1 >= 3 && leadng_1.lb1 <= 5 || leadng_1.lb1 >= 25 && leadng_1.lb1 <= 28))
				{
				}
				else if ((leadng_1.lb1 == -1 || leadng_1.lb1 == -2) && (dpi_1.lb2 < -5 && dpi_1.lb2 >= -13) || (dpi_1.lb2 ==
																													-1 ||
																												dpi_1.lb2 == -2) &&
																												   (leadng_1.lb1 < -5 &&
																													leadng_1.lb1 >= -13))
				{
				}
				else if ((leadng_1.lb1 == -1 || leadng_1.lb1 == -2) && (dpi_1.lb2 == -1 || dpi_1.lb2 == -2))
				{
				}
				else if (leadng_1.lb1 < -5 && leadng_1.lb1 >= -13 && (dpi_1.lb2 < -5 && dpi_1.lb2 >= -13))
				{
				}
			L2699:
				if (leadng_1.lb1 == 1 || leadng_1.lb1 == 2 || leadng_1.lb1 >= 6 && leadng_1.lb1 <= 17)
				{
					if (dpi_1.lb2 == 1 || dpi_1.lb2 == 2 || dpi_1.lb2 >= 6 && dpi_1.lb2 <= 17)
					{
						deltr0 = 2.f;
					}
				}

			L3699:
				r__1 = x1 - x2;
				r__2 = y1 - y2;
				r__3 = z1 - z2;
				rsqare = r__1 * r__1 + r__2 * r__2 + r__3 * r__3;
				r__1 = deltr0;
				if (rsqare > r__1 * r__1)
				{
					goto L400;
				}
				r__1 = x2 / gg_1.dx;
				ix2 = i_nint(&r__1);
				r__1 = y2 / gg_1.dy;
				iy2 = i_nint(&r__1);
				r__1 = z2 / gg_1.dz;
				iz2 = i_nint(&r__1);
				r__1 = px2 / gg_1.dpx;
				ipx2 = i_nint(&r__1);
				r__1 = py2 / gg_1.dpy;
				ipy2 = i_nint(&r__1);
				r__1 = pz2 / gg_1.dpz;
				ipz2 = i_nint(&r__1);
				cms_(&i1, &i2, &pcx, &pcy, &pcz, &srt);
				drmax = dr0max;
				distc0_(&drmax, &deltr0, &input1_1.dt, &ifirst, &pcx, &pcy, &pcz, &x1, &y1, &z1, &leadng_1.px1, &leadng_1.py1, &leadng_1.pz1, &leadng_1.em1, &x2, &y2, &z2, &px2, &py2, &pz2, &dpi_1.em2);
				if (ifirst == -1)
				{
					goto L400;
				}
				r__1 = srt / .04f;
				iss = i_nint(&r__1);
				if (iss > 2000)
				{
					iss = 2000;
				}
				if (abs(leadng_1.lb1) == 42 || abs(dpi_1.lb2) == 42)
				{
					ilb1 = abs(leadng_1.lb1);
					ilb2 = abs(dpi_1.lb2);
					if (leadng_1.lb1 == 0 || leadng_1.lb1 >= 3 && leadng_1.lb1 <= 5 || leadng_1.lb1 >= 25 && leadng_1.lb1 <= 28 || dpi_1.lb2 == 0 || dpi_1.lb2 >= 3 && dpi_1.lb2 <= 5 || dpi_1.lb2 >= 25 && dpi_1.lb2 <= 28)
					{
						goto L505;
					}
					else if ((ilb1 >= 1 && ilb1 <= 2 || ilb1 >= 6 && ilb1 <= 13 || ilb2 >= 1 && ilb2 <= 2 || ilb2 >= 6 && ilb2 <= 13) && leadng_1.lb1 * dpi_1.lb2 > 0)
					{
						goto L506;
					}
					else
					{
						goto L400;
					}
				}

				if ((leadng_1.lb1 == 23 || leadng_1.lb1 == 30) && (dpi_1.lb2 == -1 || dpi_1.lb2 == -2 || dpi_1.lb2 >= -13 && dpi_1.lb2 <= -6) || (dpi_1.lb2 == 23 || dpi_1.lb2 ==
																																										 30) &&
																																					 (leadng_1.lb1 == -1 || leadng_1.lb1 == -2 ||
																																					  leadng_1.lb1 >= -13 && leadng_1.lb1 <= -6))
				{
					bmass = .938f;
					if (srt <= bmass + .498f)
					{
						pkaon = 0.f;
					}
					else
					{
						r__2 = srt;
						r__3 = bmass;
						r__1 = (r__2 * r__2 - (r__3 * r__3 + .248004f)) / 2.f / bmass;
						pkaon = sqrt(r__1 * r__1 - .248004f);
					}
					sigela = (akpel_(&pkaon) + aknel_(&pkaon)) * .5f;
					sigsgm = akpsgm_(&pkaon) * 1.5f + aknsgm_(&pkaon);
					sig = sigela + sigsgm + akplam_(&pkaon);
					if (sig > 1e-7f)
					{
						icase = 3;
						brel = sigela / sig;
						brsgm = sigsgm / sig;
						brsig = sig;
						nchrg = 1;
						goto L3555;
					}
					goto L400;
				}

				if (leadng_1.lb1 >= -17 && leadng_1.lb1 <= -14 && (dpi_1.lb2 >= 3 && dpi_1.lb2 <= 5) || dpi_1.lb2 >= -17 &&
																											dpi_1.lb2 <= -14 && (leadng_1.lb1 >= 3 && leadng_1.lb1 <= 5))
				{
					nchrg = -100;
					if (leadng_1.lb1 == -15 && (dpi_1.lb2 == 5 || dpi_1.lb2 ==
																	  27) ||
						dpi_1.lb2 == -15 && (leadng_1.lb1 == 5 ||
											 leadng_1.lb1 == 27))
					{
						nchrg = -2;
						bmass = 1.232f;
						goto L110;
					}
					if (leadng_1.lb1 == -15 && (dpi_1.lb2 == 0 || dpi_1.lb2 == 4 || dpi_1.lb2 == 26 || dpi_1.lb2 == 28) ||
						dpi_1.lb2 == -15 && (leadng_1.lb1 == 0 ||
											 leadng_1.lb1 == 4 || leadng_1.lb1 == 26 ||
											 leadng_1.lb1 == 28) ||
						(leadng_1.lb1 == -14 ||
						 leadng_1.lb1 == -16) &&
							(dpi_1.lb2 == 5 ||
							 dpi_1.lb2 == 27) ||
						(dpi_1.lb2 == -14 ||
						 dpi_1.lb2 == -16) &&
							(leadng_1.lb1 == 5 ||
							 leadng_1.lb1 == 27))
					{
						nchrg = -1;
						bmass = .938f;
						goto L110;
					}
					if (leadng_1.lb1 == -15 && (dpi_1.lb2 == 3 || dpi_1.lb2 ==
																	  25) ||
						dpi_1.lb2 == -15 && (leadng_1.lb1 == 3 ||
											 leadng_1.lb1 == 25) ||
						leadng_1.lb1 == -17 && (dpi_1.lb2 == 5 || dpi_1.lb2 == 27) || dpi_1.lb2 == -17 && (leadng_1.lb1 == 5 || leadng_1.lb1 == 27) || (leadng_1.lb1 == -14 || leadng_1.lb1 == -16) && (dpi_1.lb2 == 0 || dpi_1.lb2 == 4 || dpi_1.lb2 == 26 || dpi_1.lb2 == 28) || (dpi_1.lb2 == -14 || dpi_1.lb2 == -16) && (leadng_1.lb1 == 0 || leadng_1.lb1 == 4 || leadng_1.lb1 == 26 || leadng_1.lb1 == 28))
					{
						nchrg = 0;
						bmass = .938f;
						goto L110;
					}
					if (leadng_1.lb1 == -17 && (dpi_1.lb2 == 0 || dpi_1.lb2 == 4 || dpi_1.lb2 == 26 || dpi_1.lb2 == 28) ||
						dpi_1.lb2 == -17 && (leadng_1.lb1 == 0 ||
											 leadng_1.lb1 == 4 || leadng_1.lb1 == 26 ||
											 leadng_1.lb1 == 28) ||
						(leadng_1.lb1 == -14 ||
						 leadng_1.lb1 == -16) &&
							(dpi_1.lb2 == 3 ||
							 dpi_1.lb2 == 25) ||
						(dpi_1.lb2 == -14 ||
						 dpi_1.lb2 == -16) &&
							(leadng_1.lb1 == 3 ||
							 leadng_1.lb1 == 25))
					{
						nchrg = 1;
						bmass = 1.232f;
					}

				L110:
					sig = 0.f;
					if (nchrg != -100 && srt >= bmass + .498f)
					{
						icase = 4;
						r__2 = srt;
						r__1 = (r__2 * r__2 - 1.1278479999999997f) / 2.f /
							   .938f;
						pkaon = sqrt(r__1 * r__1 - .248004f);
						if (leadng_1.lb1 == -14 || dpi_1.lb2 == -14)
						{
							if (nchrg >= 0)
							{
								sigma0 = akplam_(&pkaon);
							}
							if (nchrg < 0)
							{
								sigma0 = aknlam_(&pkaon);
							}
						}
						else
						{
							if (nchrg >= 0)
							{
								sigma0 = akpsgm_(&pkaon);
							}
							if (nchrg < 0)
							{
								sigma0 = aknsgm_(&pkaon);
							}
							sigma0 = akpsgm_(&pkaon) * 1.5f + aknsgm_(&pkaon);
						}
						r__1 = srt;
						r__2 = bmass + .498f;
						r__3 = srt;
						r__4 = .498f - bmass;
						r__5 = srt;
						r__6 = leadng_1.em1 + dpi_1.em2;
						r__7 = srt;
						r__8 = leadng_1.em1 - dpi_1.em2;
						sig = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / (r__5 * r__5 - r__6 * r__6) / (r__7 * r__7 - r__8 * r__8) * sigma0;
						if (nchrg == -2 || nchrg == 2)
						{
							sig *= 2.f;
						}
						if (leadng_1.lb1 == -14 || dpi_1.lb2 == -14)
						{
							sig *= 1.3333333333333333f;
						}
						else if (nchrg == -2 || nchrg == 2)
						{
							sig *= .88888888888888884f;
						}
						else
						{
							sig *= .44444444444444442f;
						}
					}
					icase = 4;
					sigela = 10.f;
					sig += sigela;
					brel = sigela / sig;
					brsgm = 0.f;
					brsig = sig;
					goto L3555;
				}
				if ((leadng_1.lb1 == 21 || leadng_1.lb1 == -30) && (dpi_1.lb2 >= 14 && dpi_1.lb2 <= 17) || (dpi_1.lb2 == 21 ||
																											dpi_1.lb2 == -30) &&
																											   (leadng_1.lb1 >= 14 &&
																												leadng_1.lb1 <= 17))
				{
					kp = 0;
					goto L3455;
				}
				if ((leadng_1.lb1 == 23 || leadng_1.lb1 == 30) && (dpi_1.lb2 <= -14 && dpi_1.lb2 >= -17) || (dpi_1.lb2 == 23 ||
																											 dpi_1.lb2 == 30) &&
																												(leadng_1.lb1 <= -14 &&
																												 leadng_1.lb1 >= -17))
				{
					kp = 1;
					goto L3455;
				}
				if ((leadng_1.lb1 == 21 || leadng_1.lb1 == -30) && (dpi_1.lb2 == 40 || dpi_1.lb2 == 41) || (dpi_1.lb2 == 21 ||
																											dpi_1.lb2 == -30) &&
																											   (leadng_1.lb1 == 40 ||
																												leadng_1.lb1 == 41))
				{
					kp = 0;
					goto L3455;
				}
				if ((leadng_1.lb1 == 23 || leadng_1.lb1 == 30) && (dpi_1.lb2 == -40 || dpi_1.lb2 == -41) || (dpi_1.lb2 == 23 ||
																											 dpi_1.lb2 == 30) &&
																												(leadng_1.lb1 == -40 ||
																												 leadng_1.lb1 == -41))
				{
					kp = 1;
					goto L3455;
				}
				kp = 3;
				if ((leadng_1.lb1 >= 3 && leadng_1.lb1 <= 5 || leadng_1.lb1 ==
																   0) &&
						(abs(dpi_1.lb2) == 40 || abs(dpi_1.lb2) == 41) ||
					(dpi_1.lb2 >= 3 && dpi_1.lb2 <= 5 || dpi_1.lb2 ==
															 0) &&
						(abs(leadng_1.lb1) == 40 || abs(leadng_1.lb1) ==
														41))
				{
					goto L3455;
				}
				if (leadng_1.lb1 >= 3 && leadng_1.lb1 <= 5 && abs(dpi_1.lb2) == 45 || dpi_1.lb2 >= 3 && dpi_1.lb2 <= 5 && abs(leadng_1.lb1) == 45)
				{
					goto L3455;
				}

				if (leadng_1.lb1 == 23 && (dpi_1.lb2 >= 14 && dpi_1.lb2 <= 17))
				{
					goto L5699;
				}
				if (dpi_1.lb2 == 23 && (leadng_1.lb1 >= 14 && leadng_1.lb1 <=
																  17))
				{
					goto L5699;
				}
				if (leadng_1.lb1 == 21 && (dpi_1.lb2 >= -17 && dpi_1.lb2 <=
																   -14))
				{
					goto L5699;
				}
				if (dpi_1.lb2 == 21 && (leadng_1.lb1 >= -17 && leadng_1.lb1 <=
																   -14))
				{
					goto L5699;
				}
				if ((leadng_1.lb1 == 1 || leadng_1.lb1 == 2 || leadng_1.lb1 >= 6 && leadng_1.lb1 <= 13) && (dpi_1.lb2 >= -17 &&
																											dpi_1.lb2 <= -14) ||
					(dpi_1.lb2 == 1 || dpi_1.lb2 == 2 || dpi_1.lb2 >= 6 && dpi_1.lb2 <= 13) && (leadng_1.lb1 >= -17 && leadng_1.lb1 <= -14))
				{
					goto L5999;
				}
				if ((leadng_1.lb1 == -1 || leadng_1.lb1 == -2 || leadng_1.lb1 <= -6 && leadng_1.lb1 >= -13) && (dpi_1.lb2 >= 14 &&
																												dpi_1.lb2 <= 17) ||
					(dpi_1.lb2 == -1 || dpi_1.lb2 == -2 || dpi_1.lb2 <= -6 && dpi_1.lb2 >= -13) && (leadng_1.lb1 >= 14 && leadng_1.lb1 <= 17))
				{
					goto L5999;
				}

				if (leadng_1.lb1 == 21 && dpi_1.lb2 == 23)
				{
					goto L8699;
				}
				if (dpi_1.lb2 == 21 && leadng_1.lb1 == 23)
				{
					goto L8699;
				}
				if (leadng_1.lb1 == 30 && dpi_1.lb2 == 21)
				{
					goto L8699;
				}
				if (dpi_1.lb2 == 30 && leadng_1.lb1 == 21)
				{
					goto L8699;
				}
				if (leadng_1.lb1 == -30 && dpi_1.lb2 == 23)
				{
					goto L8699;
				}
				if (dpi_1.lb2 == -30 && leadng_1.lb1 == 23)
				{
					goto L8699;
				}
				if (leadng_1.lb1 == -30 && dpi_1.lb2 == 30)
				{
					goto L8699;
				}
				if (dpi_1.lb2 == -30 && leadng_1.lb1 == 30)
				{
					goto L8699;
				}
				if ((leadng_1.lb1 == 23 || leadng_1.lb1 == 21 || abs(leadng_1.lb1) == 30) && (dpi_1.lb2 >= 25 && dpi_1.lb2 <= 28) || (dpi_1.lb2 == 23 || dpi_1.lb2 == 21 || abs(dpi_1.lb2) == 30) && (leadng_1.lb1 >= 25 &&
																																																	  leadng_1.lb1 <= 28))
				{
					goto L8799;
				}

				if (abs(leadng_1.lb1) == 30 && (dpi_1.lb2 >= 3 && dpi_1.lb2 <=
																	  5) ||
					abs(dpi_1.lb2) == 30 && (leadng_1.lb1 >= 3 &&
											 leadng_1.lb1 <= 5))
				{
					goto L8799;
				}

				if (leadng_1.lb1 == 29 && (dpi_1.lb2 == 1 || dpi_1.lb2 == 2 ||
										   dpi_1.lb2 >= 6 && dpi_1.lb2 <= 9) ||
					dpi_1.lb2 == 29 && (leadng_1.lb1 == 1 || leadng_1.lb1 == 2 ||
										leadng_1.lb1 >= 6 && leadng_1.lb1 <= 9))
				{
					goto L7222;
				}

				if (leadng_1.lb1 == 29 && (dpi_1.lb2 >= 3 && dpi_1.lb2 <= 5 ||
										   dpi_1.lb2 >= 21 && dpi_1.lb2 <= 28 || abs(dpi_1.lb2) == 30) ||
					dpi_1.lb2 == 29 && (leadng_1.lb1 >= 3 &&
											leadng_1.lb1 <= 5 ||
										leadng_1.lb1 >= 21 &&
											leadng_1.lb1 <= 28 ||
										abs(leadng_1.lb1) == 30))
				{
					goto L7444;
				}

				if ((abs(leadng_1.lb1) >= 14 && abs(leadng_1.lb1) <= 17 ||
					 abs(leadng_1.lb1) >= 40) &&
					(dpi_1.lb2 >= 25 &&
						 dpi_1.lb2 <= 29 ||
					 dpi_1.lb2 == 0))
				{
					goto L888;
				}
				if ((abs(dpi_1.lb2) >= 14 && abs(dpi_1.lb2) <= 17 || abs(
																		 dpi_1.lb2) >= 40) &&
					(leadng_1.lb1 >= 25 &&
						 leadng_1.lb1 <= 29 ||
					 leadng_1.lb1 == 0))
				{
					goto L888;
				}

				if ((leadng_1.lb1 == 23 || leadng_1.lb1 == 30) && (dpi_1.lb2 == 1 || dpi_1.lb2 == 2 || dpi_1.lb2 >= 6 && dpi_1.lb2 <= 13) || (dpi_1.lb2 == 23 || dpi_1.lb2 == 30) && (leadng_1.lb1 == 1 || leadng_1.lb1 == 2 ||
																																													  leadng_1.lb1 >= 6 && leadng_1.lb1 <= 13))
				{
					goto L888;
				}
				if ((leadng_1.lb1 == 21 || leadng_1.lb1 == -30) && (dpi_1.lb2 == -1 || dpi_1.lb2 == -2 || dpi_1.lb2 >= -13 && dpi_1.lb2 <= -6) || (dpi_1.lb2 == 21 || dpi_1.lb2 ==
																																										  -30) &&
																																					  (leadng_1.lb1 == -1 || leadng_1.lb1 == -2 ||
																																					   leadng_1.lb1 >= -13 && leadng_1.lb1 <= -6))
				{
					goto L888;
				}

				if (leadng_1.lb1 >= 14 && leadng_1.lb1 <= 17 && (dpi_1.lb2 >= 6 && dpi_1.lb2 <= 13) || dpi_1.lb2 >= 14 && dpi_1.lb2 <= 17 && (leadng_1.lb1 >= 6 && leadng_1.lb1 <= 13))
				{
					goto L7799;
				}
				if (leadng_1.lb1 <= -14 && leadng_1.lb1 >= -17 && (dpi_1.lb2 <= -6 && dpi_1.lb2 >= -13) || dpi_1.lb2 <= -14 &&
																											   dpi_1.lb2 >= -17 && (leadng_1.lb1 <= -6 && leadng_1.lb1 >= -13))
				{
					goto L7799;
				}

				if (abs(leadng_1.lb1) >= 40 || abs(dpi_1.lb2) >= 40 ||
					leadng_1.lb1 <= -14 && leadng_1.lb1 >= -17 ||
					dpi_1.lb2 <= -14 && dpi_1.lb2 >= -17)
				{
					goto L400;
				}

				if ((leadng_1.lb1 == -1 || leadng_1.lb1 == -2 || leadng_1.lb1 >= -13 && leadng_1.lb1 <= -6) && (dpi_1.lb2 == 1 ||
																												dpi_1.lb2 == 2 || dpi_1.lb2 >= 6 && dpi_1.lb2 <= 13))
				{
					goto L2799;
				}
				else if ((dpi_1.lb2 == -1 || dpi_1.lb2 == -2 || dpi_1.lb2 >= -13 && dpi_1.lb2 <= -6) && (leadng_1.lb1 == 1 ||
																										 leadng_1.lb1 == 2 || leadng_1.lb1 >= 6 && leadng_1.lb1 <= 13))
				{
					goto L2799;
				}

				inewka = irun;
				newka_(&icase, &inewka, &input1_1.iseed, &input1_1.dt, nt, &ictrl, &i1, &i2, &srt, &pcx, &pcy, &pcz, &iblock);
				if (ictrl == 1)
				{
					goto L400;
				}

				if (abs(leadng_1.lb1) >= 14 && abs(leadng_1.lb1) <= 17 || abs(
																			  dpi_1.lb2) >= 14 &&
																			  abs(dpi_1.lb2) <= 17)
				{
					goto L400;
				}
				if (leadng_1.lb1 >= 3 && leadng_1.lb1 <= 5 && (dpi_1.lb2 >= 3 && dpi_1.lb2 <= 5))
				{
					goto L777;
				}
				if (leadng_1.lb1 == 0 && (dpi_1.lb2 >= 3 && dpi_1.lb2 <= 5))
				{
					goto L777;
				}
				if (dpi_1.lb2 == 0 && (leadng_1.lb1 >= 3 && leadng_1.lb1 <= 5))
				{
					goto L777;
				}
				if (leadng_1.lb1 == 0 && dpi_1.lb2 == 0)
				{
					goto L777;
				}
				if (leadng_1.lb1 >= 25 && leadng_1.lb1 <= 28 && (dpi_1.lb2 >= 25 && dpi_1.lb2 <= 28))
				{
					goto L777;
				}
				if (leadng_1.lb1 >= 25 && leadng_1.lb1 <= 28 && (dpi_1.lb2 >= 3 && dpi_1.lb2 <= 5))
				{
					goto L777;
				}
				if (dpi_1.lb2 >= 25 && dpi_1.lb2 <= 28 && (leadng_1.lb1 >= 3 && leadng_1.lb1 <= 5))
				{
					goto L777;
				}
				if (leadng_1.lb1 >= 25 && leadng_1.lb1 <= 28 && dpi_1.lb2 == 0)
				{
					goto L777;
				}
				if (dpi_1.lb2 >= 25 && dpi_1.lb2 <= 28 && leadng_1.lb1 == 0)
				{
					goto L777;
				}

				if ((leadng_1.lb1 == 23 || leadng_1.lb1 == 21) && (dpi_1.lb2 >= 3 && dpi_1.lb2 <= 5))
				{
					goto L889;
				}
				if ((dpi_1.lb2 == 23 || dpi_1.lb2 == 21) && (leadng_1.lb1 >=
																 3 &&
															 leadng_1.lb1 <= 5))
				{
					goto L889;
				}

				if (abs(leadng_1.lb1) == 30 || abs(dpi_1.lb2) == 30)
				{
					goto L400;
				}
				if (leadng_1.lb1 == 21 || dpi_1.lb2 == 21)
				{
					goto L400;
				}
				if (leadng_1.lb1 == 23 || dpi_1.lb2 == 23)
				{
					goto L400;
				}

				if (leadng_1.lb1 >= 3 && leadng_1.lb1 <= 5 && (abs(dpi_1.lb2) == 1 || abs(dpi_1.lb2) == 2 || abs(dpi_1.lb2) >= 6 && abs(dpi_1.lb2) <= 13))
				{
					goto L3;
				}
				if (dpi_1.lb2 >= 3 && dpi_1.lb2 <= 5 && (abs(leadng_1.lb1) == 1 || abs(leadng_1.lb1) == 2 || abs(leadng_1.lb1) >= 6 && abs(leadng_1.lb1) <= 13))
				{
					goto L3;
				}

				if (leadng_1.lb1 >= 25 && leadng_1.lb1 <= 28 && (abs(dpi_1.lb2) == 1 || abs(dpi_1.lb2) == 2 || abs(dpi_1.lb2) >= 6 && abs(dpi_1.lb2) <= 13))
				{
					goto L33;
				}
				if (dpi_1.lb2 >= 25 && dpi_1.lb2 <= 28 && (abs(leadng_1.lb1) == 1 || abs(leadng_1.lb1) == 2 || abs(leadng_1.lb1) >= 6 && abs(leadng_1.lb1) <= 13))
				{
					goto L33;
				}

				if (leadng_1.lb1 == 0 && (abs(dpi_1.lb2) == 1 || abs(dpi_1.lb2) == 2 || abs(dpi_1.lb2) >= 6 && abs(dpi_1.lb2) <= 13))
				{
					goto L547;
				}
				if (dpi_1.lb2 == 0 && (abs(leadng_1.lb1) == 1 || abs(leadng_1.lb1) == 2 || abs(leadng_1.lb1) >= 6 && abs(leadng_1.lb1) <= 13))
				{
					goto L547;
				}

				if ((leadng_1.lb1 == 1 || leadng_1.lb1 == 2) && (dpi_1.lb2 >
																	 5 &&
																 dpi_1.lb2 <= 13))
				{
					goto L44;
				}
				if ((dpi_1.lb2 == 1 || dpi_1.lb2 == 2) && (leadng_1.lb1 > 5 &&
														   leadng_1.lb1 <= 13))
				{
					goto L44;
				}
				if ((leadng_1.lb1 == -1 || leadng_1.lb1 == -2) && (dpi_1.lb2 <
																	   -5 &&
																   dpi_1.lb2 >= -13))
				{
					goto L44;
				}
				if ((dpi_1.lb2 == -1 || dpi_1.lb2 == -2) && (leadng_1.lb1 <
																 -5 &&
															 leadng_1.lb1 >= -13))
				{
					goto L44;
				}

				if ((leadng_1.lb1 == 1 || leadng_1.lb1 == 2) && (dpi_1.lb2 ==
																	 1 ||
																 dpi_1.lb2 == 2))
				{
					goto L4;
				}
				if ((leadng_1.lb1 == -1 || leadng_1.lb1 == -2) && (dpi_1.lb2 == -1 || dpi_1.lb2 == -2))
				{
					goto L4;
				}

				if (leadng_1.lb1 > 5 && leadng_1.lb1 <= 13 && (dpi_1.lb2 > 5 && dpi_1.lb2 <= 13))
				{
					goto L444;
				}
				if (leadng_1.lb1 < -5 && leadng_1.lb1 >= -13 && (dpi_1.lb2 < -5 && dpi_1.lb2 >= -13))
				{
					goto L444;
				}

				if (leadng_1.lb1 < 3 && (dpi_1.lb2 >= 14 && dpi_1.lb2 <= 17))
				{
					goto L400;
				}
				if (dpi_1.lb2 < 3 && (leadng_1.lb1 >= 14 && leadng_1.lb1 <=
																17))
				{
					goto L400;
				}
				if (leadng_1.lb1 >= 14 && leadng_1.lb1 <= 17 && (dpi_1.lb2 >= 14 && dpi_1.lb2 <= 17))
				{
					goto L400;
				}

				goto L400;

			L547:
				if (leadng_1.lb1 * dpi_1.lb2 == 0)
				{
					r__1 = leadng_1.em1 + dpi_1.em2 + .02f;
					ece = r__1 * r__1;
					xkaon0 = 0.f;
					if (srt >= 1.63f && srt <= 1.7f)
					{
						xkaon0 = pnlka_(&srt);
					}
					if (srt > 1.7f)
					{
						xkaon0 = pnlka_(&srt) + pnska_(&srt);
					}
					xkaon0 *= 2.f;
					xkaon = xkaon0;
					xeta = xn1535_(&i1, &i2, &c__0);
					if ((i__4 = ee_1.lb[i1 - 1], abs(i__4)) >= 6 && (i__5 =
																		 ee_1.lb[i1 - 1],
																	 abs(i__5)) <= 13 ||
						(i__6 =
							 ee_1.lb[i2 - 1],
						 abs(i__6)) >= 6 &&
							(i__7 =
								 ee_1.lb[i2 - 1],
							 abs(i__7)) <= 13)
					{
						xeta = 0.f;
					}
					if (xeta + xkaon <= 1e-6f)
					{
						goto L400;
					}
					dse = sqrt((xeta + xkaon) / 3.1415926f);
					deltre = dse + .1f;
					px1cm = pcx;
					py1cm = pcy;
					pz1cm = pcz;
					distce_(&i1, &i2, &deltre, &dse, &input1_1.dt, &ece, &srt,
							&ic, &pcx, &pcy, &pcz);
					if (ic == -1)
					{
						goto L400;
					}
					kkk_1.ekaon[iss * 7 + 3] += 1;
					if (xkaon0 / (xkaon + xeta) > ranart_(&rndf77_1.nseed))
					{
						cren_(&px1cm, &py1cm, &pz1cm, &srt, &i1, &i2, &iblock);
						if (iblock == 7)
						{
							++(*lpn);
						}
						else if (iblock == -7)
						{
						}

						leadng_1.em1 = cc_1.e[i1 - 1];
						dpi_1.em2 = cc_1.e[i2 - 1];
						goto L440;
					}
					resona = 1.f;
					goto L98;
				}
			L3:
				px1cm = pcx;
				py1cm = pcy;
				pz1cm = pcz;
				xkaon0 = 0.f;
				if (srt >= 1.63f && srt <= 1.7f)
				{
					xkaon0 = pnlka_(&srt);
				}
				if (srt > 1.7f)
				{
					xkaon0 = pnlka_(&srt) + pnska_(&srt);
				}
				xkaon0 *= 2.f;

				xphi = 0.f;
				if ((leadng_1.lb1 >= 1 && leadng_1.lb1 <= 2 || leadng_1.lb1 >= 6 && leadng_1.lb1 <= 9 || (dpi_1.lb2 >= 1 && dpi_1.lb2 <= 2 || dpi_1.lb2 >= 6 && dpi_1.lb2 <= 9)) && srt > 1.958f)
				{
					pibphi_(&srt, &leadng_1.lb1, &dpi_1.lb2, &leadng_1.em1, &dpi_1.em2, &xphi, &xphin);
				}
				if ((i__4 = ee_1.lb[i1 - 1], abs(i__4)) >= 6 && (i__5 =
																	 ee_1.lb[i1 - 1],
																 abs(i__5)) <= 13 ||
					(i__6 = ee_1.lb[i2 - 1], abs(i__6)) >= 6 && (i__7 = ee_1.lb[i2 - 1],
																 abs(i__7)) <= 13)
				{
					goto L31;
				}
				r__1 = leadng_1.em1 + dpi_1.em2 + .02f;
				ec = r__1 * r__1;
				xkaon = 0.f;
				if (srt > 1.23f)
				{
					xkaon = (pionpp_(&srt) + pipp1_(&srt)) / 2.f;
				}
				if (leadng_1.lb1 * dpi_1.lb2 == 5 || leadng_1.lb1 * dpi_1.lb2 == 6 && (leadng_1.lb1 == 3 || dpi_1.lb2 == 3) || (leadng_1.lb1 * dpi_1.lb2 == -3 || leadng_1.lb1 * dpi_1.lb2 == -10 && (leadng_1.lb1 == 5 || dpi_1.lb2 == 5)))
				{
					xmax = 190.f;
					xmaxn = 0.f;
					xmaxn1 = 0.f;
					xdirct = dirct1_(&srt);
					goto L678;
				}
				if (leadng_1.lb1 * dpi_1.lb2 == 3 || leadng_1.lb1 * dpi_1.lb2 == 10 && (leadng_1.lb1 == 5 || dpi_1.lb2 == 5) || (leadng_1.lb1 * dpi_1.lb2 == -5 || leadng_1.lb1 * dpi_1.lb2 == -6 && (leadng_1.lb1 == 3 || dpi_1.lb2 == 3)))
				{
					xmax = 27.f;
					xmaxn = 9.9999999999999982f;
					xmaxn1 = 13.333333333333332f;
					xdirct = dirct2_(&srt);
					goto L678;
				}
				if ((leadng_1.lb1 == 4 || dpi_1.lb2 == 4) && ((i__4 =
																   leadng_1.lb1 * dpi_1.lb2,
															   abs(i__4)) == 4 ||
															  (i__5 =
																   leadng_1.lb1 * dpi_1.lb2,
															   abs(i__5)) == 8))
				{
					xmax = 50.f;
					xmaxn = 4.9999999999999991f;
					xmaxn1 = 6.6666666666666661f;
					xdirct = dirct3_(&srt);
					goto L678;
				}
			L678:
				xnpin1 = 0.f;
				xnpin = 0.f;
				xnpid = xnpi_(&i1, &i2, &c__1, &xmax);
				if (xmaxn1 != 0.f)
				{
					xnpin1 = xnpi_(&i1, &i2, &c__2, &xmaxn1);
				}
				if (xmaxn != 0.f)
				{
					xnpin = xnpi_(&i1, &i2, &c__0, &xmaxn);
				}
				xres = xnpid + xnpin + xnpin1;
				xnelas = xres + xdirct;
				icheck = 1;
				goto L34;
			L31:
				r__1 = leadng_1.em1 + dpi_1.em2 + .02f;
				ec = r__1 * r__1;
				xreab = reab_(&i1, &i2, &srt, &c__1);
				if (abs(leadng_1.lb1) >= 10 && abs(leadng_1.lb1) <= 13 || abs(
																			  dpi_1.lb2) >= 10 &&
																			  abs(dpi_1.lb2) <= 13)
				{
					xreab = 0.f;
				}
				xkaon = xkaon0 + xreab;
				if (abs(leadng_1.lb1) > 9 && abs(leadng_1.lb1) <= 13 || abs(
																			dpi_1.lb2) > 9 &&
																			abs(dpi_1.lb2) <= 13)
				{
					xnelas = 1.f;
				}
				else
				{
					xnelas = dpion_(&leadng_1.em1, &dpi_1.em2, &leadng_1.lb1,
									&dpi_1.lb2, &srt);
				}
				icheck = 2;
			L34:
				if (xnelas + xkaon + xphi <= 1e-6f)
				{
					goto L400;
				}
				ds = sqrt((xnelas + xkaon + xphi) / 3.1415926f);
				deltar = ds + .1f;
				distce_(&i1, &i2, &deltar, &ds, &input1_1.dt, &ec, &srt, &ic,
						&pcx, &pcy, &pcz);
				if (ic == -1)
				{
					goto L400;
				}
				kkk_1.ekaon[iss * 7 + 3] += 1;
				if (icheck == 2)
				{
					if (xnelas / (xnelas + xkaon + xphi) >= ranart_(&rndf77_1.nseed))
					{
						crdir_(&px1cm, &py1cm, &pz1cm, &srt, &i1, &i2, &iblock);
						goto L440;
					}
					else
					{
						goto L96;
					}
				}
				if ((xkaon + xphi) / (xkaon + xphi + xnelas) > ranart_(&rndf77_1.nseed))
				{
					goto L95;
				}
				if (xdirct / xnelas >= ranart_(&rndf77_1.nseed))
				{
					crdir_(&px1cm, &py1cm, &pz1cm, &srt, &i1, &i2, &iblock);
					goto L440;
				}
				if (leadng_1.lb1 * dpi_1.lb2 == 5 || leadng_1.lb1 * dpi_1.lb2 == 6 && (leadng_1.lb1 == 3 || dpi_1.lb2 == 3) || (leadng_1.lb1 * dpi_1.lb2 == -3 || leadng_1.lb1 * dpi_1.lb2 == -10 && (leadng_1.lb1 == 5 || dpi_1.lb2 == 5)))
				{

					goto L99;
				}
				else
				{
					xx = (xnpin + xnpin1) / xres;
					if (ranart_(&rndf77_1.nseed) < xx)
					{
						xx0 = xnpin / (xnpin + xnpin1);
						if (ranart_(&rndf77_1.nseed) < xx0)
						{
							resona = 0.f;
							goto L97;
						}
						else
						{
							resona = 1.f;
							goto L98;
						}
					}
					else
					{
						goto L99;
					}
				}
			L97:
				if (resona == 0.f)
				{
					i__ = i1;
					if (leadng_1.em1 < .6f)
					{
						i__ = i2;
					}
					if (leadng_1.lb1 * dpi_1.lb2 == 10 && (leadng_1.lb1 == 5 || dpi_1.lb2 == 5) || leadng_1.lb1 * dpi_1.lb2 ==
																										   -6 &&
																									   (leadng_1.lb1 == 3 || dpi_1.lb2 == 3))
					{
						ee_1.lb[i__ - 1] = 11;
						goto L303;
					}
					if ((i__4 = ee_1.lb[i1 - 1] * ee_1.lb[i2 - 1], abs(i__4)) == 4 && (ee_1.lb[i1 - 1] == 4 || ee_1.lb[i2 - 1] == 4))
					{
						ee_1.lb[i__ - 1] = 11;
						goto L303;
					}
					if ((i__4 = ee_1.lb[i1 - 1] * ee_1.lb[i2 - 1], abs(i__4)) == 8 && (ee_1.lb[i1 - 1] == 4 || ee_1.lb[i2 - 1] == 4))
					{
						ee_1.lb[i__ - 1] = 10;
						goto L303;
					}
					if (ee_1.lb[i1 - 1] * ee_1.lb[i2 - 1] == 3 || ee_1.lb[i1 - 1] * ee_1.lb[i2 - 1] == -5)
					{
						ee_1.lb[i__ - 1] = 10;
					}
				L303:
					dreson_(&i1, &i2);
					if (leadng_1.lb1 < 0 || dpi_1.lb2 < 0)
					{
						ee_1.lb[i__ - 1] = -ee_1.lb[i__ - 1];
					}
					++(*lres);
					goto L101;
				}
			L98:
				if (resona == 1.f)
				{
					i__ = i1;
					if (leadng_1.em1 < .6f)
					{
						i__ = i2;
					}
					if (leadng_1.lb1 * dpi_1.lb2 == 10 && (leadng_1.lb1 == 5 || dpi_1.lb2 == 5) || leadng_1.lb1 * dpi_1.lb2 ==
																										   -6 &&
																									   (leadng_1.lb1 == 3 || dpi_1.lb2 == 3))
					{
						ee_1.lb[i__ - 1] = 13;
						goto L304;
					}
					if ((i__4 = ee_1.lb[i1 - 1] * ee_1.lb[i2 - 1], abs(i__4)) == 4 && (ee_1.lb[i1 - 1] == 4 || ee_1.lb[i2 - 1] == 4))
					{
						ee_1.lb[i__ - 1] = 13;
						goto L304;
					}
					if ((i__4 = ee_1.lb[i1 - 1] * ee_1.lb[i2 - 1], abs(i__4)) == 8 && (ee_1.lb[i1 - 1] == 4 || ee_1.lb[i2 - 1] == 4))
					{
						ee_1.lb[i__ - 1] = 12;
						goto L304;
					}
					if (ee_1.lb[i1 - 1] * ee_1.lb[i2 - 1] == 3 || ee_1.lb[i1 - 1] * ee_1.lb[i2 - 1] == -5)
					{
						ee_1.lb[i__ - 1] = 12;
						goto L304;
					}
					if (ee_1.lb[i1 - 1] * ee_1.lb[i2 - 1] == 0)
					{
						if ((i__4 = ee_1.lb[i1 - 1], abs(i__4)) == 1 || (i__5 = ee_1.lb[i2 - 1], abs(i__5)) == 1)
						{
							ee_1.lb[i__ - 1] = 13;
							goto L304;
						}
						else
						{
							ee_1.lb[i__ - 1] = 12;
						}
					}
				L304:
					dreson_(&i1, &i2);
					if (leadng_1.lb1 < 0 || dpi_1.lb2 < 0)
					{
						ee_1.lb[i__ - 1] = -ee_1.lb[i__ - 1];
					}
					++(*lres);
					goto L101;
				}
			L99:
				++(*lres);
				i__ = i1;
				if (leadng_1.em1 <= .6f)
				{
					i__ = i2;
				}
				if (ee_1.lb[i1 - 1] * ee_1.lb[i2 - 1] == 5 || ee_1.lb[i1 - 1] * ee_1.lb[i2 - 1] == -3)
				{
					ee_1.lb[i__ - 1] = 9;
					goto L305;
				}
				if ((i__4 = ee_1.lb[i1 - 1] * ee_1.lb[i2 - 1], abs(i__4)) ==
						4 &&
					(ee_1.lb[i1 - 1] == 4 || ee_1.lb[i2 - 1] == 4))
				{
					ee_1.lb[i__ - 1] = 8;
					goto L305;
				}
				if (ee_1.lb[i1 - 1] * ee_1.lb[i2 - 1] == 10 && (ee_1.lb[i1 -
																		1] == 5 ||
																ee_1.lb[i2 - 1] == 5) ||
					ee_1.lb[i1 - 1] *
								ee_1.lb[i2 - 1] ==
							-6 &&
						(ee_1.lb[i1 - 1] == 3 ||
						 ee_1.lb[i2 - 1] == 3))
				{
					ee_1.lb[i__ - 1] = 8;
					goto L305;
				}
				if ((i__4 = ee_1.lb[i1 - 1] * ee_1.lb[i2 - 1], abs(i__4)) ==
						8 &&
					(ee_1.lb[i1 - 1] == 4 || ee_1.lb[i2 - 1] == 4))
				{
					ee_1.lb[i__ - 1] = 7;
					goto L305;
				}
				if (ee_1.lb[i1 - 1] * ee_1.lb[i2 - 1] == 3 || ee_1.lb[i1 - 1] * ee_1.lb[i2 - 1] == -5)
				{
					ee_1.lb[i__ - 1] = 7;
					goto L305;
				}
				if (ee_1.lb[i1 - 1] * ee_1.lb[i2 - 1] == 6 && (ee_1.lb[i1 - 1] == 3 || ee_1.lb[i2 - 1] == 3) || ee_1.lb[i1 - 1] *
																															ee_1.lb[i2 - 1] ==
																														-10 &&
																													(ee_1.lb[i1 - 1] == 5 ||
																													 ee_1.lb[i2 - 1] == 5))
				{
					ee_1.lb[i__ - 1] = 6;
				}
			L305:
				dreson_(&i1, &i2);
				if (leadng_1.lb1 < 0 || dpi_1.lb2 < 0)
				{
					ee_1.lb[i__ - 1] = -ee_1.lb[i__ - 1];
				}
				goto L101;
			L889:
				px1cm = pcx;
				py1cm = pcy;
				pz1cm = pcz;
				r__1 = leadng_1.em1 + dpi_1.em2 + .02f;
				ec = r__1 * r__1;
				r__1 = srt - .895f;
				spika = 60.f / (r__1 * r__1 * 4.f / .0025000000000000005f +
								1.f);

				crkpla_(&px1cm, &py1cm, &pz1cm, &ec, &srt, &spika, &emm1, &emm2, &lbp1, &lbp2, &i1, &i2, &icase, &c_b102);
				if (icase == 0)
				{
					iblock = 0;
					goto L400;
				}
				if (icase == 1)
				{
					ksreso_(&i1, &i2);
					iblock = 171;
					++(*lres);
					goto L101;
				}
				else if (icase == 2)
				{
					iblock = 71;
				}
				else if (abs(icase) == 5)
				{
					iblock = 88;
				}
				else
				{

					iblock = 222;
				}
				ee_1.lb[i1 - 1] = lbp1;
				ee_1.lb[i2 - 1] = lbp2;
				cc_1.e[i1 - 1] = emm1;
				cc_1.e[i2 - 1] = emm2;
				leadng_1.em1 = cc_1.e[i1 - 1];
				dpi_1.em2 = cc_1.e[i2 - 1];
				ntag = 0;
				goto L440;

			L33:
				leadng_1.em1 = cc_1.e[i1 - 1];
				dpi_1.em2 = cc_1.e[i2 - 1];
				xelstc = 0.f;
				if (leadng_1.lb1 >= 25 && leadng_1.lb1 <= 28 && (abs(dpi_1.lb2) == 1 || abs(dpi_1.lb2) == 2))
				{
					xelstc = erhon_(&leadng_1.em1, &dpi_1.em2, &leadng_1.lb1,
									&dpi_1.lb2, &srt);
				}
				if (dpi_1.lb2 >= 25 && dpi_1.lb2 <= 28 && (abs(leadng_1.lb1) == 1 || abs(leadng_1.lb1) == 2))
				{
					xelstc = erhon_(&leadng_1.em1, &dpi_1.em2, &leadng_1.lb1,
									&dpi_1.lb2, &srt);
				}
				r__1 = leadng_1.em1 + dpi_1.em2 + .02f;
				ec = r__1 * r__1;
				xkaon0 = 0.f;
				if (srt >= 1.63f && srt <= 1.7f)
				{
					xkaon0 = pnlka_(&srt);
				}
				if (srt > 1.7f)
				{
					xkaon0 = pnlka_(&srt) + pnska_(&srt);
				}
				if (xkaon0 < 0.f)
				{
					xkaon0 = 0.f;
				}
				xkaon0 *= 2.f;
				xkaon = xkaon0;
				ichann = 0;
				xphi = 0.f;
				if (((leadng_1.lb1 >= 1 && leadng_1.lb1 <= 2 || leadng_1.lb1 >= 6 && leadng_1.lb1 <= 9) && (dpi_1.lb2 >= 25 &&
																											dpi_1.lb2 <= 27) ||
					 (dpi_1.lb2 >= 1 && dpi_1.lb2 <= 2 || dpi_1.lb2 >= 6 && dpi_1.lb2 <= 9) && (leadng_1.lb1 >= 25 && leadng_1.lb1 <= 27)) &&
					srt > 1.958f)
				{
					pibphi_(&srt, &leadng_1.lb1, &dpi_1.lb2, &leadng_1.em1, &dpi_1.em2, &xphi, &xphin);
				}
				if (abs(leadng_1.lb1) >= 6 && dpi_1.lb2 >= 25 || leadng_1.lb1 >= 25 && abs(dpi_1.lb2) >= 6)
				{
					ichann = 1;
					ictrl = 2;
					if (leadng_1.lb1 == 28 || dpi_1.lb2 == 28)
					{
						ictrl = 3;
					}
					xreab = reab_(&i1, &i2, &srt, &ictrl);
					if (abs(leadng_1.lb1) >= 10 && abs(leadng_1.lb1) <= 13 ||
						abs(dpi_1.lb2) >= 10 && abs(dpi_1.lb2) <= 13)
					{
						xreab = 0.f;
					}
					if (xreab < 0.f)
					{
						xreab = 1e-6f;
					}
					xkaon = xkaon0 + xreab;
					xelstc = 1.f;
				}
				ds = sqrt((xkaon + xphi + xelstc) / 3.1415926f);

				deltar = ds + .1f;
				px1cm = pcx;
				py1cm = pcy;
				pz1cm = pcz;
				distce_(&i1, &i2, &deltar, &ds, &input1_1.dt, &ec, &srt, &ic,
						&pcx, &pcy, &pcz);
				if (ic == -1)
				{
					goto L400;
				}
				kkk_1.ekaon[iss * 7 + 3] += 1;
				if (xelstc / (xelstc + xkaon + xphi) > ranart_(&rndf77_1.nseed))
				{
					crdir_(&px1cm, &py1cm, &pz1cm, &srt, &i1, &i2, &iblock);
					goto L440;
				}
				crrd_(&px1cm, &py1cm, &pz1cm, &srt, &i1, &i2, &iblock, &xkaon0, &xkaon, &xphi, &xphin);
				if (iblock == 7)
				{
					++(*lpn);
				}
				else if (iblock == -7)
				{
				}
				if (iblock == 81)
				{
					++lrhor;
				}
				if (iblock == 82)
				{
					++lomgar;
				}
				leadng_1.em1 = cc_1.e[i1 - 1];
				dpi_1.em2 = cc_1.e[i2 - 1];
				goto L440;
			L95:
				crpn_(&px1cm, &py1cm, &pz1cm, &srt, &i1, &i2, &iblock, &xkaon0, &xkaon, &xphi, &xphin);
				if (iblock == 7)
				{
					++(*lpn);
				}
				else if (iblock == -7)
				{
				}
				if (iblock == 77)
				{
					++(*lpd);
				}
				if (iblock == 78)
				{
					++(*lrho);
				}
				if (iblock == 79)
				{
					++(*lomega);
				}
				leadng_1.em1 = cc_1.e[i1 - 1];
				dpi_1.em2 = cc_1.e[i2 - 1];
				goto L440;
			L96:
				crpd_(&px1cm, &py1cm, &pz1cm, &srt, &i1, &i2, &iblock, &xkaon0, &xkaon, &xphi, &xphin);
				if (iblock == 7)
				{
					++(*lpn);
				}
				else if (iblock == -7)
				{
				}
				if (iblock == 80)
				{
					++lpdr;
				}
				leadng_1.em1 = cc_1.e[i1 - 1];
				dpi_1.em2 = cc_1.e[i2 - 1];
				goto L440;
			L101:
				if (cc_1.e[i2 - 1] == 0.f)
				{
					goto L600;
				}
				if (cc_1.e[i1 - 1] == 0.f)
				{
					goto L800;
				}
			L44:
				cutoff = leadng_1.em1 + dpi_1.em2 + .02f;
				if (srt <= cutoff)
				{
					goto L400;
				}
				if (srt > 2.245f)
				{
					signn = pp2_(&srt);
				}
				else
				{
					signn = 35.f / ((srt - cutoff) * 100.f + 1.f) + 20.f;
				}
				xnd_(&pcx, &pcy, &pcz, &srt, &i1, &i2, &xinel, &sigk, &xsk1, &xsk2, &xsk3, &xsk4, &xsk5);
				sig = signn + xinel;
				r__1 = leadng_1.em1 + dpi_1.em2 + .02f;
				ec = r__1 * r__1;
				px1cm = pcx;
				py1cm = pcy;
				pz1cm = pcz;
				ianti = 0;
				if (ee_1.lb[i1 - 1] < 0 && ee_1.lb[i2 - 1] < 0)
				{
					ianti = 1;
				}
				sbbdm_(&srt, &sdprod, &ianti, &lbm, &xmm, &pfinal);
				sig += sdprod;
				ipdflag = 0;
				if (para8_1.idpert == 1)
				{
					ipert1 = 1;
					sigr0 = sig;
					dspert = sqrt(sigr0 / 3.1415926f / 10.f);
					dsrpert = dspert + .1f;
					distce_(&i1, &i2, &dsrpert, &dspert, &input1_1.dt, &ec, &srt, &ic, &px1cm, &py1cm, &pz1cm);
					if (ic == -1)
					{
						goto L363;
					}
					signn0 = 0.f;
					crnd_(&irun, &px1cm, &py1cm, &pz1cm, &srt, &i1, &i2, &iblock, &signn0, &sigr0, &sigk, &xsk1, &xsk2, &xsk3, &xsk4, &xsk5, nt, &ipert1);
					ipdflag = 1;
				L363:
					ipert1 = 0;
				}
				if (para8_1.idpert == 2)
				{
					ipert1 = 1;
				}

				ds = sqrt(sig / 31.415925999999999f);
				deltar = ds + .1f;
				distce_(&i1, &i2, &deltar, &ds, &input1_1.dt, &ec, &srt, &ic,
						&px1cm, &py1cm, &pz1cm);
				if (ic == -1)
				{
					if (ipdflag == 1)
					{
						iblock = 501;
					}
					goto L400;
				}
				kkk_1.ekaon[iss * 7 + 2] += 1;
				goto L361;
			L361:
				crnd_(&irun, &px1cm, &py1cm, &pz1cm, &srt, &i1, &i2, &iblock,
					  &signn, &sig, &sigk, &xsk1, &xsk2, &xsk3, &xsk4, &xsk5, nt, &ipert1);
				if (iblock == 0 && ipdflag == 1)
				{
					iblock = 501;
				}
				if (iblock == 11)
				{
					++(*lndk);
					goto L400;
				}
				else if (iblock == -11 || iblock == 501)
				{
					goto L400;
				}
				if (iblock == 222)
				{
					goto L400;
				}
				leadng_1.em1 = cc_1.e[i1 - 1];
				dpi_1.em2 = cc_1.e[i2 - 1];
				goto L440;
			L4:
				cutoff = leadng_1.em1 + dpi_1.em2 + .14f;
				if (srt > 2.245f)
				{
					sig = ppt_(&srt);
					signn = sig - pp1_(&srt);
				}
				else
				{
					sig = xpp_(&srt);
					if (zet[ee_1.lb[i1 - 1] + 45] * zet[ee_1.lb[i2 - 1] + 45] <= 0.f)
					{
						sig = xnp_(&srt);
					}
					if (zet[ee_1.lb[i1 - 1] + 45] * zet[ee_1.lb[i2 - 1] + 45] > 0.f)
					{
						sig = xpp_(&srt);
					}
					if (zet[ee_1.lb[i1 - 1] + 45] == 0.f && zet[ee_1.lb[i2 -
																		1] +
																45] == 0.f)
					{
						sig = xpp_(&srt);
					}
					if (ee_1.lb[i1 - 1] == -1 && ee_1.lb[i2 - 1] == -2 ||
						ee_1.lb[i2 - 1] == -1 && ee_1.lb[i1 - 1] == -2)
					{
						sig = xnp_(&srt);
					}
					if (srt < 1.897f)
					{
						signn = sig;
					}
					else
					{
						signn = 35.f / ((srt - 1.897f) * 100.f + 1.f) + 20.f;
					}
				}
				px1cm = pcx;
				py1cm = pcy;
				pz1cm = pcz;
				ianti = 0;
				if (ee_1.lb[i1 - 1] < 0 && ee_1.lb[i2 - 1] < 0)
				{
					ianti = 1;
				}
				sbbdm_(&srt, &sdprod, &ianti, &lbm, &xmm, &pfinal);
				sig += sdprod;

				ipdflag = 0;
				if (para8_1.idpert == 1)
				{
					ipert1 = 1;
					ec = 4.0481439999999997f;
					sigr0 = sig;
					dspert = sqrt(sigr0 / 3.1415926f / 10.f);
					dsrpert = dspert + .1f;
					distce_(&i1, &i2, &dsrpert, &dspert, &input1_1.dt, &ec, &srt, &ic, &px1cm, &py1cm, &pz1cm);
					if (ic == -1)
					{
						goto L365;
					}
					signn0 = 0.f;
					crnn_(&irun, &px1cm, &py1cm, &pz1cm, &srt, &i1, &i2, &iblock, &ntag, &signn0, &sigr0, nt, &ipert1);
					ipdflag = 1;
				L365:
					ipert1 = 0;
				}
				if (para8_1.idpert == 2)
				{
					ipert1 = 1;
				}

				if (signn <= 0.f)
				{
					if (ipdflag == 1)
					{
						iblock = 501;
					}
					goto L400;
				}

				ec = 3.59709f;
				ds = sqrt(sig / 3.1415926f / 10.f);
				dsr = ds + .1f;
				if (cc_1.e[i1 - 1] >= 1.f && cc_1.e[i2 - 1] >= 1.f)
				{
					ec = 4.75f;
				}
				distce_(&i1, &i2, &dsr, &ds, &input1_1.dt, &ec, &srt, &ic, &px1cm, &py1cm, &pz1cm);
				if (ic == -1)
				{
					if (ipdflag == 1)
					{
						iblock = 501;
					}
					goto L400;
				}

				goto L362;
			L362:
				kkk_1.ekaon[iss * 7] += 1;
				crnn_(&irun, &px1cm, &py1cm, &pz1cm, &srt, &i1, &i2, &iblock,
					  &ntag, &signn, &sig, nt, &ipert1);
				if (iblock == 0 && ipdflag == 1)
				{
					iblock = 501;
				}
				if (iblock == 4 || iblock == 9 || iblock >= 44 || iblock == -9 || iblock == 222 || iblock == 501)
				{

					++(*lcoll);
					if (iblock == 4)
					{
						++(*ldirt);
					}
					else if (iblock == 44)
					{
						++(*lddrho);
					}
					else if (iblock == 45)
					{
						++(*lnnrho);
					}
					else if (iblock == 46)
					{
						++(*lnnom);
					}
					else if (iblock == 222)
					{
					}
					else if (iblock == 9)
					{
						++(*lnnk);
					}
					else if (iblock == -9)
					{
					}
					goto L400;
				}
				leadng_1.em1 = cc_1.e[i1 - 1];
				dpi_1.em2 = cc_1.e[i2 - 1];
				goto L440;
			L505:
				ianti = 0;
				if (ee_1.lb[i1 - 1] < 0 || ee_1.lb[i2 - 1] < 0)
				{
					ianti = 1;
				}
				sdmbb_(&srt, &sdm, &ianti);
				px1cm = pcx;
				py1cm = pcy;
				pz1cm = pcz;
				ec = 4.0481439999999997f;
				ds = sqrt(sdm / 31.4f);
				dsr = ds + .1f;
				distce_(&i1, &i2, &dsr, &ds, &input1_1.dt, &ec, &srt, &ic, &px1cm, &py1cm, &pz1cm);
				if (ic == -1)
				{
					goto L400;
				}
				crdmbb_(&px1cm, &py1cm, &pz1cm, &srt, &i1, &i2, &iblock, &ntag, &sdm, nt, &ianti);
				++(*lcoll);
				goto L400;
			L506:
				ianti = 0;
				if (ee_1.lb[i1 - 1] < 0 || ee_1.lb[i2 - 1] < 0)
				{
					ianti = 1;
				}
				sdbelastic_(&srt, &sdb);
				px1cm = pcx;
				py1cm = pcy;
				pz1cm = pcz;
				ec = 4.0481439999999997f;
				ds = sqrt(sdb / 31.4f);
				dsr = ds + .1f;
				distce_(&i1, &i2, &dsr, &ds, &input1_1.dt, &ec, &srt, &ic, &px1cm, &py1cm, &pz1cm);
				if (ic == -1)
				{
					goto L400;
				}
				crdbel_(&px1cm, &py1cm, &pz1cm, &srt, &i1, &i2, &iblock, &ntag, &sdb, nt, &ianti);
				++(*lcoll);
				goto L400;
			L444:
				cutoff = leadng_1.em1 + dpi_1.em2 + .02f;
				if (srt <= cutoff)
				{
					goto L400;
				}
				if (srt > 2.245f)
				{
					signn = pp2_(&srt);
				}
				else
				{
					signn = 35.f / ((srt - cutoff) * 100.f + 1.f) + 20.f;
				}
				if (signn <= 0.f)
				{
					goto L400;
				}
				xddin_(&pcx, &pcy, &pcz, &srt, &i1, &i2, &xinel, &sigk, &xsk1,
					   &xsk2, &xsk3, &xsk4, &xsk5);
				sig = signn + xinel;
				r__1 = leadng_1.em1 + dpi_1.em2 + .02f;
				ec = r__1 * r__1;
				px1cm = pcx;
				py1cm = pcy;
				pz1cm = pcz;
				ianti = 0;
				if (ee_1.lb[i1 - 1] < 0 && ee_1.lb[i2 - 1] < 0)
				{
					ianti = 1;
				}
				sbbdm_(&srt, &sdprod, &ianti, &lbm, &xmm, &pfinal);
				sig += sdprod;
				ipdflag = 0;
				if (para8_1.idpert == 1)
				{
					ipert1 = 1;
					sigr0 = sig;
					dspert = sqrt(sigr0 / 3.1415926f / 10.f);
					dsrpert = dspert + .1f;
					distce_(&i1, &i2, &dsrpert, &dspert, &input1_1.dt, &ec, &srt, &ic, &px1cm, &py1cm, &pz1cm);
					if (ic == -1)
					{
						goto L367;
					}
					signn0 = 0.f;
					crdd_(&irun, &px1cm, &py1cm, &pz1cm, &srt, &i1, &i2, &iblock, &ntag, &signn0, &sigr0, nt, &ipert1);
					ipdflag = 1;
				L367:
					ipert1 = 0;
				}
				if (para8_1.idpert == 2)
				{
					ipert1 = 1;
				}

				ds = sqrt(sig / 31.4f);
				dsr = ds + .1f;
				distce_(&i1, &i2, &dsr, &ds, &input1_1.dt, &ec, &srt, &ic, &px1cm, &py1cm, &pz1cm);
				if (ic == -1)
				{
					if (ipdflag == 1)
					{
						iblock = 501;
					}
					goto L400;
				}
				goto L364;
			L364:
				kkk_1.ekaon[iss * 7 + 1] += 1;
				crdd_(&irun, &px1cm, &py1cm, &pz1cm, &srt, &i1, &i2, &iblock,
					  &ntag, &signn, &sig, nt, &ipert1);
				if (iblock == 0 && ipdflag == 1)
				{
					iblock = 501;
				}

				if (abs(iblock) == 10)
				{
					++(*lcoll);
					if (iblock == 10)
					{
						++(*lddk);
					}
					else if (iblock == -10)
					{
					}
					goto L400;
				}
				if (iblock == 222 || iblock == 501)
				{
					goto L400;
				}
				leadng_1.em1 = cc_1.e[i1 - 1];
				dpi_1.em2 = cc_1.e[i2 - 1];
				goto L440;
			L777:
				px1cm = pcx;
				py1cm = pcy;
				pz1cm = pcz;
				ec0 = leadng_1.em1 + dpi_1.em2 + .02f;
				if (srt <= ec0)
				{
					goto L400;
				}
				r__1 = leadng_1.em1 + dpi_1.em2 + .02f;
				ec = r__1 * r__1;
				ppel = 20.f;
				ipp = 1;
				if (leadng_1.lb1 < 3 || leadng_1.lb1 > 5 || dpi_1.lb2 < 3 ||
					dpi_1.lb2 > 5)
				{
					goto L778;
				}
				ppxs_(&leadng_1.lb1, &dpi_1.lb2, &srt, &ppsig, &spprho, &ipp);
				ppel = ppsig;
			L778:
				ppink = pipik_(&srt);
				ppink *= 2.f;
				if (leadng_1.lb1 >= 25 && dpi_1.lb2 >= 25)
				{
					ppink = .6f;
				}
				if ((leadng_1.lb1 == 0 || leadng_1.lb1 >= 3 && leadng_1.lb1 <=
																   5) &&
						(dpi_1.lb2 >= 25 && dpi_1.lb2 <= 28) ||
					(dpi_1.lb2 == 0 || dpi_1.lb2 >= 3 && dpi_1.lb2 <= 5) &&
						(leadng_1.lb1 >= 25 && leadng_1.lb1 <= 28))
				{
					ppink = 0.f;
					if (srt >= 1.393f)
					{
						ppink = .3f;
					}
				}
				spprr_(&leadng_1.lb1, &dpi_1.lb2, &srt);
				sppee_(&leadng_1.lb1, &dpi_1.lb2, &srt);
				spppe_(&leadng_1.lb1, &dpi_1.lb2, &srt);
				srpre_(&leadng_1.lb1, &dpi_1.lb2, &srt);
				sopoe_(&leadng_1.lb1, &dpi_1.lb2, &srt);
				srree_(&leadng_1.lb1, &dpi_1.lb2, &srt);
				ppb1_1.ppinnb = 0.f;
				if (srt > ppbmas_1.thresh[0])
				{
					getnst_(&srt);
					if (leadng_1.lb1 >= 3 && leadng_1.lb1 <= 5 && dpi_1.lb2 >= 3 && dpi_1.lb2 <= 5)
					{
						ppb1_1.ppinnb = ppbbar_(&srt);
					}
					else if (leadng_1.lb1 >= 3 && leadng_1.lb1 <= 5 &&
								 dpi_1.lb2 >= 25 && dpi_1.lb2 <= 27 ||
							 dpi_1.lb2 >=
									 3 &&
								 dpi_1.lb2 <= 5 && leadng_1.lb1 >= 25 &&
								 leadng_1.lb1 <= 27)
					{
						ppb1_1.ppinnb = prbbar_(&srt);
					}
					else if (leadng_1.lb1 >= 25 && leadng_1.lb1 <= 27 &&
							 dpi_1.lb2 >= 25 && dpi_1.lb2 <= 27)
					{
						ppb1_1.ppinnb = rrbbar_(&srt);
					}
					else if (leadng_1.lb1 >= 3 && leadng_1.lb1 <= 5 &&
								 dpi_1.lb2 == 28 ||
							 dpi_1.lb2 >= 3 && dpi_1.lb2 <= 5 && leadng_1.lb1 == 28)
					{
						ppb1_1.ppinnb = pobbar_(&srt);
					}
					else if (leadng_1.lb1 >= 25 && leadng_1.lb1 <= 27 &&
								 dpi_1.lb2 == 28 ||
							 dpi_1.lb2 >= 25 && dpi_1.lb2 <= 27 && leadng_1.lb1 == 28)
					{
						ppb1_1.ppinnb = robbar_(&srt);
					}
					else if (leadng_1.lb1 == 28 && dpi_1.lb2 == 28)
					{
						ppb1_1.ppinnb = oobbar_(&srt);
					}
					else
					{
						if (leadng_1.lb1 != 0 && dpi_1.lb2 != 0)
						{
							s_wsle(&io___379);
							do_lio(&c__9, &c__1, "missed MM lb1,lb2=", (ftnlen)18);
							do_lio(&c__3, &c__1, (char *)&leadng_1.lb1, (ftnlen)sizeof(int));
							do_lio(&c__3, &c__1, (char *)&dpi_1.lb2, (ftnlen)sizeof(int));
							e_wsle();
						}
					}
				}
				ppin = ppink + ppb1_1.ppinnb + ppmm_1.pprr + ppmm_1.ppee +
					   ppmm_1.pppe + ppmm_1.rpre + ppmm_1.xopoe +
					   ppmm_1.rree;
				if (ppel + ppin <= .01f)
				{
					goto L400;
				}
				dspp = sqrt((ppel + ppin) / 31.4f);
				dsppr = dspp + .1f;
				distce_(&i1, &i2, &dsppr, &dspp, &input1_1.dt, &ec, &srt, &ic,
						&px1cm, &py1cm, &pz1cm);
				if (ic == -1)
				{
					goto L400;
				}
				if (ppel == 0.f)
				{
					goto L400;
				}
				kkk_1.ekaon[iss * 7 + 4] += 1;
				crpp_(&px1cm, &py1cm, &pz1cm, &srt, &i1, &i2, &iblock, &ppel,
					  &ppin, &spprho, &ipp);
				if (iblock == 666)
				{
					goto L555;
				}
				if (iblock == 6)
				{
					++(*lpp);
				}
				if (iblock == 66)
				{
					++(*lppk);
				}
				else if (iblock == 366)
				{
					++(*lppk);
				}
				else if (iblock == 367)
				{
					++(*lppk);
				}
				leadng_1.em1 = cc_1.e[i1 - 1];
				dpi_1.em2 = cc_1.e[i2 - 1];
				goto L440;
			L2799:
				px1cm = pcx;
				py1cm = pcy;
				pz1cm = pcz;
				r__1 = leadng_1.em1 + dpi_1.em2 + .02f;
				ec = r__1 * r__1;
				dsppb = sqrt(xppbar_(&srt) / 3.1415926f / 10.f);
				dsppbr = dsppb + .1f;
				distce_(&i1, &i2, &dsppbr, &dsppb, &input1_1.dt, &ec, &srt, &ic, &px1cm, &py1cm, &pz1cm);
				if (ic == -1)
				{
					goto L400;
				}
				crppba_(&px1cm, &py1cm, &pz1cm, &srt, &i1, &i2, &iblock);
				leadng_1.em1 = cc_1.e[i1 - 1];
				dpi_1.em2 = cc_1.e[i2 - 1];
				goto L440;

			L3555:
				px1cm = pcx;
				py1cm = pcy;
				pz1cm = pcz;
				r__1 = leadng_1.em1 + dpi_1.em2 + .02f;
				ec = r__1 * r__1;
				dskk = sqrt(sig / 3.1415926f / 10.f);
				dskk0 = dskk + .1f;
				distce_(&i1, &i2, &dskk0, &dskk, &input1_1.dt, &ec, &srt, &ic,
						&px1cm, &py1cm, &pz1cm);
				if (ic == -1)
				{
					goto L400;
				}
				crlaba_(&px1cm, &py1cm, &pz1cm, &srt, &brel, &brsgm, &i1, &i2,
						nt, &iblock, &nchrg, &icase);
				leadng_1.em1 = cc_1.e[i1 - 1];
				dpi_1.em2 = cc_1.e[i2 - 1];
				goto L440;

			L3455:
				px1cm = pcx;
				py1cm = pcy;
				pz1cm = pcz;
				pertur_(&px1cm, &py1cm, &pz1cm, &srt, &irun, &i1, &i2, nt, &kp, &icontp);
				if (icontp == 0)
				{
					leadng_1.em1 = cc_1.e[i1 - 1];
					dpi_1.em2 = cc_1.e[i2 - 1];
					iblock = 727;
					goto L440;
				}
				if (cc_1.e[i1 - 1] == 0.f)
				{
					goto L800;
				}
				if (cc_1.e[i2 - 1] == 0.f)
				{
					goto L600;
				}
				goto L400;

			L7222:
				px1cm = pcx;
				py1cm = pcy;
				pz1cm = pcz;
				r__1 = leadng_1.em1 + dpi_1.em2 + .02f;
				ec = r__1 * r__1;
				xphib_(&leadng_1.lb1, &dpi_1.lb2, &leadng_1.em1, &dpi_1.em2, &srt, &xsk1, &xsk2, &xsk3, &xsk4, &xsk5, &sigp);
				dskk = sqrt(sigp / 3.1415926f / 10.f);
				dskk0 = dskk + .1f;
				distce_(&i1, &i2, &dskk0, &dskk, &input1_1.dt, &ec, &srt, &ic,
						&px1cm, &py1cm, &pz1cm);
				if (ic == -1)
				{
					goto L400;
				}
				crphib_(&px1cm, &py1cm, &pz1cm, &srt, &i1, &i2, &xsk1, &xsk2,
						&xsk3, &xsk4, &xsk5, &sigp, &iblock);
				leadng_1.em1 = cc_1.e[i1 - 1];
				dpi_1.em2 = cc_1.e[i2 - 1];
				goto L440;

			L7444:
				px1cm = pcx;
				py1cm = pcy;
				pz1cm = pcz;
				r__1 = leadng_1.em1 + dpi_1.em2 + .02f;
				ec = r__1 * r__1;
				phimes_(&i1, &i2, &srt, &xsk1, &xsk2, &xsk3, &xsk4, &xsk5, &xsk6, &xsk7, &sigphi);
				dskk = sqrt(sigphi / 3.1415926f / 10.f);
				dskk0 = dskk + .1f;
				distce_(&i1, &i2, &dskk0, &dskk, &input1_1.dt, &ec, &srt, &ic,
						&px1cm, &py1cm, &pz1cm);
				if (ic == -1)
				{
					goto L400;
				}
				pzrt = bb_1.p[i1 * 3 - 1] + bb_1.p[i2 * 3 - 1];
				r__1 = bb_1.p[i1 * 3 - 3];
				r__2 = bb_1.p[i1 * 3 - 2];
				r__3 = bb_1.p[i1 * 3 - 1];
				r__4 = cc_1.e[i1 - 1];
				er1 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
				r__1 = bb_1.p[i2 * 3 - 3];
				r__2 = bb_1.p[i2 * 3 - 2];
				r__3 = bb_1.p[i2 * 3 - 1];
				r__4 = cc_1.e[i2 - 1];
				er2 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
				ert = er1 + er2;
				yy = log((ert + pzrt) / (ert - pzrt)) * .5f;
				crphim_(&px1cm, &py1cm, &pz1cm, &srt, &i1, &i2, &xsk1, &xsk2,
						&xsk3, &xsk4, &xsk5, &xsk6, &sigphi, &ikkg, &ikkl, &iblock);
				leadng_1.em1 = cc_1.e[i1 - 1];
				dpi_1.em2 = cc_1.e[i2 - 1];
				goto L440;

			L7799:
				px1cm = pcx;
				py1cm = pcy;
				pz1cm = pcz;
				r__1 = leadng_1.em1 + dpi_1.em2 + .02f;
				ec = r__1 * r__1;
				lambar_(&i1, &i2, &srt, &siglab);
				dshn = sqrt(siglab / 3.1415926f / 10.f);
				dshnr = dshn + .1f;
				distce_(&i1, &i2, &dshnr, &dshn, &input1_1.dt, &ec, &srt, &ic,
						&px1cm, &py1cm, &pz1cm);
				if (ic == -1)
				{
					goto L400;
				}
				crhb_(&px1cm, &py1cm, &pz1cm, &srt, &i1, &i2, &iblock);
				leadng_1.em1 = cc_1.e[i1 - 1];
				dpi_1.em2 = cc_1.e[i2 - 1];
				goto L440;

			L5699:
				px1cm = pcx;
				py1cm = pcy;
				pz1cm = pcz;
				r__1 = leadng_1.em1 + dpi_1.em2 + .02f;
				ec = r__1 * r__1;
				xkhype_(&i1, &i2, &srt, &xky1, &xky2, &xky3, &xky4, &xky5, &xky6, &xky7, &xky8, &xky9, &xky10, &xky11, &xky12, &xky13, &xky14, &xky15, &xky16, &xky17, &sigk);
				dskk = sqrt(sigk / 3.1415926f);
				dskk0 = dskk + .1f;
				distce_(&i1, &i2, &dskk0, &dskk, &input1_1.dt, &ec, &srt, &ic,
						&px1cm, &py1cm, &pz1cm);
				if (ic == -1)
				{
					goto L400;
				}

				if (ee_1.lb[i1 - 1] == 23 || ee_1.lb[i2 - 1] == 23)
				{
					ikmp = 1;
				}
				else
				{
					ikmp = -1;
				}
				crkhyp_(&px1cm, &py1cm, &pz1cm, &srt, &i1, &i2, &xky1, &xky2,
						&xky3, &xky4, &xky5, &xky6, &xky7, &xky8, &xky9, &xky10, &xky11, &xky12, &xky13, &xky14, &xky15, &xky16,
						&xky17, &sigk, &ikmp, &iblock);
				leadng_1.em1 = cc_1.e[i1 - 1];
				dpi_1.em2 = cc_1.e[i2 - 1];
				goto L440;
			L5999:
				px1cm = pcx;
				py1cm = pcy;
				pz1cm = pcz;
				r__1 = leadng_1.em1 + dpi_1.em2 + .02f;
				ec = r__1 * r__1;
				sigkp = 15.f;
				dskk = sqrt(sigkp / 3.1415926f / 10.f);
				dskk0 = dskk + .1f;
				distce_(&i1, &i2, &dskk0, &dskk, &input1_1.dt, &ec, &srt, &ic,
						&px1cm, &py1cm, &pz1cm);
				if (ic == -1)
				{
					goto L400;
				}

				crlan_(&px1cm, &py1cm, &pz1cm, &srt, &i1, &i2, &iblock);
				leadng_1.em1 = cc_1.e[i1 - 1];
				dpi_1.em2 = cc_1.e[i2 - 1];
				goto L440;

			L8699:
				px1cm = pcx;
				py1cm = pcy;
				pz1cm = pcz;
				r__1 = leadng_1.em1 + dpi_1.em2 + .02f;
				ec = r__1 * r__1;
				crkphi_(&px1cm, &py1cm, &pz1cm, &ec, &srt, &iblock, &emm1, &emm2, &lbp1, &lbp2, &i1, &i2, &ikk, &icase, &c_b118, &c_b119);
				if (icase == 0)
				{
					iblock = 0;
					goto L400;
				}
				if (lbp1 == 29 || lbp2 == 29)
				{
					pzrt = bb_1.p[i1 * 3 - 1] + bb_1.p[i2 * 3 - 1];
					r__1 = bb_1.p[i1 * 3 - 3];
					r__2 = bb_1.p[i1 * 3 - 2];
					r__3 = bb_1.p[i1 * 3 - 1];
					r__4 = cc_1.e[i1 - 1];
					er1 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
					r__1 = bb_1.p[i2 * 3 - 3];
					r__2 = bb_1.p[i2 * 3 - 2];
					r__3 = bb_1.p[i2 * 3 - 1];
					r__4 = cc_1.e[i2 - 1];
					er2 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
					ert = er1 + er2;
					yy = log((ert + pzrt) / (ert - pzrt)) * .5f;
					iblock = 222;
					ntag = 0;
				}
				ee_1.lb[i1 - 1] = lbp1;
				ee_1.lb[i2 - 1] = lbp2;
				cc_1.e[i1 - 1] = emm1;
				cc_1.e[i2 - 1] = emm2;
				leadng_1.em1 = cc_1.e[i1 - 1];
				dpi_1.em2 = cc_1.e[i2 - 1];
				goto L440;
			L8799:
				px1cm = pcx;
				py1cm = pcy;
				pz1cm = pcz;
				r__1 = leadng_1.em1 + dpi_1.em2 + .02f;
				ec = r__1 * r__1;
				crksph_(&px1cm, &py1cm, &pz1cm, &ec, &srt, &emm1, &emm2, &lbp1, &lbp2, &i1, &i2, &ikkg, &ikkl, &iblock, &icase,
						&c_b102);
				if (icase == 0)
				{
					iblock = 0;
					goto L400;
				}

				if (lbp1 == 29 || lbp2 == 20)
				{
					pzrt = bb_1.p[i1 * 3 - 1] + bb_1.p[i2 * 3 - 1];
					r__1 = bb_1.p[i1 * 3 - 3];
					r__2 = bb_1.p[i1 * 3 - 2];
					r__3 = bb_1.p[i1 * 3 - 1];
					r__4 = cc_1.e[i1 - 1];
					er1 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
					r__1 = bb_1.p[i2 * 3 - 3];
					r__2 = bb_1.p[i2 * 3 - 2];
					r__3 = bb_1.p[i2 * 3 - 1];
					r__4 = cc_1.e[i2 - 1];
					er2 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
					ert = er1 + er2;
					yy = log((ert + pzrt) / (ert - pzrt)) * .5f;
				}
				ee_1.lb[i1 - 1] = lbp1;
				ee_1.lb[i2 - 1] = lbp2;
				cc_1.e[i1 - 1] = emm1;
				cc_1.e[i2 - 1] = emm2;
				leadng_1.em1 = cc_1.e[i1 - 1];
				dpi_1.em2 = cc_1.e[i2 - 1];
				goto L440;
			L888:
				px1cm = pcx;
				py1cm = pcy;
				pz1cm = pcz;
				r__1 = leadng_1.em1 + dpi_1.em2 + .02f;
				ec = r__1 * r__1;
				sig = 10.f;
				if (abs(leadng_1.lb1) == 14 || abs(dpi_1.lb2) == 14 || abs(leadng_1.lb1) == 30 || abs(dpi_1.lb2) == 30)
				{
					sig = 20.f;
				}
				if (leadng_1.lb1 == 29 || dpi_1.lb2 == 29)
				{
					sig = 5.f;
				}
				dskn = sqrt(sig / 3.1415926f / 10.f);
				dsknr = dskn + .1f;
				distce_(&i1, &i2, &dsknr, &dskn, &input1_1.dt, &ec, &srt, &ic,
						&px1cm, &py1cm, &pz1cm);
				if (ic == -1)
				{
					goto L400;
				}
				crkn_(&px1cm, &py1cm, &pz1cm, &srt, &i1, &i2, &iblock);
				leadng_1.em1 = cc_1.e[i1 - 1];
				dpi_1.em2 = cc_1.e[i2 - 1];
				goto L440;
			L440:
				if (iblock == 0)
				{
					goto L400;
				}
				++(*lcoll);
				ntag = 0;

				r__1 = leadng_1.em1;
				r__2 = px1cm;
				r__3 = py1cm;
				r__4 = pz1cm;
				e1cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
				p1beta = px1cm * bg_1.betax + py1cm * bg_1.betay + pz1cm * bg_1.betaz;
				transf = bg_1.gamma * (bg_1.gamma * p1beta / (bg_1.gamma + 1) + e1cm);
				pt1i1 = bg_1.betax * transf + px1cm;
				pt2i1 = bg_1.betay * transf + py1cm;
				pt3i1 = bg_1.betaz * transf + pz1cm;
				goto L90002;
			L90002:
				r__1 = dpi_1.em2;
				r__2 = px1cm;
				r__3 = py1cm;
				r__4 = pz1cm;
				e2cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
				transf = bg_1.gamma * (-bg_1.gamma * p1beta / (bg_1.gamma + 1.f) + e2cm);
				pt1i2 = bg_1.betax * transf - px1cm;
				pt2i2 = bg_1.betay * transf - py1cm;
				pt3i2 = bg_1.betaz * transf - pz1cm;
				goto L90003;
			L90003:
				if (iblock == 1)
				{
					++(*lcnne);
				}
				if (iblock == 5)
				{
					++(*ldd);
				}
				if (iblock == 2)
				{
					++(*lcnnd);
				}
				if (iblock == 8)
				{
					++(*lkn);
				}
				if (iblock == 43)
				{
					++(*ldou);
				}
				if (iblock == 3)
				{
					++(*lcndn);
				}
				bb_1.p[i1 * 3 - 3] = pt1i1;
				bb_1.p[i1 * 3 - 2] = pt2i1;
				bb_1.p[i1 * 3 - 1] = pt3i1;
				bb_1.p[i2 * 3 - 3] = pt1i2;
				bb_1.p[i2 * 3 - 2] = pt2i2;
				bb_1.p[i2 * 3 - 1] = pt3i2;
				leadng_1.px1 = bb_1.p[i1 * 3 - 3];
				leadng_1.py1 = bb_1.p[i1 * 3 - 2];
				leadng_1.pz1 = bb_1.p[i1 * 3 - 1];
				leadng_1.em1 = cc_1.e[i1 - 1];
				dpi_1.em2 = cc_1.e[i2 - 1];
				leadng_1.lb1 = ee_1.lb[i1 - 1];
				dpi_1.lb2 = ee_1.lb[i2 - 1];
				ee_1.id[i1 - 1] = 2;
				ee_1.id[i2 - 1] = 2;
				r__1 = leadng_1.em1;
				r__2 = leadng_1.px1;
				r__3 = leadng_1.py1;
				r__4 = leadng_1.pz1;
				leadng_1.e1 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 +
								   r__4 * r__4);
				id1 = ee_1.id[i1 - 1];
				goto L90004;
			L90004:
				am1 = leadng_1.em1;
				am2 = dpi_1.em2;
			L400:

			L555:
			L600:;
			}
		L800:;
		}
		n0 = mass + msum;
		i__2 = rr_1.massr[irun] + msum;
		for (n = n0 + 1; n <= i__2; ++n)
		{
			if (cc_1.e[n - 1] > 0.f || ee_1.lb[n - 1] > 5000)
			{
				++nn_1.nnn;
				pa_1.rpion[(nn_1.nnn + irun * 150001) * 3 - 450006] =
					aa_1.r__[n * 3 - 3];
				pa_1.rpion[(nn_1.nnn + irun * 150001) * 3 - 450005] =
					aa_1.r__[n * 3 - 2];
				pa_1.rpion[(nn_1.nnn + irun * 150001) * 3 - 450004] =
					aa_1.r__[n * 3 - 1];
				if (*nt == *ntmax)
				{
					ftpisv[nn_1.nnn + irun * 150001 - 150002] = ftmax_1.ftsv[n - 1];
					tdecay_1.tfdpi[nn_1.nnn + irun * 150001 - 150002] =
						tdecay_1.tfdcy[n - 1];
				}

				pb_1.ppion[(nn_1.nnn + irun * 150001) * 3 - 450006] = bb_1.p[n * 3 - 3];
				pb_1.ppion[(nn_1.nnn + irun * 150001) * 3 - 450005] = bb_1.p[n * 3 - 2];
				pb_1.ppion[(nn_1.nnn + irun * 150001) * 3 - 450004] = bb_1.p[n * 3 - 1];
				pc_1.epion[nn_1.nnn + irun * 150001 - 150002] = cc_1.e[n - 1];
				pd_1.lpion[nn_1.nnn + irun * 150001 - 150002] = ee_1.lb[n - 1];
				pe_1.propi[nn_1.nnn + irun * 150001 - 150002] = hh_1.proper[n - 1];
				dpert_1.dppion[nn_1.nnn + irun * 150001 - 150002] =
					dpert_1.dpertp[n - 1];
			}
		}
		massrn[irun] = nn_1.nnn + mass;
	}
	ia = 0;
	ib = 0;
	i__1 = run_1.num;
	for (irun = 1; irun <= i__1; ++irun)
	{
		ia += rr_1.massr[irun - 1];
		ib += massrn[irun - 1];
		i__2 = massrn[irun];
		for (ic = 1; ic <= i__2; ++ic)
		{
			ie = ia + ic;
			ig = ib + ic;
			if (ic <= mass)
			{
				rt[ig * 3 - 3] = aa_1.r__[ie * 3 - 3];
				rt[ig * 3 - 2] = aa_1.r__[ie * 3 - 2];
				rt[ig * 3 - 1] = aa_1.r__[ie * 3 - 1];
				if (*nt == *ntmax)
				{
					fttemp[ig - 1] = ftmax_1.ftsv[ie - 1];
					tdecay_1.tft[ig - 1] = tdecay_1.tfdcy[ie - 1];
				}

				pt[ig * 3 - 3] = bb_1.p[ie * 3 - 3];
				pt[ig * 3 - 2] = bb_1.p[ie * 3 - 2];
				pt[ig * 3 - 1] = bb_1.p[ie * 3 - 1];
				et[ig - 1] = cc_1.e[ie - 1];
				lt[ig - 1] = ee_1.lb[ie - 1];
				prot[ig - 1] = hh_1.proper[ie - 1];
				dptemp[ig - 1] = dpert_1.dpertp[ie - 1];
			}
			else
			{
				i0 = ic - mass;
				rt[ig * 3 - 3] = pa_1.rpion[(i0 + irun * 150001) * 3 - 450006];
				rt[ig * 3 - 2] = pa_1.rpion[(i0 + irun * 150001) * 3 - 450005];
				rt[ig * 3 - 1] = pa_1.rpion[(i0 + irun * 150001) * 3 - 450004];
				if (*nt == *ntmax)
				{
					fttemp[ig - 1] = ftpisv[i0 + irun * 150001 - 150002];
					tdecay_1.tft[ig - 1] = tdecay_1.tfdpi[i0 + irun * 150001 - 150002];
				}

				pt[ig * 3 - 3] = pb_1.ppion[(i0 + irun * 150001) * 3 - 450006];
				pt[ig * 3 - 2] = pb_1.ppion[(i0 + irun * 150001) * 3 - 450005];
				pt[ig * 3 - 1] = pb_1.ppion[(i0 + irun * 150001) * 3 - 450004];
				et[ig - 1] = pc_1.epion[i0 + irun * 150001 - 150002];
				lt[ig - 1] = pd_1.lpion[i0 + irun * 150001 - 150002];
				prot[ig - 1] = pe_1.propi[i0 + irun * 150001 - 150002];
				dptemp[ig - 1] = dpert_1.dppion[i0 + irun * 150001 - 150002];
			}
		}
	}

	il = 0;
	i__2 = run_1.num;
	for (irun = 1; irun <= i__2; ++irun)
	{
		rr_1.massr[irun] = massrn[irun];
		il += rr_1.massr[irun - 1];
		i__1 = rr_1.massr[irun];
		for (im = 1; im <= i__1; ++im)
		{
			in = il + im;
			aa_1.r__[in * 3 - 3] = rt[in * 3 - 3];
			aa_1.r__[in * 3 - 2] = rt[in * 3 - 2];
			aa_1.r__[in * 3 - 1] = rt[in * 3 - 1];
			if (*nt == *ntmax)
			{
				ftmax_1.ftsv[in - 1] = fttemp[in - 1];
				tdecay_1.tfdcy[in - 1] = tdecay_1.tft[in - 1];
			}
			bb_1.p[in * 3 - 3] = pt[in * 3 - 3];
			bb_1.p[in * 3 - 2] = pt[in * 3 - 2];
			bb_1.p[in * 3 - 1] = pt[in * 3 - 1];
			cc_1.e[in - 1] = et[in - 1];
			ee_1.lb[in - 1] = lt[in - 1];
			hh_1.proper[in - 1] = prot[in - 1];
			dpert_1.dpertp[in - 1] = dptemp[in - 1];
			if (ee_1.lb[in - 1] < 1 || ee_1.lb[in - 1] > 2)
			{
				ee_1.id[in - 1] = 0;
			}
		}
		hbtout_(&rr_1.massr[irun], nt, ntmax);
	}

	return 0;
}

int cms_(int *i1, int *i2, float *px1cm, float *py1cm,
		 float *pz1cm, float *srt)
{
	float r__1, r__2, r__3, r__4;

	double sqrt(double);

	static float s, e1, e2, em1, em2, px1, py1, pz1, px2, py2, pz2, p1beta,
		etotal, transf;

	px1 = bb_1.p[*i1 * 3 - 3];
	py1 = bb_1.p[*i1 * 3 - 2];
	pz1 = bb_1.p[*i1 * 3 - 1];
	px2 = bb_1.p[*i2 * 3 - 3];
	py2 = bb_1.p[*i2 * 3 - 2];
	pz2 = bb_1.p[*i2 * 3 - 1];
	em1 = cc_1.e[*i1 - 1];
	em2 = cc_1.e[*i2 - 1];
	r__1 = em1;
	r__2 = px1;
	r__3 = py1;
	r__4 = pz1;
	e1 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	r__1 = em2;
	r__2 = px2;
	r__3 = py2;
	r__4 = pz2;
	e2 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	r__1 = e1 + e2;
	r__2 = px1 + px2;
	r__3 = py1 + py2;
	r__4 = pz1 + pz2;
	s = r__1 * r__1 - r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
	*srt = sqrt(s);
	etotal = e1 + e2;
	bg_1.betax = (px1 + px2) / etotal;
	bg_1.betay = (py1 + py2) / etotal;
	bg_1.betaz = (pz1 + pz2) / etotal;
	r__1 = bg_1.betax;
	r__2 = bg_1.betay;
	r__3 = bg_1.betaz;
	bg_1.gamma = 1.f / sqrt(1.f - r__1 * r__1 - r__2 * r__2 - r__3 * r__3);
	p1beta = px1 * bg_1.betax + py1 * bg_1.betay + pz1 * bg_1.betaz;
	transf = bg_1.gamma * (bg_1.gamma * p1beta / (bg_1.gamma + 1) - e1);
	*px1cm = bg_1.betax * transf + px1;
	*py1cm = bg_1.betay * transf + py1;
	*pz1cm = bg_1.betaz * transf + pz1;
	return 0;
}

int distce_(int *i1, int *i2, float *deltar, float *ds,
			float *dt, float *ec, float *srt, int *ic, float *px1cm, float *py1cm,
			float *pz1cm)
{
	float r__1, r__2, r__3, r__4;

	double sqrt(double);

	static float s, e2, x1, y1, z1, x2, y2, z2, px2, py2, pz2, bbb, ddd, dzz,
		drcm, dxcm, dycm, dzcm, prcm, p1beta, drbeta, relvel, rsqare,
		transf;

	*ic = 0;
	x1 = aa_1.r__[*i1 * 3 - 3];
	y1 = aa_1.r__[*i1 * 3 - 2];
	z1 = aa_1.r__[*i1 * 3 - 1];
	leadng_1.px1 = bb_1.p[*i1 * 3 - 3];
	leadng_1.py1 = bb_1.p[*i1 * 3 - 2];
	leadng_1.pz1 = bb_1.p[*i1 * 3 - 1];
	x2 = aa_1.r__[*i2 * 3 - 3];
	y2 = aa_1.r__[*i2 * 3 - 2];
	z2 = aa_1.r__[*i2 * 3 - 1];
	px2 = bb_1.p[*i2 * 3 - 3];
	py2 = bb_1.p[*i2 * 3 - 2];
	pz2 = bb_1.p[*i2 * 3 - 1];
	leadng_1.em1 = cc_1.e[*i1 - 1];
	dpi_1.em2 = cc_1.e[*i2 - 1];
	r__1 = leadng_1.em1;
	r__2 = leadng_1.px1;
	r__3 = leadng_1.py1;
	r__4 = leadng_1.pz1;
	leadng_1.e1 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	r__1 = x1 - x2;
	r__2 = y1 - y2;
	r__3 = z1 - z2;
	rsqare = r__1 * r__1 + r__2 * r__2 + r__3 * r__3;
	r__1 = *deltar;
	if (rsqare > r__1 * r__1)
	{
		goto L400;
	}
	r__1 = dpi_1.em2;
	r__2 = px2;
	r__3 = py2;
	r__4 = pz2;
	e2 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	s = *srt * *srt;
	if (s < *ec)
	{
		goto L400;
	}
	p1beta = leadng_1.px1 * bg_1.betax + leadng_1.py1 * bg_1.betay +
			 leadng_1.pz1 * bg_1.betaz;
	transf = bg_1.gamma * (bg_1.gamma * p1beta / (bg_1.gamma + 1) -
						   leadng_1.e1);
	r__1 = *px1cm;
	r__2 = *py1cm;
	r__3 = *pz1cm;
	prcm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
	if (prcm <= 1e-5f)
	{
		goto L400;
	}
	drbeta = bg_1.betax * (x1 - x2) + bg_1.betay * (y1 - y2) + bg_1.betaz * (z1 - z2);
	transf = bg_1.gamma * bg_1.gamma * drbeta / (bg_1.gamma + 1);
	dxcm = bg_1.betax * transf + x1 - x2;
	dycm = bg_1.betay * transf + y1 - y2;
	dzcm = bg_1.betaz * transf + z1 - z2;
	r__1 = dxcm;
	r__2 = dycm;
	r__3 = dzcm;
	drcm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
	dzz = (*px1cm * dxcm + *py1cm * dycm + *pz1cm * dzcm) / prcm;
	r__1 = drcm;
	r__2 = dzz;
	if (r__1 * r__1 - r__2 * r__2 <= 0.f)
	{
		bbb = 0.f;
	}
	else
	{
		r__1 = drcm;
		r__2 = dzz;
		bbb = sqrt(r__1 * r__1 - r__2 * r__2);
	}
	if (bbb > *ds)
	{
		goto L400;
	}
	relvel = prcm * (1.f / leadng_1.e1 + 1.f / e2);
	ddd = relvel * *dt * .5f;
	if (dabs(ddd) < dabs(dzz))
	{
		goto L400;
	}
	*ic = 1;
	goto L500;
L400:
	*ic = -1;
L500:
	return 0;
}

int crnn_(int *irun, float *px, float *py, float *pz, float *srt, int *i1, int *i2, int *iblock, int *ntag, float *signn, float *sig, int *nt, int *ipert1)
{
	int i__1, i__2, i__3, i__4;
	float r__1, r__2, r__3, r__4, r__5;

	double sqrt(double), exp(double), log(double), atan2(double, double), cos(double), sin(double);

	extern int bbdangle_(float *, float *, float *, int *,
						 int *, int *, int *, float *, float *, int *);
	static float a, x, c1, c2, s1, t1, t2, s2, x1;
	static int ic;
	static float al;
	static int m12, n12;
	static float dm, fm, as, ta, es, pr, ss, cc1, ak0;
	extern double fd5_(float *, float *, float *);
	static float dm1, dm2, dm3, dm4;
	static int id1;
	static float ct1, ct2, pr2, st1, px3, py3, pz3, px4, py4, pz4, pz2, st2,
		ada;
	extern double fde_(float *, float *, float *);
	static float ana;
	static int lbd;
	extern double ang_(float *, int *);
	static int lbm;
	extern int n1535_(int *, int *, float *, float *);
	static float akp, ppd[30000], x1535;
	extern double fns_(float *, float *, float *);
	static float pxd, pyd, xmm, pzd;
	extern double ptr_(float *, int *);
	static float ppx, ppy, ppz, e1cm, e2cm;
	static int lbi1, lbi2;
	extern int ddp2_(float *, int *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *,
					 int *);
	static float dm1n, dm2n, sig1, sig3, sig4, sig2, eti1, eti2;
	extern double ppk0_(float *), ppk1_(float *);
	static float s4pi, pxi1;
	extern double x2pi_(float *), x3pi_(float *), x4pi_(float *);
	static float pyi1, xsk1, xsk2, xsk3, xsk4, xsk5, pzi1, pxi2, pyi2, pzi2;
	static int lbpd[10000];
	static float epcm, dmin__, dmax__, arho, sigk, pt1i1, pt2i1, pt3i1, pt1i2;
	extern double x33pi_(float *);
	static float srho, xdir, pt2i2, pt3i2;
	extern double xrho_(float *);
	static float e1dcm, xptr, t1dlk, t2dlk;
	static int icou1;
	static float t1dsk, t2dsk, t1nlk, t1nsk, t2nsk;
	static int ntry1, ntry2;
	extern int sbbdm_(float *, float *, int *, int *,
					  float *, float *);
	extern double omega_(float *), sigma_(float *, int *, int *,
											 int *);
	extern int ddrho_(float *, int *, float *, float *,
					  float *, float *, float *, float *, float *, float *, float *, float *,
					  float *, float *, int *);
	static int ianti;
	static float pmdlk, signd, amrho, pmdsk;
	extern double pplpk_(float *);
	static float pmnsk;
	extern int pprho_(float *, int *, float *, float *,
					  float *, float *, float *, float *, float *, float *, float *, float *,
					  float *, float *, int *);
	static float xmass, p1beta, p2beta, e2picm, dprob1, pmdlk2, pmdsk2, pmnsk2,
		aomega;
	extern int bbkaon_(int *, float *, float *, float *,
					   float *, float *, float *, float *, float *, float *, float *, float *,
					   float *, int *);
	static float pfinal;
	extern int rmasdd_(float *, float *, float *, float *, float *, int *, int *, float *, float *);
	static float somega, ppbeta, xfinal;
	extern int ppomga_(float *, int *, float *, float *,
					   float *, float *, float *, float *, float *, float *, float *, float *,
					   float *, int *);
	extern double ranart_(int *);
	static float sdprod, xdmass;
	extern int rotate_(float *, float *, float *, float *, float *, float *);
	static float transf;
	static int ndloop, idloop, ipertd;
	static float p1dbeta, p2pibeta;

	n12 = 0;
	m12 = 0;
	*iblock = 0;
	*ntag = 0;
	leadng_1.em1 = cc_1.e[*i1 - 1];
	dpi_1.em2 = cc_1.e[*i2 - 1];
	r__1 = *px;
	r__2 = *py;
	r__3 = *pz;
	pr = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
	c2 = *pz / pr;
	x1 = ranart_(&rndf77_1.nseed);
	ianti = 0;
	if (ee_1.lb[*i1 - 1] < 0 && ee_1.lb[*i2 - 1] < 0)
	{
		ianti = 1;
	}
	sbbdm_(srt, &sdprod, &ianti, &lbm, &xmm, &pfinal);
	if (para8_1.idpert == 1 && *ipert1 == 1)
	{
		if (*srt < 2.012f)
		{
			return 0;
		}
		if (((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 1 || (i__2 = ee_1.lb[*i1 - 1], abs(i__2)) == 2) && ((i__3 = ee_1.lb[*i2 - 1], abs(i__3)) == 1 || (i__4 = ee_1.lb[*i2 - 1], abs(i__4)) == 2))
		{
			goto L108;
		}
		else
		{
			return 0;
		}
	}

	if (x1 <= *signn / *sig)
	{
		r__1 = (*srt - 1.8766f) * 3.65f, r__1 *= r__1;
		as = r__1 * (r__1 * r__1);
		a = as * 6.f / (as + 1.f);
		r__1 = pr;
		ta = r__1 * r__1 * -2.f;
		x = ranart_(&rndf77_1.nseed);
		t1 = (float)log((double)(1.f - x) * exp((double)a * (double)ta) + (double)x) / a;
		c1 = 1.f - t1 / ta;
		t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
		*iblock = 1;
		goto L107;
	}
	else
	{
		if (*srt < 2.012f)
		{
			return 0;
		}
		i__3 = (i__1 = ee_1.lb[*i1 - 1], abs(i__1));
		i__4 = (i__2 = ee_1.lb[*i2 - 1], abs(i__2));
		n1535_(&i__3, &i__4, srt, &x1535);
		sig3 = (x3pi_(srt) + x33pi_(srt)) * 3.f;
		sig4 = x2pi_(srt) * 4.f;
		s4pi = x4pi_(srt);
		srho = xrho_(srt);
		somega = omega_(srt);
		akp = .498f;
		ak0 = .498f;
		ana = .94f;
		ada = 1.232f;
		al = 1.1157f;
		as = 1.1197f;
		xsk1 = 0.f;
		xsk2 = 0.f;
		xsk3 = 0.f;
		xsk4 = 0.f;
		xsk5 = 0.f;
		t1nlk = ana + al + akp;
		if (*srt <= t1nlk)
		{
			goto L222;
		}
		xsk1 = pplpk_(srt) * 1.5f;
		t1dlk = ada + al + akp;
		t2dlk = ada + al - akp;
		if (*srt <= t1dlk)
		{
			goto L222;
		}
		es = *srt;
		r__1 = es;
		r__2 = t1dlk;
		r__3 = es;
		r__4 = t2dlk;
		r__5 = es;
		pmdlk2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / (r__5 * r__5 * 4.f);
		pmdlk = sqrt(pmdlk2);
		xsk3 = pplpk_(srt) * 1.5f;
		t1nsk = ana + as + akp;
		t2nsk = ana + as - akp;
		if (*srt <= t1nsk)
		{
			goto L222;
		}
		r__1 = es;
		r__2 = t1nsk;
		r__3 = es;
		r__4 = t2nsk;
		r__5 = es;
		pmnsk2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / (r__5 * r__5 * 4.f);
		pmnsk = sqrt(pmnsk2);
		xsk2 = (ppk1_(srt) + ppk0_(srt)) * 1.5f;
		t1dsk = ada + as + akp;
		t2dsk = ada + as - akp;
		if (*srt <= t1dsk)
		{
			goto L222;
		}
		r__1 = es;
		r__2 = t1dsk;
		r__3 = es;
		r__4 = t2dsk;
		r__5 = es;
		pmdsk2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / (r__5 * r__5 * 4.f);
		pmdsk = sqrt(pmdsk2);
		xsk4 = (ppk1_(srt) + ppk0_(srt)) * 1.5f;
		if (*srt <= 2.898914f)
		{
			goto L222;
		}
		xsk5 = 1e-4f;
	L222:
		sigk = xsk1 + xsk2 + xsk3 + xsk4;
		xsk1 *= 2.f;
		xsk2 *= 2.f;
		xsk3 *= 2.f;
		xsk4 *= 2.f;
		sigk = sigk * 2.f + xsk5;
		leadng_1.lb1 = (i__1 = ee_1.lb[*i1 - 1], abs(i__1));
		dpi_1.lb2 = (i__1 = ee_1.lb[*i2 - 1], abs(i__1));
		if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 1 || leadng_1.lb1 <= 17 && leadng_1.lb1 >= 14 && (dpi_1.lb2 <= 17 && dpi_1.lb2 >= 14) ||
			leadng_1.lb1 <= 2 && (dpi_1.lb2 <= 17 && dpi_1.lb2 >= 14) ||
			dpi_1.lb2 <= 2 && (leadng_1.lb1 <= 17 && leadng_1.lb1 >= 14))
		{
			if (x1 <= (*signn + sdprod) / *sig)
			{
				goto L108;
			}
			sig1 = sigma_(srt, &c__1, &c__1, &c__0) + sigma_(srt, &c__1, &c__1, &c__1) * .5f;
			sig2 = sigma_(srt, &c__1, &c__1, &c__1) * 1.5f;
			signd = sig1 + sig2 + sig3 + sig4 + x1535 + sigk + s4pi + srho +
					somega;
			if (x1 > (*signn + signd + sdprod) / *sig)
			{
				return 0;
			}
			input_1.dir = sig3 / signd;
			if (ranart_(&rndf77_1.nseed) <= input_1.dir)
			{
				goto L106;
			}
			if (ranart_(&rndf77_1.nseed) <= sigk / (sigk + x1535 + sig4 +
													sig2 + sig1 + s4pi + srho + somega))
			{
				goto L306;
			}
			if (ranart_(&rndf77_1.nseed) <= s4pi / (x1535 + sig4 + sig2 +
													sig1 + s4pi + srho + somega))
			{
				goto L307;
			}
			if (ranart_(&rndf77_1.nseed) <= srho / (x1535 + sig4 + sig2 +
													sig1 + srho + somega))
			{
				goto L308;
			}
			if (ranart_(&rndf77_1.nseed) <= somega / (x1535 + sig4 + sig2 +
													  sig1 + somega))
			{
				goto L309;
			}
			if (ranart_(&rndf77_1.nseed) <= x1535 / (sig1 + sig2 + sig4 +
													 x1535))
			{
				n12 = 9;
			}
			else
			{
				if (ranart_(&rndf77_1.nseed) <= sig4 / (sig1 + sig2 + sig4))
				{
					n12 = 66;
					goto L1012;
				}
				else
				{
					n12 = 3;
					if (ranart_(&rndf77_1.nseed) > sig1 / (sig1 + sig2))
					{
						n12 = 4;
					}
				}
			}
			goto L1011;
		}
		if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2 && (i__2 = ee_1.lb[*i2 - 1], abs(i__2)) == 2)
		{
			if (x1 <= (*signn + sdprod) / *sig)
			{
				goto L108;
			}
			sig1 = sigma_(srt, &c__1, &c__1, &c__0) + sigma_(srt, &c__1, &c__1, &c__1) * .5f;
			sig2 = sigma_(srt, &c__1, &c__1, &c__1) * 1.5f;
			signd = sig1 + sig2 + x1535 + sig3 + sig4 + sigk + s4pi + srho +
					somega;
			if (x1 > (*signn + signd + sdprod) / *sig)
			{
				return 0;
			}
			input_1.dir = sig3 / signd;
			if (ranart_(&rndf77_1.nseed) <= input_1.dir)
			{
				goto L106;
			}
			if (ranart_(&rndf77_1.nseed) <= sigk / (sigk + x1535 + sig4 +
													sig2 + sig1 + s4pi + srho + somega))
			{
				goto L306;
			}
			if (ranart_(&rndf77_1.nseed) <= s4pi / (x1535 + sig4 + sig2 +
													sig1 + s4pi + srho + somega))
			{
				goto L307;
			}
			if (ranart_(&rndf77_1.nseed) <= srho / (x1535 + sig4 + sig2 +
													sig1 + srho + somega))
			{
				goto L308;
			}
			if (ranart_(&rndf77_1.nseed) <= somega / (x1535 + sig4 + sig2 +
													  sig1 + somega))
			{
				goto L309;
			}
			if (ranart_(&rndf77_1.nseed) <= x1535 / (x1535 + sig1 + sig2 +
													 sig4))
			{
				n12 = 10;
			}
			else
			{
				if (ranart_(&rndf77_1.nseed) <= sig4 / (sig1 + sig2 + sig4))
				{
					n12 = 67;
					goto L1013;
				}
				else
				{
					n12 = 6;
					if (ranart_(&rndf77_1.nseed) > sig1 / (sig1 + sig2))
					{
						n12 = 5;
					}
				}
			}
			goto L1011;
		}
		if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 2)
		{
			if (x1 <= (*signn + sdprod) / *sig)
			{
				goto L108;
			}
			sig1 = sigma_(srt, &c__1, &c__1, &c__1) * .5f + sigma_(srt, &c__1,
																   &c__1, &c__0) *
																.25f;
			if (input_1.nstar == 1)
			{
				sig2 = sigma_(srt, &c__2, &c__0, &c__1) * .75f;
			}
			else
			{
				sig2 = 0.f;
			}
			signd = (sig1 + sig2 + x1535) * 2.f + sig3 + sig4 + sigk + s4pi +
					srho + somega;
			if (x1 > (*signn + signd + sdprod) / *sig)
			{
				return 0;
			}
			input_1.dir = sig3 / signd;
			if (ranart_(&rndf77_1.nseed) <= input_1.dir)
			{
				goto L106;
			}
			if (ranart_(&rndf77_1.nseed) <= sigk / (signd - sig3))
			{
				goto L306;
			}
			if (ranart_(&rndf77_1.nseed) <= s4pi / (signd - sig3 - sigk))
			{
				goto L307;
			}
			if (ranart_(&rndf77_1.nseed) <= srho / (signd - sig3 - sigk -
													s4pi))
			{
				goto L308;
			}
			if (ranart_(&rndf77_1.nseed) <= somega / (signd - sig3 - sigk -
													  s4pi - srho))
			{
				goto L309;
			}
			if (ranart_(&rndf77_1.nseed) < x1535 / (sig1 + sig2 + x1535 +
													sig4 * .5f))
			{
				n12 = 11;
				if (ranart_(&rndf77_1.nseed) <= .5f)
				{
					n12 = 12;
				}
			}
			else
			{
				if (ranart_(&rndf77_1.nseed) <= sig4 / (sig4 + (sig1 + sig2) *
																   2.f))
				{
					n12 = 68;
					goto L1014;
				}
				else
				{
					if (ranart_(&rndf77_1.nseed) <= sig1 / (sig1 + sig2))
					{
						n12 = 2;
						if (ranart_(&rndf77_1.nseed) >= .5f)
						{
							n12 = 1;
						}
					}
					else
					{
						n12 = 8;
						if (ranart_(&rndf77_1.nseed) >= .5f)
						{
							n12 = 7;
						}
					}
				}
			}
		}
	L1011:
		*iblock = 2;
		dmax__ = *srt - .9383f - .005f;
		dmax__ = *srt - .9383f - .005f;
		dmin__ = 1.078f;
		if (n12 < 7)
		{
			if (dmax__ < 1.232f)
			{
				fm = fde_(&dmax__, srt, &c_b172);
			}
			else
			{
				xdmass = 1.232f;
				fm = fde_(&xdmass, srt, &c_b173);
			}
			if (fm == 0.f)
			{
				fm = 1e-9f;
			}
			ntry1 = 0;
		L10:
			dm = ranart_(&rndf77_1.nseed) * (dmax__ - dmin__) + dmin__;
			++ntry1;
			if (ranart_(&rndf77_1.nseed) > fde_(&dm, srt, &c_b173) / fm &&
				ntry1 <= 30)
			{
				goto L10;
			}
			if (dm > 1.47f)
			{
				goto L10;
			}
			goto L13;
		}
		if (n12 == 7 || n12 == 8)
		{
			if (dmax__ < 1.44f)
			{
				fm = fns_(&dmax__, srt, &c_b172);
			}
			else
			{
				xdmass = 1.44f;
				fm = fns_(&xdmass, srt, &c_b173);
			}
			if (fm == 0.f)
			{
				fm = 1e-9f;
			}
			ntry2 = 0;
		L11:
			dm = ranart_(&rndf77_1.nseed) * (dmax__ - dmin__) + dmin__;
			++ntry2;
			if (ranart_(&rndf77_1.nseed) > fns_(&dm, srt, &c_b173) / fm &&
				ntry2 <= 10)
			{
				goto L11;
			}
			if (dm > 2.14f)
			{
				goto L11;
			}
			goto L13;
		}
		if (n12 >= 17)
		{
			if (dmax__ < 1.535f)
			{
				fm = fd5_(&dmax__, srt, &c_b172);
			}
			else
			{
				xdmass = 1.535f;
				fm = fd5_(&xdmass, srt, &c_b173);
			}
			if (fm == 0.f)
			{
				fm = 1e-9f;
			}
			ntry1 = 0;
		L12:
			dm = ranart_(&rndf77_1.nseed) * (dmax__ - dmin__) + dmin__;
			++ntry1;
			if (ranart_(&rndf77_1.nseed) > fd5_(&dm, srt, &c_b173) / fm &&
				ntry1 <= 10)
			{
				goto L12;
			}
			if (dm > 1.84f)
			{
				goto L12;
			}
			goto L13;
		}
	L1012:
		*iblock = 43;
		rmasdd_(srt, &c_b185, &c_b185, &c_b187, &c_b187, &input1_1.iseed, &c__1, &dm1, &dm2);
		rmasdd_(srt, &c_b185, &c_b191, &c_b187, &c_b187, &input1_1.iseed, &c__3, &dm1n, &dm2n);
		if (n12 == 66)
		{
			xfinal = ranart_(&rndf77_1.nseed);
			if (xfinal <= .25f)
			{
				ee_1.lb[*i1 - 1] = 9;
				ee_1.lb[*i2 - 1] = 7;
				cc_1.e[*i1 - 1] = dm1;
				cc_1.e[*i2 - 1] = dm2;
				goto L200;
			}
			if (xfinal > .25f && xfinal <= .5f)
			{
				ee_1.lb[*i1 - 1] = 8;
				ee_1.lb[*i2 - 1] = 8;
				cc_1.e[*i1 - 1] = dm1;
				cc_1.e[*i2 - 1] = dm2;
				goto L200;
			}
			if (xfinal > .5f && xfinal <= .75f)
			{
				ee_1.lb[*i1 - 1] = 9;
				ee_1.lb[*i2 - 1] = 10;
				cc_1.e[*i1 - 1] = dm1n;
				cc_1.e[*i2 - 1] = dm2n;
				goto L200;
			}
			if (xfinal > .75f)
			{
				ee_1.lb[*i1 - 1] = 8;
				ee_1.lb[*i2 - 1] = 11;
				cc_1.e[*i1 - 1] = dm1n;
				cc_1.e[*i2 - 1] = dm2n;
				goto L200;
			}
		}
	L1013:
		*iblock = 43;
		rmasdd_(srt, &c_b185, &c_b185, &c_b187, &c_b187, &input1_1.iseed, &c__1, &dm1, &dm2);
		rmasdd_(srt, &c_b185, &c_b191, &c_b187, &c_b187, &input1_1.iseed, &c__3, &dm1n, &dm2n);
		if (n12 == 67)
		{
			xfinal = ranart_(&rndf77_1.nseed);
			if (xfinal <= .25f)
			{
				ee_1.lb[*i1 - 1] = 7;
				ee_1.lb[*i2 - 1] = 7;
				cc_1.e[*i1 - 1] = dm1;
				cc_1.e[*i2 - 1] = dm2;
				goto L200;
			}
			if (xfinal > .25f && xfinal <= .5f)
			{
				ee_1.lb[*i1 - 1] = 6;
				ee_1.lb[*i2 - 1] = 8;
				cc_1.e[*i1 - 1] = dm1;
				cc_1.e[*i2 - 1] = dm2;
				goto L200;
			}
			if (xfinal > .5f && xfinal <= .75f)
			{
				ee_1.lb[*i1 - 1] = 7;
				ee_1.lb[*i2 - 1] = 10;
				cc_1.e[*i1 - 1] = dm1n;
				cc_1.e[*i2 - 1] = dm2n;
				goto L200;
			}
			if (xfinal > .75f)
			{
				ee_1.lb[*i1 - 1] = 8;
				ee_1.lb[*i2 - 1] = 11;
				cc_1.e[*i1 - 1] = dm1n;
				cc_1.e[*i2 - 1] = dm2n;
				goto L200;
			}
		}
	L1014:
		*iblock = 43;
		rmasdd_(srt, &c_b185, &c_b185, &c_b187, &c_b187, &input1_1.iseed, &c__1, &dm1, &dm2);
		rmasdd_(srt, &c_b185, &c_b191, &c_b187, &c_b187, &input1_1.iseed, &c__3, &dm1n, &dm2n);
		if (n12 == 68)
		{
			xfinal = ranart_(&rndf77_1.nseed);
			if (xfinal <= .25f)
			{
				ee_1.lb[*i1 - 1] = 7;
				ee_1.lb[*i2 - 1] = 8;
				cc_1.e[*i1 - 1] = dm1;
				cc_1.e[*i2 - 1] = dm2;
				goto L200;
			}
			if (xfinal > .25f && xfinal <= .5f)
			{
				ee_1.lb[*i1 - 1] = 9;
				ee_1.lb[*i2 - 1] = 6;
				cc_1.e[*i1 - 1] = dm1;
				cc_1.e[*i2 - 1] = dm2;
				goto L200;
			}
			if (xfinal > .5f && xfinal <= .75f)
			{
				ee_1.lb[*i1 - 1] = 7;
				ee_1.lb[*i2 - 1] = 11;
				cc_1.e[*i1 - 1] = dm1n;
				cc_1.e[*i2 - 1] = dm2n;
				goto L200;
			}
			if (xfinal > .75f)
			{
				ee_1.lb[*i1 - 1] = 8;
				ee_1.lb[*i2 - 1] = 10;
				cc_1.e[*i1 - 1] = dm1n;
				cc_1.e[*i2 - 1] = dm2n;
				goto L200;
			}
		}
	L13:
		if (n12 == 1)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 1)
			{
				ee_1.lb[*i2 - 1] = 2;
				ee_1.lb[*i1 - 1] = 8;
				cc_1.e[*i1 - 1] = dm;
			}
			else
			{
				ee_1.lb[*i1 - 1] = 2;
				ee_1.lb[*i2 - 1] = 8;
				cc_1.e[*i2 - 1] = dm;
			}
			goto L200;
		}
		if (n12 == 2)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2)
			{
				ee_1.lb[*i2 - 1] = 1;
				ee_1.lb[*i1 - 1] = 7;
				cc_1.e[*i1 - 1] = dm;
			}
			else
			{
				ee_1.lb[*i1 - 1] = 1;
				ee_1.lb[*i2 - 1] = 7;
				cc_1.e[*i2 - 1] = dm;
			}
			goto L200;
		}
		if (n12 == 3)
		{
			ee_1.lb[*i1 - 1] = 9;
			cc_1.e[*i1 - 1] = dm;
			ee_1.lb[*i2 - 1] = 2;
			cc_1.e[*i2 - 1] = .939457f;
			goto L200;
		}
		if (n12 == 4)
		{
			ee_1.lb[*i2 - 1] = 1;
			ee_1.lb[*i1 - 1] = 8;
			cc_1.e[*i1 - 1] = dm;
			goto L200;
		}
		if (n12 == 5)
		{
			ee_1.lb[*i2 - 1] = 2;
			ee_1.lb[*i1 - 1] = 7;
			cc_1.e[*i1 - 1] = dm;
			goto L200;
		}
		if (n12 == 6)
		{
			ee_1.lb[*i1 - 1] = 6;
			cc_1.e[*i1 - 1] = dm;
			ee_1.lb[*i2 - 1] = 1;
			cc_1.e[*i2 - 1] = .93828f;
			goto L200;
		}
		if (n12 == 7)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 1)
			{
				ee_1.lb[*i1 - 1] = 1;
				ee_1.lb[*i2 - 1] = 10;
				cc_1.e[*i2 - 1] = dm;
			}
			else
			{
				ee_1.lb[*i2 - 1] = 1;
				ee_1.lb[*i1 - 1] = 10;
				cc_1.e[*i1 - 1] = dm;
			}
			goto L200;
		}
		if (n12 == 8)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 1)
			{
				ee_1.lb[*i2 - 1] = 2;
				ee_1.lb[*i1 - 1] = 11;
				cc_1.e[*i1 - 1] = dm;
			}
			else
			{
				ee_1.lb[*i1 - 1] = 2;
				ee_1.lb[*i2 - 1] = 11;
				cc_1.e[*i2 - 1] = dm;
			}
			goto L200;
		}
		if (n12 == 9)
		{
			if (ranart_(&rndf77_1.nseed) <= .5f)
			{
				ee_1.lb[*i2 - 1] = 1;
				ee_1.lb[*i1 - 1] = 13;
				cc_1.e[*i1 - 1] = dm;
			}
			else
			{
				ee_1.lb[*i1 - 1] = 1;
				ee_1.lb[*i2 - 1] = 13;
				cc_1.e[*i2 - 1] = dm;
			}
			goto L200;
		}
		if (n12 == 10)
		{
			if (ranart_(&rndf77_1.nseed) <= .5f)
			{
				ee_1.lb[*i2 - 1] = 2;
				ee_1.lb[*i1 - 1] = 12;
				cc_1.e[*i1 - 1] = dm;
			}
			else
			{
				ee_1.lb[*i1 - 1] = 2;
				ee_1.lb[*i2 - 1] = 12;
				cc_1.e[*i2 - 1] = dm;
			}
			goto L200;
		}
		if (n12 == 11)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2)
			{
				ee_1.lb[*i1 - 1] = 2;
				ee_1.lb[*i2 - 1] = 13;
				cc_1.e[*i2 - 1] = dm;
			}
			else
			{
				ee_1.lb[*i2 - 1] = 2;
				ee_1.lb[*i1 - 1] = 13;
				cc_1.e[*i1 - 1] = dm;
			}
			goto L200;
		}
		if (n12 == 12)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 1)
			{
				ee_1.lb[*i1 - 1] = 1;
				ee_1.lb[*i2 - 1] = 12;
				cc_1.e[*i2 - 1] = dm;
			}
			else
			{
				ee_1.lb[*i2 - 1] = 1;
				ee_1.lb[*i1 - 1] = 12;
				cc_1.e[*i1 - 1] = dm;
			}
		}
	}
L200:
	leadng_1.em1 = cc_1.e[*i1 - 1];
	dpi_1.em2 = cc_1.e[*i2 - 1];
	r__2 = *srt;
	r__3 = leadng_1.em1;
	r__4 = dpi_1.em2;
	r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
	r__5 = leadng_1.em1 * dpi_1.em2;
	pr2 = r__1 * r__1 - r__5 * r__5 * 4.f;
	if (pr2 <= 0.f)
	{
		pr2 = 1e-9f;
	}
	pr = sqrt(pr2) / (*srt * 2.f);
	if (*srt <= 2.14f)
	{
		c1 = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
	}
	if (*srt > 2.14f && *srt <= 2.4f)
	{
		c1 = ang_(srt, &input1_1.iseed);
	}
	if (*srt > 2.4f)
	{
		xptr = pr * .33f;
		cc1 = ptr_(&xptr, &input1_1.iseed);
		r__1 = pr;
		r__2 = cc1;
		c1 = sqrt(r__1 * r__1 - r__2 * r__2) / pr;
	}
	t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
	if (ianti == 1 && ee_1.lb[*i1 - 1] >= 1 && ee_1.lb[*i2 - 1] >= 1)
	{
		ee_1.lb[*i1 - 1] = -ee_1.lb[*i1 - 1];
		ee_1.lb[*i2 - 1] = -ee_1.lb[*i2 - 1];
	}
	goto L107;
L106:
	ntry1 = 0;
L123:
	ddp2_(srt, &input1_1.iseed, &px3, &py3, &pz3, &dm3, &px4, &py4, &pz4, &dm4, &ppx, &ppy, &ppz, &icou1);
	++ntry1;
	if (icou1 < 0 && ntry1 <= 40)
	{
		goto L123;
	}
	rotate_(px, py, pz, &px3, &py3, &pz3);
	rotate_(px, py, pz, &px4, &py4, &pz4);
	rotate_(px, py, pz, &ppx, &ppy, &ppz);
	++nn_1.nnn;
	xdir = ranart_(&rndf77_1.nseed);
	if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 1)
	{
		if (xdir <= .2f)
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 4;
			pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = .13496f;
			ee_1.lb[*i1 - 1] = 9;
			ee_1.lb[*i2 - 1] = 7;
			goto L205;
		}
		if (xdir <= .4f && xdir > .2f)
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 4;
			pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = .13496f;
			ee_1.lb[*i1 - 1] = 8;
			ee_1.lb[*i2 - 1] = 8;
			goto L205;
		}
		if (xdir <= .6f && xdir > .4f)
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 3;
			pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = .13957f;
			ee_1.lb[*i1 - 1] = 9;
			ee_1.lb[*i2 - 1] = 8;
			goto L205;
		}
		if (xdir <= .8f && xdir > .6f)
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 5;
			pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = .13957f;
			ee_1.lb[*i1 - 1] = 9;
			ee_1.lb[*i2 - 1] = 6;
			goto L205;
		}
		if (xdir > .8f)
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 5;
			pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = .13957f;
			ee_1.lb[*i1 - 1] = 7;
			ee_1.lb[*i2 - 1] = 8;
			goto L205;
		}
	}
	if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2 && (i__2 = ee_1.lb[*i2 - 1],
													  abs(i__2)) == 2)
	{
		if (xdir <= .2f)
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 4;
			pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = .13496f;
			ee_1.lb[*i1 - 1] = 6;
			ee_1.lb[*i2 - 1] = 7;
			goto L205;
		}
		if (xdir <= .4f && xdir > .2f)
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 3;
			pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = .13957f;
			ee_1.lb[*i1 - 1] = 6;
			ee_1.lb[*i2 - 1] = 9;
			goto L205;
		}
		if (xdir > .4f && xdir <= .6f)
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 5;
			pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = .13957f;
			ee_1.lb[*i1 - 1] = 9;
			ee_1.lb[*i2 - 1] = 8;
			goto L205;
		}
		if (xdir > .6f && xdir <= .8f)
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 4;
			pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = .13496f;
			ee_1.lb[*i1 - 1] = 7;
			ee_1.lb[*i2 - 1] = 7;
			goto L205;
		}
		if (xdir > .8f)
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 3;
			pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = .13957f;
			ee_1.lb[*i1 - 1] = 7;
			ee_1.lb[*i2 - 1] = 8;
			goto L205;
		}
	}
	if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 2)
	{
		if (xdir <= .17f)
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 4;
			pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = .13496f;
			ee_1.lb[*i1 - 1] = 6;
			ee_1.lb[*i2 - 1] = 9;
			goto L205;
		}
		if (xdir <= .34f && xdir > .17f)
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 3;
			pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = .13957f;
			ee_1.lb[*i1 - 1] = 7;
			ee_1.lb[*i2 - 1] = 9;
			goto L205;
		}
		if (xdir > .34f && xdir <= .51f)
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 5;
			pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = .13957f;
			ee_1.lb[*i1 - 1] = 7;
			ee_1.lb[*i2 - 1] = 8;
			goto L205;
		}
		if (xdir > .51f && xdir <= .68f)
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 3;
			pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = .13957f;
			ee_1.lb[*i1 - 1] = 8;
			ee_1.lb[*i2 - 1] = 8;
			goto L205;
		}
		if (xdir > .68f && xdir <= .85f)
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 4;
			pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = .13957f;
			ee_1.lb[*i1 - 1] = 7;
			ee_1.lb[*i2 - 1] = 8;
			goto L205;
		}
		if (xdir > .85f)
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 5;
			pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = .13957f;
			ee_1.lb[*i1 - 1] = 7;
			ee_1.lb[*i2 - 1] = 7;
		}
	}
L205:
	r__1 = dm3;
	r__2 = px3;
	r__3 = py3;
	r__4 = pz3;
	e1cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	p1beta = px3 * bg_1.betax + py3 * bg_1.betay + pz3 * bg_1.betaz;
	transf = bg_1.gamma * (bg_1.gamma * p1beta / (bg_1.gamma + 1) + e1cm);
	pt1i1 = bg_1.betax * transf + px3;
	pt2i1 = bg_1.betay * transf + py3;
	pt3i1 = bg_1.betaz * transf + pz3;
	eti1 = dm3;

	if (ianti == 1 && ee_1.lb[*i1 - 1] >= 1 && ee_1.lb[*i2 - 1] >= 1)
	{
		ee_1.lb[*i1 - 1] = -ee_1.lb[*i1 - 1];
		ee_1.lb[*i2 - 1] = -ee_1.lb[*i2 - 1];
		if (pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] == 3)
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 5;
		}
		else if (pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] == 5)
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 3;
		}
	}

	leadng_1.lb1 = ee_1.lb[*i1 - 1];
	r__1 = dm4;
	r__2 = px4;
	r__3 = py4;
	r__4 = pz4;
	e2cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	p2beta = px4 * bg_1.betax + py4 * bg_1.betay + pz4 * bg_1.betaz;
	transf = bg_1.gamma * (bg_1.gamma * p2beta / (bg_1.gamma + 1.f) + e2cm);
	pt1i2 = bg_1.betax * transf + px4;
	pt2i2 = bg_1.betay * transf + py4;
	pt3i2 = bg_1.betaz * transf + pz4;
	eti2 = dm4;
	dpi_1.lb2 = ee_1.lb[*i2 - 1];
	bb_1.p[*i1 * 3 - 3] = pt1i1;
	bb_1.p[*i1 * 3 - 2] = pt2i1;
	bb_1.p[*i1 * 3 - 1] = pt3i1;
	cc_1.e[*i1 - 1] = eti1;
	ee_1.lb[*i1 - 1] = leadng_1.lb1;
	bb_1.p[*i2 * 3 - 3] = pt1i2;
	bb_1.p[*i2 * 3 - 2] = pt2i2;
	bb_1.p[*i2 * 3 - 1] = pt3i2;
	cc_1.e[*i2 - 1] = eti2;
	ee_1.lb[*i2 - 1] = dpi_1.lb2;
	leadng_1.px1 = bb_1.p[*i1 * 3 - 3];
	leadng_1.py1 = bb_1.p[*i1 * 3 - 2];
	leadng_1.pz1 = bb_1.p[*i1 * 3 - 1];
	leadng_1.em1 = cc_1.e[*i1 - 1];
	ee_1.id[*i1 - 1] = 2;
	ee_1.id[*i2 - 1] = 2;
	id1 = ee_1.id[*i1 - 1];
	*iblock = 4;
	r__1 = pc_1.epion[nn_1.nnn + *irun * 150001 - 150002];
	r__2 = ppx;
	r__3 = ppy;
	r__4 = ppz;
	epcm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	ppbeta = ppx * bg_1.betax + ppy * bg_1.betay + ppz * bg_1.betaz;
	transf = bg_1.gamma * (bg_1.gamma * ppbeta / (bg_1.gamma + 1.f) + epcm);
	pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = bg_1.betax *
															   transf +
														   ppx;
	pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = bg_1.betay *
															   transf +
														   ppy;
	pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = bg_1.betaz *
															   transf +
														   ppz;
	dpert_1.dppion[nn_1.nnn + *irun * 150001 - 150002] = dpert_1.dpertp[*i1 -
																		1] *
														 dpert_1.dpertp[*i2 - 1];
	pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = aa_1.r__[*i1 * 3 -
																	3];
	pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = aa_1.r__[*i1 * 3 -
																	2];
	pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = aa_1.r__[*i1 * 3 -
																	1];

	goto L90005;
L108:
	if (para8_1.idpert == 1 && *ipert1 == 1 && para8_1.npertd >= 1)
	{
		ndloop = para8_1.npertd;
	}
	else if (para8_1.idpert == 2 && para8_1.npertd >= 1)
	{
		ndloop = para8_1.npertd + 1;
	}
	else
	{
		ndloop = 1;
	}

	dprob1 = sdprod / *sig / (float)para8_1.npertd;
	i__1 = ndloop;
	for (idloop = 1; idloop <= i__1; ++idloop)
	{
		bbdangle_(&pxd, &pyd, &pzd, nt, ipert1, &ianti, &idloop, &pfinal, &dprob1, &lbm);
		rotate_(px, py, pz, &pxd, &pyd, &pzd);
		xmass = 1.8756f;
		r__1 = xmass;
		r__2 = pxd;
		r__3 = pyd;
		r__4 = pzd;
		e1dcm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
		p1dbeta = pxd * bg_1.betax + pyd * bg_1.betay + pzd * bg_1.betaz;
		transf = bg_1.gamma * (bg_1.gamma * p1dbeta / (bg_1.gamma + 1.f) +
							   e1dcm);
		pxi1 = bg_1.betax * transf + pxd;
		pyi1 = bg_1.betay * transf + pyd;
		pzi1 = bg_1.betaz * transf + pzd;
		if (ianti == 0)
		{
			lbd = 42;
		}
		else
		{
			lbd = -42;
		}
		if (para8_1.idpert == 1 && *ipert1 == 1 && para8_1.npertd >= 1)
		{
			++nn_1.nnn;
			pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = pxi1;
			pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = pyi1;
			pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = pzi1;
			pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = 1.8756f;
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = lbd;
			pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = aa_1.r__[*i1 * 3 - 3];
			pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = aa_1.r__[*i1 * 3 - 2];
			pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = aa_1.r__[*i1 * 3 - 1];
			dpert_1.dppion[nn_1.nnn + *irun * 150001 - 150002] = sdprod / *sig / (float)para8_1.npertd;
		}
		else if (para8_1.idpert == 2 && idloop <= para8_1.npertd)
		{
			ppd[idloop * 3 - 3] = pxi1;
			ppd[idloop * 3 - 2] = pyi1;
			ppd[idloop * 3 - 1] = pzi1;
			lbpd[idloop - 1] = lbd;
		}
		else
		{
			cc_1.e[*i1 - 1] = xmm;
			r__1 = xmm;
			r__2 = pxd;
			r__3 = pyd;
			r__4 = pzd;
			e2picm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
			p2pibeta = -pxd * bg_1.betax - pyd * bg_1.betay - pzd * bg_1.betaz;
			transf = bg_1.gamma * (bg_1.gamma * p2pibeta / (bg_1.gamma + 1.f) + e2picm);
			pxi2 = bg_1.betax * transf - pxd;
			pyi2 = bg_1.betay * transf - pyd;
			pzi2 = bg_1.betaz * transf - pzd;
			bb_1.p[*i1 * 3 - 3] = pxi2;
			bb_1.p[*i1 * 3 - 2] = pyi2;
			bb_1.p[*i1 * 3 - 1] = pzi2;
			ee_1.lb[*i1 - 1] = lbm;
			leadng_1.px1 = bb_1.p[*i1 * 3 - 3];
			leadng_1.py1 = bb_1.p[*i1 * 3 - 2];
			leadng_1.pz1 = bb_1.p[*i1 * 3 - 1];
			leadng_1.em1 = cc_1.e[*i1 - 1];
			ee_1.id[*i1 - 1] = 2;
			id1 = ee_1.id[*i1 - 1];
			r__1 = leadng_1.em1;
			r__2 = leadng_1.px1;
			r__3 = leadng_1.py1;
			r__4 = leadng_1.pz1;
			leadng_1.e1 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
			leadng_1.lb1 = ee_1.lb[*i1 - 1];
			bb_1.p[*i2 * 3 - 3] = pxi1;
			bb_1.p[*i2 * 3 - 2] = pyi1;
			bb_1.p[*i2 * 3 - 1] = pzi1;
			ee_1.lb[*i2 - 1] = lbd;
			dpi_1.lb2 = ee_1.lb[*i2 - 1];
			cc_1.e[*i2 - 1] = 1.8756f;
			eti2 = cc_1.e[*i2 - 1];
			ee_1.id[*i2 - 1] = 2;
			if (para8_1.idpert == 2 && idloop == ndloop)
			{
				i__2 = para8_1.npertd;
				for (ipertd = 1; ipertd <= i__2; ++ipertd)
				{
					++nn_1.nnn;
					pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450006] =
						ppd[ipertd * 3 - 3];
					pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450005] =
						ppd[ipertd * 3 - 2];
					pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450004] =
						ppd[ipertd * 3 - 1];
					pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = 1.8756f;
					pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = lbpd[ipertd - 1];
					pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450006] =
						aa_1.r__[*i1 * 3 - 3];
					pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450005] =
						aa_1.r__[*i1 * 3 - 2];
					pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450004] =
						aa_1.r__[*i1 * 3 - 1];
					dpert_1.dppion[nn_1.nnn + *irun * 150001 - 150002] = 1.f /
																		 (float)para8_1.npertd;
				}
			}
		}
	}
	*iblock = 501;
	goto L90005;
L306:
	if (xsk5 / sigk > ranart_(&rndf77_1.nseed))
	{
		leadng_1.pz1 = bb_1.p[*i1 * 3 - 1];
		pz2 = bb_1.p[*i2 * 3 - 1];
		ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 2) + 1;
		ee_1.lb[*i2 - 1] = (int)(ranart_(&rndf77_1.nseed) * 2) + 1;
		++nn_1.nnn;
		pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 29;
		pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = 1.02f;
		*iblock = 222;
		goto L208;
	}

	*iblock = 9;
	if (ianti == 1)
	{
		*iblock = -9;
	}

	leadng_1.pz1 = bb_1.p[*i1 * 3 - 1];
	pz2 = bb_1.p[*i2 * 3 - 1];
	++nn_1.nnn;
	pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 23;
	pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = .498f;
	if (*srt <= 2.63f)
	{
		ic = 1;
		ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 2) + 1;
		ee_1.lb[*i2 - 1] = 14;
		goto L208;
	}
	if (*srt <= 2.74f && *srt > 2.63f)
	{
		if (xsk1 / (xsk1 + xsk2) > ranart_(&rndf77_1.nseed))
		{
			ic = 1;
			ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 2) + 1;
			ee_1.lb[*i2 - 1] = 14;
		}
		else
		{
			ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 2) + 1;
			ee_1.lb[*i2 - 1] = (int)(ranart_(&rndf77_1.nseed) * 3) + 15;
			ic = 2;
		}
		goto L208;
	}
	if (*srt <= 2.77f && *srt > 2.74f)
	{
		if (xsk1 / (xsk1 + xsk2 + xsk3) > ranart_(&rndf77_1.nseed))
		{
			ic = 1;
			ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 2) + 1;
			ee_1.lb[*i2 - 1] = 14;
			goto L208;
		}
		else
		{
			if (xsk2 / (xsk2 + xsk3) > ranart_(&rndf77_1.nseed))
			{
				ic = 2;
				ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 2) +
								   1;
				ee_1.lb[*i2 - 1] = (int)(ranart_(&rndf77_1.nseed) * 3) +
								   15;
			}
			else
			{
				ic = 3;
				ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 4) +
								   6;
				ee_1.lb[*i2 - 1] = 14;
			}
			goto L208;
		}
	}
	if (*srt > 2.77f)
	{
		if (xsk1 / (xsk1 + xsk2 + xsk3 + xsk4) > ranart_(&rndf77_1.nseed))
		{
			ic = 1;
			ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 2) + 1;
			ee_1.lb[*i2 - 1] = 14;
			goto L208;
		}
		else
		{
			if (xsk3 / (xsk2 + xsk3 + xsk4) > ranart_(&rndf77_1.nseed))
			{
				ic = 3;
				ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 4) +
								   6;
				ee_1.lb[*i2 - 1] = 14;
				goto L208;
			}
			else
			{
				if (xsk2 / (xsk2 + xsk4) > ranart_(&rndf77_1.nseed))
				{
					ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) *
												 2) +
									   1;
					ee_1.lb[*i2 - 1] = (int)(ranart_(&rndf77_1.nseed) *
												 3) +
									   15;
					ic = 2;
				}
				else
				{
					ic = 4;
					ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) *
												 4) +
									   6;
					ee_1.lb[*i2 - 1] = (int)(ranart_(&rndf77_1.nseed) *
												 3) +
									   15;
				}
				goto L208;
			}
		}
	}
L208:
	if (ianti == 1 && ee_1.lb[*i1 - 1] >= 1 && ee_1.lb[*i2 - 1] >= 1)
	{
		ee_1.lb[*i1 - 1] = -ee_1.lb[*i1 - 1];
		ee_1.lb[*i2 - 1] = -ee_1.lb[*i2 - 1];
		if (pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] == 23)
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 21;
		}
	}
	ntry1 = 0;
L127:
	bbkaon_(&ic, srt, &px3, &py3, &pz3, &dm3, &px4, &py4, &pz4, &dm4, &ppx, &ppy, &ppz, &icou1);
	++ntry1;
	if (icou1 < 0 && ntry1 <= 20)
	{
		goto L127;
	}
	rotate_(px, py, pz, &px3, &py3, &pz3);
	rotate_(px, py, pz, &px4, &py4, &pz4);
	rotate_(px, py, pz, &ppx, &ppy, &ppz);
	r__1 = dm3;
	r__2 = px3;
	r__3 = py3;
	r__4 = pz3;
	e1cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	p1beta = px3 * bg_1.betax + py3 * bg_1.betay + pz3 * bg_1.betaz;
	transf = bg_1.gamma * (bg_1.gamma * p1beta / (bg_1.gamma + 1) + e1cm);
	pt1i1 = bg_1.betax * transf + px3;
	pt2i1 = bg_1.betay * transf + py3;
	pt3i1 = bg_1.betaz * transf + pz3;
	eti1 = dm3;
	lbi1 = ee_1.lb[*i1 - 1];
	r__1 = dm4;
	r__2 = px4;
	r__3 = py4;
	r__4 = pz4;
	e2cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	p2beta = px4 * bg_1.betax + py4 * bg_1.betay + pz4 * bg_1.betaz;
	transf = bg_1.gamma * (bg_1.gamma * p2beta / (bg_1.gamma + 1.f) + e2cm);
	pt1i2 = bg_1.betax * transf + px4;
	pt2i2 = bg_1.betay * transf + py4;
	pt3i2 = bg_1.betaz * transf + pz4;
	eti2 = dm4;
	lbi2 = ee_1.lb[*i2 - 1];
	r__1 = ppx;
	r__2 = ppy;
	r__3 = ppz;
	epcm = sqrt(r__1 * r__1 + .248004f + r__2 * r__2 + r__3 * r__3);
	ppbeta = ppx * bg_1.betax + ppy * bg_1.betay + ppz * bg_1.betaz;
	transf = bg_1.gamma * (bg_1.gamma * ppbeta / (bg_1.gamma + 1.f) + epcm);
	pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = bg_1.betax *
															   transf +
														   ppx;
	pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = bg_1.betay *
															   transf +
														   ppy;
	pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = bg_1.betaz *
															   transf +
														   ppz;
	dpert_1.dppion[nn_1.nnn + *irun * 150001 - 150002] = dpert_1.dpertp[*i1 -
																		1] *
														 dpert_1.dpertp[*i2 - 1];
	pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = aa_1.r__[*i1 * 3 -
																	3];
	pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = aa_1.r__[*i1 * 3 -
																	2];
	pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = aa_1.r__[*i1 * 3 -
																	1];

	bb_1.p[*i1 * 3 - 3] = pt1i1;
	bb_1.p[*i1 * 3 - 2] = pt2i1;
	bb_1.p[*i1 * 3 - 1] = pt3i1;
	cc_1.e[*i1 - 1] = eti1;
	ee_1.lb[*i1 - 1] = lbi1;
	bb_1.p[*i2 * 3 - 3] = pt1i2;
	bb_1.p[*i2 * 3 - 2] = pt2i2;
	bb_1.p[*i2 * 3 - 1] = pt3i2;
	cc_1.e[*i2 - 1] = eti2;
	ee_1.lb[*i2 - 1] = lbi2;
	leadng_1.px1 = bb_1.p[*i1 * 3 - 3];
	leadng_1.py1 = bb_1.p[*i1 * 3 - 2];
	leadng_1.pz1 = bb_1.p[*i1 * 3 - 1];
	leadng_1.em1 = cc_1.e[*i1 - 1];
	ee_1.id[*i1 - 1] = 2;
	ee_1.id[*i2 - 1] = 2;
	id1 = ee_1.id[*i1 - 1];
	goto L90005;
L307:
	ntry1 = 0;
L125:
	ddrho_(srt, &input1_1.iseed, &px3, &py3, &pz3, &dm3, &px4, &py4, &pz4, &dm4, &ppx, &ppy, &ppz, &amrho, &icou1);
	++ntry1;
	if (icou1 < 0 && ntry1 <= 20)
	{
		goto L125;
	}
	rotate_(px, py, pz, &px3, &py3, &pz3);
	rotate_(px, py, pz, &px4, &py4, &pz4);
	rotate_(px, py, pz, &ppx, &ppy, &ppz);
	++nn_1.nnn;
	arho = amrho;
	xdir = ranart_(&rndf77_1.nseed);
	if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 1)
	{
		if (xdir <= .2f)
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 26;
			pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = arho;
			ee_1.lb[*i1 - 1] = 9;
			ee_1.lb[*i2 - 1] = 7;
			goto L2051;
		}
		if (xdir <= .4f && xdir > .2f)
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 26;
			pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = arho;
			ee_1.lb[*i1 - 1] = 8;
			ee_1.lb[*i2 - 1] = 8;
			goto L2051;
		}
		if (xdir <= .6f && xdir > .4f)
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 25;
			pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = arho;
			ee_1.lb[*i1 - 1] = 9;
			ee_1.lb[*i2 - 1] = 8;
			goto L2051;
		}
		if (xdir <= .8f && xdir > .6f)
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 27;
			pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = arho;
			ee_1.lb[*i1 - 1] = 9;
			ee_1.lb[*i2 - 1] = 6;
			goto L2051;
		}
		if (xdir > .8f)
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 27;
			pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = arho;
			ee_1.lb[*i1 - 1] = 7;
			ee_1.lb[*i2 - 1] = 8;
			goto L2051;
		}
	}
	if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2 && (i__2 = ee_1.lb[*i2 - 1],
													  abs(i__2)) == 2)
	{
		if (xdir <= .2f)
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 26;
			pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = arho;
			ee_1.lb[*i1 - 1] = 6;
			ee_1.lb[*i2 - 1] = 7;
			goto L2051;
		}
		if (xdir <= .4f && xdir > .2f)
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 25;
			pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = arho;
			ee_1.lb[*i1 - 1] = 6;
			ee_1.lb[*i2 - 1] = 9;
			goto L2051;
		}
		if (xdir > .4f && xdir <= .6f)
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 27;
			pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = arho;
			ee_1.lb[*i1 - 1] = 9;
			ee_1.lb[*i2 - 1] = 8;
			goto L2051;
		}
		if (xdir > .6f && xdir <= .8f)
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 26;
			pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = arho;
			ee_1.lb[*i1 - 1] = 7;
			ee_1.lb[*i2 - 1] = 7;
			goto L2051;
		}
		if (xdir > .8f)
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 25;
			pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = arho;
			ee_1.lb[*i1 - 1] = 7;
			ee_1.lb[*i2 - 1] = 8;
			goto L2051;
		}
	}
	if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 2)
	{
		if (xdir <= .17f)
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 25;
			pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = arho;
			ee_1.lb[*i1 - 1] = 6;
			ee_1.lb[*i2 - 1] = 9;
			goto L2051;
		}
		if (xdir <= .34f && xdir > .17f)
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 25;
			pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = arho;
			ee_1.lb[*i1 - 1] = 7;
			ee_1.lb[*i2 - 1] = 9;
			goto L2051;
		}
		if (xdir > .34f && xdir <= .51f)
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 27;
			pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = arho;
			ee_1.lb[*i1 - 1] = 7;
			ee_1.lb[*i2 - 1] = 8;
			goto L2051;
		}
		if (xdir > .51f && xdir <= .68f)
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 25;
			pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = arho;
			ee_1.lb[*i1 - 1] = 8;
			ee_1.lb[*i2 - 1] = 8;
			goto L2051;
		}
		if (xdir > .68f && xdir <= .85f)
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 26;
			pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = arho;
			ee_1.lb[*i1 - 1] = 7;
			ee_1.lb[*i2 - 1] = 8;
			goto L2051;
		}
		if (xdir > .85f)
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 27;
			pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = arho;
			ee_1.lb[*i1 - 1] = 7;
			ee_1.lb[*i2 - 1] = 7;
		}
	}
L2051:
	r__1 = dm3;
	r__2 = px3;
	r__3 = py3;
	r__4 = pz3;
	e1cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	p1beta = px3 * bg_1.betax + py3 * bg_1.betay + pz3 * bg_1.betaz;
	transf = bg_1.gamma * (bg_1.gamma * p1beta / (bg_1.gamma + 1) + e1cm);
	pt1i1 = bg_1.betax * transf + px3;
	pt2i1 = bg_1.betay * transf + py3;
	pt3i1 = bg_1.betaz * transf + pz3;
	eti1 = dm3;

	if (ianti == 1 && ee_1.lb[*i1 - 1] >= 1 && ee_1.lb[*i2 - 1] >= 1)
	{
		ee_1.lb[*i1 - 1] = -ee_1.lb[*i1 - 1];
		ee_1.lb[*i2 - 1] = -ee_1.lb[*i2 - 1];
		if (pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] == 25)
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 27;
		}
		else if (pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] == 27)
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 25;
		}
	}

	leadng_1.lb1 = ee_1.lb[*i1 - 1];
	r__1 = dm4;
	r__2 = px4;
	r__3 = py4;
	r__4 = pz4;
	e2cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	p2beta = px4 * bg_1.betax + py4 * bg_1.betay + pz4 * bg_1.betaz;
	transf = bg_1.gamma * (bg_1.gamma * p2beta / (bg_1.gamma + 1.f) + e2cm);
	pt1i2 = bg_1.betax * transf + px4;
	pt2i2 = bg_1.betay * transf + py4;
	pt3i2 = bg_1.betaz * transf + pz4;
	eti2 = dm4;
	dpi_1.lb2 = ee_1.lb[*i2 - 1];
	bb_1.p[*i1 * 3 - 3] = pt1i1;
	bb_1.p[*i1 * 3 - 2] = pt2i1;
	bb_1.p[*i1 * 3 - 1] = pt3i1;
	cc_1.e[*i1 - 1] = eti1;
	ee_1.lb[*i1 - 1] = leadng_1.lb1;
	bb_1.p[*i2 * 3 - 3] = pt1i2;
	bb_1.p[*i2 * 3 - 2] = pt2i2;
	bb_1.p[*i2 * 3 - 1] = pt3i2;
	cc_1.e[*i2 - 1] = eti2;
	ee_1.lb[*i2 - 1] = dpi_1.lb2;
	leadng_1.px1 = bb_1.p[*i1 * 3 - 3];
	leadng_1.py1 = bb_1.p[*i1 * 3 - 2];
	leadng_1.pz1 = bb_1.p[*i1 * 3 - 1];
	leadng_1.em1 = cc_1.e[*i1 - 1];
	ee_1.id[*i1 - 1] = 2;
	ee_1.id[*i2 - 1] = 2;
	id1 = ee_1.id[*i1 - 1];
	*iblock = 44;
	r__1 = pc_1.epion[nn_1.nnn + *irun * 150001 - 150002];
	r__2 = ppx;
	r__3 = ppy;
	r__4 = ppz;
	epcm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	ppbeta = ppx * bg_1.betax + ppy * bg_1.betay + ppz * bg_1.betaz;
	transf = bg_1.gamma * (bg_1.gamma * ppbeta / (bg_1.gamma + 1.f) + epcm);
	pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = bg_1.betax *
															   transf +
														   ppx;
	pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = bg_1.betay *
															   transf +
														   ppy;
	pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = bg_1.betaz *
															   transf +
														   ppz;
	dpert_1.dppion[nn_1.nnn + *irun * 150001 - 150002] = dpert_1.dpertp[*i1 -
																		1] *
														 dpert_1.dpertp[*i2 - 1];
	pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = aa_1.r__[*i1 * 3 -
																	3];
	pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = aa_1.r__[*i1 * 3 -
																	2];
	pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = aa_1.r__[*i1 * 3 -
																	1];

	goto L90005;
L308:
	ntry1 = 0;
L126:
	pprho_(srt, &input1_1.iseed, &px3, &py3, &pz3, &dm3, &px4, &py4, &pz4, &dm4, &ppx, &ppy, &ppz, &amrho, &icou1);
	++ntry1;
	if (icou1 < 0 && ntry1 <= 20)
	{
		goto L126;
	}
	rotate_(px, py, pz, &px3, &py3, &pz3);
	rotate_(px, py, pz, &px4, &py4, &pz4);
	rotate_(px, py, pz, &ppx, &ppy, &ppz);
	++nn_1.nnn;
	arho = amrho;
	xdir = ranart_(&rndf77_1.nseed);
	if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 1)
	{
		if (xdir <= .5f)
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 26;
			pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = arho;
			ee_1.lb[*i1 - 1] = 1;
			ee_1.lb[*i2 - 1] = 1;
			goto L2052;
		}
		else
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 27;
			pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = arho;
			ee_1.lb[*i1 - 1] = 1;
			ee_1.lb[*i2 - 1] = 2;
			goto L2052;
		}
	}
	if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2 && (i__2 = ee_1.lb[*i2 - 1],
													  abs(i__2)) == 2)
	{
		if (xdir <= .5f)
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 26;
			pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = arho;
			ee_1.lb[*i1 - 1] = 2;
			ee_1.lb[*i2 - 1] = 2;
			goto L2052;
		}
		else
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 25;
			pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = arho;
			ee_1.lb[*i1 - 1] = 1;
			ee_1.lb[*i2 - 1] = 2;
			goto L2052;
		}
	}
	if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 2)
	{
		if (xdir <= .33f)
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 26;
			pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = arho;
			ee_1.lb[*i1 - 1] = 1;
			ee_1.lb[*i2 - 1] = 2;
			goto L2052;
		}
		else if (xdir <= .67f && xdir > .34f)
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 25;
			pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = arho;
			ee_1.lb[*i1 - 1] = 1;
			ee_1.lb[*i2 - 1] = 1;
			goto L2052;
		}
		else
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 27;
			pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = arho;
			ee_1.lb[*i1 - 1] = 2;
			ee_1.lb[*i2 - 1] = 2;
			goto L2052;
		}
	}
L2052:
	r__1 = dm3;
	r__2 = px3;
	r__3 = py3;
	r__4 = pz3;
	e1cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	p1beta = px3 * bg_1.betax + py3 * bg_1.betay + pz3 * bg_1.betaz;
	transf = bg_1.gamma * (bg_1.gamma * p1beta / (bg_1.gamma + 1) + e1cm);
	pt1i1 = bg_1.betax * transf + px3;
	pt2i1 = bg_1.betay * transf + py3;
	pt3i1 = bg_1.betaz * transf + pz3;
	eti1 = dm3;

	if (ianti == 1 && ee_1.lb[*i1 - 1] >= 1 && ee_1.lb[*i2 - 1] >= 1)
	{
		ee_1.lb[*i1 - 1] = -ee_1.lb[*i1 - 1];
		ee_1.lb[*i2 - 1] = -ee_1.lb[*i2 - 1];
		if (pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] == 25)
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 27;
		}
		else if (pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] == 27)
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 25;
		}
	}

	leadng_1.lb1 = ee_1.lb[*i1 - 1];
	r__1 = dm4;
	r__2 = px4;
	r__3 = py4;
	r__4 = pz4;
	e2cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	p2beta = px4 * bg_1.betax + py4 * bg_1.betay + pz4 * bg_1.betaz;
	transf = bg_1.gamma * (bg_1.gamma * p2beta / (bg_1.gamma + 1.f) + e2cm);
	pt1i2 = bg_1.betax * transf + px4;
	pt2i2 = bg_1.betay * transf + py4;
	pt3i2 = bg_1.betaz * transf + pz4;
	eti2 = dm4;
	dpi_1.lb2 = ee_1.lb[*i2 - 1];
	bb_1.p[*i1 * 3 - 3] = pt1i1;
	bb_1.p[*i1 * 3 - 2] = pt2i1;
	bb_1.p[*i1 * 3 - 1] = pt3i1;
	cc_1.e[*i1 - 1] = eti1;
	ee_1.lb[*i1 - 1] = leadng_1.lb1;
	bb_1.p[*i2 * 3 - 3] = pt1i2;
	bb_1.p[*i2 * 3 - 2] = pt2i2;
	bb_1.p[*i2 * 3 - 1] = pt3i2;
	cc_1.e[*i2 - 1] = eti2;
	ee_1.lb[*i2 - 1] = dpi_1.lb2;
	leadng_1.px1 = bb_1.p[*i1 * 3 - 3];
	leadng_1.py1 = bb_1.p[*i1 * 3 - 2];
	leadng_1.pz1 = bb_1.p[*i1 * 3 - 1];
	leadng_1.em1 = cc_1.e[*i1 - 1];
	ee_1.id[*i1 - 1] = 2;
	ee_1.id[*i2 - 1] = 2;
	id1 = ee_1.id[*i1 - 1];
	*iblock = 45;
	r__1 = pc_1.epion[nn_1.nnn + *irun * 150001 - 150002];
	r__2 = ppx;
	r__3 = ppy;
	r__4 = ppz;
	epcm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	ppbeta = ppx * bg_1.betax + ppy * bg_1.betay + ppz * bg_1.betaz;
	transf = bg_1.gamma * (bg_1.gamma * ppbeta / (bg_1.gamma + 1.f) + epcm);
	pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = bg_1.betax *
															   transf +
														   ppx;
	pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = bg_1.betay *
															   transf +
														   ppy;
	pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = bg_1.betaz *
															   transf +
														   ppz;
	dpert_1.dppion[nn_1.nnn + *irun * 150001 - 150002] = dpert_1.dpertp[*i1 -
																		1] *
														 dpert_1.dpertp[*i2 - 1];
	pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = aa_1.r__[*i1 * 3 -
																	3];
	pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = aa_1.r__[*i1 * 3 -
																	2];
	pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = aa_1.r__[*i1 * 3 -
																	1];

	goto L90005;
L309:
	ntry1 = 0;
L138:
	ppomga_(srt, &input1_1.iseed, &px3, &py3, &pz3, &dm3, &px4, &py4, &pz4, &dm4, &ppx, &ppy, &ppz, &icou1);
	++ntry1;
	if (icou1 < 0 && ntry1 <= 20)
	{
		goto L138;
	}
	rotate_(px, py, pz, &px3, &py3, &pz3);
	rotate_(px, py, pz, &px4, &py4, &pz4);
	rotate_(px, py, pz, &ppx, &ppy, &ppz);
	++nn_1.nnn;
	aomega = .782f;
	if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 1)
	{
		pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 28;
		pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = aomega;
		ee_1.lb[*i1 - 1] = 1;
		ee_1.lb[*i2 - 1] = 1;
		goto L2053;
	}
	if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2 && (i__2 = ee_1.lb[*i2 - 1],
													  abs(i__2)) == 2)
	{
		pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 28;
		pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = aomega;
		ee_1.lb[*i1 - 1] = 2;
		ee_1.lb[*i2 - 1] = 2;
		goto L2053;
	}
	if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 2)
	{
		pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 28;
		pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = aomega;
		ee_1.lb[*i1 - 1] = 1;
		ee_1.lb[*i2 - 1] = 2;
		goto L2053;
	}
L2053:
	r__1 = dm3;
	r__2 = px3;
	r__3 = py3;
	r__4 = pz3;
	e1cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	p1beta = px3 * bg_1.betax + py3 * bg_1.betay + pz3 * bg_1.betaz;
	transf = bg_1.gamma * (bg_1.gamma * p1beta / (bg_1.gamma + 1) + e1cm);
	pt1i1 = bg_1.betax * transf + px3;
	pt2i1 = bg_1.betay * transf + py3;
	pt3i1 = bg_1.betaz * transf + pz3;
	eti1 = dm3;
	if (ianti == 1 && ee_1.lb[*i1 - 1] >= 1 && ee_1.lb[*i2 - 1] >= 1)
	{
		ee_1.lb[*i1 - 1] = -ee_1.lb[*i1 - 1];
		ee_1.lb[*i2 - 1] = -ee_1.lb[*i2 - 1];
	}
	leadng_1.lb1 = ee_1.lb[*i1 - 1];
	r__1 = dm4;
	r__2 = px4;
	r__3 = py4;
	r__4 = pz4;
	e2cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	p2beta = px4 * bg_1.betax + py4 * bg_1.betay + pz4 * bg_1.betaz;
	transf = bg_1.gamma * (bg_1.gamma * p2beta / (bg_1.gamma + 1.f) + e2cm);
	pt1i2 = bg_1.betax * transf + px4;
	pt2i2 = bg_1.betay * transf + py4;
	pt3i2 = bg_1.betaz * transf + pz4;
	eti2 = dm4;
	dpi_1.lb2 = ee_1.lb[*i2 - 1];
	bb_1.p[*i1 * 3 - 3] = pt1i1;
	bb_1.p[*i1 * 3 - 2] = pt2i1;
	bb_1.p[*i1 * 3 - 1] = pt3i1;
	cc_1.e[*i1 - 1] = eti1;
	ee_1.lb[*i1 - 1] = leadng_1.lb1;
	bb_1.p[*i2 * 3 - 3] = pt1i2;
	bb_1.p[*i2 * 3 - 2] = pt2i2;
	bb_1.p[*i2 * 3 - 1] = pt3i2;
	cc_1.e[*i2 - 1] = eti2;
	ee_1.lb[*i2 - 1] = dpi_1.lb2;
	leadng_1.px1 = bb_1.p[*i1 * 3 - 3];
	leadng_1.py1 = bb_1.p[*i1 * 3 - 2];
	leadng_1.pz1 = bb_1.p[*i1 * 3 - 1];
	leadng_1.em1 = cc_1.e[*i1 - 1];
	ee_1.id[*i1 - 1] = 2;
	ee_1.id[*i2 - 1] = 2;
	id1 = ee_1.id[*i1 - 1];
	*iblock = 46;
	r__1 = pc_1.epion[nn_1.nnn + *irun * 150001 - 150002];
	r__2 = ppx;
	r__3 = ppy;
	r__4 = ppz;
	epcm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	ppbeta = ppx * bg_1.betax + ppy * bg_1.betay + ppz * bg_1.betaz;
	transf = bg_1.gamma * (bg_1.gamma * ppbeta / (bg_1.gamma + 1.f) + epcm);
	pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = bg_1.betax *
															   transf +
														   ppx;
	pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = bg_1.betay *
															   transf +
														   ppy;
	pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = bg_1.betaz *
															   transf +
														   ppz;
	dpert_1.dppion[nn_1.nnn + *irun * 150001 - 150002] = dpert_1.dpertp[*i1 -
																		1] *
														 dpert_1.dpertp[*i2 - 1];
	pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = aa_1.r__[*i1 * 3 -
																	3];
	pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = aa_1.r__[*i1 * 3 -
																	2];
	pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = aa_1.r__[*i1 * 3 -
																	1];

	goto L90005;
L90005:
	return 0;
L107:
	if (*px == 0.f && *py == 0.f)
	{
		t2 = 0.f;
	}
	else
	{
		t2 = atan2(*py, *px);
	}
	r__1 = c1;
	s1 = 1.f - r__1 * r__1;
	if (s1 <= 0.f)
	{
		s1 = 0.f;
	}
	s1 = sqrt(s1);
	r__1 = c2;
	s2 = sqrt(1.f - r__1 * r__1);
	ct1 = cos(t1);
	st1 = sin(t1);
	ct2 = cos(t2);
	st2 = sin(t2);
	*pz = pr * (c1 * c2 - s1 * s2 * ct1);
	ss = c2 * s1 * ct1 + s2 * c1;
	*px = pr * (ss * ct2 - s1 * st1 * st2);
	*py = pr * (ss * st2 + s1 * st1 * ct2);
	return 0;
}

int crpp_(float *px, float *py, float *pz, float *srt, int *i1, int *i2, int *iblock, float *ppel, float *ppin, float *spprho, int *ipp)
{
	float r__1, r__2, r__3, r__4, r__5;

	double sqrt(double), cos(double), sin(double);

	static float c1, s1, t1, pr, ei1, ei2;
	static int lb1, lb2;
	static float em1, em2, ct1, pr2, px0, py0, pz0, st1;
	static int lbb1, lbb2, lb1i, lb2i, ntag;
	extern int pi2et2_(int *, int *, int *,
					   int *, float *, float *, int *, int *);
	static float ranpi;
	extern int pi2ro2_(int *, int *, int *,
					   int *, float *, float *, int *, int *),
		ro2et2_(int
					*,
				int *, int *, int *, float *, float *, int *,
				int *),
		pi3eta_(int *, int *, int *, int *,
				float *, float *, int *, int *),
		bbarfs_(int *, int *, float *, float *, int *, int *);
	extern double ranart_(int *);
	extern int rotate_(float *, float *, float *, float *, float *, float *), opioet_(int *, int *, int *, int *, float *, float *, int *, int *), rhores_(int *, int *), rpiret_(int *, int *, int *, int *, float *, float *, int *, int *);

	lb1i = ee_1.lb[*i1 - 1];
	lb2i = ee_1.lb[*i2 - 1];
	px0 = *px;
	py0 = *py;
	pz0 = *pz;
	*iblock = 1;
	if (*srt > .996f && *ppin / (*ppin + *ppel) > ranart_(&rndf77_1.nseed))
	{
		ranpi = ranart_(&rndf77_1.nseed);
		if (ppmm_1.pprr / *ppin >= ranpi)
		{
			pi2ro2_(i1, i2, &lbb1, &lbb2, &ei1, &ei2, iblock, &input1_1.iseed);
		}
		else if ((ppmm_1.pprr + ppmm_1.ppee) / *ppin >= ranpi)
		{
			pi2et2_(i1, i2, &lbb1, &lbb2, &ei1, &ei2, iblock, &input1_1.iseed);
		}
		else if ((ppmm_1.pprr + ppmm_1.ppee + ppmm_1.pppe) / *ppin >= ranpi)
		{
			pi3eta_(i1, i2, &lbb1, &lbb2, &ei1, &ei2, iblock, &input1_1.iseed);
		}
		else if ((ppmm_1.pprr + ppmm_1.ppee + ppmm_1.pppe + ppmm_1.rpre) / *ppin >= ranpi)
		{
			rpiret_(i1, i2, &lbb1, &lbb2, &ei1, &ei2, iblock, &input1_1.iseed);
		}
		else if ((ppmm_1.pprr + ppmm_1.ppee + ppmm_1.pppe + ppmm_1.rpre +
				  ppmm_1.xopoe) /
					 *ppin >=
				 ranpi)
		{
			opioet_(i1, i2, &lbb1, &lbb2, &ei1, &ei2, iblock, &input1_1.iseed);
		}
		else if ((ppmm_1.pprr + ppmm_1.ppee + ppmm_1.pppe + ppmm_1.rpre +
				  ppmm_1.xopoe + ppmm_1.rree) /
					 *ppin >=
				 ranpi)
		{
			ro2et2_(i1, i2, &lbb1, &lbb2, &ei1, &ei2, iblock, &input1_1.iseed);
		}
		else if ((ppmm_1.pprr + ppmm_1.ppee + ppmm_1.pppe + ppmm_1.rpre +
				  ppmm_1.xopoe + ppmm_1.rree + ppb1_1.ppinnb) /
					 *ppin >=
				 ranpi)
		{
			bbarfs_(&lbb1, &lbb2, &ei1, &ei2, iblock, &input1_1.iseed);
		}
		else
		{
			*iblock = 66;
			ei1 = .498f;
			ei2 = .498f;
			lbb1 = 21;
			lbb2 = 23;
			lb1 = ee_1.lb[*i1 - 1];
			lb2 = ee_1.lb[*i2 - 1];
			if ((lb1 == 0 || lb1 >= 3 && lb1 <= 5) && (lb2 >= 25 && lb2 <= 28) || (lb2 == 0 || lb2 >= 3 && lb2 <= 5) && (lb1 >= 25 &&
																														 lb1 <= 28))
			{
				ei1 = .895f;
				ei2 = .498f;
				if (ranart_(&rndf77_1.nseed) >= .5f)
				{
					*iblock = 366;
					lbb1 = 30;
					lbb2 = 21;
				}
				else
				{
					*iblock = 367;
					lbb1 = -30;
					lbb2 = 23;
				}
			}
		}
		cc_1.e[*i1 - 1] = ei1;
		cc_1.e[*i2 - 1] = ei2;
		ee_1.lb[*i1 - 1] = lbb1;
		ee_1.lb[*i2 - 1] = lbb2;
	}
	else
	{
		if ((ee_1.lb[*i1 - 1] < 3 || ee_1.lb[*i1 - 1] > 5) && (ee_1.lb[*i2 -
																	   1] < 3 ||
															   ee_1.lb[*i2 - 1] > 5))
		{
			return 0;
		}
		*iblock = 6;
		if (*ipp == 1 || *ipp == 4 || *ipp == 6)
		{
			goto L10;
		}
		if (*spprho / *ppel > ranart_(&rndf77_1.nseed))
		{
			goto L20;
		}
	}
L10:
	ntag = 0;
	em1 = cc_1.e[*i1 - 1];
	em2 = cc_1.e[*i2 - 1];
	r__2 = *srt;
	r__3 = em1;
	r__4 = em2;
	r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
	r__5 = em1 * em2;
	pr2 = r__1 * r__1 - r__5 * r__5 * 4.f;
	if (pr2 <= 0.f)
	{
		pr2 = 1e-9f;
	}
	pr = sqrt(pr2) / (*srt * 2.f);
	c1 = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
	t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
	r__1 = c1;
	s1 = sqrt(1.f - r__1 * r__1);
	ct1 = cos(t1);
	st1 = sin(t1);
	*pz = pr * c1;
	*px = pr * s1 * ct1;
	*py = pr * s1 * st1;
	rotate_(&px0, &py0, &pz0, px, py, pz);
	return 0;
L20:
	*iblock = 666;
	rhores_(i1, i2);
	if (*ipp == 2)
	{
		ee_1.lb[*i1 - 1] = 27;
	}
	if (*ipp == 3)
	{
		ee_1.lb[*i1 - 1] = 26;
	}
	if (*ipp == 5)
	{
		ee_1.lb[*i1 - 1] = 25;
	}
	return 0;
}

int crnd_(int *irun, float *px, float *py, float *pz, float *srt, int *i1, int *i2, int *iblock, float *signn, float *sig, float *sigk, float *xsk1, float *xsk2, float *xsk3, float *xsk4, float *xsk5, int *nt, int *ipert1)
{
	int i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8;
	float r__1, r__2, r__3, r__4;

	double sqrt(double), exp(double), log(double), atan2(double, double), cos(double), sin(double);

	extern int bbdangle_(float *, float *, float *, int *,
						 int *, int *, int *, float *, float *, int *);
	static float a, x, c1, c2, s1, t1, t2, s2, x1;
	static int ic, m12, n12;
	static float dm, fm, as, ta, pr, ss, cc1;
	extern double fd5_(float *, float *, float *);
	static float dm3, dm4, pf2, ct1, ct2;
	static int id1;
	static float am1, am2, st1, st2, pz2, px3, py3, pz3, px4, py4, pz4;
	static int lbd;
	extern double ang_(float *, int *);
	static int lbm;
	extern int m1535_(int *, int *, float *, float *);
	static float x1440, ppd[30000], x1535;
	extern double fns_(float *, float *, float *);
	static float prf, xmm, pxd, pyd;
	extern double ptr_(float *, int *);
	static float ppx, ppy, ppz, pzd, e1cm, e2cm;
	static int lbi1, lbi2;
	static float eti1, eti2, pxi1, pyi1, pzi1, pxi2, pyi2, pzi2;
	static int lbpd[10000];
	static float epcm, dmin__;
	static int ntag;
	static float dmax__, pt1i1, pt2i1, pt3i1, pt1i2, pt2i2, pt3i2, e1dcm, xptr;
	static int icou1, ntry1, ntry2;
	extern int sbbdm_(float *, float *, int *, int *,
					  float *, float *);
	extern double sigma_(float *, int *, int *, int *), denom_(
																		   float *, float *);
	static int ianti;
	static float signd, sigdn, renom, xmass, p1beta, p2beta, e2picm, dprob1,
		renom1;
	extern int bbkaon_(int *, float *, float *, float *,
					   float *, float *, float *, float *, float *, float *, float *, float *,
					   float *, int *);
	static float deltam, pfinal, ppbeta;
	extern double ranart_(int *);
	static float sdprod, renomn, xdmass;
	extern int rotate_(float *, float *, float *, float *, float *, float *);
	static float transf;
	static int ndloop, idloop, ipertd;
	static float p1dbeta, p2pibeta;

	n12 = 0;
	m12 = 0;
	*iblock = 0;
	ntag = 0;
	leadng_1.em1 = cc_1.e[*i1 - 1];
	dpi_1.em2 = cc_1.e[*i2 - 1];
	r__1 = *px;
	r__2 = *py;
	r__3 = *pz;
	pr = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
	c2 = *pz / pr;
	x1 = ranart_(&rndf77_1.nseed);
	ianti = 0;
	if (ee_1.lb[*i1 - 1] < 0 && ee_1.lb[*i2 - 1] < 0)
	{
		ianti = 1;
	}
	sbbdm_(srt, &sdprod, &ianti, &lbm, &xmm, &pfinal);
	if (para8_1.idpert == 1 && *ipert1 == 1)
	{
		if (*srt < 2.012f)
		{
			return 0;
		}
		if (((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 1 || (i__2 = ee_1.lb[*i1 - 1], abs(i__2)) == 2) && ((i__3 = ee_1.lb[*i2 - 1], abs(i__3)) >= 6 && (i__4 = ee_1.lb[*i2 - 1], abs(i__4)) <= 13))
		{
			goto L108;
		}
		else if (((i__1 = ee_1.lb[*i2 - 1], abs(i__1)) == 1 || (i__2 =
																	ee_1.lb[*i2 - 1],
																abs(i__2)) == 2) &&
				 ((i__3 = ee_1.lb[*i1 -
								  1],
				   abs(i__3)) >= 6 &&
				  (i__4 = ee_1.lb[*i1 - 1], abs(i__4)) <=
					  13))
		{
			goto L108;
		}
		else
		{
			return 0;
		}
	}
	if (x1 <= *signn / *sig)
	{
		r__1 = (*srt - 1.8766f) * 3.65f, r__1 *= r__1;
		as = r__1 * (r__1 * r__1);
		a = as * 6.f / (as + 1.f);
		r__1 = pr;
		ta = r__1 * r__1 * -2.f;
		x = ranart_(&rndf77_1.nseed);
		t1 = (float)log((double)(1.f - x) * exp((double)a * (double)ta) + (double)x) / a;
		c1 = 1.f - t1 / ta;
		t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
		*iblock = 1;
		goto L107;
	}
	else
	{
		if (*srt < 2.04f)
		{
			return 0;
		}
		if (((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2 || (i__2 = ee_1.lb[*i2 - 1], abs(i__2)) == 2) && ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 20 || ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 13)
		{
			if (x1 <= (*signn + sdprod) / *sig)
			{
				goto L108;
			}
		}

		r__1 = *srt;
		prf = sqrt(r__1 * r__1 * .25f - .88040689000000005f);
		if (leadng_1.em1 > 1.f)
		{
			deltam = leadng_1.em1;
		}
		else
		{
			deltam = dpi_1.em2;
		}
		r__1 = prf;
		renom = deltam * (r__1 * r__1) / denom_(srt, &c_b173) / pr;
		r__1 = prf;
		renomn = deltam * (r__1 * r__1) / denom_(srt, &c_b234) / pr;
		r__1 = prf;
		renom1 = deltam * (r__1 * r__1) / denom_(srt, &c_b235) / pr;
		if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2 && (i__2 = ee_1.lb[*i2 - 1], abs(i__2)) == 6)
		{
			renom = 0.f;
		}
		if ((i__1 = ee_1.lb[*i2 - 1], abs(i__1)) == 2 && (i__2 = ee_1.lb[*i1 - 1], abs(i__2)) == 6)
		{
			renom = 0.f;
		}
		if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 1 && (i__2 = ee_1.lb[*i2 - 1], abs(i__2)) == 9)
		{
			renom = 0.f;
		}
		if ((i__1 = ee_1.lb[*i2 - 1], abs(i__1)) == 1 && (i__2 = ee_1.lb[*i1 - 1], abs(i__2)) == 9)
		{
			renom = 0.f;
		}
		i__3 = (i__1 = ee_1.lb[*i1 - 1], abs(i__1));
		i__4 = (i__2 = ee_1.lb[*i2 - 1], abs(i__2));
		m1535_(&i__3, &i__4, srt, &x1535);
		x1440 = sigma_(srt, &c__2, &c__0, &c__1) * .75f;
		if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2 && (i__2 = ee_1.lb[*i2 - 1], abs(i__2)) == 6 || (i__3 = ee_1.lb[*i2 - 1], abs(i__3)) == 2 && (i__4 = ee_1.lb[*i1 - 1], abs(i__4)) == 6 || (i__5 = ee_1.lb[*i1 - 1], abs(i__5)) == 1 && (i__6 = ee_1.lb[*i2 - 1], abs(i__6)) == 9 || (i__7 = ee_1.lb[*i2 - 1], abs(i__7)) == 1 && (i__8 = ee_1.lb[*i1 - 1], abs(i__8)) == 9)
		{
			if (x1 <= (*signn + sdprod) / *sig)
			{
				goto L108;
			}
			if ((*sigk + *signn + sdprod) / *sig >= x1)
			{
				goto L306;
			}
		}
		if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 18 && ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2 || (i__2 = ee_1.lb[*i2 - 1], abs(i__2)) == 2))
		{
			signd = sigma_(srt, &c__1, &c__1, &c__0) + sigma_(srt, &c__1, &c__1, &c__1) * .5f;
			sigdn = signd * .25f * renom;
			if (x1 <= (*signn + sdprod) / *sig)
			{
				goto L108;
			}
			if (x1 > (*signn + sigdn + x1440 + x1535 + *sigk + sdprod) / *sig)
			{
				return 0;
			}

			if (*sigk / (*sigk + sigdn + x1440 + x1535) > ranart_(&rndf77_1.nseed))
			{
				goto L306;
			}
			if (ranart_(&rndf77_1.nseed) < sigdn / (sigdn + x1440 + x1535))
			{
				m12 = 3;
				goto L206;
			}
			else
			{
				if (ranart_(&rndf77_1.nseed) < x1440 / (x1440 + x1535))
				{
					m12 = 37;
				}
				else
				{
					return 0;
				}
				goto L204;
			}
		}
		if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 6 && ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 1 || (i__2 = ee_1.lb[*i2 - 1], abs(i__2)) == 1))
		{
			signd = sigma_(srt, &c__1, &c__1, &c__0) + sigma_(srt, &c__1, &c__1, &c__1) * .5f;
			sigdn = signd * .25f * renom;
			if (x1 <= (*signn + sdprod) / *sig)
			{
				goto L108;
			}
			if (x1 > (*signn + sigdn + x1440 + x1535 + *sigk + sdprod) / *sig)
			{
				return 0;
			}

			if (*sigk / (*sigk + sigdn + x1440 + x1535) > ranart_(&rndf77_1.nseed))
			{
				goto L306;
			}
			if (ranart_(&rndf77_1.nseed) < sigdn / (sigdn + x1440 + x1535))
			{
				m12 = 6;
				goto L206;
			}
			else
			{
				if (ranart_(&rndf77_1.nseed) < x1440 / (x1440 + x1535))
				{
					m12 = 47;
				}
				else
				{
					return 0;
				}
				goto L204;
			}
		}
		if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 8 && ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 1 || (i__2 = ee_1.lb[*i2 - 1], abs(i__2)) == 1))
		{
			signd = sigma_(srt, &c__1, &c__1, &c__1) * 1.5f;
			sigdn = signd * .25f * renom;
			if (x1 <= (*signn + sdprod) / *sig)
			{
				goto L108;
			}
			if (x1 > (*signn + sigdn + x1440 + x1535 + *sigk + sdprod) / *sig)
			{
				return 0;
			}

			if (*sigk / (*sigk + sigdn + x1440 + x1535) > ranart_(&rndf77_1.nseed))
			{
				goto L306;
			}
			if (ranart_(&rndf77_1.nseed) < sigdn / (sigdn + x1440 + x1535))
			{
				m12 = 4;
				goto L206;
			}
			else
			{
				if (ranart_(&rndf77_1.nseed) < x1440 / (x1440 + x1535))
				{
					m12 = 39;
				}
				else
				{
					m12 = 40;
				}
				goto L204;
			}
		}
		if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 14 && ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2 || (i__2 = ee_1.lb[*i2 - 1], abs(i__2)) == 2))
		{
			signd = sigma_(srt, &c__1, &c__1, &c__1) * 1.5f;
			sigdn = signd * .25f * renom;
			if (x1 <= (*signn + sdprod) / *sig)
			{
				goto L108;
			}
			if (x1 > (*signn + sigdn + x1440 + x1535 + *sigk + sdprod) / *sig)
			{
				return 0;
			}

			if (*sigk / (*sigk + sigdn + x1440 + x1535) > ranart_(&rndf77_1.nseed))
			{
				goto L306;
			}
			if (ranart_(&rndf77_1.nseed) < sigdn / (sigdn + x1440 + x1535))
			{
				m12 = 5;
				goto L206;
			}
			else
			{
				if (ranart_(&rndf77_1.nseed) < x1440 / (x1440 + x1535))
				{
					m12 = 48;
				}
				else
				{
					m12 = 49;
				}
				goto L204;
			}
		}
		if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 16 && ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2 || (i__2 = ee_1.lb[*i2 - 1], abs(i__2)) == 2))
		{
			signd = sigma_(srt, &c__1, &c__1, &c__1) * .5f + sigma_(srt, &c__1, &c__1, &c__0) * .25f;
			sigdn = signd * .5f * renom;
			if (x1 <= (*signn + sdprod) / *sig)
			{
				goto L108;
			}
			if (x1 > (*signn + sigdn + x1440 * 2.f + x1535 * 2.f + *sigk +
					  sdprod) /
						 *sig)
			{
				return 0;
			}

			if (*sigk / (*sigk + sigdn + x1440 * 2 + x1535 * 2) > ranart_(&rndf77_1.nseed))
			{
				goto L306;
			}
			if (ranart_(&rndf77_1.nseed) < sigdn / (sigdn + x1440 * 2.f +
													x1535 * 2.f))
			{
				m12 = 1;
				goto L206;
			}
			else
			{
				if (ranart_(&rndf77_1.nseed) < x1440 / (x1440 + x1535))
				{
					m12 = 41;
					if (ranart_(&rndf77_1.nseed) <= .5f)
					{
						m12 = 43;
					}
				}
				else
				{
					m12 = 42;
					if (ranart_(&rndf77_1.nseed) <= .5f)
					{
						m12 = 44;
					}
				}
				goto L204;
			}
		}
		if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 7)
		{
			signd = sigma_(srt, &c__1, &c__1, &c__1) * .5f + sigma_(srt, &c__1, &c__1, &c__0) * .25f;
			sigdn = signd * .5f * renom;
			if (x1 <= (*signn + sdprod) / *sig)
			{
				goto L108;
			}
			if (x1 > (*signn + sigdn + x1440 * 2.f + x1535 * 2.f + *sigk +
					  sdprod) /
						 *sig)
			{
				return 0;
			}

			if (*sigk / (*sigk + sigdn + x1440 * 2 + x1535 * 2) > ranart_(&rndf77_1.nseed))
			{
				goto L306;
			}
			if (ranart_(&rndf77_1.nseed) < sigdn / (sigdn + x1440 * 2.f +
													x1535 * 2.f))
			{
				m12 = 2;
				goto L206;
			}
			else
			{
				if (ranart_(&rndf77_1.nseed) < x1440 / (x1440 + x1535))
				{
					m12 = 50;
					if (ranart_(&rndf77_1.nseed) <= .5f)
					{
						m12 = 51;
					}
				}
				else
				{
					m12 = 52;
					if (ranart_(&rndf77_1.nseed) <= .5f)
					{
						m12 = 53;
					}
				}
				goto L204;
			}
		}
		if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 10 && ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 1 || (i__2 = ee_1.lb[*i2 - 1], abs(i__2)) == 1))
		{
			signd = sigma_(srt, &c__2, &c__0, &c__1) * .75f;
			sigdn = signd * renomn;
			if (x1 <= (*signn + sdprod) / *sig)
			{
				goto L108;
			}
			if (x1 > (*signn + sigdn + x1535 + *sigk + sdprod) / *sig)
			{
				return 0;
			}

			if (*sigk / (*sigk + sigdn + x1535) > ranart_(&rndf77_1.nseed))
			{
				goto L306;
			}
			if (ranart_(&rndf77_1.nseed) < sigdn / (sigdn + x1535))
			{
				m12 = 7;
				goto L206;
			}
			else
			{
				m12 = 54;
				if (ranart_(&rndf77_1.nseed) <= .5f)
				{
					m12 = 55;
				}
			}
			goto L204;
		}
		if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 22 && ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2 || (i__2 = ee_1.lb[*i2 - 1], abs(i__2)) == 2))
		{
			signd = sigma_(srt, &c__2, &c__0, &c__1) * .75f;
			sigdn = signd * renomn;
			if (x1 <= (*signn + sdprod) / *sig)
			{
				goto L108;
			}
			if (x1 > (*signn + sigdn + x1535 + *sigk + sdprod) / *sig)
			{
				return 0;
			}

			if (*sigk / (*sigk + sigdn + x1535) > ranart_(&rndf77_1.nseed))
			{
				goto L306;
			}
			if (ranart_(&rndf77_1.nseed) < sigdn / (sigdn + x1535))
			{
				m12 = 8;
				goto L206;
			}
			else
			{
				m12 = 56;
				if (ranart_(&rndf77_1.nseed) <= .5f)
				{
					m12 = 57;
				}
			}
			goto L204;
		}
		if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 12 || (i__2 = ee_1.lb[*i1 - 1], abs(i__2)) == 13 || (i__3 = ee_1.lb[*i2 - 1], abs(i__3)) == 12 || (i__4 = ee_1.lb[*i2 - 1], abs(i__4)) == 13)
		{
			signd = x1535;
			sigdn = signd * renom1;
			if (x1 <= (*signn + sdprod) / *sig)
			{
				goto L108;
			}
			if (x1 > (*signn + sigdn + *sigk + sdprod) / *sig)
			{
				return 0;
			}

			if (*sigk / (*sigk + sigdn) > ranart_(&rndf77_1.nseed))
			{
				goto L306;
			}
			if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 24)
			{
				m12 = 10;
			}
			if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 12)
			{
				m12 = 12;
			}
			if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 26)
			{
				m12 = 11;
			}
			if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 13)
			{
				m12 = 9;
			}
			goto L206;
		}
	L204:
		dmax__ = *srt - .9383f - .005f;
		dmin__ = 1.078f;
		if (m12 == 37 || m12 == 39 || m12 == 41 || m12 == 43 || m12 == 46 ||
			m12 == 48 || m12 == 50 || m12 == 51)
		{
			if (dmax__ < 1.44f)
			{
				fm = fns_(&dmax__, srt, &c_b172);
			}
			else
			{
				xdmass = 1.44f;
				fm = fns_(&xdmass, srt, &c_b173);
			}
			if (fm == 0.f)
			{
				fm = 1e-9f;
			}
			ntry2 = 0;
		L11:
			dm = ranart_(&rndf77_1.nseed) * (dmax__ - dmin__) + dmin__;
			++ntry2;
			if (ranart_(&rndf77_1.nseed) > fns_(&dm, srt, &c_b173) / fm &&
				ntry2 <= 10)
			{
				goto L11;
			}
			if (dm > 2.14f)
			{
				goto L11;
			}
			goto L13;
		}
		else
		{
			if (dmax__ < 1.535f)
			{
				fm = fd5_(&dmax__, srt, &c_b172);
			}
			else
			{
				xdmass = 1.535f;
				fm = fd5_(&xdmass, srt, &c_b173);
			}
			if (fm == 0.f)
			{
				fm = 1e-9f;
			}
			ntry1 = 0;
		L12:
			dm = ranart_(&rndf77_1.nseed) * (dmax__ - dmin__) + dmin__;
			++ntry1;
			if (ranart_(&rndf77_1.nseed) > fd5_(&dm, srt, &c_b173) / fm &&
				ntry1 <= 10)
			{
				goto L12;
			}
			if (dm > 1.84f)
			{
				goto L12;
			}
		}
	L13:
		prf = 0.f;
		r__2 = *srt;
		r__3 = dm;
		r__1 = (r__2 * r__2 - r__3 * r__3 + .88040689000000005f) / (*srt *
																	2.f);
		pf2 = r__1 * r__1 - .88040689000000005f;
		if (pf2 > 0.f)
		{
			prf = sqrt(pf2);
		}
		if (m12 == 37)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 9)
			{
				ee_1.lb[*i1 - 1] = 1;
				cc_1.e[*i1 - 1] = .93828f;
				ee_1.lb[*i2 - 1] = 11;
				cc_1.e[*i2 - 1] = dm;
			}
			else
			{
				ee_1.lb[*i2 - 1] = 1;
				cc_1.e[*i2 - 1] = .93828f;
				ee_1.lb[*i1 - 1] = 11;
				cc_1.e[*i1 - 1] = dm;
			}
			goto L207;
		}
		if (m12 == 38)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 9)
			{
				ee_1.lb[*i1 - 1] = 1;
				cc_1.e[*i1 - 1] = .93828f;
				ee_1.lb[*i2 - 1] = 13;
				cc_1.e[*i2 - 1] = dm;
			}
			else
			{
				ee_1.lb[*i2 - 1] = 1;
				cc_1.e[*i2 - 1] = .93828f;
				ee_1.lb[*i1 - 1] = 13;
				cc_1.e[*i1 - 1] = dm;
			}
			goto L207;
		}
		if (m12 == 39)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 8)
			{
				ee_1.lb[*i1 - 1] = 1;
				cc_1.e[*i1 - 1] = .93828f;
				ee_1.lb[*i2 - 1] = 11;
				cc_1.e[*i2 - 1] = dm;
			}
			else
			{
				ee_1.lb[*i2 - 1] = 1;
				cc_1.e[*i2 - 1] = .93828f;
				ee_1.lb[*i1 - 1] = 11;
				cc_1.e[*i1 - 1] = dm;
			}
			goto L207;
		}
		if (m12 == 40)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 8)
			{
				ee_1.lb[*i1 - 1] = 1;
				cc_1.e[*i1 - 1] = .93828f;
				ee_1.lb[*i2 - 1] = 13;
				cc_1.e[*i2 - 1] = dm;
			}
			else
			{
				ee_1.lb[*i2 - 1] = 1;
				cc_1.e[*i2 - 1] = .93828f;
				ee_1.lb[*i1 - 1] = 13;
				cc_1.e[*i1 - 1] = dm;
			}
			goto L207;
		}
		if (m12 == 41)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 8)
			{
				ee_1.lb[*i1 - 1] = 2;
				cc_1.e[*i1 - 1] = .939457f;
				ee_1.lb[*i2 - 1] = 11;
				cc_1.e[*i2 - 1] = dm;
			}
			else
			{
				ee_1.lb[*i2 - 1] = 2;
				cc_1.e[*i2 - 1] = .939457f;
				ee_1.lb[*i1 - 1] = 11;
				cc_1.e[*i1 - 1] = dm;
			}
			goto L207;
		}
		if (m12 == 42)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 8)
			{
				ee_1.lb[*i1 - 1] = 2;
				cc_1.e[*i1 - 1] = .939457f;
				ee_1.lb[*i2 - 1] = 13;
				cc_1.e[*i2 - 1] = dm;
			}
			else
			{
				ee_1.lb[*i2 - 1] = 2;
				cc_1.e[*i2 - 1] = .939457f;
				ee_1.lb[*i1 - 1] = 13;
				cc_1.e[*i1 - 1] = dm;
			}
			goto L207;
		}
		if (m12 == 43)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 8)
			{
				ee_1.lb[*i1 - 1] = 1;
				cc_1.e[*i1 - 1] = .93828f;
				ee_1.lb[*i2 - 1] = 10;
				cc_1.e[*i2 - 1] = dm;
			}
			else
			{
				ee_1.lb[*i2 - 1] = 1;
				cc_1.e[*i2 - 1] = .93828f;
				ee_1.lb[*i1 - 1] = 10;
				cc_1.e[*i1 - 1] = dm;
			}
			goto L207;
		}
		if (m12 == 44)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 8)
			{
				ee_1.lb[*i1 - 1] = 1;
				cc_1.e[*i1 - 1] = .93828f;
				ee_1.lb[*i2 - 1] = 12;
				cc_1.e[*i2 - 1] = dm;
			}
			else
			{
				ee_1.lb[*i2 - 1] = 1;
				cc_1.e[*i2 - 1] = .93828f;
				ee_1.lb[*i1 - 1] = 12;
				cc_1.e[*i1 - 1] = dm;
			}
			goto L207;
		}
		if (m12 == 46)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 6)
			{
				ee_1.lb[*i1 - 1] = 2;
				cc_1.e[*i1 - 1] = .939457f;
				ee_1.lb[*i2 - 1] = 10;
				cc_1.e[*i2 - 1] = dm;
			}
			else
			{
				ee_1.lb[*i2 - 1] = 2;
				cc_1.e[*i2 - 1] = .939457f;
				ee_1.lb[*i1 - 1] = 10;
				cc_1.e[*i1 - 1] = dm;
			}
			goto L207;
		}
		if (m12 == 47)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 6)
			{
				ee_1.lb[*i1 - 1] = 2;
				cc_1.e[*i1 - 1] = .939457f;
				ee_1.lb[*i2 - 1] = 12;
				cc_1.e[*i2 - 1] = dm;
			}
			else
			{
				ee_1.lb[*i2 - 1] = 2;
				cc_1.e[*i2 - 1] = .939457f;
				ee_1.lb[*i1 - 1] = 12;
				cc_1.e[*i1 - 1] = dm;
			}
			goto L207;
		}
		if (m12 == 48)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 7)
			{
				ee_1.lb[*i1 - 1] = 2;
				cc_1.e[*i1 - 1] = .939457f;
				ee_1.lb[*i2 - 1] = 11;
				cc_1.e[*i2 - 1] = dm;
			}
			else
			{
				ee_1.lb[*i2 - 1] = 2;
				cc_1.e[*i2 - 1] = .939457f;
				ee_1.lb[*i1 - 1] = 11;
				cc_1.e[*i1 - 1] = dm;
			}
			goto L207;
		}
		if (m12 == 49)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 7)
			{
				ee_1.lb[*i1 - 1] = 2;
				cc_1.e[*i1 - 1] = .939457f;
				ee_1.lb[*i2 - 1] = 12;
				cc_1.e[*i2 - 1] = dm;
			}
			else
			{
				ee_1.lb[*i2 - 1] = 2;
				cc_1.e[*i2 - 1] = .939457f;
				ee_1.lb[*i1 - 1] = 12;
				cc_1.e[*i1 - 1] = dm;
			}
			goto L207;
		}
		if (m12 == 50)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 7)
			{
				ee_1.lb[*i1 - 1] = 1;
				cc_1.e[*i1 - 1] = .93828f;
				ee_1.lb[*i2 - 1] = 10;
				cc_1.e[*i2 - 1] = dm;
			}
			else
			{
				ee_1.lb[*i2 - 1] = 1;
				cc_1.e[*i2 - 1] = .93828f;
				ee_1.lb[*i1 - 1] = 10;
				cc_1.e[*i1 - 1] = dm;
			}
			goto L207;
		}
		if (m12 == 51)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 7)
			{
				ee_1.lb[*i1 - 1] = 2;
				cc_1.e[*i1 - 1] = .939457f;
				ee_1.lb[*i2 - 1] = 11;
				cc_1.e[*i2 - 1] = dm;
			}
			else
			{
				ee_1.lb[*i2 - 1] = 2;
				cc_1.e[*i2 - 1] = .939457f;
				ee_1.lb[*i1 - 1] = 11;
				cc_1.e[*i1 - 1] = dm;
			}
			goto L207;
		}
		if (m12 == 52)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 7)
			{
				ee_1.lb[*i1 - 1] = 1;
				cc_1.e[*i1 - 1] = .93828f;
				ee_1.lb[*i2 - 1] = 12;
				cc_1.e[*i2 - 1] = dm;
			}
			else
			{
				ee_1.lb[*i2 - 1] = 1;
				cc_1.e[*i2 - 1] = .93828f;
				ee_1.lb[*i1 - 1] = 12;
				cc_1.e[*i1 - 1] = dm;
			}
			goto L207;
		}
		if (m12 == 53)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 7)
			{
				ee_1.lb[*i1 - 1] = 2;
				cc_1.e[*i1 - 1] = .939457f;
				ee_1.lb[*i2 - 1] = 13;
				cc_1.e[*i2 - 1] = dm;
			}
			else
			{
				ee_1.lb[*i2 - 1] = 2;
				cc_1.e[*i2 - 1] = .939457f;
				ee_1.lb[*i1 - 1] = 13;
				cc_1.e[*i1 - 1] = dm;
			}
			goto L207;
		}
		if (m12 == 54)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 10)
			{
				ee_1.lb[*i1 - 1] = 2;
				cc_1.e[*i1 - 1] = .939457f;
				ee_1.lb[*i2 - 1] = 13;
				cc_1.e[*i2 - 1] = dm;
			}
			else
			{
				ee_1.lb[*i2 - 1] = 2;
				cc_1.e[*i2 - 1] = .939457f;
				ee_1.lb[*i1 - 1] = 13;
				cc_1.e[*i1 - 1] = dm;
			}
			goto L207;
		}
		if (m12 == 55)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 10)
			{
				ee_1.lb[*i1 - 1] = 1;
				cc_1.e[*i1 - 1] = .93828f;
				ee_1.lb[*i2 - 1] = 12;
				cc_1.e[*i2 - 1] = dm;
			}
			else
			{
				ee_1.lb[*i2 - 1] = 1;
				cc_1.e[*i2 - 1] = .93828f;
				ee_1.lb[*i1 - 1] = 12;
				cc_1.e[*i1 - 1] = dm;
			}
			goto L207;
		}
		if (m12 == 56)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 11)
			{
				ee_1.lb[*i1 - 1] = 2;
				cc_1.e[*i1 - 1] = .939457f;
				ee_1.lb[*i2 - 1] = 13;
				cc_1.e[*i2 - 1] = dm;
			}
			else
			{
				ee_1.lb[*i2 - 1] = 2;
				cc_1.e[*i2 - 1] = .939457f;
				ee_1.lb[*i1 - 1] = 13;
				cc_1.e[*i1 - 1] = dm;
			}
			goto L207;
		}
		if (m12 == 57)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 11)
			{
				ee_1.lb[*i1 - 1] = 1;
				cc_1.e[*i1 - 1] = .93828f;
				ee_1.lb[*i2 - 1] = 12;
				cc_1.e[*i2 - 1] = dm;
			}
			else
			{
				ee_1.lb[*i2 - 1] = 1;
				cc_1.e[*i2 - 1] = .93828f;
				ee_1.lb[*i1 - 1] = 12;
				cc_1.e[*i1 - 1] = dm;
			}
		}
		goto L207;
	L206:
		if (m12 == 1)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 8)
			{
				ee_1.lb[*i2 - 1] = 2;
				ee_1.lb[*i1 - 1] = 1;
				cc_1.e[*i1 - 1] = .93828f;
			}
			else
			{
				ee_1.lb[*i1 - 1] = 2;
				ee_1.lb[*i2 - 1] = 1;
				cc_1.e[*i2 - 1] = .93828f;
			}
			goto L207;
		}
		if (m12 == 2)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 7)
			{
				ee_1.lb[*i2 - 1] = 1;
				ee_1.lb[*i1 - 1] = 2;
				cc_1.e[*i1 - 1] = .939457f;
			}
			else
			{
				ee_1.lb[*i1 - 1] = 1;
				ee_1.lb[*i2 - 1] = 2;
				cc_1.e[*i2 - 1] = .939457f;
			}
			goto L207;
		}
		if (m12 == 3)
		{
			ee_1.lb[*i1 - 1] = 1;
			ee_1.lb[*i2 - 1] = 1;
			cc_1.e[*i1 - 1] = .93828f;
			cc_1.e[*i2 - 1] = .93828f;
			goto L207;
		}
		if (m12 == 4)
		{
			ee_1.lb[*i1 - 1] = 1;
			ee_1.lb[*i2 - 1] = 1;
			cc_1.e[*i1 - 1] = .93828f;
			cc_1.e[*i2 - 1] = .93828f;
			goto L207;
		}
		if (m12 == 5)
		{
			ee_1.lb[*i1 - 1] = 2;
			ee_1.lb[*i2 - 1] = 2;
			cc_1.e[*i1 - 1] = .939457f;
			cc_1.e[*i2 - 1] = .939457f;
			goto L207;
		}
		if (m12 == 6)
		{
			ee_1.lb[*i1 - 1] = 2;
			ee_1.lb[*i2 - 1] = 2;
			cc_1.e[*i1 - 1] = .939457f;
			cc_1.e[*i2 - 1] = .939457f;
			goto L207;
		}
		if (m12 == 7)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 1)
			{
				ee_1.lb[*i1 - 1] = 1;
				ee_1.lb[*i2 - 1] = 2;
				cc_1.e[*i1 - 1] = .93828f;
				cc_1.e[*i2 - 1] = .939457f;
			}
			else
			{
				ee_1.lb[*i1 - 1] = 2;
				ee_1.lb[*i2 - 1] = 1;
				cc_1.e[*i1 - 1] = .939457f;
				cc_1.e[*i2 - 1] = .93828f;
			}
			goto L207;
		}
		if (m12 == 8)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2)
			{
				ee_1.lb[*i1 - 1] = 2;
				ee_1.lb[*i2 - 1] = 1;
				cc_1.e[*i1 - 1] = .939457f;
				cc_1.e[*i2 - 1] = .93828f;
			}
			else
			{
				ee_1.lb[*i1 - 1] = 1;
				ee_1.lb[*i2 - 1] = 2;
				cc_1.e[*i1 - 1] = .93828f;
				cc_1.e[*i2 - 1] = .939457f;
			}
			goto L207;
		}
		if (m12 == 9)
		{
			ee_1.lb[*i1 - 1] = 1;
			ee_1.lb[*i2 - 1] = 1;
			cc_1.e[*i1 - 1] = .93828f;
			cc_1.e[*i2 - 1] = .93828f;
			goto L207;
		}
		if (m12 == 12)
		{
			ee_1.lb[*i1 - 1] = 2;
			ee_1.lb[*i2 - 1] = 1;
			cc_1.e[*i1 - 1] = .939457f;
			cc_1.e[*i2 - 1] = .93828f;
			goto L207;
		}
		if (m12 == 11)
		{
			ee_1.lb[*i1 - 1] = 2;
			ee_1.lb[*i2 - 1] = 1;
			cc_1.e[*i1 - 1] = .939457f;
			cc_1.e[*i2 - 1] = .93828f;
			goto L207;
		}
		if (m12 == 12)
		{
			ee_1.lb[*i1 - 1] = 1;
			ee_1.lb[*i2 - 1] = 2;
			cc_1.e[*i1 - 1] = .93828f;
			cc_1.e[*i2 - 1] = .939457f;
		}
	L207:
		pr = prf;
		c1 = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
		if (*srt <= 2.14f)
		{
			c1 = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
		}
		if (*srt > 2.14f && *srt <= 2.4f)
		{
			c1 = ang_(srt, &input1_1.iseed);
		}
		if (*srt > 2.4f)
		{
			xptr = pr * .33f;
			cc1 = ptr_(&xptr, &input1_1.iseed);
			r__1 = pr;
			r__2 = cc1;
			c1 = sqrt(r__1 * r__1 - r__2 * r__2) / pr;
		}
		t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
		*iblock = 3;
	}
	if (ianti == 1 && ee_1.lb[*i1 - 1] >= 1 && ee_1.lb[*i2 - 1] >= 1)
	{
		ee_1.lb[*i1 - 1] = -ee_1.lb[*i1 - 1];
		ee_1.lb[*i2 - 1] = -ee_1.lb[*i2 - 1];
	}
L107:
	if (*px == 0.f && *py == 0.f)
	{
		t2 = 0.f;
	}
	else
	{
		t2 = atan2(*py, *px);
	}
	r__1 = c1;
	s1 = sqrt(1.f - r__1 * r__1);
	r__1 = c2;
	s2 = sqrt(1.f - r__1 * r__1);
	ct1 = cos(t1);
	st1 = sin(t1);
	ct2 = cos(t2);
	st2 = sin(t2);
	*pz = pr * (c1 * c2 - s1 * s2 * ct1);
	ss = c2 * s1 * ct1 + s2 * c1;
	*px = pr * (ss * ct2 - s1 * st1 * st2);
	*py = pr * (ss * st2 + s1 * st1 * ct2);
	return 0;
L306:
	if (*xsk5 / *sigk > ranart_(&rndf77_1.nseed))
	{
		leadng_1.pz1 = bb_1.p[*i1 * 3 - 1];
		pz2 = bb_1.p[*i2 * 3 - 1];
		ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 2) + 1;
		ee_1.lb[*i2 - 1] = (int)(ranart_(&rndf77_1.nseed) * 2) + 1;
		++nn_1.nnn;
		pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 29;
		pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = 1.02f;
		*iblock = 222;
		goto L208;
	}
	*iblock = 11;
	if (ianti == 1)
	{
		*iblock = -11;
	}

	leadng_1.pz1 = bb_1.p[*i1 * 3 - 1];
	pz2 = bb_1.p[*i2 * 3 - 1];
	++nn_1.nnn;
	pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 23;
	pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = .498f;
	if (*srt <= 2.63f)
	{
		ic = 1;
		ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 2) + 1;
		ee_1.lb[*i2 - 1] = 14;
		goto L208;
	}
	if (*srt <= 2.74f && *srt > 2.63f)
	{
		if (*xsk1 / (*xsk1 + *xsk2) > ranart_(&rndf77_1.nseed))
		{
			ic = 1;
			ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 2) + 1;
			ee_1.lb[*i2 - 1] = 14;
		}
		else
		{
			ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 2) + 1;
			ee_1.lb[*i2 - 1] = (int)(ranart_(&rndf77_1.nseed) * 3) + 15;
			ic = 2;
		}
		goto L208;
	}
	if (*srt <= 2.77f && *srt > 2.74f)
	{
		if (*xsk1 / (*xsk1 + *xsk2 + *xsk3) > ranart_(&rndf77_1.nseed))
		{
			ic = 1;
			ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 2) + 1;
			ee_1.lb[*i2 - 1] = 14;
			goto L208;
		}
		else
		{
			if (*xsk2 / (*xsk2 + *xsk3) > ranart_(&rndf77_1.nseed))
			{
				ic = 2;
				ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 2) +
								   1;
				ee_1.lb[*i2 - 1] = (int)(ranart_(&rndf77_1.nseed) * 3) +
								   15;
			}
			else
			{
				ic = 3;
				ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 4) +
								   6;
				ee_1.lb[*i2 - 1] = 14;
			}
			goto L208;
		}
	}
	if (*srt > 2.77f)
	{
		if (*xsk1 / (*xsk1 + *xsk2 + *xsk3 + *xsk4) > ranart_(&rndf77_1.nseed))
		{
			ic = 1;
			ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 2) + 1;
			ee_1.lb[*i2 - 1] = 14;
			goto L208;
		}
		else
		{
			if (*xsk3 / (*xsk2 + *xsk3 + *xsk4) > ranart_(&rndf77_1.nseed))
			{
				ic = 3;
				ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 4) +
								   6;
				ee_1.lb[*i2 - 1] = 14;
				goto L208;
			}
			else
			{
				if (*xsk2 / (*xsk2 + *xsk4) > ranart_(&rndf77_1.nseed))
				{
					ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) *
												 2) +
									   1;
					ee_1.lb[*i2 - 1] = (int)(ranart_(&rndf77_1.nseed) *
												 3) +
									   15;
					ic = 2;
				}
				else
				{
					ic = 4;
					ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) *
												 4) +
									   6;
					ee_1.lb[*i2 - 1] = (int)(ranart_(&rndf77_1.nseed) *
												 3) +
									   15;
				}
				goto L208;
			}
		}
	}
L208:
	if (ianti == 1 && ee_1.lb[*i1 - 1] >= 1 && ee_1.lb[*i2 - 1] >= 1)
	{
		ee_1.lb[*i1 - 1] = -ee_1.lb[*i1 - 1];
		ee_1.lb[*i2 - 1] = -ee_1.lb[*i2 - 1];
		if (pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] == 23)
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 21;
		}
	}
	lbi1 = ee_1.lb[*i1 - 1];
	lbi2 = ee_1.lb[*i2 - 1];
	ntry1 = 0;
L128:
	bbkaon_(&ic, srt, &px3, &py3, &pz3, &dm3, &px4, &py4, &pz4, &dm4, &ppx, &ppy, &ppz, &icou1);
	++ntry1;
	if (icou1 < 0 && ntry1 <= 20)
	{
		goto L128;
	}
	rotate_(px, py, pz, &px3, &py3, &pz3);
	rotate_(px, py, pz, &px4, &py4, &pz4);
	rotate_(px, py, pz, &ppx, &ppy, &ppz);
	r__1 = dm3;
	r__2 = px3;
	r__3 = py3;
	r__4 = pz3;
	e1cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	p1beta = px3 * bg_1.betax + py3 * bg_1.betay + pz3 * bg_1.betaz;
	transf = bg_1.gamma * (bg_1.gamma * p1beta / (bg_1.gamma + 1) + e1cm);
	pt1i1 = bg_1.betax * transf + px3;
	pt2i1 = bg_1.betay * transf + py3;
	pt3i1 = bg_1.betaz * transf + pz3;
	eti1 = dm3;
	r__1 = dm4;
	r__2 = px4;
	r__3 = py4;
	r__4 = pz4;
	e2cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	p2beta = px4 * bg_1.betax + py4 * bg_1.betay + pz4 * bg_1.betaz;
	transf = bg_1.gamma * (bg_1.gamma * p2beta / (bg_1.gamma + 1.f) + e2cm);
	pt1i2 = bg_1.betax * transf + px4;
	pt2i2 = bg_1.betay * transf + py4;
	pt3i2 = bg_1.betaz * transf + pz4;
	eti2 = dm4;
	r__1 = ppx;
	r__2 = ppy;
	r__3 = ppz;
	epcm = sqrt(r__1 * r__1 + .248004f + r__2 * r__2 + r__3 * r__3);
	ppbeta = ppx * bg_1.betax + ppy * bg_1.betay + ppz * bg_1.betaz;
	transf = bg_1.gamma * (bg_1.gamma * ppbeta / (bg_1.gamma + 1.f) + epcm);
	pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = bg_1.betax *
															   transf +
														   ppx;
	pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = bg_1.betay *
															   transf +
														   ppy;
	pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = bg_1.betaz *
															   transf +
														   ppz;
	dpert_1.dppion[nn_1.nnn + *irun * 150001 - 150002] = dpert_1.dpertp[*i1 -
																		1] *
														 dpert_1.dpertp[*i2 - 1];
	pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = aa_1.r__[*i1 * 3 -
																	3];
	pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = aa_1.r__[*i1 * 3 -
																	2];
	pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = aa_1.r__[*i1 * 3 -
																	1];

	bb_1.p[*i1 * 3 - 3] = pt1i1;
	bb_1.p[*i1 * 3 - 2] = pt2i1;
	bb_1.p[*i1 * 3 - 1] = pt3i1;
	cc_1.e[*i1 - 1] = eti1;
	ee_1.lb[*i1 - 1] = lbi1;
	bb_1.p[*i2 * 3 - 3] = pt1i2;
	bb_1.p[*i2 * 3 - 2] = pt2i2;
	bb_1.p[*i2 * 3 - 1] = pt3i2;
	cc_1.e[*i2 - 1] = eti2;
	ee_1.lb[*i2 - 1] = lbi2;
	leadng_1.px1 = bb_1.p[*i1 * 3 - 3];
	leadng_1.py1 = bb_1.p[*i1 * 3 - 2];
	leadng_1.pz1 = bb_1.p[*i1 * 3 - 1];
	leadng_1.em1 = cc_1.e[*i1 - 1];
	ee_1.id[*i1 - 1] = 2;
	ee_1.id[*i2 - 1] = 2;
	id1 = ee_1.id[*i1 - 1];
	if (pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] != 29)
	{
		*iblock = 11;
	}
	leadng_1.lb1 = ee_1.lb[*i1 - 1];
	dpi_1.lb2 = ee_1.lb[*i2 - 1];
	am1 = leadng_1.em1;
	am2 = dpi_1.em2;
	r__1 = leadng_1.em1;
	r__2 = leadng_1.px1;
	r__3 = leadng_1.py1;
	r__4 = leadng_1.pz1;
	leadng_1.e1 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	return 0;
L108:
	if (para8_1.idpert == 1 && *ipert1 == 1 && para8_1.npertd >= 1)
	{
		ndloop = para8_1.npertd;
	}
	else if (para8_1.idpert == 2 && para8_1.npertd >= 1)
	{
		ndloop = para8_1.npertd + 1;
	}
	else
	{
		ndloop = 1;
	}

	dprob1 = sdprod / *sig / (float)para8_1.npertd;
	i__1 = ndloop;
	for (idloop = 1; idloop <= i__1; ++idloop)
	{
		bbdangle_(&pxd, &pyd, &pzd, nt, ipert1, &ianti, &idloop, &pfinal, &dprob1, &lbm);
		rotate_(px, py, pz, &pxd, &pyd, &pzd);
		xmass = 1.8756f;
		r__1 = xmass;
		r__2 = pxd;
		r__3 = pyd;
		r__4 = pzd;
		e1dcm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
		p1dbeta = pxd * bg_1.betax + pyd * bg_1.betay + pzd * bg_1.betaz;
		transf = bg_1.gamma * (bg_1.gamma * p1dbeta / (bg_1.gamma + 1.f) +
							   e1dcm);
		pxi1 = bg_1.betax * transf + pxd;
		pyi1 = bg_1.betay * transf + pyd;
		pzi1 = bg_1.betaz * transf + pzd;
		if (ianti == 0)
		{
			lbd = 42;
		}
		else
		{
			lbd = -42;
		}
		if (para8_1.idpert == 1 && *ipert1 == 1 && para8_1.npertd >= 1)
		{
			++nn_1.nnn;
			pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = pxi1;
			pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = pyi1;
			pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = pzi1;
			pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = 1.8756f;
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = lbd;
			pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = aa_1.r__[*i1 * 3 - 3];
			pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = aa_1.r__[*i1 * 3 - 2];
			pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = aa_1.r__[*i1 * 3 - 1];
			dpert_1.dppion[nn_1.nnn + *irun * 150001 - 150002] = sdprod / *sig / (float)para8_1.npertd;
		}
		else if (para8_1.idpert == 2 && idloop <= para8_1.npertd)
		{
			ppd[idloop * 3 - 3] = pxi1;
			ppd[idloop * 3 - 2] = pyi1;
			ppd[idloop * 3 - 1] = pzi1;
			lbpd[idloop - 1] = lbd;
		}
		else
		{
			cc_1.e[*i1 - 1] = xmm;
			r__1 = xmm;
			r__2 = pxd;
			r__3 = pyd;
			r__4 = pzd;
			e2picm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
			p2pibeta = -pxd * bg_1.betax - pyd * bg_1.betay - pzd * bg_1.betaz;
			transf = bg_1.gamma * (bg_1.gamma * p2pibeta / (bg_1.gamma + 1.f) + e2picm);
			pxi2 = bg_1.betax * transf - pxd;
			pyi2 = bg_1.betay * transf - pyd;
			pzi2 = bg_1.betaz * transf - pzd;
			bb_1.p[*i1 * 3 - 3] = pxi2;
			bb_1.p[*i1 * 3 - 2] = pyi2;
			bb_1.p[*i1 * 3 - 1] = pzi2;
			ee_1.lb[*i1 - 1] = lbm;
			leadng_1.px1 = bb_1.p[*i1 * 3 - 3];
			leadng_1.py1 = bb_1.p[*i1 * 3 - 2];
			leadng_1.pz1 = bb_1.p[*i1 * 3 - 1];
			leadng_1.em1 = cc_1.e[*i1 - 1];
			ee_1.id[*i1 - 1] = 2;
			id1 = ee_1.id[*i1 - 1];
			r__1 = leadng_1.em1;
			r__2 = leadng_1.px1;
			r__3 = leadng_1.py1;
			r__4 = leadng_1.pz1;
			leadng_1.e1 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
			leadng_1.lb1 = ee_1.lb[*i1 - 1];
			bb_1.p[*i2 * 3 - 3] = pxi1;
			bb_1.p[*i2 * 3 - 2] = pyi1;
			bb_1.p[*i2 * 3 - 1] = pzi1;
			ee_1.lb[*i2 - 1] = lbd;
			dpi_1.lb2 = ee_1.lb[*i2 - 1];
			cc_1.e[*i2 - 1] = 1.8756f;
			eti2 = cc_1.e[*i2 - 1];
			ee_1.id[*i2 - 1] = 2;
			if (para8_1.idpert == 2 && idloop == ndloop)
			{
				i__2 = para8_1.npertd;
				for (ipertd = 1; ipertd <= i__2; ++ipertd)
				{
					++nn_1.nnn;
					pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450006] =
						ppd[ipertd * 3 - 3];
					pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450005] =
						ppd[ipertd * 3 - 2];
					pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450004] =
						ppd[ipertd * 3 - 1];
					pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = 1.8756f;
					pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = lbpd[ipertd - 1];
					pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450006] =
						aa_1.r__[*i1 * 3 - 3];
					pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450005] =
						aa_1.r__[*i1 * 3 - 2];
					pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450004] =
						aa_1.r__[*i1 * 3 - 1];
					dpert_1.dppion[nn_1.nnn + *irun * 150001 - 150002] = 1.f /
																		 (float)para8_1.npertd;
				}
			}
		}
	}
	*iblock = 501;
	return 0;
}

int crdd_(int *irun, float *px, float *py, float *pz, float *srt, int *i1, int *i2, int *iblock, int *ntag, float *signn, float *sig, int *nt, int *ipert1)
{
	int i__1, i__2, i__3, i__4, i__5, i__6;
	float r__1, r__2, r__3, r__4, r__5;

	double sqrt(double), atan2(double, double), exp(double),
		log(double), cos(double), sin(double);

	extern int bbdangle_(float *, float *, float *, int *,
						 int *, int *, int *, float *, float *, int *);
	static float a, x, c1, c2, s1, t1, t2, s2, x1;
	static int ic;
	static float al;
	static int m12, n12;
	static float dm, fm, as, ta, es, pr, ss, cc1, ak0;
	extern double fd5_(float *, float *, float *);
	static float dm3, dm4, ct1, s2d, ct2;
	static int id1;
	static float am1, am2, pr2, st1, st2, pz2, px3, py3, pz3, px4, py4, pz4,
		ada, ana;
	static int idd, lbd, ich;
	extern double ang_(float *, int *);
	static int lbm;
	extern int n1535_(int *, int *, float *, float *);
	static float akp, ppd[30000], x1535;
	extern double fns_(float *, float *, float *);
	static float pxd, xmm, pyd, pzd;
	extern double ptr_(float *, int *);
	static float ppx, ppy, ppz, e1cm, e2cm;
	static int lbi1, lbi2;
	static float eti1, sig2, eti2;
	extern double ppk0_(float *), ppk1_(float *);
	static float pxi1, pyi1, pzi1, pxi2, pyi2, xsk1, xsk2, xsk3, xsk4, xsk5,
		pzi2;
	static int lbpd[10000];
	static float epcm, dmin__, dmax__, sigk, pt1i1, pt2i1, pt3i1, pt1i2, pt2i2,
		pt3i2, e1dcm, xptr, t1dlk, t2dlk;
	static int icou1;
	static float t1dsk, t2dsk, t1nlk, t1nsk, t2nsk;
	static int ntry1, ntry2;
	extern int sbbdm_(float *, float *, int *, int *,
					  float *, float *);
	extern double sigma_(float *, int *, int *, int *);
	static int ianti;
	static float pmdlk, signd, pmdsk;
	extern double pplpk_(float *);
	static float pmnsk, xmass;
	extern double reab2d_(int *, int *, float *);
	static float p1beta, p2beta, e2picm, dprob1, pmdlk2, pmdsk2, pmnsk2;
	extern int bbkaon_(int *, float *, float *, float *,
					   float *, float *, float *, float *, float *, float *, float *, float *,
					   float *, int *);
	static float pfinal, ppbeta;
	extern double ranart_(int *);
	static float sdprod, xdmass;
	extern int rotate_(float *, float *, float *, float *, float *, float *);
	static float transf;
	static int ndloop, idloop, ipertd;
	static float p1dbeta, p2pibeta;

	n12 = 0;
	m12 = 0;
	*iblock = 0;
	*ntag = 0;
	leadng_1.em1 = cc_1.e[*i1 - 1];
	dpi_1.em2 = cc_1.e[*i2 - 1];
	r__1 = *px;
	r__2 = *py;
	r__3 = *pz;
	pr = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
	c2 = *pz / pr;
	if (*px == 0.f && *py == 0.f)
	{
		t2 = 0.f;
	}
	else
	{
		t2 = atan2(*py, *px);
	}
	x1 = ranart_(&rndf77_1.nseed);
	ianti = 0;
	if (ee_1.lb[*i1 - 1] < 0 && ee_1.lb[*i2 - 1] < 0)
	{
		ianti = 1;
	}
	sbbdm_(srt, &sdprod, &ianti, &lbm, &xmm, &pfinal);
	if (para8_1.idpert == 1 && *ipert1 == 1)
	{
		if (*srt < 2.012f)
		{
			return 0;
		}
		if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) >= 6 && (i__2 = ee_1.lb[*i1 - 1], abs(i__2)) <= 13 && ((i__3 = ee_1.lb[*i2 - 1], abs(i__3)) >= 6 && (i__4 = ee_1.lb[*i2 - 1], abs(i__4)) <= 13))
		{
			goto L108;
		}
		else
		{
			return 0;
		}
	}
	if (x1 <= *signn / *sig)
	{
		r__1 = (*srt - 1.8766f) * 3.65f, r__1 *= r__1;
		as = r__1 * (r__1 * r__1);
		a = as * 6.f / (as + 1.f);
		r__1 = pr;
		ta = r__1 * r__1 * -2.f;
		x = ranart_(&rndf77_1.nseed);
		t1 = (float)log((double)(1.f - x) * exp((double)a * (double)ta) + (double)x) / a;
		c1 = 1.f - t1 / ta;
		t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
		*iblock = 20;
		goto L107;
	}
	else
	{
		if (*srt < 2.15f)
		{
			return 0;
		}
		i__3 = (i__1 = ee_1.lb[*i1 - 1], abs(i__1));
		i__4 = (i__2 = ee_1.lb[*i2 - 1], abs(i__2));
		n1535_(&i__3, &i__4, srt, &x1535);
		akp = .498f;
		ak0 = .498f;
		ana = .938f;
		ada = 1.232f;
		al = 1.1157f;
		as = 1.1197f;
		xsk1 = 0.f;
		xsk2 = 0.f;
		xsk3 = 0.f;
		xsk4 = 0.f;
		xsk5 = 0.f;
		t1nlk = ana + al + akp;
		if (*srt <= t1nlk)
		{
			goto L222;
		}
		xsk1 = pplpk_(srt) * 1.5f;
		t1dlk = ada + al + akp;
		t2dlk = ada + al - akp;
		if (*srt <= t1dlk)
		{
			goto L222;
		}
		es = *srt;
		r__1 = es;
		r__2 = t1dlk;
		r__3 = es;
		r__4 = t2dlk;
		r__5 = es;
		pmdlk2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / (r__5 * r__5 * 4.f);
		pmdlk = sqrt(pmdlk2);
		xsk3 = pplpk_(srt) * 1.5f;
		t1nsk = ana + as + akp;
		t2nsk = ana + as - akp;
		if (*srt <= t1nsk)
		{
			goto L222;
		}
		r__1 = es;
		r__2 = t1nsk;
		r__3 = es;
		r__4 = t2nsk;
		r__5 = es;
		pmnsk2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / (r__5 * r__5 * 4.f);
		pmnsk = sqrt(pmnsk2);
		xsk2 = (ppk1_(srt) + ppk0_(srt)) * 1.5f;
		t1dsk = ada + as + akp;
		t2dsk = ada + as - akp;
		if (*srt <= t1dsk)
		{
			goto L222;
		}
		r__1 = es;
		r__2 = t1dsk;
		r__3 = es;
		r__4 = t2dsk;
		r__5 = es;
		pmdsk2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / (r__5 * r__5 * 4.f);
		pmdsk = sqrt(pmdsk2);
		xsk4 = (ppk1_(srt) + ppk0_(srt)) * 1.5f;
		if (*srt <= 2.898914f)
		{
			goto L222;
		}
		xsk5 = 1e-4f;
	L222:
		sigk = xsk1 + xsk2 + xsk3 + xsk4;
		xsk1 *= 2.f;
		xsk2 *= 2.f;
		xsk3 *= 2.f;
		xsk4 *= 2.f;
		sigk = sigk * 2.f + xsk5;
		s2d = reab2d_(i1, i2, srt);
		s2d = 0.f;
		if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) >= 12 && (i__2 = ee_1.lb[*i2 - 1], abs(i__2)) >= 12 || (i__3 = ee_1.lb[*i1 - 1], abs(i__3)) >= 12 && (i__4 = ee_1.lb[*i2 - 1], abs(i__4)) >= 6 || (i__5 = ee_1.lb[*i2 - 1], abs(i__5)) >= 12 && (i__6 = ee_1.lb[*i1 - 1], abs(i__6)) >= 6)
		{
			signd = sigk + s2d;
			if (x1 <= (*signn + sdprod) / *sig)
			{
				goto L108;
			}
			if (x1 > (signd + *signn + sdprod) / *sig)
			{
				return 0;
			}

			if ((sigk + sdprod) / *sig >= ranart_(&rndf77_1.nseed))
			{
				goto L306;
			}

			goto L1012;
		}
		idd = (i__1 = ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1], abs(i__1));
		if (idd == 63 || idd == 64 || idd == 48 || idd == 49 || idd == 121 ||
			idd == 100 || idd == 88 || idd == 66 || idd == 90 || idd == 70)
		{
			signd = x1535 + sigk + s2d;
			if (x1 <= (*signn + sdprod) / *sig)
			{
				goto L108;
			}
			if (x1 > (*signn + signd + sdprod) / *sig)
			{
				return 0;
			}

			if (sigk / signd > ranart_(&rndf77_1.nseed))
			{
				goto L306;
			}
			if (s2d / (x1535 + s2d) > ranart_(&rndf77_1.nseed))
			{
				goto L1012;
			}
			if (idd == 63)
			{
				n12 = 17;
			}
			if (idd == 64)
			{
				n12 = 20;
			}
			if (idd == 48)
			{
				n12 = 23;
			}
			if (idd == 49)
			{
				n12 = 24;
			}
			if (idd == 121)
			{
				n12 = 25;
			}
			if (idd == 100)
			{
				n12 = 26;
			}
			if (idd == 88)
			{
				n12 = 29;
			}
			if (idd == 66)
			{
				n12 = 31;
			}
			if (idd == 90)
			{
				n12 = 32;
			}
			if (idd == 70)
			{
				n12 = 35;
			}
			goto L1011;
		}
		if (idd == 110 || idd == 77 || idd == 80)
		{
			if (x1 <= (*signn + sdprod) / *sig)
			{
				goto L108;
			}
			if (x1 > (*signn + x1535 + sigk + s2d + sdprod) / *sig)
			{
				return 0;
			}

			if (sigk / (x1535 + sigk + s2d) > ranart_(&rndf77_1.nseed))
			{
				goto L306;
			}
			if (s2d / (x1535 + s2d) > ranart_(&rndf77_1.nseed))
			{
				goto L1012;
			}
			if (idd == 77)
			{
				n12 = 30;
			}
			if (idd == 77 && ranart_(&rndf77_1.nseed) <= .5f)
			{
				n12 = 36;
			}
			if (idd == 80)
			{
				n12 = 34;
			}
			if (idd == 80 && ranart_(&rndf77_1.nseed) <= .5f)
			{
				n12 = 35;
			}
			if (idd == 110)
			{
				n12 = 27;
			}
			if (idd == 110 && ranart_(&rndf77_1.nseed) <= .5f)
			{
				n12 = 28;
			}
			goto L1011;
		}
		if (idd == 54 || idd == 56)
		{
			sig2 = sigma_(srt, &c__2, &c__0, &c__1) * .75f;
			signd = (sig2 + x1535) * 2.f + sigk + s2d;
			if (x1 <= (*signn + sdprod) / *sig)
			{
				goto L108;
			}
			if (x1 > (*signn + signd + sdprod) / *sig)
			{
				return 0;
			}

			if (sigk / signd > ranart_(&rndf77_1.nseed))
			{
				goto L306;
			}
			if (s2d / ((sig2 + x1535) * 2.f + s2d) > ranart_(&rndf77_1.nseed))
			{
				goto L1012;
			}
			if (ranart_(&rndf77_1.nseed) < x1535 / (sig2 + x1535))
			{
				if (idd == 54)
				{
					n12 = 18;
				}
				if (idd == 54 && ranart_(&rndf77_1.nseed) <= .5f)
				{
					n12 = 19;
				}
				if (idd == 56)
				{
					n12 = 21;
				}
				if (idd == 56 && ranart_(&rndf77_1.nseed) <= .5f)
				{
					n12 = 22;
				}
			}
			else
			{
				if (idd == 54)
				{
					n12 = 13;
				}
				if (idd == 54 && ranart_(&rndf77_1.nseed) <= .5f)
				{
					n12 = 14;
				}
				if (idd == 56)
				{
					n12 = 15;
				}
				if (idd == 56 && ranart_(&rndf77_1.nseed) <= .5f)
				{
					n12 = 16;
				}
			}
		}
	L1011:
		*iblock = 5;
		dmax__ = *srt - .9383f - .005f;
		dmin__ = 1.078f;
		if (n12 >= 13 && n12 <= 16)
		{
			if (dmax__ < 1.44f)
			{
				fm = fns_(&dmax__, srt, &c_b172);
			}
			else
			{
				xdmass = 1.44f;
				fm = fns_(&xdmass, srt, &c_b173);
			}
			if (fm == 0.f)
			{
				fm = 1e-9f;
			}
			ntry2 = 0;
		L11:
			dm = ranart_(&rndf77_1.nseed) * (dmax__ - dmin__) + dmin__;
			++ntry2;
			if (ranart_(&rndf77_1.nseed) > fns_(&dm, srt, &c_b173) / fm &&
				ntry2 <= 10)
			{
				goto L11;
			}
			if (dm > 2.14f)
			{
				goto L11;
			}
			goto L13;
		}
		if (n12 >= 17 && n12 <= 36)
		{
			if (dmax__ < 1.535f)
			{
				fm = fd5_(&dmax__, srt, &c_b172);
			}
			else
			{
				xdmass = 1.535f;
				fm = fd5_(&xdmass, srt, &c_b173);
			}
			if (fm == 0.f)
			{
				fm = 1e-9f;
			}
			ntry1 = 0;
		L12:
			dm = ranart_(&rndf77_1.nseed) * (dmax__ - dmin__) + dmin__;
			++ntry1;
			if (ranart_(&rndf77_1.nseed) > fd5_(&dm, srt, &c_b173) / fm &&
				ntry1 <= 10)
			{
				goto L12;
			}
			if (dm > 1.84f)
			{
				goto L12;
			}
		}
	L13:
		if (n12 == 13)
		{
			if (ranart_(&rndf77_1.nseed) <= .5f)
			{
				ee_1.lb[*i2 - 1] = 11;
				cc_1.e[*i2 - 1] = dm;
				ee_1.lb[*i1 - 1] = 2;
				cc_1.e[*i1 - 1] = .939457f;
			}
			else
			{
				ee_1.lb[*i1 - 1] = 11;
				cc_1.e[*i1 - 1] = dm;
				ee_1.lb[*i2 - 1] = 2;
				cc_1.e[*i2 - 1] = .939457f;
			}
			goto L200;
		}
		if (n12 == 14)
		{
			if (ranart_(&rndf77_1.nseed) <= .5f)
			{
				ee_1.lb[*i2 - 1] = 10;
				cc_1.e[*i2 - 1] = dm;
				ee_1.lb[*i1 - 1] = 1;
				cc_1.e[*i1 - 1] = .93828f;
			}
			else
			{
				ee_1.lb[*i1 - 1] = 10;
				cc_1.e[*i1 - 1] = dm;
				ee_1.lb[*i2 - 1] = 1;
				cc_1.e[*i2 - 1] = .93828f;
			}
			goto L200;
		}
		if (n12 == 15)
		{
			if (ranart_(&rndf77_1.nseed) <= .5f)
			{
				ee_1.lb[*i2 - 1] = 11;
				cc_1.e[*i2 - 1] = dm;
				ee_1.lb[*i1 - 1] = 2;
				cc_1.e[*i1 - 1] = .939457f;
			}
			else
			{
				ee_1.lb[*i1 - 1] = 11;
				cc_1.e[*i1 - 1] = dm;
				ee_1.lb[*i2 - 1] = 2;
				cc_1.e[*i2 - 1] = .939457f;
			}
			goto L200;
		}
		if (n12 == 16)
		{
			if (ranart_(&rndf77_1.nseed) <= .5f)
			{
				ee_1.lb[*i2 - 1] = 10;
				cc_1.e[*i2 - 1] = dm;
				ee_1.lb[*i1 - 1] = 1;
				cc_1.e[*i1 - 1] = .93828f;
			}
			else
			{
				ee_1.lb[*i1 - 1] = 10;
				cc_1.e[*i1 - 1] = dm;
				ee_1.lb[*i2 - 1] = 1;
				cc_1.e[*i2 - 1] = .93828f;
			}
			goto L200;
		}
		if (n12 == 17)
		{
			ee_1.lb[*i2 - 1] = 13;
			cc_1.e[*i2 - 1] = dm;
			ee_1.lb[*i1 - 1] = 1;
			cc_1.e[*i1 - 1] = .93828f;
			goto L200;
		}
		if (n12 == 18)
		{
			if (ranart_(&rndf77_1.nseed) <= .5f)
			{
				ee_1.lb[*i2 - 1] = 12;
				cc_1.e[*i2 - 1] = dm;
				ee_1.lb[*i1 - 1] = 1;
				cc_1.e[*i1 - 1] = .93828f;
			}
			else
			{
				ee_1.lb[*i1 - 1] = 12;
				cc_1.e[*i1 - 1] = dm;
				ee_1.lb[*i2 - 1] = 1;
				cc_1.e[*i2 - 1] = .93828f;
			}
			goto L200;
		}
		if (n12 == 19)
		{
			if (ranart_(&rndf77_1.nseed) <= .5f)
			{
				ee_1.lb[*i2 - 1] = 13;
				cc_1.e[*i2 - 1] = dm;
				ee_1.lb[*i1 - 1] = 2;
				cc_1.e[*i1 - 1] = .939457f;
			}
			else
			{
				ee_1.lb[*i1 - 1] = 13;
				cc_1.e[*i1 - 1] = dm;
				ee_1.lb[*i2 - 1] = 2;
				cc_1.e[*i2 - 1] = .939457f;
			}
			goto L200;
		}
		if (n12 == 20)
		{
			if (ranart_(&rndf77_1.nseed) <= .5f)
			{
				ee_1.lb[*i2 - 1] = 13;
				cc_1.e[*i2 - 1] = dm;
				ee_1.lb[*i1 - 1] = 1;
				cc_1.e[*i1 - 1] = .93828f;
			}
			else
			{
				ee_1.lb[*i1 - 1] = 13;
				cc_1.e[*i1 - 1] = dm;
				ee_1.lb[*i2 - 1] = 1;
				cc_1.e[*i2 - 1] = .93828f;
			}
			goto L200;
		}
		if (n12 == 21)
		{
			if (ranart_(&rndf77_1.nseed) <= .5f)
			{
				ee_1.lb[*i2 - 1] = 13;
				cc_1.e[*i2 - 1] = dm;
				ee_1.lb[*i1 - 1] = 2;
				cc_1.e[*i1 - 1] = .939457f;
			}
			else
			{
				ee_1.lb[*i1 - 1] = 13;
				cc_1.e[*i1 - 1] = dm;
				ee_1.lb[*i2 - 1] = 2;
				cc_1.e[*i2 - 1] = .939457f;
			}
			goto L200;
		}
		if (n12 == 22)
		{
			if (ranart_(&rndf77_1.nseed) <= .5f)
			{
				ee_1.lb[*i2 - 1] = 12;
				cc_1.e[*i2 - 1] = dm;
				ee_1.lb[*i1 - 1] = 1;
				cc_1.e[*i1 - 1] = .93828f;
			}
			else
			{
				ee_1.lb[*i1 - 1] = 12;
				cc_1.e[*i1 - 1] = dm;
				ee_1.lb[*i2 - 1] = 1;
				cc_1.e[*i2 - 1] = .93828f;
			}
			goto L200;
		}
		if (n12 == 23)
		{
			if (ranart_(&rndf77_1.nseed) <= .5f)
			{
				ee_1.lb[*i2 - 1] = 12;
				cc_1.e[*i2 - 1] = dm;
				ee_1.lb[*i1 - 1] = 2;
				cc_1.e[*i1 - 1] = .939457f;
			}
			else
			{
				ee_1.lb[*i1 - 1] = 12;
				cc_1.e[*i1 - 1] = dm;
				ee_1.lb[*i2 - 1] = 2;
				cc_1.e[*i2 - 1] = .939457f;
			}
			goto L200;
		}
		if (n12 == 24)
		{
			ee_1.lb[*i2 - 1] = 12;
			cc_1.e[*i2 - 1] = dm;
			ee_1.lb[*i1 - 1] = 2;
			cc_1.e[*i1 - 1] = .939457f;
			goto L200;
		}
		if (n12 == 25)
		{
			ee_1.lb[*i2 - 1] = 12;
			cc_1.e[*i2 - 1] = dm;
			ee_1.lb[*i1 - 1] = 1;
			cc_1.e[*i1 - 1] = .93828f;
			goto L200;
		}
		if (n12 == 26)
		{
			ee_1.lb[*i2 - 1] = 12;
			cc_1.e[*i2 - 1] = dm;
			ee_1.lb[*i1 - 1] = 2;
			cc_1.e[*i1 - 1] = .939457f;
			goto L200;
		}
		if (n12 == 27)
		{
			if (ranart_(&rndf77_1.nseed) <= .5f)
			{
				ee_1.lb[*i2 - 1] = 13;
				cc_1.e[*i2 - 1] = dm;
				ee_1.lb[*i1 - 1] = 2;
				cc_1.e[*i1 - 1] = .939457f;
			}
			else
			{
				ee_1.lb[*i1 - 1] = 13;
				cc_1.e[*i1 - 1] = dm;
				ee_1.lb[*i2 - 1] = 2;
				cc_1.e[*i2 - 1] = .939457f;
			}
			goto L200;
		}
		if (n12 == 28)
		{
			if (ranart_(&rndf77_1.nseed) <= .5f)
			{
				ee_1.lb[*i2 - 1] = 12;
				cc_1.e[*i2 - 1] = dm;
				ee_1.lb[*i1 - 1] = 1;
				cc_1.e[*i1 - 1] = .93828f;
			}
			else
			{
				ee_1.lb[*i1 - 1] = 12;
				cc_1.e[*i1 - 1] = dm;
				ee_1.lb[*i2 - 1] = 1;
				cc_1.e[*i2 - 1] = .93828f;
			}
			goto L200;
		}
		if (n12 == 27)
		{
			if (ranart_(&rndf77_1.nseed) <= .5f)
			{
				ee_1.lb[*i2 - 1] = 13;
				cc_1.e[*i2 - 1] = dm;
				ee_1.lb[*i1 - 1] = 2;
				cc_1.e[*i1 - 1] = .939457f;
			}
			else
			{
				ee_1.lb[*i1 - 1] = 13;
				cc_1.e[*i1 - 1] = dm;
				ee_1.lb[*i2 - 1] = 2;
				cc_1.e[*i2 - 1] = .939457f;
			}
			goto L200;
		}
		if (n12 == 29)
		{
			if (ranart_(&rndf77_1.nseed) <= .5f)
			{
				ee_1.lb[*i2 - 1] = 13;
				cc_1.e[*i2 - 1] = dm;
				ee_1.lb[*i1 - 1] = 1;
				cc_1.e[*i1 - 1] = .93828f;
			}
			else
			{
				ee_1.lb[*i1 - 1] = 13;
				cc_1.e[*i1 - 1] = dm;
				ee_1.lb[*i2 - 1] = 1;
				cc_1.e[*i2 - 1] = .93828f;
			}
			goto L200;
		}
		if (n12 == 30)
		{
			if (ranart_(&rndf77_1.nseed) <= .5f)
			{
				ee_1.lb[*i2 - 1] = 13;
				cc_1.e[*i2 - 1] = dm;
				ee_1.lb[*i1 - 1] = 2;
				cc_1.e[*i1 - 1] = .939457f;
			}
			else
			{
				ee_1.lb[*i1 - 1] = 13;
				cc_1.e[*i1 - 1] = dm;
				ee_1.lb[*i2 - 1] = 2;
				cc_1.e[*i2 - 1] = .939457f;
			}
			goto L200;
		}
		if (n12 == 31)
		{
			if (ranart_(&rndf77_1.nseed) <= .5f)
			{
				ee_1.lb[*i2 - 1] = 12;
				cc_1.e[*i2 - 1] = dm;
				ee_1.lb[*i1 - 1] = 2;
				cc_1.e[*i1 - 1] = .939457f;
			}
			else
			{
				ee_1.lb[*i1 - 1] = 12;
				cc_1.e[*i1 - 1] = dm;
				ee_1.lb[*i2 - 1] = 2;
				cc_1.e[*i2 - 1] = .939457f;
			}
			goto L200;
		}
		if (n12 == 32)
		{
			if (ranart_(&rndf77_1.nseed) <= .5f)
			{
				ee_1.lb[*i2 - 1] = 13;
				cc_1.e[*i2 - 1] = dm;
				ee_1.lb[*i1 - 1] = 1;
				cc_1.e[*i1 - 1] = .93828f;
			}
			else
			{
				ee_1.lb[*i1 - 1] = 13;
				cc_1.e[*i1 - 1] = dm;
				ee_1.lb[*i2 - 1] = 1;
				cc_1.e[*i2 - 1] = .93828f;
			}
			goto L200;
		}
		if (n12 == 33)
		{
			if (ranart_(&rndf77_1.nseed) <= .5f)
			{
				ee_1.lb[*i2 - 1] = 13;
				cc_1.e[*i2 - 1] = dm;
				ee_1.lb[*i1 - 1] = 2;
				cc_1.e[*i1 - 1] = .939457f;
			}
			else
			{
				ee_1.lb[*i1 - 1] = 13;
				cc_1.e[*i1 - 1] = dm;
				ee_1.lb[*i2 - 1] = 2;
				cc_1.e[*i2 - 1] = .939457f;
			}
			goto L200;
		}
		if (n12 == 34)
		{
			if (ranart_(&rndf77_1.nseed) <= .5f)
			{
				ee_1.lb[*i2 - 1] = 12;
				cc_1.e[*i2 - 1] = dm;
				ee_1.lb[*i1 - 1] = 1;
				cc_1.e[*i1 - 1] = .93828f;
			}
			else
			{
				ee_1.lb[*i1 - 1] = 12;
				cc_1.e[*i1 - 1] = dm;
				ee_1.lb[*i2 - 1] = 1;
				cc_1.e[*i2 - 1] = .93828f;
			}
			goto L200;
		}
		if (n12 == 35)
		{
			if (ranart_(&rndf77_1.nseed) <= .5f)
			{
				ee_1.lb[*i2 - 1] = 12;
				cc_1.e[*i2 - 1] = dm;
				ee_1.lb[*i1 - 1] = 2;
				cc_1.e[*i1 - 1] = .939457f;
			}
			else
			{
				ee_1.lb[*i1 - 1] = 12;
				cc_1.e[*i1 - 1] = dm;
				ee_1.lb[*i2 - 1] = 2;
				cc_1.e[*i2 - 1] = .939457f;
			}
			goto L200;
		}
		if (n12 == 36)
		{
			if (ranart_(&rndf77_1.nseed) <= .5f)
			{
				ee_1.lb[*i2 - 1] = 12;
				cc_1.e[*i2 - 1] = dm;
				ee_1.lb[*i1 - 1] = 1;
				cc_1.e[*i1 - 1] = .93828f;
			}
			else
			{
				ee_1.lb[*i1 - 1] = 12;
				cc_1.e[*i1 - 1] = dm;
				ee_1.lb[*i2 - 1] = 1;
				cc_1.e[*i2 - 1] = .93828f;
			}
			goto L200;
		}
	L1012:
		*iblock = 55;
		leadng_1.lb1 = ee_1.lb[*i1 - 1];
		dpi_1.lb2 = ee_1.lb[*i2 - 1];
		ich = (i__1 = leadng_1.lb1 * dpi_1.lb2, abs(i__1));
		if (ich == 54)
		{
			if (ranart_(&rndf77_1.nseed) <= .5f)
			{
				ee_1.lb[*i2 - 1] = 1;
				cc_1.e[*i2 - 1] = .93828f;
				ee_1.lb[*i1 - 1] = 2;
				cc_1.e[*i1 - 1] = .939457f;
			}
			else
			{
				ee_1.lb[*i1 - 1] = 1;
				cc_1.e[*i1 - 1] = .93828f;
				ee_1.lb[*i2 - 1] = 2;
				cc_1.e[*i2 - 1] = .939457f;
			}
			goto L200;
		}
		if (ich == 56)
		{
			if (ranart_(&rndf77_1.nseed) <= .5f)
			{
				ee_1.lb[*i2 - 1] = 1;
				cc_1.e[*i2 - 1] = .93828f;
				ee_1.lb[*i1 - 1] = 2;
				cc_1.e[*i1 - 1] = .939457f;
			}
			else
			{
				ee_1.lb[*i1 - 1] = 1;
				cc_1.e[*i1 - 1] = .93828f;
				ee_1.lb[*i2 - 1] = 2;
				cc_1.e[*i2 - 1] = .939457f;
			}
			goto L200;
		}
		if (ich == 63)
		{
			ee_1.lb[*i2 - 1] = 1;
			cc_1.e[*i2 - 1] = .93828f;
			ee_1.lb[*i1 - 1] = 1;
			cc_1.e[*i1 - 1] = .93828f;
			goto L200;
		}
		if (ich == 64)
		{
			ee_1.lb[*i2 - 1] = 1;
			cc_1.e[*i2 - 1] = .93828f;
			ee_1.lb[*i1 - 1] = 1;
			cc_1.e[*i1 - 1] = .93828f;
			goto L200;
		}
		if (ich == 48)
		{
			ee_1.lb[*i2 - 1] = 2;
			cc_1.e[*i2 - 1] = .939457f;
			ee_1.lb[*i1 - 1] = 2;
			cc_1.e[*i1 - 1] = .939457f;
			goto L200;
		}
		if (ich == 36)
		{
			ee_1.lb[*i2 - 1] = 2;
			cc_1.e[*i2 - 1] = .939457f;
			ee_1.lb[*i1 - 1] = 2;
			cc_1.e[*i1 - 1] = .939457f;
			goto L200;
		}
		if (ich == 121 || ich == 169 || ich == 143)
		{
			ee_1.lb[*i2 - 1] = 1;
			cc_1.e[*i2 - 1] = .93828f;
			ee_1.lb[*i1 - 1] = 1;
			cc_1.e[*i1 - 1] = .93828f;
			goto L200;
		}
		if (ich == 100 || ich == 144 || ich == 120)
		{
			ee_1.lb[*i2 - 1] = 2;
			cc_1.e[*i2 - 1] = .939457f;
			ee_1.lb[*i1 - 1] = 2;
			cc_1.e[*i1 - 1] = .939457f;
			goto L200;
		}
		if (ich == 110 || ich == 156 || ich == 130 || ich == 132)
		{
			if (ranart_(&rndf77_1.nseed) <= .5f)
			{
				ee_1.lb[*i2 - 1] = 1;
				cc_1.e[*i2 - 1] = .93828f;
				ee_1.lb[*i1 - 1] = 2;
				cc_1.e[*i1 - 1] = .939457f;
			}
			else
			{
				ee_1.lb[*i1 - 1] = 1;
				cc_1.e[*i1 - 1] = .93828f;
				ee_1.lb[*i2 - 1] = 2;
				cc_1.e[*i2 - 1] = .939457f;
			}
			goto L200;
		}
		if (ich == 88 || ich == 104)
		{
			ee_1.lb[*i2 - 1] = 1;
			cc_1.e[*i2 - 1] = .93828f;
			ee_1.lb[*i1 - 1] = 1;
			cc_1.e[*i1 - 1] = .93828f;
			goto L200;
		}
		if (ich == 77 || ich == 91)
		{
			if (ranart_(&rndf77_1.nseed) <= .5f)
			{
				ee_1.lb[*i2 - 1] = 1;
				cc_1.e[*i2 - 1] = .93828f;
				ee_1.lb[*i1 - 1] = 2;
				cc_1.e[*i1 - 1] = .939457f;
			}
			else
			{
				ee_1.lb[*i1 - 1] = 1;
				cc_1.e[*i1 - 1] = .93828f;
				ee_1.lb[*i2 - 1] = 2;
				cc_1.e[*i2 - 1] = .939457f;
			}
			goto L200;
		}
		if (ich == 66 || ich == 78)
		{
			ee_1.lb[*i2 - 1] = 2;
			cc_1.e[*i2 - 1] = .939457f;
			ee_1.lb[*i1 - 1] = 2;
			cc_1.e[*i1 - 1] = .939457f;
			goto L200;
		}
		if (ich == 90 || ich == 108)
		{
			ee_1.lb[*i2 - 1] = 1;
			cc_1.e[*i2 - 1] = .93828f;
			ee_1.lb[*i1 - 1] = 1;
			cc_1.e[*i1 - 1] = .93828f;
			goto L200;
		}
		if (ich == 70 || ich == 84)
		{
			ee_1.lb[*i2 - 1] = 2;
			cc_1.e[*i2 - 1] = .939457f;
			ee_1.lb[*i1 - 1] = 2;
			cc_1.e[*i1 - 1] = .939457f;
			goto L200;
		}
		if (ich == 80 || ich == 96)
		{
			if (ranart_(&rndf77_1.nseed) <= .5f)
			{
				ee_1.lb[*i2 - 1] = 2;
				cc_1.e[*i2 - 1] = .939457f;
				ee_1.lb[*i1 - 1] = 1;
				cc_1.e[*i1 - 1] = .93828f;
			}
			else
			{
				ee_1.lb[*i1 - 1] = 2;
				cc_1.e[*i1 - 1] = .939457f;
				ee_1.lb[*i2 - 1] = 1;
				cc_1.e[*i2 - 1] = .93828f;
			}
			goto L200;
		}
		ee_1.lb[*i1 - 1] = 1;
		cc_1.e[*i1 - 1] = .93828f;
		ee_1.lb[*i2 - 1] = 2;
		cc_1.e[*i2 - 1] = .939457f;
	L200:
		leadng_1.em1 = cc_1.e[*i1 - 1];
		dpi_1.em2 = cc_1.e[*i2 - 1];
		r__2 = *srt;
		r__3 = leadng_1.em1;
		r__4 = dpi_1.em2;
		r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
		r__5 = leadng_1.em1 * dpi_1.em2;
		pr2 = r__1 * r__1 - r__5 * r__5 * 4.f;
		if (pr2 <= 0.f)
		{
			pr2 = 1e-9f;
		}
		pr = sqrt(pr2) / (*srt * 2.f);
		if (*srt <= 2.14f)
		{
			c1 = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
		}
		if (*srt > 2.14f && *srt <= 2.4f)
		{
			c1 = ang_(srt, &input1_1.iseed);
		}
		if (*srt > 2.4f)
		{
			xptr = pr * .33f;
			cc1 = ptr_(&xptr, &input1_1.iseed);
			r__1 = pr;
			r__2 = cc1;
			c1 = sqrt(r__1 * r__1 - r__2 * r__2) / pr;
		}
		t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
		if (ianti == 1 && ee_1.lb[*i1 - 1] >= 1 && ee_1.lb[*i2 - 1] >= 1)
		{
			ee_1.lb[*i1 - 1] = -ee_1.lb[*i1 - 1];
			ee_1.lb[*i2 - 1] = -ee_1.lb[*i2 - 1];
		}
	}
L107:
	r__1 = c1;
	s1 = sqrt(1.f - r__1 * r__1);
	r__1 = c2;
	s2 = sqrt(1.f - r__1 * r__1);
	ct1 = cos(t1);
	st1 = sin(t1);
	ct2 = cos(t2);
	st2 = sin(t2);
	*pz = pr * (c1 * c2 - s1 * s2 * ct1);
	ss = c2 * s1 * ct1 + s2 * c1;
	*px = pr * (ss * ct2 - s1 * st1 * st2);
	*py = pr * (ss * st2 + s1 * st1 * ct2);
	return 0;
L306:
	if (xsk5 / sigk > ranart_(&rndf77_1.nseed))
	{
		leadng_1.pz1 = bb_1.p[*i1 * 3 - 1];
		pz2 = bb_1.p[*i2 * 3 - 1];
		ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 2) + 1;
		ee_1.lb[*i2 - 1] = (int)(ranart_(&rndf77_1.nseed) * 2) + 1;
		++nn_1.nnn;
		pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 29;
		pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = 1.02f;
		*iblock = 222;
		goto L208;
	}
	*iblock = 10;
	if (ianti == 1)
	{
		*iblock = -10;
	}
	leadng_1.pz1 = bb_1.p[*i1 * 3 - 1];
	pz2 = bb_1.p[*i2 * 3 - 1];
	++nn_1.nnn;
	pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 23;
	pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = .498f;
	if (*srt <= 2.63f)
	{
		ic = 1;
		ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 2) + 1;
		ee_1.lb[*i2 - 1] = 14;
		goto L208;
	}
	if (*srt <= 2.74f && *srt > 2.63f)
	{
		if (xsk1 / (xsk1 + xsk2) > ranart_(&rndf77_1.nseed))
		{
			ic = 1;
			ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 2) + 1;
			ee_1.lb[*i2 - 1] = 14;
		}
		else
		{
			ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 2) + 1;
			ee_1.lb[*i2 - 1] = (int)(ranart_(&rndf77_1.nseed) * 3) + 15;
			ic = 2;
		}
		goto L208;
	}
	if (*srt <= 2.77f && *srt > 2.74f)
	{
		if (xsk1 / (xsk1 + xsk2 + xsk3) > ranart_(&rndf77_1.nseed))
		{
			ic = 1;
			ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 2) + 1;
			ee_1.lb[*i2 - 1] = 14;
			goto L208;
		}
		else
		{
			if (xsk2 / (xsk2 + xsk3) > ranart_(&rndf77_1.nseed))
			{
				ic = 2;
				ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 2) +
								   1;
				ee_1.lb[*i2 - 1] = (int)(ranart_(&rndf77_1.nseed) * 3) +
								   15;
			}
			else
			{
				ic = 3;
				ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 4) +
								   6;
				ee_1.lb[*i2 - 1] = 14;
			}
			goto L208;
		}
	}
	if (*srt > 2.77f)
	{
		if (xsk1 / (xsk1 + xsk2 + xsk3 + xsk4) > ranart_(&rndf77_1.nseed))
		{
			ic = 1;
			ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 2) + 1;
			ee_1.lb[*i2 - 1] = 14;
			goto L208;
		}
		else
		{
			if (xsk3 / (xsk2 + xsk3 + xsk4) > ranart_(&rndf77_1.nseed))
			{
				ic = 3;
				ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 4) +
								   6;
				ee_1.lb[*i2 - 1] = 14;
				goto L208;
			}
			else
			{
				if (xsk2 / (xsk2 + xsk4) > ranart_(&rndf77_1.nseed))
				{
					ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) *
												 2) +
									   1;
					ee_1.lb[*i2 - 1] = (int)(ranart_(&rndf77_1.nseed) *
												 3) +
									   15;
					ic = 2;
				}
				else
				{
					ic = 4;
					ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) *
												 4) +
									   6;
					ee_1.lb[*i2 - 1] = (int)(ranart_(&rndf77_1.nseed) *
												 3) +
									   15;
				}
				goto L208;
			}
		}
	}
L208:
	if (ianti == 1 && ee_1.lb[*i1 - 1] >= 1 && ee_1.lb[*i2 - 1] >= 1)
	{
		ee_1.lb[*i1 - 1] = -ee_1.lb[*i1 - 1];
		ee_1.lb[*i2 - 1] = -ee_1.lb[*i2 - 1];
		if (pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] == 23)
		{
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 21;
		}
	}
	lbi1 = ee_1.lb[*i1 - 1];
	lbi2 = ee_1.lb[*i2 - 1];
	ntry1 = 0;
L129:
	bbkaon_(&ic, srt, &px3, &py3, &pz3, &dm3, &px4, &py4, &pz4, &dm4, &ppx, &ppy, &ppz, &icou1);
	++ntry1;
	if (icou1 < 0 && ntry1 <= 20)
	{
		goto L129;
	}
	rotate_(px, py, pz, &px3, &py3, &pz3);
	rotate_(px, py, pz, &px4, &py4, &pz4);
	rotate_(px, py, pz, &ppx, &ppy, &ppz);
	r__1 = dm3;
	r__2 = px3;
	r__3 = py3;
	r__4 = pz3;
	e1cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	p1beta = px3 * bg_1.betax + py3 * bg_1.betay + pz3 * bg_1.betaz;
	transf = bg_1.gamma * (bg_1.gamma * p1beta / (bg_1.gamma + 1) + e1cm);
	pt1i1 = bg_1.betax * transf + px3;
	pt2i1 = bg_1.betay * transf + py3;
	pt3i1 = bg_1.betaz * transf + pz3;
	eti1 = dm3;
	r__1 = dm4;
	r__2 = px4;
	r__3 = py4;
	r__4 = pz4;
	e2cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	p2beta = px4 * bg_1.betax + py4 * bg_1.betay + pz4 * bg_1.betaz;
	transf = bg_1.gamma * (bg_1.gamma * p2beta / (bg_1.gamma + 1.f) + e2cm);
	pt1i2 = bg_1.betax * transf + px4;
	pt2i2 = bg_1.betay * transf + py4;
	pt3i2 = bg_1.betaz * transf + pz4;
	eti2 = dm4;
	r__1 = ppx;
	r__2 = ppy;
	r__3 = ppz;
	epcm = sqrt(r__1 * r__1 + .248004f + r__2 * r__2 + r__3 * r__3);
	ppbeta = ppx * bg_1.betax + ppy * bg_1.betay + ppz * bg_1.betaz;
	transf = bg_1.gamma * (bg_1.gamma * ppbeta / (bg_1.gamma + 1.f) + epcm);
	pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = bg_1.betax *
															   transf +
														   ppx;
	pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = bg_1.betay *
															   transf +
														   ppy;
	pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = bg_1.betaz *
															   transf +
														   ppz;
	dpert_1.dppion[nn_1.nnn + *irun * 150001 - 150002] = dpert_1.dpertp[*i1 -
																		1] *
														 dpert_1.dpertp[*i2 - 1];
	pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = aa_1.r__[*i1 * 3 -
																	3];
	pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = aa_1.r__[*i1 * 3 -
																	2];
	pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = aa_1.r__[*i1 * 3 -
																	1];

	bb_1.p[*i1 * 3 - 3] = pt1i1;
	bb_1.p[*i1 * 3 - 2] = pt2i1;
	bb_1.p[*i1 * 3 - 1] = pt3i1;
	cc_1.e[*i1 - 1] = eti1;
	ee_1.lb[*i1 - 1] = lbi1;
	bb_1.p[*i2 * 3 - 3] = pt1i2;
	bb_1.p[*i2 * 3 - 2] = pt2i2;
	bb_1.p[*i2 * 3 - 1] = pt3i2;
	cc_1.e[*i2 - 1] = eti2;
	ee_1.lb[*i2 - 1] = lbi2;
	leadng_1.px1 = bb_1.p[*i1 * 3 - 3];
	leadng_1.py1 = bb_1.p[*i1 * 3 - 2];
	leadng_1.pz1 = bb_1.p[*i1 * 3 - 1];
	leadng_1.em1 = cc_1.e[*i1 - 1];
	ee_1.id[*i1 - 1] = 2;
	ee_1.id[*i2 - 1] = 2;
	id1 = ee_1.id[*i1 - 1];
	leadng_1.lb1 = ee_1.lb[*i1 - 1];
	dpi_1.lb2 = ee_1.lb[*i2 - 1];
	am1 = leadng_1.em1;
	am2 = dpi_1.em2;
	r__1 = leadng_1.em1;
	r__2 = leadng_1.px1;
	r__3 = leadng_1.py1;
	r__4 = leadng_1.pz1;
	leadng_1.e1 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	return 0;
L108:
	if (para8_1.idpert == 1 && *ipert1 == 1 && para8_1.npertd >= 1)
	{
		ndloop = para8_1.npertd;
	}
	else if (para8_1.idpert == 2 && para8_1.npertd >= 1)
	{
		ndloop = para8_1.npertd + 1;
	}
	else
	{
		ndloop = 1;
	}

	dprob1 = sdprod / *sig / (float)para8_1.npertd;
	i__1 = ndloop;
	for (idloop = 1; idloop <= i__1; ++idloop)
	{
		bbdangle_(&pxd, &pyd, &pzd, nt, ipert1, &ianti, &idloop, &pfinal, &dprob1, &lbm);
		rotate_(px, py, pz, &pxd, &pyd, &pzd);
		xmass = 1.8756f;
		r__1 = xmass;
		r__2 = pxd;
		r__3 = pyd;
		r__4 = pzd;
		e1dcm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
		p1dbeta = pxd * bg_1.betax + pyd * bg_1.betay + pzd * bg_1.betaz;
		transf = bg_1.gamma * (bg_1.gamma * p1dbeta / (bg_1.gamma + 1.f) +
							   e1dcm);
		pxi1 = bg_1.betax * transf + pxd;
		pyi1 = bg_1.betay * transf + pyd;
		pzi1 = bg_1.betaz * transf + pzd;
		if (ianti == 0)
		{
			lbd = 42;
		}
		else
		{
			lbd = -42;
		}
		if (para8_1.idpert == 1 && *ipert1 == 1 && para8_1.npertd >= 1)
		{
			++nn_1.nnn;
			pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = pxi1;
			pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = pyi1;
			pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = pzi1;
			pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = 1.8756f;
			pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = lbd;
			pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = aa_1.r__[*i1 * 3 - 3];
			pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = aa_1.r__[*i1 * 3 - 2];
			pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = aa_1.r__[*i1 * 3 - 1];
			dpert_1.dppion[nn_1.nnn + *irun * 150001 - 150002] = sdprod / *sig / (float)para8_1.npertd;
		}
		else if (para8_1.idpert == 2 && idloop <= para8_1.npertd)
		{
			ppd[idloop * 3 - 3] = pxi1;
			ppd[idloop * 3 - 2] = pyi1;
			ppd[idloop * 3 - 1] = pzi1;
			lbpd[idloop - 1] = lbd;
		}
		else
		{
			cc_1.e[*i1 - 1] = xmm;
			r__1 = xmm;
			r__2 = pxd;
			r__3 = pyd;
			r__4 = pzd;
			e2picm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
			p2pibeta = -pxd * bg_1.betax - pyd * bg_1.betay - pzd * bg_1.betaz;
			transf = bg_1.gamma * (bg_1.gamma * p2pibeta / (bg_1.gamma + 1.f) + e2picm);
			pxi2 = bg_1.betax * transf - pxd;
			pyi2 = bg_1.betay * transf - pyd;
			pzi2 = bg_1.betaz * transf - pzd;
			bb_1.p[*i1 * 3 - 3] = pxi2;
			bb_1.p[*i1 * 3 - 2] = pyi2;
			bb_1.p[*i1 * 3 - 1] = pzi2;
			ee_1.lb[*i1 - 1] = lbm;
			leadng_1.px1 = bb_1.p[*i1 * 3 - 3];
			leadng_1.py1 = bb_1.p[*i1 * 3 - 2];
			leadng_1.pz1 = bb_1.p[*i1 * 3 - 1];
			leadng_1.em1 = cc_1.e[*i1 - 1];
			ee_1.id[*i1 - 1] = 2;
			id1 = ee_1.id[*i1 - 1];
			r__1 = leadng_1.em1;
			r__2 = leadng_1.px1;
			r__3 = leadng_1.py1;
			r__4 = leadng_1.pz1;
			leadng_1.e1 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
			leadng_1.lb1 = ee_1.lb[*i1 - 1];
			bb_1.p[*i2 * 3 - 3] = pxi1;
			bb_1.p[*i2 * 3 - 2] = pyi1;
			bb_1.p[*i2 * 3 - 1] = pzi1;
			ee_1.lb[*i2 - 1] = lbd;
			dpi_1.lb2 = ee_1.lb[*i2 - 1];
			cc_1.e[*i2 - 1] = 1.8756f;
			eti2 = cc_1.e[*i2 - 1];
			ee_1.id[*i2 - 1] = 2;
			if (para8_1.idpert == 2 && idloop == ndloop)
			{
				i__2 = para8_1.npertd;
				for (ipertd = 1; ipertd <= i__2; ++ipertd)
				{
					++nn_1.nnn;
					pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450006] =
						ppd[ipertd * 3 - 3];
					pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450005] =
						ppd[ipertd * 3 - 2];
					pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450004] =
						ppd[ipertd * 3 - 1];
					pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = 1.8756f;
					pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = lbpd[ipertd - 1];
					pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450006] =
						aa_1.r__[*i1 * 3 - 3];
					pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450005] =
						aa_1.r__[*i1 * 3 - 2];
					pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450004] =
						aa_1.r__[*i1 * 3 - 1];
					dpert_1.dppion[nn_1.nnn + *irun * 150001 - 150002] = 1.f /
																		 (float)para8_1.npertd;
				}
			}
		}
	}
	*iblock = 501;
	return 0;
}

int init_(int *minnum, int *maxnum, int *num,
		  float *radius, float *x0, float *z0, float *p0, float *gamma, int *iseed, int *mass, int *iopt)
{
	int i__1, i__2;
	float r__1, r__2, r__3;
	double d__1;

	double sqrt(double), exp(double), pow_dd(double *, double *);

	static int i__;
	static float x, y, z__, px, py, pz, beta;
	static int idir;
	static float sign;
	static int irun;
	static float ptot[3], rhow0, epart;
	static int idnum, npart;
	static float rdist, rhows, pfermi;
	extern double ranart_(int *);

	if (*p0 != 0.f)
	{
		sign = *p0 / dabs(*p0);
	}
	else
	{
		sign = 0.f;
	}
	r__1 = *gamma;
	beta = sign * sqrt(r__1 * r__1 - 1.f) / *gamma;
	if (*minnum == 1)
	{
		idnum = 1;
	}
	else
	{
		idnum = -1;
	}
	i__1 = *num;
	for (irun = 1; irun <= i__1; ++irun)
	{
		i__2 = *maxnum + (irun - 1) * *mass;
		for (i__ = *minnum + (irun - 1) * *mass; i__ <= i__2; ++i__)
		{
			ee_1.id[i__ - 1] = idnum;
			cc_1.e[i__ - 1] = .9383f;
		}
		i__2 = *maxnum + (irun - 1) * *mass;
		for (i__ = *minnum + (irun - 1) * *mass; i__ <= i__2; ++i__)
		{
		L200:
			x = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
			y = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
			z__ = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
			if (x * x + y * y + z__ * z__ > 1.f)
			{
				goto L200;
			}
			aa_1.r__[i__ * 3 - 3] = x * *radius;
			aa_1.r__[i__ * 3 - 2] = y * *radius;
			aa_1.r__[i__ * 3 - 1] = z__ * *radius;
		}
	}
	if (*iopt != 3)
	{
		rhow0 = .168f;
		i__1 = *num;
		for (irun = 1; irun <= i__1; ++irun)
		{
			i__2 = *maxnum + (irun - 1) * *mass;
			for (i__ = *minnum + (irun - 1) * *mass; i__ <= i__2; ++i__)
			{
			L500:
				px = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
				py = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
				pz = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
				if (px * px + py * py + pz * pz > 1.f)
				{
					goto L500;
				}
				r__1 = aa_1.r__[i__ * 3 - 3];
				r__2 = aa_1.r__[i__ * 3 - 2];
				r__3 = aa_1.r__[i__ * 3 - 1];
				rdist = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
				rhows = rhow0 / (exp((rdist - *radius) / .55f) + 1.f);
				d__1 = (double)(rhows * 14.804406096562142f);
				pfermi = pow_dd(&d__1, &c_b5) * .197f;
				if (*iopt == 2)
				{
					pfermi = .27f;
				}
				if (*iopt == 4)
				{
					pfermi = 0.f;
				}
				bb_1.p[i__ * 3 - 3] = pfermi * px;
				bb_1.p[i__ * 3 - 2] = pfermi * py;
				bb_1.p[i__ * 3 - 1] = pfermi * pz;
			}

			for (idir = 1; idir <= 3; ++idir)
			{
				ptot[idir - 1] = 0.f;
			}
			npart = 0;
			i__2 = *maxnum + (irun - 1) * *mass;
			for (i__ = *minnum + (irun - 1) * *mass; i__ <= i__2; ++i__)
			{
				++npart;
				for (idir = 1; idir <= 3; ++idir)
				{
					ptot[idir - 1] += bb_1.p[idir + i__ * 3 - 4];
				}
			}
			i__2 = *maxnum + (irun - 1) * *mass;
			for (i__ = *minnum + (irun - 1) * *mass; i__ <= i__2; ++i__)
			{
				for (idir = 1; idir <= 3; ++idir)
				{
					bb_1.p[idir + i__ * 3 - 4] -= ptot[idir - 1] / (float)
																	   npart;
				}
				if (*iopt == 1 || *iopt == 2)
				{
					r__1 = bb_1.p[i__ * 3 - 3];
					r__2 = bb_1.p[i__ * 3 - 2];
					r__3 = bb_1.p[i__ * 3 - 1];
					epart = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 +
								 .88040689000000005f);
					bb_1.p[i__ * 3 - 1] = *gamma * (bb_1.p[i__ * 3 - 1] +
													beta * epart);
				}
				else
				{
					bb_1.p[i__ * 3 - 1] += *p0;
				}
			}
		}
	}
	else
	{
		i__1 = *num;
		for (irun = 1; irun <= i__1; ++irun)
		{
			i__2 = *maxnum + (irun - 1) * *mass;
			for (i__ = *minnum + (irun - 1) * *mass; i__ <= i__2; ++i__)
			{
				bb_1.p[i__ * 3 - 3] = 0.f;
				bb_1.p[i__ * 3 - 2] = 0.f;
				bb_1.p[i__ * 3 - 1] = *p0;
			}
		}
	}
	i__1 = *num;
	for (irun = 1; irun <= i__1; ++irun)
	{
		i__2 = *maxnum + (irun - 1) * *mass;
		for (i__ = *minnum + (irun - 1) * *mass; i__ <= i__2; ++i__)
		{
			aa_1.r__[i__ * 3 - 3] += *x0;
			aa_1.r__[i__ * 3 - 1] = (aa_1.r__[i__ * 3 - 1] + *z0) / *gamma;
		}
	}

	return 0;
}

int dens_(int *ipot, int *mass, int *num,
		  int *nesc)
{
	static float zet[91] = {1.f, 0.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
						   0.f, 0.f, 0.f, -1.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
						   -1.f, 0.f, 1.f, 0.f, -1.f, 0.f, -1.f, 0.f, -2.f, -1.f, 0.f, 1.f, 0.f, 0.f, 0.f,
						   0.f, -1.f, 0.f, 1.f, 0.f, -1.f, 0.f, 1.f, -1.f, 0.f, 1.f, 2.f, 0.f, 1.f, 0.f,
						   1.f, 0.f, -1.f, 0.f, 1.f, 0.f, 0.f, 0.f, -1.f, 0.f, 1.f, 0.f, -1.f, 0.f, 1.f,
						   0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, -1.f, 0.f, 0.f, 0.f,
						   0.f, -1.f};

	int i__1, i__2;
	float r__1, r__2, r__3, r__4;
	double d__1, d__2;

	int i_nint(float *);
	double sqrt(double), pow_dd(double *, double *);

	static float a, b;
	static int i__, j;
	static float s, u;
	static int ix, iy, iz;
	static float big, pxl[82369], pyl[82369],
		pzl[82369], rho0,
		denr;
	static int irun, msum;
	static float gamma, small, smass, smass2;

	for (iz = -24; iz <= 24; ++iz)
	{
		for (iy = -20; iy <= 20; ++iy)
		{
			for (ix = -20; ix <= 20; ++ix)
			{
				dd_1.rho[ix + (iy + iz * 41) * 41 + 41184] = 0.f;
				dd_1.rhon[ix + (iy + iz * 41) * 41 + 41184] = 0.f;
				dd_1.rhop[ix + (iy + iz * 41) * 41 + 41184] = 0.f;
				ddpi_1.pirho[ix + (iy + iz * 41) * 41 + 41184] = 0.f;
				pxl[ix + (iy + iz * 41) * 41 + 41184] = 0.f;
				pyl[ix + (iy + iz * 41) * 41 + 41184] = 0.f;
				pzl[ix + (iy + iz * 41) * 41 + 41184] = 0.f;
				tt_1.pel[ix + (iy + iz * 41) * 41 + 41184] = 0.f;
				bbb_1.bxx[ix + (iy + iz * 41) * 41 + 41184] = 0.f;
				bbb_1.byy[ix + (iy + iz * 41) * 41 + 41184] = 0.f;
				bbb_1.bzz[ix + (iy + iz * 41) * 41 + 41184] = 0.f;
			}
		}
	}

	*nesc = 0;
	big = 1.f / ((float)(*num) * 3.f);
	small = 1.f / ((float)(*num) * 9.f);

	msum = 0;
	i__1 = *num;
	for (irun = 1; irun <= i__1; ++irun)
	{
		msum += rr_1.massr[irun - 1];
		i__2 = rr_1.massr[irun];
		for (j = 1; j <= i__2; ++j)
		{
			i__ = j + msum;
			ix = i_nint(&aa_1.r__[i__ * 3 - 3]);
			iy = i_nint(&aa_1.r__[i__ * 3 - 2]);
			iz = i_nint(&aa_1.r__[i__ * 3 - 1]);
			if (ix <= -20 || ix >= 20 || iy <= -20 || iy >= 20 || iz <= -24 ||
				iz >= 24)
			{
				++(*nesc);
			}
			else
			{

				if (j > *mass)
				{
					goto L30;
				}
				dd_1.rho[ix + (iy + iz * 41) * 41 + 41184] += big;
				dd_1.rho[ix + 1 + (iy + iz * 41) * 41 + 41184] += small;
				dd_1.rho[ix - 1 + (iy + iz * 41) * 41 + 41184] += small;
				dd_1.rho[ix + (iy + 1 + iz * 41) * 41 + 41184] += small;
				dd_1.rho[ix + (iy - 1 + iz * 41) * 41 + 41184] += small;
				dd_1.rho[ix + (iy + (iz + 1) * 41) * 41 + 41184] += small;
				dd_1.rho[ix + (iy + (iz - 1) * 41) * 41 + 41184] += small;
				if (zet[ee_1.lb[i__ - 1] + 45] != 0.f)
				{
					dd_1.rhop[ix + (iy + iz * 41) * 41 + 41184] += big;
					dd_1.rhop[ix + 1 + (iy + iz * 41) * 41 + 41184] += small;
					dd_1.rhop[ix - 1 + (iy + iz * 41) * 41 + 41184] += small;
					dd_1.rhop[ix + (iy + 1 + iz * 41) * 41 + 41184] += small;
					dd_1.rhop[ix + (iy - 1 + iz * 41) * 41 + 41184] += small;
					dd_1.rhop[ix + (iy + (iz + 1) * 41) * 41 + 41184] +=
						small;
					dd_1.rhop[ix + (iy + (iz - 1) * 41) * 41 + 41184] +=
						small;
					goto L40;
				}
				if (zet[ee_1.lb[i__ - 1] + 45] == 0.f)
				{
					dd_1.rhon[ix + (iy + iz * 41) * 41 + 41184] += big;
					dd_1.rhon[ix + 1 + (iy + iz * 41) * 41 + 41184] += small;
					dd_1.rhon[ix - 1 + (iy + iz * 41) * 41 + 41184] += small;
					dd_1.rhon[ix + (iy + 1 + iz * 41) * 41 + 41184] += small;
					dd_1.rhon[ix + (iy - 1 + iz * 41) * 41 + 41184] += small;
					dd_1.rhon[ix + (iy + (iz + 1) * 41) * 41 + 41184] +=
						small;
					dd_1.rhon[ix + (iy + (iz - 1) * 41) * 41 + 41184] +=
						small;
					goto L40;
				}
			L30:
				ddpi_1.pirho[ix + (iy + iz * 41) * 41 + 41184] += big;
				ddpi_1.pirho[ix + 1 + (iy + iz * 41) * 41 + 41184] += small;
				ddpi_1.pirho[ix - 1 + (iy + iz * 41) * 41 + 41184] += small;
				ddpi_1.pirho[ix + (iy + 1 + iz * 41) * 41 + 41184] += small;
				ddpi_1.pirho[ix + (iy - 1 + iz * 41) * 41 + 41184] += small;
				ddpi_1.pirho[ix + (iy + (iz + 1) * 41) * 41 + 41184] += small;
				ddpi_1.pirho[ix + (iy + (iz - 1) * 41) * 41 + 41184] += small;
			L40:
				pxl[ix + (iy + iz * 41) * 41 + 41184] += bb_1.p[i__ * 3 - 3] *
														 big;
				pxl[ix + 1 + (iy + iz * 41) * 41 + 41184] += bb_1.p[i__ * 3 -
																	3] *
															 small;
				pxl[ix - 1 + (iy + iz * 41) * 41 + 41184] += bb_1.p[i__ * 3 -
																	3] *
															 small;
				pxl[ix + (iy + 1 + iz * 41) * 41 + 41184] += bb_1.p[i__ * 3 -
																	3] *
															 small;
				pxl[ix + (iy - 1 + iz * 41) * 41 + 41184] += bb_1.p[i__ * 3 -
																	3] *
															 small;
				pxl[ix + (iy + (iz + 1) * 41) * 41 + 41184] += bb_1.p[i__ * 3 - 3] * small;
				pxl[ix + (iy + (iz - 1) * 41) * 41 + 41184] += bb_1.p[i__ * 3 - 3] * small;
				pyl[ix + (iy + iz * 41) * 41 + 41184] += bb_1.p[i__ * 3 - 2] *
														 big;
				pyl[ix + 1 + (iy + iz * 41) * 41 + 41184] += bb_1.p[i__ * 3 -
																	2] *
															 small;
				pyl[ix - 1 + (iy + iz * 41) * 41 + 41184] += bb_1.p[i__ * 3 -
																	2] *
															 small;
				pyl[ix + (iy + 1 + iz * 41) * 41 + 41184] += bb_1.p[i__ * 3 -
																	2] *
															 small;
				pyl[ix + (iy - 1 + iz * 41) * 41 + 41184] += bb_1.p[i__ * 3 -
																	2] *
															 small;
				pyl[ix + (iy + (iz + 1) * 41) * 41 + 41184] += bb_1.p[i__ * 3 - 2] * small;
				pyl[ix + (iy + (iz - 1) * 41) * 41 + 41184] += bb_1.p[i__ * 3 - 2] * small;
				pzl[ix + (iy + iz * 41) * 41 + 41184] += bb_1.p[i__ * 3 - 1] *
														 big;
				pzl[ix + 1 + (iy + iz * 41) * 41 + 41184] += bb_1.p[i__ * 3 -
																	1] *
															 small;
				pzl[ix - 1 + (iy + iz * 41) * 41 + 41184] += bb_1.p[i__ * 3 -
																	1] *
															 small;
				pzl[ix + (iy + 1 + iz * 41) * 41 + 41184] += bb_1.p[i__ * 3 -
																	1] *
															 small;
				pzl[ix + (iy - 1 + iz * 41) * 41 + 41184] += bb_1.p[i__ * 3 -
																	1] *
															 small;
				pzl[ix + (iy + (iz + 1) * 41) * 41 + 41184] += bb_1.p[i__ * 3 - 1] * small;
				pzl[ix + (iy + (iz - 1) * 41) * 41 + 41184] += bb_1.p[i__ * 3 - 1] * small;
				r__1 = cc_1.e[i__ - 1];
				r__2 = bb_1.p[i__ * 3 - 3];
				r__3 = bb_1.p[i__ * 3 - 2];
				r__4 = bb_1.p[i__ * 3 - 1];
				tt_1.pel[ix + (iy + iz * 41) * 41 + 41184] += sqrt(r__1 *
																	   r__1 +
																   r__2 * r__2 + r__3 * r__3 + r__4 * r__4) *
															  big;
				r__1 = cc_1.e[i__ - 1];
				r__2 = bb_1.p[i__ * 3 - 3];
				r__3 = bb_1.p[i__ * 3 - 2];
				r__4 = bb_1.p[i__ * 3 - 1];
				tt_1.pel[ix + 1 + (iy + iz * 41) * 41 + 41184] += sqrt(r__1 *
																		   r__1 +
																	   r__2 * r__2 + r__3 * r__3 + r__4 * r__4) *
																  small;
				r__1 = cc_1.e[i__ - 1];
				r__2 = bb_1.p[i__ * 3 - 3];
				r__3 = bb_1.p[i__ * 3 - 2];
				r__4 = bb_1.p[i__ * 3 - 1];
				tt_1.pel[ix - 1 + (iy + iz * 41) * 41 + 41184] += sqrt(r__1 *
																		   r__1 +
																	   r__2 * r__2 + r__3 * r__3 + r__4 * r__4) *
																  small;
				r__1 = cc_1.e[i__ - 1];
				r__2 = bb_1.p[i__ * 3 - 3];
				r__3 = bb_1.p[i__ * 3 - 2];
				r__4 = bb_1.p[i__ * 3 - 1];
				tt_1.pel[ix + (iy + 1 + iz * 41) * 41 + 41184] += sqrt(r__1 *
																		   r__1 +
																	   r__2 * r__2 + r__3 * r__3 + r__4 * r__4) *
																  small;
				r__1 = cc_1.e[i__ - 1];
				r__2 = bb_1.p[i__ * 3 - 3];
				r__3 = bb_1.p[i__ * 3 - 2];
				r__4 = bb_1.p[i__ * 3 - 1];
				tt_1.pel[ix + (iy - 1 + iz * 41) * 41 + 41184] += sqrt(r__1 *
																		   r__1 +
																	   r__2 * r__2 + r__3 * r__3 + r__4 * r__4) *
																  small;
				r__1 = cc_1.e[i__ - 1];
				r__2 = bb_1.p[i__ * 3 - 3];
				r__3 = bb_1.p[i__ * 3 - 2];
				r__4 = bb_1.p[i__ * 3 - 1];
				tt_1.pel[ix + (iy + (iz + 1) * 41) * 41 + 41184] += sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4) *
																	small;
				r__1 = cc_1.e[i__ - 1];
				r__2 = bb_1.p[i__ * 3 - 3];
				r__3 = bb_1.p[i__ * 3 - 2];
				r__4 = bb_1.p[i__ * 3 - 1];
				tt_1.pel[ix + (iy + (iz - 1) * 41) * 41 + 41184] += sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4) *
																	small;
			}
		}
	}

	for (iz = -24; iz <= 24; ++iz)
	{
		for (iy = -20; iy <= 20; ++iy)
		{
			for (ix = -20; ix <= 20; ++ix)
			{
				if (dd_1.rho[ix + (iy + iz * 41) * 41 + 41184] == 0.f ||
					tt_1.pel[ix + (iy + iz * 41) * 41 + 41184] == 0.f)
				{
					goto L101;
				}
				r__1 = tt_1.pel[ix + (iy + iz * 41) * 41 + 41184];
				r__2 = pxl[ix + (iy + iz * 41) * 41 + 41184];
				r__3 = pyl[ix + (iy + iz * 41) * 41 + 41184];
				r__4 = pzl[ix + (iy + iz * 41) * 41 + 41184];
				smass2 = r__1 * r__1 - r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
				if (smass2 <= 0.f)
				{
					smass2 = 1e-6f;
				}
				smass = sqrt(smass2);
				if (smass == 0.f)
				{
					smass = 1e-6f;
				}
				gamma = tt_1.pel[ix + (iy + iz * 41) * 41 + 41184] / smass;
				if (gamma == 0.f)
				{
					goto L101;
				}
				bbb_1.bxx[ix + (iy + iz * 41) * 41 + 41184] = pxl[ix + (iy + iz * 41) * 41 + 41184] / tt_1.pel[ix + (iy + iz * 41) * 41 + 41184];
				bbb_1.byy[ix + (iy + iz * 41) * 41 + 41184] = pyl[ix + (iy + iz * 41) * 41 + 41184] / tt_1.pel[ix + (iy + iz * 41) * 41 + 41184];
				bbb_1.bzz[ix + (iy + iz * 41) * 41 + 41184] = pzl[ix + (iy + iz * 41) * 41 + 41184] / tt_1.pel[ix + (iy + iz * 41) * 41 + 41184];
				dd_1.rho[ix + (iy + iz * 41) * 41 + 41184] /= gamma;
				dd_1.rhon[ix + (iy + iz * 41) * 41 + 41184] /= gamma;
				dd_1.rhop[ix + (iy + iz * 41) * 41 + 41184] /= gamma;
				ddpi_1.pirho[ix + (iy + iz * 41) * 41 + 41184] /= gamma;
				r__1 = gamma;
				tt_1.pel[ix + (iy + iz * 41) * 41 + 41184] /= r__1 * r__1;
				rho0 = .163f;
				if (*ipot == 0)
				{
					u = 0.f;
					goto L70;
				}
				if (*ipot == 1 || *ipot == 6)
				{
					a = -.1236f;
					b = .0704f;
					s = 2.f;
					goto L60;
				}
				if (*ipot == 2 || *ipot == 7)
				{
					a = -.218f;
					b = .164f;
					s = 1.3333333333333333f;
				}
				if (*ipot == 3)
				{
					a = -.3581f;
					b = .3048f;
					s = 1.167f;
					goto L60;
				}
				if (*ipot == 4)
				{
					denr = dd_1.rho[ix + (iy + iz * 41) * 41 + 41184] / rho0;
					b = .3048f;
					s = 1.167f;
					if (denr <= 4.f || denr > 7.f)
					{
						a = -.3581f;
					}
					else
					{
						d__1 = (double)denr;
						d__2 = (double)denr;
						a = -b * pow_dd(&d__1, &c_b342) - pow_dd(&d__2, &c_b343) * .023999999999999997f;
					}
					goto L60;
				}
			L60:
				r__1 = dd_1.rho[ix + (iy + iz * 41) * 41 + 41184];
				d__1 = (double)(dd_1.rho[ix + (iy + iz * 41) * 41 +
											 41184] /
									rho0);
				d__2 = (double)s;
				u = a * .5f * (r__1 * r__1) / rho0 + b / (s + 1) * pow_dd(&d__1, &d__2) * dd_1.rho[ix + (iy + iz * 41) * 41 + 41184];
			L70:
				tt_1.pel[ix + (iy + iz * 41) * 41 + 41184] += u;
			L101:;
			}
		}
	}
	return 0;
}

int gradu_(int *iopt, int *ix, int *iy, int *iz, float *gradx, float *grady, float *gradz)
{
	float r__1, r__2;
	double d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8;

	double pow_dd(double *, double *);

	static float ef, eh, cf0, den0, ene0, denr, eqgp, acoef, expnt, acoef1,
		acoef2, expnt2, rxmins, rymins, rzmins, rxplus, ryplus, rzplus;

	rxplus = dd_1.rho[*ix + 1 + (*iy + *iz * 41) * 41 + 41184] / .167f;
	rxmins = dd_1.rho[*ix - 1 + (*iy + *iz * 41) * 41 + 41184] / .167f;
	ryplus = dd_1.rho[*ix + (*iy + 1 + *iz * 41) * 41 + 41184] / .167f;
	rymins = dd_1.rho[*ix + (*iy - 1 + *iz * 41) * 41 + 41184] / .167f;
	rzplus = dd_1.rho[*ix + (*iy + (*iz + 1) * 41) * 41 + 41184] / .167f;
	rzmins = dd_1.rho[*ix + (*iy + (*iz - 1) * 41) * 41 + 41184] / .167f;
	den0 = dd_1.rho[*ix + (*iy + *iz * 41) * 41 + 41184] / .167f;
	ene0 = tt_1.pel[*ix + (*iy + *iz * 41) * 41 + 41184];
	switch (*iopt)
	{
	case 1:
		goto L1;
	case 2:
		goto L2;
	case 3:
		goto L3;
	case 4:
		goto L4;
	case 5:
		goto L5;
	}
	if (*iopt == 6)
	{
		goto L6;
	}
	if (*iopt == 7)
	{
		goto L7;
	}

L1:
	r__1 = rxplus;
	r__2 = rxmins;
	*gradx = (rxplus - rxmins) * -.062f + (r__1 * r__1 - r__2 * r__2) *
											  .03525f;
	r__1 = ryplus;
	r__2 = rymins;
	*grady = (ryplus - rymins) * -.062f + (r__1 * r__1 - r__2 * r__2) *
											  .03525f;
	r__1 = rzplus;
	r__2 = rzmins;
	*gradz = (rzplus - rzmins) * -.062f + (r__1 * r__1 - r__2 * r__2) *
											  .03525f;
	return 0;

L2:
	expnt = 1.3333333f;
	d__1 = (double)rxplus;
	d__2 = (double)expnt;
	d__3 = (double)rxmins;
	d__4 = (double)expnt;
	*gradx = (rxplus - rxmins) * -.109f + (pow_dd(&d__1, &d__2) - pow_dd(&d__3, &d__4)) * .082f;
	d__1 = (double)ryplus;
	d__2 = (double)expnt;
	d__3 = (double)rymins;
	d__4 = (double)expnt;
	*grady = (ryplus - rymins) * -.109f + (pow_dd(&d__1, &d__2) - pow_dd(&d__3, &d__4)) * .082f;
	d__1 = (double)rzplus;
	d__2 = (double)expnt;
	d__3 = (double)rzmins;
	d__4 = (double)expnt;
	*gradz = (rzplus - rzmins) * -.109f + (pow_dd(&d__1, &d__2) - pow_dd(&d__3, &d__4)) * .082f;
	return 0;

L3:
	expnt = 1.1666667f;
	acoef = .178f;
	d__1 = (double)rxplus;
	d__2 = (double)expnt;
	d__3 = (double)rxmins;
	d__4 = (double)expnt;
	*gradx = -acoef * (rxplus - rxmins) + (pow_dd(&d__1, &d__2) - pow_dd(&d__3, &d__4)) * .1515f;
	d__1 = (double)ryplus;
	d__2 = (double)expnt;
	d__3 = (double)rymins;
	d__4 = (double)expnt;
	*grady = -acoef * (ryplus - rymins) + (pow_dd(&d__1, &d__2) - pow_dd(&d__3, &d__4)) * .1515f;
	d__1 = (double)rzplus;
	d__2 = (double)expnt;
	d__3 = (double)rzmins;
	d__4 = (double)expnt;
	*gradz = -acoef * (rzplus - rzmins) + (pow_dd(&d__1, &d__2) - pow_dd(&d__3, &d__4)) * .1515f;
	return 0;

L4:
	eh = 4.f;
	eqgp = 7.f;
	acoef = .178f;
	expnt = 1.1666667f;
	denr = dd_1.rho[*ix + (*iy + *iz * 41) * 41 + 41184] / .167f;
	if (denr <= eh || denr >= eqgp)
	{
		d__1 = (double)rxplus;
		d__2 = (double)expnt;
		d__3 = (double)rxmins;
		d__4 = (double)expnt;
		*gradx = -acoef * (rxplus - rxmins) + (pow_dd(&d__1, &d__2) - pow_dd(&d__3, &d__4)) * .1515f;
		d__1 = (double)ryplus;
		d__2 = (double)expnt;
		d__3 = (double)rymins;
		d__4 = (double)expnt;
		*grady = -acoef * (ryplus - rymins) + (pow_dd(&d__1, &d__2) - pow_dd(&d__3, &d__4)) * .1515f;
		d__1 = (double)rzplus;
		d__2 = (double)expnt;
		d__3 = (double)rzmins;
		d__4 = (double)expnt;
		*gradz = -acoef * (rzplus - rzmins) + (pow_dd(&d__1, &d__2) - pow_dd(&d__3, &d__4)) * .1515f;
	}
	else
	{
		acoef1 = .178f;
		acoef2 = 0.f;
		expnt2 = .66666666666666663f;
		d__1 = (double)rxplus;
		d__2 = (double)expnt;
		d__3 = (double)rxmins;
		d__4 = (double)expnt;
		d__5 = (double)rxplus;
		d__6 = (double)expnt2;
		d__7 = (double)rxmins;
		d__8 = (double)expnt2;
		*gradx = -acoef1 * (pow_dd(&d__1, &d__2) - pow_dd(&d__3, &d__4)) -
				 acoef2 * (pow_dd(&d__5, &d__6) - pow_dd(&d__7, &d__8));
		d__1 = (double)ryplus;
		d__2 = (double)expnt;
		d__3 = (double)rymins;
		d__4 = (double)expnt;
		d__5 = (double)ryplus;
		d__6 = (double)expnt2;
		d__7 = (double)rymins;
		d__8 = (double)expnt2;
		*grady = -acoef1 * (pow_dd(&d__1, &d__2) - pow_dd(&d__3, &d__4)) -
				 acoef2 * (pow_dd(&d__5, &d__6) - pow_dd(&d__7, &d__8));
		d__1 = (double)rzplus;
		d__2 = (double)expnt;
		d__3 = (double)rzmins;
		d__4 = (double)expnt;
		d__5 = (double)rzplus;
		d__6 = (double)expnt2;
		d__7 = (double)rzmins;
		d__8 = (double)expnt2;
		*gradz = -acoef1 * (pow_dd(&d__1, &d__2) - pow_dd(&d__3, &d__4)) -
				 acoef2 * (pow_dd(&d__5, &d__6) - pow_dd(&d__7, &d__8));
	}
	return 0;

L5:
	expnt = 2.77f;
	d__1 = (double)rxplus;
	d__2 = (double)expnt;
	d__3 = (double)rxmins;
	d__4 = (double)expnt;
	*gradx = (rxplus - rxmins) * -.0516f + (pow_dd(&d__1, &d__2) - pow_dd(&d__3, &d__4)) * .02498f;
	d__1 = (double)ryplus;
	d__2 = (double)expnt;
	d__3 = (double)rymins;
	d__4 = (double)expnt;
	*grady = (ryplus - rymins) * -.0516f + (pow_dd(&d__1, &d__2) - pow_dd(&d__3, &d__4)) * .02498f;
	d__1 = (double)rzplus;
	d__2 = (double)expnt;
	d__3 = (double)rzmins;
	d__4 = (double)expnt;
	*gradz = (rzplus - rzmins) * -.0516f + (pow_dd(&d__1, &d__2) - pow_dd(&d__3, &d__4)) * .02498f;
	return 0;

L6:
	if (ene0 <= .5f)
	{
		r__1 = rxplus;
		r__2 = rxmins;
		*gradx = (rxplus - rxmins) * -.062f + (r__1 * r__1 - r__2 * r__2) *
												  .03525f;
		r__1 = ryplus;
		r__2 = rymins;
		*grady = (ryplus - rymins) * -.062f + (r__1 * r__1 - r__2 * r__2) *
												  .03525f;
		r__1 = rzplus;
		r__2 = rzmins;
		*gradz = (rzplus - rzmins) * -.062f + (r__1 * r__1 - r__2 * r__2) *
												  .03525f;
		return 0;
	}
	if (ene0 > .5f && ene0 <= 1.5f)
	{
		ef = .035999999999999997f;
		d__1 = (double)rxplus;
		d__2 = (double)rxmins;
		*gradx = ef * -.5f * (pow_dd(&d__1, &c_b352) - pow_dd(&d__2, &c_b352));
		d__1 = (double)ryplus;
		d__2 = (double)rymins;
		*grady = ef * -.5f * (pow_dd(&d__1, &c_b352) - pow_dd(&d__2, &c_b352));
		d__1 = (double)rzplus;
		d__2 = (double)rzmins;
		*gradz = ef * -.5f * (pow_dd(&d__1, &c_b352) - pow_dd(&d__2, &c_b352));
		return 0;
	}
	if (ene0 > 1.5f)
	{
		ef = .035999999999999997f;
		cf0 = .8f;
		d__1 = (double)rxplus;
		d__2 = (double)rxmins;
		d__3 = (double)rxplus;
		d__4 = (double)rxmins;
		*gradx = cf0 * .5f * (pow_dd(&d__1, &c_b358) - pow_dd(&d__2, &c_b358)) - ef * .5f * (pow_dd(&d__3, &c_b352) - pow_dd(&d__4, &c_b352));
		d__1 = (double)ryplus;
		d__2 = (double)rymins;
		d__3 = (double)ryplus;
		d__4 = (double)rymins;
		*grady = cf0 * .5f * (pow_dd(&d__1, &c_b358) - pow_dd(&d__2, &c_b358)) - ef * .5f * (pow_dd(&d__3, &c_b352) - pow_dd(&d__4, &c_b352));
		d__1 = (double)rzplus;
		d__2 = (double)rzmins;
		d__3 = (double)rzplus;
		d__4 = (double)rzmins;
		*gradz = cf0 * .5f * (pow_dd(&d__1, &c_b358) - pow_dd(&d__2, &c_b358)) - ef * .5f * (pow_dd(&d__3, &c_b352) - pow_dd(&d__4, &c_b352));
		return 0;
	}

L7:
	if (den0 <= 4.5f)
	{
		expnt = 1.1666667f;
		acoef = .178f;
		d__1 = (double)rxplus;
		d__2 = (double)expnt;
		d__3 = (double)rxmins;
		d__4 = (double)expnt;
		*gradx = -acoef * (rxplus - rxmins) + (pow_dd(&d__1, &d__2) - pow_dd(&d__3, &d__4)) * .1515f;
		d__1 = (double)ryplus;
		d__2 = (double)expnt;
		d__3 = (double)rymins;
		d__4 = (double)expnt;
		*grady = -acoef * (ryplus - rymins) + (pow_dd(&d__1, &d__2) - pow_dd(&d__3, &d__4)) * .1515f;
		d__1 = (double)rzplus;
		d__2 = (double)expnt;
		d__3 = (double)rzmins;
		d__4 = (double)expnt;
		*gradz = -acoef * (rzplus - rzmins) + (pow_dd(&d__1, &d__2) - pow_dd(&d__3, &d__4)) * .1515f;
		return 0;
	}
	if (den0 > 4.5f && den0 <= 5.1f)
	{
		ef = .035999999999999997f;
		d__1 = (double)rxplus;
		d__2 = (double)rxmins;
		*gradx = ef * -.5f * (pow_dd(&d__1, &c_b352) - pow_dd(&d__2, &c_b352));
		d__1 = (double)ryplus;
		d__2 = (double)rymins;
		*grady = ef * -.5f * (pow_dd(&d__1, &c_b352) - pow_dd(&d__2, &c_b352));
		d__1 = (double)rzplus;
		d__2 = (double)rzmins;
		*gradz = ef * -.5f * (pow_dd(&d__1, &c_b352) - pow_dd(&d__2, &c_b352));
		return 0;
	}
	if (den0 > 5.1f)
	{
		ef = .035999999999999997f;
		cf0 = .8f;
		d__1 = (double)rxplus;
		d__2 = (double)rxmins;
		d__3 = (double)rxplus;
		d__4 = (double)rxmins;
		*gradx = cf0 * .5f * (pow_dd(&d__1, &c_b358) - pow_dd(&d__2, &c_b358)) - ef * .5f * (pow_dd(&d__3, &c_b352) - pow_dd(&d__4, &c_b352));
		d__1 = (double)ryplus;
		d__2 = (double)rymins;
		d__3 = (double)ryplus;
		d__4 = (double)rymins;
		*grady = cf0 * .5f * (pow_dd(&d__1, &c_b358) - pow_dd(&d__2, &c_b358)) - ef * .5f * (pow_dd(&d__3, &c_b352) - pow_dd(&d__4, &c_b352));
		d__1 = (double)rzplus;
		d__2 = (double)rzmins;
		d__3 = (double)rzplus;
		d__4 = (double)rzmins;
		*gradz = cf0 * .5f * (pow_dd(&d__1, &c_b358) - pow_dd(&d__2, &c_b358)) - ef * .5f * (pow_dd(&d__3, &c_b352) - pow_dd(&d__4, &c_b352));
		return 0;
	}
	return 0;
}

int graduk_(int *ix, int *iy, int *iz, float *gradxk, float *gradyk, float *gradzk)
{
	static float rxmins, rymins, rzmins, rxplus, ryplus, rzplus;

	rxplus = dd_1.rho[*ix + 1 + (*iy + *iz * 41) * 41 + 41184];
	rxmins = dd_1.rho[*ix - 1 + (*iy + *iz * 41) * 41 + 41184];
	ryplus = dd_1.rho[*ix + (*iy + 1 + *iz * 41) * 41 + 41184];
	rymins = dd_1.rho[*ix + (*iy - 1 + *iz * 41) * 41 + 41184];
	rzplus = dd_1.rho[*ix + (*iy + (*iz + 1) * 41) * 41 + 41184];
	rzmins = dd_1.rho[*ix + (*iy + (*iz - 1) * 41) * 41 + 41184];
	*gradxk = (rxplus - rxmins) / 2.f;
	*gradyk = (ryplus - rymins) / 2.f;
	*gradzk = (rzplus - rzmins) / 2.f;
	return 0;
}

int gradup_(int *ix, int *iy, int *iz, float *gradxp, float *gradyp, float *gradzp)
{
	static float rxmins, rymins, rzmins, rxplus, ryplus, rzplus;

	rxplus = dd_1.rhop[*ix + 1 + (*iy + *iz * 41) * 41 + 41184] / .168f;
	rxmins = dd_1.rhop[*ix - 1 + (*iy + *iz * 41) * 41 + 41184] / .168f;
	ryplus = dd_1.rhop[*ix + (*iy + 1 + *iz * 41) * 41 + 41184] / .168f;
	rymins = dd_1.rhop[*ix + (*iy - 1 + *iz * 41) * 41 + 41184] / .168f;
	rzplus = dd_1.rhop[*ix + (*iy + (*iz + 1) * 41) * 41 + 41184] / .168f;
	rzmins = dd_1.rhop[*ix + (*iy + (*iz - 1) * 41) * 41 + 41184] / .168f;
	*gradxp = (rxplus - rxmins) / 2.f;
	*gradyp = (ryplus - rymins) / 2.f;
	*gradzp = (rzplus - rzmins) / 2.f;
	return 0;
}

int gradun_(int *ix, int *iy, int *iz, float *gradxn, float *gradyn, float *gradzn)
{
	static float rxmins, rymins, rzmins, rxplus, ryplus, rzplus;

	rxplus = dd_1.rhon[*ix + 1 + (*iy + *iz * 41) * 41 + 41184] / .168f;
	rxmins = dd_1.rhon[*ix - 1 + (*iy + *iz * 41) * 41 + 41184] / .168f;
	ryplus = dd_1.rhon[*ix + (*iy + 1 + *iz * 41) * 41 + 41184] / .168f;
	rymins = dd_1.rhon[*ix + (*iy - 1 + *iz * 41) * 41 + 41184] / .168f;
	rzplus = dd_1.rhon[*ix + (*iy + (*iz + 1) * 41) * 41 + 41184] / .168f;
	rzmins = dd_1.rhon[*ix + (*iy + (*iz - 1) * 41) * 41 + 41184] / .168f;
	*gradxn = (rxplus - rxmins) / 2.f;
	*gradyn = (ryplus - rymins) / 2.f;
	*gradzn = (rzplus - rzmins) / 2.f;
	return 0;
}

double fde_(float *dmass, float *srt, float *con)
{
	float ret_val, r__1, r__2, r__3, r__4, r__5, r__6;

	double sqrt(double);

	static float p1, fd, p11, am0, amn, avpi;
	extern double width_(float *);

	amn = .938869f;
	avpi = .13803333f;
	am0 = 1.232f;
	r__1 = am0;
	r__3 = *dmass;
	r__2 = r__3 * r__3 - 1.5178240000000001f;
	r__4 = am0;
	r__5 = width_(dmass);
	fd = r__1 * r__1 * 4.f * width_(dmass) / (r__2 * r__2 + r__4 * r__4 * (r__5 * r__5));
	if (*con == 1.f)
	{
		r__2 = *srt;
		r__3 = *dmass;
		r__4 = amn;
		r__1 = r__2 * r__2 + r__3 * r__3 - r__4 * r__4;
		r__5 = *srt;
		r__6 = *dmass;
		p11 = r__1 * r__1 / (r__5 * r__5 * 4.f) - r__6 * r__6;
		if (p11 <= 0.f)
		{
			p11 = 1e-6f;
		}
		p1 = sqrt(p11);
	}
	else
	{
		*dmass = amn + avpi;
		r__2 = *srt;
		r__3 = *dmass;
		r__4 = amn;
		r__1 = r__2 * r__2 + r__3 * r__3 - r__4 * r__4;
		r__5 = *srt;
		r__6 = *dmass;
		p11 = r__1 * r__1 / (r__5 * r__5 * 4.f) - r__6 * r__6;
		if (p11 <= 0.f)
		{
			p11 = 1e-6f;
		}
		p1 = sqrt(p11);
	}
	ret_val = fd * p1 * *dmass;
	return ret_val;
}

double fd5_(float *dmass, float *srt, float *con)
{
	float ret_val, r__1, r__2, r__3, r__4, r__5, r__6;

	double sqrt(double);

	static float p1, fd, am0, amn;
	extern double w1535_(float *);
	static float avpi;

	amn = .938869f;
	avpi = .13803333f;
	am0 = 1.535f;
	r__1 = am0;
	r__3 = *dmass;
	r__2 = r__3 * r__3 - 2.3562249999999998f;
	r__4 = am0;
	r__5 = w1535_(dmass);
	fd = r__1 * r__1 * 4.f * w1535_(dmass) / (r__2 * r__2 + r__4 * r__4 * (r__5 * r__5));
	if (*con == 1.f)
	{
		r__2 = *srt;
		r__3 = *dmass;
		r__4 = amn;
		r__1 = r__2 * r__2 + r__3 * r__3 - r__4 * r__4;
		r__5 = *srt;
		r__6 = *dmass;
		p1 = sqrt(r__1 * r__1 / (r__5 * r__5 * 4.f) - r__6 * r__6);
	}
	else
	{
		*dmass = amn + avpi;
		r__2 = *srt;
		r__3 = *dmass;
		r__4 = amn;
		r__1 = r__2 * r__2 + r__3 * r__3 - r__4 * r__4;
		r__5 = *srt;
		r__6 = *dmass;
		p1 = sqrt(r__1 * r__1 / (r__5 * r__5 * 4.f) - r__6 * r__6);
	}
	ret_val = fd * p1 * *dmass;
	return ret_val;
}

double fns_(float *dmass, float *srt, float *con)
{
	float ret_val, r__1, r__2, r__3, r__4, r__5, r__6;

	double sqrt(double);

	static float p1, fn, an0, amn, avpi, width;

	width = .2f;
	amn = .938869f;
	avpi = .13803333f;
	an0 = 1.43f;
	r__1 = an0;
	r__3 = *dmass;
	r__2 = r__3 * r__3 - 2.0735999999999999f;
	r__4 = an0;
	r__5 = width;
	fn = r__1 * r__1 * 4.f * width / (r__2 * r__2 + r__4 * r__4 * (r__5 * r__5));
	if (*con == 1.f)
	{
		r__2 = *srt;
		r__3 = *dmass;
		r__4 = amn;
		r__1 = r__2 * r__2 + r__3 * r__3 - r__4 * r__4;
		r__5 = *srt;
		r__6 = *dmass;
		p1 = sqrt(r__1 * r__1 / (r__5 * r__5 * 4.f) - r__6 * r__6);
	}
	else
	{
		*dmass = amn + avpi;
		r__2 = *srt;
		r__3 = *dmass;
		r__4 = amn;
		r__1 = r__2 * r__2 + r__3 * r__3 - r__4 * r__4;
		r__5 = *srt;
		r__6 = *dmass;
		p1 = sqrt(r__1 * r__1 / (r__5 * r__5 * 4.f) - r__6 * r__6);
	}
	ret_val = fn * p1 * *dmass;
	return ret_val;
}

int decay_(int *irun, int *i__, int *nnn,
		   int *iseed, float *wid, int *nt)
{
	int i__1;

	static float x3, x4, x5, x6, x8, dm;
	static int lbi, lbm, nlab, nalb;
	static float ctrl;
	extern int dkine_(int *, int *, int *,
					  int *, int *, float *, int *);
	static int lbanti, lbsave;
	static float dpsave;
	extern double ranart_(int *);
	static float xmsave, pxsave, pysave, pzsave;

	lbanti = ee_1.lb[*i__ - 1];

	dm = cc_1.e[*i__ - 1];
	if ((i__1 = ee_1.lb[*i__ - 1], abs(i__1)) == 11)
	{
		x3 = ranart_(&rndf77_1.nseed);
		if (x3 > .33333333333333331f)
		{
			ee_1.lb[*i__ - 1] = 2;
			nlab = 2;
			pd_1.lpion[*nnn + *irun * 150001 - 150002] = 5;
			pc_1.epion[*nnn + *irun * 150001 - 150002] = .13957f;
		}
		else
		{
			ee_1.lb[*i__ - 1] = 1;
			nlab = 1;
			pd_1.lpion[*nnn + *irun * 150001 - 150002] = 4;
			pc_1.epion[*nnn + *irun * 150001 - 150002] = .13496f;
		}
	}
	else if ((i__1 = ee_1.lb[*i__ - 1], abs(i__1)) == 10)
	{
		x4 = ranart_(&rndf77_1.nseed);
		if (x4 > .33333333333333331f)
		{
			ee_1.lb[*i__ - 1] = 1;
			nlab = 1;
			pd_1.lpion[*nnn + *irun * 150001 - 150002] = 3;
			pc_1.epion[*nnn + *irun * 150001 - 150002] = .13957f;
		}
		else
		{
			ee_1.lb[*i__ - 1] = 2;
			nalb = 2;
			pd_1.lpion[*nnn + *irun * 150001 - 150002] = 4;
			pc_1.epion[*nnn + *irun * 150001 - 150002] = .13496f;
		}
	}
	else if ((i__1 = ee_1.lb[*i__ - 1], abs(i__1)) == 12)
	{
		ctrl = .65f;
		if (dm <= 1.49f)
		{
			ctrl = -1.f;
		}
		x5 = ranart_(&rndf77_1.nseed);
		if (x5 >= ctrl)
		{
			x6 = ranart_(&rndf77_1.nseed);
			if (x6 > .33333333333333331f)
			{
				ee_1.lb[*i__ - 1] = 1;
				nlab = 1;
				pd_1.lpion[*nnn + *irun * 150001 - 150002] = 3;
				pc_1.epion[*nnn + *irun * 150001 - 150002] = .13957f;
			}
			else
			{
				ee_1.lb[*i__ - 1] = 2;
				nalb = 2;
				pd_1.lpion[*nnn + *irun * 150001 - 150002] = 4;
				pc_1.epion[*nnn + *irun * 150001 - 150002] = .13496f;
			}
		}
		else
		{
			ee_1.lb[*i__ - 1] = 2;
			nlab = 2;
			pd_1.lpion[*nnn + *irun * 150001 - 150002] = 0;
			pc_1.epion[*nnn + *irun * 150001 - 150002] = .5475f;
		}
	}
	else if ((i__1 = ee_1.lb[*i__ - 1], abs(i__1)) == 13)
	{
		ctrl = .65f;
		if (dm <= 1.49f)
		{
			ctrl = -1.f;
		}
		x5 = ranart_(&rndf77_1.nseed);
		if (x5 >= ctrl)
		{
			x8 = ranart_(&rndf77_1.nseed);
			if (x8 > .33333333333333331f)
			{
				ee_1.lb[*i__ - 1] = 2;
				nlab = 2;
				pd_1.lpion[*nnn + *irun * 150001 - 150002] = 5;
				pc_1.epion[*nnn + *irun * 150001 - 150002] = .13957f;
			}
			else
			{
				ee_1.lb[*i__ - 1] = 1;
				nlab = 1;
				pd_1.lpion[*nnn + *irun * 150001 - 150002] = 4;
				pc_1.epion[*nnn + *irun * 150001 - 150002] = .13496f;
			}
		}
		else
		{
			ee_1.lb[*i__ - 1] = 1;
			nlab = 1;
			pd_1.lpion[*nnn + *irun * 150001 - 150002] = 0;
			pc_1.epion[*nnn + *irun * 150001 - 150002] = .5475f;
		}
	}

	dkine_(irun, i__, nnn, &nlab, iseed, wid, nt);

	if (lbanti < 0)
	{
		lbi = ee_1.lb[*i__ - 1];
		if (lbi == 1 || lbi == 2)
		{
			lbi = -lbi;
		}
		else if (lbi == 3)
		{
			lbi = 5;
		}
		else if (lbi == 5)
		{
			lbi = 3;
		}
		ee_1.lb[*i__ - 1] = lbi;

		lbi = pd_1.lpion[*nnn + *irun * 150001 - 150002];
		if (lbi == 3)
		{
			lbi = 5;
		}
		else if (lbi == 5)
		{
			lbi = 3;
		}
		else if (lbi == 1 || lbi == 2)
		{
			lbi = -lbi;
		}
		pd_1.lpion[*nnn + *irun * 150001 - 150002] = lbi;
	}

	if (*nt == input2_1.ntmax)
	{
		lbm = pd_1.lpion[*nnn + *irun * 150001 - 150002];
		if (lbm == 0 || lbm == 25 || lbm == 26 || lbm == 27)
		{
			lbsave = lbm;
			xmsave = pc_1.epion[*nnn + *irun * 150001 - 150002];
			pxsave = pb_1.ppion[(*nnn + *irun * 150001) * 3 - 450006];
			pysave = pb_1.ppion[(*nnn + *irun * 150001) * 3 - 450005];
			pzsave = pb_1.ppion[(*nnn + *irun * 150001) * 3 - 450004];
			dpsave = dpert_1.dppion[*nnn + *irun * 150001 - 150002];
			pd_1.lpion[*nnn + *irun * 150001 - 150002] = ee_1.lb[*i__ - 1];
			pc_1.epion[*nnn + *irun * 150001 - 150002] = cc_1.e[*i__ - 1];
			pb_1.ppion[(*nnn + *irun * 150001) * 3 - 450006] = bb_1.p[*i__ *
																		  3 -
																	  3];
			pb_1.ppion[(*nnn + *irun * 150001) * 3 - 450005] = bb_1.p[*i__ *
																		  3 -
																	  2];
			pb_1.ppion[(*nnn + *irun * 150001) * 3 - 450004] = bb_1.p[*i__ *
																		  3 -
																	  1];
			dpert_1.dppion[*nnn + *irun * 150001 - 150002] = dpert_1.dpertp[*i__ - 1];
			ee_1.lb[*i__ - 1] = lbsave;
			cc_1.e[*i__ - 1] = xmsave;
			bb_1.p[*i__ * 3 - 3] = pxsave;
			bb_1.p[*i__ * 3 - 2] = pysave;
			bb_1.p[*i__ * 3 - 1] = pzsave;
			dpert_1.dpertp[*i__ - 1] = dpsave;
		}
	}
	return 0;
}

int dkine_(int *irun, int *i__, int *nnn,
		   int *nlab, int *iseed, float *wid, int *nt)
{
	float r__1, r__2, r__3, r__4, r__5, r__6, r__7, r__8;

	double sqrt(double), log(double);

	static float q, q2, gd, am, dm, en, ep, pm, qs, px, py, pz, rx, ry, rz, qx,
		qy, qz, fgd, bdx, bdy, bdz, bpp, bpn, pxn, pyn, pxp, pyp, pzp,
		pzn, tau0, devio, edelta;
	extern double ranart_(int *);
	static float taudcy;

	px = bb_1.p[*i__ * 3 - 3];
	py = bb_1.p[*i__ * 3 - 2];
	pz = bb_1.p[*i__ * 3 - 1];
	rx = aa_1.r__[*i__ * 3 - 3];
	ry = aa_1.r__[*i__ * 3 - 2];
	rz = aa_1.r__[*i__ * 3 - 1];
	dm = cc_1.e[*i__ - 1];
	r__1 = dm;
	r__2 = px;
	r__3 = py;
	r__4 = pz;
	edelta = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	pm = pc_1.epion[*nnn + *irun * 150001 - 150002];
	am = .93828f;
	if (*nlab == 2)
	{
		am = .939457f;
	}
	r__2 = dm;
	r__3 = am;
	r__4 = pm;
	r__1 = (r__2 * r__2 - r__3 * r__3 + r__4 * r__4) / (dm * 2.f);
	r__5 = pm;
	q2 = r__1 * r__1 - r__5 * r__5;
	if (q2 <= 0.f)
	{
		q2 = 1e-9f;
	}
	q = sqrt(q2);
L11:
	qx = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
	qy = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
	qz = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
	r__1 = qx;
	r__2 = qy;
	r__3 = qz;
	qs = r__1 * r__1 + r__2 * r__2 + r__3 * r__3;
	if (qs > 1.f)
	{
		goto L11;
	}
	pxp = q * qx / sqrt(qs);
	pyp = q * qy / sqrt(qs);
	pzp = q * qz / sqrt(qs);
	r__1 = q;
	r__2 = pm;
	ep = sqrt(r__1 * r__1 + r__2 * r__2);
	pxn = -pxp;
	pyn = -pyp;
	pzn = -pzp;
	r__1 = q;
	r__2 = am;
	en = sqrt(r__1 * r__1 + r__2 * r__2);
	gd = edelta / dm;
	fgd = gd / (gd + 1.f);
	bdx = px / edelta;
	bdy = py / edelta;
	bdz = pz / edelta;
	bpp = bdx * pxp + bdy * pyp + bdz * pzp;
	bpn = bdx * pxn + bdy * pyn + bdz * pzn;
	bb_1.p[*i__ * 3 - 3] = pxn + bdx * gd * (fgd * bpn + en);
	bb_1.p[*i__ * 3 - 2] = pyn + bdy * gd * (fgd * bpn + en);
	bb_1.p[*i__ * 3 - 1] = pzn + bdz * gd * (fgd * bpn + en);
	cc_1.e[*i__ - 1] = am;
	pb_1.ppion[(*nnn + *irun * 150001) * 3 - 450006] = pxp + bdx * gd * (fgd * bpp + ep);
	pb_1.ppion[(*nnn + *irun * 150001) * 3 - 450005] = pyp + bdy * gd * (fgd * bpp + ep);
	pb_1.ppion[(*nnn + *irun * 150001) * 3 - 450004] = pzp + bdz * gd * (fgd * bpp + ep);
	dpert_1.dppion[*nnn + *irun * 150001 - 150002] = dpert_1.dpertp[*i__ - 1];
	pa_1.rpion[(*nnn + *irun * 150001) * 3 - 450006] = aa_1.r__[*i__ * 3 - 3];
	pa_1.rpion[(*nnn + *irun * 150001) * 3 - 450005] = aa_1.r__[*i__ * 3 - 2];
	pa_1.rpion[(*nnn + *irun * 150001) * 3 - 450004] = aa_1.r__[*i__ * 3 - 1];

	r__1 = pc_1.epion[*nnn + *irun * 150001 - 150002];
	r__2 = pb_1.ppion[(*nnn + *irun * 150001) * 3 - 450006];
	r__3 = pb_1.ppion[(*nnn + *irun * 150001) * 3 - 450005];
	r__4 = pb_1.ppion[(*nnn + *irun * 150001) * 3 - 450004];
	r__5 = cc_1.e[*i__ - 1];
	r__6 = bb_1.p[*i__ * 3 - 3];
	r__7 = bb_1.p[*i__ * 3 - 2];
	r__8 = bb_1.p[*i__ * 3 - 1];
	devio = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4) +
			sqrt(r__5 * r__5 + r__6 * r__6 + r__7 * r__7 + r__8 * r__8) -
			leadng_1.e1;
	if (*nt == input2_1.ntmax)
	{
		tau0 = .19733f / *wid;
		taudcy = tau0 * -1.f * log(1.f - ranart_(&rndf77_1.nseed));
		taudcy = taudcy * leadng_1.e1 / leadng_1.em1;
		leadng_1.tfnl += taudcy;
		leadng_1.xfnl += leadng_1.px1 / leadng_1.e1 * taudcy;
		leadng_1.yfnl += leadng_1.py1 / leadng_1.e1 * taudcy;
		leadng_1.zfnl += leadng_1.pz1 / leadng_1.e1 * taudcy;
		aa_1.r__[*i__ * 3 - 3] = leadng_1.xfnl;
		aa_1.r__[*i__ * 3 - 2] = leadng_1.yfnl;
		aa_1.r__[*i__ * 3 - 1] = leadng_1.zfnl;
		tdecay_1.tfdcy[*i__ - 1] = leadng_1.tfnl;
		pa_1.rpion[(*nnn + *irun * 150001) * 3 - 450006] = leadng_1.xfnl;
		pa_1.rpion[(*nnn + *irun * 150001) * 3 - 450005] = leadng_1.yfnl;
		pa_1.rpion[(*nnn + *irun * 150001) * 3 - 450004] = leadng_1.zfnl;
		tdecay_1.tfdpi[*nnn + *irun * 150001 - 150002] = leadng_1.tfnl;
	}
	return 0;
}

int decay2_(int *irun, int *i__, int *nnn,
			int *iseed, float *wid, int *nt)
{
	int i__1;

	static float x3, dm;
	static int lbi, nlab;
	extern int dkine2_(int *, int *, int *,
					   int *, int *, float *, int *);
	static int lbanti;
	extern double ranart_(int *);

	lbanti = ee_1.lb[*i__ - 1];

	dm = cc_1.e[*i__ - 1];
	if ((i__1 = ee_1.lb[*i__ - 1], abs(i__1)) == 11)
	{
		x3 = ranart_(&rndf77_1.nseed);
		if (x3 < .33333333333333331f)
		{
			ee_1.lb[*i__ - 1] = 2;
			nlab = 2;
			pd_1.lpion[*nnn + *irun * 150001 - 150002] = 5;
			pc_1.epion[*nnn + *irun * 150001 - 150002] = .13957f;
			pd_1.lpion[*nnn + 1 + *irun * 150001 - 150002] = 4;
			pc_1.epion[*nnn + 1 + *irun * 150001 - 150002] = .13496f;
		}
		else if (x3 < .66666666666666663f && x3 > .33333333333333331f)
		{
			ee_1.lb[*i__ - 1] = 1;
			nlab = 1;
			pd_1.lpion[*nnn + *irun * 150001 - 150002] = 5;
			pc_1.epion[*nnn + *irun * 150001 - 150002] = .13957f;
			pd_1.lpion[*nnn + 1 + *irun * 150001 - 150002] = 3;
			pc_1.epion[*nnn + 1 + *irun * 150001 - 150002] = .13957f;
		}
		else
		{
			ee_1.lb[*i__ - 1] = 1;
			nlab = 1;
			pd_1.lpion[*nnn + *irun * 150001 - 150002] = 4;
			pc_1.epion[*nnn + *irun * 150001 - 150002] = .13496f;
			pd_1.lpion[*nnn + 1 + *irun * 150001 - 150002] = 4;
			pc_1.epion[*nnn + 1 + *irun * 150001 - 150002] = .13496f;
		}
	}
	else if ((i__1 = ee_1.lb[*i__ - 1], abs(i__1)) == 10)
	{
		x3 = ranart_(&rndf77_1.nseed);
		if (x3 < .33333333333333331f)
		{
			ee_1.lb[*i__ - 1] = 2;
			nlab = 2;
			pd_1.lpion[*nnn + *irun * 150001 - 150002] = 4;
			pc_1.epion[*nnn + *irun * 150001 - 150002] = .13496f;
			pd_1.lpion[*nnn + 1 + *irun * 150001 - 150002] = 4;
			pc_1.epion[*nnn + 1 + *irun * 150001 - 150002] = .13496f;
		}
		else if (x3 < .66666666666666663f && x3 > .33333333333333331f)
		{
			ee_1.lb[*i__ - 1] = 1;
			nlab = 1;
			pd_1.lpion[*nnn + *irun * 150001 - 150002] = 3;
			pc_1.epion[*nnn + *irun * 150001 - 150002] = .13957f;
			pd_1.lpion[*nnn + 1 + *irun * 150001 - 150002] = 4;
			pc_1.epion[*nnn + 1 + *irun * 150001 - 150002] = .13496f;
		}
		else
		{
			ee_1.lb[*i__ - 1] = 2;
			nlab = 2;
			pd_1.lpion[*nnn + *irun * 150001 - 150002] = 5;
			pc_1.epion[*nnn + *irun * 150001 - 150002] = .13957f;
			pd_1.lpion[*nnn + 1 + *irun * 150001 - 150002] = 3;
			pc_1.epion[*nnn + 1 + *irun * 150001 - 150002] = .13957f;
		}
	}
	dkine2_(irun, i__, nnn, &nlab, iseed, wid, nt);

	if (lbanti < 0)
	{
		lbi = ee_1.lb[*i__ - 1];
		if (lbi == 1 || lbi == 2)
		{
			lbi = -lbi;
		}
		else if (lbi == 3)
		{
			lbi = 5;
		}
		else if (lbi == 5)
		{
			lbi = 3;
		}
		ee_1.lb[*i__ - 1] = lbi;

		lbi = pd_1.lpion[*nnn + *irun * 150001 - 150002];
		if (lbi == 3)
		{
			lbi = 5;
		}
		else if (lbi == 5)
		{
			lbi = 3;
		}
		else if (lbi == 1 || lbi == 2)
		{
			lbi = -lbi;
		}
		pd_1.lpion[*nnn + *irun * 150001 - 150002] = lbi;

		lbi = pd_1.lpion[*nnn + 1 + *irun * 150001 - 150002];
		if (lbi == 3)
		{
			lbi = 5;
		}
		else if (lbi == 5)
		{
			lbi = 3;
		}
		else if (lbi == 1 || lbi == 2)
		{
			lbi = -lbi;
		}
		pd_1.lpion[*nnn + 1 + *irun * 150001 - 150002] = lbi;
	}

	return 0;
}

int dkine2_(int *irun, int *i__, int *nnn,
			int *nlab, int *iseed, float *wid, int *nt)
{
	float r__1, r__2, r__3, r__4, r__5, r__6, r__7, r__8, r__9, r__10, r__11,
		r__12;

	double sqrt(double), cos(double), sin(double), log(double);

	static float q, q2, gd, am, dm, en, ep, qs, px, py, pz, rx, ry, rz, qx, qy,
		qz, gd1, ep0, bp0, pm1, pm2, p1m, p2m, p3m, p1p, p2p, p3p, px0,
		py0, pz0, fai, fgd, bdx, bdy, bdz, bpp, epn, epp, bpn, css, pxn,
		pyn, pxp, sss, pyp, pzp, pzn, fgd1, bpn1, bpp1, tau0, pmax, pmax2,
		betax, betay, betaz, enucl, devio, epion1, epion2, edelta;
	extern double ranart_(int *);
	static float taudcy;

	px = bb_1.p[*i__ * 3 - 3];
	py = bb_1.p[*i__ * 3 - 2];
	pz = bb_1.p[*i__ * 3 - 1];
	rx = aa_1.r__[*i__ * 3 - 3];
	ry = aa_1.r__[*i__ * 3 - 2];
	rz = aa_1.r__[*i__ * 3 - 1];
	dm = cc_1.e[*i__ - 1];
	r__1 = dm;
	r__2 = px;
	r__3 = py;
	r__4 = pz;
	edelta = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	pm1 = pc_1.epion[*nnn + *irun * 150001 - 150002];
	pm2 = pc_1.epion[*nnn + 1 + *irun * 150001 - 150002];
	am = .939457f;
	if (*nlab == 1)
	{
		am = .93828f;
	}
	r__1 = dm;
	r__2 = am + pm1 + pm2;
	r__3 = dm;
	r__4 = am - pm1 - pm2;
	r__5 = dm;
	pmax2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / 4 / (r__5 * r__5);
	pmax = sqrt(pmax2);
	css = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
	r__1 = css;
	sss = sqrt(1 - r__1 * r__1);
	fai = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
	px0 = pmax * sss * cos(fai);
	py0 = pmax * sss * sin(fai);
	pz0 = pmax * css;
	r__1 = px0;
	r__2 = py0;
	r__3 = pz0;
	r__4 = am;
	ep0 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	betax = -px0 / (dm - ep0);
	betay = -py0 / (dm - ep0);
	betaz = -pz0 / (dm - ep0);
	r__1 = betax;
	r__2 = betay;
	r__3 = betaz;
	gd1 = 1.f / sqrt(1 - r__1 * r__1 - r__2 * r__2 - r__3 * r__3);
	fgd1 = gd1 / (gd1 + 1);
	r__1 = (dm - ep0) / (gd1 * 2.f);
	r__2 = pm1;
	q2 = r__1 * r__1 - r__2 * r__2;
	if (q2 <= 0.f)
	{
		q2 = 1e-9f;
	}
	q = sqrt(q2);
L11:
	qx = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
	qy = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
	qz = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
	r__1 = qx;
	r__2 = qy;
	r__3 = qz;
	qs = r__1 * r__1 + r__2 * r__2 + r__3 * r__3;
	if (qs > 1.f)
	{
		goto L11;
	}
	pxp = q * qx / sqrt(qs);
	pyp = q * qy / sqrt(qs);
	pzp = q * qz / sqrt(qs);
	r__1 = q;
	r__2 = pm1;
	ep = sqrt(r__1 * r__1 + r__2 * r__2);
	pxn = -pxp;
	pyn = -pyp;
	pzn = -pzp;
	r__1 = q;
	r__2 = pm2;
	en = sqrt(r__1 * r__1 + r__2 * r__2);
	bpp1 = betax * pxp + betay * pyp + betaz * pzp;
	bpn1 = betax * pxn + betay * pyn + betaz * pzn;
	p1m = pxn + betax * gd1 * (fgd1 * bpn1 + en);
	p2m = pyn + betay * gd1 * (fgd1 * bpn1 + en);
	p3m = pzn + betaz * gd1 * (fgd1 * bpn1 + en);
	r__1 = p1m;
	r__2 = p2m;
	r__3 = p3m;
	r__4 = pm2;
	epn = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	p1p = pxp + betax * gd1 * (fgd1 * bpp1 + ep);
	p2p = pyp + betay * gd1 * (fgd1 * bpp1 + ep);
	p3p = pzp + betaz * gd1 * (fgd1 * bpp1 + ep);
	r__1 = p1p;
	r__2 = p2p;
	r__3 = p3p;
	r__4 = pm1;
	epp = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	gd = edelta / dm;
	fgd = gd / (gd + 1.f);
	bdx = px / edelta;
	bdy = py / edelta;
	bdz = pz / edelta;
	bp0 = bdx * px0 + bdy * py0 + bdz * pz0;
	bpp = bdx * p1p + bdy * p2p + bdz * p3p;
	bpn = bdx * p1m + bdy * p2m + bdz * p3m;
	bb_1.p[*i__ * 3 - 3] = px0 + bdx * gd * (fgd * bp0 + ep0);
	bb_1.p[*i__ * 3 - 2] = py0 + bdy * gd * (fgd * bp0 + ep0);
	bb_1.p[*i__ * 3 - 1] = pz0 + bdz * gd * (fgd * bp0 + ep0);
	cc_1.e[*i__ - 1] = am;
	ee_1.id[*i__ - 1] = 0;
	r__1 = bb_1.p[*i__ * 3 - 3];
	r__2 = bb_1.p[*i__ * 3 - 2];
	r__3 = bb_1.p[*i__ * 3 - 1];
	r__4 = cc_1.e[*i__ - 1];
	enucl = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	pb_1.ppion[(*nnn + *irun * 150001) * 3 - 450006] = p1p + bdx * gd * (fgd * bpp + epp);
	pb_1.ppion[(*nnn + *irun * 150001) * 3 - 450005] = p2p + bdy * gd * (fgd * bpp + epp);
	pb_1.ppion[(*nnn + *irun * 150001) * 3 - 450004] = p3p + bdz * gd * (fgd * bpp + epp);
	r__1 = pb_1.ppion[(*nnn + *irun * 150001) * 3 - 450006];
	r__2 = pb_1.ppion[(*nnn + *irun * 150001) * 3 - 450005];
	r__3 = pb_1.ppion[(*nnn + *irun * 150001) * 3 - 450004];
	r__4 = pc_1.epion[*nnn + *irun * 150001 - 150002];
	epion1 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	pa_1.rpion[(*nnn + *irun * 150001) * 3 - 450006] = aa_1.r__[*i__ * 3 - 3];
	pa_1.rpion[(*nnn + *irun * 150001) * 3 - 450005] = aa_1.r__[*i__ * 3 - 2];
	pa_1.rpion[(*nnn + *irun * 150001) * 3 - 450004] = aa_1.r__[*i__ * 3 - 1];
	pb_1.ppion[(*nnn + 1 + *irun * 150001) * 3 - 450006] = p1m + bdx * gd * (fgd * bpn + epn);
	pb_1.ppion[(*nnn + 1 + *irun * 150001) * 3 - 450005] = p2m + bdy * gd * (fgd * bpn + epn);
	pb_1.ppion[(*nnn + 1 + *irun * 150001) * 3 - 450004] = p3m + bdz * gd * (fgd * bpn + epn);
	dpert_1.dppion[*nnn + *irun * 150001 - 150002] = dpert_1.dpertp[*i__ - 1];
	dpert_1.dppion[*nnn + 1 + *irun * 150001 - 150002] = dpert_1.dpertp[*i__ - 1];

	r__1 = pb_1.ppion[(*nnn + 1 + *irun * 150001) * 3 - 450006];
	r__2 = pb_1.ppion[(*nnn + 1 + *irun * 150001) * 3 - 450005];
	r__3 = pb_1.ppion[(*nnn + 1 + *irun * 150001) * 3 - 450004];
	r__4 = pc_1.epion[*nnn + 1 + *irun * 150001 - 150002];
	epion2 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	pa_1.rpion[(*nnn + 1 + *irun * 150001) * 3 - 450006] = aa_1.r__[*i__ * 3 - 3];
	pa_1.rpion[(*nnn + 1 + *irun * 150001) * 3 - 450005] = aa_1.r__[*i__ * 3 - 2];
	pa_1.rpion[(*nnn + 1 + *irun * 150001) * 3 - 450004] = aa_1.r__[*i__ * 3 - 1];

	r__1 = pc_1.epion[*nnn + *irun * 150001 - 150002];
	r__2 = pb_1.ppion[(*nnn + *irun * 150001) * 3 - 450006];
	r__3 = pb_1.ppion[(*nnn + *irun * 150001) * 3 - 450005];
	r__4 = pb_1.ppion[(*nnn + *irun * 150001) * 3 - 450004];
	r__5 = cc_1.e[*i__ - 1];
	r__6 = bb_1.p[*i__ * 3 - 3];
	r__7 = bb_1.p[*i__ * 3 - 2];
	r__8 = bb_1.p[*i__ * 3 - 1];
	r__9 = pc_1.epion[*nnn + 1 + *irun * 150001 - 150002];
	r__10 = pb_1.ppion[(*nnn + 1 + *irun * 150001) * 3 - 450006];
	r__11 = pb_1.ppion[(*nnn + 1 + *irun * 150001) * 3 - 450005];
	r__12 = pb_1.ppion[(*nnn + 1 + *irun * 150001) * 3 - 450004];
	devio = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4) +
			sqrt(r__5 * r__5 + r__6 * r__6 + r__7 * r__7 + r__8 * r__8) +
			sqrt(r__9 * r__9 + r__10 * r__10 + r__11 * r__11 + r__12 * r__12) - leadng_1.e1;
	if (*nt == input2_1.ntmax)
	{
		tau0 = .19733f / *wid;
		taudcy = tau0 * -1.f * log(1.f - ranart_(&rndf77_1.nseed));
		taudcy = taudcy * leadng_1.e1 / leadng_1.em1;
		leadng_1.tfnl += taudcy;
		leadng_1.xfnl += leadng_1.px1 / leadng_1.e1 * taudcy;
		leadng_1.yfnl += leadng_1.py1 / leadng_1.e1 * taudcy;
		leadng_1.zfnl += leadng_1.pz1 / leadng_1.e1 * taudcy;
		aa_1.r__[*i__ * 3 - 3] = leadng_1.xfnl;
		aa_1.r__[*i__ * 3 - 2] = leadng_1.yfnl;
		aa_1.r__[*i__ * 3 - 1] = leadng_1.zfnl;
		tdecay_1.tfdcy[*i__ - 1] = leadng_1.tfnl;
		pa_1.rpion[(*nnn + *irun * 150001) * 3 - 450006] = leadng_1.xfnl;
		pa_1.rpion[(*nnn + *irun * 150001) * 3 - 450005] = leadng_1.yfnl;
		pa_1.rpion[(*nnn + *irun * 150001) * 3 - 450004] = leadng_1.zfnl;
		tdecay_1.tfdpi[*nnn + *irun * 150001 - 150002] = leadng_1.tfnl;
		pa_1.rpion[(*nnn + 1 + *irun * 150001) * 3 - 450006] = leadng_1.xfnl;
		pa_1.rpion[(*nnn + 1 + *irun * 150001) * 3 - 450005] = leadng_1.yfnl;
		pa_1.rpion[(*nnn + 1 + *irun * 150001) * 3 - 450004] = leadng_1.zfnl;
		tdecay_1.tfdpi[*nnn + 1 + *irun * 150001 - 150002] = leadng_1.tfnl;
	}
	return 0;
}

int dreson_(int *i1, int *i2)
{
	int i__1, i__2, i__3, i__4;
	float r__1, r__2, r__3, r__4;

	double sqrt(double);

	static int i__;
	static float e10, e20, dm;

	r__1 = cc_1.e[*i1 - 1];
	r__2 = bb_1.p[*i1 * 3 - 3];
	r__3 = bb_1.p[*i1 * 3 - 2];
	r__4 = bb_1.p[*i1 * 3 - 1];
	e10 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	r__1 = cc_1.e[*i2 - 1];
	r__2 = bb_1.p[*i2 * 3 - 3];
	r__3 = bb_1.p[*i2 * 3 - 2];
	r__4 = bb_1.p[*i2 * 3 - 1];
	e20 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	if ((i__1 = ee_1.lb[*i2 - 1], abs(i__1)) == 1 || (i__2 = ee_1.lb[*i2 - 1], abs(i__2)) == 2 || (i__3 = ee_1.lb[*i2 - 1], abs(i__3)) >= 6 && (i__4 = ee_1.lb[*i2 - 1], abs(i__4)) <= 17)
	{
		cc_1.e[*i1 - 1] = 0.f;
		i__ = *i2;
	}
	else
	{
		cc_1.e[*i2 - 1] = 0.f;
		i__ = *i1;
	}
	bb_1.p[i__ * 3 - 3] = bb_1.p[*i1 * 3 - 3] + bb_1.p[*i2 * 3 - 3];
	bb_1.p[i__ * 3 - 2] = bb_1.p[*i1 * 3 - 2] + bb_1.p[*i2 * 3 - 2];
	bb_1.p[i__ * 3 - 1] = bb_1.p[*i1 * 3 - 1] + bb_1.p[*i2 * 3 - 1];
	r__1 = e10 + e20;
	r__2 = bb_1.p[i__ * 3 - 3];
	r__3 = bb_1.p[i__ * 3 - 2];
	r__4 = bb_1.p[i__ * 3 - 1];
	dm = sqrt(r__1 * r__1 - r__2 * r__2 - r__3 * r__3 - r__4 * r__4);
	cc_1.e[i__ - 1] = dm;
	return 0;
}

int rhores_(int *i1, int *i2)
{
	float r__1, r__2, r__3, r__4;

	double sqrt(double);

	static float e10, e20, dm;

	r__1 = cc_1.e[*i1 - 1];
	r__2 = bb_1.p[*i1 * 3 - 3];
	r__3 = bb_1.p[*i1 * 3 - 2];
	r__4 = bb_1.p[*i1 * 3 - 1];
	e10 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	r__1 = cc_1.e[*i2 - 1];
	r__2 = bb_1.p[*i2 * 3 - 3];
	r__3 = bb_1.p[*i2 * 3 - 2];
	r__4 = bb_1.p[*i2 * 3 - 1];
	e20 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	bb_1.p[*i1 * 3 - 3] += bb_1.p[*i2 * 3 - 3];
	bb_1.p[*i1 * 3 - 2] += bb_1.p[*i2 * 3 - 2];
	bb_1.p[*i1 * 3 - 1] += bb_1.p[*i2 * 3 - 1];
	r__1 = e10 + e20;
	r__2 = bb_1.p[*i1 * 3 - 3];
	r__3 = bb_1.p[*i1 * 3 - 2];
	r__4 = bb_1.p[*i1 * 3 - 1];
	dm = sqrt(r__1 * r__1 - r__2 * r__2 - r__3 * r__3 - r__4 * r__4);
	cc_1.e[*i1 - 1] = dm;
	cc_1.e[*i2 - 1] = 0.f;
	return 0;
}

double xnpi_(int *i1, int *i2, int *la, float *xmax)
{
	float ret_val, r__1, r__2, r__3, r__4, r__5;

	double sqrt(double);

	static float f1, p1, p2, p3, e10, e20, dm, gam;
	extern double w1440_(float *), w1535_(float *);
	static float avpi;
	extern double width_(float *);
	static float pdelt2, pstar2, avmass;

	avmass = .93886849999999999f;
	avpi = .13803333333333334f;
	r__1 = cc_1.e[*i1 - 1];
	r__2 = bb_1.p[*i1 * 3 - 3];
	r__3 = bb_1.p[*i1 * 3 - 2];
	r__4 = bb_1.p[*i1 * 3 - 1];
	e10 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	r__1 = cc_1.e[*i2 - 1];
	r__2 = bb_1.p[*i2 * 3 - 3];
	r__3 = bb_1.p[*i2 * 3 - 2];
	r__4 = bb_1.p[*i2 * 3 - 1];
	e20 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	p1 = bb_1.p[*i1 * 3 - 3] + bb_1.p[*i2 * 3 - 3];
	p2 = bb_1.p[*i1 * 3 - 2] + bb_1.p[*i2 * 3 - 2];
	p3 = bb_1.p[*i1 * 3 - 1] + bb_1.p[*i2 * 3 - 1];
	r__1 = e10 + e20;
	r__2 = p1;
	r__3 = p2;
	r__4 = p3;
	dm = sqrt(r__1 * r__1 - r__2 * r__2 - r__3 * r__3 - r__4 * r__4);
	if (dm <= 1.1f)
	{
		ret_val = 1e-9f;
		return ret_val;
	}
	if (*la == 1)
	{
		gam = width_(&dm);
		r__1 = gam;
		r__2 = gam;
		r__3 = dm - 1.232f;
		f1 = r__1 * r__1 * .25f / (r__2 * r__2 * .25f + r__3 * r__3);
		pdelt2 = .051622f;
		goto L10;
	}
	if (*la == 0)
	{
		gam = w1440_(&dm);
		r__1 = gam;
		r__2 = gam;
		r__3 = dm - 1.44f;
		f1 = r__1 * r__1 * .25f / (r__2 * r__2 * .25f + r__3 * r__3);
		pdelt2 = .157897f;
		goto L10;
	}
	if (*la == 2)
	{
		gam = w1535_(&dm);
		r__1 = gam;
		r__2 = gam;
		r__3 = dm - 1.535f;
		f1 = r__1 * r__1 * .25f / (r__2 * r__2 * .25f + r__3 * r__3);
		pdelt2 = .2181f;
	}
L10:
	r__2 = dm;
	r__3 = avmass;
	r__4 = avpi;
	r__1 = (r__2 * r__2 - r__3 * r__3 + r__4 * r__4) / (dm * 2.f);
	r__5 = avpi;
	pstar2 = r__1 * r__1 - r__5 * r__5;
	if (pstar2 <= 0.f)
	{
		ret_val = 1e-9f;
	}
	else
	{
		ret_val = f1 * (pdelt2 / pstar2) * *xmax / 10.f;
	}
	return ret_val;
}

double sigma_(float *srt, int *id, int *ioi, int *iof)
{
	float ret_val, r__1, r__2, r__3, r__4, r__5, r__6, r__7;
	double d__1, d__2;

	double atan(double), log(double), sqrt(double), pow_dd(double *, double *);

	static float q, s, t, p0, q0, p2, s0, t0, q2, p02, q02, pr, ss, am0, pr2,
		ss0, alfa, beta, deln, amass, zplus, amass0, zminus;

	if (*id == 1)
	{
		amass0 = 1.22f;
		t0 = .12f;
	}
	else
	{
		amass0 = 1.43f;
		t0 = .2f;
	}
	if (*ioi == 1 && *iof == 1)
	{
		alfa = 3.772f;
		beta = 1.262f;
		am0 = 1.188f;
		t = .09902f;
	}
	if (*ioi == 1 && *iof == 0)
	{
		alfa = 15.28f;
		beta = 0.f;
		am0 = 1.245f;
		t = .1374f;
	}
	if (*ioi == 0 && *iof == 1)
	{
		alfa = 146.3f;
		beta = 0.f;
		am0 = 1.472f;
		t = .02649f;
	}
	zplus = (*srt - .9383f - amass0) * 2.f / t0;
	zminus = (1.0767f - amass0) * 2.f / t0;
	deln = atan(zplus) - atan(zminus);
	if (deln == 0.f)
	{
		deln = 1e-6f;
	}
	r__1 = zplus;
	r__2 = zminus;
	amass = amass0 + t0 / 4.f * log((r__1 * r__1 + 1.f) / (r__2 * r__2 + 1.f)) / deln;
	r__1 = *srt;
	s = r__1 * r__1;
	p2 = s / 4.f - .88040689000000005f;
	r__1 = am0 + .9383f;
	s0 = r__1 * r__1;
	p02 = s0 / 4.f - .88040689000000005f;
	p0 = sqrt(p02);
	r__1 = .9383f - amass;
	r__2 = amass + .9383f;
	pr2 = (s - r__1 * r__1) * (s - r__2 * r__2) / (s * 4.f);
	if (pr2 > 1e-6f)
	{
		pr = sqrt(pr2);
	}
	else
	{
		pr = 0.f;
		ret_val = 1e-6f;
		return ret_val;
	}
	r__1 = amass;
	ss = r__1 * r__1;
	q2 = (ss - .63984001000000013f) * (ss - 1.1592828900000001f) / (ss * 4.f);
	if (q2 > 1e-6f)
	{
		q = sqrt(q2);
	}
	else
	{
		q = 0.f;
		ret_val = 1e-6f;
		return ret_val;
	}
	r__1 = am0;
	ss0 = r__1 * r__1;
	q02 = (ss0 - .63984001000000013f) * (ss0 - 1.1592828900000001f) / (ss0 * 4.f);
	q0 = sqrt(q02);
	d__1 = (double)(pr / p0);
	d__2 = (double)beta;
	r__1 = am0;
	r__2 = t;
	r__3 = q / q0;
	r__5 = am0;
	r__4 = ss - r__5 * r__5;
	r__6 = am0;
	r__7 = t;
	ret_val = .12233087920268615f / (p2 * 2.f) * alfa * pow_dd(&d__1, &d__2) *
			  (r__1 * r__1) * (r__2 * r__2) * (r__3 * (r__3 * r__3)) / (r__4 * r__4 + r__6 * r__6 * (r__7 * r__7));
	ret_val *= 10.f;
	if (ret_val == 0.f)
	{
		ret_val = 1e-6f;
	}
	return ret_val;
}

double denom_(float *srt, float *con)
{
	int i__1;
	float ret_val, r__1, r__2, r__3, r__4, r__5, r__6;

	double sqrt(double);

	static float f;
	static int i__;
	static float q, s, a1, p0, p1, q2, dm, tq, am0, amn, amp, sum, amin, amax,
		avpi;
	static int nmax;
	static float dmass;

	avpi = .13803333333333334f;
	am0 = 1.232f;
	amn = .9383f;
	amp = avpi;
	amax = *srt - .9383f;
	amin = avpi + .9383f;
	nmax = 200;
	dmass = (amax - amin) / (float)nmax;
	sum = 0.f;
	i__1 = nmax + 1;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		dm = amin + (float)(i__ - 1) * dmass;
		if (*con == 1.f)
		{
			r__2 = dm;
			r__3 = amn;
			r__4 = amp;
			r__1 = (r__2 * r__2 - r__3 * r__3 + r__4 * r__4) / (dm * 2.f);
			r__5 = amp;
			q2 = r__1 * r__1 - r__5 * r__5;
			if (q2 > 0.f)
			{
				q = sqrt(q2);
			}
			else
			{
				q = 1e-6f;
			}
			r__1 = q;
			r__2 = amp;
			r__3 = q / amp;
			tq = r__1 * (r__1 * r__1) * .47f / (r__2 * r__2 * (r__3 * r__3 * .6f + 1.f));
		}
		else if (*con == 2.f)
		{
			tq = .2f;
			am0 = 1.44f;
		}
		else if (*con == -1.f)
		{
			tq = .1f;
			am0 = 1.535f;
		}
		r__1 = am0;
		r__2 = am0;
		r__3 = tq;
		r__5 = dm;
		r__6 = am0;
		r__4 = r__5 * r__5 - r__6 * r__6;
		a1 = tq * 4.f * (r__1 * r__1) / (r__2 * r__2 * (r__3 * r__3) + r__4 * r__4);
		r__1 = *srt;
		s = r__1 * r__1;
		r__2 = dm;
		r__3 = amn;
		r__1 = s + r__2 * r__2 - r__3 * r__3;
		r__4 = dm;
		p0 = r__1 * r__1 / (s * 4.f) - r__4 * r__4;
		if (p0 <= 0.f)
		{
			p1 = 1e-6f;
		}
		else
		{
			p1 = sqrt(p0);
		}
		f = dm * a1 * p1;
		if (i__ == 1 || i__ == nmax + 1)
		{
			sum += f * .5f;
		}
		else
		{
			sum += f;
		}
	}
	ret_val = sum * dmass / 6.2831852000000001f;
	return ret_val;
}

double ang_(float *srt, int *iseed)
{
	float ret_val, r__1, r__2;
	double d__1;

	double sqrt(double), pow_dd(double *, double *);

	static float p, q, x, b1s, b2s, ang1, ang2;
	extern double ranart_(int *);

	if (*srt > 2.14f && *srt <= 2.4f)
	{
		r__1 = *srt;
		b1s = 29.03f - *srt * 23.75f + r__1 * r__1 * 4.865f;
		r__1 = *srt;
		b2s = *srt * 25.53f - 30.33f - r__1 * r__1 * 5.301f;
	}
	if (*srt > 2.4f)
	{
		b1s = .06f;
		b2s = .4f;
	}
	x = ranart_(&rndf77_1.nseed);
	p = b1s / b2s;
	q = (x * 2.f - 1.f) * (b1s + b2s) / b2s;
	r__1 = q / 2.f;
	r__2 = p / 3.f;
	if (-q / 2.f + sqrt(r__1 * r__1 + r__2 * (r__2 * r__2)) >= 0.f)
	{
		r__1 = q / 2.f;
		r__2 = p / 3.f;
		d__1 = (double)(-q / 2.f + sqrt(r__1 * r__1 + r__2 * (r__2 *
																  r__2)));
		ang1 = pow_dd(&d__1, &c_b5);
	}
	else
	{
		r__1 = q / 2.f;
		r__2 = p / 3.f;
		d__1 = (double)(q / 2.f - sqrt(r__1 * r__1 + r__2 * (r__2 * r__2)));
		ang1 = -pow_dd(&d__1, &c_b5);
	}
	r__1 = q / 2.f;
	r__2 = p / 3.f;
	if (-q / 2.f - sqrt(r__1 * r__1 + r__2 * (r__2 * r__2)) >= 0.f)
	{
		r__1 = q / 2.f;
		r__2 = p / 3.f;
		d__1 = (double)(-q / 2.f - sqrt(r__1 * r__1 + r__2 * (r__2 *
																  r__2)));
		ang2 = pow_dd(&d__1, &c_b5);
	}
	else
	{
		r__1 = q / 2.f;
		r__2 = p / 3.f;
		d__1 = (double)(q / 2.f + sqrt(r__1 * r__1 + r__2 * (r__2 * r__2)));
		ang2 = -pow_dd(&d__1, &c_b5);
	}
	ret_val = ang1 + ang2;
	return ret_val;
}

double pnlka_(float *srt)
{
	float ret_val;

	static float t1, aka, ala, ana, sbbk;

	ala = 1.116f;
	aka = .498f;
	ana = .939f;
	t1 = ala + aka;
	if (*srt <= t1)
	{
		ret_val = 0.f;
	}
	else
	{
		if (*srt < 1.7f)
		{
			sbbk = (*srt - t1) * 9.8901098901098905f;
		}
		if (*srt >= 1.7f)
		{
			sbbk = .09f / (*srt - 1.6f);
		}
		ret_val = sbbk * .25f;
		ret_val /= 10.f;
	}
	return ret_val;
}

double pnska_(float *srt)
{
	float ret_val;

	static float t1, aka, ala, ana, asa, sbb1, sbb2;

	if (*srt > 3.f)
	{
		ret_val = 0.f;
		return ret_val;
	}
	ala = 1.116f;
	aka = .498f;
	ana = .939f;
	asa = 1.197f;
	t1 = asa + aka;
	if (*srt <= t1)
	{
		ret_val = 0.f;
		return ret_val;
	}
	if (*srt < 1.9f)
	{
		sbb1 = (*srt - t1) * 3.2110091743119265f;
	}
	if (*srt >= 1.9f)
	{
		sbb1 = .14f / (*srt - 1.7f);
	}
	sbb2 = 0.f;
	if (*srt > 1.682f)
	{
		sbb2 = (1.f - (*srt - 1.682f) * .75f) * .5f;
	}
	ret_val = (sbb1 + sbb2) * .25f;
	ret_val /= 10.f;
	return ret_val;
}

double fkaon_(float *p, float *pmax)
{
	float ret_val, r__1;

	static float fmax;

	fmax = .148f;
	if (*pmax == 0.f)
	{
		*pmax = 1e-6f;
	}
	r__1 = *p / *pmax;
	ret_val = (1.f - *p / *pmax) * (r__1 * r__1);
	if (ret_val > fmax)
	{
		ret_val = fmax;
	}
	ret_val /= fmax;
	return ret_val;
}

int m1535_(int *lb1, int *lb2, float *srt, float *x1535)
{
	float r__1;

	static float s0, sigma;

	s0 = 2.424f;
	*x1535 = 0.f;
	if (*srt <= s0)
	{
		return 0;
	}
	r__1 = *srt - s0;
	sigma = (*srt - s0) * .20399999999999999f / (r__1 * r__1 + .058f);
	if (*lb1 * *lb2 == 18 && (*lb1 == 2 || *lb2 == 2) || *lb1 * *lb2 == 6 && (*lb1 == 1 || *lb2 == 1) || *lb1 * *lb2 == 8 && (*lb1 == 1 || *lb2 == 1))
	{
		*x1535 = sigma;
		return 0;
	}
	if (*lb1 * *lb2 == 7)
	{
		*x1535 = sigma * 3.f;
		return 0;
	}
	if (*lb1 * *lb2 == 11 || *lb1 * *lb2 == 20 && (*lb1 == 2 || *lb2 == 2))
	{
		*x1535 = sigma;
		return 0;
	}
	if (*lb1 * *lb2 == 10 && (*lb1 == 1 || *lb2 == 1) || *lb1 * *lb2 == 22 &&
															 (*lb1 == 2 || *lb2 == 2))
	{
		*x1535 = sigma * 3.f;
	}
	return 0;
}

int n1535_(int *lb1, int *lb2, float *srt, float *x1535)
{
	float r__1;

	static float s0, sigma;

	s0 = 2.424f;
	*x1535 = 0.f;
	if (*srt <= s0)
	{
		return 0;
	}
	r__1 = *srt - s0;
	sigma = (*srt - s0) * .20399999999999999f / (r__1 * r__1 + .058f);
	if (*lb1 * *lb2 == 1 || *lb1 == 2 && *lb2 == 2)
	{
		*x1535 = sigma;
		return 0;
	}
	if (*lb1 * *lb2 == 2)
	{
		*x1535 = sigma * 3.f;
		return 0;
	}
	if (*lb1 * *lb2 == 63 && (*lb1 == 7 || *lb2 == 7) || *lb1 * *lb2 == 64 && (*lb1 == 8 || *lb2 == 8) || *lb1 * *lb2 == 48 && (*lb1 == 6 || *lb2 == 6) || *lb1 * *lb2 == 49 && (*lb1 == 7 || *lb2 == 7))
	{
		*x1535 = sigma;
		return 0;
	}
	if (*lb1 * *lb2 == 54 && (*lb1 == 6 || *lb2 == 6) || *lb1 * *lb2 == 56 &&
															 (*lb1 == 7 || *lb2 == 7))
	{
		*x1535 = sigma * 3.f;
		return 0;
	}
	if (*lb1 == 10 && *lb2 == 10 || *lb1 == 11 && *lb2 == 11)
	{
		*x1535 = sigma;
	}
	if (*lb1 * *lb2 == 110 && (*lb1 == 10 || *lb2 == 10))
	{
		*x1535 = sigma * 3.f;
	}
	return 0;
}

int wida1_(float *dmass, float *rhomp, float *wa1, int *iseed)
{
	float r__1, r__2, r__3, r__4;

	double sqrt(double);

	static float epi, qqp, qqp2, erho, coupa, epirho;
	extern double rhomas_(float *, int *);
	static float pimass, rhomax;
	static int icount;

	pimass = .137265f;
	coupa = 14.8f;

	rhomax = *dmass - pimass - .02f;
	if (rhomax <= 0.f)
	{
		*rhomp = 0.f;
		*wa1 = -10.f;
	}
	icount = 0;
L711:
	*rhomp = rhomas_(&rhomax, iseed);
	++icount;
	if (*dmass <= pimass + *rhomp)
	{
		if (icount <= 100)
		{
			goto L711;
		}
		else
		{
			*rhomp = 0.f;
			*wa1 = -10.f;
			return 0;
		}
	}
	r__1 = *dmass;
	r__2 = *rhomp + pimass;
	r__3 = *dmass;
	r__4 = *rhomp - pimass;
	qqp2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4);
	qqp = sqrt(qqp2) / (*dmass * 2.f);
	r__1 = pimass;
	r__2 = qqp;
	epi = sqrt(r__1 * r__1 + r__2 * r__2);
	r__1 = *rhomp;
	r__2 = qqp;
	erho = sqrt(r__1 * r__1 + r__2 * r__2);
	r__2 = qqp;
	r__1 = epi * erho + r__2 * r__2;
	r__3 = *rhomp;
	r__4 = epi;
	epirho = r__1 * r__1 * 2.f + r__3 * r__3 * (r__4 * r__4);
	r__1 = coupa;
	r__2 = *dmass;
	*wa1 = r__1 * r__1 * qqp * epirho / (r__2 * r__2 * 75.398399999999995f);
	return 0;
}

double w1535_(float *dmass)
{
	float ret_val, r__1, r__2, r__3, r__4, r__5;

	double sqrt(double);

	static float aux, qavail, avmass, pimass;

	avmass = .938868f;
	pimass = .137265f;
	r__2 = *dmass;
	r__3 = avmass;
	r__4 = pimass;
	r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
	r__5 = avmass * pimass;
	aux = r__1 * r__1 * .25f - r__5 * r__5;
	if (aux > 0.f)
	{
		r__1 = *dmass;
		qavail = sqrt(aux / (r__1 * r__1));
	}
	else
	{
		qavail = 1e-6f;
	}
	ret_val = qavail * .15f / .467f;
	return ret_val;
}

double w1440_(float *dmass)
{
	float ret_val, r__1, r__2, r__3, r__4, r__5;

	double sqrt(double);

	static float aux, qavail, avmass, pimass;

	avmass = .938868f;
	pimass = .137265f;
	r__2 = *dmass;
	r__3 = avmass;
	r__4 = pimass;
	r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
	r__5 = avmass * pimass;
	aux = r__1 * r__1 * .25f - r__5 * r__5;
	if (aux > 0.f)
	{
		qavail = sqrt(aux) / *dmass;
	}
	else
	{
		qavail = 1e-6f;
	}
	r__1 = qavail / .397f;
	ret_val = r__1 * (r__1 * r__1) * .2f;
	return ret_val;
}

double xn1535_(int *i1, int *i2, int *la)
{
	float ret_val, r__1, r__2, r__3, r__4;

	double sqrt(double);

	static float f1, p1, p2, p3, e10, e20, dm, gam;
	extern double w1535_(float *);
	static float gam0, avpi, xmax, avmass;

	avmass = .93886849999999999f;
	avpi = .13803333333333334f;
	r__1 = cc_1.e[*i1 - 1];
	r__2 = bb_1.p[*i1 * 3 - 3];
	r__3 = bb_1.p[*i1 * 3 - 2];
	r__4 = bb_1.p[*i1 * 3 - 1];
	e10 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	r__1 = cc_1.e[*i2 - 1];
	r__2 = bb_1.p[*i2 * 3 - 3];
	r__3 = bb_1.p[*i2 * 3 - 2];
	r__4 = bb_1.p[*i2 * 3 - 1];
	e20 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	p1 = bb_1.p[*i1 * 3 - 3] + bb_1.p[*i2 * 3 - 3];
	p2 = bb_1.p[*i1 * 3 - 2] + bb_1.p[*i2 * 3 - 2];
	p3 = bb_1.p[*i1 * 3 - 1] + bb_1.p[*i2 * 3 - 1];
	r__1 = e10 + e20;
	r__2 = p1;
	r__3 = p2;
	r__4 = p3;
	dm = sqrt(r__1 * r__1 - r__2 * r__2 - r__3 * r__3 - r__4 * r__4);
	if (dm <= 1.1f)
	{
		ret_val = 1e-6f;
		return ret_val;
	}
	gam = w1535_(&dm);
	gam0 = .15f;
	r__1 = gam0;
	r__2 = gam;
	r__3 = dm - 1.535f;
	f1 = r__1 * r__1 * .25f / (r__2 * r__2 * .25f + r__3 * r__3);
	if (*la == 1)
	{
		xmax = 11.3f;
	}
	else
	{
		xmax = 74.f;
	}
	ret_val = f1 * xmax / 10.f;
	return ret_val;
}

double fdelta_(float *dmass)
{
	float ret_val, r__1, r__2, r__3;

	static float fd, am0, amn, avpi;
	extern double width_(float *);

	amn = .938869f;
	avpi = .13803333f;
	am0 = 1.232f;
	r__1 = width_(dmass);
	r__2 = *dmass - 1.232f;
	r__3 = width_(dmass);
	fd = r__1 * r__1 * .25f / (r__2 * r__2 + r__3 * r__3 * .25f);
	ret_val = fd;
	return ret_val;
}

double width_(float *dmass)
{
	float ret_val, r__1, r__2, r__3, r__4, r__5;

	double sqrt(double);

	static float aux, qavail, avmass, pimass;

	avmass = .938868f;
	pimass = .137265f;
	r__2 = *dmass;
	r__3 = avmass;
	r__4 = pimass;
	r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
	r__5 = avmass * pimass;
	aux = r__1 * r__1 * .25f - r__5 * r__5;
	if (aux > 0.f)
	{
		r__1 = *dmass;
		qavail = sqrt(aux / (r__1 * r__1));
	}
	else
	{
		qavail = 1e-6f;
	}
	r__1 = qavail;
	r__2 = pimass;
	r__3 = qavail / pimass;
	ret_val = r__1 * (r__1 * r__1) * .47f / (r__2 * r__2 * (r__3 * r__3 * .6f + 1.f));
	return ret_val;
}

int ddp2_(float *srt, int *iseed, float *px, float *py,
		  float *pz, float *dm1, float *pnx, float *pny, float *pnz, float *dm2, float *ppx, float *ppy, float *ppz, int *icou1)
{
	float r__1, r__2, r__3, r__4, r__5, r__6;

	double sqrt(double), cos(double), sin(double);

	static float v, w, x, f00, ga, ek, en, ep, x00, pi, bx, by, bz, pn, pt,
		pn2, fai, amn, amp, eln, sig, pnt;
	extern double ptr_(float *, int *);
	static float srt1, fain, elnc, xmax, xptr, fmax00, pbeta, ptmax, pzmax;
	static int ntrym, ntryx;
	static float trans0, ptmax2, pzmax2;
	extern int rmasdd_(float *, float *, float *, float *, float *, int *, int *, float *, float *);
	extern double ranart_(int *);
	static float xratio;

	*icou1 = 0;
	pi = 3.1415926f;
	amn = .93892500000000001f;
	amp = .137265f;
	srt1 = *srt - amp - .02f;
	ntrym = 0;
L8:
	rmasdd_(&srt1, &c_b185, &c_b185, &c_b187, &c_b187, iseed, &c__1, dm1, dm2);
	++ntrym;
	v = .43f;
	w = -.84f;
	r__1 = *srt;
	r__2 = *dm1 + *dm2 + amp;
	r__3 = *srt;
	r__4 = *dm1 - amp - *dm2;
	r__5 = *srt;
	ptmax2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / 4.f /
			 (r__5 * r__5);
	if (ptmax2 <= 0.f)
	{
		goto L8;
	}
	ptmax = sqrt(ptmax2) * 1.f / 3.f;
L7:
	pt = ptr_(&ptmax, iseed);
	r__1 = *srt;
	r__2 = *dm1 + *dm2 + amp;
	r__3 = *srt;
	r__4 = *dm1 - amp - *dm2;
	r__5 = *srt;
	r__6 = pt;
	pzmax2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / 4.f /
				 (r__5 * r__5) -
			 r__6 * r__6;
	if (pzmax2 < 0.f && ntrym <= 100)
	{
		goto L7;
	}
	else
	{
		pzmax2 = 1e-9f;
	}
	pzmax = sqrt(pzmax2);
	xmax = pzmax * 2.f / *srt;
	ntryx = 0;
	fmax00 = 1.056f;
	x00 = .26f;
	if (dabs(xmax) > .26f)
	{
		f00 = fmax00;
	}
	else
	{
		r__1 = xmax;
		f00 = v * dabs(xmax) + 1.f + w * (r__1 * r__1);
	}
L9:
	x = xmax * (1.f - ranart_(&rndf77_1.nseed) * 2.f);
	++ntryx;
	r__1 = x;
	xratio = (v * dabs(x) + 1.f + w * (r__1 * r__1)) / f00;
	if (xratio < ranart_(&rndf77_1.nseed) && ntryx <= 50)
	{
		goto L9;
	}
	*pz = *srt * .5f * x;
	fai = pi * 2.f * ranart_(&rndf77_1.nseed);
	*px = pt * cos(fai);
	*py = pt * sin(fai);
	r__1 = *dm1;
	r__2 = pt;
	r__3 = *pz;
	ek = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
	eln = *srt - ek;
	if (eln <= 0.f)
	{
		*icou1 = -1;
		return 0;
	}
	bx = -(*px) / eln;
	by = -(*py) / eln;
	bz = -(*pz) / eln;
	r__1 = bx;
	r__2 = by;
	r__3 = bz;
	ga = 1.f / sqrt(1.f - r__1 * r__1 - r__2 * r__2 - r__3 * r__3);
	elnc = eln / ga;
	r__2 = elnc;
	r__3 = *dm2;
	r__4 = amp;
	r__1 = (r__2 * r__2 + r__3 * r__3 - r__4 * r__4) / (elnc * 2.f);
	r__5 = *dm2;
	pn2 = r__1 * r__1 - r__5 * r__5;
	if (pn2 <= 0.f)
	{
		*icou1 = -1;
		return 0;
	}
	pn = sqrt(pn2);
	xptr = pn * .33f;
	pnt = ptr_(&xptr, iseed);
	fain = pi * 2.f * ranart_(&rndf77_1.nseed);
	*pnx = pnt * cos(fain);
	*pny = pnt * sin(fain);
	sig = 1.f;
	if (x > 0.f)
	{
		sig = -1.f;
	}
	r__1 = pn;
	r__2 = pnt;
	*pnz = sig * sqrt(r__1 * r__1 - r__2 * r__2);
	r__1 = *dm2;
	r__2 = *pnx;
	r__3 = *pny;
	r__4 = *pnz;
	en = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	*ppx = -(*pnx);
	*ppy = -(*pny);
	*ppz = -(*pnz);
	r__1 = amp;
	r__2 = *ppx;
	r__3 = *ppy;
	r__4 = *ppz;
	ep = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	pbeta = *pnx * bx + *pny * by + *pnz * bz;
	trans0 = ga * (ga * pbeta / (ga + 1.f) + en);
	*pnx = bx * trans0 + *pnx;
	*pny = by * trans0 + *pny;
	*pnz = bz * trans0 + *pnz;
	if (ep == 0.f)
	{
		ep = 1e-9f;
	}
	pbeta = *ppx * bx + *ppy * by + *ppz * bz;
	trans0 = ga * (ga * pbeta / (ga + 1.f) + ep);
	*ppx = bx * trans0 + *ppx;
	*ppy = by * trans0 + *ppy;
	*ppz = bz * trans0 + *ppz;
	return 0;
}

int ddrho_(float *srt, int *iseed, float *px, float *py,
		   float *pz, float *dm1, float *pnx, float *pny, float *pnz, float *dm2, float *ppx, float *ppy, float *ppz, float *amp, int *icou1)
{
	float r__1, r__2, r__3, r__4, r__5, r__6;

	double sqrt(double), cos(double), sin(double);

	static float v, w, x, f00, ga, ek, en, ep, x00, pi, bx, by, bz, pn, pt,
		pn2, fai, amn, eln, sig, pnt;
	extern double ptr_(float *, int *);
	static float srt1, fain, elnc, xmax, xptr, fmax00, pbeta, ptmax, pzmax;
	static int ntrym, ntryx;
	static float trans0, ptmax2, pzmax2;
	extern int rmasdd_(float *, float *, float *, float *, float *, int *, int *, float *, float *);
	extern double ranart_(int *), rhomas_(float *, int *);
	static float rhomax, xratio;

	*icou1 = 0;
	pi = 3.1415926f;
	amn = .93892500000000001f;
	*amp = .77000000000000002f;
	srt1 = *srt - *amp - .02f;
	ntrym = 0;
L8:
	rmasdd_(&srt1, &c_b185, &c_b185, &c_b187, &c_b187, iseed, &c__1, dm1, dm2);
	++ntrym;
	rhomax = *srt - *dm1 - *dm2 - .02f;
	if (rhomax <= 0.f && ntrym <= 20)
	{
		goto L8;
	}
	*amp = rhomas_(&rhomax, iseed);
	v = .43f;
	w = -.84f;
	r__1 = *srt;
	r__2 = *dm1 + *dm2 + *amp;
	r__3 = *srt;
	r__4 = *dm1 - *amp - *dm2;
	r__5 = *srt;
	ptmax2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / 4.f /
			 (r__5 * r__5);
	ptmax = sqrt(ptmax2) * 1.f / 3.f;
L7:
	pt = ptr_(&ptmax, iseed);
	r__1 = *srt;
	r__2 = *dm1 + *dm2 + *amp;
	r__3 = *srt;
	r__4 = *dm1 - *amp - *dm2;
	r__5 = *srt;
	r__6 = pt;
	pzmax2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / 4.f /
				 (r__5 * r__5) -
			 r__6 * r__6;
	if (pzmax2 < 0.f && ntrym <= 100)
	{
		goto L7;
	}
	else
	{
		pzmax2 = 1e-6f;
	}
	pzmax = sqrt(pzmax2);
	xmax = pzmax * 2.f / *srt;
	ntryx = 0;
	fmax00 = 1.056f;
	x00 = .26f;
	if (dabs(xmax) > .26f)
	{
		f00 = fmax00;
	}
	else
	{
		r__1 = xmax;
		f00 = v * dabs(xmax) + 1.f + w * (r__1 * r__1);
	}
L9:
	x = xmax * (1.f - ranart_(&rndf77_1.nseed) * 2.f);
	++ntryx;
	r__1 = x;
	xratio = (v * dabs(x) + 1.f + w * (r__1 * r__1)) / f00;
	if (xratio < ranart_(&rndf77_1.nseed) && ntryx <= 50)
	{
		goto L9;
	}
	*pz = *srt * .5f * x;
	fai = pi * 2.f * ranart_(&rndf77_1.nseed);
	*px = pt * cos(fai);
	*py = pt * sin(fai);
	r__1 = *dm1;
	r__2 = pt;
	r__3 = *pz;
	ek = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
	eln = *srt - ek;
	if (eln <= 0.f)
	{
		*icou1 = -1;
		return 0;
	}
	bx = -(*px) / eln;
	by = -(*py) / eln;
	bz = -(*pz) / eln;
	r__1 = bx;
	r__2 = by;
	r__3 = bz;
	ga = 1.f / sqrt(1.f - r__1 * r__1 - r__2 * r__2 - r__3 * r__3);
	elnc = eln / ga;
	r__2 = elnc;
	r__3 = *dm2;
	r__4 = *amp;
	r__1 = (r__2 * r__2 + r__3 * r__3 - r__4 * r__4) / (elnc * 2.f);
	r__5 = *dm2;
	pn2 = r__1 * r__1 - r__5 * r__5;
	if (pn2 <= 0.f)
	{
		*icou1 = -1;
		return 0;
	}
	pn = sqrt(pn2);
	xptr = pn * .33f;
	pnt = ptr_(&xptr, iseed);
	fain = pi * 2.f * ranart_(&rndf77_1.nseed);
	*pnx = pnt * cos(fain);
	*pny = pnt * sin(fain);
	sig = 1.f;
	if (x > 0.f)
	{
		sig = -1.f;
	}
	r__1 = pn;
	r__2 = pnt;
	*pnz = sig * sqrt(r__1 * r__1 - r__2 * r__2);
	r__1 = *dm2;
	r__2 = *pnx;
	r__3 = *pny;
	r__4 = *pnz;
	en = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	*ppx = -(*pnx);
	*ppy = -(*pny);
	*ppz = -(*pnz);
	r__1 = *amp;
	r__2 = *ppx;
	r__3 = *ppy;
	r__4 = *ppz;
	ep = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	pbeta = *pnx * bx + *pny * by + *pnz * bz;
	trans0 = ga * (ga * pbeta / (ga + 1.f) + en);
	*pnx = bx * trans0 + *pnx;
	*pny = by * trans0 + *pny;
	*pnz = bz * trans0 + *pnz;
	if (ep == 0.f)
	{
		ep = 1e-9f;
	}
	pbeta = *ppx * bx + *ppy * by + *ppz * bz;
	trans0 = ga * (ga * pbeta / (ga + 1.f) + ep);
	*ppx = bx * trans0 + *ppx;
	*ppy = by * trans0 + *ppy;
	*ppz = bz * trans0 + *ppz;
	return 0;
}

int pprho_(float *srt, int *iseed, float *px, float *py,
		   float *pz, float *dm1, float *pnx, float *pny, float *pnz, float *dm2, float *ppx, float *ppy, float *ppz, float *amp, int *icou1)
{
	float r__1, r__2, r__3, r__4, r__5, r__6;

	double sqrt(double), cos(double), sin(double);

	static float v, w, x, f00, ga, ek, en, ep, x00, pi, bx, by, bz, pn, pt,
		pn2, fai, amn, eln, sig, pnt;
	extern double ptr_(float *, int *);
	static float fain, elnc;
	static int icou;
	static float xmax, xptr, fmax00, pbeta, ptmax, pzmax;
	static int ntrym, ntryx;
	static float trans0, ptmax2, pzmax2;
	extern double ranart_(int *), rhomas_(float *, int *);
	static float rhomax, xratio;

	ntrym = 0;
	*icou1 = 0;
	pi = 3.1415926f;
	amn = .93892500000000001f;
	*dm1 = amn;
	*dm2 = amn;
	rhomax = *srt - *dm1 - *dm2 - .02f;
	if (rhomax <= 0.f)
	{
		icou = -1;
		return 0;
	}
	*amp = rhomas_(&rhomax, iseed);
	v = .43f;
	w = -.84f;
	r__1 = *srt;
	r__2 = *dm1 + *dm2 + *amp;
	r__3 = *srt;
	r__4 = *dm1 - *amp - *dm2;
	r__5 = *srt;
	ptmax2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / 4.f /
			 (r__5 * r__5);
	ptmax = sqrt(ptmax2) * 1.f / 3.f;
L7:
	pt = ptr_(&ptmax, iseed);
	r__1 = *srt;
	r__2 = *dm1 + *dm2 + *amp;
	r__3 = *srt;
	r__4 = *dm1 - *amp - *dm2;
	r__5 = *srt;
	r__6 = pt;
	pzmax2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / 4.f /
				 (r__5 * r__5) -
			 r__6 * r__6;
	++ntrym;
	if (pzmax2 < 0.f && ntrym <= 100)
	{
		goto L7;
	}
	else
	{
		pzmax2 = 1e-6f;
	}
	pzmax = sqrt(pzmax2);
	xmax = pzmax * 2.f / *srt;
	ntryx = 0;
	fmax00 = 1.056f;
	x00 = .26f;
	if (dabs(xmax) > .26f)
	{
		f00 = fmax00;
	}
	else
	{
		r__1 = xmax;
		f00 = v * dabs(xmax) + 1.f + w * (r__1 * r__1);
	}
L9:
	x = xmax * (1.f - ranart_(&rndf77_1.nseed) * 2.f);
	++ntryx;
	r__1 = x;
	xratio = (v * dabs(x) + 1.f + w * (r__1 * r__1)) / f00;
	if (xratio < ranart_(&rndf77_1.nseed) && ntryx <= 50)
	{
		goto L9;
	}
	*pz = *srt * .5f * x;
	fai = pi * 2.f * ranart_(&rndf77_1.nseed);
	*px = pt * cos(fai);
	*py = pt * sin(fai);
	r__1 = *dm1;
	r__2 = pt;
	r__3 = *pz;
	ek = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
	eln = *srt - ek;
	if (eln <= 0.f)
	{
		*icou1 = -1;
		return 0;
	}
	bx = -(*px) / eln;
	by = -(*py) / eln;
	bz = -(*pz) / eln;
	r__1 = bx;
	r__2 = by;
	r__3 = bz;
	ga = 1.f / sqrt(1.f - r__1 * r__1 - r__2 * r__2 - r__3 * r__3);
	elnc = eln / ga;
	r__2 = elnc;
	r__3 = *dm2;
	r__4 = *amp;
	r__1 = (r__2 * r__2 + r__3 * r__3 - r__4 * r__4) / (elnc * 2.f);
	r__5 = *dm2;
	pn2 = r__1 * r__1 - r__5 * r__5;
	if (pn2 <= 0.f)
	{
		*icou1 = -1;
		return 0;
	}
	pn = sqrt(pn2);
	xptr = pn * .33f;
	pnt = ptr_(&xptr, iseed);
	fain = pi * 2.f * ranart_(&rndf77_1.nseed);
	*pnx = pnt * cos(fain);
	*pny = pnt * sin(fain);
	sig = 1.f;
	if (x > 0.f)
	{
		sig = -1.f;
	}
	r__1 = pn;
	r__2 = pnt;
	*pnz = sig * sqrt(r__1 * r__1 - r__2 * r__2);
	r__1 = *dm2;
	r__2 = *pnx;
	r__3 = *pny;
	r__4 = *pnz;
	en = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	*ppx = -(*pnx);
	*ppy = -(*pny);
	*ppz = -(*pnz);
	r__1 = *amp;
	r__2 = *ppx;
	r__3 = *ppy;
	r__4 = *ppz;
	ep = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	pbeta = *pnx * bx + *pny * by + *pnz * bz;
	trans0 = ga * (ga * pbeta / (ga + 1.f) + en);
	*pnx = bx * trans0 + *pnx;
	*pny = by * trans0 + *pny;
	*pnz = bz * trans0 + *pnz;
	if (ep == 0.f)
	{
		ep = 1e-9f;
	}
	pbeta = *ppx * bx + *ppy * by + *ppz * bz;
	trans0 = ga * (ga * pbeta / (ga + 1.f) + ep);
	*ppx = bx * trans0 + *ppx;
	*ppy = by * trans0 + *ppy;
	*ppz = bz * trans0 + *ppz;
	return 0;
}

int ppomga_(float *srt, int *iseed, float *px, float *py,
			float *pz, float *dm1, float *pnx, float *pny, float *pnz, float *dm2, float *ppx, float *ppy, float *ppz, int *icou1)
{
	float r__1, r__2, r__3, r__4, r__5, r__6;

	double sqrt(double), cos(double), sin(double);

	static float v, w, x, f00, ga, ek, en, ep, x00, pi, bx, by, bz, pn, pt,
		pn2, fai, amn, amp, eln, sig, pnt;
	extern double ptr_(float *, int *);
	static float fain, elnc, xmax, xptr, fmax00, pbeta, ptmax, pzmax;
	static int ntrym, ntryx;
	static float trans0, ptmax2, pzmax2;
	extern double ranart_(int *);
	static float xratio;

	ntrym = 0;
	*icou1 = 0;
	pi = 3.1415926f;
	amn = .93892500000000001f;
	amp = .78200000000000003f;
	*dm1 = amn;
	*dm2 = amn;
	v = .43f;
	w = -.84f;
	r__1 = *srt;
	r__2 = *dm1 + *dm2 + amp;
	r__3 = *srt;
	r__4 = *dm1 - amp - *dm2;
	r__5 = *srt;
	ptmax2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / 4.f /
			 (r__5 * r__5);
	ptmax = sqrt(ptmax2) * 1.f / 3.f;
L7:
	pt = ptr_(&ptmax, iseed);
	r__1 = *srt;
	r__2 = *dm1 + *dm2 + amp;
	r__3 = *srt;
	r__4 = *dm1 - amp - *dm2;
	r__5 = *srt;
	r__6 = pt;
	pzmax2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / 4.f /
				 (r__5 * r__5) -
			 r__6 * r__6;
	++ntrym;
	if (pzmax2 < 0.f && ntrym <= 100)
	{
		goto L7;
	}
	else
	{
		pzmax2 = 1e-9f;
	}
	pzmax = sqrt(pzmax2);
	xmax = pzmax * 2.f / *srt;
	ntryx = 0;
	fmax00 = 1.056f;
	x00 = .26f;
	if (dabs(xmax) > .26f)
	{
		f00 = fmax00;
	}
	else
	{
		r__1 = xmax;
		f00 = v * dabs(xmax) + 1.f + w * (r__1 * r__1);
	}
L9:
	x = xmax * (1.f - ranart_(&rndf77_1.nseed) * 2.f);
	++ntryx;
	r__1 = x;
	xratio = (v * dabs(x) + 1.f + w * (r__1 * r__1)) / f00;
	if (xratio < ranart_(&rndf77_1.nseed) && ntryx <= 50)
	{
		goto L9;
	}
	*pz = *srt * .5f * x;
	fai = pi * 2.f * ranart_(&rndf77_1.nseed);
	*px = pt * cos(fai);
	*py = pt * sin(fai);
	r__1 = *dm1;
	r__2 = pt;
	r__3 = *pz;
	ek = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
	eln = *srt - ek;
	if (eln <= 0.f)
	{
		*icou1 = -1;
		return 0;
	}
	bx = -(*px) / eln;
	by = -(*py) / eln;
	bz = -(*pz) / eln;
	r__1 = bx;
	r__2 = by;
	r__3 = bz;
	ga = 1.f / sqrt(1.f - r__1 * r__1 - r__2 * r__2 - r__3 * r__3);
	elnc = eln / ga;
	r__2 = elnc;
	r__3 = *dm2;
	r__4 = amp;
	r__1 = (r__2 * r__2 + r__3 * r__3 - r__4 * r__4) / (elnc * 2.f);
	r__5 = *dm2;
	pn2 = r__1 * r__1 - r__5 * r__5;
	if (pn2 <= 0.f)
	{
		*icou1 = -1;
		return 0;
	}
	pn = sqrt(pn2);
	xptr = pn * .33f;
	pnt = ptr_(&xptr, iseed);
	fain = pi * 2.f * ranart_(&rndf77_1.nseed);
	*pnx = pnt * cos(fain);
	*pny = pnt * sin(fain);
	sig = 1.f;
	if (x > 0.f)
	{
		sig = -1.f;
	}
	r__1 = pn;
	r__2 = pnt;
	*pnz = sig * sqrt(r__1 * r__1 - r__2 * r__2);
	r__1 = *dm2;
	r__2 = *pnx;
	r__3 = *pny;
	r__4 = *pnz;
	en = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	*ppx = -(*pnx);
	*ppy = -(*pny);
	*ppz = -(*pnz);
	r__1 = amp;
	r__2 = *ppx;
	r__3 = *ppy;
	r__4 = *ppz;
	ep = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	pbeta = *pnx * bx + *pny * by + *pnz * bz;
	trans0 = ga * (ga * pbeta / (ga + 1.f) + en);
	*pnx = bx * trans0 + *pnx;
	*pny = by * trans0 + *pny;
	*pnz = bz * trans0 + *pnz;
	if (ep == 0.f)
	{
		ep = 1e-9f;
	}
	pbeta = *ppx * bx + *ppy * by + *ppz * bz;
	trans0 = ga * (ga * pbeta / (ga + 1.f) + ep);
	*ppx = bx * trans0 + *ppx;
	*ppy = by * trans0 + *ppy;
	*ppz = bz * trans0 + *ppz;
	return 0;
}

double rmass_(float *dmax__, int *iseed)
{
	float ret_val;

	static float dm, fm, dmin__;
	static int ntry1;
	extern double fdelta_(float *), ranart_(int *);

	dmin__ = 1.078f;
	if (*dmax__ < 1.232f)
	{
		fm = fdelta_(dmax__);
	}
	else
	{
		fm = 1.f;
	}
	if (fm == 0.f)
	{
		fm = 1e-6f;
	}
	ntry1 = 0;
L10:
	dm = ranart_(&rndf77_1.nseed) * (*dmax__ - dmin__) + dmin__;
	++ntry1;
	if (ranart_(&rndf77_1.nseed) > fdelta_(&dm) / fm && ntry1 <= 10)
	{
		goto L10;
	}
	if (dm > 1.47f)
	{
		goto L10;
	}
	ret_val = dm;
	return ret_val;
}

double frho_(float *dmass)
{
	float ret_val, r__1, r__2, r__3;

	static float fd, am0, wid;

	am0 = .77f;
	wid = .153f;
	r__1 = wid;
	r__2 = *dmass - am0;
	r__3 = wid;
	fd = r__1 * r__1 * .25f / (r__2 * r__2 + r__3 * r__3 * .25f);
	ret_val = fd;
	return ret_val;
}

double rhomas_(float *dmax__, int *iseed)
{
	float ret_val;

	static float dm, fm, dmin__;
	extern double frho_(float *);
	static int ntry1;
	extern double ranart_(int *);

	dmin__ = .28f;
	if (*dmax__ < .77f)
	{
		fm = frho_(dmax__);
	}
	else
	{
		fm = 1.f;
	}
	if (fm == 0.f)
	{
		fm = 1e-6f;
	}
	ntry1 = 0;
L10:
	dm = ranart_(&rndf77_1.nseed) * (*dmax__ - dmin__) + dmin__;
	++ntry1;
	if (ranart_(&rndf77_1.nseed) > frho_(&dm) / fm && ntry1 <= 10)
	{
		goto L10;
	}
	if (dm > 1.07f)
	{
		goto L10;
	}
	ret_val = dm;
	return ret_val;
}

double x2pi_(float *srt)
{
	static float earray[15] = {2.23f, 2.81f, 3.67f, 4.f, 4.95f, 5.52f, 5.97f, 6.04f,
							  6.6f, 6.9f, 7.87f, 8.11f, 10.01f, 16.f, 19.f};
	static float xarray[15] = {1.22f, 2.51f, 2.67f, 2.95f, 2.96f, 2.84f, 2.8f, 3.2f,
							  2.7f, 3.f, 2.54f, 2.46f, 2.4f, 1.66f, 1.5f};

	float ret_val, r__1, r__2, r__3, r__4;

	double sqrt(double), log(double), exp(double);

	static int ie;
	static float plab, xmin, ymin, xmax, ymax, pmass;

	pmass = .9383f;
	ret_val = 1e-6f;
	if (*srt <= 2.2f)
	{
		return ret_val;
	}
	r__2 = *srt;
	r__3 = pmass;
	r__1 = (r__2 * r__2 - r__3 * r__3 * 2.f) / (pmass * 2.f);
	r__4 = pmass;
	plab = sqrt(r__1 * r__1 - r__4 * r__4);
	if (plab < earray[0])
	{
		ret_val = xarray[0];
		return ret_val;
	}

	for (ie = 1; ie <= 15; ++ie)
	{
		if (earray[ie - 1] == plab)
		{
			ret_val = xarray[ie - 1];
			return ret_val;
		}
		else if (earray[ie - 1] > plab)
		{
			ymin = log(xarray[ie - 2]);
			ymax = log(xarray[ie - 1]);
			xmin = log(earray[ie - 2]);
			xmax = log(earray[ie - 1]);
			ret_val = exp(ymin + (log(plab) - xmin) * (ymax - ymin) / (xmax - xmin));
			return ret_val;
		}
	}
	return ret_val;
}

double x3pi_(float *srt)
{
	static float xarray[12] = {.02f, .4f, 1.15f, 1.6f, 2.19f, 2.85f, 2.3f, 3.1f,
							  2.47f, 2.6f, 2.4f, 1.7f};
	static float earray[12] = {2.23f, 2.81f, 3.67f, 4.f, 4.95f, 5.52f, 5.97f, 6.04f,
							  6.6f, 6.9f, 10.01f, 19.f};

	float ret_val, r__1, r__2, r__3, r__4;

	double sqrt(double), log(double), exp(double);

	static int ie;
	static float plab, xmin, ymin, xmax, ymax, pmass;

	pmass = .9383f;
	ret_val = 1e-6f;
	if (*srt <= 2.3f)
	{
		return ret_val;
	}
	r__2 = *srt;
	r__3 = pmass;
	r__1 = (r__2 * r__2 - r__3 * r__3 * 2.f) / (pmass * 2.f);
	r__4 = pmass;
	plab = sqrt(r__1 * r__1 - r__4 * r__4);
	if (plab < earray[0])
	{
		ret_val = xarray[0];
		return ret_val;
	}

	for (ie = 1; ie <= 12; ++ie)
	{
		if (earray[ie - 1] == plab)
		{
			ret_val = xarray[ie - 1];
			return ret_val;
		}
		else if (earray[ie - 1] > plab)
		{
			ymin = log(xarray[ie - 2]);
			ymax = log(xarray[ie - 1]);
			xmin = log(earray[ie - 2]);
			xmax = log(earray[ie - 1]);
			ret_val = exp(ymin + (log(plab) - xmin) * (ymax - ymin) / (xmax - xmin));
			return ret_val;
		}
	}
	return ret_val;
}

double x33pi_(float *srt)
{
	static float xarray[12] = {.02f, .22f, .74f, 1.1f, 1.76f, 1.84f, 2.2f, 2.4f,
							  2.15f, 2.6f, 2.3f, 1.7f};
	static float earray[12] = {2.23f, 2.81f, 3.67f, 4.f, 4.95f, 5.52f, 5.97f, 6.04f,
							  6.6f, 6.9f, 10.01f, 19.f};

	float ret_val, r__1, r__2, r__3, r__4;

	double sqrt(double), log(double), exp(double);

	static int ie;
	static float plab, xmin, ymin, xmax, ymax, pmass;

	pmass = .9383f;
	ret_val = 1e-6f;
	if (*srt <= 2.3f)
	{
		return ret_val;
	}
	r__2 = *srt;
	r__3 = pmass;
	r__1 = (r__2 * r__2 - r__3 * r__3 * 2.f) / (pmass * 2.f);
	r__4 = pmass;
	plab = sqrt(r__1 * r__1 - r__4 * r__4);
	if (plab < earray[0])
	{
		ret_val = xarray[0];
		return ret_val;
	}

	for (ie = 1; ie <= 12; ++ie)
	{
		if (earray[ie - 1] == plab)
		{
			ret_val = xarray[ie - 1];
			return ret_val;
		}
		else if (earray[ie - 1] > plab)
		{
			ymin = log(xarray[ie - 2]);
			ymax = log(xarray[ie - 1]);
			xmin = log(earray[ie - 2]);
			xmax = log(earray[ie - 1]);
			ret_val = exp(ymin + (log(plab) - xmin) * (ymax - ymin) / (xmax - xmin));
			return ret_val;
		}
	}
	return ret_val;
}

double x4pi_(float *srt)
{
	float ret_val, r__1, r__2, r__3, r__4, r__5;

	double sqrt(double);

	static float al, as, es, ak0;
	extern double pp1_(float *);
	static float xk1, xk2, xk3, xk4, ada, ana, akp;
	extern double s1535_(float *), ppk0_(float *), ppk1_(float *), x2pi_(float *), x3pi_(float *);
	static float pps1, pps2;
	extern double x33pi_(float *);
	static float t1dlk, t2dlk, t1dsk, t2dsk, t1nlk, t2nlk, t1nsk, t2nsk;
	extern double sigma_(float *, int *, int *, int *);
	static float pmdlk, pmdsk, xkaon, pmnlk, pmass;
	extern double pplpk_(float *);
	static float pmnsk, pmdlk2, pmdsk2, pmnlk2, pmnsk2, xpp2pi, xpp3pi, ppsngl;

	akp = .498f;
	ak0 = .498f;
	ana = .94f;
	ada = 1.232f;
	al = 1.1157f;
	as = 1.1197f;
	pmass = .9383f;
	es = *srt;
	if (es <= 4.f)
	{
		ret_val = 0.f;
	}
	else
	{
		xpp2pi = x2pi_(&es) * 4.f;
		xpp3pi = (x3pi_(&es) + x33pi_(&es)) * 3.f;
		pps1 = sigma_(&es, &c__1, &c__1, &c__0) + sigma_(&es, &c__1, &c__1, &c__1) * .5f;
		pps2 = sigma_(&es, &c__1, &c__1, &c__1) * 1.5f;
		ppsngl = pps1 + pps2 + s1535_(&es);
		xk1 = 0.f;
		xk2 = 0.f;
		xk3 = 0.f;
		xk4 = 0.f;
		t1nlk = ana + al + akp;
		t2nlk = ana + al - akp;
		if (es <= t1nlk)
		{
			goto L333;
		}
		r__1 = es;
		r__2 = t1nlk;
		r__3 = es;
		r__4 = t2nlk;
		r__5 = es;
		pmnlk2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / (r__5 * r__5 * 4.f);
		pmnlk = sqrt(pmnlk2);
		xk1 = pplpk_(&es);
		t1dlk = ada + al + akp;
		t2dlk = ada + al - akp;
		if (es <= t1dlk)
		{
			goto L333;
		}
		r__1 = es;
		r__2 = t1dlk;
		r__3 = es;
		r__4 = t2dlk;
		r__5 = es;
		pmdlk2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / (r__5 * r__5 * 4.f);
		pmdlk = sqrt(pmdlk2);
		xk3 = pplpk_(&es);
		t1nsk = ana + as + akp;
		t2nsk = ana + as - akp;
		if (es <= t1nsk)
		{
			goto L333;
		}
		r__1 = es;
		r__2 = t1nsk;
		r__3 = es;
		r__4 = t2nsk;
		r__5 = es;
		pmnsk2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / (r__5 * r__5 * 4.f);
		pmnsk = sqrt(pmnsk2);
		xk2 = ppk1_(&es) + ppk0_(&es);
		t1dsk = ada + as + akp;
		t2dsk = ada + as - akp;
		if (es <= t1dsk)
		{
			goto L333;
		}
		r__1 = es;
		r__2 = t1dsk;
		r__3 = es;
		r__4 = t2dsk;
		r__5 = es;
		pmdsk2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / (r__5 * r__5 * 4.f);
		pmdsk = sqrt(pmdsk2);
		xk4 = ppk1_(&es) + ppk0_(&es);
	L333:
		xkaon = (xk1 + xk2 + xk3 + xk4) * 3.f;
		ret_val = pp1_(&es) - ppsngl - xpp2pi - xpp3pi - xkaon;
		if (ret_val <= 0.f)
		{
			ret_val = 1e-6f;
		}
	}
	return ret_val;
}

double pp1_(float *srt)
{
	float ret_val, r__1, r__2, r__3, r__4;
	double d__1, d__2;

	double sqrt(double), log(double), pow_dd(double *, double *);

	static float a, b, c__, d__, an, plab, pmin, pmax, plab2, pmass;

	pmass = .9383f;
	ret_val = 0.f;
	r__2 = *srt;
	r__3 = pmass;
	r__1 = (r__2 * r__2 - r__3 * r__3 * 2.f) / (pmass * 2.f);
	r__4 = pmass;
	plab2 = r__1 * r__1 - r__4 * r__4;
	if (plab2 <= 0.f)
	{
		return ret_val;
	}
	plab = sqrt(plab2);
	pmin = .968f;
	pmax = 2080.f;
	if (plab < pmin || plab > pmax)
	{
		ret_val = 0.f;
		return ret_val;
	}
	a = 30.9f;
	b = -28.9f;
	c__ = .192f;
	d__ = -.835f;
	an = -2.46f;
	d__1 = (double)plab;
	d__2 = (double)an;
	r__1 = log(plab);
	ret_val = a + b * pow_dd(&d__1, &d__2) + c__ * (r__1 * r__1);
	if (ret_val <= 0.f)
	{
		ret_val = 0.f;
	}
	return ret_val;
}

double pp2_(float *srt)
{
	float ret_val, r__1, r__2, r__3, r__4;
	double d__1, d__2;

	double sqrt(double), log(double), pow_dd(double *, double *);

	static float a, b, c__, d__, an, plab, pmin, pmax, pmass;

	pmass = .9383f;
	r__2 = *srt;
	r__3 = pmass;
	r__1 = (r__2 * r__2 - r__3 * r__3 * 2.f) / (pmass * 2.f);
	r__4 = pmass;
	plab = sqrt(r__1 * r__1 - r__4 * r__4);
	pmin = 2.f;
	pmax = 2050.f;
	if (plab > pmax)
	{
		ret_val = 8.f;
		return ret_val;
	}
	if (plab < pmin)
	{
		ret_val = 25.f;
		return ret_val;
	}
	a = 11.2f;
	b = 25.5f;
	c__ = .151f;
	d__ = -1.62f;
	an = -1.12f;
	d__1 = (double)plab;
	d__2 = (double)an;
	r__1 = log(plab);
	ret_val = a + b * pow_dd(&d__1, &d__2) + c__ * (r__1 * r__1) + d__ * log(plab);
	if (ret_val <= 0.f)
	{
		ret_val = 0.f;
	}
	return ret_val;
}

double ppt_(float *srt)
{
	float ret_val, r__1, r__2, r__3, r__4;
	double d__1, d__2;

	double sqrt(double), log(double), pow_dd(double *, double *);

	static float a, b, c__, d__, an, plab, pmin, pmax, pmass;

	pmass = .9383f;
	r__2 = *srt;
	r__3 = pmass;
	r__1 = (r__2 * r__2 - r__3 * r__3 * 2.f) / (pmass * 2.f);
	r__4 = pmass;
	plab = sqrt(r__1 * r__1 - r__4 * r__4);
	pmin = 3.f;
	pmax = 2100.f;
	if (plab < pmin || plab > pmax)
	{
		ret_val = 55.f;
		return ret_val;
	}
	a = 45.6f;
	b = 219.f;
	c__ = .41f;
	d__ = -3.41f;
	an = -4.23f;
	d__1 = (double)plab;
	d__2 = (double)an;
	r__1 = log(plab);
	ret_val = a + b * pow_dd(&d__1, &d__2) + c__ * (r__1 * r__1) + d__ * log(plab);
	if (ret_val <= 0.f)
	{
		ret_val = 0.f;
	}
	return ret_val;
}

double s1535_(float *srt)
{
	float ret_val, r__1;

	static float s0;

	s0 = 2.424f;
	ret_val = 0.f;
	if (*srt <= s0)
	{
		return ret_val;
	}
	r__1 = *srt - s0;
	ret_val = (*srt - s0) * .20399999999999999f / (r__1 * r__1 + .058f);
	return ret_val;
}

int tablem_(void)
{
	static int l;
	static float x, rr, anorm;
	extern double ptdis_(float *);
	static float ptmax;

	ptmax = 2.01f;
	anorm = ptdis_(&ptmax);
	for (l = 0; l <= 200; ++l)
	{
		x = (float)(l + 1) * .01f;
		rr = ptdis_(&x) / anorm;
		table_1.earray[l] = rr;
		table_1.xarray[l] = x;
	}
	return 0;
}

double ptdis_(float *x)
{
	float ret_val, r__1, r__2;

	double exp(double);

	static float b, c__, d__;

	b = 3.78f;
	c__ = .47f;
	d__ = 3.6f;
	r__1 = *x;
	r__2 = d__;
	ret_val = 1.f / (b * 2.f) * (1.f - exp(-b * (r__1 * r__1))) - c__ / d__ * *x * exp(-d__ * *x) - c__ / (r__2 * r__2) * (exp(-d__ * *x) - 1.f);
	return ret_val;
}

int ppxs_(int *lb1, int *lb2, float *srt, float *ppsig,
		  float *spprho, int *ipp)
{
	float r__1, r__2;

	double sqrt(double), atan(double), sin(double);

	static float q, s0, s1, s2, d00, d11, d20, erh, esi, erho, trho, esigma,
		tsigma;

	*ppsig = 0.f;
	*spprho = 0.f;
	*ipp = 0;
	if (*srt <= .3f)
	{
		return 0;
	}
	r__1 = *srt / 2;
	q = sqrt(r__1 * r__1 - .019600000000000003f);
	esigma = .81200000000000006f;
	tsigma = q * 2.06f;
	erho = .77f;
	r__2 = q / erho;
	r__1 = q / .14f / (r__2 * r__2 + 1.f);
	trho = q * .095f * (r__1 * r__1);
	esi = esigma - *srt;
	if (esi == 0.f)
	{
		d00 = 1.5707963f;
		goto L10;
	}
	d00 = atan(tsigma / 2.f / esi);
L10:
	erh = erho - *srt;
	if (erh == 0.f)
	{
		d11 = 1.5707963f;
		goto L20;
	}
	d11 = atan(trho / 2.f / erh);
L20:
	d20 = q * -.12f / .14f;
	r__1 = sin(d00);
	r__2 = q;
	s0 = r__1 * r__1 * 25.132740800000001f / (r__2 * r__2);
	r__1 = sin(d11);
	r__2 = q;
	s1 = r__1 * r__1 * 75.398222400000009f / (r__2 * r__2);
	r__1 = sin(d20);
	r__2 = q;
	s2 = r__1 * r__1 * 125.663704f / (r__2 * r__2);
	s0 = s0 * .038809000000000003f * 10.f;
	s1 = s1 * .038809000000000003f * 10.f;
	s2 = s2 * .038809000000000003f * 10.f;
	*spprho = s1 / 2.f;
	if (*lb1 == 5 && *lb2 == 5)
	{
		*ipp = 1;
		*ppsig = s2;
		return 0;
	}
	if (*lb1 == 5 && *lb2 == 4 || *lb1 == 4 && *lb2 == 5)
	{
		*ipp = 2;
		*ppsig = s2 / 2.f + s1 / 2.f;
		return 0;
	}
	if (*lb1 == 5 && *lb2 == 3 || *lb1 == 3 && *lb2 == 5)
	{
		*ipp = 3;
		*ppsig = s2 / 6.f + s1 / 2.f + s0 / 3.f;
		return 0;
	}
	if (*lb1 == 4 && *lb2 == 4)
	{
		*ipp = 4;
		*ppsig = s2 * 2 / 3.f + s0 / 3.f;
		return 0;
	}
	if (*lb1 == 4 && *lb2 == 3 || *lb1 == 3 && *lb2 == 4)
	{
		*ipp = 5;
		*ppsig = s2 / 2.f + s1 / 2.f;
		return 0;
	}
	if (*lb1 == 3 && *lb2 == 3)
	{
		*ipp = 6;
		*ppsig = s2;
	}
	return 0;
}

double pplpk_(float *srt)
{
	float ret_val, r__1, r__2, r__3, r__4;
	double d__1, d__2;

	double sqrt(double), log(double), pow_dd(double *, double *);

	static float a, b, c__, an, plab, pmin, pmax, pmass;

	pmass = .9383f;
	ret_val = 0.f;
	r__2 = *srt;
	r__3 = pmass;
	r__1 = (r__2 * r__2 - r__3 * r__3 * 2.f) / (pmass * 2.f);
	r__4 = pmass;
	plab = sqrt(r__1 * r__1 - r__4 * r__4);
	pmin = 2.82f;
	pmax = 25.f;
	if (plab > pmax)
	{
		ret_val = .036f;
		return ret_val;
	}
	if (plab < pmin)
	{
		ret_val = 0.f;
		return ret_val;
	}
	a = .0654f;
	b = -3.16f;
	c__ = -.0029f;
	an = -4.14f;
	d__1 = (double)plab;
	d__2 = (double)an;
	r__1 = log(plab);
	ret_val = a + b * pow_dd(&d__1, &d__2) + c__ * (r__1 * r__1);
	if (ret_val <= 0.f)
	{
		ret_val = 0.f;
	}
	return ret_val;
}

double ppk0_(float *srt)
{
	static float xarray[7] = {.03f, .025f, .025f, .026f, .02f, .014f, .06f};
	static float earray[7] = {3.67f, 4.95f, 5.52f, 6.05f, 6.92f, 7.87f, 10.f};

	float ret_val, r__1, r__2, r__3, r__4;

	double sqrt(double), log(double), exp(double);

	static int ie;
	static float plab, xmin, ymin, xmax, ymax, pmass;

	pmass = .9383f;
	ret_val = 0.f;
	if (*srt <= 2.63f)
	{
		return ret_val;
	}
	if (*srt > 4.54f)
	{
		ret_val = .037f;
		return ret_val;
	}
	r__2 = *srt;
	r__3 = pmass;
	r__1 = (r__2 * r__2 - r__3 * r__3 * 2.f) / (pmass * 2.f);
	r__4 = pmass;
	plab = sqrt(r__1 * r__1 - r__4 * r__4);
	if (plab < earray[0])
	{
		ret_val = xarray[0];
		return ret_val;
	}

	for (ie = 1; ie <= 7; ++ie)
	{
		if (earray[ie - 1] == plab)
		{
			ret_val = xarray[ie - 1];
			goto L10;
		}
		else if (earray[ie - 1] > plab)
		{
			ymin = log(xarray[ie - 2]);
			ymax = log(xarray[ie - 1]);
			xmin = log(earray[ie - 2]);
			xmax = log(earray[ie - 1]);
			ret_val = exp(ymin + (log(plab) - xmin) * (ymax - ymin) / (xmax - xmin));
			goto L10;
		}
	}
L10:
	return ret_val;
}

double ppk1_(float *srt)
{
	static float xarray[7] = {.013f, .025f, .016f, .012f, .017f, .029f, .025f};
	static float earray[7] = {3.67f, 4.95f, 5.52f, 5.97f, 6.05f, 6.92f, 7.87f};

	float ret_val, r__1, r__2, r__3, r__4;

	double sqrt(double), log(double), exp(double);

	static int ie;
	static float plab, xmin, ymin, xmax, ymax, pmass;

	pmass = .9383f;
	ret_val = 0.f;
	if (*srt <= 2.63f)
	{
		return ret_val;
	}
	if (*srt > 4.08f)
	{
		ret_val = .025f;
		return ret_val;
	}
	r__2 = *srt;
	r__3 = pmass;
	r__1 = (r__2 * r__2 - r__3 * r__3 * 2.f) / (pmass * 2.f);
	r__4 = pmass;
	plab = sqrt(r__1 * r__1 - r__4 * r__4);
	if (plab < earray[0])
	{
		ret_val = xarray[0];
		return ret_val;
	}

	for (ie = 1; ie <= 7; ++ie)
	{
		if (earray[ie - 1] == plab)
		{
			ret_val = xarray[ie - 1];
			goto L10;
		}
		else if (earray[ie - 1] > plab)
		{
			ymin = log(xarray[ie - 2]);
			ymax = log(xarray[ie - 1]);
			xmin = log(earray[ie - 2]);
			xmax = log(earray[ie - 1]);
			ret_val = exp(ymin + (log(plab) - xmin) * (ymax - ymin) / (xmax - xmin));
			goto L10;
		}
	}
L10:
	return ret_val;
}

int crpn_(float *px, float *py, float *pz, float *srt, int *i1, int *i2, int *iblock, float *xkaon0, float *xkaon, float *xphi, float *xphin)
{
	int i__1, i__2;
	float r__1, r__2, r__3, r__4, r__5;

	double sqrt(double), cos(double), sin(double);

	static float c1, s1, t1, x1, x2, dm;
	static int ii, jj;
	static float pr, cc1, em1, em2, ct1, pr2, px0, py0, pz0, st1;
	static int ipi;
	extern double ptr_(float *, int *);
	static int ntag;
	static float dmax__, arho, xptr;
	static int kaonc, ianti;
	extern double pnlka_(float *), pnska_(float *), rmass_(float *, int *), twopi_(float *);
	static float aomega;
	extern int rmasdd_(float *, float *, float *, float *, float *, int *, int *, float *, float *);
	static float ameson;
	extern double ranart_(int *), threpi_(float *);
	extern int rotate_(float *, float *, float *, float *, float *, float *);
	extern double fourpi_(float *);

	px0 = *px;
	py0 = *py;
	pz0 = *pz;
	*iblock = 1;
	x1 = ranart_(&rndf77_1.nseed);
	ianti = 0;
	if (ee_1.lb[*i1 - 1] < 0 || ee_1.lb[*i2 - 1] < 0)
	{
		ianti = 1;
	}
	if (*xkaon0 / (*xkaon + *xphi) >= x1)
	{
		*iblock = 7;
		if (ianti == 1)
		{
			*iblock = -7;
		}
		ntag = 0;
		kaonc = 0;
		if (pnlka_(srt) / (pnlka_(srt) + pnska_(srt)) > ranart_(&rndf77_1.nseed))
		{
			kaonc = 1;
		}
		if (cc_1.e[*i1 - 1] <= .2f)
		{
			ee_1.lb[*i1 - 1] = 23;
			cc_1.e[*i1 - 1] = .498f;
			if (kaonc == 1)
			{
				ee_1.lb[*i2 - 1] = 14;
				cc_1.e[*i2 - 1] = 1.1157f;
			}
			else
			{
				ee_1.lb[*i2 - 1] = (int)(ranart_(&rndf77_1.nseed) * 3) +
								   15;
				cc_1.e[*i2 - 1] = 1.1974f;
			}
			if (ianti == 1)
			{
				ee_1.lb[*i1 - 1] = 21;
				ee_1.lb[*i2 - 1] = -ee_1.lb[*i2 - 1];
			}
		}
		else
		{
			ee_1.lb[*i2 - 1] = 23;
			cc_1.e[*i2 - 1] = .498f;
			if (kaonc == 1)
			{
				ee_1.lb[*i1 - 1] = 14;
				cc_1.e[*i1 - 1] = 1.1157f;
			}
			else
			{
				ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 3) +
								   15;
				cc_1.e[*i1 - 1] = 1.1974f;
			}
			if (ianti == 1)
			{
				ee_1.lb[*i2 - 1] = 21;
				ee_1.lb[*i1 - 1] = -ee_1.lb[*i1 - 1];
			}
		}
		em1 = cc_1.e[*i1 - 1];
		em2 = cc_1.e[*i2 - 1];
		goto L50;
	}
	else if (*xphi / (*xkaon + *xphi) >= x1)
	{
		*iblock = 222;
		if (*xphin / *xphi >= ranart_(&rndf77_1.nseed))
		{
			ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 2) + 1;
			cc_1.e[*i1 - 1] = .939457f;
		}
		else
		{
			ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 4) + 6;
			cc_1.e[*i1 - 1] = 1.232f;
		}
		if (ianti == 1)
		{
			ee_1.lb[*i1 - 1] = -ee_1.lb[*i1 - 1];
		}
		ee_1.lb[*i2 - 1] = 29;
		cc_1.e[*i2 - 1] = 1.02f;
		em1 = cc_1.e[*i1 - 1];
		em2 = cc_1.e[*i2 - 1];
		goto L50;
	}
	else
	{
		if (ranart_(&rndf77_1.nseed) <= twopi_(srt) / (twopi_(srt) + threpi_(srt) + fourpi_(srt)))
		{
			*iblock = 77;
		}
		else
		{
			if (threpi_(srt) / (threpi_(srt) + fourpi_(srt)) > ranart_(&rndf77_1.nseed))
			{
				*iblock = 78;
			}
			else
			{
				*iblock = 79;
			}
		}
		ntag = 0;
		x2 = ranart_(&rndf77_1.nseed);
		if (*iblock == 77)
		{
			dmax__ = *srt - .13496f - .02f;
			dm = rmass_(&dmax__, &input1_1.iseed);
			if (ee_1.lb[*i1 - 1] == 1 && ee_1.lb[*i2 - 1] == 5 || ee_1.lb[*i1 - 1] == 5 && ee_1.lb[*i2 - 1] == 1 || (ee_1.lb[*i1 - 1] == -1 && ee_1.lb[*i2 - 1] == 3 || ee_1.lb[*i1 - 1] == 3 && ee_1.lb[*i2 - 1] == -1))
			{
				if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 1)
				{
					ii = *i1;
					if (x2 <= .5f)
					{
						ee_1.lb[*i1 - 1] = 8;
						cc_1.e[*i1 - 1] = dm;
						ee_1.lb[*i2 - 1] = 5;
						cc_1.e[*i2 - 1] = .13496f;
						goto L40;
					}
					else
					{
						ee_1.lb[*i1 - 1] = 9;
						cc_1.e[*i1 - 1] = dm;
						ee_1.lb[*i2 - 1] = 4;
						ipi = 4;
						cc_1.e[*i2 - 1] = .13496f;
						goto L40;
					}
				}
				else
				{
					ii = *i2;
					if (x2 <= .5f)
					{
						ee_1.lb[*i2 - 1] = 8;
						cc_1.e[*i2 - 1] = dm;
						ee_1.lb[*i1 - 1] = 5;
						cc_1.e[*i1 - 1] = .13496f;
						goto L40;
					}
					else
					{
						ee_1.lb[*i2 - 1] = 9;
						cc_1.e[*i2 - 1] = dm;
						ee_1.lb[*i1 - 1] = 4;
						cc_1.e[*i1 - 1] = .13496f;
						goto L40;
					}
				}
			}
			if (ee_1.lb[*i1 - 1] == 1 && ee_1.lb[*i2 - 1] == 3 || ee_1.lb[*i1 - 1] == 3 && ee_1.lb[*i2 - 1] == 1 || (ee_1.lb[*i1 - 1] == -1 && ee_1.lb[*i2 - 1] == 5 || ee_1.lb[*i1 - 1] == 5 && ee_1.lb[*i2 - 1] == -1))
			{
				if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 1)
				{
					ii = *i1;
					if (x2 <= .33f)
					{
						ee_1.lb[*i1 - 1] = 6;
						cc_1.e[*i1 - 1] = dm;
						ee_1.lb[*i2 - 1] = 5;
						cc_1.e[*i2 - 1] = .13496f;
						goto L40;
					}
					if (x2 > .33f && x2 <= .67f)
					{
						ee_1.lb[*i1 - 1] = 7;
						cc_1.e[*i1 - 1] = dm;
						ee_1.lb[*i2 - 1] = 4;
						cc_1.e[*i2 - 1] = .13496f;
						goto L40;
					}
					if (x2 > .67f)
					{
						ee_1.lb[*i1 - 1] = 8;
						cc_1.e[*i1 - 1] = dm;
						ee_1.lb[*i2 - 1] = 3;
						cc_1.e[*i2 - 1] = .13496f;
						goto L40;
					}
				}
				else
				{
					ii = *i2;
					if (x2 <= .33f)
					{
						ee_1.lb[*i2 - 1] = 6;
						cc_1.e[*i2 - 1] = dm;
						ee_1.lb[*i1 - 1] = 5;
						cc_1.e[*i1 - 1] = .13496f;
						goto L40;
					}
					if (x2 > .33f && x2 <= .67f)
					{
						ee_1.lb[*i2 - 1] = 7;
						cc_1.e[*i2 - 1] = dm;
						ee_1.lb[*i1 - 1] = 4;
						cc_1.e[*i1 - 1] = .13496f;
						goto L40;
					}
					if (x2 > .67f)
					{
						ee_1.lb[*i2 - 1] = 8;
						cc_1.e[*i2 - 1] = dm;
						ee_1.lb[*i1 - 1] = 3;
						cc_1.e[*i1 - 1] = .13496f;
						goto L40;
					}
				}
			}
			if (ee_1.lb[*i1 - 1] == 2 && ee_1.lb[*i2 - 1] == 5 || ee_1.lb[*i1 - 1] == 5 && ee_1.lb[*i2 - 1] == 2 || (ee_1.lb[*i1 - 1] == -2 && ee_1.lb[*i2 - 1] == 3 || ee_1.lb[*i1 - 1] == 3 && ee_1.lb[*i2 - 1] == -2))
			{
				if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2)
				{
					ii = *i1;
					if (x2 <= .33f)
					{
						ee_1.lb[*i1 - 1] = 8;
						cc_1.e[*i1 - 1] = dm;
						ee_1.lb[*i2 - 1] = 4;
						cc_1.e[*i2 - 1] = .13496f;
						goto L40;
					}
					if (x2 > .33f && x2 <= .67f)
					{
						ee_1.lb[*i1 - 1] = 7;
						cc_1.e[*i1 - 1] = dm;
						ee_1.lb[*i2 - 1] = 5;
						cc_1.e[*i2 - 1] = .13496f;
						goto L40;
					}
					if (x2 > .67f)
					{
						ee_1.lb[*i1 - 1] = 9;
						cc_1.e[*i1 - 1] = dm;
						ee_1.lb[*i2 - 1] = 3;
						cc_1.e[*i2 - 1] = .13496f;
						goto L40;
					}
				}
				else
				{
					ii = *i2;
					if (x2 <= .33f)
					{
						ee_1.lb[*i2 - 1] = 8;
						cc_1.e[*i2 - 1] = dm;
						ee_1.lb[*i1 - 1] = 4;
						cc_1.e[*i1 - 1] = .13496f;
						goto L40;
					}
					if (x2 > .33f && x2 <= .67f)
					{
						ee_1.lb[*i2 - 1] = 7;
						cc_1.e[*i2 - 1] = dm;
						ee_1.lb[*i1 - 1] = 5;
						cc_1.e[*i1 - 1] = .13496f;
						goto L40;
					}
					if (x2 > .67f)
					{
						ee_1.lb[*i2 - 1] = 9;
						cc_1.e[*i2 - 1] = dm;
						ee_1.lb[*i1 - 1] = 3;
						cc_1.e[*i1 - 1] = .13496f;
						goto L40;
					}
				}
			}
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 1 && ee_1.lb[*i2 - 1] == 4 || ee_1.lb[*i1 - 1] == 4 && (i__2 = ee_1.lb[*i2 - 1],
																												abs(i__2)) == 1)
			{
				if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 1)
				{
					ii = *i1;
					if (x2 <= .33f)
					{
						ee_1.lb[*i1 - 1] = 8;
						cc_1.e[*i1 - 1] = dm;
						ee_1.lb[*i2 - 1] = 4;
						cc_1.e[*i2 - 1] = .13496f;
						goto L40;
					}
					if (x2 > .33f && x2 <= .67f)
					{
						ee_1.lb[*i1 - 1] = 7;
						cc_1.e[*i1 - 1] = dm;
						ee_1.lb[*i2 - 1] = 5;
						cc_1.e[*i2 - 1] = .13496f;
						goto L40;
					}
					if (x2 > .67f)
					{
						ee_1.lb[*i1 - 1] = 9;
						cc_1.e[*i1 - 1] = dm;
						ee_1.lb[*i2 - 1] = 3;
						cc_1.e[*i2 - 1] = .13496f;
						goto L40;
					}
				}
				else
				{
					ii = *i2;
					if (x2 <= .33f)
					{
						ee_1.lb[*i2 - 1] = 8;
						cc_1.e[*i2 - 1] = dm;
						ee_1.lb[*i1 - 1] = 4;
						cc_1.e[*i1 - 1] = .13496f;
						goto L40;
					}
					if (x2 > .33f && x2 <= .67f)
					{
						ee_1.lb[*i2 - 1] = 7;
						cc_1.e[*i2 - 1] = dm;
						ee_1.lb[*i1 - 1] = 5;
						cc_1.e[*i1 - 1] = .13496f;
						goto L40;
					}
					if (x2 > .67f)
					{
						ee_1.lb[*i2 - 1] = 9;
						cc_1.e[*i2 - 1] = dm;
						ee_1.lb[*i1 - 1] = 3;
						cc_1.e[*i1 - 1] = .13496f;
						goto L40;
					}
				}
			}
			if (ee_1.lb[*i1 - 1] == 2 && ee_1.lb[*i2 - 1] == 3 || ee_1.lb[*i1 - 1] == 3 && ee_1.lb[*i2 - 1] == 2 || (ee_1.lb[*i1 - 1] == -2 && ee_1.lb[*i2 - 1] == 5 || ee_1.lb[*i1 - 1] == 5 && ee_1.lb[*i2 - 1] == -2))
			{
				if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2)
				{
					ii = *i1;
					if (x2 <= .5f)
					{
						ee_1.lb[*i1 - 1] = 6;
						cc_1.e[*i1 - 1] = dm;
						ee_1.lb[*i2 - 1] = 4;
						cc_1.e[*i2 - 1] = .13496f;
						goto L40;
					}
					else
					{
						ee_1.lb[*i1 - 1] = 7;
						cc_1.e[*i1 - 1] = dm;
						ee_1.lb[*i2 - 1] = 3;
						cc_1.e[*i2 - 1] = .13496f;
						goto L40;
					}
				}
				else
				{
					ii = *i2;
					if (x2 <= .5f)
					{
						ee_1.lb[*i2 - 1] = 6;
						cc_1.e[*i2 - 1] = dm;
						ee_1.lb[*i1 - 1] = 4;
						cc_1.e[*i1 - 1] = .13496f;
						goto L40;
					}
					else
					{
						ee_1.lb[*i2 - 1] = 7;
						cc_1.e[*i2 - 1] = dm;
						ee_1.lb[*i1 - 1] = 3;
						cc_1.e[*i1 - 1] = .13496f;
						goto L40;
					}
				}
			}
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2 && ee_1.lb[*i2 - 1] == 4 || ee_1.lb[*i1 - 1] == 4 && (i__2 = ee_1.lb[*i2 - 1],
																												abs(i__2)) == 2)
			{
				if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2)
				{
					ii = *i1;
					if (x2 <= .33f)
					{
						ee_1.lb[*i1 - 1] = 7;
						cc_1.e[*i1 - 1] = dm;
						ee_1.lb[*i2 - 1] = 4;
						cc_1.e[*i2 - 1] = .13496f;
						goto L40;
					}
					if (x2 <= .67f && x2 > .33f)
					{
						ee_1.lb[*i1 - 1] = 6;
						cc_1.e[*i1 - 1] = dm;
						ee_1.lb[*i2 - 1] = 5;
						cc_1.e[*i2 - 1] = .13496f;
						goto L40;
					}
					if (x2 > .67f)
					{
						ee_1.lb[*i1 - 1] = 8;
						cc_1.e[*i1 - 1] = dm;
						ee_1.lb[*i2 - 1] = 3;
						cc_1.e[*i2 - 1] = .13496f;
						goto L40;
					}
				}
				else
				{
					ii = *i2;
					if (x2 <= .33f)
					{
						ee_1.lb[*i2 - 1] = 7;
						cc_1.e[*i2 - 1] = dm;
						ee_1.lb[*i1 - 1] = 4;
						cc_1.e[*i1 - 1] = .13496f;
						goto L40;
					}
					if (x2 <= .67f && x2 > .33f)
					{
						ee_1.lb[*i2 - 1] = 6;
						cc_1.e[*i2 - 1] = dm;
						ee_1.lb[*i1 - 1] = 5;
						cc_1.e[*i1 - 1] = .13496f;
						goto L40;
					}
					if (x2 > .67f)
					{
						ee_1.lb[*i2 - 1] = 8;
						cc_1.e[*i2 - 1] = dm;
						ee_1.lb[*i1 - 1] = 3;
						cc_1.e[*i1 - 1] = .13496f;
						goto L40;
					}
				}
			}
		}
		if (*iblock == 78)
		{
			rmasdd_(srt, &c_b185, &c_b503, &c_b187, &c_b505, &input1_1.iseed,
					&c__4, &dm, &ameson);
			arho = ameson;
			if (ee_1.lb[*i1 - 1] == 1 && ee_1.lb[*i2 - 1] == 5 || ee_1.lb[*i1 - 1] == 5 && ee_1.lb[*i2 - 1] == 1 || (ee_1.lb[*i1 - 1] == -1 && ee_1.lb[*i2 - 1] == 3 || ee_1.lb[*i1 - 1] == 3 && ee_1.lb[*i2 - 1] == -1))
			{
				if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 1)
				{
					ii = *i1;
					if (x2 <= .5f)
					{
						ee_1.lb[*i1 - 1] = 8;
						cc_1.e[*i1 - 1] = dm;
						ee_1.lb[*i2 - 1] = 27;
						cc_1.e[*i2 - 1] = arho;
						goto L40;
					}
					else
					{
						ee_1.lb[*i1 - 1] = 9;
						cc_1.e[*i1 - 1] = dm;
						ee_1.lb[*i2 - 1] = 26;
						cc_1.e[*i2 - 1] = arho;
						goto L40;
					}
				}
				else
				{
					ii = *i2;
					if (x2 <= .5f)
					{
						ee_1.lb[*i2 - 1] = 8;
						cc_1.e[*i2 - 1] = dm;
						ee_1.lb[*i1 - 1] = 27;
						cc_1.e[*i1 - 1] = arho;
						goto L40;
					}
					else
					{
						ee_1.lb[*i2 - 1] = 9;
						cc_1.e[*i2 - 1] = dm;
						ee_1.lb[*i1 - 1] = 26;
						cc_1.e[*i1 - 1] = arho;
						goto L40;
					}
				}
			}
			if (ee_1.lb[*i1 - 1] == 1 && ee_1.lb[*i2 - 1] == 3 || ee_1.lb[*i1 - 1] == 3 && ee_1.lb[*i2 - 1] == 1 || (ee_1.lb[*i1 - 1] == -1 && ee_1.lb[*i2 - 1] == 5 || ee_1.lb[*i1 - 1] == 5 && ee_1.lb[*i2 - 1] == -1))
			{
				if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 1)
				{
					ii = *i1;
					if (x2 <= .33f)
					{
						ee_1.lb[*i1 - 1] = 6;
						cc_1.e[*i1 - 1] = dm;
						ee_1.lb[*i2 - 1] = 27;
						cc_1.e[*i2 - 1] = arho;
						goto L40;
					}
					if (x2 > .33f && x2 <= .67f)
					{
						ee_1.lb[*i1 - 1] = 7;
						cc_1.e[*i1 - 1] = dm;
						ee_1.lb[*i2 - 1] = 26;
						cc_1.e[*i2 - 1] = arho;
						goto L40;
					}
					if (x2 > .67f)
					{
						ee_1.lb[*i1 - 1] = 8;
						cc_1.e[*i1 - 1] = dm;
						ee_1.lb[*i2 - 1] = 25;
						cc_1.e[*i2 - 1] = arho;
						goto L40;
					}
				}
				else
				{
					ii = *i2;
					if (x2 <= .33f)
					{
						ee_1.lb[*i2 - 1] = 6;
						cc_1.e[*i2 - 1] = dm;
						ee_1.lb[*i1 - 1] = 27;
						cc_1.e[*i1 - 1] = arho;
						goto L40;
					}
					if (x2 > .33f && x2 <= .67f)
					{
						ee_1.lb[*i2 - 1] = 7;
						cc_1.e[*i2 - 1] = dm;
						ee_1.lb[*i1 - 1] = 26;
						cc_1.e[*i1 - 1] = arho;
						goto L40;
					}
					if (x2 > .67f)
					{
						ee_1.lb[*i2 - 1] = 8;
						cc_1.e[*i2 - 1] = dm;
						ee_1.lb[*i1 - 1] = 25;
						cc_1.e[*i1 - 1] = arho;
						goto L40;
					}
				}
			}
			if (ee_1.lb[*i1 - 1] == 2 && ee_1.lb[*i2 - 1] == 5 || ee_1.lb[*i1 - 1] == 5 && ee_1.lb[*i2 - 1] == 2 || (ee_1.lb[*i1 - 1] == -2 && ee_1.lb[*i2 - 1] == 3 || ee_1.lb[*i1 - 1] == 3 && ee_1.lb[*i2 - 1] == -2))
			{
				if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2)
				{
					ii = *i1;
					if (x2 <= .33f)
					{
						ee_1.lb[*i1 - 1] = 8;
						cc_1.e[*i1 - 1] = dm;
						ee_1.lb[*i2 - 1] = 26;
						cc_1.e[*i2 - 1] = arho;
						goto L40;
					}
					if (x2 > .33f && x2 <= .67f)
					{
						ee_1.lb[*i1 - 1] = 7;
						cc_1.e[*i1 - 1] = dm;
						ee_1.lb[*i2 - 1] = 27;
						cc_1.e[*i2 - 1] = arho;
						goto L40;
					}
					if (x2 > .67f)
					{
						ee_1.lb[*i1 - 1] = 9;
						cc_1.e[*i1 - 1] = dm;
						ee_1.lb[*i2 - 1] = 25;
						cc_1.e[*i2 - 1] = arho;
						goto L40;
					}
				}
				else
				{
					ii = *i2;
					if (x2 <= .33f)
					{
						ee_1.lb[*i2 - 1] = 8;
						cc_1.e[*i2 - 1] = dm;
						ee_1.lb[*i1 - 1] = 26;
						cc_1.e[*i1 - 1] = arho;
						goto L40;
					}
					if (x2 > .33f && x2 <= .67f)
					{
						ee_1.lb[*i2 - 1] = 7;
						cc_1.e[*i2 - 1] = dm;
						ee_1.lb[*i1 - 1] = 27;
						cc_1.e[*i1 - 1] = arho;
						goto L40;
					}
					if (x2 > .67f)
					{
						ee_1.lb[*i2 - 1] = 9;
						cc_1.e[*i2 - 1] = dm;
						ee_1.lb[*i1 - 1] = 25;
						cc_1.e[*i1 - 1] = arho;
						goto L40;
					}
				}
			}
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 1 && ee_1.lb[*i2 - 1] == 4 || ee_1.lb[*i1 - 1] == 4 && (i__2 = ee_1.lb[*i2 - 1],
																												abs(i__2)) == 1)
			{
				if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 1)
				{
					ii = *i1;
					if (x2 <= .33f)
					{
						ee_1.lb[*i1 - 1] = 7;
						cc_1.e[*i1 - 1] = dm;
						ee_1.lb[*i2 - 1] = 27;
						cc_1.e[*i2 - 1] = arho;
						goto L40;
					}
					if (x2 > .33f && x2 <= .67f)
					{
						ee_1.lb[*i1 - 1] = 8;
						cc_1.e[*i1 - 1] = dm;
						ee_1.lb[*i2 - 1] = 26;
						cc_1.e[*i2 - 1] = arho;
						goto L40;
					}
					if (x2 > .67f)
					{
						ee_1.lb[*i1 - 1] = 9;
						cc_1.e[*i1 - 1] = dm;
						ee_1.lb[*i2 - 1] = 25;
						cc_1.e[*i2 - 1] = arho;
						goto L40;
					}
				}
				else
				{
					ii = *i2;
					if (x2 <= .33f)
					{
						ee_1.lb[*i2 - 1] = 7;
						cc_1.e[*i2 - 1] = dm;
						ee_1.lb[*i1 - 1] = 27;
						cc_1.e[*i1 - 1] = arho;
						goto L40;
					}
					if (x2 > .33f && x2 <= .67f)
					{
						ee_1.lb[*i2 - 1] = 8;
						cc_1.e[*i2 - 1] = dm;
						ee_1.lb[*i1 - 1] = 26;
						cc_1.e[*i1 - 1] = arho;
						goto L40;
					}
					if (x2 > .67f)
					{
						ee_1.lb[*i2 - 1] = 9;
						cc_1.e[*i2 - 1] = dm;
						ee_1.lb[*i1 - 1] = 25;
						cc_1.e[*i1 - 1] = arho;
						goto L40;
					}
				}
			}
			if (ee_1.lb[*i1 - 1] == 2 && ee_1.lb[*i2 - 1] == 3 || ee_1.lb[*i1 - 1] == 3 && ee_1.lb[*i2 - 1] == 2 || (ee_1.lb[*i1 - 1] == -2 && ee_1.lb[*i2 - 1] == 5 || ee_1.lb[*i1 - 1] == 5 && ee_1.lb[*i2 - 1] == -2))
			{
				if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2)
				{
					ii = *i1;
					if (x2 <= .5f)
					{
						ee_1.lb[*i1 - 1] = 6;
						cc_1.e[*i1 - 1] = dm;
						ee_1.lb[*i2 - 1] = 26;
						cc_1.e[*i2 - 1] = arho;
						goto L40;
					}
					else
					{
						ee_1.lb[*i1 - 1] = 7;
						cc_1.e[*i1 - 1] = dm;
						ee_1.lb[*i2 - 1] = 25;
						cc_1.e[*i2 - 1] = arho;
						goto L40;
					}
				}
				else
				{
					ii = *i2;
					if (x2 <= .5f)
					{
						ee_1.lb[*i2 - 1] = 6;
						cc_1.e[*i2 - 1] = dm;
						ee_1.lb[*i1 - 1] = 26;
						cc_1.e[*i1 - 1] = arho;
						goto L40;
					}
					else
					{
						ee_1.lb[*i2 - 1] = 7;
						cc_1.e[*i2 - 1] = dm;
						ee_1.lb[*i1 - 1] = 25;
						cc_1.e[*i1 - 1] = arho;
						goto L40;
					}
				}
			}
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2 && ee_1.lb[*i2 - 1] == 4 || ee_1.lb[*i1 - 1] == 4 && (i__2 = ee_1.lb[*i2 - 1],
																												abs(i__2)) == 2)
			{
				if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2)
				{
					ii = *i1;
					if (x2 <= .33f)
					{
						ee_1.lb[*i1 - 1] = 7;
						cc_1.e[*i1 - 1] = dm;
						ee_1.lb[*i2 - 1] = 26;
						cc_1.e[*i2 - 1] = arho;
						goto L40;
					}
					if (x2 > .33f && x2 <= .67f)
					{
						ee_1.lb[*i1 - 1] = 6;
						cc_1.e[*i1 - 1] = dm;
						ee_1.lb[*i2 - 1] = 27;
						cc_1.e[*i2 - 1] = arho;
						goto L40;
					}
					if (x2 > .67f)
					{
						ee_1.lb[*i1 - 1] = 8;
						cc_1.e[*i1 - 1] = dm;
						ee_1.lb[*i2 - 1] = 25;
						cc_1.e[*i2 - 1] = arho;
					}
				}
				else
				{
					ii = *i2;
					if (x2 <= .33f)
					{
						ee_1.lb[*i2 - 1] = 7;
						cc_1.e[*i2 - 1] = dm;
						ee_1.lb[*i1 - 1] = 26;
						cc_1.e[*i1 - 1] = arho;
						goto L40;
					}
					if (x2 <= .67f && x2 > .33f)
					{
						ee_1.lb[*i2 - 1] = 6;
						cc_1.e[*i2 - 1] = dm;
						ee_1.lb[*i1 - 1] = 27;
						cc_1.e[*i1 - 1] = arho;
						goto L40;
					}
					if (x2 > .67f)
					{
						ee_1.lb[*i2 - 1] = 8;
						cc_1.e[*i2 - 1] = dm;
						ee_1.lb[*i1 - 1] = 25;
						cc_1.e[*i1 - 1] = arho;
					}
				}
			}
		}
		if (*iblock == 79)
		{
			aomega = .782f;
			dmax__ = *srt - .782f - .02f;
			dm = rmass_(&dmax__, &input1_1.iseed);
			if (ee_1.lb[*i1 - 1] == 1 && ee_1.lb[*i2 - 1] == 5 || ee_1.lb[*i1 - 1] == 5 && ee_1.lb[*i2 - 1] == 1 || (ee_1.lb[*i1 - 1] == -1 && ee_1.lb[*i2 - 1] == 3 || ee_1.lb[*i1 - 1] == 3 && ee_1.lb[*i2 - 1] == -1))
			{
				if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 1)
				{
					ii = *i1;
					ee_1.lb[*i1 - 1] = 9;
					cc_1.e[*i1 - 1] = dm;
					ee_1.lb[*i2 - 1] = 28;
					cc_1.e[*i2 - 1] = aomega;
					goto L40;
				}
				else
				{
					ii = *i2;
					ee_1.lb[*i2 - 1] = 9;
					cc_1.e[*i2 - 1] = dm;
					ee_1.lb[*i1 - 1] = 28;
					cc_1.e[*i1 - 1] = aomega;
					goto L40;
				}
			}
			if (ee_1.lb[*i1 - 1] == 1 && ee_1.lb[*i2 - 1] == 3 || ee_1.lb[*i1 - 1] == 3 && ee_1.lb[*i2 - 1] == 1 || (ee_1.lb[*i1 - 1] == -1 && ee_1.lb[*i2 - 1] == 5 || ee_1.lb[*i1 - 1] == 5 && ee_1.lb[*i2 - 1] == -1))
			{
				if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 1)
				{
					ii = *i1;
					ee_1.lb[*i1 - 1] = 7;
					cc_1.e[*i1 - 1] = dm;
					ee_1.lb[*i2 - 1] = 28;
					cc_1.e[*i2 - 1] = aomega;
					goto L40;
				}
				else
				{
					ii = *i2;
					ee_1.lb[*i2 - 1] = 7;
					cc_1.e[*i2 - 1] = dm;
					ee_1.lb[*i1 - 1] = 28;
					cc_1.e[*i1 - 1] = aomega;
					goto L40;
				}
			}
			if (ee_1.lb[*i1 - 1] == 2 && ee_1.lb[*i2 - 1] == 5 || ee_1.lb[*i1 - 1] == 5 && ee_1.lb[*i2 - 1] == 2 || (ee_1.lb[*i1 - 1] == -2 && ee_1.lb[*i2 - 1] == 3 || ee_1.lb[*i1 - 1] == 3 && ee_1.lb[*i2 - 1] == -2))
			{
				if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2)
				{
					ii = *i1;
					ee_1.lb[*i1 - 1] = 8;
					cc_1.e[*i1 - 1] = dm;
					ee_1.lb[*i2 - 1] = 28;
					cc_1.e[*i2 - 1] = aomega;
					goto L40;
				}
				else
				{
					ii = *i2;
					ee_1.lb[*i2 - 1] = 8;
					cc_1.e[*i2 - 1] = dm;
					ee_1.lb[*i1 - 1] = 28;
					cc_1.e[*i1 - 1] = aomega;
					goto L40;
				}
			}
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 1 && ee_1.lb[*i2 - 1] == 4 || ee_1.lb[*i1 - 1] == 4 && (i__2 = ee_1.lb[*i2 - 1],
																												abs(i__2)) == 1)
			{
				if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 1)
				{
					ii = *i1;
					ee_1.lb[*i1 - 1] = 8;
					cc_1.e[*i1 - 1] = dm;
					ee_1.lb[*i2 - 1] = 28;
					cc_1.e[*i2 - 1] = aomega;
					goto L40;
				}
				else
				{
					ii = *i2;
					ee_1.lb[*i2 - 1] = 8;
					cc_1.e[*i2 - 1] = dm;
					ee_1.lb[*i1 - 1] = 28;
					cc_1.e[*i1 - 1] = aomega;
					goto L40;
				}
			}
			if (ee_1.lb[*i1 - 1] == 2 && ee_1.lb[*i2 - 1] == 3 || ee_1.lb[*i1 - 1] == 3 && ee_1.lb[*i2 - 1] == 2 || (ee_1.lb[*i1 - 1] == -2 && ee_1.lb[*i2 - 1] == 5 || ee_1.lb[*i1 - 1] == 5 && ee_1.lb[*i2 - 1] == -2))
			{
				if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2)
				{
					ii = *i1;
					ee_1.lb[*i1 - 1] = 6;
					cc_1.e[*i1 - 1] = dm;
					ee_1.lb[*i2 - 1] = 28;
					cc_1.e[*i2 - 1] = aomega;
					goto L40;
				}
				else
				{
					ii = *i2;
					ee_1.lb[*i2 - 1] = 6;
					cc_1.e[*i2 - 1] = dm;
					ee_1.lb[*i1 - 1] = 28;
					cc_1.e[*i1 - 1] = aomega;
				}
			}
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2 && ee_1.lb[*i2 - 1] == 4 || ee_1.lb[*i1 - 1] == 4 && (i__2 = ee_1.lb[*i2 - 1],
																												abs(i__2)) == 2)
			{
				if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2)
				{
					ii = *i1;
					ee_1.lb[*i1 - 1] = 7;
					cc_1.e[*i1 - 1] = dm;
					ee_1.lb[*i2 - 1] = 28;
					cc_1.e[*i2 - 1] = aomega;
					goto L40;
				}
				else
				{
					ii = *i2;
					ee_1.lb[*i2 - 1] = 7;
					cc_1.e[*i2 - 1] = dm;
					ee_1.lb[*i1 - 1] = 26;
					cc_1.e[*i1 - 1] = arho;
					goto L40;
				}
			}
		}
	L40:
		em1 = cc_1.e[*i1 - 1];
		em2 = cc_1.e[*i2 - 1];
		if (ianti == 1 && ee_1.lb[*i1 - 1] >= 1 && ee_1.lb[*i2 - 1] >= 1)
		{
			ee_1.lb[ii - 1] = -ee_1.lb[ii - 1];
			jj = *i2;
			if (ii == *i2)
			{
				jj = *i1;
			}
			if (*iblock == 77)
			{
				if (ee_1.lb[jj - 1] == 3)
				{
					ee_1.lb[jj - 1] = 5;
				}
				else if (ee_1.lb[jj - 1] == 5)
				{
					ee_1.lb[jj - 1] = 3;
				}
			}
			else if (*iblock == 78)
			{
				if (ee_1.lb[jj - 1] == 25)
				{
					ee_1.lb[jj - 1] = 27;
				}
				else if (ee_1.lb[jj - 1] == 27)
				{
					ee_1.lb[jj - 1] = 25;
				}
			}
		}
	}
L50:
	r__2 = *srt;
	r__3 = em1;
	r__4 = em2;
	r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
	r__5 = em1 * em2;
	pr2 = r__1 * r__1 - r__5 * r__5 * 4.f;
	if (pr2 <= 0.f)
	{
		pr2 = 1e-8f;
	}
	pr = sqrt(pr2) / (*srt * 2.f);
	xptr = pr * .33f;
	cc1 = ptr_(&xptr, &input1_1.iseed);
	r__1 = pr;
	r__2 = cc1;
	c1 = sqrt(r__1 * r__1 - r__2 * r__2) / pr;
	t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
	r__1 = c1;
	s1 = sqrt(1.f - r__1 * r__1);
	ct1 = cos(t1);
	st1 = sin(t1);
	*pz = pr * c1;
	*px = pr * s1 * ct1;
	*py = pr * s1 * st1;
	rotate_(&px0, &py0, &pz0, px, py, pz);
	return 0;
}

int cren_(float *px, float *py, float *pz, float *srt, int *i1, int *i2, int *iblock)
{
	float r__1, r__2, r__3, r__4, r__5;

	double sqrt(double), cos(double), sin(double);

	static float c1, s1, t1, pr, em1, em2, ct1, pr2, px0, py0, pz0, st1;
	static int ntag, kaonc, ianti;
	extern double pnlka_(float *), pnska_(float *), ranart_(int *);

	px0 = *px;
	py0 = *py;
	pz0 = *pz;
	ntag = 0;
	*iblock = 7;
	ianti = 0;
	if (ee_1.lb[*i1 - 1] < 0 || ee_1.lb[*i2 - 1] < 0)
	{
		ianti = 1;
		*iblock = -7;
	}
	kaonc = 0;
	if (pnlka_(srt) / (pnlka_(srt) + pnska_(srt)) > ranart_(&rndf77_1.nseed))
	{
		kaonc = 1;
	}
	if (cc_1.e[*i1 - 1] <= .6f)
	{
		ee_1.lb[*i1 - 1] = 23;
		cc_1.e[*i1 - 1] = .498f;
		if (kaonc == 1)
		{
			ee_1.lb[*i2 - 1] = 14;
			cc_1.e[*i2 - 1] = 1.1157f;
		}
		else
		{
			ee_1.lb[*i2 - 1] = (int)(ranart_(&rndf77_1.nseed) * 3) + 15;
			cc_1.e[*i2 - 1] = 1.1974f;
		}
		if (ianti == 1)
		{
			ee_1.lb[*i1 - 1] = 21;
			ee_1.lb[*i2 - 1] = -ee_1.lb[*i2 - 1];
		}
	}
	else
	{
		ee_1.lb[*i2 - 1] = 23;
		cc_1.e[*i2 - 1] = .498f;
		if (kaonc == 1)
		{
			ee_1.lb[*i1 - 1] = 14;
			cc_1.e[*i1 - 1] = 1.1157f;
		}
		else
		{
			ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 3) + 15;
			cc_1.e[*i1 - 1] = 1.1974f;
		}
		if (ianti == 1)
		{
			ee_1.lb[*i2 - 1] = 21;
			ee_1.lb[*i1 - 1] = -ee_1.lb[*i1 - 1];
		}
	}
	em1 = cc_1.e[*i1 - 1];
	em2 = cc_1.e[*i2 - 1];
	r__2 = *srt;
	r__3 = em1;
	r__4 = em2;
	r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
	r__5 = em1 * em2;
	pr2 = r__1 * r__1 - r__5 * r__5 * 4.f;
	if (pr2 <= 0.f)
	{
		pr2 = 1e-9f;
	}
	pr = sqrt(pr2) / (*srt * 2.f);
	c1 = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
	t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
	r__1 = c1;
	s1 = sqrt(1.f - r__1 * r__1);
	ct1 = cos(t1);
	st1 = sin(t1);
	*pz = pr * c1;
	*px = pr * s1 * ct1;
	*py = pr * s1 * st1;
	return 0;
}

int crdir_(float *px, float *py, float *pz, float *srt, int *i1, int *i2, int *iblock)
{
	float r__1, r__2, r__3, r__4, r__5;

	double sqrt(double), cos(double), sin(double);

	static float c1, s1, t1, pr, cc1, em1, em2, ct1, pr2, px0, py0, pz0, st1;
	extern double ptr_(float *, int *);
	static int ntag;
	static float xptr;
	extern double ranart_(int *);
	extern int rotate_(float *, float *, float *, float *, float *, float *);

	px0 = *px;
	py0 = *py;
	pz0 = *pz;
	*iblock = 999;
	ntag = 0;
	em1 = cc_1.e[*i1 - 1];
	em2 = cc_1.e[*i2 - 1];
	r__2 = *srt;
	r__3 = em1;
	r__4 = em2;
	r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
	r__5 = em1 * em2;
	pr2 = r__1 * r__1 - r__5 * r__5 * 4.f;
	if (pr2 <= 0.f)
	{
		pr2 = 1e-9f;
	}
	pr = sqrt(pr2) / (*srt * 2.f);
	xptr = pr * .33f;
	cc1 = ptr_(&xptr, &input1_1.iseed);
	r__1 = pr;
	r__2 = cc1;
	c1 = sqrt(r__1 * r__1 - r__2 * r__2) / pr;
	t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
	r__1 = c1;
	s1 = sqrt(1.f - r__1 * r__1);
	ct1 = cos(t1);
	st1 = sin(t1);
	*pz = pr * c1;
	*px = pr * s1 * ct1;
	*py = pr * s1 * st1;
	rotate_(&px0, &py0, &pz0, px, py, pz);
	return 0;
}

int crpd_(float *px, float *py, float *pz, float *srt, int *i1, int *i2, int *iblock, float *xkaon0, float *xkaon, float *xphi, float *xphin)
{
	int i__1, i__2, i__3, i__4;
	float r__1, r__2, r__3, r__4, r__5;

	double sqrt(double), cos(double), sin(double);

	static float c1, s1, t1, x1, x2;
	static int ii, jj;
	static float pr, cc1, em1, em2, ct1, pr2, px0, py0, pz0, st1;
	extern double ptr_(float *, int *);
	static int ntag;
	static float xptr;
	static int kaonc, ianti;
	extern double pnlka_(float *), pnska_(float *), ranart_(int *);
	extern int rotate_(float *, float *, float *, float *, float *, float *);

	px0 = *px;
	py0 = *py;
	pz0 = *pz;
	*iblock = 1;
	x1 = ranart_(&rndf77_1.nseed);
	ianti = 0;
	if (ee_1.lb[*i1 - 1] < 0 || ee_1.lb[*i2 - 1] < 0)
	{
		ianti = 1;
	}
	if (*xkaon0 / (*xkaon + *xphi) >= x1)
	{
		*iblock = 7;
		if (ianti == 1)
		{
			*iblock = -7;
		}
		ntag = 0;
		kaonc = 0;
		if (pnlka_(srt) / (pnlka_(srt) + pnska_(srt)) > ranart_(&rndf77_1.nseed))
		{
			kaonc = 1;
		}
		if (cc_1.e[*i1 - 1] <= .2f)
		{
			ee_1.lb[*i1 - 1] = 23;
			cc_1.e[*i1 - 1] = .498f;
			if (kaonc == 1)
			{
				ee_1.lb[*i2 - 1] = 14;
				cc_1.e[*i2 - 1] = 1.1157f;
			}
			else
			{
				ee_1.lb[*i2 - 1] = (int)(ranart_(&rndf77_1.nseed) * 3) +
								   15;
				cc_1.e[*i2 - 1] = 1.1974f;
			}
			if (ianti == 1)
			{
				ee_1.lb[*i1 - 1] = 21;
				ee_1.lb[*i2 - 1] = -ee_1.lb[*i2 - 1];
			}
		}
		else
		{
			ee_1.lb[*i2 - 1] = 23;
			cc_1.e[*i2 - 1] = .498f;
			if (kaonc == 1)
			{
				ee_1.lb[*i1 - 1] = 14;
				cc_1.e[*i1 - 1] = 1.1157f;
			}
			else
			{
				ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 3) +
								   15;
				cc_1.e[*i1 - 1] = 1.1974f;
			}
			if (ianti == 1)
			{
				ee_1.lb[*i2 - 1] = 21;
				ee_1.lb[*i1 - 1] = -ee_1.lb[*i1 - 1];
			}
		}
		em1 = cc_1.e[*i1 - 1];
		em2 = cc_1.e[*i2 - 1];
		goto L50;
	}
	else if (*xphi / (*xkaon + *xphi) >= x1)
	{
		*iblock = 222;
		if (*xphin / *xphi >= ranart_(&rndf77_1.nseed))
		{
			ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 2) + 1;
			cc_1.e[*i1 - 1] = .939457f;
		}
		else
		{
			ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 4) + 6;
			cc_1.e[*i1 - 1] = 1.232f;
		}
		if (ianti == 1)
		{
			ee_1.lb[*i1 - 1] = -ee_1.lb[*i1 - 1];
		}
		ee_1.lb[*i2 - 1] = 29;
		cc_1.e[*i2 - 1] = 1.02f;
		em1 = cc_1.e[*i1 - 1];
		em2 = cc_1.e[*i2 - 1];
		goto L50;
	}
	else
	{
		x2 = ranart_(&rndf77_1.nseed);
		*iblock = 80;
		ntag = 0;
		if (ee_1.lb[*i1 - 1] == 8 && ee_1.lb[*i2 - 1] == 5 || ee_1.lb[*i1 - 1] == 5 && ee_1.lb[*i2 - 1] == 8 || (ee_1.lb[*i1 - 1] == -8 && ee_1.lb[*i2 - 1] == 3 || ee_1.lb[*i1 - 1] == 3 && ee_1.lb[*i2 - 1] == -8))
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 8)
			{
				ii = *i1;
				ee_1.lb[*i1 - 1] = 1;
				cc_1.e[*i1 - 1] = .939457f;
				ee_1.lb[*i2 - 1] = 5;
				cc_1.e[*i2 - 1] = .13496f;
				goto L40;
			}
			else
			{
				ii = *i2;
				ee_1.lb[*i2 - 1] = 1;
				cc_1.e[*i2 - 1] = .939457f;
				ee_1.lb[*i1 - 1] = 5;
				cc_1.e[*i1 - 1] = .13496f;
				goto L40;
			}
		}

		if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 7 && ee_1.lb[*i2 - 1] ==
															 4 ||
			ee_1.lb[*i1 - 1] == 4 && (i__2 = ee_1.lb[*i2 - 1], abs(
																   i__2)) == 7)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 7)
			{
				ii = *i1;
				if (x2 <= .5f)
				{
					ee_1.lb[*i1 - 1] = 2;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 4;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i1 - 1] = 1;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 3;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
			}
			else
			{
				ii = *i2;
				if (x2 <= .5f)
				{
					ee_1.lb[*i2 - 1] = 2;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 4;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i2 - 1] = 1;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 3;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
			}
		}
		if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 8 && ee_1.lb[*i2 - 1] ==
															 4 ||
			ee_1.lb[*i1 - 1] == 4 && (i__2 = ee_1.lb[*i2 - 1], abs(
																   i__2)) == 8)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 8)
			{
				ii = *i1;
				if (x2 <= .5f)
				{
					ee_1.lb[*i1 - 1] = 2;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 5;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i1 - 1] = 1;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 4;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
			}
			else
			{
				ii = *i2;
				if (x2 <= .5f)
				{
					ee_1.lb[*i2 - 1] = 2;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 5;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i2 - 1] = 1;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 4;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
			}
		}
		if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 6 && ee_1.lb[*i2 - 1] ==
															 4 ||
			ee_1.lb[*i1 - 1] == 4 && (i__2 = ee_1.lb[*i2 - 1], abs(
																   i__2)) == 6)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 6)
			{
				ii = *i1;
				ee_1.lb[*i1 - 1] = 2;
				cc_1.e[*i1 - 1] = .939457f;
				ee_1.lb[*i2 - 1] = 3;
				cc_1.e[*i2 - 1] = .13496f;
				goto L40;
			}
			else
			{
				ii = *i2;
				ee_1.lb[*i2 - 1] = 2;
				cc_1.e[*i2 - 1] = .939457f;
				ee_1.lb[*i1 - 1] = 3;
				cc_1.e[*i1 - 1] = .13496f;
				goto L40;
			}
		}
		if (ee_1.lb[*i1 - 1] == 8 && ee_1.lb[*i2 - 1] == 3 || ee_1.lb[*i1 - 1] == 3 && ee_1.lb[*i2 - 1] == 8 || (ee_1.lb[*i1 - 1] == -8 && ee_1.lb[*i2 - 1] == 5 || ee_1.lb[*i1 - 1] == 5 && ee_1.lb[*i2 - 1] == -8))
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 8)
			{
				ii = *i1;
				if (x2 <= .5f)
				{
					ee_1.lb[*i1 - 1] = 2;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 4;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i1 - 1] = 1;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 3;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
			}
			else
			{
				ii = *i2;
				if (x2 <= .5f)
				{
					ee_1.lb[*i2 - 1] = 2;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 4;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i2 - 1] = 1;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 3;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
			}
		}
		if (ee_1.lb[*i1 - 1] == 7 && ee_1.lb[*i2 - 1] == 5 || ee_1.lb[*i1 - 1] == 5 && ee_1.lb[*i2 - 1] == 7 || (ee_1.lb[*i1 - 1] == -7 && ee_1.lb[*i2 - 1] == 3 || ee_1.lb[*i1 - 1] == 3 && ee_1.lb[*i2 - 1] == -7))
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 7)
			{
				ii = *i1;
				if (x2 <= .5f)
				{
					ee_1.lb[*i1 - 1] = 2;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 5;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i1 - 1] = 1;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 4;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
			}
			else
			{
				ii = *i2;
				if (x2 <= .5f)
				{
					ee_1.lb[*i2 - 1] = 2;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 5;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i2 - 1] = 1;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 4;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
			}
		}
		if (ee_1.lb[*i1 - 1] == 7 && ee_1.lb[*i2 - 1] == 3 || ee_1.lb[*i1 - 1] == 3 && ee_1.lb[*i2 - 1] == 7 || (ee_1.lb[*i1 - 1] == -7 && ee_1.lb[*i2 - 1] == 5 || ee_1.lb[*i1 - 1] == 5 && ee_1.lb[*i2 - 1] == -7))
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 7)
			{
				ii = *i1;
				ee_1.lb[*i1 - 1] = 2;
				cc_1.e[*i1 - 1] = .939457f;
				ee_1.lb[*i2 - 1] = 3;
				cc_1.e[*i2 - 1] = .13496f;
				goto L40;
			}
			else
			{
				ii = *i2;
				ee_1.lb[*i2 - 1] = 2;
				cc_1.e[*i2 - 1] = .939457f;
				ee_1.lb[*i1 - 1] = 3;
				cc_1.e[*i1 - 1] = .13496f;
				goto L40;
			}
		}
		if (ee_1.lb[*i1 - 1] == 6 && ee_1.lb[*i2 - 1] == 5 || ee_1.lb[*i1 - 1] == 5 && ee_1.lb[*i2 - 1] == 6 || (ee_1.lb[*i1 - 1] == -6 && ee_1.lb[*i2 - 1] == 3 || ee_1.lb[*i1 - 1] == 3 && ee_1.lb[*i2 - 1] == -6))
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 6)
			{
				ii = *i1;
				if (x2 <= .5f)
				{
					ee_1.lb[*i1 - 1] = 2;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 4;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i1 - 1] = 1;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 3;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
			}
			else
			{
				ii = *i2;
				if (x2 <= .5f)
				{
					ee_1.lb[*i2 - 1] = 2;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 4;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i2 - 1] = 1;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 3;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
			}
		}

		if (ee_1.lb[*i1 - 1] == 9 && ee_1.lb[*i2 - 1] == 3 || ee_1.lb[*i1 - 1] == 3 && ee_1.lb[*i2 - 1] == 9 || (ee_1.lb[*i1 - 1] == -9 && ee_1.lb[*i2 - 1] == 5 || ee_1.lb[*i1 - 1] == 5 && ee_1.lb[*i2 - 1] == -9))
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 9)
			{
				ii = *i1;
				if (x2 <= .5f)
				{
					ee_1.lb[*i1 - 1] = 2;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 5;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i1 - 1] = 1;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 4;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
			}
			else
			{
				ii = *i2;
				if (x2 <= .5f)
				{
					ee_1.lb[*i2 - 1] = 2;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 5;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i2 - 1] = 1;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 4;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
			}
		}
		if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 9 && ee_1.lb[*i2 - 1] ==
															 4 ||
			ee_1.lb[*i1 - 1] == 4 && (i__2 = ee_1.lb[*i2 - 1], abs(
																   i__2)) == 9)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 9)
			{
				ii = *i1;
				ee_1.lb[*i1 - 1] = 1;
				cc_1.e[*i1 - 1] = .939457f;
				ee_1.lb[*i2 - 1] = 5;
				cc_1.e[*i2 - 1] = .13496f;
				goto L40;
			}
			else
			{
				ii = *i2;
				ee_1.lb[*i2 - 1] = 1;
				cc_1.e[*i2 - 1] = .939457f;
				ee_1.lb[*i1 - 1] = 5;
				cc_1.e[*i1 - 1] = .13496f;
				goto L40;
			}
		}
		if (ee_1.lb[*i1 - 1] == 11 && ee_1.lb[*i2 - 1] == 5 || ee_1.lb[*i1 - 1] == 5 && ee_1.lb[*i2 - 1] == 11 || ee_1.lb[*i1 - 1] == 13 && ee_1.lb[*i2 - 1] == 5 || ee_1.lb[*i1 - 1] == 5 && ee_1.lb[*i2 - 1] == 13 || (ee_1.lb[*i1 - 1] == -11 && ee_1.lb[*i2 - 1] == 3 || ee_1.lb[*i1 - 1] == 3 && ee_1.lb[*i2 - 1] == -11 || ee_1.lb[*i1 - 1] == -13 && ee_1.lb[*i2 - 1] == 3 || ee_1.lb[*i1 - 1] == 3 && ee_1.lb[*i2 - 1] == -13))
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 11 || (i__2 = ee_1.lb[*i1 - 1], abs(i__2)) == 13)
			{
				ii = *i1;
				ee_1.lb[*i1 - 1] = 1;
				cc_1.e[*i1 - 1] = .939457f;
				ee_1.lb[*i2 - 1] = 5;
				cc_1.e[*i2 - 1] = .13496f;
				goto L40;
			}
			else
			{
				ii = *i2;
				ee_1.lb[*i2 - 1] = 1;
				cc_1.e[*i2 - 1] = .939457f;
				ee_1.lb[*i1 - 1] = 5;
				cc_1.e[*i1 - 1] = .13496f;
				goto L40;
			}
		}
		if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 10 && ee_1.lb[*i2 - 1] ==
															  4 ||
			ee_1.lb[*i1 - 1] == 4 && (i__2 = ee_1.lb[*i2 - 1], abs(
																   i__2)) == 10 ||
			ee_1.lb[*i1 - 1] == 4 && (i__3 = ee_1.lb[*i2 - 1], abs(i__3)) == 12 || ee_1.lb[*i2 - 1] == 4 && (i__4 = ee_1.lb[*i1 - 1], abs(i__4)) == 12)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 10 || (i__2 = ee_1.lb[*i1 - 1], abs(i__2)) == 12)
			{
				ii = *i1;
				if (x2 <= .5f)
				{
					ee_1.lb[*i1 - 1] = 2;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 4;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i1 - 1] = 1;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 3;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
			}
			else
			{
				ii = *i2;
				if (x2 <= .5f)
				{
					ee_1.lb[*i2 - 1] = 2;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 4;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i2 - 1] = 1;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 3;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
			}
		}
		if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 11 && ee_1.lb[*i2 - 1] ==
															  4 ||
			ee_1.lb[*i1 - 1] == 4 && (i__2 = ee_1.lb[*i2 - 1], abs(
																   i__2)) == 11 ||
			ee_1.lb[*i1 - 1] == 4 && (i__3 = ee_1.lb[*i2 - 1], abs(i__3)) == 13 || ee_1.lb[*i2 - 1] == 4 && (i__4 = ee_1.lb[*i1 - 1], abs(i__4)) == 13)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 11 || (i__2 = ee_1.lb[*i1 - 1], abs(i__2)) == 13)
			{
				ii = *i1;
				if (x2 <= .5f)
				{
					ee_1.lb[*i1 - 1] = 2;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 5;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i1 - 1] = 1;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 4;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
			}
			else
			{
				ii = *i2;
				if (x2 <= .5f)
				{
					ee_1.lb[*i2 - 1] = 2;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 5;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i2 - 1] = 1;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 4;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
			}
		}
		if (ee_1.lb[*i1 - 1] == 11 && ee_1.lb[*i2 - 1] == 3 || ee_1.lb[*i1 - 1] == 3 && ee_1.lb[*i2 - 1] == 11 || ee_1.lb[*i1 - 1] == 3 && ee_1.lb[*i2 - 1] == 13 || ee_1.lb[*i2 - 1] == 3 && ee_1.lb[*i1 - 1] == 13 || (ee_1.lb[*i1 - 1] == -11 && ee_1.lb[*i2 - 1] == 5 || ee_1.lb[*i1 - 1] == 5 && ee_1.lb[*i2 - 1] == -11 || ee_1.lb[*i1 - 1] == 5 && ee_1.lb[*i2 - 1] == -13 || ee_1.lb[*i2 - 1] == 5 && ee_1.lb[*i1 - 1] == -13))
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 11 || (i__2 = ee_1.lb[*i1 - 1], abs(i__2)) == 13)
			{
				ii = *i1;
				if (x2 <= .5f)
				{
					ee_1.lb[*i1 - 1] = 2;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 4;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i1 - 1] = 1;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 3;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
			}
			else
			{
				ii = *i2;
				if (x2 <= .5f)
				{
					ee_1.lb[*i2 - 1] = 2;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 4;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i2 - 1] = 1;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 3;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
			}
		}
		if (ee_1.lb[*i1 - 1] == 10 && ee_1.lb[*i2 - 1] == 5 || ee_1.lb[*i1 - 1] == 5 && ee_1.lb[*i2 - 1] == 10 || ee_1.lb[*i1 - 1] == 12 && ee_1.lb[*i2 - 1] == 5 || ee_1.lb[*i1 - 1] == 5 && ee_1.lb[*i2 - 1] == 12 || (ee_1.lb[*i1 - 1] == -10 && ee_1.lb[*i2 - 1] == 3 || ee_1.lb[*i1 - 1] == 3 && ee_1.lb[*i2 - 1] == -10 || ee_1.lb[*i1 - 1] == -12 && ee_1.lb[*i2 - 1] == 3 || ee_1.lb[*i1 - 1] == 3 && ee_1.lb[*i2 - 1] == -12))
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 10 || (i__2 = ee_1.lb[*i1 - 1], abs(i__2)) == 12)
			{
				ii = *i1;
				if (x2 <= .5f)
				{
					ee_1.lb[*i1 - 1] = 2;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 5;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i1 - 1] = 1;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 4;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
			}
			else
			{
				ii = *i2;
				if (x2 <= .5f)
				{
					ee_1.lb[*i2 - 1] = 2;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 5;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i2 - 1] = 1;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 4;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
			}
		}
		if (ee_1.lb[*i1 - 1] == 10 && ee_1.lb[*i2 - 1] == 3 || ee_1.lb[*i1 - 1] == 3 && ee_1.lb[*i2 - 1] == 10 || ee_1.lb[*i1 - 1] == 3 && ee_1.lb[*i2 - 1] == 12 || ee_1.lb[*i1 - 1] == 12 && ee_1.lb[*i2 - 1] == 3 || (ee_1.lb[*i1 - 1] == -10 && ee_1.lb[*i2 - 1] == 5 || ee_1.lb[*i1 - 1] == 5 && ee_1.lb[*i2 - 1] == -10 || ee_1.lb[*i1 - 1] == 5 && ee_1.lb[*i2 - 1] == -12 || ee_1.lb[*i1 - 1] == -12 && ee_1.lb[*i2 - 1] == 5))
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 10 || (i__2 = ee_1.lb[*i1 - 1], abs(i__2)) == 12)
			{
				ii = *i1;
				ee_1.lb[*i1 - 1] = 2;
				cc_1.e[*i1 - 1] = .939457f;
				ee_1.lb[*i2 - 1] = 3;
				cc_1.e[*i2 - 1] = .13496f;
				goto L40;
			}
			else
			{
				ii = *i2;
				ee_1.lb[*i2 - 1] = 2;
				cc_1.e[*i2 - 1] = .939457f;
				ee_1.lb[*i1 - 1] = 3;
				cc_1.e[*i1 - 1] = .13496f;
				goto L40;
			}
		}
	L40:
		em1 = cc_1.e[*i1 - 1];
		em2 = cc_1.e[*i2 - 1];
		if (ianti == 1 && ee_1.lb[*i1 - 1] >= 1 && ee_1.lb[*i2 - 1] >= 1)
		{
			ee_1.lb[ii - 1] = -ee_1.lb[ii - 1];
			jj = *i2;
			if (ii == *i2)
			{
				jj = *i1;
			}
			if (ee_1.lb[jj - 1] == 3)
			{
				ee_1.lb[jj - 1] = 5;
			}
			else if (ee_1.lb[jj - 1] == 5)
			{
				ee_1.lb[jj - 1] = 3;
			}
		}
	}
L50:
	r__2 = *srt;
	r__3 = em1;
	r__4 = em2;
	r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
	r__5 = em1 * em2;
	pr2 = r__1 * r__1 - r__5 * r__5 * 4.f;
	if (pr2 <= 0.f)
	{
		pr2 = 1e-9f;
	}
	pr = sqrt(pr2) / (*srt * 2.f);
	xptr = pr * .33f;
	cc1 = ptr_(&xptr, &input1_1.iseed);
	r__1 = pr;
	r__2 = cc1;
	c1 = sqrt(r__1 * r__1 - r__2 * r__2) / pr;
	t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
	r__1 = c1;
	s1 = sqrt(1.f - r__1 * r__1);
	ct1 = cos(t1);
	st1 = sin(t1);
	*pz = pr * c1;
	*px = pr * s1 * ct1;
	*py = pr * s1 * st1;
	rotate_(&px0, &py0, &pz0, px, py, pz);
	return 0;
}

int crrd_(float *px, float *py, float *pz, float *srt, int *i1, int *i2, int *iblock, float *xkaon0, float *xkaon, float *xphi, float *xphin)
{
	int i__1, i__2, i__3, i__4;
	float r__1, r__2, r__3, r__4, r__5;

	double sqrt(double), cos(double), sin(double);

	static float c1, s1, t1, x1, x2;
	static int ii, jj;
	static float pr, cc1, em1, em2, ct1, pr2, px0, py0, pz0, st1;
	extern double ptr_(float *, int *);
	static int ntag;
	static float xptr;
	static int kaonc, ianti;
	extern double pnlka_(float *), pnska_(float *), ranart_(int *);
	extern int rotate_(float *, float *, float *, float *, float *, float *);

	px0 = *px;
	py0 = *py;
	pz0 = *pz;
	*iblock = 1;
	ianti = 0;
	if (ee_1.lb[*i1 - 1] < 0 || ee_1.lb[*i2 - 1] < 0)
	{
		ianti = 1;
	}
	x1 = ranart_(&rndf77_1.nseed);
	if (*xkaon0 / (*xkaon + *xphi) >= x1)
	{
		*iblock = 7;
		if (ianti == 1)
		{
			*iblock = -7;
		}
		ntag = 0;
		kaonc = 0;
		if (pnlka_(srt) / (pnlka_(srt) + pnska_(srt)) > ranart_(&rndf77_1.nseed))
		{
			kaonc = 1;
		}
		if (cc_1.e[*i1 - 1] <= .92f)
		{
			ee_1.lb[*i1 - 1] = 23;
			cc_1.e[*i1 - 1] = .498f;
			if (kaonc == 1)
			{
				ee_1.lb[*i2 - 1] = 14;
				cc_1.e[*i2 - 1] = 1.1157f;
			}
			else
			{
				ee_1.lb[*i2 - 1] = (int)(ranart_(&rndf77_1.nseed) * 3) +
								   15;
				cc_1.e[*i2 - 1] = 1.1974f;
			}
			if (ianti == 1)
			{
				ee_1.lb[*i1 - 1] = 21;
				ee_1.lb[*i2 - 1] = -ee_1.lb[*i2 - 1];
			}
		}
		else
		{
			ee_1.lb[*i2 - 1] = 23;
			cc_1.e[*i2 - 1] = .498f;
			if (kaonc == 1)
			{
				ee_1.lb[*i1 - 1] = 14;
				cc_1.e[*i1 - 1] = 1.1157f;
			}
			else
			{
				ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 3) +
								   15;
				cc_1.e[*i1 - 1] = 1.1974f;
			}
			if (ianti == 1)
			{
				ee_1.lb[*i2 - 1] = 21;
				ee_1.lb[*i1 - 1] = -ee_1.lb[*i1 - 1];
			}
		}
		em1 = cc_1.e[*i1 - 1];
		em2 = cc_1.e[*i2 - 1];
		goto L50;
	}
	else if (*xphi / (*xkaon + *xphi) >= x1)
	{
		*iblock = 222;
		if (*xphin / *xphi >= ranart_(&rndf77_1.nseed))
		{
			ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 2) + 1;
			cc_1.e[*i1 - 1] = .939457f;
		}
		else
		{
			ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 4) + 6;
			cc_1.e[*i1 - 1] = 1.232f;
		}
		if (ianti == 1)
		{
			ee_1.lb[*i1 - 1] = -ee_1.lb[*i1 - 1];
		}
		ee_1.lb[*i2 - 1] = 29;
		cc_1.e[*i2 - 1] = 1.02f;
		em1 = cc_1.e[*i1 - 1];
		em2 = cc_1.e[*i2 - 1];
		goto L50;
	}
	else
	{
		x2 = ranart_(&rndf77_1.nseed);
		*iblock = 81;
		ntag = 0;
		if (ee_1.lb[*i1 - 1] == 28 || ee_1.lb[*i2 - 1] == 28)
		{
			goto L60;
		}
		if (ee_1.lb[*i1 - 1] == 8 && ee_1.lb[*i2 - 1] == 27 || ee_1.lb[*i1 - 1] == 27 && ee_1.lb[*i2 - 1] == 8 || (ee_1.lb[*i1 - 1] == -8 && ee_1.lb[*i2 - 1] == 25 || ee_1.lb[*i1 - 1] == 25 && ee_1.lb[*i2 - 1] == -8))
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 8)
			{
				ii = *i1;
				ee_1.lb[*i1 - 1] = 1;
				cc_1.e[*i1 - 1] = .939457f;
				ee_1.lb[*i2 - 1] = 5;
				cc_1.e[*i2 - 1] = .13496f;
				goto L40;
			}
			else
			{
				ii = *i2;
				ee_1.lb[*i2 - 1] = 1;
				cc_1.e[*i2 - 1] = .939457f;
				ee_1.lb[*i1 - 1] = 5;
				cc_1.e[*i1 - 1] = .13496f;
				goto L40;
			}
		}
		if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 7 && ee_1.lb[*i2 - 1] ==
															 26 ||
			ee_1.lb[*i1 - 1] == 26 && (i__2 = ee_1.lb[*i2 - 1], abs(
																	i__2)) == 7)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 7)
			{
				ii = *i1;
				if (x2 <= .5f)
				{
					ee_1.lb[*i1 - 1] = 2;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 4;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i1 - 1] = 1;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 3;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
			}
			else
			{
				ii = *i2;
				if (x2 <= .5f)
				{
					ee_1.lb[*i2 - 1] = 2;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 4;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i2 - 1] = 1;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 3;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
			}
		}
		if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 8 && ee_1.lb[*i2 - 1] ==
															 26 ||
			ee_1.lb[*i1 - 1] == 26 && (i__2 = ee_1.lb[*i2 - 1], abs(
																	i__2)) == 8)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 8)
			{
				ii = *i1;
				if (x2 <= .5f)
				{
					ee_1.lb[*i1 - 1] = 2;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 5;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i1 - 1] = 1;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 4;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
			}
			else
			{
				ii = *i2;
				if (x2 <= .5f)
				{
					ee_1.lb[*i2 - 1] = 2;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 5;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i2 - 1] = 1;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 4;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
			}
		}
		if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 6 && ee_1.lb[*i2 - 1] ==
															 26 ||
			ee_1.lb[*i1 - 1] == 26 && (i__2 = ee_1.lb[*i2 - 1], abs(
																	i__2)) == 6)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 6)
			{
				ii = *i1;
				ee_1.lb[*i1 - 1] = 2;
				cc_1.e[*i1 - 1] = .939457f;
				ee_1.lb[*i2 - 1] = 3;
				cc_1.e[*i2 - 1] = .13496f;
				goto L40;
			}
			else
			{
				ii = *i2;
				ee_1.lb[*i2 - 1] = 2;
				cc_1.e[*i2 - 1] = .939457f;
				ee_1.lb[*i1 - 1] = 3;
				cc_1.e[*i1 - 1] = .13496f;
				goto L40;
			}
		}
		if (ee_1.lb[*i1 - 1] == 8 && ee_1.lb[*i2 - 1] == 25 || ee_1.lb[*i1 - 1] == 25 && ee_1.lb[*i2 - 1] == 8 || (ee_1.lb[*i1 - 1] == -8 && ee_1.lb[*i2 - 1] == 27 || ee_1.lb[*i1 - 1] == 27 && ee_1.lb[*i2 - 1] == -8))
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 8)
			{
				ii = *i1;
				if (x2 <= .5f)
				{
					ee_1.lb[*i1 - 1] = 2;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 4;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i1 - 1] = 1;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 3;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
			}
			else
			{
				ii = *i2;
				if (x2 <= .5f)
				{
					ee_1.lb[*i2 - 1] = 2;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 4;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i2 - 1] = 1;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 3;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
			}
		}
		if (ee_1.lb[*i1 - 1] == 7 && ee_1.lb[*i2 - 1] == 27 || ee_1.lb[*i1 - 1] == 27 && ee_1.lb[*i2 - 1] == 7 || (ee_1.lb[*i1 - 1] == -7 && ee_1.lb[*i2 - 1] == 25 || ee_1.lb[*i1 - 1] == 25 && ee_1.lb[*i2 - 1] == -7))
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 7)
			{
				ii = *i1;
				if (x2 <= .5f)
				{
					ee_1.lb[*i1 - 1] = 2;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 5;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i1 - 1] = 1;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 4;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
			}
			else
			{
				ii = *i2;
				if (x2 <= .5f)
				{
					ee_1.lb[*i2 - 1] = 2;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 5;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i2 - 1] = 1;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 4;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
			}
		}
		if (ee_1.lb[*i1 - 1] == 7 && ee_1.lb[*i2 - 1] == 25 || ee_1.lb[*i1 - 1] == 25 && ee_1.lb[*i2 - 1] == 7 || (ee_1.lb[*i1 - 1] == -7 && ee_1.lb[*i2 - 1] == 27 || ee_1.lb[*i1 - 1] == 27 && ee_1.lb[*i2 - 1] == -7))
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 7)
			{
				ii = *i1;
				ee_1.lb[*i1 - 1] = 2;
				cc_1.e[*i1 - 1] = .939457f;
				ee_1.lb[*i2 - 1] = 3;
				cc_1.e[*i2 - 1] = .13496f;
				goto L40;
			}
			else
			{
				ii = *i2;
				ee_1.lb[*i2 - 1] = 2;
				cc_1.e[*i2 - 1] = .939457f;
				ee_1.lb[*i1 - 1] = 3;
				cc_1.e[*i1 - 1] = .13496f;
				goto L40;
			}
		}
		if (ee_1.lb[*i1 - 1] == 6 && ee_1.lb[*i2 - 1] == 27 || ee_1.lb[*i1 - 1] == 27 && ee_1.lb[*i2 - 1] == 6 || (ee_1.lb[*i1 - 1] == -6 && ee_1.lb[*i2 - 1] == 25 || ee_1.lb[*i1 - 1] == 25 && ee_1.lb[*i2 - 1] == -6))
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 6)
			{
				ii = *i1;
				if (x2 <= .5f)
				{
					ee_1.lb[*i1 - 1] = 2;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 4;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i1 - 1] = 1;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 3;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
			}
			else
			{
				ii = *i2;
				if (x2 <= .5f)
				{
					ee_1.lb[*i2 - 1] = 2;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 4;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i2 - 1] = 1;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 3;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
			}
		}
		if (ee_1.lb[*i1 - 1] == 9 && ee_1.lb[*i2 - 1] == 25 || ee_1.lb[*i1 - 1] == 25 && ee_1.lb[*i2 - 1] == 9 || (ee_1.lb[*i1 - 1] == -9 && ee_1.lb[*i2 - 1] == 27 || ee_1.lb[*i1 - 1] == 27 && ee_1.lb[*i2 - 1] == -9))
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 9)
			{
				ii = *i1;
				if (x2 <= .5f)
				{
					ee_1.lb[*i1 - 1] = 2;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 5;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i1 - 1] = 1;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 4;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
			}
			else
			{
				ii = *i2;
				if (x2 <= .5f)
				{
					ee_1.lb[*i2 - 1] = 2;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 5;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i2 - 1] = 1;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 4;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
			}
		}
		if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 9 && ee_1.lb[*i2 - 1] ==
															 26 ||
			ee_1.lb[*i1 - 1] == 26 && (i__2 = ee_1.lb[*i2 - 1], abs(
																	i__2)) == 9)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 9)
			{
				ii = *i1;
				ee_1.lb[*i1 - 1] = 1;
				cc_1.e[*i1 - 1] = .939457f;
				ee_1.lb[*i2 - 1] = 5;
				cc_1.e[*i2 - 1] = .13496f;
				goto L40;
			}
			else
			{
				ii = *i2;
				ee_1.lb[*i2 - 1] = 1;
				cc_1.e[*i2 - 1] = .939457f;
				ee_1.lb[*i1 - 1] = 5;
				cc_1.e[*i1 - 1] = .13496f;
				goto L40;
			}
		}
		if (ee_1.lb[*i1 - 1] == 11 && ee_1.lb[*i2 - 1] == 27 || ee_1.lb[*i1 - 1] == 27 && ee_1.lb[*i2 - 1] == 11 || ee_1.lb[*i1 - 1] == 13 && ee_1.lb[*i2 - 1] == 27 || ee_1.lb[*i1 - 1] == 27 && ee_1.lb[*i2 - 1] == 13 || (ee_1.lb[*i1 - 1] == -11 && ee_1.lb[*i2 - 1] == 25 || ee_1.lb[*i1 - 1] == 25 && ee_1.lb[*i2 - 1] == -11 || ee_1.lb[*i1 - 1] == -13 && ee_1.lb[*i2 - 1] == 25 || ee_1.lb[*i1 - 1] == 25 && ee_1.lb[*i2 - 1] == -13))
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 11 || (i__2 = ee_1.lb[*i1 - 1], abs(i__2)) == 13)
			{
				ii = *i1;
				ee_1.lb[*i1 - 1] = 1;
				cc_1.e[*i1 - 1] = .939457f;
				ee_1.lb[*i2 - 1] = 5;
				cc_1.e[*i2 - 1] = .13496f;
				goto L40;
			}
			else
			{
				ii = *i2;
				ee_1.lb[*i2 - 1] = 1;
				cc_1.e[*i2 - 1] = .939457f;
				ee_1.lb[*i1 - 1] = 5;
				cc_1.e[*i1 - 1] = .13496f;
				goto L40;
			}
		}
		if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 10 && ee_1.lb[*i2 - 1] ==
															  26 ||
			ee_1.lb[*i1 - 1] == 26 && (i__2 = ee_1.lb[*i2 - 1], abs(
																	i__2)) == 10 ||
			ee_1.lb[*i1 - 1] == 26 && (i__3 = ee_1.lb[*i2 - 1], abs(i__3)) == 12 || ee_1.lb[*i2 - 1] == 26 && (i__4 = ee_1.lb[*i1 - 1], abs(i__4)) == 12)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 10 || (i__2 = ee_1.lb[*i1 - 1], abs(i__2)) == 12)
			{
				ii = *i1;
				if (x2 <= .5f)
				{
					ee_1.lb[*i1 - 1] = 2;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 4;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i1 - 1] = 1;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 3;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
			}
			else
			{
				ii = *i2;
				if (x2 <= .5f)
				{
					ee_1.lb[*i2 - 1] = 2;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 4;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i2 - 1] = 1;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 3;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
			}
		}
		if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 11 && ee_1.lb[*i2 - 1] ==
															  26 ||
			ee_1.lb[*i1 - 1] == 26 && (i__2 = ee_1.lb[*i2 - 1], abs(
																	i__2)) == 11 ||
			ee_1.lb[*i1 - 1] == 26 && (i__3 = ee_1.lb[*i2 - 1], abs(i__3)) == 13 || ee_1.lb[*i2 - 1] == 26 && (i__4 = ee_1.lb[*i1 - 1], abs(i__4)) == 13)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 11 || (i__2 = ee_1.lb[*i1 - 1], abs(i__2)) == 13)
			{
				ii = *i1;
				if (x2 <= .5f)
				{
					ee_1.lb[*i1 - 1] = 2;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 5;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i1 - 1] = 1;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 4;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
			}
			else
			{
				ii = *i2;
				if (x2 <= .5f)
				{
					ee_1.lb[*i2 - 1] = 2;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 5;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i2 - 1] = 1;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 4;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
			}
		}
		if (ee_1.lb[*i1 - 1] == 11 && ee_1.lb[*i2 - 1] == 25 || ee_1.lb[*i1 - 1] == 25 && ee_1.lb[*i2 - 1] == 11 || ee_1.lb[*i1 - 1] == 25 && ee_1.lb[*i2 - 1] == 13 || ee_1.lb[*i2 - 1] == 25 && ee_1.lb[*i1 - 1] == 13 || (ee_1.lb[*i1 - 1] == -11 && ee_1.lb[*i2 - 1] == 27 || ee_1.lb[*i1 - 1] == 27 && ee_1.lb[*i2 - 1] == -11 || ee_1.lb[*i1 - 1] == 27 && ee_1.lb[*i2 - 1] == -13 || ee_1.lb[*i2 - 1] == 27 && ee_1.lb[*i1 - 1] == -13))
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 11 || (i__2 = ee_1.lb[*i1 - 1], abs(i__2)) == 13)
			{
				ii = *i1;
				if (x2 <= .5f)
				{
					ee_1.lb[*i1 - 1] = 2;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 4;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i1 - 1] = 1;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 3;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
			}
			else
			{
				ii = *i2;
				if (x2 <= .5f)
				{
					ee_1.lb[*i2 - 1] = 2;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 4;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i2 - 1] = 1;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 3;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
			}
		}
		if (ee_1.lb[*i1 - 1] == 10 && ee_1.lb[*i2 - 1] == 27 || ee_1.lb[*i1 - 1] == 27 && ee_1.lb[*i2 - 1] == 10 || ee_1.lb[*i1 - 1] == 12 && ee_1.lb[*i2 - 1] == 27 || ee_1.lb[*i1 - 1] == 27 && ee_1.lb[*i2 - 1] == 12 || (ee_1.lb[*i1 - 1] == -10 && ee_1.lb[*i2 - 1] == 25 || ee_1.lb[*i1 - 1] == 25 && ee_1.lb[*i2 - 1] == -10 || ee_1.lb[*i1 - 1] == -12 && ee_1.lb[*i2 - 1] == 25 || ee_1.lb[*i1 - 1] == 25 && ee_1.lb[*i2 - 1] == -12))
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 10 || (i__2 = ee_1.lb[*i1 - 1], abs(i__2)) == 12)
			{
				ii = *i1;
				if (x2 <= .5f)
				{
					ee_1.lb[*i1 - 1] = 2;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 5;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i1 - 1] = 1;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 4;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
			}
			else
			{
				ii = *i2;
				if (x2 <= .5f)
				{
					ee_1.lb[*i2 - 1] = 2;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 5;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i2 - 1] = 1;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 4;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
			}
		}
		if (ee_1.lb[*i1 - 1] == 10 && ee_1.lb[*i2 - 1] == 25 || ee_1.lb[*i1 - 1] == 25 && ee_1.lb[*i2 - 1] == 10 || ee_1.lb[*i1 - 1] == 25 && ee_1.lb[*i2 - 1] == 12 || ee_1.lb[*i1 - 1] == 12 && ee_1.lb[*i2 - 1] == 25 || (ee_1.lb[*i1 - 1] == -10 && ee_1.lb[*i2 - 1] == 27 || ee_1.lb[*i1 - 1] == 27 && ee_1.lb[*i2 - 1] == -10 || ee_1.lb[*i1 - 1] == 27 && ee_1.lb[*i2 - 1] == -12 || ee_1.lb[*i1 - 1] == -12 && ee_1.lb[*i2 - 1] == 27))
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 10 || (i__2 = ee_1.lb[*i1 - 1], abs(i__2)) == 12)
			{
				ii = *i1;
				ee_1.lb[*i1 - 1] = 2;
				cc_1.e[*i1 - 1] = .939457f;
				ee_1.lb[*i2 - 1] = 3;
				cc_1.e[*i2 - 1] = .13496f;
				goto L40;
			}
			else
			{
				ii = *i2;
				ee_1.lb[*i2 - 1] = 2;
				cc_1.e[*i2 - 1] = .939457f;
				ee_1.lb[*i1 - 1] = 3;
				cc_1.e[*i1 - 1] = .13496f;
				goto L40;
			}
		}
	L60:
		*iblock = 82;
		if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 7 && ee_1.lb[*i2 - 1] ==
															 28 ||
			ee_1.lb[*i1 - 1] == 28 && (i__2 = ee_1.lb[*i2 - 1], abs(
																	i__2)) == 7)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 7)
			{
				ii = *i1;
				if (x2 <= .5f)
				{
					ee_1.lb[*i1 - 1] = 2;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 4;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i1 - 1] = 1;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 3;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
			}
			else
			{
				ii = *i2;
				if (x2 <= .5f)
				{
					ee_1.lb[*i2 - 1] = 2;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 4;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i2 - 1] = 1;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 3;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
			}
		}
		if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 8 && ee_1.lb[*i2 - 1] ==
															 28 ||
			ee_1.lb[*i1 - 1] == 28 && (i__2 = ee_1.lb[*i2 - 1], abs(
																	i__2)) == 8)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 8)
			{
				ii = *i1;
				if (x2 <= .5f)
				{
					ee_1.lb[*i1 - 1] = 2;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 5;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i1 - 1] = 1;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 4;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
			}
			else
			{
				ii = *i2;
				if (x2 <= .5f)
				{
					ee_1.lb[*i2 - 1] = 2;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 5;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i2 - 1] = 1;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 4;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
			}
		}
		if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 6 && ee_1.lb[*i2 - 1] ==
															 28 ||
			ee_1.lb[*i1 - 1] == 28 && (i__2 = ee_1.lb[*i2 - 1], abs(
																	i__2)) == 6)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 6)
			{
				ii = *i1;
				ee_1.lb[*i1 - 1] = 2;
				cc_1.e[*i1 - 1] = .939457f;
				ee_1.lb[*i2 - 1] = 3;
				cc_1.e[*i2 - 1] = .13496f;
				goto L40;
			}
			else
			{
				ii = *i2;
				ee_1.lb[*i2 - 1] = 2;
				cc_1.e[*i2 - 1] = .939457f;
				ee_1.lb[*i1 - 1] = 3;
				cc_1.e[*i1 - 1] = .13496f;
				goto L40;
			}
		}
		if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 9 && ee_1.lb[*i2 - 1] ==
															 28 ||
			ee_1.lb[*i1 - 1] == 28 && (i__2 = ee_1.lb[*i2 - 1], abs(
																	i__2)) == 9)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 9)
			{
				ii = *i1;
				ee_1.lb[*i1 - 1] = 1;
				cc_1.e[*i1 - 1] = .939457f;
				ee_1.lb[*i2 - 1] = 5;
				cc_1.e[*i2 - 1] = .13496f;
				goto L40;
			}
			else
			{
				ii = *i2;
				ee_1.lb[*i2 - 1] = 1;
				cc_1.e[*i2 - 1] = .939457f;
				ee_1.lb[*i1 - 1] = 5;
				cc_1.e[*i1 - 1] = .13496f;
				goto L40;
			}
		}
		if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 10 && ee_1.lb[*i2 - 1] ==
															  28 ||
			ee_1.lb[*i1 - 1] == 28 && (i__2 = ee_1.lb[*i2 - 1], abs(
																	i__2)) == 10 ||
			ee_1.lb[*i1 - 1] == 28 && (i__3 = ee_1.lb[*i2 - 1], abs(i__3)) == 12 || ee_1.lb[*i2 - 1] == 28 && (i__4 = ee_1.lb[*i1 - 1], abs(i__4)) == 12)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 10 || (i__2 = ee_1.lb[*i1 - 1], abs(i__2)) == 12)
			{
				ii = *i1;
				if (x2 <= .5f)
				{
					ee_1.lb[*i1 - 1] = 2;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 4;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i1 - 1] = 1;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 3;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
			}
			else
			{
				ii = *i2;
				if (x2 <= .5f)
				{
					ee_1.lb[*i2 - 1] = 2;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 4;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i2 - 1] = 1;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 3;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
			}
		}
		if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 11 && ee_1.lb[*i2 - 1] ==
															  28 ||
			ee_1.lb[*i1 - 1] == 28 && (i__2 = ee_1.lb[*i2 - 1], abs(
																	i__2)) == 11 ||
			ee_1.lb[*i1 - 1] == 28 && (i__3 = ee_1.lb[*i2 - 1], abs(i__3)) == 13 || ee_1.lb[*i2 - 1] == 28 && (i__4 = ee_1.lb[*i1 - 1], abs(i__4)) == 13)
		{
			if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 11 || (i__2 = ee_1.lb[*i1 - 1], abs(i__2)) == 13)
			{
				ii = *i1;
				if (x2 <= .5f)
				{
					ee_1.lb[*i1 - 1] = 2;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 5;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i1 - 1] = 1;
					cc_1.e[*i1 - 1] = .939457f;
					ee_1.lb[*i2 - 1] = 4;
					cc_1.e[*i2 - 1] = .13496f;
					goto L40;
				}
			}
			else
			{
				ii = *i2;
				if (x2 <= .5f)
				{
					ee_1.lb[*i2 - 1] = 2;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 5;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
				else
				{
					ee_1.lb[*i2 - 1] = 1;
					cc_1.e[*i2 - 1] = .939457f;
					ee_1.lb[*i1 - 1] = 4;
					cc_1.e[*i1 - 1] = .13496f;
					goto L40;
				}
			}
		}
	L40:
		em1 = cc_1.e[*i1 - 1];
		em2 = cc_1.e[*i2 - 1];
		if (ianti == 1 && ee_1.lb[*i1 - 1] >= 1 && ee_1.lb[*i2 - 1] >= 1)
		{
			ee_1.lb[ii - 1] = -ee_1.lb[ii - 1];
			jj = *i2;
			if (ii == *i2)
			{
				jj = *i1;
			}
			if (ee_1.lb[jj - 1] == 3)
			{
				ee_1.lb[jj - 1] = 5;
			}
			else if (ee_1.lb[jj - 1] == 5)
			{
				ee_1.lb[jj - 1] = 3;
			}
		}
	}
L50:
	r__2 = *srt;
	r__3 = em1;
	r__4 = em2;
	r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
	r__5 = em1 * em2;
	pr2 = r__1 * r__1 - r__5 * r__5 * 4.f;
	if (pr2 <= 0.f)
	{
		pr2 = 1e-9f;
	}
	pr = sqrt(pr2) / (*srt * 2.f);
	xptr = pr * .33f;
	cc1 = ptr_(&xptr, &input1_1.iseed);
	r__1 = pr;
	r__2 = cc1;
	c1 = sqrt(r__1 * r__1 - r__2 * r__2) / pr;
	t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
	r__1 = c1;
	s1 = sqrt(1.f - r__1 * r__1);
	ct1 = cos(t1);
	st1 = sin(t1);
	*pz = pr * c1;
	*px = pr * s1 * ct1;
	*py = pr * s1 * st1;
	rotate_(&px0, &py0, &pz0, px, py, pz);
	return 0;
}

int crlaba_(float *px, float *py, float *pz, float *srt, float *brel, float *brsgm, int *i1, int *i2, int *nt, int *iblock, int *nchrg, int *icase)
{
	float r__1, r__2, r__3, r__4, r__5;

	double sqrt(double), cos(double), sin(double);

	static float c1, s1, t1, pr, em1, em2, ct1, pr2, px0, py0, pz0, st1, rrr;
	extern double ranart_(int *);
	extern int rotate_(float *, float *, float *, float *, float *, float *);

	px0 = *px;
	py0 = *py;
	pz0 = *pz;

	if (*icase == 3)
	{
		rrr = ranart_(&rndf77_1.nseed);
		if (rrr < *brel)
		{
			*iblock = 8;
		}
		else
		{
			*iblock = 100;
			if (rrr < *brel + *brsgm)
			{
				ee_1.lb[*i1 - 1] = -15 - (int)(ranart_(&rndf77_1.nseed) *
												   3);
				cc_1.e[*i1 - 1] = 1.1974f;
			}
			else
			{
				ee_1.lb[*i1 - 1] = -14;
				cc_1.e[*i1 - 1] = 1.1157f;
			}
			ee_1.lb[*i2 - 1] = (int)(ranart_(&rndf77_1.nseed) * 3) + 3;
			cc_1.e[*i2 - 1] = .138f;
		}
	}

	if (*icase == 4)
	{
		rrr = ranart_(&rndf77_1.nseed);
		if (rrr < *brel)
		{
			*iblock = 8;
		}
		else
		{
			*iblock = 102;
			ee_1.lb[*i1 - 1] = 23;
			ee_1.lb[*i2 - 1] = -1 - (int)(ranart_(&rndf77_1.nseed) * 2);
			if (*nchrg == -2)
			{
				ee_1.lb[*i2 - 1] = -6;
			}
			if (*nchrg == 1)
			{
				ee_1.lb[*i2 - 1] = -9;
			}
			cc_1.e[*i1 - 1] = .498f;
			cc_1.e[*i2 - 1] = .938f;
			if (*nchrg == -2 || *nchrg == 1)
			{
				cc_1.e[*i2 - 1] = 1.232f;
			}
		}
	}

	em1 = cc_1.e[*i1 - 1];
	em2 = cc_1.e[*i2 - 1];
	r__2 = *srt;
	r__3 = em1;
	r__4 = em2;
	r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
	r__5 = em1 * em2;
	pr2 = r__1 * r__1 - r__5 * r__5 * 4.f;
	if (pr2 <= 0.f)
	{
		pr2 = 1e-9f;
	}
	pr = sqrt(pr2) / (*srt * 2.f);
	c1 = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
	t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
	r__1 = c1;
	s1 = sqrt(1.f - r__1 * r__1);
	ct1 = cos(t1);
	st1 = sin(t1);
	*pz = pr * c1;
	*px = pr * s1 * ct1;
	*py = pr * s1 * st1;
	rotate_(&px0, &py0, &pz0, px, py, pz);
	return 0;
}

int crkn_(float *px, float *py, float *pz, float *srt, int *i1, int *i2, int *iblock)
{
	float r__1, r__2, r__3, r__4, r__5;

	double sqrt(double), cos(double), sin(double);

	static float c1, s1, t1, pr, em1, em2, ct1, pr2, px0, py0, pz0, st1;
	static int ntag;
	extern double ranart_(int *);

	px0 = *px;
	py0 = *py;
	pz0 = *pz;
	*iblock = 8;
	ntag = 0;
	em1 = cc_1.e[*i1 - 1];
	em2 = cc_1.e[*i2 - 1];
	r__2 = *srt;
	r__3 = em1;
	r__4 = em2;
	r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
	r__5 = em1 * em2;
	pr2 = r__1 * r__1 - r__5 * r__5 * 4.f;
	if (pr2 <= 0.f)
	{
		pr2 = 1e-9f;
	}
	pr = sqrt(pr2) / (*srt * 2.f);
	c1 = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
	t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
	r__1 = c1;
	s1 = sqrt(1.f - r__1 * r__1);
	ct1 = cos(t1);
	st1 = sin(t1);
	*pz = pr * c1;
	*px = pr * s1 * ct1;
	*py = pr * s1 * st1;
	return 0;
}

int crppba_(float *px, float *py, float *pz, float *srt, int *i1, int *i2, int *iblock)
{
	float r__1, r__2, r__3, r__4, r__5;

	double sqrt(double), cos(double), sin(double);

	static float c1, s1, t1, pr, em1, em2, ct1, pr2, px0, py0, pz0, st1;
	static int npion, nchrg1, nchrg2;
	static float pmass1, pmass2;
	extern int pbarfs_(float *, int *, int *);
	extern double ranart_(int *);
	extern int rotate_(float *, float *, float *, float *, float *, float *);

	px0 = *px;
	py0 = *py;
	pz0 = *pz;
	pbarfs_(srt, &npion, &input1_1.iseed);
	nchrg1 = (int)(ranart_(&rndf77_1.nseed) * 3) + 3;
	nchrg2 = (int)(ranart_(&rndf77_1.nseed) * 3) + 3;
	pmass1 = .13496f;
	pmass2 = .13496f;
	if (nchrg1 == 3 || nchrg1 == 5)
	{
		pmass1 = .13957f;
	}
	if (nchrg2 == 3 || nchrg2 == 5)
	{
		pmass2 = .13957f;
	}
	if (npion == 2)
	{
		*iblock = 1902;
		ee_1.lb[*i1 - 1] = nchrg1;
		cc_1.e[*i1 - 1] = pmass1;
		ee_1.lb[*i2 - 1] = nchrg2;
		cc_1.e[*i2 - 1] = pmass2;
		goto L50;
	}
	if (npion == 3)
	{
		*iblock = 1903;
		ee_1.lb[*i1 - 1] = nchrg1;
		cc_1.e[*i1 - 1] = pmass1;
		ee_1.lb[*i2 - 1] = nchrg2 + 22;
		cc_1.e[*i2 - 1] = .769f;
		goto L50;
	}
	if (npion == 4)
	{
		*iblock = 1904;
		if (ranart_(&rndf77_1.nseed) >= .5f)
		{
			ee_1.lb[*i1 - 1] = nchrg1 + 22;
			cc_1.e[*i1 - 1] = .769f;
			ee_1.lb[*i2 - 1] = nchrg2 + 22;
			cc_1.e[*i2 - 1] = .769f;
		}
		else
		{
			ee_1.lb[*i1 - 1] = nchrg1;
			cc_1.e[*i1 - 1] = pmass1;
			ee_1.lb[*i2 - 1] = 28;
			cc_1.e[*i2 - 1] = .782f;
		}
		goto L50;
	}
	if (npion == 5)
	{
		*iblock = 1905;
		ee_1.lb[*i1 - 1] = nchrg1 + 22;
		cc_1.e[*i1 - 1] = .769f;
		ee_1.lb[*i2 - 1] = 28;
		cc_1.e[*i2 - 1] = .782f;
		goto L50;
	}
	if (npion == 6)
	{
		*iblock = 1906;
		ee_1.lb[*i1 - 1] = 28;
		cc_1.e[*i1 - 1] = .782f;
		ee_1.lb[*i2 - 1] = 28;
		cc_1.e[*i2 - 1] = .782f;
	}
L50:
	em1 = cc_1.e[*i1 - 1];
	em2 = cc_1.e[*i2 - 1];
	r__2 = *srt;
	r__3 = em1;
	r__4 = em2;
	r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
	r__5 = em1 * em2;
	pr2 = r__1 * r__1 - r__5 * r__5 * 4.f;
	if (pr2 <= 0.f)
	{
		pr2 = 1e-8f;
	}
	pr = sqrt(pr2) / (*srt * 2.f);
	c1 = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
	t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
	r__1 = c1;
	s1 = sqrt(1.f - r__1 * r__1);
	ct1 = cos(t1);
	st1 = sin(t1);
	*pz = pr * c1;
	*px = pr * s1 * ct1;
	*py = pr * s1 * st1;
	rotate_(&px0, &py0, &pz0, px, py, pz);
	return 0;
}

int crkkpi_(int *i1, int *i2, float *xsk1, float *xsk2,
			float *xsk3, float *xsk4, float *xsk5, float *xsk6, float *xsk7, float *xsk8, float *xsk9, float *xsk10, float *xsk11, float *sigk, int *iblock, int *lbp1, int *lbp2, float *emm1, float *emm2)
{
	static float x1;
	extern double ranart_(int *);
	extern int rhores_(int *, int *);

	*iblock = 1907;
	x1 = ranart_(&rndf77_1.nseed) * *sigk;
	*xsk2 = *xsk1 + *xsk2;
	*xsk3 = *xsk2 + *xsk3;
	*xsk4 = *xsk3 + *xsk4;
	*xsk5 = *xsk4 + *xsk5;
	*xsk6 = *xsk5 + *xsk6;
	*xsk7 = *xsk6 + *xsk7;
	*xsk8 = *xsk7 + *xsk8;
	*xsk9 = *xsk8 + *xsk9;
	*xsk10 = *xsk9 + *xsk10;
	if (x1 <= *xsk1)
	{
		ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 3) + 3;
		ee_1.lb[*i2 - 1] = (int)(ranart_(&rndf77_1.nseed) * 3) + 3;
		cc_1.e[*i1 - 1] = .13957f;
		cc_1.e[*i2 - 1] = .13957f;
		goto L100;
	}
	else if (x1 <= *xsk2)
	{
		ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 3) + 3;
		ee_1.lb[*i2 - 1] = (int)(ranart_(&rndf77_1.nseed) * 3) + 25;
		cc_1.e[*i1 - 1] = .13957f;
		cc_1.e[*i2 - 1] = .769f;
		goto L100;
	}
	else if (x1 <= *xsk3)
	{
		ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 3) + 3;
		ee_1.lb[*i2 - 1] = 28;
		cc_1.e[*i1 - 1] = .13957f;
		cc_1.e[*i2 - 1] = .782f;
		goto L100;
	}
	else if (x1 <= *xsk4)
	{
		ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 3) + 3;
		ee_1.lb[*i2 - 1] = 0;
		cc_1.e[*i1 - 1] = .13957f;
		cc_1.e[*i2 - 1] = .5473f;
		goto L100;
	}
	else if (x1 <= *xsk5)
	{
		ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 3) + 25;
		ee_1.lb[*i2 - 1] = (int)(ranart_(&rndf77_1.nseed) * 3) + 25;
		cc_1.e[*i1 - 1] = .769f;
		cc_1.e[*i2 - 1] = .769f;
		goto L100;
	}
	else if (x1 <= *xsk6)
	{
		ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 3) + 25;
		ee_1.lb[*i2 - 1] = 28;
		cc_1.e[*i1 - 1] = .769f;
		cc_1.e[*i2 - 1] = .782f;
		goto L100;
	}
	else if (x1 <= *xsk7)
	{
		ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 3) + 25;
		ee_1.lb[*i2 - 1] = 0;
		cc_1.e[*i1 - 1] = .769f;
		cc_1.e[*i2 - 1] = .5473f;
		goto L100;
	}
	else if (x1 <= *xsk8)
	{
		ee_1.lb[*i1 - 1] = 28;
		ee_1.lb[*i2 - 1] = 28;
		cc_1.e[*i1 - 1] = .782f;
		cc_1.e[*i2 - 1] = .782f;
		goto L100;
	}
	else if (x1 <= *xsk9)
	{
		ee_1.lb[*i1 - 1] = 28;
		ee_1.lb[*i2 - 1] = 0;
		cc_1.e[*i1 - 1] = .782f;
		cc_1.e[*i2 - 1] = .5473f;
		goto L100;
	}
	else if (x1 <= *xsk10)
	{
		ee_1.lb[*i1 - 1] = 0;
		ee_1.lb[*i2 - 1] = 0;
		cc_1.e[*i1 - 1] = .5473f;
		cc_1.e[*i2 - 1] = .5473f;
	}
	else
	{
		*iblock = 222;
		rhores_(i1, i2);
		ee_1.lb[*i1 - 1] = 29;
		cc_1.e[*i2 - 1] = 0.f;
	}
L100:
	*lbp1 = ee_1.lb[*i1 - 1];
	*lbp2 = ee_1.lb[*i2 - 1];
	*emm1 = cc_1.e[*i1 - 1];
	*emm2 = cc_1.e[*i2 - 1];
	return 0;
}

int crkhyp_(float *px, float *py, float *pz, float *srt, int *i1, int *i2, float *xky1, float *xky2, float *xky3, float *xky4,
			float *xky5, float *xky6, float *xky7, float *xky8, float *xky9, float *xky10, float *xky11, float *xky12, float *xky13, float *xky14, float *xky15, float *xky16, float *xky17, float *sigk, int *ikmp, int *iblock)
{
	float r__1, r__2, r__3, r__4, r__5;

	double sqrt(double), cos(double), sin(double);

	static float c1, s1, t1, x1, pr, em1, em2, ct1, pr2, px0, py0, pz0, st1;
	extern double ranart_(int *);
	extern int rotate_(float *, float *, float *, float *, float *, float *);

	px0 = *px;
	py0 = *py;
	pz0 = *pz;
	*iblock = 1908;

	x1 = ranart_(&rndf77_1.nseed) * *sigk;
	*xky2 = *xky1 + *xky2;
	*xky3 = *xky2 + *xky3;
	*xky4 = *xky3 + *xky4;
	*xky5 = *xky4 + *xky5;
	*xky6 = *xky5 + *xky6;
	*xky7 = *xky6 + *xky7;
	*xky8 = *xky7 + *xky8;
	*xky9 = *xky8 + *xky9;
	*xky10 = *xky9 + *xky10;
	*xky11 = *xky10 + *xky11;
	*xky12 = *xky11 + *xky12;
	*xky13 = *xky12 + *xky13;
	*xky14 = *xky13 + *xky14;
	*xky15 = *xky14 + *xky15;
	*xky16 = *xky15 + *xky16;
	if (x1 <= *xky1)
	{
		ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 3) + 3;
		ee_1.lb[*i2 - 1] = (int)(ranart_(&rndf77_1.nseed) * 2) + 1;
		cc_1.e[*i1 - 1] = .14f;
		cc_1.e[*i2 - 1] = .93828f;
		goto L100;
	}
	else if (x1 <= *xky2)
	{
		ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 3) + 3;
		ee_1.lb[*i2 - 1] = (int)(ranart_(&rndf77_1.nseed) * 4) + 6;
		cc_1.e[*i1 - 1] = .14f;
		cc_1.e[*i2 - 1] = 1.232f;
		goto L100;
	}
	else if (x1 <= *xky3)
	{
		ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 3) + 3;
		ee_1.lb[*i2 - 1] = (int)(ranart_(&rndf77_1.nseed) * 2) + 10;
		cc_1.e[*i1 - 1] = .14f;
		cc_1.e[*i2 - 1] = 1.44f;
		goto L100;
	}
	else if (x1 <= *xky4)
	{
		ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 3) + 3;
		ee_1.lb[*i2 - 1] = (int)(ranart_(&rndf77_1.nseed) * 2) + 12;
		cc_1.e[*i1 - 1] = .14f;
		cc_1.e[*i2 - 1] = 1.535f;
		goto L100;
	}
	else if (x1 <= *xky5)
	{
		ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 3) + 25;
		ee_1.lb[*i2 - 1] = (int)(ranart_(&rndf77_1.nseed) * 2) + 1;
		cc_1.e[*i1 - 1] = .769f;
		cc_1.e[*i2 - 1] = .93828f;
		goto L100;
	}
	else if (x1 <= *xky6)
	{
		ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 3) + 25;
		ee_1.lb[*i2 - 1] = (int)(ranart_(&rndf77_1.nseed) * 4) + 6;
		cc_1.e[*i1 - 1] = .769f;
		cc_1.e[*i2 - 1] = 1.232f;
		goto L100;
	}
	else if (x1 <= *xky7)
	{
		ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 3) + 25;
		ee_1.lb[*i2 - 1] = (int)(ranart_(&rndf77_1.nseed) * 2) + 10;
		cc_1.e[*i1 - 1] = .769f;
		cc_1.e[*i2 - 1] = 1.44f;
		goto L100;
	}
	else if (x1 <= *xky8)
	{
		ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 3) + 25;
		ee_1.lb[*i2 - 1] = (int)(ranart_(&rndf77_1.nseed) * 2) + 12;
		cc_1.e[*i1 - 1] = .769f;
		cc_1.e[*i2 - 1] = 1.535f;
		goto L100;
	}
	else if (x1 <= *xky9)
	{
		ee_1.lb[*i1 - 1] = 28;
		ee_1.lb[*i2 - 1] = (int)(ranart_(&rndf77_1.nseed) * 2) + 1;
		cc_1.e[*i1 - 1] = .782f;
		cc_1.e[*i2 - 1] = .93828f;
		goto L100;
	}
	else if (x1 <= *xky10)
	{
		ee_1.lb[*i1 - 1] = 28;
		ee_1.lb[*i2 - 1] = (int)(ranart_(&rndf77_1.nseed) * 4) + 6;
		cc_1.e[*i1 - 1] = .782f;
		cc_1.e[*i2 - 1] = 1.232f;
		goto L100;
	}
	else if (x1 <= *xky11)
	{
		ee_1.lb[*i1 - 1] = 28;
		ee_1.lb[*i2 - 1] = (int)(ranart_(&rndf77_1.nseed) * 2) + 10;
		cc_1.e[*i1 - 1] = .782f;
		cc_1.e[*i2 - 1] = 1.44f;
		goto L100;
	}
	else if (x1 <= *xky12)
	{
		ee_1.lb[*i1 - 1] = 28;
		ee_1.lb[*i2 - 1] = (int)(ranart_(&rndf77_1.nseed) * 2) + 12;
		cc_1.e[*i1 - 1] = .782f;
		cc_1.e[*i2 - 1] = 1.535f;
		goto L100;
	}
	else if (x1 <= *xky13)
	{
		ee_1.lb[*i1 - 1] = 0;
		ee_1.lb[*i2 - 1] = (int)(ranart_(&rndf77_1.nseed) * 2) + 1;
		cc_1.e[*i1 - 1] = .5473f;
		cc_1.e[*i2 - 1] = .93828f;
		goto L100;
	}
	else if (x1 <= *xky14)
	{
		ee_1.lb[*i1 - 1] = 0;
		ee_1.lb[*i2 - 1] = (int)(ranart_(&rndf77_1.nseed) * 4) + 6;
		cc_1.e[*i1 - 1] = .5473f;
		cc_1.e[*i2 - 1] = 1.232f;
		goto L100;
	}
	else if (x1 <= *xky15)
	{
		ee_1.lb[*i1 - 1] = 0;
		ee_1.lb[*i2 - 1] = (int)(ranart_(&rndf77_1.nseed) * 2) + 10;
		cc_1.e[*i1 - 1] = .5473f;
		cc_1.e[*i2 - 1] = 1.44f;
		goto L100;
	}
	else if (x1 <= *xky16)
	{
		ee_1.lb[*i1 - 1] = 0;
		ee_1.lb[*i2 - 1] = (int)(ranart_(&rndf77_1.nseed) * 2) + 12;
		cc_1.e[*i1 - 1] = .5473f;
		cc_1.e[*i2 - 1] = 1.535f;
		goto L100;
	}
	else
	{
		ee_1.lb[*i1 - 1] = 29;
		ee_1.lb[*i2 - 1] = (int)(ranart_(&rndf77_1.nseed) * 2) + 1;
		cc_1.e[*i1 - 1] = 1.02f;
		cc_1.e[*i2 - 1] = .939457f;
		*iblock = 222;
		goto L100;
	}
L100:
	if (*ikmp == -1)
	{
		ee_1.lb[*i2 - 1] = -ee_1.lb[*i2 - 1];
	}
	em1 = cc_1.e[*i1 - 1];
	em2 = cc_1.e[*i2 - 1];
	r__2 = *srt;
	r__3 = em1;
	r__4 = em2;
	r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
	r__5 = em1 * em2;
	pr2 = r__1 * r__1 - r__5 * r__5 * 4.f;
	if (pr2 <= 0.f)
	{
		pr2 = 1e-8f;
	}
	pr = sqrt(pr2) / (*srt * 2.f);
	c1 = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
	t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
	r__1 = c1;
	s1 = sqrt(1.f - r__1 * r__1);
	ct1 = cos(t1);
	st1 = sin(t1);
	*pz = pr * c1;
	*px = pr * s1 * ct1;
	*py = pr * s1 * st1;
	rotate_(&px0, &py0, &pz0, px, py, pz);
	return 0;
}

int crlan_(float *px, float *py, float *pz, float *srt, int *i1, int *i2, int *iblock)
{
	float r__1, r__2, r__3, r__4, r__5;

	double sqrt(double), cos(double), sin(double);

	static float c1, s1, t1, pr, em1, em2, ct1, pr2, px0, py0, pz0, st1;
	static int ntag;
	extern double ranart_(int *);

	px0 = *px;
	py0 = *py;
	pz0 = *pz;
	*iblock = 71;
	ntag = 0;
	if (ee_1.lb[*i1 - 1] >= 14 && ee_1.lb[*i1 - 1] <= 17 || ee_1.lb[*i2 - 1] >= 14 && ee_1.lb[*i2 - 1] <= 17)
	{
		ee_1.lb[*i1 - 1] = 21;
	}
	else
	{
		ee_1.lb[*i1 - 1] = 23;
	}
	ee_1.lb[*i2 - 1] = (int)(ranart_(&rndf77_1.nseed) * 3) + 3;
	cc_1.e[*i1 - 1] = .498f;
	cc_1.e[*i2 - 1] = .138f;
	em1 = cc_1.e[*i1 - 1];
	em2 = cc_1.e[*i2 - 1];
	r__2 = *srt;
	r__3 = em1;
	r__4 = em2;
	r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
	r__5 = em1 * em2;
	pr2 = r__1 * r__1 - r__5 * r__5 * 4.f;
	if (pr2 <= 0.f)
	{
		pr2 = 1e-9f;
	}
	pr = sqrt(pr2) / (*srt * 2.f);
	c1 = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
	t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
	r__1 = c1;
	s1 = sqrt(1.f - r__1 * r__1);
	ct1 = cos(t1);
	st1 = sin(t1);
	*pz = pr * c1;
	*px = pr * s1 * ct1;
	*py = pr * s1 * st1;
	return 0;
}

int crkpla_(float *px, float *py, float *pz, float *ec, float *srt, float *spika, float *emm1, float *emm2, int *lbp1, int *lbp2, int *i1, int *i2, int *icase, float *srhoks)
{
	float r__1, r__2, r__3, r__4, r__5, r__6;

	double sqrt(double), cos(double), sin(double);

	static float c1, s1, t1;
	static int ic;
	static float pr, ct1, pr2, px0, py0, pz0, st1, pdd, pff, sig1, sig2, xkp0,
		xkp1, xkp2, xkp3, xkp4, xkp5, xkp6, xkp7, xkp8, xkp9, sigm, dskn,
		xkp10, randu, sigkp, dsknr;
	extern int distce_(int *, int *, float *, float *,
					   float *, float *, float *, int *, float *, float *, float *);
	static float sigpik;
	extern double ranart_(int *);

	*emm1 = 0.f;
	*emm2 = 0.f;
	*lbp1 = 0;
	*lbp2 = 0;
	xkp0 = *spika;
	xkp1 = 0.f;
	xkp2 = 0.f;
	xkp3 = 0.f;
	xkp4 = 0.f;
	xkp5 = 0.f;
	xkp6 = 0.f;
	xkp7 = 0.f;
	xkp8 = 0.f;
	xkp9 = 0.f;
	xkp10 = 0.f;
	sigm = 15.f;
	r__1 = *srt;
	r__2 = *srt;
	pdd = (r__1 * r__1 - .40063836159999994f) * (r__2 * r__2 -
												 .13179804160000003f);

	if (*srt < 2.0551569999999999f)
	{
		goto L70;
	}
	r__1 = *srt;
	r__2 = *srt;
	xkp1 = sigm * 1.3333333333333333f * (r__1 * r__1 - 4.2236702946489997f) *
		   (r__2 * r__2 - .031061595048999975f) / pdd;
	if (*srt > 2.3476999999999997f)
	{
		r__1 = *srt;
		r__2 = *srt;
		xkp2 = sigm * 5.333333333333333f * (r__1 * r__1 - 5.5116952899999987f) * (r__2 * r__2 - .013525690000000016f) / pdd;
	}
	if (*srt > 2.5556999999999999f)
	{
		r__1 = *srt;
		r__2 = *srt;
		xkp3 = sigm * 1.3333333333333333f * (r__1 * r__1 - 6.5316024899999992f) * (r__2 * r__2 - .10517049000000002f) /
			   pdd;
	}
	if (*srt > 2.6506999999999996f)
	{
		r__1 = *srt;
		r__2 = *srt;
		xkp4 = sigm * 1.3333333333333333f * (r__1 * r__1 - 7.0262104899999978f) * (r__2 * r__2 - .17581249000000002f) /
			   pdd;
	}

	if (*srt > 2.136857f)
	{
		r__1 = *srt;
		r__2 = *srt;
		xkp5 = sigm * 4.f * (r__1 * r__1 - 4.5661578384489996f) * (r__2 * r__2 - .066534591249000019f) / pdd;
	}
	if (*srt > 2.4294000000000002f)
	{
		r__1 = *srt;
		r__2 = *srt;
		xkp6 = sigm * 16.f * (r__1 * r__1 - 5.901984360000001f) * (r__2 * r__2 - .0011971599999999975f) / pdd;
	}
	if (*srt > 2.6374f)
	{
		r__1 = *srt;
		r__2 = *srt;
		xkp7 = sigm * 4.f * (r__1 * r__1 - 6.95587876f) * (r__2 * r__2 - .058854759999999964f) / pdd;
	}
	if (*srt > 2.7324000000000002f)
	{
		r__1 = *srt;
		r__2 = *srt;
		xkp8 = sigm * 4.f * (r__1 * r__1 - 7.4660097600000013f) * (r__2 * r__2 - .11397375999999994f) / pdd;
	}
L70:
	sig1 = 195.639f;
	sig2 = 372.378f;
	if (*srt > 1.518f)
	{
		r__1 = *srt;
		r__2 = *srt;
		pff = sqrt((r__1 * r__1 - 2.3043240000000003f) * (r__2 * r__2 -
														  .272484f));
		r__1 = *srt;
		xkp9 = sig1 * pff / sqrt(pdd) * 1.f / 32.f / 3.1415926f / (r__1 * r__1);
		if (*srt > 1.915f)
		{
			r__1 = *srt;
			r__2 = *srt;
			pff = sqrt((r__1 * r__1 - 3.6672250000000002f) * (r__2 * r__2 -
															  .015625f));
			r__1 = *srt;
			xkp10 = sig2 * pff / sqrt(pdd) * 3.f / 32.f / 3.1415926f / (r__1 * r__1);
		}
	}
	sigpik = 0.f;
	if (*srt > 1.6640000000000001f)
	{
		r__1 = *srt;
		r__2 = *srt;
		r__3 = *srt;
		r__4 = *px;
		r__5 = *py;
		r__6 = *pz;
		sigpik = *srhoks * 9.f * (r__1 * r__1 - .015625f) * (r__2 * r__2 - 2.7722250000000002f) / 4 / (r__3 * r__3) / (r__4 * r__4 + r__5 * r__5 + r__6 * r__6);
		if (*srt > 1.677f)
		{
			sigpik = sigpik * 12.f / 9.f;
		}
	}

	sigkp = xkp0 + xkp1 + xkp2 + xkp3 + xkp4 + xkp5 + xkp6 + xkp7 + xkp8 +
			xkp9 + xkp10 + sigpik;
	*icase = 0;
	dskn = sqrt(sigkp / 3.1415926f / 10.f);
	dsknr = dskn + .1f;
	distce_(i1, i2, &dsknr, &dskn, &input1_1.dt, ec, srt, &ic, px, py, pz);
	if (ic == -1)
	{
		return 0;
	}

	randu = ranart_(&rndf77_1.nseed) * sigkp;
	xkp1 = xkp0 + xkp1;
	xkp2 = xkp1 + xkp2;
	xkp3 = xkp2 + xkp3;
	xkp4 = xkp3 + xkp4;
	xkp5 = xkp4 + xkp5;
	xkp6 = xkp5 + xkp6;
	xkp7 = xkp6 + xkp7;
	xkp8 = xkp7 + xkp8;
	xkp9 = xkp8 + xkp9;
	xkp10 = xkp9 + xkp10;

	if (randu <= xkp0)
	{
		*icase = 1;
		return 0;
	}
	else
	{
		*icase = 2;
		if (randu <= xkp1)
		{
			*lbp1 = -14;
			*lbp2 = (int)(ranart_(&rndf77_1.nseed) * 2) + 1;
			*emm1 = 1.1157f;
			*emm2 = .939457f;
			goto L60;
		}
		else if (randu <= xkp2)
		{
			*lbp1 = -14;
			*lbp2 = (int)(ranart_(&rndf77_1.nseed) * 4) + 6;
			*emm1 = 1.1157f;
			*emm2 = 1.232f;
			goto L60;
		}
		else if (randu <= xkp3)
		{
			*lbp1 = -14;
			*lbp2 = (int)(ranart_(&rndf77_1.nseed) * 2) + 10;
			*emm1 = 1.1157f;
			*emm2 = 1.44f;
			goto L60;
		}
		else if (randu <= xkp4)
		{
			*lbp1 = -14;
			*lbp2 = (int)(ranart_(&rndf77_1.nseed) * 2) + 12;
			*emm1 = 1.1157f;
			*emm2 = 1.535f;
			goto L60;
		}
		else if (randu <= xkp5)
		{
			*lbp1 = -15 - (int)(ranart_(&rndf77_1.nseed) * 3);
			*lbp2 = (int)(ranart_(&rndf77_1.nseed) * 2) + 1;
			*emm1 = 1.1974f;
			*emm2 = .939457f;
			goto L60;
		}
		else if (randu <= xkp6)
		{
			*lbp1 = -15 - (int)(ranart_(&rndf77_1.nseed) * 3);
			*lbp2 = (int)(ranart_(&rndf77_1.nseed) * 4) + 6;
			*emm1 = 1.1974f;
			*emm2 = 1.232f;
			goto L60;
		}
		else if (randu < xkp7)
		{
			*lbp1 = -15 - (int)(ranart_(&rndf77_1.nseed) * 3);
			*lbp2 = (int)(ranart_(&rndf77_1.nseed) * 2) + 10;
			*emm1 = 1.1974f;
			*emm2 = 1.44f;
			goto L60;
		}
		else if (randu < xkp8)
		{
			*lbp1 = -15 - (int)(ranart_(&rndf77_1.nseed) * 3);
			*lbp2 = (int)(ranart_(&rndf77_1.nseed) * 2) + 12;
			*emm1 = 1.1974f;
			*emm2 = 1.535f;
			goto L60;
		}
		else if (randu < xkp9)
		{
			*icase = 3;
			*lbp1 = 29;
			*lbp2 = 23;
			*emm1 = 1.02f;
			*emm2 = .498f;
			if (ee_1.lb[*i1 - 1] == 21 || ee_1.lb[*i2 - 1] == 21)
			{
				*lbp2 = 21;
				*icase = -3;
			}
			goto L60;
		}
		else if (randu < xkp10)
		{
			*icase = 4;
			*lbp1 = 29;
			*lbp2 = 30;
			*emm1 = 1.02f;
			*emm2 = .895f;
			if (ee_1.lb[*i1 - 1] == 21 || ee_1.lb[*i2 - 1] == 21)
			{
				*lbp2 = -30;
				*icase = -4;
			}
			goto L60;
		}
		else
		{
			*icase = 5;
			*lbp1 = (int)(ranart_(&rndf77_1.nseed) * 3) + 25;
			*lbp2 = 30;
			*emm1 = .769f;
			*emm2 = .895f;
			if (*srt > 1.677f && ranart_(&rndf77_1.nseed) < .25f)
			{
				*lbp1 = 28;
				*emm1 = .782f;
			}
			if (ee_1.lb[*i1 - 1] == 21 || ee_1.lb[*i2 - 1] == 21)
			{
				*lbp2 = -30;
				*icase = -5;
			}
		}
	}

L60:
	if (*icase == 2 && (ee_1.lb[*i1 - 1] == 21 || ee_1.lb[*i2 - 1] == 21))
	{
		*lbp1 = -(*lbp1);
		*lbp2 = -(*lbp2);
	}
	px0 = *px;
	py0 = *py;
	pz0 = *pz;
	r__2 = *srt;
	r__3 = *emm1;
	r__4 = *emm2;
	r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
	r__5 = *emm1 * *emm2;
	pr2 = r__1 * r__1 - r__5 * r__5 * 4.f;
	if (pr2 <= 0.f)
	{
		pr2 = 1e-9f;
	}
	pr = sqrt(pr2) / (*srt * 2.f);
	c1 = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
	t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
	r__1 = c1;
	s1 = sqrt(1.f - r__1 * r__1);
	ct1 = cos(t1);
	st1 = sin(t1);
	*pz = pr * c1;
	*px = pr * s1 * ct1;
	*py = pr * s1 * st1;
	return 0;
}

int crkphi_(float *px, float *py, float *pz, float *ec, float *srt, int *iblock, float *emm1, float *emm2, int *lbp1, int *lbp2, int *i1, int *i2, int *ikk, int *icase, float *rrkk, float *prkk)
{
	float r__1, r__2, r__3, r__4, r__5;
	double d__1;

	double pow_dd(double *, double *), sqrt(double), cos(double), sin(double);

	static float c1, s1, t1;
	static int ic;
	static float pr;
	static int lb1, lb2;
	static float ct1, pr2, px0, py0, pz0, st1, pii, sig, dnr, sig1, sig2, sig3,
		xsk1, srr1, srr2, srr3, xsk2, xsk3, xsk4, xsk5, xsk6, xsk7, xsk8,
		xsk9, sigm, dskn, xsk10, xsk11, ranx, srri, srrt, sigm0, prkk0,
		rrkk0, sigks, dsknr, sigks1, sigks2, sigks3, sigks4;
	extern int distce_(int *, int *, float *, float *,
					   float *, float *, float *, int *, float *, float *, float *),
		crkkpi_(int *, int *, float *, float *, float *, float *,
				float *, float *, float *, float *, float *, float *, float *, float *,
				int *, int *, int *, float *, float *);
	extern double ranart_(int *);
	extern int xkkann_(float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *, float *,
					   float *),
		crkspi_(int *, int *, float *, float *, float *,
				float *, float *, int *, int *, int *, float *, float *),
		xkksan_(int *, int *, float *, float *, float *, float *,
				float *, float *, float *);

	lb1 = ee_1.lb[*i1 - 1];
	lb2 = ee_1.lb[*i2 - 1];
	*icase = 0;
	if (*srt < 1.15496f)
	{
		sig1 = 0.f;
		sig2 = 0.f;
		sig3 = 0.f;
	}
	else
	{

		if (lb1 == 23 && lb2 == 21 || lb2 == 23 && lb1 == 21)
		{
			dnr = 4.f;
			*ikk = 2;
		}
		else if (lb1 == 21 && lb2 == 30 || lb2 == 21 && lb1 == 30 || lb1 == 23 && lb2 == -30 || lb2 == 23 && lb1 == -30)
		{
			dnr = 12.f;
			*ikk = 1;
		}
		else
		{
			dnr = 36.f;
			*ikk = 0;
		}
		sig1 = 0.f;
		sig2 = 0.f;
		sig3 = 0.f;
		srri = cc_1.e[*i1 - 1] + cc_1.e[*i2 - 1];
		srr1 = 1.15496f;
		srr2 = 1.8019000000000001f;
		srr3 = 1.79f;

		r__1 = *srt;
		r__2 = cc_1.e[*i1 - 1] + cc_1.e[*i2 - 1];
		r__3 = *srt;
		r__4 = cc_1.e[*i1 - 1] - cc_1.e[*i2 - 1];
		pii = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4);
		srrt = *srt - dmax(srri, srr1);
		if (srrt < .3f && srrt > .01f)
		{
			d__1 = (double)srrt;
			sig = 1.69f / (pow_dd(&d__1, &c_b529) - .407f);
		}
		else
		{
			d__1 = (double)srrt;
			sig = pow_dd(&d__1, &c_b530) * .008f + 3.74f;
		}
		r__1 = *srt;
		r__2 = *srt;
		sig1 = sig * (9.f / dnr) * (r__1 * r__1 - 1.3339326015999999f) * (r__2 * r__2 - .78329580160000012f) / pii;
		if (*srt > 1.8019000000000001f)
		{
			srrt = *srt - dmax(srri, srr2);
			if (srrt < .3f && srrt > .01f)
			{
				d__1 = (double)srrt;
				sig = 1.69f / (pow_dd(&d__1, &c_b529) - .407f);
			}
			else
			{
				d__1 = (double)srrt;
				sig = pow_dd(&d__1, &c_b530) * .008f + 3.74f;
			}
			r__1 = *srt;
			r__2 = *srt;
			sig2 = sig * (9.f / dnr) * (r__1 * r__1 - 3.24684361f) * (r__2 * r__2 - .05669160999999999f) / pii;
		}
		if (*srt > 1.79f)
		{
			srrt = *srt - dmax(srri, srr3);
			if (srrt < .3f && srrt > .01f)
			{
				d__1 = (double)srrt;
				sig = 1.69f / (pow_dd(&d__1, &c_b529) - .407f);
			}
			else
			{
				d__1 = (double)srrt;
				sig = pow_dd(&d__1, &c_b530) * .008f + 3.74f;
			}
			r__1 = *srt;
			r__2 = *srt;
			sig3 = sig * (27.f / dnr) * (r__1 * r__1 - 3.2040999999999999f) *
				   (r__2 * r__2 - .0625f) / pii;
		}
	}
	rrkk0 = *rrkk;
	prkk0 = *prkk;
	sigm = 0.f;
	if (lb1 == 23 && lb2 == 21 || lb2 == 23 && lb1 == 21)
	{
		xkkann_(srt, &xsk1, &xsk2, &xsk3, &xsk4, &xsk5, &xsk6, &xsk7, &xsk8, &xsk9, &xsk10, &xsk11, &sigm, &rrkk0);
	}
	else if (lb1 == 21 && lb2 == 30 || lb2 == 21 && lb1 == 30 || lb1 == 23 && lb2 == -30 || lb2 == 23 && lb1 == -30)
	{
		xkksan_(i1, i2, srt, &sigks1, &sigks2, &sigks3, &sigks4, &sigm, &prkk0);
	}
	else
	{
	}

	sigm0 = sigm;
	sigks = sig1 + sig2 + sig3 + sigm;
	dskn = sqrt(sigks / 3.1415926f / 10.f);
	dsknr = dskn + .1f;
	distce_(i1, i2, &dsknr, &dskn, &input1_1.dt, ec, srt, &ic, px, py, pz);
	if (ic == -1)
	{
		return 0;
	}
	*icase = 1;
	ranx = ranart_(&rndf77_1.nseed);
	*lbp1 = 29;
	*emm1 = 1.02f;
	if (ranx <= sig1 / sigks)
	{
		*lbp2 = (int)(ranart_(&rndf77_1.nseed) * 3) + 3;
		*emm2 = .13496f;
	}
	else if (ranx <= (sig1 + sig2) / sigks)
	{
		*lbp2 = 28;
		*emm2 = .7819f;
	}
	else if (ranx <= (sig1 + sig2 + sig3) / sigks)
	{
		*lbp2 = (int)(ranart_(&rndf77_1.nseed) * 3) + 25;
		*emm2 = .77f;
	}
	else
	{
		if (lb1 == 23 && lb2 == 21 || lb2 == 23 && lb1 == 21)
		{
			crkkpi_(i1, i2, &xsk1, &xsk2, &xsk3, &xsk4, &xsk5, &xsk6, &xsk7, &xsk8, &xsk9, &xsk10, &xsk11, &sigm0, iblock, lbp1, lbp2,
					emm1, emm2);
		}
		else if (lb1 == 21 && lb2 == 30 || lb2 == 21 && lb1 == 30 || lb1 == 23 && lb2 == -30 || lb2 == 23 && lb1 == -30)
		{
			crkspi_(i1, i2, &sigks1, &sigks2, &sigks3, &sigks4, &sigm0,
					iblock, lbp1, lbp2, emm1, emm2);
		}
		else
		{
		}
	}

	px0 = *px;
	py0 = *py;
	pz0 = *pz;
	r__2 = *srt;
	r__3 = *emm1;
	r__4 = *emm2;
	r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
	r__5 = *emm1 * *emm2;
	pr2 = r__1 * r__1 - r__5 * r__5 * 4.f;
	if (pr2 <= 0.f)
	{
		pr2 = 1e-9f;
	}
	pr = sqrt(pr2) / (*srt * 2.f);
	c1 = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
	t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
	r__1 = c1;
	s1 = sqrt(1.f - r__1 * r__1);
	ct1 = cos(t1);
	st1 = sin(t1);
	*pz = pr * c1;
	*px = pr * s1 * ct1;
	*py = pr * s1 * st1;
	return 0;
}

int crksph_(float *px, float *py, float *pz, float *ec, float *srt, float *emm1, float *emm2, int *lbp1, int *lbp2, int *i1, int *i2, int *ikkg, int *ikkl, int *iblock,
			int *icase, float *srhoks)
{
	float r__1, r__2, r__3, r__4, r__5;

	double sqrt(double), cos(double), sin(double);

	static float c1, s1, t1;
	static int ic;
	static float pr;
	static int lb1, lb2;
	static float ct1, pr2, px0, py0, pz0, st1, pff, pii, dnr, sig1, sig2,
		sig11, sig22, dskn, ranx, sigkm, sigks, dsknr, sigela;
	extern int distce_(int *, int *, float *, float *,
					   float *, float *, float *, int *, float *, float *, float *);
	extern double ranart_(int *);

	lb1 = ee_1.lb[*i1 - 1];
	lb2 = ee_1.lb[*i2 - 1];
	*icase = 0;
	sigela = 10.f;
	sigkm = 0.f;
	if (lb1 >= 25 && lb1 <= 28 || lb2 >= 25 && lb2 <= 28)
	{
		if (abs(lb1) == 30 || abs(lb2) == 30)
		{
			sigkm = *srhoks;
		}
		else if ((lb1 == 23 || lb1 == 21 || lb2 == 23 || lb2 == 21) && *srt > 1.03457f)
		{
			sigkm = *srhoks;
		}
	}
	if (*srt < 1.518f)
	{
		sig11 = 0.f;
		sig22 = 0.f;
	}
	else
	{
		if (abs(lb1) == 30 && (lb2 >= 3 && lb2 <= 5) || abs(lb2) == 30 && (lb1 >= 3 && lb1 <= 5))
		{
			dnr = 18.f;
			*ikkl = 0;
			*iblock = 225;
			sig1 = 2047.042f;
			sig2 = 1496.692f;
		}
		else if (lb1 == 23 || lb1 == 21 && (lb2 >= 25 && lb2 <= 27) || (lb2 == 23 || lb2 == 21 && (lb1 >= 25 && lb1 <= 27)))
		{
			dnr = 18.f;
			*ikkl = 1;
			*iblock = 224;
			sig1 = 526.702f;
			sig2 = 1313.96f;
		}
		else if (abs(lb1) == 30 && (lb2 >= 25 && lb2 <= 27) || abs(lb2) ==
																	   30 &&
																   (lb1 >= 25 && lb1 <= 27))
		{
			dnr = 54.f;
			*ikkl = 0;
			*iblock = 225;
			sig1 = 1371.257f;
			sig2 = 6999.84f;
		}
		else if ((lb1 == 23 || lb1 == 21) && lb2 == 28 || (lb2 == 23 || lb2 == 21) && lb1 == 28)
		{
			dnr = 6.f;
			*ikkl = 1;
			*iblock = 224;
			sig1 = 355.429f;
			sig2 = 440.558f;
		}
		else
		{
			dnr = 18.f;
			*ikkl = 0;
			*iblock = 225;
			sig1 = 482.292f;
			sig2 = 1698.903f;
		}
		sig11 = 0.f;
		sig22 = 0.f;
		r__1 = *srt;
		r__2 = cc_1.e[*i1 - 1] + cc_1.e[*i2 - 1];
		r__3 = *srt;
		r__4 = cc_1.e[*i1 - 1] - cc_1.e[*i2 - 1];
		pii = sqrt((r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4));
		r__1 = *srt;
		r__2 = *srt;
		pff = sqrt((r__1 * r__1 - 2.3043240000000003f) * (r__2 * r__2 -
														  .272484f));
		r__1 = *srt;
		sig11 = sig1 * pff / pii * 6.f / dnr / 32.f / 3.1415926f / (r__1 * r__1);

		if (*srt > 1.915f)
		{
			r__1 = *srt;
			r__2 = *srt;
			pff = sqrt((r__1 * r__1 - 3.6672250000000002f) * (r__2 * r__2 -
															  .015625f));
			r__1 = *srt;
			sig22 = sig2 * pff / pii * 18.f / dnr / 32.f / 3.1415926f / (r__1 * r__1);
		}
	}
	sigks = sig11 + sig22 + sigela + sigkm;

	dskn = sqrt(sigks / 3.1415926f / 10.f);
	dsknr = dskn + .1f;
	distce_(i1, i2, &dsknr, &dskn, &input1_1.dt, ec, srt, &ic, px, py, pz);
	if (ic == -1)
	{
		return 0;
	}
	*icase = 1;
	ranx = ranart_(&rndf77_1.nseed);
	if (ranx <= sigela / sigks)
	{
		*lbp1 = lb1;
		*emm1 = cc_1.e[*i1 - 1];
		*lbp2 = lb2;
		*emm2 = cc_1.e[*i2 - 1];
		*iblock = 111;
	}
	else if (ranx <= (sigela + sigkm) / sigks)
	{
		*lbp1 = (int)(ranart_(&rndf77_1.nseed) * 3) + 3;
		*emm1 = .14f;
		if (lb1 == 23 || lb2 == 23)
		{
			*lbp2 = 30;
			*emm2 = .895f;
		}
		else if (lb1 == 21 || lb2 == 21)
		{
			*lbp2 = -30;
			*emm2 = .895f;
		}
		else if (lb1 == 30 || lb2 == 30)
		{
			*lbp2 = 23;
			*emm2 = .498f;
		}
		else
		{
			*lbp2 = 21;
			*emm2 = .498f;
		}
		*iblock = 112;
	}
	else if (ranx <= (sigela + sigkm + sig11) / sigks)
	{
		*lbp2 = 23;
		*emm2 = .498f;
		*ikkg = 1;
		if (lb1 == 21 || lb2 == 21 || lb1 == -30 || lb2 == -30)
		{
			*lbp2 = 21;
			*iblock += -100;
		}
		*lbp1 = 29;
		*emm1 = 1.02f;
	}
	else
	{
		*lbp2 = 30;
		*emm2 = .895f;
		*ikkg = 0;
		*iblock += 2;
		if (lb1 == 21 || lb2 == 21 || lb1 == -30 || lb2 == -30)
		{
			*lbp2 = -30;
			*iblock += -100;
		}
		*lbp1 = 29;
		*emm1 = 1.02f;
	}

	px0 = *px;
	py0 = *py;
	pz0 = *pz;
	r__2 = *srt;
	r__3 = *emm1;
	r__4 = *emm2;
	r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
	r__5 = *emm1 * *emm2;
	pr2 = r__1 * r__1 - r__5 * r__5 * 4.f;
	if (pr2 <= 0.f)
	{
		pr2 = 1e-9f;
	}
	pr = sqrt(pr2) / (*srt * 2.f);
	c1 = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
	t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
	r__1 = c1;
	s1 = sqrt(1.f - r__1 * r__1);
	ct1 = cos(t1);
	st1 = sin(t1);
	*pz = pr * c1;
	*px = pr * s1 * ct1;
	*py = pr * s1 * st1;
	return 0;
}

int bbkaon_(int *ic, float *srt, float *px, float *py, float *pz, float *ana, float *plx, float *ply, float *plz, float *ala, float *pkx,
			float *pky, float *pkz, int *icou1)
{
	float r__1, r__2, r__3, r__4, r__5;

	double sqrt(double), cos(double), sin(double);

	static float t1, t2, ga, ek, el, en, cs, pi, bx, pk, by, bz, pn, ss, dm1,
		pn2, aka, fai, eln, csn, ssn, fain, elnc, dmax__, prob, pmax;
	static int ntry;
	static float pbeta;
	extern double fkaon_(float *, float *), rmass_(float *, int *);
	static float trans0;
	extern double ranart_(int *);

	pi = 3.1415962f;
	*icou1 = 0;
	aka = .498f;
	*ala = 1.116f;
	if (*ic == 2 || *ic == 4)
	{
		*ala = 1.197f;
	}
	*ana = .939f;
	if (*ic > 2)
	{
		dmax__ = *srt - aka - *ala - .02f;
		dm1 = rmass_(&dmax__, &input1_1.iseed);
		*ana = dm1;
	}
	t1 = aka + *ana + *ala;
	t2 = *ana + *ala - aka;
	if (*srt <= t1)
	{
		*icou1 = -1;
		return 0;
	}
	r__1 = *srt;
	r__2 = t1;
	r__3 = *srt;
	r__4 = t2;
	pmax = sqrt((r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4)) / (*srt * 2.f);
	if (pmax == 0.f)
	{
		pmax = 1e-9f;
	}
	ntry = 0;
L1:
	pk = pmax * ranart_(&rndf77_1.nseed);
	++ntry;
	prob = fkaon_(&pk, &pmax);
	if (prob < ranart_(&rndf77_1.nseed) && ntry <= 40)
	{
		goto L1;
	}
	cs = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
	r__1 = cs;
	ss = sqrt(1.f - r__1 * r__1);
	fai = ranart_(&rndf77_1.nseed) * 6.2800000000000002f;
	*pkx = pk * ss * cos(fai);
	*pky = pk * ss * sin(fai);
	*pkz = pk * cs;
	r__1 = aka;
	r__2 = pk;
	ek = sqrt(r__1 * r__1 + r__2 * r__2);
	eln = *srt - ek;
	if (eln <= 0.f)
	{
		*icou1 = -1;
		return 0;
	}
	bx = -(*pkx) / eln;
	by = -(*pky) / eln;
	bz = -(*pkz) / eln;
	r__1 = bx;
	r__2 = by;
	r__3 = bz;
	ga = 1.f / sqrt(1.f - r__1 * r__1 - r__2 * r__2 - r__3 * r__3);
	elnc = eln / ga;
	r__2 = elnc;
	r__3 = *ana;
	r__4 = *ala;
	r__1 = (r__2 * r__2 + r__3 * r__3 - r__4 * r__4) / (elnc * 2.f);
	r__5 = *ana;
	pn2 = r__1 * r__1 - r__5 * r__5;
	if (pn2 <= 0.f)
	{
		pn2 = 1e-9f;
	}
	pn = sqrt(pn2);
	csn = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
	r__1 = csn;
	ssn = sqrt(1.f - r__1 * r__1);
	fain = ranart_(&rndf77_1.nseed) * 6.2800000000000002f;
	*px = pn * ssn * cos(fain);
	*py = pn * ssn * sin(fain);
	*pz = pn * csn;
	r__1 = *ana;
	en = sqrt(r__1 * r__1 + pn2);
	*plx = -(*px);
	*ply = -(*py);
	*plz = -(*pz);
	pbeta = *px * bx + *py * by + *pz * bz;
	trans0 = ga * (ga * pbeta / (ga + 1.f) + en);
	*px = bx * trans0 + *px;
	*py = by * trans0 + *py;
	*pz = bz * trans0 + *pz;
	r__1 = *ala;
	r__2 = *plx;
	r__3 = *ply;
	r__4 = *plz;
	el = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	pbeta = *plx * bx + *ply * by + *plz * bz;
	trans0 = ga * (ga * pbeta / (ga + 1.f) + el);
	*plx = bx * trans0 + *plx;
	*ply = by * trans0 + *ply;
	*plz = bz * trans0 + *plz;
	return 0;
}

double pipik_(float *srt)
{
	static float xarray[5] = {.001f, .7f, 1.5f, 1.7f, 2.f};
	static float earray[5] = {1.f, 1.2f, 1.6f, 2.f, 2.4f};

	float ret_val;

	double log(double), exp(double);

	static int ie;
	static float xmin, ymin, xmax, ymax, pmass;

	pmass = .9383f;
	ret_val = 0.f;
	if (*srt <= 1.f)
	{
		return ret_val;
	}
	if (*srt > 2.4f)
	{
		ret_val = 1.f;
		return ret_val;
	}
	if (*srt < earray[0])
	{
		ret_val = xarray[0] / 2.f;
		return ret_val;
	}

	for (ie = 1; ie <= 5; ++ie)
	{
		if (earray[ie - 1] == *srt)
		{
			ret_val = xarray[ie - 1];
			goto L10;
		}
		else if (earray[ie - 1] > *srt)
		{
			ymin = log(xarray[ie - 2]);
			ymax = log(xarray[ie - 1]);
			xmin = log(earray[ie - 2]);
			xmax = log(earray[ie - 1]);
			ret_val = exp(ymin + (log(*srt) - xmin) * (ymax - ymin) / (xmax - xmin));
			goto L10;
		}
	}
L10:
	ret_val /= 2.f;
	return ret_val;
}

double pionpp_(float *srt)
{
	float ret_val, r__1, r__2, r__3, r__4, r__5;
	double d__1, d__2;

	double sqrt(double), log(double), pow_dd(double *, double *);

	static float a, b, c__, d__, an, plab, pmin, pmax, pmass, pmass1;

	pmass = .14f;
	pmass1 = .938f;
	ret_val = 1e-5f;
	if (*srt <= 1.22f)
	{
		return ret_val;
	}
	r__2 = *srt;
	r__3 = pmass;
	r__4 = pmass1;
	r__1 = (r__2 * r__2 - r__3 * r__3 - r__4 * r__4) / (pmass1 * 2.f);
	r__5 = pmass;
	plab = sqrt(r__1 * r__1 - r__5 * r__5);
	pmin = .3f;
	pmax = 25.f;
	if (plab > pmax)
	{
		ret_val = 2.f;
		return ret_val;
	}
	if (plab < pmin)
	{
		ret_val = 0.f;
		return ret_val;
	}
	a = 24.3f;
	b = -12.3f;
	c__ = .324f;
	an = -1.91f;
	d__ = -2.44f;
	d__1 = (double)plab;
	d__2 = (double)an;
	r__1 = log(plab);
	ret_val = a + b * pow_dd(&d__1, &d__2) + c__ * (r__1 * r__1) + d__ * log(plab);
	if (ret_val <= 0.f)
	{
		ret_val = 0.f;
	}
	ret_val /= 10.f;
	return ret_val;
}

double pipp1_(float *srt)
{
	float ret_val, r__1, r__2, r__3, r__4, r__5;
	double d__1, d__2;

	double sqrt(double), log(double), pow_dd(double *, double *);

	static float a, b, c__, d__, an, plab, pmin, pmax, pmass, pmass1;

	pmass = .14f;
	pmass1 = .938f;
	ret_val = 1e-4f;
	if (*srt <= 1.22f)
	{
		return ret_val;
	}
	r__2 = *srt;
	r__3 = pmass;
	r__4 = pmass1;
	r__1 = (r__2 * r__2 - r__3 * r__3 - r__4 * r__4) / (pmass1 * 2.f);
	r__5 = pmass;
	plab = sqrt(r__1 * r__1 - r__5 * r__5);
	pmin = .3f;
	pmax = 25.f;
	if (plab > pmax)
	{
		ret_val = 2.f;
		return ret_val;
	}
	if (plab < pmin)
	{
		ret_val = 0.f;
		return ret_val;
	}
	a = 26.6f;
	b = -7.18f;
	c__ = .327f;
	an = -1.86f;
	d__ = -2.81f;
	d__1 = (double)plab;
	d__2 = (double)an;
	r__1 = log(plab);
	ret_val = a + b * pow_dd(&d__1, &d__2) + c__ * (r__1 * r__1) + d__ * log(plab);
	if (ret_val <= 0.f)
	{
		ret_val = 0.f;
	}
	ret_val /= 10.f;
	return ret_val;
}

double xrho_(float *srt)
{
	float ret_val, r__1;

	static float es, trho, xrho0, esmin, pmass, rmass;

	pmass = .9383f;
	rmass = .77f;
	trho = .151f;
	ret_val = 1e-9f;
	if (*srt <= 2.67f)
	{
		return ret_val;
	}
	esmin = rmass + 1.8766f - trho / 2.f;
	es = *srt;
	r__1 = es - esmin;
	xrho0 = (es - esmin) * .24f / (r__1 * r__1 + 1.4f);
	ret_val = xrho0 * 3.f;
	return ret_val;
}

double omega_(float *srt)
{
	float ret_val, r__1;

	static float es, esmin, omass, pmass, tomega;

	pmass = .9383f;
	omass = .782f;
	tomega = .0084f;
	ret_val = 1e-8f;
	if (*srt <= 2.68f)
	{
		return ret_val;
	}
	esmin = omass + 1.8766f - tomega / 2.f;
	es = *srt;
	r__1 = es - esmin;
	ret_val = (es - esmin) * .36f / (r__1 * r__1 + 1.25f);
	return ret_val;
}

double twopi_(float *srt)
{
	static float xarray[19] = {3e-6f, 1.87f, 11.f, 14.9f, 9.35f, 7.65f, 4.62f, 3.45f,
							  2.41f, 1.85f, 1.65f, 1.5f, 1.32f, 1.17f, 1.16f, 1.f, .856f, .745f, 3e-6f};
	static float earray[19] = {1.22f, 1.47f, 1.72f, 1.97f, 2.22f, 2.47f, 2.72f,
							  2.97f, 3.22f, 3.47f, 3.72f, 3.97f, 4.22f, 4.47f, 4.72f, 4.97f, 5.22f, 5.47f,
							  5.72f};

	float ret_val;

	double log(double), exp(double);

	static int ie;
	static float plab, xmin, ymin, xmax, ymax, pmass, pmass1;

	pmass = .14f;
	pmass1 = .938f;
	ret_val = 1e-6f;
	if (*srt <= 1.22f)
	{
		return ret_val;
	}
	plab = *srt;
	if (plab < earray[0])
	{
		ret_val = 1e-5f;
		return ret_val;
	}

	for (ie = 1; ie <= 19; ++ie)
	{
		if (earray[ie - 1] == plab)
		{
			ret_val = xarray[ie - 1];
			return ret_val;
		}
		else if (earray[ie - 1] > plab)
		{
			ymin = log(xarray[ie - 2]);
			ymax = log(xarray[ie - 1]);
			xmin = log(earray[ie - 2]);
			xmax = log(earray[ie - 1]);
			ret_val = exp(ymin + (log(plab) - xmin) * (ymax - ymin) / (xmax - xmin));
			return ret_val;
		}
	}
	return ret_val;
}

double threpi_(float *srt)
{
	static float xarray[15] = {8e-6f, 6.1999999e-5f, 1.88194f, 5.02569f,
							  11.80154f, 13.92114f, 15.07308f, 11.79571f, 11.53772f, 10.01197f,
							  9.792673f, 9.465264f, 8.97049f, 7.944254f, 6.88632f};
	static float earray[15] = {1.22f, 1.47f, 1.72f, 1.97f, 2.22f, 2.47f, 2.72f,
							  2.97f, 3.22f, 3.47f, 3.72f, 3.97f, 4.22f, 4.47f, 4.72f};

	float ret_val;

	double log(double), exp(double);

	static int ie;
	static float plab, xmin, ymin, xmax, ymax, pmass, pmass1;

	pmass = .14f;
	pmass1 = .938f;
	ret_val = 1e-6f;
	if (*srt <= 1.36f)
	{
		return ret_val;
	}
	plab = *srt;
	if (plab < earray[0])
	{
		ret_val = 1e-5f;
		return ret_val;
	}

	for (ie = 1; ie <= 15; ++ie)
	{
		if (earray[ie - 1] == plab)
		{
			ret_val = xarray[ie - 1];
			return ret_val;
		}
		else if (earray[ie - 1] > plab)
		{
			ymin = log(xarray[ie - 2]);
			ymax = log(xarray[ie - 1]);
			xmin = log(earray[ie - 2]);
			xmax = log(earray[ie - 1]);
			ret_val = exp(ymin + (log(plab) - xmin) * (ymax - ymin) / (xmax - xmin));
			return ret_val;
		}
	}
	return ret_val;
}

double fourpi_(float *srt)
{
	static float xarray[10] = {1e-4f, 1.986597f, 6.411932f, 7.636956f, 9.598362f,
							  9.88974f, 10.24317f, 10.80138f, 11.86988f, 12.83925f};
	static float earray[10] = {2.468f, 2.718f, 2.968f, 3.22f, 3.47f, 3.72f, 3.97f,
							  4.22f, 4.47f, 4.72f};

	float ret_val;

	double log(double), exp(double);

	static int ie;
	static float plab, xmin, ymin, xmax, ymax, pmass, pmass1;

	pmass = .14f;
	pmass1 = .938f;
	ret_val = 1e-6f;
	if (*srt <= 1.52f)
	{
		return ret_val;
	}
	plab = *srt;
	if (plab < earray[0])
	{
		ret_val = 1e-5f;
		return ret_val;
	}

	for (ie = 1; ie <= 10; ++ie)
	{
		if (earray[ie - 1] == plab)
		{
			ret_val = xarray[ie - 1];
			return ret_val;
		}
		else if (earray[ie - 1] > plab)
		{
			ymin = log(xarray[ie - 2]);
			ymax = log(xarray[ie - 1]);
			xmin = log(earray[ie - 2]);
			xmax = log(earray[ie - 1]);
			ret_val = exp(ymin + (log(plab) - xmin) * (ymax - ymin) / (xmax - xmin));
			return ret_val;
		}
	}
	return ret_val;
}

double reab_(int *i1, int *i2, float *srt, int *ictrl)
{
	float ret_val, r__1, r__2, r__3, r__4, r__5;

	static float ed;
	static int lb1, lb2;
	static float pin2, xpro, arho1, pout2;
	extern double twopi_(float *);
	static float factor;
	extern double threpi_(float *), fourpi_(float *);

	lb1 = ee_1.lb[*i1 - 1];
	lb2 = ee_1.lb[*i2 - 1];
	ret_val = 0.f;
	if (*ictrl == 1 && *srt <= 1.238f)
	{
		return ret_val;
	}
	if (*ictrl == 3 && *srt <= 1.8799999999999999f)
	{
		return ret_val;
	}
	r__2 = *srt;
	r__1 = (r__2 * r__2 + .019600000000000003f - .87984399999999985f) / (*srt * 2.f);
	pin2 = r__1 * r__1 - .019600000000000003f;
	if (pin2 <= 0.f)
	{
		return ret_val;
	}
	if (*ictrl == 1)
	{
		if (cc_1.e[*i1 - 1] > 1.f)
		{
			ed = cc_1.e[*i1 - 1];
		}
		else
		{
			ed = cc_1.e[*i2 - 1];
		}
		r__2 = *srt;
		r__3 = ed;
		r__1 = (r__2 * r__2 + .019600000000000003f - r__3 * r__3) / (*srt *
																	 2.f);
		pout2 = r__1 * r__1 - .019600000000000003f;
		if (pout2 <= 0.f)
		{
			return ret_val;
		}
		xpro = twopi_(srt) / 10.f;
		factor = .33333333333333331f;
		if (lb1 == 8 && lb2 == 5 || lb1 == 5 && lb2 == 8 || (lb1 == -8 && lb2 == 3 || lb1 == 3 && lb2 == -8))
		{
			factor = .25f;
		}
		if (abs(lb1) >= 10 && abs(lb1) <= 13 || abs(lb2) >= 10 && abs(lb2) <=
																	  13)
		{
			factor = 1.f;
		}
		ret_val = factor * pin2 / pout2 * xpro;
		return ret_val;
	}
	if (*ictrl == 2)
	{
		if (ee_1.lb[*i2 - 1] >= 25)
		{
			ed = cc_1.e[*i1 - 1];
			arho1 = cc_1.e[*i2 - 1];
		}
		else
		{
			ed = cc_1.e[*i2 - 1];
			arho1 = cc_1.e[*i1 - 1];
		}
		if (*srt <= arho1 + 1.0779999999999998f + .02f)
		{
			return ret_val;
		}
		r__2 = *srt;
		r__3 = arho1;
		r__4 = ed;
		r__1 = (r__2 * r__2 + r__3 * r__3 - r__4 * r__4) / (*srt * 2.f);
		r__5 = arho1;
		pout2 = r__1 * r__1 - r__5 * r__5;
		if (pout2 <= 0.f)
		{
			return ret_val;
		}
		xpro = threpi_(srt) / 10.f;
		factor = .33333333333333331f;
		if (lb1 == 8 && lb2 == 27 || lb1 == 27 && lb2 == 8 || (lb1 == -8 && lb2 == 25 || lb1 == 25 && lb2 == -8))
		{
			factor = .25f;
		}
		if (abs(lb1) >= 10 && abs(lb1) <= 13 || abs(lb2) >= 10 && abs(lb2) <=
																	  13)
		{
			factor = 1.f;
		}
		ret_val = factor * pin2 / pout2 * xpro;
		return ret_val;
	}
	if (*ictrl == 3)
	{
		if (cc_1.e[*i1 - 1] > 1.f)
		{
			ed = cc_1.e[*i1 - 1];
		}
		if (cc_1.e[*i2 - 1] > 1.f)
		{
			ed = cc_1.e[*i2 - 1];
		}
		r__2 = *srt;
		r__3 = ed;
		r__1 = (r__2 * r__2 + .61152400000000007f - r__3 * r__3) / (*srt *
																	2.f);
		pout2 = r__1 * r__1 - .61152400000000007f;
		if (pout2 <= 0.f)
		{
			return ret_val;
		}
		xpro = fourpi_(srt) / 10.f;
		factor = .16666666666666666f;
		if (abs(lb1) >= 10 && abs(lb1) <= 13 || abs(lb2) >= 10 && abs(lb2) <=
																	  13)
		{
			factor = .33333333333333331f;
		}
		ret_val = factor * pin2 / pout2 * xpro;
	}
	return ret_val;
}

double reab2d_(int *i1, int *i2, float *srt)
{
	int i__1;
	float ret_val, r__1, r__2, r__3, r__4, r__5;

	static float ed1, ed2;
	static int lb1, lb2;
	static float pin2;
	extern double x2pi_(float *);
	static float xpro, pout2, factor;

	ret_val = 0.f;
	lb1 = (i__1 = ee_1.lb[*i1 - 1], abs(i__1));
	lb2 = (i__1 = ee_1.lb[*i2 - 1], abs(i__1));
	ed1 = cc_1.e[*i1 - 1];
	ed2 = cc_1.e[*i2 - 1];
	r__1 = *srt / 2.f;
	pin2 = r__1 * r__1 - .87984399999999985f;
	r__2 = *srt;
	r__3 = ed1;
	r__4 = ed2;
	r__1 = (r__2 * r__2 + r__3 * r__3 - r__4 * r__4) / (*srt * 2.f);
	r__5 = ed1;
	pout2 = r__1 * r__1 - r__5 * r__5;
	if (pout2 <= 0.f)
	{
		return ret_val;
	}
	xpro = x2pi_(srt);
	factor = .25f;
	if (lb1 >= 10 && lb1 <= 13 && (lb2 >= 10 && lb2 <= 13))
	{
		factor = 1.f;
	}
	if (lb1 >= 6 && lb1 <= 9 && (lb2 > 10 && lb2 <= 13))
	{
		factor = .5f;
	}
	if (lb2 >= 6 && lb2 <= 9 && (lb1 > 10 && lb1 <= 13))
	{
		factor = .5f;
	}
	ret_val = factor * pin2 / pout2 * xpro;
	return ret_val;
}

int rotate_(float *px0, float *py0, float *pz0, float *px, float *py, float *pz)
{
	float r__1, r__2, r__3;

	double sqrt(double), atan2(double, double), cos(double),
		sin(double);

	static float c1, c2, s1, s2, t2, t1, pr, ss, ct1, ct2, pr0, st2, st1;

	r__1 = *px0;
	r__2 = *py0;
	r__3 = *pz0;
	pr0 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
	if (pr0 == 0.f)
	{
		pr0 = 1e-8f;
	}
	c2 = *pz0 / pr0;
	if (*px0 == 0.f && *py0 == 0.f)
	{
		t2 = 0.f;
	}
	else
	{
		t2 = atan2(*py0, *px0);
	}
	r__1 = c2;
	s2 = sqrt(1.f - r__1 * r__1);
	ct2 = cos(t2);
	st2 = sin(t2);
	r__1 = *px;
	r__2 = *py;
	r__3 = *pz;
	pr = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
	if (pr == 0.f)
	{
		pr = 1e-7f;
	}
	c1 = *pz / pr;
	if (*px == 0.f && *py == 0.f)
	{
		t1 = 0.f;
	}
	else
	{
		t1 = atan2(*py, *px);
	}
	r__1 = c1;
	s1 = sqrt(1.f - r__1 * r__1);
	ct1 = cos(t1);
	st1 = sin(t1);
	ss = c2 * s1 * ct1 + s2 * c1;
	*px = pr * (ss * ct2 - s1 * st1 * st2);
	*py = pr * (ss * st2 + s1 * st1 * ct2);
	*pz = pr * (c1 * c2 - s1 * s2 * ct1);
	return 0;
}

double xpp_(float *srt)
{
	static float earray[14] = {20.f, 30.f, 40.f, 60.f, 80.f, 100.f, 170.f, 250.f,
							  310.f, 350.f, 460.f, 560.f, 660.f, 800.f};
	static float xarray[14] = {150.f, 90.f, 80.6f, 48.f, 36.6f, 31.6f, 25.9f, 24.f,
							  23.1f, 24.f, 28.3f, 33.6f, 41.5f, 47.f};

	float ret_val, r__1;

	double log(double), exp(double);

	static int ie;
	static float ekin, xmin, ymin, xmax, ymax, pmass;

	pmass = .9383f;
	r__1 = *srt / (pmass * 2.f);
	ekin = pmass * 2e3f * (r__1 * r__1 - 1.f);
	if (ekin < earray[0])
	{
		ret_val = xarray[0];
		if (ret_val > 55.f)
		{
			ret_val = 55.f;
		}
		return ret_val;
	}
	if (ekin > earray[13])
	{
		ret_val = xarray[13];
		return ret_val;
	}

	for (ie = 1; ie <= 14; ++ie)
	{
		if (earray[ie - 1] == ekin)
		{
			ret_val = xarray[ie - 1];
			if (ret_val > 55.f)
			{
				ret_val = 55.f;
			}
			return ret_val;
		}
		if (earray[ie - 1] > ekin)
		{
			ymin = log(xarray[ie - 2]);
			ymax = log(xarray[ie - 1]);
			xmin = log(earray[ie - 2]);
			xmax = log(earray[ie - 1]);
			ret_val = exp(ymin + (log(ekin) - xmin) * (ymax - ymin) / (xmax - xmin));
			if (ret_val > 55.f)
			{
				ret_val = 55.f;
			}
			goto L50;
		}
	}
L50:
	return ret_val;
}

double xnp_(float *srt)
{
	static float earray[11] = {20.f, 30.f, 40.f, 60.f, 90.f, 135.f, 200.f, 300.f,
							  400.f, 600.f, 800.f};
	static float xarray[11] = {410.f, 270.f, 214.5f, 130.f, 78.f, 53.5f, 41.6f,
							  35.9f, 34.2f, 34.3f, 34.9f};

	float ret_val, r__1;

	double log(double), exp(double);

	static int ie;
	static float ekin, xmin, ymin, xmax, ymax, pmass;

	pmass = .9383f;
	r__1 = *srt / (pmass * 2.f);
	ekin = pmass * 2e3f * (r__1 * r__1 - 1.f);
	if (ekin < earray[0])
	{
		ret_val = xarray[0];
		if (ret_val > 55.f)
		{
			ret_val = 55.f;
		}
		return ret_val;
	}
	if (ekin > earray[10])
	{
		ret_val = xarray[10];
		return ret_val;
	}

	for (ie = 1; ie <= 11; ++ie)
	{
		if (earray[ie - 1] == ekin)
		{
			ret_val = xarray[ie - 1];
			if (ret_val > 55.f)
			{
				ret_val = 55.f;
			}
			return ret_val;
		}
		if (earray[ie - 1] > ekin)
		{
			ymin = log(xarray[ie - 2]);
			ymax = log(xarray[ie - 1]);
			xmin = log(earray[ie - 2]);
			xmax = log(earray[ie - 1]);
			ret_val = exp(ymin + (log(ekin) - xmin) * (ymax - ymin) / (xmax - xmin));
			if (ret_val > 55.f)
			{
				ret_val = 55.f;
			}
			goto L50;
		}
	}
L50:
	return ret_val;
}

double ptr_(float *ptmax, int *iseed)
{
	float ret_val;

	double log(double), exp(double);

	static int ie;
	static float xt, xmin, ymin, xmax, ymax;
	extern double ptdis_(float *), ranart_(int *);
	static float tryial;

	ret_val = 0.f;
	if (*ptmax <= .01f)
	{
		ret_val = *ptmax;
		return ret_val;
	}
	if (*ptmax > 2.01f)
	{
		*ptmax = 2.01f;
	}
	tryial = ptdis_(ptmax) / ptdis_(&c_b561);
	xt = ranart_(&rndf77_1.nseed) * tryial;
	for (ie = 1; ie <= 200; ++ie)
	{
		if (table_1.earray[ie] == xt)
		{
			ret_val = table_1.xarray[ie];
			return ret_val;
		}
		if (table_1.xarray[ie - 1] <= 1e-5f)
		{
			goto L50;
		}
		if (table_1.xarray[ie] <= 1e-5f)
		{
			goto L50;
		}
		if (table_1.earray[ie - 1] <= 1e-5f)
		{
			goto L50;
		}
		if (table_1.earray[ie] <= 1e-5f)
		{
			goto L50;
		}
		if (table_1.earray[ie] > xt)
		{
			ymin = log(table_1.xarray[ie - 1]);
			ymax = log(table_1.xarray[ie]);
			xmin = log(table_1.earray[ie - 1]);
			xmax = log(table_1.earray[ie]);
			ret_val = exp(ymin + (log(xt) - xmin) * (ymax - ymin) / (xmax - xmin));
			if (ret_val > *ptmax)
			{
				ret_val = *ptmax;
			}
			return ret_val;
		}
	L50:;
	}
	return ret_val;
}

int xnd_(float *px, float *py, float *pz, float *srt, int *i1, int *i2, float *xinel, float *sigk, float *xsk1, float *xsk2,
		 float *xsk3, float *xsk4, float *xsk5)
{
	int i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8;
	float r__1, r__2, r__3, r__4, r__5;

	double sqrt(double);

	static float al, as, es, pr, ak0, em1, em2, ada, ana;
	extern int m1535_(int *, int *, float *, float *);
	static float akp, x1440, x1535, prf;
	extern double ppk0_(float *), ppk1_(float *);
	static float t1dlk, t2dlk, t1dsk, t2dsk, t1nlk, t1nsk, t2nsk;
	extern double sigma_(float *, int *, int *, int *), denom_(
																		   float *, float *);
	static float signd, pmdlk, sigdn, pmdsk, renom;
	extern double pplpk_(float *);
	static float pmnsk, pmdlk2, pmdsk2, renom1, pmnsk2, deltam, renomn;

	*xinel = 0.f;
	*sigk = 0.f;
	*xsk1 = 0.f;
	*xsk2 = 0.f;
	*xsk3 = 0.f;
	*xsk4 = 0.f;
	*xsk5 = 0.f;
	em1 = cc_1.e[*i1 - 1];
	em2 = cc_1.e[*i2 - 1];
	r__1 = *px;
	r__2 = *py;
	r__3 = *pz;
	pr = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
	if (*srt < 2.04f)
	{
		return 0;
	}
	r__1 = *srt;
	prf = sqrt(r__1 * r__1 * .25f - .88040689000000005f);
	if (em1 > 1.f)
	{
		deltam = em1;
	}
	else
	{
		deltam = em2;
	}
	r__1 = prf;
	renom = deltam * (r__1 * r__1) / denom_(srt, &c_b173) / pr;
	r__1 = prf;
	renomn = deltam * (r__1 * r__1) / denom_(srt, &c_b234) / pr;
	r__1 = prf;
	renom1 = deltam * (r__1 * r__1) / denom_(srt, &c_b235) / pr;
	if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2 && (i__2 = ee_1.lb[*i2 - 1],
													  abs(i__2)) == 6)
	{
		renom = 0.f;
	}
	if ((i__1 = ee_1.lb[*i2 - 1], abs(i__1)) == 2 && (i__2 = ee_1.lb[*i1 - 1],
													  abs(i__2)) == 6)
	{
		renom = 0.f;
	}
	if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 1 && (i__2 = ee_1.lb[*i2 - 1],
													  abs(i__2)) == 9)
	{
		renom = 0.f;
	}
	if ((i__1 = ee_1.lb[*i2 - 1], abs(i__1)) == 1 && (i__2 = ee_1.lb[*i1 - 1],
													  abs(i__2)) == 9)
	{
		renom = 0.f;
	}
	i__3 = (i__1 = ee_1.lb[*i1 - 1], abs(i__1));
	i__4 = (i__2 = ee_1.lb[*i2 - 1], abs(i__2));
	m1535_(&i__3, &i__4, srt, &x1535);
	x1440 = sigma_(srt, &c__2, &c__0, &c__1) * .75f;
	akp = .498f;
	ak0 = .498f;
	ana = .94f;
	ada = 1.232f;
	al = 1.1157f;
	as = 1.1197f;
	*xsk1 = 0.f;
	*xsk2 = 0.f;
	*xsk3 = 0.f;
	*xsk4 = 0.f;
	*xsk5 = 0.f;
	t1nlk = ana + al + akp;
	if (*srt <= t1nlk)
	{
		goto L222;
	}
	*xsk1 = pplpk_(srt) * 1.5f;
	t1dlk = ada + al + akp;
	t2dlk = ada + al - akp;
	if (*srt <= t1dlk)
	{
		goto L222;
	}
	es = *srt;
	r__1 = es;
	r__2 = t1dlk;
	r__3 = es;
	r__4 = t2dlk;
	r__5 = es;
	pmdlk2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / (r__5 * r__5 * 4.f);
	pmdlk = sqrt(pmdlk2);
	*xsk3 = pplpk_(srt) * 1.5f;
	t1nsk = ana + as + akp;
	t2nsk = ana + as - akp;
	if (*srt <= t1nsk)
	{
		goto L222;
	}
	r__1 = es;
	r__2 = t1nsk;
	r__3 = es;
	r__4 = t2nsk;
	r__5 = es;
	pmnsk2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / (r__5 * r__5 * 4.f);
	pmnsk = sqrt(pmnsk2);
	*xsk2 = (ppk1_(srt) + ppk0_(srt)) * 1.5f;
	t1dsk = ada + as + akp;
	t2dsk = ada + as - akp;
	if (*srt <= t1dsk)
	{
		goto L222;
	}
	r__1 = es;
	r__2 = t1dsk;
	r__3 = es;
	r__4 = t2dsk;
	r__5 = es;
	pmdsk2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / (r__5 * r__5 * 4.f);
	pmdsk = sqrt(pmdsk2);
	*xsk4 = (ppk1_(srt) + ppk0_(srt)) * 1.5f;
	if (*srt <= 2.898914f)
	{
		goto L222;
	}
	*xsk5 = 1e-4f;
L222:
	*sigk = *xsk1 + *xsk2 + *xsk3 + *xsk4;
	*xsk1 *= 2.f;
	*xsk2 *= 2.f;
	*xsk3 *= 2.f;
	*xsk4 *= 2.f;
	*sigk = *sigk * 2.f + *xsk5;
	if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2 && (i__2 = ee_1.lb[*i2 - 1],
													  abs(i__2)) == 6 ||
		(i__3 = ee_1.lb[*i2 - 1], abs(i__3)) == 2 && (i__4 = ee_1.lb[*i1 - 1], abs(i__4)) == 6 || (i__5 = ee_1.lb[*i1 - 1], abs(i__5)) == 1 && (i__6 = ee_1.lb[*i2 - 1], abs(i__6)) == 9 || (i__7 = ee_1.lb[*i2 - 1], abs(i__7)) == 1 && (i__8 = ee_1.lb[*i1 - 1], abs(i__8)) == 9)
	{
		*xinel = *sigk;
		return 0;
	}
	if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 18 && ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2 || (i__2 = ee_1.lb[*i2 - 1], abs(i__2)) == 2))
	{
		signd = sigma_(srt, &c__1, &c__1, &c__0) + sigma_(srt, &c__1, &c__1, &c__1) * .5f;
		sigdn = signd * .25f * renom;
		*xinel = sigdn + x1440 + x1535 + *sigk;
		return 0;
	}
	if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 6 && ((i__1 = ee_1.lb[*i1 - 1],
													  abs(i__1)) == 1 ||
													 (i__2 = ee_1.lb[*i2 - 1], abs(i__2)) == 1))
	{
		signd = sigma_(srt, &c__1, &c__1, &c__0) + sigma_(srt, &c__1, &c__1, &c__1) * .5f;
		sigdn = signd * .25f * renom;
		*xinel = sigdn + x1440 + x1535 + *sigk;
		return 0;
	}
	if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 8 && ((i__1 = ee_1.lb[*i1 - 1],
													  abs(i__1)) == 1 ||
													 (i__2 = ee_1.lb[*i2 - 1], abs(i__2)) == 1))
	{
		signd = sigma_(srt, &c__1, &c__1, &c__1) * 1.5f;
		sigdn = signd * .25f * renom;
		*xinel = sigdn + x1440 + x1535 + *sigk;
		return 0;
	}
	if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 14 && ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2 && (i__2 = ee_1.lb[*i2 - 1], abs(i__2)) == 2))
	{
		signd = sigma_(srt, &c__1, &c__1, &c__1) * 1.5f;
		sigdn = signd * .25f * renom;
		*xinel = sigdn + x1440 + x1535 + *sigk;
		return 0;
	}
	if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 16 && ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2 || (i__2 = ee_1.lb[*i2 - 1], abs(i__2)) == 2))
	{
		signd = sigma_(srt, &c__1, &c__1, &c__1) * .5f + sigma_(srt, &c__1, &c__1, &c__0) * .25f;
		sigdn = signd * .5f * renom;
		*xinel = sigdn + x1440 * 2.f + x1535 * 2.f + *sigk;
		return 0;
	}
	if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 7)
	{
		signd = sigma_(srt, &c__1, &c__1, &c__1) * .5f + sigma_(srt, &c__1, &c__1, &c__0) * .25f;
		sigdn = signd * .5f * renom;
		*xinel = sigdn + x1440 * 2.f + x1535 * 2.f + *sigk;
		return 0;
	}
	if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 10 && ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 1 || (i__2 = ee_1.lb[*i2 - 1], abs(i__2)) == 1))
	{
		signd = sigma_(srt, &c__2, &c__0, &c__1) * .75f;
		sigdn = signd * renomn;
		*xinel = sigdn + x1535 + *sigk;
		return 0;
	}
	if (ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1] == 22 && ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 2 || (i__2 = ee_1.lb[*i2 - 1], abs(i__2)) == 2))
	{
		signd = sigma_(srt, &c__2, &c__0, &c__1) * .75f;
		sigdn = signd * renomn;
		*xinel = sigdn + x1535 + *sigk;
		return 0;
	}
	if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) == 12 || (i__2 = ee_1.lb[*i1 - 1], abs(i__2)) == 13 || (i__3 = ee_1.lb[*i2 - 1], abs(i__3)) == 12 || (i__4 = ee_1.lb[*i2 - 1], abs(i__4)) == 13)
	{
		signd = x1535;
		sigdn = signd * renom1;
		*xinel = sigdn + *sigk;
		return 0;
	}
	return 0;
}

int xddin_(float *px, float *py, float *pz, float *srt, int *i1, int *i2, float *xinel, float *sigk, float *xsk1, float *xsk2,
		   float *xsk3, float *xsk4, float *xsk5)
{
	int i__1, i__2, i__3, i__4, i__5, i__6;
	float r__1, r__2, r__3, r__4, r__5;

	double sqrt(double);

	static float al, as, es, pr, ak0, em1, em2, s2d, ada, ana;
	static int idd;
	extern int n1535_(int *, int *, float *, float *);
	static float akp, x1535, sig2;
	extern double ppk0_(float *), ppk1_(float *);
	static float t1dlk, t2dlk, t1dsk, t2dsk, t1nlk, t1nsk, t2nsk;
	extern double sigma_(float *, int *, int *, int *);
	static float pmdlk, pmdsk;
	extern double pplpk_(float *);
	static float pmnsk;
	extern double reab2d_(int *, int *, float *);
	static float pmdlk2, pmdsk2, pmnsk2;

	*xinel = 0.f;
	*sigk = 0.f;
	*xsk1 = 0.f;
	*xsk2 = 0.f;
	*xsk3 = 0.f;
	*xsk4 = 0.f;
	*xsk5 = 0.f;
	em1 = cc_1.e[*i1 - 1];
	em2 = cc_1.e[*i2 - 1];
	r__1 = *px;
	r__2 = *py;
	r__3 = *pz;
	pr = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
	i__3 = (i__1 = ee_1.lb[*i1 - 1], abs(i__1));
	i__4 = (i__2 = ee_1.lb[*i2 - 1], abs(i__2));
	n1535_(&i__3, &i__4, srt, &x1535);

	akp = .498f;
	ak0 = .498f;
	ana = .94f;
	ada = 1.232f;
	al = 1.1157f;
	as = 1.1197f;
	*xsk1 = 0.f;
	*xsk2 = 0.f;
	*xsk3 = 0.f;
	*xsk4 = 0.f;
	t1nlk = ana + al + akp;
	if (*srt <= t1nlk)
	{
		goto L222;
	}
	*xsk1 = pplpk_(srt) * 1.5f;
	t1dlk = ada + al + akp;
	t2dlk = ada + al - akp;
	if (*srt <= t1dlk)
	{
		goto L222;
	}
	es = *srt;
	r__1 = es;
	r__2 = t1dlk;
	r__3 = es;
	r__4 = t2dlk;
	r__5 = es;
	pmdlk2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / (r__5 * r__5 * 4.f);
	pmdlk = sqrt(pmdlk2);
	*xsk3 = pplpk_(srt) * 1.5f;
	t1nsk = ana + as + akp;
	t2nsk = ana + as - akp;
	if (*srt <= t1nsk)
	{
		goto L222;
	}
	r__1 = es;
	r__2 = t1nsk;
	r__3 = es;
	r__4 = t2nsk;
	r__5 = es;
	pmnsk2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / (r__5 * r__5 * 4.f);
	pmnsk = sqrt(pmnsk2);
	*xsk2 = (ppk1_(srt) + ppk0_(srt)) * 1.5f;
	t1dsk = ada + as + akp;
	t2dsk = ada + as - akp;
	if (*srt <= t1dsk)
	{
		goto L222;
	}
	r__1 = es;
	r__2 = t1dsk;
	r__3 = es;
	r__4 = t2dsk;
	r__5 = es;
	pmdsk2 = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / (r__5 * r__5 * 4.f);
	pmdsk = sqrt(pmdsk2);
	*xsk4 = (ppk1_(srt) + ppk0_(srt)) * 1.5f;
	if (*srt <= 2.898914f)
	{
		goto L222;
	}
	*xsk5 = 1e-4f;
L222:
	*sigk = *xsk1 + *xsk2 + *xsk3 + *xsk4;
	*xsk1 *= 2.f;
	*xsk2 *= 2.f;
	*xsk3 *= 2.f;
	*xsk4 *= 2.f;
	*sigk = *sigk * 2.f + *xsk5;
	idd = (i__1 = ee_1.lb[*i1 - 1] * ee_1.lb[*i2 - 1], abs(i__1));
	s2d = reab2d_(i1, i2, srt);
	s2d = 0.f;
	if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) >= 12 && (i__2 = ee_1.lb[*i2 - 1], abs(i__2)) >= 12 || (i__3 = ee_1.lb[*i1 - 1], abs(i__3)) >= 12 && (i__4 = ee_1.lb[*i2 - 1], abs(i__4)) >= 6 || (i__5 = ee_1.lb[*i2 - 1], abs(i__5)) >= 12 && (i__6 = ee_1.lb[*i1 - 1], abs(i__6)) >= 6)
	{
		*xinel = *sigk + s2d;
		return 0;
	}
	if (idd == 63 || idd == 64 || idd == 48 || idd == 49 || idd == 121 || idd == 100 || idd == 88 || idd == 66 || idd == 90 || idd == 70)
	{
		*xinel = x1535 + *sigk + s2d;
		return 0;
	}
	if (idd == 110 || idd == 77 || idd == 80)
	{
		*xinel = x1535 + *sigk + s2d;
		return 0;
	}
	if (idd == 54 || idd == 56)
	{
		sig2 = sigma_(srt, &c__2, &c__0, &c__1) * .75f;
		*xinel = (sig2 + x1535) * 2.f + *sigk + s2d;
		return 0;
	}
	return 0;
}

double dirct1_(float *srt)
{
	static float earray[122] = {1.5683f, 1.5783f, 1.5883f, 1.5983f, 1.6083f,
							   1.6183f, 1.6283f, 1.6383f, 1.6483f, 1.6583f, 1.6683f, 1.6783f, 1.6883f,
							   1.6983f, 1.7083f, 1.7183f, 1.7283f, 1.7383f, 1.7483f, 1.7583f, 1.7683f,
							   1.7783f, 1.7883f, 1.7983f, 1.8083f, 1.8183f, 1.8283f, 1.8383f, 1.8483f,
							   1.8583f, 1.8683f, 1.8783f, 1.8883f, 1.8983f, 1.9083f, 1.9183f, 1.9283f,
							   1.9383f, 1.9483f, 1.9583f, 1.9683f, 1.9783f, 1.9883f, 1.9983f, 2.0083f,
							   2.0183f, 2.0283f, 2.0383f, 2.0483f, 2.0583f, 2.0683f, 2.0783f, 2.0883f,
							   2.0983f, 2.1083f, 2.1183f, 2.1283f, 2.1383f, 2.1483f, 2.1583f, 2.1683f,
							   2.1783f, 2.1883f, 2.1983f, 2.2083f, 2.2183f, 2.2283f, 2.2383f, 2.2483f,
							   2.2583f, 2.2683f, 2.2783f, 2.2883f, 2.2983f, 2.3083f, 2.3183f, 2.3283f,
							   2.3383f, 2.3483f, 2.3583f, 2.3683f, 2.3783f, 2.3883f, 2.3983f, 2.4083f,
							   2.4183f, 2.4283f, 2.4383f, 2.4483f, 2.4583f, 2.4683f, 2.4783f, 2.4883f,
							   2.4983f, 2.5083f, 2.5183f, 2.5283f, 2.5383f, 2.5483f, 2.5583f, 2.5683f,
							   2.5783f, 2.5883f, 2.5983f, 2.6083f, 2.6183f, 2.6283f, 2.6383f, 2.6483f,
							   2.6583f, 2.6683f, 2.6783f, 2.6883f, 2.6983f, 2.7083f, 2.7183f, 2.7283f,
							   2.7383f, 2.7483f, 2.7583f, 2.7683f, 2.7783f};
	static float xarray[122] = {.017764091f, .5643668f, .8150568f, 1.045565f,
							   2.133695f, 3.327922f, 4.206488f, 3.471242f, 4.486876f, 5.542213f,
							   6.800052f, 7.192446f, 6.829848f, 6.580306f, 6.86841f, 8.527946f,
							   10.1572f, 9.716511f, 9.298335f, 8.90131f, 10.31213f, 10.52185f,
							   11.1763f, 11.61639f, 12.05577f, 12.71596f, 13.46036f, 14.2206f,
							   14.65449f, 14.94775f, 14.9331f, 15.32907f, 16.56481f, 16.29422f,
							   15.18548f, 14.12658f, 13.72544f, 13.24488f, 13.31003f, 14.4268f,
							   12.84423f, 12.49025f, 12.14858f, 11.8187f, 11.18993f, 11.35816f,
							   11.09447f, 10.83873f, 10.61592f, 10.53754f, 9.425521f, 8.195912f,
							   9.661075f, 9.696192f, 9.200142f, 8.953734f, 8.715461f, 8.484999f,
							   8.320765f, 8.255512f, 8.190969f, 8.127125f, 8.079508f, 8.073004f,
							   8.010611f, 7.948909f, 7.887895f, 7.761005f, 7.62629f, 7.494696f,
							   7.366132f, 7.530178f, 8.392097f, 9.046881f, 8.962544f, 8.879403f,
							   8.797427f, 8.716601f, 8.636904f, 8.558312f, 8.404368f, 8.328978f,
							   8.254617f, 8.181265f, 8.108907f, 8.037527f, 7.9671f, 7.897617f,
							   7.829057f, 7.761405f, 7.694647f, 7.628764f, 7.563742f, 7.49957f,
							   7.387562f, 7.273281f, 7.161334f, 6.973375f, 6.529592f, 6.280323f,
							   6.293136f, 6.305725f, 6.318097f, 6.330258f, 6.342214f, 6.353968f,
							   6.365528f, 6.376895f, 6.388079f, 6.399081f, 6.409906f, 6.42056f,
							   6.431045f, 6.441367f, 6.451529f, 6.461533f, 6.471386f, 6.481091f,
							   6.49065f, 6.476413f, 6.297259f, 6.097826f};

	float ret_val;

	double log(double), exp(double);

	static int ie;
	static float xmin, ymin, xmax, ymax;

	if (*srt < earray[0])
	{
		ret_val = 1e-5f;
		return ret_val;
	}
	if (*srt > earray[121])
	{
		ret_val = xarray[121];
		ret_val /= 10.f;
		return ret_val;
	}

	for (ie = 1; ie <= 122; ++ie)
	{
		if (earray[ie - 1] == *srt)
		{
			ret_val = xarray[ie - 1];
			ret_val /= 10.f;
			return ret_val;
		}
		if (earray[ie - 1] > *srt)
		{
			ymin = log(xarray[ie - 2]);
			ymax = log(xarray[ie - 1]);
			xmin = log(earray[ie - 2]);
			xmax = log(earray[ie - 1]);
			ret_val = exp(ymin + (log(*srt) - xmin) * (ymax - ymin) / (xmax - xmin));
			ret_val /= 10.f;
			goto L50;
		}
	}
L50:
	return ret_val;
}

double dirct2_(float *srt)
{
	static float earray[122] = {1.5683f, 1.5783f, 1.5883f, 1.5983f, 1.6083f,
							   1.6183f, 1.6283f, 1.6383f, 1.6483f, 1.6583f, 1.6683f, 1.6783f, 1.6883f,
							   1.6983f, 1.7083f, 1.7183f, 1.7283f, 1.7383f, 1.7483f, 1.7583f, 1.7683f,
							   1.7783f, 1.7883f, 1.7983f, 1.8083f, 1.8183f, 1.8283f, 1.8383f, 1.8483f,
							   1.8583f, 1.8683f, 1.8783f, 1.8883f, 1.8983f, 1.9083f, 1.9183f, 1.9283f,
							   1.9383f, 1.9483f, 1.9583f, 1.9683f, 1.9783f, 1.9883f, 1.9983f, 2.0083f,
							   2.0183f, 2.0283f, 2.0383f, 2.0483f, 2.0583f, 2.0683f, 2.0783f, 2.0883f,
							   2.0983f, 2.1083f, 2.1183f, 2.1283f, 2.1383f, 2.1483f, 2.1583f, 2.1683f,
							   2.1783f, 2.1883f, 2.1983f, 2.2083f, 2.2183f, 2.2283f, 2.2383f, 2.2483f,
							   2.2583f, 2.2683f, 2.2783f, 2.2883f, 2.2983f, 2.3083f, 2.3183f, 2.3283f,
							   2.3383f, 2.3483f, 2.3583f, 2.3683f, 2.3783f, 2.3883f, 2.3983f, 2.4083f,
							   2.4183f, 2.4283f, 2.4383f, 2.4483f, 2.4583f, 2.4683f, 2.4783f, 2.4883f,
							   2.4983f, 2.5083f, 2.5183f, 2.5283f, 2.5383f, 2.5483f, 2.5583f, 2.5683f,
							   2.5783f, 2.5883f, 2.5983f, 2.6083f, 2.6183f, 2.6283f, 2.6383f, 2.6483f,
							   2.6583f, 2.6683f, 2.6783f, 2.6883f, 2.6983f, 2.7083f, 2.7183f, 2.7283f,
							   2.7383f, 2.7483f, 2.7583f, 2.7683f, 2.7783f};
	static float xarray[122] = {.5773182f, 1.404156f, 2.578629f, 3.832013f,
							   4.906011f, 9.076963f, 13.10492f, 10.65975f, 15.31156f, 19.77611f,
							   19.92874f, 18.68979f, 19.80114f, 18.39536f, 14.34269f, 13.35353f,
							   13.58822f, 14.57031f, 10.24686f, 11.23386f, 9.764803f, 10.35652f,
							   10.53539f, 10.07524f, 9.582198f, 9.596469f, 9.818489f, 9.012848f,
							   9.378012f, 9.529244f, 9.529698f, 8.835624f, 6.671396f, 8.797758f,
							   8.133437f, 7.866227f, 7.823946f, 7.808504f, 7.791755f, 7.502062f,
							   7.417275f, 7.592349f, 7.752028f, 7.910585f, 8.068122f, 8.224736f,
							   8.075289f, 7.895902f, 7.721359f, 7.551512f, 7.386224f, 7.225343f,
							   7.068739f, 6.916284f, 6.767842f, 6.623294f, 6.48252f, 6.345404f,
							   6.211833f, 7.33951f, 7.531462f, 7.724824f, 7.91962f, 7.848021f,
							   7.639856f, 7.571083f, 7.508881f, 7.447474f, 7.386855f, 7.327011f,
							   7.164454f, 7.001266f, 6.842526f, 6.688094f, 6.537823f, 6.391583f,
							   6.249249f, 6.110689f, 5.97579f, 5.8942f, 5.959503f, 6.024602f,
							   6.089505f, 6.154224f, 6.21876f, 6.283128f, 6.347331f, 6.297411f,
							   6.120248f, 5.948606f, 6.494864f, 6.357106f, 6.222824f, 6.09191f,
							   5.964267f, 5.839795f, 5.718402f, 5.599994f, 5.499146f, 5.451325f,
							   5.404156f, 5.357625f, 5.311721f, 5.266435f, 5.301964f, 5.343963f,
							   5.385833f, 5.427577f, 5.4692f, 5.510702f, 5.552088f, 5.593359f,
							   5.63452f, 5.67557f, 5.716515f, 5.757356f, 5.798093f, 5.838732f,
							   5.879272f, 5.919717f, 5.960068f, 5.980941f};

	float ret_val;

	double log(double), exp(double);

	static int ie;
	static float xmin, ymin, xmax, ymax;

	if (*srt < earray[0])
	{
		ret_val = 1e-5f;
		return ret_val;
	}
	if (*srt > earray[121])
	{
		ret_val = xarray[121];
		ret_val /= 10.f;
		return ret_val;
	}

	for (ie = 1; ie <= 122; ++ie)
	{
		if (earray[ie - 1] == *srt)
		{
			ret_val = xarray[ie - 1];
			ret_val /= 10.f;
			return ret_val;
		}
		if (earray[ie - 1] > *srt)
		{
			ymin = log(xarray[ie - 2]);
			ymax = log(xarray[ie - 1]);
			xmin = log(earray[ie - 2]);
			xmax = log(earray[ie - 1]);
			ret_val = exp(ymin + (log(*srt) - xmin) * (ymax - ymin) / (xmax - xmin));
			ret_val /= 10.f;
			goto L50;
		}
	}
L50:
	return ret_val;
}

double erhon_(float *em1, float *em2, int *lb1, int *lb2, float *srt)
{
	static float arrayj[19] = {.5f, 1.5f, .5f, .5f, 2.5f, 2.5f, 1.5f, .5f, 1.5f, 3.5f,
							  1.5f, .5f, 1.5f, .5f, 2.5f, .5f, 1.5f, 2.5f, 3.5f};
	static float arrayl[19] = {1.f, 2.f, 0.f, 0.f, 2.f, 3.f, 2.f, 1.f, 1.f, 3.f, 1.f,
							  0.f, 2.f, 0.f, 3.f, 1.f, 1.f, 2.f, 3.f};
	static float arraym[19] = {1.44f, 1.52f, 1.535f, 1.65f, 1.675f, 1.68f, 1.7f,
							  1.71f, 1.72f, 1.99f, 1.6f, 1.62f, 1.7f, 1.9f, 1.905f, 1.91f, 1.86f, 1.93f,
							  1.95f};
	static float arrayw[19] = {.2f, .125f, .15f, .15f, .155f, .125f, .1f, .11f, .2f,
							  .29f, .25f, .16f, .28f, .15f, .3f, .22f, .25f, .25f, .24f};
	static float arrayb[19] = {.15f, .2f, .05f, .175f, .025f, .125f, .1f, .2f, .53f,
							  .34f, .05f, .07f, .15f, .45f, .45f, .058f, .08f, .12f, .08f};

	int i__1, i__2;
	float ret_val;

	static float pi;
	static int ir;
	static float xs, xs0;
	extern double fdr_(float *, float *, float *, float *, float *, float *,
						   float *, float *);
	static float branch;

	pi = 3.1415926f;
	xs = 0.f;
	for (ir = 1; ir <= 19; ++ir)
	{
		if (ir <= 8)
		{
			if (*lb1 * *lb2 == 27 && (*lb1 == 1 || *lb2 == 1) || *lb1 * *lb2 == 50 && (*lb1 == 2 || *lb2 == 2) || (*lb1 * *lb2 == -25 && (*lb1 == -1 || *lb2 == -1) || *lb1 * *lb2 == -54 && (*lb1 == -2 || *lb2 == -2)))
			{
				branch = 0.f;
			}
			if ((i__1 = *lb1 * *lb2, abs(i__1)) == 26 && (abs(*lb1) == 1 ||
														  abs(*lb2) == 1) ||
				(i__2 = *lb1 * *lb2, abs(i__2)) == 52 && (abs(*lb1) == 2 || abs(*lb2) == 2))
			{
				branch = .33333333333333331f;
			}
			if (*lb1 * *lb2 == 54 && (*lb1 == 2 || *lb2 == 2) || *lb1 * *lb2 == 25 && (*lb1 == 1 || *lb2 == 1) || (*lb1 * *lb2 == -50 && (*lb1 == -2 || *lb2 == -2) || *lb1 * *lb2 == -27 && (*lb1 == -1 || *lb2 == -1)))
			{
				branch = .66666666666666663f;
			}
		}
		else
		{
			if (*lb1 * *lb2 == 27 && (*lb1 == 1 || *lb2 == 1) || *lb1 * *lb2 == 50 && (*lb1 == 2 || *lb2 == 2) || (*lb1 * *lb2 == -25 && (*lb1 == -1 || *lb2 == -1) || *lb1 * *lb2 == -54 && (*lb1 == -2 || *lb2 == -2)))
			{
				branch = 1.f;
			}
			if ((i__1 = *lb1 * *lb2, abs(i__1)) == 26 && (abs(*lb1) == 1 ||
														  abs(*lb2) == 1) ||
				(i__2 = *lb1 * *lb2, abs(i__2)) == 52 && (abs(*lb1) == 2 || abs(*lb2) == 2))
			{
				branch = .66666666666666663f;
			}
			if (*lb1 * *lb2 == 54 && (*lb1 == 2 || *lb2 == 2) || *lb1 * *lb2 == 25 && (*lb1 == 1 || *lb2 == 1) || (*lb1 * *lb2 == -50 && (*lb1 == -2 || *lb2 == -2) || *lb1 * *lb2 == -27 && (*lb1 == -1 || *lb2 == -1)))
			{
				branch = .33333333333333331f;
			}
		}
		xs0 = fdr_(&arraym[ir - 1], &arrayj[ir - 1], &arrayl[ir - 1], &arrayw[ir - 1], &arrayb[ir - 1], srt, em1, em2);
		xs += pi * 1.3f * branch * xs0 * .038927290000000003f;
	}
	ret_val = xs;
	return ret_val;
}

double fdr_(float *dmass, float *aj, float *al, float *width, float *widb0,
				float *srt, float *em1, float *em2)
{
	float ret_val, r__1, r__2, r__3, r__4, r__5;
	double d__1, d__2, d__3, d__4;

	double sqrt(double), pow_dd(double *, double *);

	static float b, q, q0, ak2, ak02, amd, amp;

	amd = *em1;
	amp = *em2;
	r__2 = *dmass;
	r__3 = amd;
	r__4 = amp;
	r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
	r__5 = amp * amd;
	ak02 = r__1 * r__1 * .25f - r__5 * r__5;
	if (ak02 > 0.f)
	{
		q0 = sqrt(ak02 / *dmass);
	}
	else
	{
		q0 = 0.f;
		ret_val = 0.f;
		return ret_val;
	}
	r__2 = *srt;
	r__3 = amd;
	r__4 = amp;
	r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
	r__5 = amp * amd;
	ak2 = r__1 * r__1 * .25f - r__5 * r__5;
	if (ak2 > 0.f)
	{
		q = sqrt(ak2 / *dmass);
	}
	else
	{
		q = 0.f;
		ret_val = 0.f;
		return ret_val;
	}
	d__1 = (double)(q / q0);
	d__2 = (double)(*al * 2.f + 1);
	d__3 = (double)(q / q0);
	d__4 = (double)(*al * 2);
	b = *widb0 * 1.2f * *dmass / *srt * pow_dd(&d__1, &d__2) / (pow_dd(&d__3, &d__4) * .2f + 1.f);
	r__1 = *width;
	r__2 = *srt - *dmass;
	r__3 = *width;
	r__4 = q;
	ret_val = (*aj * 2.f + 1) * (r__1 * r__1) * b / (r__2 * r__2 + r__3 * r__3 * .25f) / (r__4 * r__4 * 6.f);
	return ret_val;
}

double dirct3_(float *srt)
{
	static float arrayj[17] = {1.5f, .5f, 2.5f, 2.5f, 1.5f, .5f, 1.5f, 3.5f, 1.5f, .5f,
							  1.5f, .5f, 2.5f, .5f, 1.5f, 2.5f, 3.5f};
	static float arrayl[17] = {2.f, 0.f, 2.f, 3.f, 2.f, 1.f, 1.f, 3.f, 1.f, 0.f, 2.f,
							  0.f, 3.f, 1.f, 1.f, 2.f, 3.f};
	static float arraym[17] = {1.52f, 1.65f, 1.675f, 1.68f, 1.7f, 1.71f, 1.72f,
							  1.99f, 1.6f, 1.62f, 1.7f, 1.9f, 1.905f, 1.91f, 1.86f, 1.93f, 1.95f};
	static float arrayw[17] = {.125f, .15f, .155f, .125f, .1f, .11f, .2f, .29f, .25f,
							  .16f, .28f, .15f, .3f, .22f, .25f, .25f, .24f};
	static float arrayb[17] = {.55f, .6f, .375f, .6f, .1f, .15f, .15f, .05f, .35f, .3f,
							  .15f, .1f, .1f, .22f, .2f, .09f, .4f};

	float ret_val;

	static float pi;
	static int ir;
	static float xs;
	extern double fd1_(float *, float *, float *, float *, float *, float *);
	static float xs0, amn, amp, branch;

	pi = 3.1415926f;
	amn = .938f;
	amp = .138f;
	xs = 0.f;
	branch = .33333333333333331f;
	for (ir = 1; ir <= 17; ++ir)
	{
		if (ir > 8)
		{
			branch = .66666666666666663f;
		}
		xs0 = fd1_(&arraym[ir - 1], &arrayj[ir - 1], &arrayl[ir - 1], &arrayw[ir - 1], &arrayb[ir - 1], srt);
		xs += pi * 1.3f * branch * xs0 * .038927290000000003f;
	}
	ret_val = xs;
	return ret_val;
}

double fd1_(float *dmass, float *aj, float *al, float *width, float *widb0,
				float *srt)
{
	float ret_val, r__1, r__2, r__3, r__4, r__5;
	double d__1, d__2, d__3, d__4;

	double sqrt(double), pow_dd(double *, double *);

	static float b, q, q0, ak2, ak02, amd, amn, amp;

	amn = .938f;
	amp = .138f;
	amd = amn;
	r__2 = *dmass;
	r__3 = amd;
	r__4 = amp;
	r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
	r__5 = amp * amd;
	ak02 = r__1 * r__1 * .25f - r__5 * r__5;
	if (ak02 > 0.f)
	{
		q0 = sqrt(ak02 / *dmass);
	}
	else
	{
		q0 = 0.f;
		ret_val = 0.f;
		return ret_val;
	}
	r__2 = *srt;
	r__3 = amd;
	r__4 = amp;
	r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
	r__5 = amp * amd;
	ak2 = r__1 * r__1 * .25f - r__5 * r__5;
	if (ak2 > 0.f)
	{
		q = sqrt(ak2 / *dmass);
	}
	else
	{
		q = 0.f;
		ret_val = 0.f;
		return ret_val;
	}
	d__1 = (double)(q / q0);
	d__2 = (double)(*al * 2.f + 1);
	d__3 = (double)(q / q0);
	d__4 = (double)(*al * 2);
	b = *widb0 * 1.2f * *dmass / *srt * pow_dd(&d__1, &d__2) / (pow_dd(&d__3, &d__4) * .2f + 1.f);
	r__1 = *width;
	r__2 = *srt - *dmass;
	r__3 = *width;
	r__4 = q;
	ret_val = (*aj * 2.f + 1) * (r__1 * r__1) * b / (r__2 * r__2 + r__3 * r__3 * .25f) / (r__4 * r__4 * 2.f);
	return ret_val;
}

double dpion_(float *em1, float *em2, int *lb1, int *lb2, float *srt)
{
	static float arrayj[19] = {.5f, 1.5f, .5f, .5f, 2.5f, 2.5f, 1.5f, .5f, 1.5f, 3.5f,
							  1.5f, .5f, 1.5f, .5f, 2.5f, .5f, 1.5f, 2.5f, 3.5f};
	static float arrayl[19] = {1.f, 2.f, 0.f, 0.f, 2.f, 3.f, 2.f, 1.f, 1.f, 3.f, 1.f,
							  0.f, 2.f, 0.f, 3.f, 1.f, 1.f, 2.f, 3.f};
	static float arraym[19] = {1.44f, 1.52f, 1.535f, 1.65f, 1.675f, 1.68f, 1.7f,
							  1.71f, 1.72f, 1.99f, 1.6f, 1.62f, 1.7f, 1.9f, 1.905f, 1.91f, 1.86f, 1.93f,
							  1.95f};
	static float arrayw[19] = {.2f, .125f, .15f, .15f, .155f, .125f, .1f, .11f, .2f,
							  .29f, .25f, .16f, .28f, .15f, .3f, .22f, .25f, .25f, .24f};
	static float arrayb[19] = {.15f, .25f, 0.f, .05f, .575f, .125f, .379f, .1f, .1f,
							  .062f, .45f, .6f, .6984f, .05f, .25f, .089f, .19f, .2f, .13f};

	int i__1, i__2;
	float ret_val;

	static float pi;
	static int ir;
	static float xs;
	extern double fd2_(float *, float *, float *, float *, float *, float *,
						   float *, float *);
	static float xs0, amn, amp, branch;

	pi = 3.1415926f;
	amn = .94f;
	amp = .14f;
	xs = 0.f;
	for (ir = 1; ir <= 19; ++ir)
	{
		branch = 0.f;
		if (ir <= 8)
		{
			if (*lb1 * *lb2 == 35 && (*lb1 == 5 || *lb2 == 5) || *lb1 * *lb2 == 24 && (*lb1 == 3 || *lb2 == 3) || (*lb1 * *lb2 == -21 && (*lb1 == 3 || *lb2 == 3) || *lb1 * *lb2 == -40 && (*lb1 == 5 || *lb2 == 5)))
			{
				branch = .16666666666666666f;
			}
			if ((i__1 = *lb1 * *lb2, abs(i__1)) == 28 && (*lb1 == 4 || *lb2 ==
																		   4) ||
				(i__2 = *lb1 * *lb2, abs(i__2)) == 32 && (*lb1 ==
															  4 ||
														  *lb2 == 4))
			{
				branch = .33333333333333331f;
			}
			if (*lb1 * *lb2 == 30 && (*lb1 == 5 || *lb2 == 5) || *lb1 * *lb2 == 27 && (*lb1 == 3 || *lb2 == 3) || (*lb1 * *lb2 == -18 && (*lb1 == 3 || *lb2 == 3) || *lb1 * *lb2 == -45 && (*lb1 == 5 || *lb2 == 5)))
			{
				branch = .5f;
			}
		}
		else
		{
			if (*lb1 * *lb2 == 40 && (*lb1 == 5 || *lb2 == 5) || *lb1 * *lb2 == 30 && (*lb1 == 5 || *lb2 == 5) || (*lb1 * *lb2 == -24 && (*lb1 == 3 || *lb2 == 3) || *lb1 * *lb2 == -18 && (*lb1 == 3 || *lb2 == 3)))
			{
				branch = .40000000000000002f;
			}
			if (*lb1 * *lb2 == 27 && (*lb1 == 3 || *lb2 == 3) || *lb1 * *lb2 == 21 && (*lb1 == 3 || *lb2 == 3) || (*lb1 * *lb2 == -45 && (*lb1 == 5 || *lb2 == 5) || *lb1 * *lb2 == -35 && (*lb1 == 5 || *lb2 == 5)))
			{
				branch = .40000000000000002f;
			}
			if (*lb1 * *lb2 == 35 && (*lb1 == 5 || *lb2 == 5) || *lb1 * *lb2 == 24 && (*lb1 == 3 || *lb2 == 3) || (*lb1 * *lb2 == -21 && (*lb1 == 3 || *lb2 == 3) || *lb1 * *lb2 == -40 && (*lb1 == 5 || *lb2 == 5)))
			{
				branch = .53333333333333333f;
			}
			if ((i__1 = *lb1 * *lb2, abs(i__1)) == 28 && (*lb1 == 4 || *lb2 ==
																		   4) ||
				(i__2 = *lb1 * *lb2, abs(i__2)) == 32 && (*lb1 ==
															  4 ||
														  *lb2 == 4))
			{
				branch = .066666666666666666f;
			}
			if ((i__1 = *lb1 * *lb2, abs(i__1)) == 36 && (*lb1 == 4 || *lb2 ==
																		   4) ||
				(i__2 = *lb1 * *lb2, abs(i__2)) == 24 && (*lb1 ==
															  4 ||
														  *lb2 == 4))
			{
				branch = .59999999999999998f;
			}
		}
		xs0 = fd2_(&arraym[ir - 1], &arrayj[ir - 1], &arrayl[ir - 1], &arrayw[ir - 1], &arrayb[ir - 1], em1, em2, srt);
		xs += pi * 1.3f * branch * xs0 * .038927290000000003f;
	}
	ret_val = xs;
	return ret_val;
}

double fd2_(float *dmass, float *aj, float *al, float *width, float *widb0,
				float *em1, float *em2, float *srt)
{
	float ret_val, r__1, r__2, r__3, r__4, r__5;
	double d__1, d__2, d__3, d__4;

	double sqrt(double), pow_dd(double *, double *);

	static float b, q, q0, ak2, ak02, amd, amp;

	amp = *em1;
	amd = *em2;
	r__2 = *dmass;
	r__3 = amd;
	r__4 = amp;
	r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
	r__5 = amp * amd;
	ak02 = r__1 * r__1 * .25f - r__5 * r__5;
	if (ak02 > 0.f)
	{
		q0 = sqrt(ak02 / *dmass);
	}
	else
	{
		q0 = 0.f;
		ret_val = 0.f;
		return ret_val;
	}
	r__2 = *srt;
	r__3 = amd;
	r__4 = amp;
	r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
	r__5 = amp * amd;
	ak2 = r__1 * r__1 * .25f - r__5 * r__5;
	if (ak2 > 0.f)
	{
		q = sqrt(ak2 / *dmass);
	}
	else
	{
		q = 0.f;
		ret_val = 0.f;
		return ret_val;
	}
	d__1 = (double)(q / q0);
	d__2 = (double)(*al * 2.f + 1);
	d__3 = (double)(q / q0);
	d__4 = (double)(*al * 2);
	b = *widb0 * 1.2f * *dmass / *srt * pow_dd(&d__1, &d__2) / (pow_dd(&d__3, &d__4) * .2f + 1.f);
	r__1 = *width;
	r__2 = *srt - *dmass;
	r__3 = *width;
	r__4 = q;
	ret_val = (*aj * 2.f + 1) * (r__1 * r__1) * b / (r__2 * r__2 + r__3 * r__3 * .25f) / (r__4 * r__4 * 4.f);
	return ret_val;
}

int rmasdd_(float *srt, float *am10, float *am20, float *dmin1,
			float *dmin2, int *iseed, int *ic, float *dm1, float *dm2)
{
	float r__1, r__2, r__3, r__4, r__5, r__6;

	static float q2, fm1, fm2, fff, amn, amp, prob;
	static int ntry;
	static float dmax1, dmax2, prob0;
	static int ntry1, ntry2, ictrl;
	extern double fmassd_(float *), ranart_(int *), fmassn_(float *),
		fmassr_(float *);

	amn = .94f;
	amp = .14f;
	dmax1 = *srt - *dmin2;
L5:
	ntry1 = 0;
	ntry2 = 0;
	ntry = 0;
	ictrl = 0;
L10:
	*dm1 = ranart_(&rndf77_1.nseed) * (dmax1 - *dmin1) + *dmin1;
	++ntry1;
	if (ictrl == 0)
	{
		dmax2 = *srt - *dm1;
	}
L20:
	*dm2 = ranart_(&rndf77_1.nseed) * (dmax2 - *dmin2) + *dmin2;
	++ntry2;
	r__2 = *srt;
	r__3 = *dm1;
	r__4 = *dm2;
	r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
	r__5 = *dm1;
	r__6 = *dm2;
	q2 = r__1 * r__1 - r__5 * r__5 * 4.f * (r__6 * r__6);
	if (q2 <= 0.f)
	{
		dmax2 = *dm2 - .01f;
		ictrl = 1;
		goto L20;
	}
	if (dmax1 < *am10)
	{
		if (*ic == 1)
		{
			fm1 = fmassd_(&dmax1);
		}
		if (*ic == 2)
		{
			fm1 = fmassn_(&dmax1);
		}
		if (*ic == 3)
		{
			fm1 = fmassd_(&dmax1);
		}
		if (*ic == 4)
		{
			fm1 = fmassd_(&dmax1);
		}
	}
	else
	{
		if (*ic == 1)
		{
			fm1 = fmassd_(am10);
		}
		if (*ic == 2)
		{
			fm1 = fmassn_(am10);
		}
		if (*ic == 3)
		{
			fm1 = fmassd_(am10);
		}
		if (*ic == 4)
		{
			fm1 = fmassd_(am10);
		}
	}
	if (dmax2 < *am20)
	{
		if (*ic == 1)
		{
			fm2 = fmassd_(&dmax2);
		}
		if (*ic == 2)
		{
			fm2 = fmassn_(&dmax2);
		}
		if (*ic == 3)
		{
			fm2 = fmassn_(&dmax2);
		}
		if (*ic == 4)
		{
			fm2 = fmassr_(&dmax2);
		}
	}
	else
	{
		if (*ic == 1)
		{
			fm2 = fmassd_(am20);
		}
		if (*ic == 2)
		{
			fm2 = fmassn_(am20);
		}
		if (*ic == 3)
		{
			fm2 = fmassn_(am20);
		}
		if (*ic == 4)
		{
			fm2 = fmassr_(am20);
		}
	}
	if (fm1 == 0.f)
	{
		fm1 = 1e-4f;
	}
	if (fm2 == 0.f)
	{
		fm2 = 1e-4f;
	}
	prob0 = fm1 * fm2;
	if (*ic == 1)
	{
		prob = fmassd_(dm1) * fmassd_(dm2);
	}
	if (*ic == 2)
	{
		prob = fmassn_(dm1) * fmassn_(dm2);
	}
	if (*ic == 3)
	{
		prob = fmassd_(dm1) * fmassn_(dm2);
	}
	if (*ic == 4)
	{
		prob = fmassd_(dm1) * fmassr_(dm2);
	}
	if (prob <= 1e-6f)
	{
		prob = 1e-6f;
	}
	fff = prob / prob0;
	++ntry;
	if (ranart_(&rndf77_1.nseed) > fff && ntry <= 20)
	{
		goto L10;
	}
	if ((r__1 = *am10 - .77f, dabs(r__1)) <= .01f && *dm1 > 1.07f || (r__2 = *am10 - 1.232f, dabs(r__2)) <= .01f && *dm1 > 1.47f || (r__3 = *am10 - 1.44f, dabs(r__3)) <= .01f && *dm1 > 2.14f)
	{
		goto L5;
	}
	if ((r__1 = *am20 - .77f, dabs(r__1)) <= .01f && *dm2 > 1.07f || (r__2 = *am20 - 1.232f, dabs(r__2)) <= .01f && *dm2 > 1.47f || (r__3 = *am20 - 1.44f, dabs(r__3)) <= .01f && *dm2 > 2.14f)
	{
		goto L5;
	}
	return 0;
}

double fmassd_(float *dmass)
{
	float ret_val, r__1, r__2, r__3, r__4, r__5;

	static float am0;
	extern double width_(float *);

	am0 = 1.232f;
	r__2 = *dmass;
	r__3 = am0;
	r__1 = r__2 * r__2 - r__3 * r__3;
	r__4 = am0;
	r__5 = width_(dmass);
	ret_val = am0 * width_(dmass) / (r__1 * r__1 + r__4 * r__4 * (r__5 * r__5));
	return ret_val;
}

double fmassn_(float *dmass)
{
	float ret_val, r__1, r__2, r__3, r__4, r__5;

	static float am0;
	extern double w1440_(float *);

	am0 = 1.44f;
	r__2 = *dmass;
	r__3 = am0;
	r__1 = r__2 * r__2 - r__3 * r__3;
	r__4 = am0;
	r__5 = w1440_(dmass);
	ret_val = am0 * w1440_(dmass) / (r__1 * r__1 + r__4 * r__4 * (r__5 * r__5));
	return ret_val;
}

double fmassr_(float *dmass)
{
	float ret_val, r__1, r__2, r__3, r__4, r__5;

	static float am0, wid;

	am0 = .77f;
	wid = .153f;
	r__2 = *dmass;
	r__3 = am0;
	r__1 = r__2 * r__2 - r__3 * r__3;
	r__4 = am0;
	r__5 = wid;
	ret_val = am0 * wid / (r__1 * r__1 + r__4 * r__4 * (r__5 * r__5));
	return ret_val;
}

int flow_(int *nt)
{
	int i__1, i__2;
	float r__1, r__2, r__3, r__4;

	int i_nint(float *);
	double sqrt(double), log(double);

	static int i__, j, m;
	static float e00;
	static int kk;
	static float y00;
	static int is;
	static float dy;
	static int iy, ly, npr, npt;
	static float ypr[161], dnuc, dypr;
	static int nrun;
	static float ycut1, ycut2, dnuck;
	static int nkaon;
	static float dnucp, ykaon[161];
	static int npion;
	static float ypion[161], pxpro[161], dykaon, pxkaon[161], dypion, pxpion[161];

	ycut1 = -2.6f;
	ycut2 = 2.6f;
	dy = .2f;
	r__1 = (ycut2 - ycut1) / dy;
	ly = i_nint(&r__1);
	for (kk = -80; kk <= 80; ++kk)
	{
		pxpion[kk + 80] = 0.f;
		pxpro[kk + 80] = 0.f;
		pxkaon[kk + 80] = 0.f;
	}
	i__1 = ly;
	for (j = -ly; j <= i__1; ++j)
	{
		ypion[j + 80] = 0.f;
		ykaon[j + 80] = 0.f;
		ypr[j + 80] = 0.f;
	}
	nkaon = 0;
	npr = 0;
	npion = 0;
	is = 0;
	i__1 = run_1.num;
	for (nrun = 1; nrun <= i__1; ++nrun)
	{
		is += rr_1.massr[nrun - 1];
		i__2 = rr_1.massr[nrun];
		for (j = 1; j <= i__2; ++j)
		{
			i__ = j + is;
			r__1 = bb_1.p[i__ * 3 - 3];
			r__2 = bb_1.p[i__ * 3 - 2];
			r__3 = bb_1.p[i__ * 3 - 1];
			r__4 = cc_1.e[i__ - 1];
			e00 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
			y00 = log((e00 + bb_1.p[i__ * 3 - 1]) / (e00 - bb_1.p[i__ * 3 - 1])) * .5f;
			if (dabs(y00) >= ycut2)
			{
				goto L20;
			}
			r__1 = y00 / dy;
			iy = i_nint(&r__1);
			if (abs(iy) >= 80)
			{
				goto L20;
			}
			if (cc_1.e[i__ - 1] == 0.f)
			{
				goto L20;
			}
			if (ee_1.lb[i__ - 1] >= 25)
			{
				goto L20;
			}
			if (ee_1.lb[i__ - 1] <= 5 && ee_1.lb[i__ - 1] >= 3)
			{
				goto L50;
			}
			if (ee_1.lb[i__ - 1] == 1 || ee_1.lb[i__ - 1] == 2)
			{
				goto L200;
			}
			if (ee_1.lb[i__ - 1] >= 6 && ee_1.lb[i__ - 1] <= 17)
			{
				goto L200;
			}
			if (ee_1.lb[i__ - 1] == 23)
			{
				goto L400;
			}
			goto L20;
		L50:
			++npion;
			ypion[iy + 80] += 1;
			pxpion[iy + 80] += bb_1.p[i__ * 3 - 3] / cc_1.e[i__ - 1];
			goto L20;
		L200:
			++npr;
			pxpro[iy + 80] += bb_1.p[i__ * 3 - 3] / cc_1.e[i__ - 1];
			ypr[iy + 80] += 1.f;
			goto L20;
		L400:
			++nkaon;
			ykaon[iy + 80] += 1.f;
			pxkaon[iy + 80] += bb_1.p[i__ * 3 - 3] / cc_1.e[i__ - 1];
		L20:;
		}
	}
	for (npt = -10; npt <= 10; ++npt)
	{
		if (ypr[npt + 80] == 0.f)
		{
			goto L101;
		}
		pxpro[npt + 80] = -pxpro[npt + 80] / ypr[npt + 80];
		dnuc = pxpro[npt + 80] / sqrt(ypr[npt + 80]);
	L101:
		if (ypion[npt + 80] == 0.f)
		{
			goto L102;
		}
		pxpion[npt + 80] = -pxpion[npt + 80] / ypion[npt + 80];
		dnucp = pxpion[npt + 80] / sqrt(ypion[npt + 80]);
	L102:
		if (ykaon[npt + 80] == 0.f)
		{
			goto L3;
		}
		pxkaon[npt + 80] = -pxkaon[npt + 80] / ykaon[npt + 80];
		dnuck = pxkaon[npt + 80] / sqrt(ykaon[npt + 80]);
	L3:;
	}
	i__2 = ly;
	for (m = -ly; m <= i__2; ++m)
	{
		dypr = 0.f;
		if (ypr[m + 80] != 0.f)
		{
			dypr = sqrt(ypr[m + 80]) / (float)nrun / dy;
		}
		ypr[m + 80] = ypr[m + 80] / (float)nrun / dy;
		dypion = 0.f;
		if (ypion[m + 80] != 0.f)
		{
			dypion = sqrt(ypion[m + 80]) / (float)nrun / dy;
		}
		ypion[m + 80] = ypion[m + 80] / (float)nrun / dy;
		dykaon = 0.f;
		if (ykaon[m + 80] != 0.f)
		{
			dykaon = sqrt(ykaon[m + 80]) / (float)nrun / dy;
		}
		ykaon[m + 80] = ykaon[m + 80] / (float)nrun / dy;
	}
	return 0;
}

double xppbar_(float *srt)
{
	float ret_val, r__1, r__2;
	double d__1;

	double sqrt(double), pow_dd(double *, double *);

	static float plab, plab2;

	ret_val = 1e-6f;
	r__2 = *srt;
	r__1 = r__2 * r__2 / 1.8766f - .9383f;
	plab2 = r__1 * r__1 - .88040689000000005f;
	if (plab2 > 0.f)
	{
		plab = sqrt(plab2);
		d__1 = (double)plab;
		ret_val = 67.f / pow_dd(&d__1, &c_b646);
		if (ret_val > 400.f)
		{
			ret_val = 400.f;
		}
	}
	return ret_val;
}

int pbarfs_(float *srt, int *npion, int *iseed)
{
	float r__1;

	double pow_ri(float *, int *);

	static int n;
	static float ene, pmax, pnpi[6];
	static int ntry;
	static float thisp, factor[6];
	extern double ranart_(int *);

	factor[1] = 1.f;
	factor[2] = .117f;
	factor[3] = .00327f;
	factor[4] = 3.58e-5f;
	factor[5] = 1.93e-7f;
	r__1 = *srt / .14f;
	ene = r__1 * (r__1 * r__1) / 59.217624386248559f;
	for (n = 2; n <= 6; ++n)
	{
		pnpi[n - 1] = pow_ri(&ene, &n) * factor[n - 1];
	}
	r__1 = max(pnpi[1], pnpi[2]), r__1 = max(r__1, pnpi[3]), r__1 = max(r__1, pnpi[4]);
	pmax = dmax(r__1, pnpi[5]);
	ntry = 0;
L10:
	*npion = (int)(ranart_(&rndf77_1.nseed) * 5) + 2;
	if (*npion > 6)
	{
		goto L10;
	}
	thisp = pnpi[*npion - 1] / pmax;
	++ntry;
	if (thisp < ranart_(&rndf77_1.nseed) && ntry <= 20)
	{
		goto L10;
	}
	return 0;
}

int xkkann_(float *srt, float *xsk1, float *xsk2, float *xsk3,
			float *xsk4, float *xsk5, float *xsk6, float *xsk7, float *xsk8, float *xsk9, float *xsk10, float *xsk11, float *sigk, float *rrkk)
{
	float r__1, r__2, r__3, r__4, r__5;

	double pow_dd(double *, double *), sqrt(double);

	static float s, pf2, pi2, xm1, xm2, fwdp, pkaon;
	extern double pipik_(float *);
	static float xpion0;

	r__1 = *srt;
	s = r__1 * r__1;
	*sigk = 1e-8f;
	*xsk1 = 0.f;
	*xsk2 = 0.f;
	*xsk3 = 0.f;
	*xsk4 = 0.f;
	*xsk5 = 0.f;
	*xsk6 = 0.f;
	*xsk7 = 0.f;
	*xsk8 = 0.f;
	*xsk9 = 0.f;
	*xsk10 = 0.f;
	*xsk11 = 0.f;
	xpion0 = pipik_(srt);
	xpion0 *= 2.f;
	pi2 = s * (s - .99201600000000001f);
	if (pi2 <= 0.f)
	{
		return 0;
	}
	xm1 = .14f;
	xm2 = .14f;
	r__1 = xm1 + xm2;
	r__2 = xm1 - xm2;
	pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
	if (pf2 > 0.f)
	{
		*xsk1 = pf2 * 2.25f / pi2 * xpion0;
	}
	xm1 = .14f;
	xm2 = .5473f;
	r__1 = xm1 + xm2;
	r__2 = xm1 - xm2;
	pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
	if (pf2 > 0.f)
	{
		*xsk4 = pf2 * .75f / pi2 * xpion0;
	}
	xm1 = .5473f;
	xm2 = .5473f;
	r__1 = xm1 + xm2;
	r__2 = xm1 - xm2;
	pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
	if (pf2 > 0.f)
	{
		*xsk10 = pf2 * .25f / pi2 * xpion0;
	}
	xpion0 = *rrkk;
	xm1 = .77f;
	xm2 = .77f;
	r__1 = xm1 + xm2;
	r__2 = xm1 - xm2;
	pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
	if (pf2 > 0.f)
	{
		*xsk5 = pf2 * 20.25f / pi2 * xpion0;
	}
	xm1 = .77f;
	xm2 = .7819f;
	r__1 = xm1 + xm2;
	r__2 = xm1 - xm2;
	pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
	if (pf2 > 0.f)
	{
		*xsk6 = pf2 * 6.75f / pi2 * xpion0;
	}
	xm1 = .7819f;
	xm2 = .7819f;
	r__1 = xm1 + xm2;
	r__2 = xm1 - xm2;
	pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
	if (pf2 > 0.f)
	{
		*xsk8 = pf2 * 2.25f / pi2 * xpion0;
	}
	fwdp = pow_dd(&c_b651, &c_b652) * 1.68f / 6.f / 1.02f / 1.02f;
	r__1 = *srt;
	pkaon = sqrt(r__1 * r__1 - .99201600000000001f) * .5f;
	r__1 = fwdp * 1.02f;
	r__3 = *srt;
	r__2 = r__3 * r__3 - 1.0404f;
	r__4 = fwdp * 1.02f;
	r__5 = pkaon;
	*xsk11 = r__1 * r__1 * 3.6688075497330002f / (r__2 * r__2 + r__4 * r__4) /
			 (r__5 * r__5);

	*sigk = *xsk1 + *xsk2 + *xsk3 + *xsk4 + *xsk5 + *xsk6 + *xsk7 + *xsk8 + *xsk9 + *xsk10 + *xsk11;
	return 0;
}

int xphib_(int *lb1, int *lb2, float *em1, float *em2,
		   float *srt, float *xsk1, float *xsk2, float *xsk3, float *xsk4, float *xsk5,
		   float *sigp)
{
	float r__1;
	double d__1;

	double pow_dd(double *, double *);

	static float xsk6, srrt;

	*sigp = 1e-8f;
	*xsk1 = 0.f;
	*xsk2 = 0.f;
	*xsk3 = 0.f;
	*xsk4 = 0.f;
	*xsk5 = 0.f;
	xsk6 = 0.f;
	srrt = *srt - (*em1 + *em2);
	*xsk1 = 8.f;

	if (*srt > 1.074417f)
	{
		d__1 = (double)srrt;
		*xsk2 = pow_dd(&d__1, &c_b654) * .0235f;
	}

	if (*srt > 1.36696f)
	{
		if (srrt < .7f)
		{
			d__1 = (double)srrt;
			*xsk3 = pow_dd(&d__1, &c_b655) * .0119f;
		}
		else
		{
			d__1 = (double)srrt;
			*xsk3 = pow_dd(&d__1, &c_b656) * .013f;
		}
	}

	if (*srt > 1.709457f)
	{
		if (srrt < .7f)
		{
			d__1 = (double)srrt;
			*xsk4 = pow_dd(&d__1, &c_b657) * .0166f;
		}
		else
		{
			d__1 = (double)srrt;
			*xsk4 = pow_dd(&d__1, &c_b658) * .0189f;
		}
	}

	if (*srt > 2.0019999999999998f)
	{
		if (srrt < .7f)
		{
			d__1 = (double)srrt;
			*xsk5 = pow_dd(&d__1, &c_b655) * .0119f;
		}
		else
		{
			d__1 = (double)srrt;
			*xsk5 = pow_dd(&d__1, &c_b656) * .013f;
		}
	}

	if (*lb1 >= 1 && *lb1 <= 2 || *lb2 >= 1 && *lb2 <= 2)
	{
		if (*srt > 1.6136999999999999f)
		{
			r__1 = srrt + 3.508f;
			xsk6 = 1.715f / (r__1 * r__1 - 12.138f);
		}
	}
	*sigp = *xsk1 + *xsk2 + *xsk3 + *xsk4 + *xsk5 + xsk6;
	return 0;
}

int crphib_(float *px, float *py, float *pz, float *srt, int *i1, int *i2, float *xsk1, float *xsk2, float *xsk3, float *xsk4,
			float *xsk5, float *sigp, int *iblock)
{
	float r__1, r__2, r__3, r__4, r__5;

	double sqrt(double), cos(double), sin(double);

	static float c1, s1, t1, x1, pr, em1, em2, ct1, pr2, px0, py0, pz0, st1;
	extern double ranart_(int *);
	extern int rotate_(float *, float *, float *, float *, float *, float *);

	px0 = *px;
	py0 = *py;
	pz0 = *pz;
	*iblock = 223;

	x1 = ranart_(&rndf77_1.nseed) * *sigp;
	*xsk2 = *xsk1 + *xsk2;
	*xsk3 = *xsk2 + *xsk3;
	*xsk4 = *xsk3 + *xsk4;
	*xsk5 = *xsk4 + *xsk5;

	if (x1 <= *xsk1)
	{
		*iblock = 20;
		goto L100;
	}
	else if (x1 <= *xsk2)
	{
		ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 3) + 3;
		ee_1.lb[*i2 - 1] = (int)(ranart_(&rndf77_1.nseed) * 2) + 1;
		cc_1.e[*i1 - 1] = .13496f;
		cc_1.e[*i2 - 1] = .939457f;
		goto L100;
	}
	else if (x1 <= *xsk3)
	{
		ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 3) + 3;
		ee_1.lb[*i2 - 1] = (int)(ranart_(&rndf77_1.nseed) * 4) + 6;
		cc_1.e[*i1 - 1] = .13496f;
		cc_1.e[*i2 - 1] = 1.232f;
		goto L100;
	}
	else if (x1 <= *xsk4)
	{
		ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 3) + 25;
		ee_1.lb[*i2 - 1] = (int)(ranart_(&rndf77_1.nseed) * 2) + 1;
		cc_1.e[*i1 - 1] = .77f;
		cc_1.e[*i2 - 1] = .939457f;
		goto L100;
	}
	else if (x1 <= *xsk5)
	{
		ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 3) + 25;
		ee_1.lb[*i2 - 1] = (int)(ranart_(&rndf77_1.nseed) * 4) + 6;
		cc_1.e[*i1 - 1] = .77f;
		cc_1.e[*i2 - 1] = 1.232f;
		goto L100;
	}
	else
	{
		ee_1.lb[*i1 - 1] = 23;
		ee_1.lb[*i2 - 1] = 14;
		cc_1.e[*i1 - 1] = .498f;
		cc_1.e[*i2 - 1] = 1.1157f;
		*iblock = 221;
	}
L100:
	em1 = cc_1.e[*i1 - 1];
	em2 = cc_1.e[*i2 - 1];
	r__2 = *srt;
	r__3 = em1;
	r__4 = em2;
	r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
	r__5 = em1 * em2;
	pr2 = r__1 * r__1 - r__5 * r__5 * 4.f;
	if (pr2 <= 0.f)
	{
		pr2 = 1e-8f;
	}
	pr = sqrt(pr2) / (*srt * 2.f);
	c1 = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
	t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
	r__1 = c1;
	s1 = sqrt(1.f - r__1 * r__1);
	ct1 = cos(t1);
	st1 = sin(t1);
	*pz = pr * c1;
	*px = pr * s1 * ct1;
	*py = pr * s1 * st1;
	rotate_(&px0, &py0, &pz0, px, py, pz);
	return 0;
}

int pibphi_(float *srt, int *lb1, int *lb2, float *em1,
			float *em2, float *xphi, float *xphin)
{
	float r__1, r__2, r__3, r__4, r__5, r__6;
	double d__1;

	double pow_dd(double *, double *);

	static float sig, srrt, xphid;

	*xphi = 0.f;
	*xphin = 0.f;
	xphid = 0.f;

	if (*lb1 >= 3 && *lb1 <= 5 || *lb2 >= 3 && *lb2 <= 5)
	{

		if (abs(*lb1) >= 1 && abs(*lb1) <= 2 || abs(*lb2) >= 1 && abs(*lb2) <=
																	  2)
		{
			if (*srt > 1.959457f)
			{
				srrt = *srt - 1.959457f;
				d__1 = (double)srrt;
				sig = pow_dd(&d__1, &c_b654) * .0235f;
				r__1 = *srt;
				r__2 = *srt;
				r__3 = *srt;
				r__4 = *em1 + *em2;
				r__5 = *srt;
				r__6 = *em1 - *em2;
				*xphin = sig * 1.f * (r__1 * r__1 - 3.839471734849f) * (r__2 * r__2 - .0064871748490000049f) / (r__3 * r__3 - r__4 * r__4) / (r__5 * r__5 - r__6 * r__6);
			}
			if (*srt > 2.2519999999999998f)
			{
				srrt = *srt - 2.2519999999999998f;
				d__1 = (double)srrt;
				sig = pow_dd(&d__1, &c_b654) * .0235f;
				r__1 = *srt;
				r__2 = *srt;
				r__3 = *srt;
				r__4 = *em1 + *em2;
				r__5 = *srt;
				r__6 = *em1 - *em2;
				xphid = sig * 4.f * (r__1 * r__1 - 5.0715039999999991f) * (r__2 * r__2 - .044943999999999984f) / (r__3 * r__3 - r__4 * r__4) / (r__5 * r__5 - r__6 * r__6);
			}
		}
		else
		{
			if (*srt > 1.959457f)
			{
				srrt = *srt - 1.959457f;
				if (srrt < .7f)
				{
					d__1 = (double)srrt;
					sig = pow_dd(&d__1, &c_b655) * .0119f;
				}
				else
				{
					d__1 = (double)srrt;
					sig = pow_dd(&d__1, &c_b656) * .013f;
				}
				r__1 = *srt;
				r__2 = *srt;
				r__3 = *srt;
				r__4 = *em1 + *em2;
				r__5 = *srt;
				r__6 = *em1 - *em2;
				*xphin = sig * .25f * (r__1 * r__1 - 3.839471734849f) * (r__2 * r__2 - .0064871748490000049f) / (r__3 * r__3 - r__4 * r__4) / (r__5 * r__5 - r__6 * r__6);
			}
			if (*srt > 2.2519999999999998f)
			{
				srrt = *srt - 2.2519999999999998f;
				if (srrt < .7f)
				{
					d__1 = (double)srrt;
					sig = pow_dd(&d__1, &c_b655) * .0119f;
				}
				else
				{
					d__1 = (double)srrt;
					sig = pow_dd(&d__1, &c_b656) * .013f;
				}
				r__1 = *srt;
				r__2 = *srt;
				r__3 = *srt;
				r__4 = *em1 + *em2;
				r__5 = *srt;
				r__6 = *em1 - *em2;
				xphid = sig * 1.f * (r__1 * r__1 - 5.0715039999999991f) * (r__2 * r__2 - .044943999999999984f) / (r__3 * r__3 - r__4 * r__4) / (r__5 * r__5 - r__6 * r__6);
			}
		}
	}
	else
	{

		if (abs(*lb1) >= 1 && abs(*lb1) <= 2 || abs(*lb2) >= 1 && abs(*lb2) <=
																	  2)
		{

			if (*srt > 1.959457f)
			{
				srrt = *srt - 1.959457f;
				if (srrt < .7f)
				{
					d__1 = (double)srrt;
					sig = pow_dd(&d__1, &c_b657) * .0166f;
				}
				else
				{
					d__1 = (double)srrt;
					sig = pow_dd(&d__1, &c_b658) * .0189f;
				}
				r__1 = *srt;
				r__2 = *srt;
				r__3 = *srt;
				r__4 = *em1 + *em2;
				r__5 = *srt;
				r__6 = *em1 - *em2;
				*xphin = sig * .33333333333333331f * (r__1 * r__1 - 3.839471734849f) * (r__2 * r__2 - .0064871748490000049f) / (r__3 * r__3 - r__4 * r__4) /
						 (r__5 * r__5 - r__6 * r__6);
			}
			if (*srt > 2.2519999999999998f)
			{
				srrt = *srt - 2.2519999999999998f;
				if (srrt < .7f)
				{
					d__1 = (double)srrt;
					sig = pow_dd(&d__1, &c_b657) * .0166f;
				}
				else
				{
					d__1 = (double)srrt;
					sig = pow_dd(&d__1, &c_b658) * .0189f;
				}
				r__1 = *srt;
				r__2 = *srt;
				r__3 = *srt;
				r__4 = *em1 + *em2;
				r__5 = *srt;
				r__6 = *em1 - *em2;
				xphid = sig * 1.3333333333333333f * (r__1 * r__1 - 5.0715039999999991f) * (r__2 * r__2 - .044943999999999984f) / (r__3 * r__3 - r__4 * r__4) /
						(r__5 * r__5 - r__6 * r__6);
			}
		}
		else
		{
			if (*srt > 1.959457f)
			{
				srrt = *srt - 1.959457f;
				if (srrt < .7f)
				{
					d__1 = (double)srrt;
					sig = pow_dd(&d__1, &c_b655) * .0119f;
				}
				else
				{
					d__1 = (double)srrt;
					sig = pow_dd(&d__1, &c_b656) * .013f;
				}
				r__1 = *srt;
				r__2 = *srt;
				r__3 = *srt;
				r__4 = *em1 + *em2;
				r__5 = *srt;
				r__6 = *em1 - *em2;
				*xphin = sig * .083333333333333329f * (r__1 * r__1 - 3.839471734849f) * (r__2 * r__2 - .0064871748490000049f) / (r__3 * r__3 - r__4 * r__4) /
						 (r__5 * r__5 - r__6 * r__6);
			}
			if (*srt > 2.2519999999999998f)
			{
				srrt = *srt - 2.2519999999999998f;
				if (srrt < .7f)
				{
					d__1 = (double)srrt;
					sig = pow_dd(&d__1, &c_b655) * .0119f;
				}
				else
				{
					d__1 = (double)srrt;
					sig = pow_dd(&d__1, &c_b656) * .013f;
				}
				r__1 = *srt;
				r__2 = *srt;
				r__3 = *srt;
				r__4 = *em1 + *em2;
				r__5 = *srt;
				r__6 = *em1 - *em2;
				xphid = sig * .33333333333333331f * (r__1 * r__1 - 5.0715039999999991f) * (r__2 * r__2 - .044943999999999984f) / (r__3 * r__3 - r__4 * r__4) /
						(r__5 * r__5 - r__6 * r__6);
			}
		}
	}
	*xphin /= 10.f;
	xphid /= 10.f;
	*xphi = *xphin + xphid;
	return 0;
}

int phimes_(int *i1, int *i2, float *srt, float *xsk1,
			float *xsk2, float *xsk3, float *xsk4, float *xsk5, float *xsk6, float *xsk7, float *sigphi)
{
	float r__1, r__2;
	double d__1;

	double sqrt(double), pow_dd(double *, double *);

	static float s;
	static int lb1, lb2;
	static float em1, em2, pff, pii, srr, srr1, srr2, akap, srrt;

	r__1 = *srt;
	s = r__1 * r__1;
	*sigphi = 1e-8f;
	*xsk1 = 0.f;
	*xsk2 = 0.f;
	*xsk3 = 0.f;
	*xsk4 = 0.f;
	*xsk5 = 0.f;
	*xsk6 = 0.f;
	*xsk7 = 0.f;
	em1 = cc_1.e[*i1 - 1];
	em2 = cc_1.e[*i2 - 1];
	lb1 = ee_1.lb[*i1 - 1];
	lb2 = ee_1.lb[*i2 - 1];
	akap = .498f;
	*xsk1 = 5.f;
	r__1 = em1 + em2;
	r__2 = em1 - em2;
	pii = sqrt((s - r__1 * r__1) * (s - r__2 * r__2));
	if (lb1 == 23 || lb2 == 23 || lb1 == 21 || lb2 == 21)
	{
		if (*srt > akap + .13496f)
		{
			r__1 = akap + .13496f;
			r__2 = .13496f - akap;
			pff = sqrt((s - r__1 * r__1) * (s - r__2 * r__2));
			*xsk2 = pff * 195.639f / pii / 32.f / 3.1415926f / s;
		}
		if (*srt > akap + .77f)
		{
			r__1 = akap + .77f;
			r__2 = .77f - akap;
			pff = sqrt((s - r__1 * r__1) * (s - r__2 * r__2));
			*xsk3 = pff * 526.702f / pii / 32.f / 3.1415926f / s;
		}
		if (*srt > akap + .7819f)
		{
			r__1 = akap + .7819f;
			r__2 = .7819f - akap;
			pff = sqrt((s - r__1 * r__1) * (s - r__2 * r__2));
			*xsk4 = pff * 355.429f / pii / 32.f / 3.1415926f / s;
		}
		if (*srt > 1.02996f)
		{
			pff = sqrt((s - 1.0608176015999999f) * (s - .57766080160000011f));
			*xsk5 = pff * 2047.042f / pii / 32.f / 3.1415926f / s;
		}
		if (*srt > 1.665f)
		{
			pff = sqrt((s - 2.7722250000000002f) * (s - .015625f));
			*xsk6 = pff * 1371.257f / pii / 32.f / 3.1415926f / s;
		}
		if (*srt > 1.6769000000000001f)
		{
			pff = sqrt((s - 2.81199361f) * (s - .012791609999999995f));
			*xsk7 = pff * 482.292f / pii / 32.f / 3.1415926f / s;
		}
	}
	else if (abs(lb1) == 30 || abs(lb2) == 30)
	{
		if (*srt > akap + .13496f)
		{
			r__1 = akap + .13496f;
			r__2 = .13496f - akap;
			pff = sqrt((s - r__1 * r__1) * (s - r__2 * r__2));
			*xsk2 = pff * 372.378f / pii / 32.f / 3.1415926f / s;
		}
		if (*srt > akap + .77f)
		{
			r__1 = akap + .77f;
			r__2 = .77f - akap;
			pff = sqrt((s - r__1 * r__1) * (s - r__2 * r__2));
			*xsk3 = pff * 1313.96f / pii / 32.f / 3.1415926f / s;
		}
		if (*srt > akap + .7819f)
		{
			r__1 = akap + .7819f;
			r__2 = .7819f - akap;
			pff = sqrt((s - r__1 * r__1) * (s - r__2 * r__2));
			*xsk4 = pff * 440.558f / pii / 32.f / 3.1415926f / s;
		}
		if (*srt > 1.02996f)
		{
			pff = sqrt((s - 1.0608176015999999f) * (s - .57766080160000011f));
			*xsk5 = pff * 1496.692f / pii / 32.f / 3.1415926f / s;
		}
		if (*srt > 1.665f)
		{
			pff = sqrt((s - 2.7722250000000002f) * (s - .015625f));
			*xsk6 = pff * 6999.84f / pii / 32.f / 3.1415926f / s;
		}
		if (*srt > 1.6769000000000001f)
		{
			pff = sqrt((s - 2.81199361f) * (s - .012791609999999995f));
			*xsk7 = pff * 1698.903f / pii / 32.f / 3.1415926f / s;
		}
	}
	else
	{

		srr1 = em1 + em2;
		if (*srt > akap + akap)
		{
			srrt = *srt - srr1;
			if (srrt < .3f && srrt > .01f)
			{
				d__1 = (double)srrt;
				*xsk2 = 1.69f / (pow_dd(&d__1, &c_b529) - .407f);
			}
			else
			{
				d__1 = (double)srrt;
				*xsk2 = pow_dd(&d__1, &c_b530) * .008f + 3.74f;
			}
		}
		if (*srt > akap + .895f)
		{
			srr2 = akap + .895f;
			srr = dmax(srr1, srr2);
			srrt = *srt - srr;
			if (srrt < .3f && srrt > .01f)
			{
				d__1 = (double)srrt;
				*xsk3 = 1.69f / (pow_dd(&d__1, &c_b529) - .407f);
			}
			else
			{
				d__1 = (double)srrt;
				*xsk3 = pow_dd(&d__1, &c_b530) * .008f + 3.74f;
			}
		}
		if (*srt > 1.79f)
		{
			srr2 = 1.79f;
			srr = dmax(srr1, srr2);
			srrt = *srt - srr;
			if (srrt < .3f && srrt > .01f)
			{
				d__1 = (double)srrt;
				*xsk4 = 1.69f / (pow_dd(&d__1, &c_b529) - .407f);
			}
			else
			{
				d__1 = (double)srrt;
				*xsk4 = pow_dd(&d__1, &c_b530) * .008f + 3.74f;
			}
		}
	}
	*sigphi = *xsk1 + *xsk2 + *xsk3 + *xsk4 + *xsk5 + *xsk6 + *xsk7;
	return 0;
}

int crphim_(float *px, float *py, float *pz, float *srt, int *i1, int *i2, float *xsk1, float *xsk2, float *xsk3, float *xsk4,
			float *xsk5, float *xsk6, float *sigphi, int *ikkg, int *ikkl,
			int *iblock)
{
	float r__1, r__2, r__3, r__4, r__5;

	double sqrt(double), cos(double), sin(double);

	static float c1, s1, t1, x1, pr;
	static int lb1, lb2;
	static float em1, em2, ct1, pr2, px0, py0, pz0, st1;
	static int iad1, iad2;
	extern double ranart_(int *);
	extern int rotate_(float *, float *, float *, float *, float *, float *);

	px0 = *px;
	py0 = *py;
	pz0 = *pz;
	lb1 = ee_1.lb[*i1 - 1];
	lb2 = ee_1.lb[*i2 - 1];
	x1 = ranart_(&rndf77_1.nseed) * *sigphi;
	*xsk2 = *xsk1 + *xsk2;
	*xsk3 = *xsk2 + *xsk3;
	*xsk4 = *xsk3 + *xsk4;
	*xsk5 = *xsk4 + *xsk5;
	*xsk6 = *xsk5 + *xsk6;
	if (x1 <= *xsk1)
	{
		*iblock = 20;
		goto L100;
	}
	else
	{

		if (lb1 == 23 || lb1 == 21 || abs(lb1) == 30 || lb2 == 23 || lb2 == 21 || abs(lb2) == 30)
		{

			if (lb1 == 23 || lb2 == 23)
			{
				*ikkl = 1;
				*iblock = 224;
				iad1 = 23;
				iad2 = 30;
			}
			else if (lb1 == 30 || lb2 == 30)
			{
				*ikkl = 0;
				*iblock = 226;
				iad1 = 23;
				iad2 = 30;
			}
			else if (lb1 == 21 || lb2 == 21)
			{
				*ikkl = 1;
				*iblock = 124;
				iad1 = 21;
				iad2 = -30;
			}
			else
			{
				*ikkl = 0;
				*iblock = 126;
				iad1 = 21;
				iad2 = -30;
			}
			if (x1 <= *xsk2)
			{
				ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 3) +
								   3;
				ee_1.lb[*i2 - 1] = iad1;
				cc_1.e[*i1 - 1] = .13496f;
				cc_1.e[*i2 - 1] = .498f;
				*ikkg = 1;
				goto L100;
			}
			else if (x1 <= *xsk3)
			{
				ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 3) +
								   25;
				ee_1.lb[*i2 - 1] = iad1;
				cc_1.e[*i1 - 1] = .77f;
				cc_1.e[*i2 - 1] = .498f;
				*ikkg = 1;
				goto L100;
			}
			else if (x1 <= *xsk4)
			{
				ee_1.lb[*i1 - 1] = 28;
				ee_1.lb[*i2 - 1] = iad1;
				cc_1.e[*i1 - 1] = .7819f;
				cc_1.e[*i2 - 1] = .498f;
				*ikkg = 1;
				goto L100;
			}
			else if (x1 <= *xsk5)
			{
				ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 3) +
								   3;
				ee_1.lb[*i2 - 1] = iad2;
				cc_1.e[*i1 - 1] = .13496f;
				cc_1.e[*i2 - 1] = .895f;
				*ikkg = 0;
				++(*iblock);
				goto L100;
			}
			else if (x1 <= *xsk6)
			{
				ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 3) +
								   25;
				ee_1.lb[*i2 - 1] = iad2;
				cc_1.e[*i1 - 1] = .77f;
				cc_1.e[*i2 - 1] = .895f;
				*ikkg = 0;
				++(*iblock);
				goto L100;
			}
			else
			{
				ee_1.lb[*i1 - 1] = 28;
				ee_1.lb[*i2 - 1] = iad2;
				cc_1.e[*i1 - 1] = .7819f;
				cc_1.e[*i2 - 1] = .895f;
				*ikkg = 0;
				++(*iblock);
				goto L100;
			}
		}
		else
		{
			*iblock = 223;
			if (x1 <= *xsk2)
			{
				ee_1.lb[*i1 - 1] = 23;
				ee_1.lb[*i2 - 1] = 21;
				cc_1.e[*i1 - 1] = .498f;
				cc_1.e[*i2 - 1] = .498f;
				*ikkg = 2;
				*ikkl = 0;
				goto L100;
			}
			else if (x1 <= *xsk3)
			{
				ee_1.lb[*i1 - 1] = 23;
				ee_1.lb[*i2 - 1] = -30;
				if (ranart_(&rndf77_1.nseed) <= .5f)
				{
					ee_1.lb[*i1 - 1] = 21;
					ee_1.lb[*i2 - 1] = 30;
				}
				cc_1.e[*i1 - 1] = .498f;
				cc_1.e[*i2 - 1] = .895f;
				*ikkg = 1;
				*ikkl = 0;
				goto L100;
			}
			else if (x1 <= *xsk4)
			{
				ee_1.lb[*i1 - 1] = 30;
				ee_1.lb[*i2 - 1] = -30;
				cc_1.e[*i1 - 1] = .895f;
				cc_1.e[*i2 - 1] = .895f;
				*ikkg = 0;
				*ikkl = 0;
				goto L100;
			}
		}
	}

L100:
	em1 = cc_1.e[*i1 - 1];
	em2 = cc_1.e[*i2 - 1];
	r__2 = *srt;
	r__3 = em1;
	r__4 = em2;
	r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
	r__5 = em1 * em2;
	pr2 = r__1 * r__1 - r__5 * r__5 * 4.f;
	if (pr2 <= 0.f)
	{
		pr2 = 1e-8f;
	}
	pr = sqrt(pr2) / (*srt * 2.f);
	c1 = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
	t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
	r__1 = c1;
	s1 = sqrt(1.f - r__1 * r__1);
	ct1 = cos(t1);
	st1 = sin(t1);
	*pz = pr * c1;
	*px = pr * s1 * ct1;
	*py = pr * s1 * st1;
	rotate_(&px0, &py0, &pz0, px, py, pz);
	return 0;
}

int xkhype_(int *i1, int *i2, float *srt, float *xky1,
			float *xky2, float *xky3, float *xky4, float *xky5, float *xky6, float *xky7, float *xky8, float *xky9, float *xky10, float *xky11, float *xky12,
			float *xky13, float *xky14, float *xky15, float *xky16, float *xky17, float *sigk)
{
	float r__1, r__2;

	static float s;
	static int lb1, lb2;
	static float pf2, pi2, xm1, xm2, ddf, sig, srrt;
	extern double pnlka_(float *), pnska_(float *);
	static float xkaon0;

	r__1 = *srt;
	s = r__1 * r__1;
	*sigk = 1e-8f;
	*xky1 = 0.f;
	*xky2 = 0.f;
	*xky3 = 0.f;
	*xky4 = 0.f;
	*xky5 = 0.f;
	*xky6 = 0.f;
	*xky7 = 0.f;
	*xky8 = 0.f;
	*xky9 = 0.f;
	*xky10 = 0.f;
	*xky11 = 0.f;
	*xky12 = 0.f;
	*xky13 = 0.f;
	*xky14 = 0.f;
	*xky15 = 0.f;
	*xky16 = 0.f;
	*xky17 = 0.f;
	lb1 = ee_1.lb[*i1 - 1];
	lb2 = ee_1.lb[*i2 - 1];
	if (abs(lb1) == 14 || abs(lb2) == 14)
	{
		xkaon0 = pnlka_(srt);
		xkaon0 *= 2.f;
		pi2 = (s - 2.6049960000000003f) * (s - .38192400000000015f);
	}
	else
	{
		xkaon0 = pnska_(srt);
		xkaon0 *= 2.f;
		pi2 = (s - 2.8594810000000002f) * (s - .48302500000000009f);
	}
	if (pi2 <= 0.f)
	{
		return 0;
	}
	xm1 = .14f;
	xm2 = .93828f;
	r__1 = xm1 + xm2;
	r__2 = xm1 - xm2;
	pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
	if (pf2 > 0.f)
	{
		*xky1 = pf2 * 3.f / pi2 * xkaon0;
	}
	xm1 = .14f;
	xm2 = 1.232f;
	r__1 = xm1 + xm2;
	r__2 = xm1 - xm2;
	pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
	if (pf2 > 0.f)
	{
		*xky2 = pf2 * 12.f / pi2 * xkaon0;
	}
	xm1 = .14f;
	xm2 = 1.44f;
	r__1 = xm1 + xm2;
	r__2 = xm1 - xm2;
	pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
	if (pf2 > 0.f)
	{
		*xky3 = pf2 * 3.f / pi2 * xkaon0;
	}
	xm1 = .14f;
	xm2 = 1.535f;
	r__1 = xm1 + xm2;
	r__2 = xm1 - xm2;
	pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
	if (pf2 > 0.f)
	{
		*xky4 = pf2 * 3.f / pi2 * xkaon0;
	}
	xm1 = .769f;
	xm2 = .93828f;
	r__1 = xm1 + xm2;
	r__2 = xm1 - xm2;
	pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
	if (pf2 > 0.f)
	{
		*xky5 = pf2 * 9.f / pi2 * xkaon0;
	}
	xm1 = .769f;
	xm2 = 1.232f;
	r__1 = xm1 + xm2;
	r__2 = xm1 - xm2;
	pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
	if (pf2 > 0.f)
	{
		*xky6 = pf2 * 36.f / pi2 * xkaon0;
	}
	xm1 = .769f;
	xm2 = 1.44f;
	r__1 = xm1 + xm2;
	r__2 = xm1 - xm2;
	pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
	if (pf2 > 0.f)
	{
		*xky7 = pf2 * 9.f / pi2 * xkaon0;
	}
	xm1 = .769f;
	xm2 = 1.535f;
	r__1 = xm1 + xm2;
	r__2 = xm1 - xm2;
	pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
	if (pf2 > 0.f)
	{
		*xky8 = pf2 * 9.f / pi2 * xkaon0;
	}
	xm1 = .782f;
	xm2 = .93828f;
	r__1 = xm1 + xm2;
	r__2 = xm1 - xm2;
	pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
	if (pf2 > 0.f)
	{
		*xky9 = pf2 * 3.f / pi2 * xkaon0;
	}
	xm1 = .782f;
	xm2 = 1.232f;
	r__1 = xm1 + xm2;
	r__2 = xm1 - xm2;
	pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
	if (pf2 > 0.f)
	{
		*xky10 = pf2 * 12.f / pi2 * xkaon0;
	}
	xm1 = .782f;
	xm2 = 1.44f;
	r__1 = xm1 + xm2;
	r__2 = xm1 - xm2;
	pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
	if (pf2 > 0.f)
	{
		*xky11 = pf2 * 3.f / pi2 * xkaon0;
	}
	xm1 = .782f;
	xm2 = 1.535f;
	r__1 = xm1 + xm2;
	r__2 = xm1 - xm2;
	pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
	if (pf2 > 0.f)
	{
		*xky12 = pf2 * 3.f / pi2 * xkaon0;
	}
	xm1 = .5473f;
	xm2 = .93828f;
	r__1 = xm1 + xm2;
	r__2 = xm1 - xm2;
	pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
	if (pf2 > 0.f)
	{
		*xky13 = pf2 * 1.f / pi2 * xkaon0;
	}
	xm1 = .5473f;
	xm2 = 1.232f;
	r__1 = xm1 + xm2;
	r__2 = xm1 - xm2;
	pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
	if (pf2 > 0.f)
	{
		*xky14 = pf2 * 4.f / pi2 * xkaon0;
	}
	xm1 = .5473f;
	xm2 = 1.44f;
	r__1 = xm1 + xm2;
	r__2 = xm1 - xm2;
	pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
	if (pf2 > 0.f)
	{
		*xky15 = pf2 * 1.f / pi2 * xkaon0;
	}
	xm1 = .5473f;
	xm2 = 1.535f;
	r__1 = xm1 + xm2;
	r__2 = xm1 - xm2;
	pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
	if (pf2 > 0.f)
	{
		*xky16 = pf2 * 1.f / pi2 * xkaon0;
	}
	if (lb1 == 14 || lb2 == 14)
	{
		if (*srt > 1.959457f)
		{
			srrt = *srt - 1.959457f;
			r__1 = srrt + 3.508f;
			sig = 1.715f / (r__1 * r__1 - 12.138f);
			xm1 = .939457f;
			xm2 = 1.02f;
			r__1 = xm1 + xm2;
			r__2 = xm1 - xm2;
			pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
			*xky17 = pf2 * 3.f / pi2 * sig / 10.f;
		}
	}
	if (abs(lb1) >= 15 && abs(lb1) <= 17 || abs(lb2) >= 15 && abs(lb2) <= 17)
	{
		ddf = 3.f;
		*xky1 /= ddf;
		*xky2 /= ddf;
		*xky3 /= ddf;
		*xky4 /= ddf;
		*xky5 /= ddf;
		*xky6 /= ddf;
		*xky7 /= ddf;
		*xky8 /= ddf;
		*xky9 /= ddf;
		*xky10 /= ddf;
		*xky11 /= ddf;
		*xky12 /= ddf;
		*xky13 /= ddf;
		*xky14 /= ddf;
		*xky15 /= ddf;
		*xky16 /= ddf;
	}
	*sigk = *xky1 + *xky2 + *xky3 + *xky4 + *xky5 + *xky6 + *xky7 + *xky8 + *xky9 + *xky10 + *xky11 + *xky12 + *xky13 + *xky14 + *xky15 + *xky16 + *xky17;
	return 0;
}

int ppbdat_(void)
{
	return 0;
}

int getnst_(float *srt)
{
	int i__1;
	float r__1, r__2, r__3;

	static int i__;
	static float pf2;

	r__1 = *srt;
	ppb1_1.s = r__1 * r__1;
	ppbmas_1.nstate = 0;
	ppb1_1.wtot = 0.f;
	if (*srt <= ppbmas_1.thresh[0])
	{
		return 0;
	}
	for (i__ = 1; i__ <= 15; ++i__)
	{
		ppbmas_1.weight[i__ - 1] = 0.f;
		if (*srt > ppbmas_1.thresh[i__ - 1])
		{
			ppbmas_1.nstate = i__;
		}
	}
	i__1 = ppbmas_1.nstate;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		r__1 = ppbmas_1.ppbm[i__ - 1] + ppbmas_1.ppbm[i__ + 14];
		r__2 = ppbmas_1.ppbm[i__ - 1] - ppbmas_1.ppbm[i__ + 14];
		pf2 = (ppb1_1.s - r__1 * r__1) * (ppb1_1.s - r__2 * r__2) / 4 /
			  ppb1_1.s;
		ppbmas_1.weight[i__ - 1] = pf2 * ppbmas_1.niso[i__ - 1];
		ppb1_1.wtot += ppbmas_1.weight[i__ - 1];
	}
	r__1 = *srt / .14f;
	ppb1_1.ene = r__1 * (r__1 * r__1) / 59.217624386248559f;
	r__1 = ppb1_1.ene;
	r__2 = ppb1_1.ene;
	r__3 = ppb1_1.ene, r__3 *= r__3;
	ppb1_1.fsum = ppb1_1.factr2[1] + ppb1_1.factr2[2] * ppb1_1.ene +
				  ppb1_1.factr2[3] * (r__1 * r__1) + ppb1_1.factr2[4] * (r__2 * (r__2 * r__2)) + ppb1_1.factr2[5] * (r__3 * r__3);
	return 0;
}

double ppbbar_(float *srt)
{
	float ret_val;

	static float pi2, sppb2p;
	extern double xppbar_(float *);

	sppb2p = xppbar_(srt) * ppb1_1.factr2[1] / ppb1_1.fsum;
	pi2 = (ppb1_1.s - .078400000000000011f) / 4;
	ret_val = sppb2p * .44444444444444442f / pi2 * ppb1_1.wtot;
	return ret_val;
}

double prbbar_(float *srt)
{
	float ret_val;

	static float pi2, sppb3p;
	extern double xppbar_(float *);

	sppb3p = xppbar_(srt) * ppb1_1.factr2[2] * ppb1_1.ene / ppb1_1.fsum;
	pi2 = (ppb1_1.s - .82810000000000006f) * (ppb1_1.s - .39690000000000003f) / 4 / ppb1_1.s;
	ret_val = sppb3p * .14814814814814814f / pi2 * ppb1_1.wtot;
	return ret_val;
}

double rrbbar_(float *srt)
{
	float ret_val, r__1;

	static float pi2, sppb4p;
	extern double xppbar_(float *);

	r__1 = ppb1_1.ene;
	sppb4p = xppbar_(srt) * ppb1_1.factr2[3] * (r__1 * r__1) / ppb1_1.fsum;
	pi2 = (ppb1_1.s - 2.3715999999999999f) / 4;
	ret_val = sppb4p / 2 * .049382716049382713f / pi2 * ppb1_1.wtot;
	return ret_val;
}

double pobbar_(float *srt)
{
	float ret_val, r__1;

	static float pi2, sppb4p;
	extern double xppbar_(float *);

	r__1 = ppb1_1.ene;
	sppb4p = xppbar_(srt) * ppb1_1.factr2[3] * (r__1 * r__1) / ppb1_1.fsum;
	pi2 = (ppb1_1.s - .85008400000000006f) * (ppb1_1.s - .41216400000000003f) / 4 / ppb1_1.s;
	ret_val = sppb4p / 2 * .44444444444444442f / pi2 * ppb1_1.wtot;
	return ret_val;
}

double robbar_(float *srt)
{
	float ret_val, r__1;

	static float pi2, sppb5p;
	extern double xppbar_(float *);

	r__1 = ppb1_1.ene;
	sppb5p = xppbar_(srt) * ppb1_1.factr2[4] * (r__1 * (r__1 * r__1)) /
			 ppb1_1.fsum;
	pi2 = (ppb1_1.s - 2.4087040000000002f) * (ppb1_1.s - 1.4400000000000025e-4f) / 4 / ppb1_1.s;
	ret_val = sppb5p * .14814814814814814f / pi2 * ppb1_1.wtot;
	return ret_val;
}

double oobbar_(float *srt)
{
	float ret_val, r__1;

	static float pi2, sppb6p;
	extern double xppbar_(float *);

	r__1 = ppb1_1.ene, r__1 *= r__1;
	sppb6p = xppbar_(srt) * ppb1_1.factr2[5] * (r__1 * r__1) / ppb1_1.fsum;
	pi2 = (ppb1_1.s - 2.4460960000000003f) / 4;
	ret_val = sppb6p * .44444444444444442f / pi2 * ppb1_1.wtot;
	return ret_val;
}

int bbarfs_(int *lbb1, int *lbb2, float *ei1, float *ei2, int *iblock, int *iseed)
{
	int i__1;

	static int i__;
	static float rd, rd1, rd2;
	static int ifs;
	static float wsum;
	extern double ranart_(int *);

	rd = ranart_(&rndf77_1.nseed);
	wsum = 0.f;
	i__1 = ppbmas_1.nstate;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		wsum += ppbmas_1.weight[i__ - 1];
		if (rd <= wsum / ppb1_1.wtot)
		{
			ifs = i__;
			*ei1 = ppbmas_1.ppbm[i__ - 1];
			*ei2 = ppbmas_1.ppbm[i__ + 14];
			goto L10;
		}
	}
L10:
	if (ifs == 1)
	{
		*iblock = 1801;
		*lbb1 = -1;
		*lbb2 = 1;
	}
	else if (ifs == 2)
	{
		if (ranart_(&rndf77_1.nseed) <= .5f)
		{
			*iblock = 18021;
			*lbb1 = -1;
			*lbb2 = 2;
		}
		else
		{
			*iblock = 18022;
			*lbb1 = 1;
			*lbb2 = -2;
		}
	}
	else if (ifs == 3)
	{
		*iblock = 1803;
		*lbb1 = -2;
		*lbb2 = 2;
	}
	else if (ifs == 4 || ifs == 5)
	{
		rd = ranart_(&rndf77_1.nseed);
		if (rd <= .5f)
		{
			if (ifs == 4)
			{
				*iblock = 18041;
				*lbb1 = -1;
			}
			else
			{
				*iblock = 18051;
				*lbb1 = -2;
			}
			rd2 = ranart_(&rndf77_1.nseed);
			if (rd2 <= .25f)
			{
				*lbb2 = 6;
			}
			else if (rd2 <= .5f)
			{
				*lbb2 = 7;
			}
			else if (rd2 <= .75f)
			{
				*lbb2 = 8;
			}
			else
			{
				*lbb2 = 9;
			}
		}
		else
		{
			if (ifs == 4)
			{
				*iblock = 18042;
				*lbb1 = 1;
			}
			else
			{
				*iblock = 18052;
				*lbb1 = 2;
			}
			rd2 = ranart_(&rndf77_1.nseed);
			if (rd2 <= .25f)
			{
				*lbb2 = -6;
			}
			else if (rd2 <= .5f)
			{
				*lbb2 = -7;
			}
			else if (rd2 <= .75f)
			{
				*lbb2 = -8;
			}
			else
			{
				*lbb2 = -9;
			}
		}
	}
	else if (ifs == 6 || ifs == 7)
	{
		rd = ranart_(&rndf77_1.nseed);
		if (rd <= .5f)
		{
			if (ifs == 6)
			{
				*iblock = 18061;
				*lbb1 = -1;
			}
			else
			{
				*iblock = 18071;
				*lbb1 = -2;
			}
			rd2 = ranart_(&rndf77_1.nseed);
			if (rd2 <= .5f)
			{
				*lbb2 = 10;
			}
			else
			{
				*lbb2 = 11;
			}
		}
		else
		{
			if (ifs == 6)
			{
				*iblock = 18062;
				*lbb1 = 1;
			}
			else
			{
				*iblock = 18072;
				*lbb1 = 2;
			}
			rd2 = ranart_(&rndf77_1.nseed);
			if (rd2 <= .5f)
			{
				*lbb2 = -10;
			}
			else
			{
				*lbb2 = -11;
			}
		}
	}
	else if (ifs == 8)
	{
		*iblock = 1808;
		rd1 = ranart_(&rndf77_1.nseed);
		if (rd1 <= .25f)
		{
			*lbb1 = 6;
		}
		else if (rd1 <= .5f)
		{
			*lbb1 = 7;
		}
		else if (rd1 <= .75f)
		{
			*lbb1 = 8;
		}
		else
		{
			*lbb1 = 9;
		}
		rd2 = ranart_(&rndf77_1.nseed);
		if (rd2 <= .25f)
		{
			*lbb2 = -6;
		}
		else if (rd2 <= .5f)
		{
			*lbb2 = -7;
		}
		else if (rd2 <= .75f)
		{
			*lbb2 = -8;
		}
		else
		{
			*lbb2 = -9;
		}
	}
	else if (ifs == 9 || ifs == 10)
	{
		rd = ranart_(&rndf77_1.nseed);
		if (rd <= .5f)
		{
			if (ifs == 9)
			{
				*iblock = 18091;
				*lbb1 = -1;
			}
			else
			{
				*iblock = 18101;
				*lbb1 = -2;
			}
			rd2 = ranart_(&rndf77_1.nseed);
			if (rd2 <= .5f)
			{
				*lbb2 = 12;
			}
			else
			{
				*lbb2 = 13;
			}
		}
		else
		{
			if (ifs == 9)
			{
				*iblock = 18092;
				*lbb1 = 1;
			}
			else
			{
				*iblock = 18102;
				*lbb1 = 2;
			}
			rd2 = ranart_(&rndf77_1.nseed);
			if (rd2 <= .5f)
			{
				*lbb2 = -12;
			}
			else
			{
				*lbb2 = -13;
			}
		}
	}
	else if (ifs == 11 || ifs == 12)
	{
		rd = ranart_(&rndf77_1.nseed);
		if (rd <= .5f)
		{
			rd1 = ranart_(&rndf77_1.nseed);
			if (rd1 <= .25f)
			{
				*lbb1 = -6;
			}
			else if (rd1 <= .5f)
			{
				*lbb1 = -7;
			}
			else if (rd1 <= .75f)
			{
				*lbb1 = -8;
			}
			else
			{
				*lbb1 = -9;
			}
			if (ifs == 11)
			{
				*iblock = 18111;
				rd2 = ranart_(&rndf77_1.nseed);
				if (rd2 <= .5f)
				{
					*lbb2 = 10;
				}
				else
				{
					*lbb2 = 11;
				}
			}
			else
			{
				*iblock = 18121;
				rd2 = ranart_(&rndf77_1.nseed);
				if (rd2 <= .5f)
				{
					*lbb2 = 12;
				}
				else
				{
					*lbb2 = 13;
				}
			}
		}
		else
		{
			rd1 = ranart_(&rndf77_1.nseed);
			if (rd1 <= .25f)
			{
				*lbb1 = 6;
			}
			else if (rd1 <= .5f)
			{
				*lbb1 = 7;
			}
			else if (rd1 <= .75f)
			{
				*lbb1 = 8;
			}
			else
			{
				*lbb1 = 9;
			}
			if (ifs == 11)
			{
				*iblock = 18112;
				rd2 = ranart_(&rndf77_1.nseed);
				if (rd2 <= .5f)
				{
					*lbb2 = -10;
				}
				else
				{
					*lbb2 = -11;
				}
			}
			else
			{
				*iblock = 18122;
				rd2 = ranart_(&rndf77_1.nseed);
				if (rd2 <= .5f)
				{
					*lbb2 = -12;
				}
				else
				{
					*lbb2 = -13;
				}
			}
		}
	}
	else if (ifs == 13)
	{
		*iblock = 1813;
		rd1 = ranart_(&rndf77_1.nseed);
		if (rd1 <= .5f)
		{
			*lbb1 = 10;
		}
		else
		{
			*lbb1 = 11;
		}
		rd2 = ranart_(&rndf77_1.nseed);
		if (rd2 <= .5f)
		{
			*lbb2 = -10;
		}
		else
		{
			*lbb2 = -11;
		}
	}
	else if (ifs == 14)
	{
		rd = ranart_(&rndf77_1.nseed);
		if (rd <= .5f)
		{
			*iblock = 18141;
			rd1 = ranart_(&rndf77_1.nseed);
			if (rd1 <= .5f)
			{
				*lbb1 = -10;
			}
			else
			{
				*lbb1 = -11;
			}
			rd2 = ranart_(&rndf77_1.nseed);
			if (rd2 <= .5f)
			{
				*lbb2 = 12;
			}
			else
			{
				*lbb2 = 13;
			}
		}
		else
		{
			*iblock = 18142;
			rd1 = ranart_(&rndf77_1.nseed);
			if (rd1 <= .5f)
			{
				*lbb1 = 10;
			}
			else
			{
				*lbb1 = 11;
			}
			rd2 = ranart_(&rndf77_1.nseed);
			if (rd2 <= .5f)
			{
				*lbb2 = -12;
			}
			else
			{
				*lbb2 = -13;
			}
		}
	}
	else if (ifs == 15)
	{
		*iblock = 1815;
		rd1 = ranart_(&rndf77_1.nseed);
		if (rd1 <= .5f)
		{
			*lbb1 = 12;
		}
		else
		{
			*lbb1 = 13;
		}
		rd2 = ranart_(&rndf77_1.nseed);
		if (rd2 <= .5f)
		{
			*lbb2 = -12;
		}
		else
		{
			*lbb2 = -13;
		}
	}
	else
	{
	}
	return 0;
}

int spprr_(int *lb1, int *lb2, float *srt)
{
	extern double ptor_(float *), rtop_(float *);

	ppmm_1.pprr = 0.f;
	if (*lb1 >= 3 && *lb1 <= 5 && (*lb2 >= 3 && *lb2 <= 5))
	{
		if (*srt > 1.54f)
		{
			ppmm_1.pprr = ptor_(srt);
		}
	}
	else if (*lb1 >= 25 && *lb1 <= 27 && (*lb2 >= 25 && *lb2 <= 27))
	{
		ppmm_1.pprr = rtop_(srt);
	}

	return 0;
}

double ptor_(float *srt)
{
	float ret_val, r__1;

	static float s2;
	extern double rtop_(float *);

	r__1 = *srt;
	s2 = r__1 * r__1;
	ret_val = (s2 - 2.3715999999999999f) * 9 / (s2 - .078400000000000011f) *
			  rtop_(srt);
	return ret_val;
}

double rtop_(float *srt)
{
	float ret_val;

	ret_val = 5.f;
	return ret_val;
}

int pi2ro2_(int *i1, int *i2, int *lbb1, int *lbb2, float *ei1, float *ei2, int *iblock, int *iseed)
{
	extern double ranart_(int *);

	if (ee_1.lb[*i1 - 1] >= 3 && ee_1.lb[*i1 - 1] <= 5 && (ee_1.lb[*i2 - 1] >= 3 && ee_1.lb[*i2 - 1] <= 5))
	{
		*iblock = 1850;
		*ei1 = .77f;
		*ei2 = .77f;
		*lbb1 = (int)(ranart_(&rndf77_1.nseed) * 3) + 25;
		*lbb2 = (int)(ranart_(&rndf77_1.nseed) * 3) + 25;
	}
	else if (ee_1.lb[*i1 - 1] >= 25 && ee_1.lb[*i1 - 1] <= 27 && (ee_1.lb[*i2 - 1] >= 25 && ee_1.lb[*i2 - 1] <= 27))
	{
		*iblock = 1851;
		*lbb1 = (int)(ranart_(&rndf77_1.nseed) * 3) + 3;
		*lbb2 = (int)(ranart_(&rndf77_1.nseed) * 3) + 3;
		*ei1 = .13957f;
		*ei2 = .13957f;
		if (*lbb1 == 4)
		{
			*ei1 = .13496f;
		}
		if (*lbb2 == 4)
		{
			*ei2 = .13496f;
		}
	}
	return 0;
}

int sppee_(int *lb1, int *lb2, float *srt)
{
	extern double ptoe_(float *), etop_(float *);

	ppmm_1.ppee = 0.f;
	if (*lb1 >= 3 && *lb1 <= 5 && (*lb2 >= 3 && *lb2 <= 5))
	{
		if (*srt > 1.095f)
		{
			ppmm_1.ppee = ptoe_(srt);
		}
	}
	else if (*lb1 == 0 && *lb2 == 0)
	{
		ppmm_1.ppee = etop_(srt);
	}
	return 0;
}

double ptoe_(float *srt)
{
	float ret_val, r__1;

	static float s2;
	extern double etop_(float *);

	r__1 = *srt;
	s2 = r__1 * r__1;
	ret_val = (s2 - 1.199025f) * .1111111111111111f / (s2 - .078400000000000011f) * etop_(srt);
	return ret_val;
}

double etop_(float *srt)
{
	float ret_val;

	ret_val = 5.f;
	return ret_val;
}

int pi2et2_(int *i1, int *i2, int *lbb1, int *lbb2, float *ei1, float *ei2, int *iblock, int *iseed)
{
	extern double ranart_(int *);

	if (ee_1.lb[*i1 - 1] >= 3 && ee_1.lb[*i1 - 1] <= 5 && (ee_1.lb[*i2 - 1] >= 3 && ee_1.lb[*i2 - 1] <= 5))
	{
		*iblock = 1860;
		*ei1 = .5475f;
		*ei2 = .5475f;
		*lbb1 = 0;
		*lbb2 = 0;
	}
	else if (ee_1.lb[*i1 - 1] == 0 && ee_1.lb[*i2 - 1] == 0)
	{
		*iblock = 1861;
		*lbb1 = (int)(ranart_(&rndf77_1.nseed) * 3) + 3;
		*lbb2 = (int)(ranart_(&rndf77_1.nseed) * 3) + 3;
		*ei1 = .13957f;
		*ei2 = .13957f;
		if (*lbb1 == 4)
		{
			*ei1 = .13496f;
		}
		if (*lbb2 == 4)
		{
			*ei2 = .13496f;
		}
	}
	return 0;
}

int spppe_(int *lb1, int *lb2, float *srt)
{
	extern double pptope_(float *), petopp_(float *);

	ppmm_1.pppe = 0.f;
	if (*lb1 >= 3 && *lb1 <= 5 && (*lb2 >= 3 && *lb2 <= 5))
	{
		if (*srt > .6875f)
		{
			ppmm_1.pppe = pptope_(srt);
		}
	}
	else if (*lb1 >= 3 && *lb1 <= 5 && *lb2 == 0)
	{
		ppmm_1.pppe = petopp_(srt);
	}
	else if (*lb2 >= 3 && *lb2 <= 5 && *lb1 == 0)
	{
		ppmm_1.pppe = petopp_(srt);
	}
	return 0;
}

double pptope_(float *srt)
{
	float ret_val, r__1;

	double sqrt(double);

	static float s2, pf2, pi2;
	extern double petopp_(float *);

	r__1 = *srt;
	s2 = r__1 * r__1;
	pf2 = (s2 - .47265625f) * (s2 - .16605624999999999f) / 2 / sqrt(s2);
	pi2 = (s2 - .078400000000000011f) * s2 / 2 / sqrt(s2);
	ret_val = pf2 * .33333333333333331f / pi2 * petopp_(srt);
	return ret_val;
}

double petopp_(float *srt)
{
	float ret_val;

	ret_val = 5.f;
	return ret_val;
}

int pi3eta_(int *i1, int *i2, int *lbb1, int *lbb2, float *ei1, float *ei2, int *iblock, int *iseed)
{
	extern double ranart_(int *);

	if (ee_1.lb[*i1 - 1] >= 3 && ee_1.lb[*i1 - 1] <= 5 && (ee_1.lb[*i2 - 1] >= 3 && ee_1.lb[*i2 - 1] <= 5))
	{
		*iblock = 1870;
		*ei1 = .13957f;
		*ei2 = .5475f;
		*lbb1 = (int)(ranart_(&rndf77_1.nseed) * 3) + 3;
		if (*lbb1 == 4)
		{
			*ei1 = .13496f;
		}
		*lbb2 = 0;
	}
	else if (ee_1.lb[*i1 - 1] >= 3 && ee_1.lb[*i1 - 1] <= 5 && ee_1.lb[*i2 - 1] == 0 || ee_1.lb[*i2 - 1] >= 3 && ee_1.lb[*i2 - 1] <= 5 &&
																							ee_1.lb[*i1 - 1] == 0)
	{
		*iblock = 1871;
		*lbb1 = (int)(ranart_(&rndf77_1.nseed) * 3) + 3;
		*lbb2 = (int)(ranart_(&rndf77_1.nseed) * 3) + 3;
		*ei1 = .13957f;
		*ei2 = .13957f;
		if (*lbb1 == 4)
		{
			*ei1 = .13496f;
		}
		if (*lbb2 == 4)
		{
			*ei2 = .13496f;
		}
	}
	return 0;
}

int srpre_(int *lb1, int *lb2, float *srt)
{
	extern double rptore_(float *), retorp_(float *);

	ppmm_1.rpre = 0.f;
	if (*lb1 >= 25 && *lb1 <= 27 && *lb2 >= 3 && *lb2 <= 5)
	{
		if (*srt > 1.3174999999999999f)
		{
			ppmm_1.rpre = rptore_(srt);
		}
	}
	else if (*lb2 >= 25 && *lb2 <= 27 && *lb1 >= 3 && *lb1 <= 5)
	{
		if (*srt > 1.3174999999999999f)
		{
			ppmm_1.rpre = rptore_(srt);
		}
	}
	else if (*lb1 >= 25 && *lb1 <= 27 && *lb2 == 0)
	{
		if (*srt > .91000000000000003f)
		{
			ppmm_1.rpre = retorp_(srt);
		}
	}
	else if (*lb2 >= 25 && *lb2 <= 27 && *lb1 == 0)
	{
		if (*srt > .91000000000000003f)
		{
			ppmm_1.rpre = retorp_(srt);
		}
	}
	return 0;
}

double rptore_(float *srt)
{
	float ret_val, r__1;

	double sqrt(double);

	static float s2, pf2, pi2;
	extern double retorp_(float *);

	r__1 = *srt;
	s2 = r__1 * r__1;
	pf2 = (s2 - 1.7358062499999998f) * (s2 - .049506250000000016f) / 2 / sqrt(s2);
	pi2 = (s2 - .82810000000000006f) * (s2 - .39690000000000003f) / 2 / sqrt(s2);
	ret_val = pf2 * .33333333333333331f / pi2 * retorp_(srt);
	return ret_val;
}

double retorp_(float *srt)
{
	float ret_val;

	ret_val = 5.f;
	return ret_val;
}

int rpiret_(int *i1, int *i2, int *lbb1, int *lbb2, float *ei1, float *ei2, int *iblock, int *iseed)
{
	extern double ranart_(int *);

	if (ee_1.lb[*i1 - 1] >= 25 && ee_1.lb[*i1 - 1] <= 27 && ee_1.lb[*i2 - 1] >= 3 && ee_1.lb[*i2 - 1] <= 5 || ee_1.lb[*i1 - 1] >= 3 && ee_1.lb[*i1 - 1] <= 5 && ee_1.lb[*i2 - 1] >= 25 && ee_1.lb[*i2 - 1] <= 27)
	{
		*iblock = 1880;
		*ei1 = .77f;
		*ei2 = .5475f;
		*lbb1 = (int)(ranart_(&rndf77_1.nseed) * 3) + 25;
		*lbb2 = 0;
	}
	else if (ee_1.lb[*i1 - 1] >= 25 && ee_1.lb[*i1 - 1] <= 27 && ee_1.lb[*i2 - 1] == 0 || ee_1.lb[*i2 - 1] >= 25 && ee_1.lb[*i2 - 1] <= 27 && ee_1.lb[*i1 - 1] == 0)
	{
		*iblock = 1881;
		*lbb1 = (int)(ranart_(&rndf77_1.nseed) * 3) + 25;
		*lbb2 = (int)(ranart_(&rndf77_1.nseed) * 3) + 3;
		*ei1 = .77f;
		*ei2 = .13957f;
		if (*lbb2 == 4)
		{
			*ei2 = .13496f;
		}
	}
	return 0;
}

int sopoe_(int *lb1, int *lb2, float *srt)
{
	extern double xop2oe_(float *), xoe2op_(float *);

	ppmm_1.xopoe = 0.f;
	if (*lb1 == 28 && *lb2 >= 3 && *lb2 <= 5 || *lb2 == 28 && *lb1 >= 3 && *lb1 <= 5)
	{
		if (*srt > 1.3294999999999999f)
		{
			ppmm_1.xopoe = xop2oe_(srt);
		}
	}
	else if (*lb1 == 28 && *lb2 == 0 || *lb1 == 0 && *lb2 == 28)
	{
		if (*srt > 1.3294999999999999f)
		{
			ppmm_1.xopoe = xoe2op_(srt);
		}
	}
	return 0;
}

double xop2oe_(float *srt)
{
	float ret_val, r__1;

	double sqrt(double);

	static float s2, pf2, pi2;
	extern double xoe2op_(float *);

	r__1 = *srt;
	s2 = r__1 * r__1;
	pf2 = (s2 - 1.7675702499999997f) * (s2 - .054990250000000018f) / 2 / sqrt(s2);
	pi2 = (s2 - .85008400000000006f) * (s2 - .41216400000000003f) / 2 / sqrt(s2);
	ret_val = pf2 * .33333333333333331f / pi2 * xoe2op_(srt);
	return ret_val;
}

double xoe2op_(float *srt)
{
	float ret_val;

	ret_val = 5.f;
	return ret_val;
}

int opioet_(int *i1, int *i2, int *lbb1, int *lbb2, float *ei1, float *ei2, int *iblock, int *iseed)
{
	extern double ranart_(int *);

	if (ee_1.lb[*i1 - 1] >= 3 && ee_1.lb[*i1 - 1] <= 5 && ee_1.lb[*i2 - 1] == 28 || ee_1.lb[*i2 - 1] >= 3 && ee_1.lb[*i2 - 1] <= 5 && ee_1.lb[*i1 - 1] == 28)
	{
		*iblock = 1890;
		*ei1 = .782f;
		*ei2 = .5475f;
		*lbb1 = 28;
		*lbb2 = 0;
	}
	else if (ee_1.lb[*i1 - 1] == 28 && ee_1.lb[*i2 - 1] == 0 || ee_1.lb[*i1 - 1] == 0 && ee_1.lb[*i2 - 1] == 28)
	{
		*iblock = 1891;
		*lbb1 = 28;
		*lbb2 = (int)(ranart_(&rndf77_1.nseed) * 3) + 3;
		*ei1 = .782f;
		*ei2 = .13957f;
		if (*lbb2 == 4)
		{
			*ei2 = .13496f;
		}
	}
	return 0;
}

int srree_(int *lb1, int *lb2, float *srt)
{
	extern double rrtoee_(float *), eetorr_(float *);

	ppmm_1.rree = 0.f;
	if (*lb1 >= 25 && *lb1 <= 27 && *lb2 >= 25 && *lb2 <= 27)
	{
		if (*srt > 1.095f)
		{
			ppmm_1.rree = rrtoee_(srt);
		}
	}
	else if (*lb1 == 0 && *lb2 == 0)
	{
		if (*srt > 1.54f)
		{
			ppmm_1.rree = eetorr_(srt);
		}
	}
	return 0;
}

double eetorr_(float *srt)
{
	float ret_val, r__1;

	static float s2;
	extern double rrtoee_(float *);

	r__1 = *srt;
	s2 = r__1 * r__1;
	ret_val = (s2 - 2.3715999999999999f) * 81.f / (s2 - 1.199025f) * rrtoee_(srt);
	return ret_val;
}

double rrtoee_(float *srt)
{
	float ret_val;

	ret_val = 5.f;
	return ret_val;
}

int ro2et2_(int *i1, int *i2, int *lbb1, int *lbb2, float *ei1, float *ei2, int *iblock, int *iseed)
{
	extern double ranart_(int *);

	if (ee_1.lb[*i1 - 1] >= 25 && ee_1.lb[*i1 - 1] <= 27 && ee_1.lb[*i2 - 1] >= 25 && ee_1.lb[*i2 - 1] <= 27)
	{
		*iblock = 1895;
		*ei1 = .5475f;
		*ei2 = .5475f;
		*lbb1 = 0;
		*lbb2 = 0;
	}
	else if (ee_1.lb[*i1 - 1] == 0 && ee_1.lb[*i2 - 1] == 0)
	{
		*iblock = 1896;
		*lbb1 = (int)(ranart_(&rndf77_1.nseed) * 3) + 25;
		*lbb2 = (int)(ranart_(&rndf77_1.nseed) * 3) + 25;
		*ei1 = .77f;
		*ei2 = .77f;
	}
	return 0;
}

int xkksan_(int *i1, int *i2, float *srt, float *sigks1, float *sigks2, float *sigks3, float *sigks4, float *sigk, float *prkk)
{
	float r__1, r__2;

	static float s, pf2, pi2, xm1, xm2, xpion0;

	r__1 = *srt;
	s = r__1 * r__1;
	*sigks1 = 1e-8f;
	*sigks2 = 1e-8f;
	*sigks3 = 1e-8f;
	*sigks4 = 1e-8f;
	xpion0 = *prkk;
	xpion0 /= 2;
	r__1 = cc_1.e[*i1 - 1] + cc_1.e[*i2 - 1];
	r__2 = cc_1.e[*i1 - 1] - cc_1.e[*i2 - 1];
	pi2 = (s - r__1 * r__1) * (s - r__2 * r__2);
	*sigk = 1e-8f;
	if (pi2 <= 0.f)
	{
		return 0;
	}
	xm1 = .14f;
	xm2 = .77f;
	r__1 = xm1 + xm2;
	r__2 = xm1 - xm2;
	pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
	if (pi2 > 0.f && pf2 > 0.f)
	{
		*sigks1 = pf2 * 6.75f / pi2 * xpion0;
	}
	xm1 = .14f;
	xm2 = .7819f;
	r__1 = xm1 + xm2;
	r__2 = xm1 - xm2;
	pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
	if (pi2 > 0.f && pf2 > 0.f)
	{
		*sigks2 = pf2 * 2.25f / pi2 * xpion0;
	}
	xm1 = .77f;
	xm2 = .5473f;
	r__1 = xm1 + xm2;
	r__2 = xm1 - xm2;
	pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
	if (pf2 > 0.f)
	{
		*sigks3 = pf2 * 2.25f / pi2 * xpion0;
	}
	xm1 = .7819f;
	xm2 = .5473f;
	r__1 = xm1 + xm2;
	r__2 = xm1 - xm2;
	pf2 = (s - r__1 * r__1) * (s - r__2 * r__2);
	if (pf2 > 0.f)
	{
		*sigks4 = pf2 * .75f / pi2 * xpion0;
	}
	*sigk = *sigks1 + *sigks2 + *sigks3 + *sigks4;
	return 0;
}

int crkspi_(int *i1, int *i2, float *xsk1, float *xsk2,
			float *xsk3, float *xsk4, float *sigk, int *iblock, int *lbp1,
			int *lbp2, float *emm1, float *emm2)
{
	static float x1;
	extern double ranart_(int *);

	*iblock = 466;
	x1 = ranart_(&rndf77_1.nseed) * *sigk;
	*xsk2 = *xsk1 + *xsk2;
	*xsk3 = *xsk2 + *xsk3;
	*xsk4 = *xsk3 + *xsk4;
	if (x1 <= *xsk1)
	{
		ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 3) + 3;
		ee_1.lb[*i2 - 1] = (int)(ranart_(&rndf77_1.nseed) * 3) + 25;
		cc_1.e[*i1 - 1] = .13957f;
		cc_1.e[*i2 - 1] = .77f;
	}
	else if (x1 <= *xsk2)
	{
		ee_1.lb[*i1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 3) + 3;
		ee_1.lb[*i2 - 1] = 28;
		cc_1.e[*i1 - 1] = .13957f;
		cc_1.e[*i2 - 1] = .782f;
	}
	else if (x1 <= *xsk3)
	{
		ee_1.lb[*i1 - 1] = 0;
		ee_1.lb[*i2 - 1] = (int)(ranart_(&rndf77_1.nseed) * 3) + 25;
		cc_1.e[*i1 - 1] = .548f;
		cc_1.e[*i2 - 1] = .77f;
	}
	else
	{
		ee_1.lb[*i1 - 1] = 0;
		ee_1.lb[*i2 - 1] = 28;
		cc_1.e[*i1 - 1] = .548f;
		cc_1.e[*i2 - 1] = .782f;
	}
	if (ee_1.lb[*i1 - 1] == 4)
	{
		cc_1.e[*i1 - 1] = .13496f;
	}
	*lbp1 = ee_1.lb[*i1 - 1];
	*lbp2 = ee_1.lb[*i2 - 1];
	*emm1 = cc_1.e[*i1 - 1];
	*emm2 = cc_1.e[*i2 - 1];
	return 0;
}

int ksreso_(int *i1, int *i2)
{
	float r__1, r__2, r__3, r__4;

	double sqrt(double);

	static int i__;
	static float e10, e20, dm;

	r__1 = cc_1.e[*i1 - 1];
	r__2 = bb_1.p[*i1 * 3 - 3];
	r__3 = bb_1.p[*i1 * 3 - 2];
	r__4 = bb_1.p[*i1 * 3 - 1];
	e10 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	r__1 = cc_1.e[*i2 - 1];
	r__2 = bb_1.p[*i2 * 3 - 3];
	r__3 = bb_1.p[*i2 * 3 - 2];
	r__4 = bb_1.p[*i2 * 3 - 1];
	e20 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	if (ee_1.lb[*i2 - 1] == 21 || ee_1.lb[*i2 - 1] == 23)
	{
		cc_1.e[*i1 - 1] = 0.f;
		i__ = *i2;
	}
	else
	{
		cc_1.e[*i2 - 1] = 0.f;
		i__ = *i1;
	}
	if (ee_1.lb[i__ - 1] == 23)
	{
		ee_1.lb[i__ - 1] = 30;
	}
	else if (ee_1.lb[i__ - 1] == 21)
	{
		ee_1.lb[i__ - 1] = -30;
	}
	bb_1.p[i__ * 3 - 3] = bb_1.p[*i1 * 3 - 3] + bb_1.p[*i2 * 3 - 3];
	bb_1.p[i__ * 3 - 2] = bb_1.p[*i1 * 3 - 2] + bb_1.p[*i2 * 3 - 2];
	bb_1.p[i__ * 3 - 1] = bb_1.p[*i1 * 3 - 1] + bb_1.p[*i2 * 3 - 1];
	r__1 = e10 + e20;
	r__2 = bb_1.p[i__ * 3 - 3];
	r__3 = bb_1.p[i__ * 3 - 2];
	r__4 = bb_1.p[i__ * 3 - 1];
	dm = sqrt(r__1 * r__1 - r__2 * r__2 - r__3 * r__3 - r__4 * r__4);
	cc_1.e[i__ - 1] = dm;
	return 0;
}

int pertur_(float *px, float *py, float *pz, float *srt, int *irun, int *i1, int *i2, int *nt, int *kp, int *icont)
{
	float r__1, r__2, r__3, r__4, r__5, r__6, r__7, r__8;

	double sqrt(double), cos(double), sin(double);

	static float c1, s1, t1, x1, y1, z1, x2, y2, z2, ec;
	static int ic;
	static float ds, pr;
	static int lb1, lb2;
	static float em1, em2, ct1, pr2, px0, py0, pz0, st1;
	static int idn, idp;
	static float app, pii, pff, sig, dfr, dsr, xpt, ypt, zpt, e1cm, e2cm, acap,
		akal, akap, alas, asap, cmat, ames, aomp, brpp, ppt11, ppt12,
		ppt13, ppt21, ppt22, ppt23, srrt;
	static int lbpp1, lbpp2;
	static float empp1, prob1, prob2, empp2, sigca, pkaon, sigpe, sigpi, xrand,
		sigom, p1beta;
	static int icsbel;
	static float sigcal, sigcas, sigeta;
	extern int distce_(int *, int *, float *, float *,
					   float *, float *, float *, int *, float *, float *, float *);
	extern double aknpsg_(float *), ranart_(int *);
	static float sigomm, transf;
	extern int rotate_(float *, float *, float *, float *, float *, float *);

	px0 = *px;
	py0 = *py;
	pz0 = *pz;
	lb1 = ee_1.lb[*i1 - 1];
	em1 = cc_1.e[*i1 - 1];
	x1 = aa_1.r__[*i1 * 3 - 3];
	y1 = aa_1.r__[*i1 * 3 - 2];
	z1 = aa_1.r__[*i1 * 3 - 1];
	prob1 = hh_1.proper[*i1 - 1];

	lb2 = ee_1.lb[*i2 - 1];
	em2 = cc_1.e[*i2 - 1];
	x2 = aa_1.r__[*i2 * 3 - 3];
	y2 = aa_1.r__[*i2 * 3 - 2];
	z2 = aa_1.r__[*i2 * 3 - 1];
	prob2 = hh_1.proper[*i2 - 1];

	*icont = 1;
	icsbel = -1;
	if ((lb1 == 21 || lb1 == 23 || abs(lb1) == 30) && (abs(lb2) >= 14 && abs(
																			 lb2) <= 17))
	{
		goto L60;
	}
	if ((lb2 == 21 || lb2 == 23 || abs(lb2) == 30) && (abs(lb1) >= 14 && abs(
																			 lb1) <= 17))
	{
		goto L60;
	}
	if ((lb1 == 21 || lb1 == 23 || abs(lb1) == 30) && (abs(lb2) == 40 || abs(
																			 lb2) == 41))
	{
		goto L70;
	}
	if ((lb2 == 21 || lb2 == 23 || abs(lb2) == 30) && (abs(lb1) == 40 || abs(
																			 lb1) == 41))
	{
		goto L70;
	}

	if ((lb1 >= 3 && lb1 <= 5 || lb1 == 0) && (abs(lb2) == 40 || abs(lb2) ==
																	 41) ||
		(lb2 >= 3 && lb2 <= 5 || lb2 == 0) && (abs(lb1) == 40 ||
											   abs(lb1) == 41))
	{
		goto L90;
	}
	if (lb1 >= 3 && lb1 <= 5 && abs(lb2) == 45 || lb2 >= 3 && lb2 <= 5 && abs(lb1) == 45)
	{
		goto L110;
	}

L60:
	if (abs(lb1) >= 14 && abs(lb1) <= 17)
	{
		asap = cc_1.e[*i1 - 1];
		akap = cc_1.e[*i2 - 1];
		idp = *i1;
	}
	else
	{
		asap = cc_1.e[*i2 - 1];
		akap = cc_1.e[*i1 - 1];
		idp = *i2;
	}
	app = .138f;
	if (*srt < app + 1.3213f)
	{
		return 0;
	}
	srrt = *srt - (app + 1.3213f) + (akap + .939457f);
	r__2 = srrt;
	r__3 = akap;
	r__1 = (r__2 * r__2 - (r__3 * r__3 + .88257945484900002f)) / 2.f /
		   .939457f;
	r__4 = akap;
	pkaon = sqrt(r__1 * r__1 - r__4 * r__4);
	sigca = (aknpsg_(&pkaon) + aknpsg_(&pkaon)) * 1.5f;
	r__1 = *srt;
	r__2 = akap + .939457f;
	r__3 = *srt;
	r__4 = .939457f - akap;
	pii = sqrt((r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4));
	r__1 = *srt;
	r__2 = asap + app;
	r__3 = *srt;
	r__4 = asap - app;
	pff = sqrt((r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4));
	cmat = sigca * pii / pff;
	r__1 = *srt;
	r__2 = app + 1.3213f;
	r__3 = *srt;
	r__4 = 1.3213f - app;
	r__5 = *srt;
	r__6 = asap + akap;
	r__7 = *srt;
	r__8 = asap - akap;
	sigpi = cmat * sqrt((r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4)) / sqrt((r__5 * r__5 - r__6 * r__6) * (r__7 * r__7 - r__8 * r__8));

	sigeta = 0.f;
	if (*srt > 1.8693f)
	{
		srrt = *srt - 1.8693f + (akap + .939457f);
		r__2 = srrt;
		r__3 = akap;
		r__1 = (r__2 * r__2 - (r__3 * r__3 + .88257945484900002f)) / 2.f /
			   .939457f;
		r__4 = akap;
		pkaon = sqrt(r__1 * r__1 - r__4 * r__4);
		sigca = (aknpsg_(&pkaon) + aknpsg_(&pkaon)) * 1.5f;
		cmat = sigca * pii / pff;
		r__1 = *srt;
		r__2 = *srt;
		r__3 = *srt;
		r__4 = asap + akap;
		r__5 = *srt;
		r__6 = asap - akap;
		sigeta = cmat * sqrt((r__1 * r__1 - 3.4942824899999998f) * (r__2 * r__2 - .59799288999999978f)) / sqrt((r__3 * r__3 - r__4 * r__4) * (r__5 * r__5 - r__6 * r__6));
	}

	sigca = sigpi + sigeta;
	sigpe = 0.f;
	sig = dmax(sigpe, sigca);
	ds = sqrt(sig / 31.4f);
	dsr = ds + .1f;
	r__1 = em1 + em2 + .02f;
	ec = r__1 * r__1;
	distce_(i1, i2, &dsr, &ds, &input1_1.dt, &ec, srt, &ic, px, py, pz);
	if (ic == -1)
	{
		return 0;
	}
	brpp = sigca / sig;

	if (lb1 >= 14 && lb1 <= 17 || lb2 >= 14 && lb2 <= 17)
	{
		lbpp1 = (int)(ranart_(&rndf77_1.nseed) * 2) + 40;
	}
	else
	{
		lbpp1 = -40 - (int)(ranart_(&rndf77_1.nseed) * 2);
	}
	empp1 = 1.3213f;
	if (ranart_(&rndf77_1.nseed) < sigpi / sigca)
	{
		lbpp2 = (int)(ranart_(&rndf77_1.nseed) * 3) + 3;
		empp2 = .138f;
	}
	else
	{
		lbpp2 = 0;
		empp2 = .548f;
	}
	if (ranart_(&rndf77_1.nseed) < brpp)
	{
		*icont = 0;
		ee_1.lb[*i1 - 1] = lbpp1;
		cc_1.e[*i1 - 1] = empp1;
		hh_1.proper[*i1 - 1] = brpp;
		ee_1.lb[*i2 - 1] = lbpp2;
		cc_1.e[*i2 - 1] = empp2;
		hh_1.proper[*i2 - 1] = 1.f;
	}
	goto L700;
L70:
	if (abs(lb1) == 40 || abs(lb1) == 41)
	{
		acap = cc_1.e[*i1 - 1];
		akap = cc_1.e[*i2 - 1];
		idp = *i1;
	}
	else
	{
		acap = cc_1.e[*i2 - 1];
		akap = cc_1.e[*i1 - 1];
		idp = *i2;
	}
	app = .138f;
	ames = .138f;
	if (*srt < ames + 1.6724f)
	{
		return 0;
	}
	srrt = *srt - (ames + 1.6724f) + (akap + .939457f);
	r__2 = srrt;
	r__3 = akap;
	r__1 = (r__2 * r__2 - (r__3 * r__3 + .88257945484900002f)) / 2.f /
		   .939457f;
	r__4 = akap;
	pkaon = sqrt(r__1 * r__1 - r__4 * r__4);
	sigomm = (aknpsg_(&pkaon) + aknpsg_(&pkaon)) * 1.5f;
	r__1 = *srt;
	r__2 = akap + .939457f;
	r__3 = *srt;
	r__4 = .939457f - akap;
	r__5 = *srt;
	r__6 = app + 1.1974f;
	r__7 = *srt;
	r__8 = 1.1974f - app;
	cmat = sigomm * sqrt((r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4)) / sqrt((r__5 * r__5 - r__6 * r__6) * (r__7 * r__7 - r__8 * r__8));
	r__1 = *srt;
	r__2 = ames + 1.6724f;
	r__3 = *srt;
	r__4 = 1.6724f - ames;
	r__5 = *srt;
	r__6 = acap + akap;
	r__7 = *srt;
	r__8 = acap - akap;
	sigom = cmat * sqrt((r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4)) / sqrt((r__5 * r__5 - r__6 * r__6) * (r__7 * r__7 - r__8 * r__8));
	sigpe = 0.f;
	sig = dmax(sigpe, sigom);
	ds = sqrt(sig / 31.4f);
	dsr = ds + .1f;
	r__1 = em1 + em2 + .02f;
	ec = r__1 * r__1;
	distce_(i1, i2, &dsr, &ds, &input1_1.dt, &ec, srt, &ic, px, py, pz);
	if (ic == -1)
	{
		return 0;
	}
	brpp = sigom / sig;

	if (lb1 >= 40 && lb1 <= 41 || lb2 >= 40 && lb2 <= 41)
	{
		lbpp1 = 45;
	}
	else
	{
		lbpp1 = -45;
	}
	empp1 = 1.6724f;
	lbpp2 = (int)(ranart_(&rndf77_1.nseed) * 3) + 3;
	empp2 = ames;

	xrand = ranart_(&rndf77_1.nseed);
	if (xrand < hh_1.proper[idp - 1] * brpp)
	{
		*icont = 0;
		ee_1.lb[*i1 - 1] = lbpp1;
		cc_1.e[*i1 - 1] = empp1;
		hh_1.proper[*i1 - 1] = hh_1.proper[idp - 1] * brpp;
		ee_1.lb[*i2 - 1] = lbpp2;
		cc_1.e[*i2 - 1] = empp2;
		hh_1.proper[*i2 - 1] = 1.f;
	}
	else if (xrand < brpp)
	{
		cc_1.e[idp - 1] = 0.f;
	}
	goto L700;
L90:
	if (abs(lb1) == 40 || abs(lb1) == 41)
	{
		acap = cc_1.e[*i1 - 1];
		app = cc_1.e[*i2 - 1];
		idp = *i1;
		idn = *i2;
	}
	else
	{
		acap = cc_1.e[*i2 - 1];
		app = cc_1.e[*i1 - 1];
		idp = *i2;
		idn = *i1;
	}
	akal = .498f;

	alas = 1.1157f;
	if (*srt <= alas + .498f)
	{
		return 0;
	}
	srrt = *srt - (acap + app) + 1.437457f;
	r__2 = srrt;
	r__1 = (r__2 * r__2 - 1.1305834548489999f) / 2.f / .939457f;
	pkaon = sqrt(r__1 * r__1 - .248004f);
	sigca = (aknpsg_(&pkaon) + aknpsg_(&pkaon)) * 1.5f;
	r__1 = *srt;
	r__2 = *srt;
	r__3 = *srt;
	r__4 = alas + .138f;
	r__5 = *srt;
	r__6 = alas - .138f;
	cmat = sigca * sqrt((r__1 * r__1 - 2.066282626849f) * (r__2 * r__2 - .194884282849f)) / sqrt((r__3 * r__3 - r__4 * r__4) * (r__5 * r__5 - r__6 * r__6));
	r__1 = *srt;
	r__2 = acap + app;
	r__3 = *srt;
	r__4 = acap - app;
	r__5 = *srt;
	r__6 = alas + .498f;
	r__7 = *srt;
	r__8 = alas - .498f;
	sigca = cmat * sqrt((r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4)) / sqrt((r__5 * r__5 - r__6 * r__6) * (r__7 * r__7 - r__8 * r__8));
	dfr = .33333333333333331f;
	if (ee_1.lb[idn - 1] == 0)
	{
		dfr = 1.f;
	}
	r__1 = *srt;
	r__2 = alas + .498f;
	r__3 = *srt;
	r__4 = alas - .498f;
	r__5 = *srt;
	r__6 = acap + app;
	r__7 = *srt;
	r__8 = acap - app;
	sigcal = sigca * dfr * (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / (r__5 * r__5 - r__6 * r__6) / (r__7 * r__7 - r__8 * r__8);

	alas = 1.1974f;
	if (*srt <= alas + .498f)
	{
		sigcas = 0.f;
	}
	else
	{
		srrt = *srt - (acap + app) + 1.437457f;
		r__2 = srrt;
		r__1 = (r__2 * r__2 - 1.1305834548489999f) / 2.f / .939457f;
		pkaon = sqrt(r__1 * r__1 - .248004f);
		sigca = (aknpsg_(&pkaon) + aknpsg_(&pkaon)) * 1.5f;
		r__1 = *srt;
		r__2 = *srt;
		r__3 = *srt;
		r__4 = alas + .138f;
		r__5 = *srt;
		r__6 = alas - .138f;
		cmat = sigca * sqrt((r__1 * r__1 - 2.066282626849f) * (r__2 * r__2 - .194884282849f)) / sqrt((r__3 * r__3 - r__4 * r__4) * (r__5 * r__5 - r__6 * r__6));
		r__1 = *srt;
		r__2 = acap + app;
		r__3 = *srt;
		r__4 = acap - app;
		r__5 = *srt;
		r__6 = alas + .498f;
		r__7 = *srt;
		r__8 = alas - .498f;
		sigca = cmat * sqrt((r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4)) / sqrt((r__5 * r__5 - r__6 * r__6) * (r__7 * r__7 - r__8 * r__8));
		dfr = 1.f;
		if (ee_1.lb[idn - 1] == 0)
		{
			dfr = 3.f;
		}
		r__1 = *srt;
		r__2 = alas + .498f;
		r__3 = *srt;
		r__4 = alas - .498f;
		r__5 = *srt;
		r__6 = acap + app;
		r__7 = *srt;
		r__8 = acap - app;
		sigcas = sigca * dfr * (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) / (r__5 * r__5 - r__6 * r__6) / (r__7 * r__7 - r__8 * r__8);
	}

	sig = sigcal + sigcas;
	brpp = 1.f;
	ds = sqrt(sig / 31.4f);
	dsr = ds + .1f;
	r__1 = em1 + em2 + .02f;
	ec = r__1 * r__1;
	distce_(i1, i2, &dsr, &ds, &input1_1.dt, &ec, srt, &ic, px, py, pz);

	if (ic == -1)
	{
		ds = sqrt(.63694267515923575f);
		dsr = ds + .1f;
		distce_(i1, i2, &dsr, &ds, &input1_1.dt, &ec, srt, &icsbel, px, py,
				pz);
		if (icsbel == -1)
		{
			return 0;
		}
		empp1 = em1;
		empp2 = em2;
		goto L700;
	}

	if (sigcal / sig > ranart_(&rndf77_1.nseed))
	{
		if (lb1 == 40 || lb1 == 41 || lb2 == 40 || lb2 == 41)
		{
			lbpp1 = 21;
			lbpp2 = 14;
		}
		else
		{
			lbpp1 = 23;
			lbpp2 = -14;
		}
		alas = 1.1157f;
	}
	else
	{
		if (lb1 == 40 || lb1 == 41 || lb2 == 40 || lb2 == 41)
		{
			lbpp1 = 21;
			lbpp2 = (int)(ranart_(&rndf77_1.nseed) * 3) + 15;
		}
		else
		{
			lbpp1 = 23;
			lbpp2 = -15 - (int)(ranart_(&rndf77_1.nseed) * 3);
		}
		alas = 1.1974f;
	}
	empp1 = .498f;
	empp2 = alas;

	if (ranart_(&rndf77_1.nseed) < hh_1.proper[idp - 1])
	{
		*icont = 0;
		ee_1.lb[*i1 - 1] = lbpp1;
		cc_1.e[*i1 - 1] = empp1;
		hh_1.proper[*i1 - 1] = 1.f;
		ee_1.lb[*i2 - 1] = lbpp2;
		cc_1.e[*i2 - 1] = empp2;
		hh_1.proper[*i2 - 1] = 1.f;
		goto L700;
	}
	else
	{
		cc_1.e[idp - 1] = 0.f;
	}
	return 0;

L110:
	if (lb1 == 45 || lb1 == -45)
	{
		aomp = cc_1.e[*i1 - 1];
		app = cc_1.e[*i2 - 1];
		idp = *i1;
		idn = *i2;
	}
	else
	{
		aomp = cc_1.e[*i2 - 1];
		app = cc_1.e[*i1 - 1];
		idp = *i2;
		idn = *i1;
	}
	akal = .498f;
	if (*srt <= 1.8192999999999999f)
	{
		return 0;
	}
	srrt = *srt - (app + 1.6724f) + 1.437457f;
	r__2 = srrt;
	r__1 = (r__2 * r__2 - 1.1305834548489999f) / 2.f / .939457f;
	pkaon = sqrt(r__1 * r__1 - .248004f);
	sigca = (aknpsg_(&pkaon) + aknpsg_(&pkaon)) * 1.5f;
	r__1 = *srt;
	r__2 = *srt;
	r__3 = *srt;
	r__4 = *srt;
	cmat = sigca * sqrt((r__1 * r__1 - 2.066282626849f) * (r__2 * r__2 - .194884282849f)) / sqrt((r__3 * r__3 - 1.7832931599999997f) * (r__4 * r__4 - 1.1223283600000002f));
	r__1 = *srt;
	r__2 = aomp + app;
	r__3 = *srt;
	r__4 = aomp - app;
	r__5 = *srt;
	r__6 = *srt;
	sigom = cmat * sqrt((r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4)) / sqrt((r__5 * r__5 - 3.3098524899999995f) * (r__6 * r__6 - .67782288999999984f));
	dfr = .66666666666666663f;
	r__1 = *srt;
	r__2 = *srt;
	r__3 = *srt;
	r__4 = aomp + app;
	r__5 = *srt;
	r__6 = aomp - app;
	sigom = sigom * dfr * (r__1 * r__1 - 3.3098524899999995f) * (r__2 * r__2 - .67782288999999984f) / (r__3 * r__3 - r__4 * r__4) / (r__5 * r__5 - r__6 * r__6);

	brpp = 1.f;
	ds = sqrt(sigom / 31.4f);
	dsr = ds + .1f;
	r__1 = em1 + em2 + .02f;
	ec = r__1 * r__1;
	distce_(i1, i2, &dsr, &ds, &input1_1.dt, &ec, srt, &ic, px, py, pz);

	if (ic == -1)
	{
		ds = sqrt(.63694267515923575f);
		dsr = ds + .1f;
		distce_(i1, i2, &dsr, &ds, &input1_1.dt, &ec, srt, &icsbel, px, py,
				pz);
		if (icsbel == -1)
		{
			return 0;
		}
		empp1 = em1;
		empp2 = em2;
		goto L700;
	}

	if (lb1 == 45 || lb2 == 45)
	{
		lbpp1 = (int)(ranart_(&rndf77_1.nseed) * 2) + 40;
		lbpp2 = 21;
	}
	else
	{
		lbpp1 = -40 - (int)(ranart_(&rndf77_1.nseed) * 2);
		lbpp2 = 23;
	}
	empp1 = 1.3213f;
	empp2 = .498f;

	if (ranart_(&rndf77_1.nseed) < hh_1.proper[idp - 1])
	{
		*icont = 0;
		ee_1.lb[*i1 - 1] = lbpp1;
		cc_1.e[*i1 - 1] = empp1;
		hh_1.proper[*i1 - 1] = hh_1.proper[idp - 1];
		ee_1.lb[*i2 - 1] = lbpp2;
		cc_1.e[*i2 - 1] = empp2;
		hh_1.proper[*i2 - 1] = 1.f;
	}
	else
	{
		cc_1.e[idp - 1] = 0.f;
	}
	goto L700;

L700:
	r__2 = *srt;
	r__3 = empp1;
	r__4 = empp2;
	r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
	r__5 = empp1 * empp2;
	pr2 = r__1 * r__1 - r__5 * r__5 * 4.f;
	if (pr2 <= 0.f)
	{
		pr2 = 1e-8f;
	}
	pr = sqrt(pr2) / (*srt * 2.f);
	c1 = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
	t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
	r__1 = c1;
	s1 = sqrt(1.f - r__1 * r__1);
	ct1 = cos(t1);
	st1 = sin(t1);
	*pz = pr * c1;
	*px = pr * s1 * ct1;
	*py = pr * s1 * st1;
	rotate_(&px0, &py0, &pz0, px, py, pz);
	if (*icont == 0)
	{
		return 0;
	}

	r__1 = empp1;
	r__2 = *px;
	r__3 = *py;
	r__4 = *pz;
	e1cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	p1beta = *px * bg_1.betax + *py * bg_1.betay + *pz * bg_1.betaz;
	transf = bg_1.gamma * (bg_1.gamma * p1beta / (bg_1.gamma + 1) + e1cm);
	ppt11 = bg_1.betax * transf + *px;
	ppt12 = bg_1.betay * transf + *py;
	ppt13 = bg_1.betaz * transf + *pz;

	if (icsbel != -1)
	{
		bb_1.p[*i1 * 3 - 3] = ppt11;
		bb_1.p[*i1 * 3 - 2] = ppt12;
		bb_1.p[*i1 * 3 - 1] = ppt13;
		r__1 = empp2;
		r__2 = *px;
		r__3 = *py;
		r__4 = *pz;
		e2cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
		transf = bg_1.gamma * (-bg_1.gamma * p1beta / (bg_1.gamma + 1) + e2cm);
		ppt21 = bg_1.betax * transf - *px;
		ppt22 = bg_1.betay * transf - *py;
		ppt23 = bg_1.betaz * transf - *pz;
		bb_1.p[*i2 * 3 - 3] = ppt21;
		bb_1.p[*i2 * 3 - 2] = ppt22;
		bb_1.p[*i2 * 3 - 1] = ppt23;
		return 0;
	}
	xpt = x1;
	ypt = y1;
	zpt = z1;

	++nn_1.nnn;
	pe_1.propi[nn_1.nnn + *irun * 150001 - 150002] = hh_1.proper[idp - 1] *
													 brpp;
	pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = lbpp1;
	pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = empp1;
	pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = xpt;
	pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = ypt;
	pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = zpt;
	pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = ppt11;
	pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = ppt12;
	pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = ppt13;
	dpert_1.dppion[nn_1.nnn + *irun * 150001 - 150002] = dpert_1.dpertp[*i1 -
																		1] *
														 dpert_1.dpertp[*i2 - 1];
	return 0;
}

int crhb_(float *px, float *py, float *pz, float *srt, int *i1, int *i2, int *iblock)
{
	float r__1, r__2, r__3, r__4, r__5;

	double sqrt(double), cos(double), sin(double);

	static float c1, s1, t1, pr, em1, em2, ct1, pr2, px0, py0, pz0, st1;
	static int ntag;
	extern double ranart_(int *);

	px0 = *px;
	py0 = *py;
	pz0 = *pz;
	*iblock = 144;
	ntag = 0;
	em1 = cc_1.e[*i1 - 1];
	em2 = cc_1.e[*i2 - 1];
	r__2 = *srt;
	r__3 = em1;
	r__4 = em2;
	r__1 = r__2 * r__2 - r__3 * r__3 - r__4 * r__4;
	r__5 = em1 * em2;
	pr2 = r__1 * r__1 - r__5 * r__5 * 4.f;
	if (pr2 <= 0.f)
	{
		pr2 = 1e-9f;
	}
	pr = sqrt(pr2) / (*srt * 2.f);
	c1 = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
	t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
	r__1 = c1;
	s1 = sqrt(1.f - r__1 * r__1);
	ct1 = cos(t1);
	st1 = sin(t1);
	*pz = pr * c1;
	*px = pr * s1 * ct1;
	*py = pr * s1 * st1;
	return 0;
}

int lambar_(int *i1, int *i2, float *srt, float *siglab)
{
	int i__1, i__2;
	float r__1, r__2, r__3;
	double d__1;

	double sqrt(double), pow_dd(double *, double *);

	static float emb, eml, plab, pthr, plab2;

	*siglab = 1e-6f;
	if ((i__1 = ee_1.lb[*i1 - 1], abs(i__1)) >= 14 && (i__2 = ee_1.lb[*i1 - 1], abs(i__2)) <= 17)
	{
		eml = cc_1.e[*i1 - 1];
		emb = cc_1.e[*i2 - 1];
	}
	else
	{
		eml = cc_1.e[*i2 - 1];
		emb = cc_1.e[*i1 - 1];
	}
	r__1 = *srt;
	r__2 = eml;
	r__3 = emb;
	pthr = r__1 * r__1 - r__2 * r__2 - r__3 * r__3;
	if (pthr > 0.f)
	{
		r__1 = pthr / 2.f / emb;
		r__2 = eml;
		plab2 = r__1 * r__1 - r__2 * r__2;
		if (plab2 > 0.f)
		{
			plab = sqrt(plab2);
			d__1 = (double)plab;
			*siglab = .43f / pow_dd(&d__1, &c_b736) + 12.f;
			if (*siglab > 200.f)
			{
				*siglab = 200.f;
			}
		}
	}
	return 0;
}

int distc0_(float *drmax, float *deltr0, float *dt, int *ifirst, float *px1cm, float *py1cm, float *pz1cm, float *x1, float *y1,
			float *z1, float *px1, float *py1, float *pz1, float *em1, float *x2, float *y2, float *z2, float *px2, float *py2, float *pz2, float *em2)
{
	float r__1, r__2, r__3, r__4;

	double sqrt(double);

	static float e1, e2, bbb, ddd, dzz, drcm, dxcm, dycm, dzcm, prcm, p1beta,
		drbeta, relvel, transf;

	*ifirst = -1;
	r__1 = *em1;
	r__2 = *px1;
	r__3 = *py1;
	r__4 = *pz1;
	e1 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	r__1 = *em2;
	r__2 = *px2;
	r__3 = *py2;
	r__4 = *pz2;
	e2 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	p1beta = *px1 * bg_1.betax + *py1 * bg_1.betay + *pz1 * bg_1.betaz;
	transf = bg_1.gamma * (bg_1.gamma * p1beta / (bg_1.gamma + 1) - e1);
	r__1 = *px1cm;
	r__2 = *py1cm;
	r__3 = *pz1cm;
	prcm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
	if (prcm <= 1e-5f)
	{
		return 0;
	}
	drbeta = bg_1.betax * (*x1 - *x2) + bg_1.betay * (*y1 - *y2) + bg_1.betaz * (*z1 - *z2);
	transf = bg_1.gamma * bg_1.gamma * drbeta / (bg_1.gamma + 1);
	dxcm = bg_1.betax * transf + *x1 - *x2;
	dycm = bg_1.betay * transf + *y1 - *y2;
	dzcm = bg_1.betaz * transf + *z1 - *z2;
	r__1 = dxcm;
	r__2 = dycm;
	r__3 = dzcm;
	drcm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
	dzz = (*px1cm * dxcm + *py1cm * dycm + *pz1cm * dzcm) / prcm;
	r__1 = drcm;
	r__2 = dzz;
	if (r__1 * r__1 - r__2 * r__2 <= 0.f)
	{
		bbb = 0.f;
	}
	else
	{
		r__1 = drcm;
		r__2 = dzz;
		bbb = sqrt(r__1 * r__1 - r__2 * r__2);
	}
	if (bbb > *drmax)
	{
		return 0;
	}
	relvel = prcm * (1.f / e1 + 1.f / e2);
	ddd = relvel * *dt * .5f;
	if (dabs(ddd) < dabs(dzz))
	{
		return 0;
	}
	*ifirst = 1;
	return 0;
}

int sbbdm_(float *srt, float *sdprod, int *ianti, int *lbm, float *xmm, float *pfinal)
{
	float r__1, r__2;

	double sqrt(double);

	static float pifactor, pinitial, s, sbbdomega, x1, threshold, fs;
	static int ilb1, ilb2;
	static float snew, sbbdpi;
	extern double fnndpi_(float *), ranart_(int *);
	static float sbbdeta, sbbdrho;

	*sdprod = 0.f;
	sbbdpi = 0.f;
	sbbdrho = 0.f;
	sbbdomega = 0.f;
	sbbdeta = 0.f;
	return 0;

	if (*srt <= leadng_1.em1 + dpi_1.em2)
	{
		return 0;
	}

	ilb1 = abs(leadng_1.lb1);
	ilb2 = abs(dpi_1.lb2);
	r__1 = *srt;
	s = r__1 * r__1;
	r__1 = leadng_1.em1 + dpi_1.em2;
	r__2 = leadng_1.em1 - dpi_1.em2;
	pinitial = sqrt((s - r__1 * r__1) * (s - r__2 * r__2)) / 2.f / *srt;
	fs = fnndpi_(&s);
	if (para8_1.idxsec == 1 || para8_1.idxsec == 2)
	{
	}
	else
	{
		if (ilb1 >= 1 && ilb1 <= 2 && ilb2 >= 1 && ilb2 <= 2)
		{
			pifactor = 1.125f;
		}
		else if (ilb1 >= 1 && ilb1 <= 2 && ilb2 >= 6 && ilb2 <= 9 || ilb2 >=
																			 1 &&
																		 ilb2 <= 2 && ilb1 >= 6 && ilb1 <= 9)
		{
			pifactor = .140625f;
		}
		else if (ilb1 >= 1 && ilb1 <= 2 && ilb2 >= 10 && ilb2 <= 13 || ilb2 >= 1 && ilb2 <= 2 && ilb1 >= 10 && ilb1 <= 13)
		{
			pifactor = .5625f;
		}
		else if (ilb1 >= 6 && ilb1 <= 9 && ilb2 >= 6 && ilb2 <= 9)
		{
			pifactor = .0703125f;
		}
		else if (ilb1 >= 6 && ilb1 <= 9 && ilb2 >= 10 && ilb2 <= 13 || ilb2 >= 6 && ilb2 <= 9 && ilb1 >= 10 && ilb1 <= 13)
		{
			pifactor = .140625f;
		}
		else if (ilb1 >= 10 && ilb1 <= 11 && ilb2 >= 10 && ilb2 <= 11 ||
				 ilb2 >= 12 && ilb2 <= 13 && ilb1 >= 12 && ilb1 <= 13)
		{
			pifactor = 1.125f;
		}
		else if (ilb1 >= 10 && ilb1 <= 11 && ilb2 >= 12 && ilb2 <= 13 ||
				 ilb2 >= 10 && ilb2 <= 11 && ilb1 >= 12 && ilb1 <= 13)
		{
			pifactor = .5625f;
		}
	}
	if (ilb1 * ilb2 == 1)
	{
		*lbm = 5;
		if (*ianti == 1)
		{
			*lbm = 3;
		}
		*xmm = .13957f;
	}
	else if (ilb1 == 2 && ilb2 == 2)
	{
		*lbm = 3;
		if (*ianti == 1)
		{
			*lbm = 5;
		}
		*xmm = .13957f;
	}
	else if (ilb1 * ilb2 == 2)
	{
		*lbm = 4;
		*xmm = .13496f;
	}
	else
	{
		*lbm = (int)(ranart_(&rndf77_1.nseed) * 3) + 3;
		if (*lbm == 4)
		{
			*xmm = .13496f;
		}
		else
		{
			*xmm = .13957f;
		}
	}

	if (*srt >= *xmm + 1.8756f)
	{
		r__1 = *xmm + 1.8756f;
		r__2 = 1.8756f - *xmm;
		*pfinal = sqrt((s - r__1 * r__1) * (s - r__2 * r__2)) / 2.f / *srt;
		if (ilb1 == 1 && ilb2 == 1 || ilb1 == 2 && ilb2 == 2)
		{
			sbbdpi = fs * *pfinal / pinitial / 4.f;
		}
		else if (ilb1 == 1 && ilb2 == 2 || ilb1 == 2 && ilb2 == 1)
		{
			sbbdpi = fs * *pfinal / pinitial / 4.f / 2.f;
		}
		else
		{
			if (para8_1.idxsec == 1)
			{
				sbbdpi = fs * *pfinal / pinitial * 3.f / 16.f;
			}
			else if (para8_1.idxsec == 2 || para8_1.idxsec == 4)
			{
				r__1 = *xmm + 1.8756f, r__2 = leadng_1.em1 + dpi_1.em2;
				threshold = dmax(r__1, r__2);
				r__1 = *srt - threshold + 2.012f;
				snew = r__1 * r__1;
				if (para8_1.idxsec == 2)
				{
					sbbdpi = fnndpi_(&snew) * *pfinal / pinitial * 3.f / 16.f;
				}
				else if (para8_1.idxsec == 4)
				{
					sbbdpi = fnndpi_(&snew) * *pfinal / pinitial / 6.f *
							 pifactor;
				}
			}
			else if (para8_1.idxsec == 3)
			{
				sbbdpi = fs * *pfinal / pinitial / 6.f * pifactor;
			}
		}
	}

	if (*srt > 2.6456f)
	{
		*pfinal = sqrt((s - 6.9991993599999995f) * (s - 1.2223513599999998f)) / 2.f / *srt;
		if (para8_1.idxsec == 1)
		{
			sbbdrho = fs * *pfinal / pinitial * 3.f / 16.f;
		}
		else if (para8_1.idxsec == 2 || para8_1.idxsec == 4)
		{
			r__1 = 2.6456f, r__2 = leadng_1.em1 + dpi_1.em2;
			threshold = dmax(r__1, r__2);
			r__1 = *srt - threshold + 2.012f;
			snew = r__1 * r__1;
			if (para8_1.idxsec == 2)
			{
				sbbdrho = fnndpi_(&snew) * *pfinal / pinitial * 3.f / 16.f;
			}
			else if (para8_1.idxsec == 4)
			{
				sbbdrho = fnndpi_(&snew) * *pfinal / pinitial / 6.f * (pifactor * 3.f);
			}
		}
		else if (para8_1.idxsec == 3)
		{
			sbbdrho = fs * *pfinal / pinitial / 6.f * (pifactor * 3.f);
		}
	}

	if (*srt > 2.6576f)
	{
		*pfinal = sqrt((s - 7.0628377599999999f) * (s - 1.1959609599999999f)) / 2.f / *srt;
		if (para8_1.idxsec == 1)
		{
			sbbdomega = fs * *pfinal / pinitial * 3.f / 16.f;
		}
		else if (para8_1.idxsec == 2 || para8_1.idxsec == 4)
		{
			r__1 = 2.6576f, r__2 = leadng_1.em1 + dpi_1.em2;
			threshold = dmax(r__1, r__2);
			r__1 = *srt - threshold + 2.012f;
			snew = r__1 * r__1;
			if (para8_1.idxsec == 2)
			{
				sbbdomega = fnndpi_(&snew) * *pfinal / pinitial * 3.f / 16.f;
			}
			else if (para8_1.idxsec == 4)
			{
				sbbdomega = fnndpi_(&snew) * *pfinal / pinitial / 6.f *
							pifactor;
			}
		}
		else if (para8_1.idxsec == 3)
		{
			sbbdomega = fs * *pfinal / pinitial / 6.f * pifactor;
		}
	}

	if (*srt > 2.4236f)
	{
		*pfinal = sqrt((s - 5.8738369600000002f) * (s - 1.7625217599999996f)) / 2.f / *srt;
		if (para8_1.idxsec == 1)
		{
			sbbdeta = fs * *pfinal / pinitial * 3.f / 16.f;
		}
		else if (para8_1.idxsec == 2 || para8_1.idxsec == 4)
		{
			r__1 = 2.4236f, r__2 = leadng_1.em1 + dpi_1.em2;
			threshold = dmax(r__1, r__2);
			r__1 = *srt - threshold + 2.012f;
			snew = r__1 * r__1;
			if (para8_1.idxsec == 2)
			{
				sbbdeta = fnndpi_(&snew) * *pfinal / pinitial * 3.f / 16.f;
			}
			else if (para8_1.idxsec == 4)
			{
				sbbdeta = fnndpi_(&snew) * *pfinal / pinitial / 6.f * (pifactor / 3.f);
			}
		}
		else if (para8_1.idxsec == 3)
		{
			sbbdeta = fs * *pfinal / pinitial / 6.f * (pifactor / 3.f);
		}
	}

	*sdprod = sbbdpi + sbbdrho + sbbdomega + sbbdeta;
	if (*sdprod <= 0.f)
	{
		return 0;
	}

	x1 = ranart_(&rndf77_1.nseed);
	if (x1 <= sbbdpi / *sdprod)
	{
	}
	else if (x1 <= (sbbdpi + sbbdrho) / *sdprod)
	{
		*lbm = (int)(ranart_(&rndf77_1.nseed) * 3) + 25;
		*xmm = .77f;
	}
	else if (x1 <= (sbbdpi + sbbdrho + sbbdomega) / *sdprod)
	{
		*lbm = 28;
		*xmm = .782f;
	}
	else
	{
		*lbm = 0;
		*xmm = .548f;
	}

	return 0;
}

int bbdangle_(float *pxd, float *pyd, float *pzd, int *nt,
			  int *ipert1, int *ianti, int *idloop, float *pfinal, float *dprob1, int *lbm)
{
	float r__1;

	double sqrt(double), cos(double), sin(double);
	int s_wsle(cilist *), do_lio(int *, int *, char *, ftnlen),
		e_wsle(void);

	static float c1, s1, t1, ct1, st1, dprob;
	extern double ranart_(int *);

	static cilist io___2331 = {0, 91, 0, 0, 0};
	static cilist io___2332 = {0, 91, 0, 0, 0};
	static cilist io___2333 = {0, 91, 0, 0, 0};
	static cilist io___2334 = {0, 91, 0, 0, 0};

	c1 = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
	t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
	r__1 = c1;
	s1 = sqrt(1.f - r__1 * r__1);
	ct1 = cos(t1);
	st1 = sin(t1);
	*pzd = *pfinal * c1;
	*pxd = *pfinal * s1 * ct1;
	*pyd = *pfinal * s1 * st1;
	if (para8_1.idpert == 1 && para8_1.npertd >= 1)
	{
		dprob = *dprob1;
	}
	else if (para8_1.idpert == 2 && para8_1.npertd >= 1)
	{
		dprob = 1.f / (float)para8_1.npertd;
	}
	if (*ianti == 0)
	{
		if (para8_1.idpert == 0 || para8_1.idpert == 1 && *ipert1 == 0 ||
			para8_1.idpert == 2 && *idloop == para8_1.npertd + 1)
		{
			s_wsle(&io___2331);
			do_lio(&c__3, &c__1, (char *)&leadng_1.lb1, (ftnlen)sizeof(int));
			do_lio(&c__9, &c__1, " *", (ftnlen)2);
			do_lio(&c__3, &c__1, (char *)&dpi_1.lb2, (ftnlen)sizeof(int));
			do_lio(&c__9, &c__1, " ->d+", (ftnlen)5);
			do_lio(&c__3, &c__1, (char *)&(*lbm), (ftnlen)sizeof(int));
			do_lio(&c__9, &c__1, " (regular d prodn)    @evt#", (ftnlen)27);
			do_lio(&c__3, &c__1, (char *)&arevt_1.iaevt, (ftnlen)sizeof(int));
			do_lio(&c__9, &c__1, " @nt=", (ftnlen)5);
			do_lio(&c__3, &c__1, (char *)&(*nt), (ftnlen)sizeof(int));
			e_wsle();
		}
		else if ((para8_1.idpert == 1 || para8_1.idpert == 2) && *idloop ==
																	 para8_1.npertd)
		{
			s_wsle(&io___2332);
			do_lio(&c__3, &c__1, (char *)&leadng_1.lb1, (ftnlen)sizeof(int));
			do_lio(&c__9, &c__1, " *", (ftnlen)2);
			do_lio(&c__3, &c__1, (char *)&dpi_1.lb2, (ftnlen)sizeof(int));
			do_lio(&c__9, &c__1, " ->d+", (ftnlen)5);
			do_lio(&c__3, &c__1, (char *)&(*lbm), (ftnlen)sizeof(int));
			do_lio(&c__9, &c__1, " (pert d prodn)       @evt#", (ftnlen)27);
			do_lio(&c__3, &c__1, (char *)&arevt_1.iaevt, (ftnlen)sizeof(int));
			do_lio(&c__9, &c__1, " @nt=", (ftnlen)5);
			do_lio(&c__3, &c__1, (char *)&(*nt), (ftnlen)sizeof(int));
			do_lio(&c__9, &c__1, " @prob=", (ftnlen)7);
			do_lio(&c__4, &c__1, (char *)&dprob, (ftnlen)sizeof(float));
			e_wsle();
		}
	}
	else
	{
		if (para8_1.idpert == 0 || para8_1.idpert == 1 && *ipert1 == 0 ||
			para8_1.idpert == 2 && *idloop == para8_1.npertd + 1)
		{
			s_wsle(&io___2333);
			do_lio(&c__3, &c__1, (char *)&leadng_1.lb1, (ftnlen)sizeof(int));
			do_lio(&c__9, &c__1, " *", (ftnlen)2);
			do_lio(&c__3, &c__1, (char *)&dpi_1.lb2, (ftnlen)sizeof(int));
			do_lio(&c__9, &c__1, " ->d+", (ftnlen)5);
			do_lio(&c__3, &c__1, (char *)&(*lbm), (ftnlen)sizeof(int));
			do_lio(&c__9, &c__1, " (regular dbar prodn) @evt#", (ftnlen)27);
			do_lio(&c__3, &c__1, (char *)&arevt_1.iaevt, (ftnlen)sizeof(int));
			do_lio(&c__9, &c__1, " @nt=", (ftnlen)5);
			do_lio(&c__3, &c__1, (char *)&(*nt), (ftnlen)sizeof(int));
			e_wsle();
		}
		else if ((para8_1.idpert == 1 || para8_1.idpert == 2) && *idloop ==
																	 para8_1.npertd)
		{
			s_wsle(&io___2334);
			do_lio(&c__3, &c__1, (char *)&leadng_1.lb1, (ftnlen)sizeof(int));
			do_lio(&c__9, &c__1, " *", (ftnlen)2);
			do_lio(&c__3, &c__1, (char *)&dpi_1.lb2, (ftnlen)sizeof(int));
			do_lio(&c__9, &c__1, " ->d+", (ftnlen)5);
			do_lio(&c__3, &c__1, (char *)&(*lbm), (ftnlen)sizeof(int));
			do_lio(&c__9, &c__1, " (pert dbar prodn)    @evt#", (ftnlen)27);
			do_lio(&c__3, &c__1, (char *)&arevt_1.iaevt, (ftnlen)sizeof(int));
			do_lio(&c__9, &c__1, " @nt=", (ftnlen)5);
			do_lio(&c__3, &c__1, (char *)&(*nt), (ftnlen)sizeof(int));
			do_lio(&c__9, &c__1, " @prob=", (ftnlen)7);
			do_lio(&c__4, &c__1, (char *)&dprob, (ftnlen)sizeof(float));
			e_wsle();
		}
	}

	return 0;
}

int sdmbb_(float *srt, float *sdm, int *ianti)
{
	float r__1, r__2;

	double sqrt(double);

	static float pinitial, s, threshold, xnnfactor, fs, snew;
	extern double fdpiel_(float *);
	static float pfinal;
	extern double fnndpi_(float *), ranart_(int *);

	*sdm = 0.f;
	dpisig_1.sdmel = 0.f;
	dpisig_1.sdmnn = 0.f;
	dpisig_1.sdmnd = 0.f;
	dpisig_1.sdmns = 0.f;
	dpisig_1.sdmnp = 0.f;
	dpisig_1.sdmdd = 0.f;
	dpisig_1.sdmds = 0.f;
	dpisig_1.sdmdp = 0.f;
	dpisig_1.sdmss = 0.f;
	dpisig_1.sdmsp = 0.f;
	dpisig_1.sdmpp = 0.f;
	return 0;

	if (*srt <= leadng_1.em1 + dpi_1.em2)
	{
		return 0;
	}
	r__1 = *srt;
	s = r__1 * r__1;
	r__1 = leadng_1.em1 + dpi_1.em2;
	r__2 = leadng_1.em1 - dpi_1.em2;
	pinitial = sqrt((s - r__1 * r__1) * (s - r__2 * r__2)) / 2.f / *srt;
	fs = fnndpi_(&s);
	if (para8_1.idxsec == 1 || para8_1.idxsec == 2)
	{
		if (leadng_1.lb1 >= 3 && leadng_1.lb1 <= 5 || dpi_1.lb2 >= 3 &&
														  dpi_1.lb2 <= 5)
		{
			xnnfactor = .88888888888888884f;
		}
		else if (leadng_1.lb1 >= 25 && leadng_1.lb1 <= 27 || dpi_1.lb2 >=
																	 25 &&
																 dpi_1.lb2 <= 27)
		{
			xnnfactor = .29629629629629628f;
		}
		else if (leadng_1.lb1 == 28 || dpi_1.lb2 == 28)
		{
			xnnfactor = .88888888888888884f;
		}
		else if (leadng_1.lb1 == 0 || dpi_1.lb2 == 0)
		{
			xnnfactor = 2.6666666666666665f;
		}
	}
	else
	{
	}
	if (para8_1.idxsec == 1 || para8_1.idxsec == 3)
	{
		dpisig_1.sdmel = fdpiel_(&s);
	}
	else if (para8_1.idxsec == 2 || para8_1.idxsec == 4)
	{
		threshold = leadng_1.em1 + dpi_1.em2;
		r__1 = *srt - threshold + 2.012f;
		snew = r__1 * r__1;
		dpisig_1.sdmel = fdpiel_(&snew);
	}

	if ((leadng_1.lb1 == 5 || dpi_1.lb2 == 5 || leadng_1.lb1 == 27 ||
		 dpi_1.lb2 == 27) &&
			*ianti == 0 ||
		(leadng_1.lb1 == 3 ||
		 dpi_1.lb2 == 3 || leadng_1.lb1 == 25 || dpi_1.lb2 == 25) &&
			*ianti == 1)
	{
		dpifsl_1.lbnn1 = 1;
		dpifsl_1.lbnn2 = 1;
		dpifsm_1.xmnn1 = .93828f;
		dpifsm_1.xmnn2 = .93828f;
	}
	else if (leadng_1.lb1 == 3 || dpi_1.lb2 == 3 || leadng_1.lb1 == 26 ||
			 dpi_1.lb2 == 26 || leadng_1.lb1 == 28 || dpi_1.lb2 == 28 ||
			 leadng_1.lb1 == 0 || dpi_1.lb2 == 0)
	{
		dpifsl_1.lbnn1 = 2;
		dpifsl_1.lbnn2 = 1;
		dpifsm_1.xmnn1 = .939457f;
		dpifsm_1.xmnn2 = .93828f;
	}
	else
	{
		dpifsl_1.lbnn1 = 2;
		dpifsl_1.lbnn2 = 2;
		dpifsm_1.xmnn1 = .939457f;
		dpifsm_1.xmnn2 = .939457f;
	}
	if (*srt > dpifsm_1.xmnn1 + dpifsm_1.xmnn2)
	{
		r__1 = dpifsm_1.xmnn1 + dpifsm_1.xmnn2;
		r__2 = dpifsm_1.xmnn1 - dpifsm_1.xmnn2;
		pfinal = sqrt((s - r__1 * r__1) * (s - r__2 * r__2)) / 2.f / *srt;
		if (para8_1.idxsec == 1)
		{
			dpisig_1.sdmnn = fs * pfinal / pinitial * 3.f / 16.f * xnnfactor;
		}
		else if (para8_1.idxsec == 2 || para8_1.idxsec == 4)
		{
			r__1 = dpifsm_1.xmnn1 + dpifsm_1.xmnn2, r__2 = leadng_1.em1 +
														   dpi_1.em2;
			threshold = dmax(r__1, r__2);
			r__1 = *srt - threshold + 2.012f;
			snew = r__1 * r__1;
			if (para8_1.idxsec == 2)
			{
				dpisig_1.sdmnn = fnndpi_(&snew) * pfinal / pinitial * 3.f /
								 16.f * xnnfactor;
			}
			else if (para8_1.idxsec == 4)
			{
				dpisig_1.sdmnn = fnndpi_(&snew) * pfinal / pinitial / 6.f;
			}
		}
		else if (para8_1.idxsec == 3)
		{
			dpisig_1.sdmnn = fs * pfinal / pinitial / 6.f;
		}
	}

	dpifsl_1.lbnd1 = (int)(ranart_(&rndf77_1.nseed) * 2) + 1;
	dpifsl_1.lbnd2 = (int)(ranart_(&rndf77_1.nseed) * 4) + 6;
	if (dpifsl_1.lbnd1 == 1)
	{
		dpifsm_1.xmnd1 = .93828f;
	}
	else if (dpifsl_1.lbnd1 == 2)
	{
		dpifsm_1.xmnd1 = .939457f;
	}
	dpifsm_1.xmnd2 = 1.232f;
	if (*srt > dpifsm_1.xmnd1 + dpifsm_1.xmnd2)
	{
		r__1 = dpifsm_1.xmnd1 + dpifsm_1.xmnd2;
		r__2 = dpifsm_1.xmnd1 - dpifsm_1.xmnd2;
		pfinal = sqrt((s - r__1 * r__1) * (s - r__2 * r__2)) / 2.f / *srt;
		if (para8_1.idxsec == 1)
		{
			dpisig_1.sdmnd = fs * pfinal / pinitial * 3.f / 16.f * (xnnfactor * 8.f);
		}
		else if (para8_1.idxsec == 2 || para8_1.idxsec == 4)
		{
			r__1 = dpifsm_1.xmnd1 + dpifsm_1.xmnd2, r__2 = leadng_1.em1 +
														   dpi_1.em2;
			threshold = dmax(r__1, r__2);
			r__1 = *srt - threshold + 2.012f;
			snew = r__1 * r__1;
			if (para8_1.idxsec == 2)
			{
				dpisig_1.sdmnd = fnndpi_(&snew) * pfinal / pinitial * 3.f /
								 16.f * (xnnfactor * 8.f);
			}
			else if (para8_1.idxsec == 4)
			{
				dpisig_1.sdmnd = fnndpi_(&snew) * pfinal / pinitial / 6.f;
			}
		}
		else if (para8_1.idxsec == 3)
		{
			dpisig_1.sdmnd = fs * pfinal / pinitial / 6.f;
		}
	}

	dpifsl_1.lbns1 = (int)(ranart_(&rndf77_1.nseed) * 2) + 1;
	dpifsl_1.lbns2 = (int)(ranart_(&rndf77_1.nseed) * 2) + 10;
	if (dpifsl_1.lbns1 == 1)
	{
		dpifsm_1.xmns1 = .93828f;
	}
	else if (dpifsl_1.lbns1 == 2)
	{
		dpifsm_1.xmns1 = .939457f;
	}
	dpifsm_1.xmns2 = 1.44f;
	if (*srt > dpifsm_1.xmns1 + dpifsm_1.xmns2)
	{
		r__1 = dpifsm_1.xmns1 + dpifsm_1.xmns2;
		r__2 = dpifsm_1.xmns1 - dpifsm_1.xmns2;
		pfinal = sqrt((s - r__1 * r__1) * (s - r__2 * r__2)) / 2.f / *srt;
		if (para8_1.idxsec == 1)
		{
			dpisig_1.sdmns = fs * pfinal / pinitial * 3.f / 16.f * (xnnfactor * 2.f);
		}
		else if (para8_1.idxsec == 2 || para8_1.idxsec == 4)
		{
			r__1 = dpifsm_1.xmns1 + dpifsm_1.xmns2, r__2 = leadng_1.em1 +
														   dpi_1.em2;
			threshold = dmax(r__1, r__2);
			r__1 = *srt - threshold + 2.012f;
			snew = r__1 * r__1;
			if (para8_1.idxsec == 2)
			{
				dpisig_1.sdmns = fnndpi_(&snew) * pfinal / pinitial * 3.f /
								 16.f * (xnnfactor * 2.f);
			}
			else if (para8_1.idxsec == 4)
			{
				dpisig_1.sdmns = fnndpi_(&snew) * pfinal / pinitial / 6.f;
			}
		}
		else if (para8_1.idxsec == 3)
		{
			dpisig_1.sdmns = fs * pfinal / pinitial / 6.f;
		}
	}

	dpifsl_1.lbnp1 = (int)(ranart_(&rndf77_1.nseed) * 2) + 1;
	dpifsl_1.lbnp2 = (int)(ranart_(&rndf77_1.nseed) * 2) + 12;
	if (dpifsl_1.lbnp1 == 1)
	{
		dpifsm_1.xmnp1 = .93828f;
	}
	else if (dpifsl_1.lbnp1 == 2)
	{
		dpifsm_1.xmnp1 = .939457f;
	}
	dpifsm_1.xmnp2 = 1.535f;
	if (*srt > dpifsm_1.xmnp1 + dpifsm_1.xmnp2)
	{
		r__1 = dpifsm_1.xmnp1 + dpifsm_1.xmnp2;
		r__2 = dpifsm_1.xmnp1 - dpifsm_1.xmnp2;
		pfinal = sqrt((s - r__1 * r__1) * (s - r__2 * r__2)) / 2.f / *srt;
		if (para8_1.idxsec == 1)
		{
			dpisig_1.sdmnp = fs * pfinal / pinitial * 3.f / 16.f * (xnnfactor * 2.f);
		}
		else if (para8_1.idxsec == 2 || para8_1.idxsec == 4)
		{
			r__1 = dpifsm_1.xmnp1 + dpifsm_1.xmnp2, r__2 = leadng_1.em1 +
														   dpi_1.em2;
			threshold = dmax(r__1, r__2);
			r__1 = *srt - threshold + 2.012f;
			snew = r__1 * r__1;
			if (para8_1.idxsec == 2)
			{
				dpisig_1.sdmnp = fnndpi_(&snew) * pfinal / pinitial * 3.f /
								 16.f * (xnnfactor * 2.f);
			}
			else if (para8_1.idxsec == 4)
			{
				dpisig_1.sdmnp = fnndpi_(&snew) * pfinal / pinitial / 6.f;
			}
		}
		else if (para8_1.idxsec == 3)
		{
			dpisig_1.sdmnp = fs * pfinal / pinitial / 6.f;
		}
	}

	dpifsl_1.lbdd1 = (int)(ranart_(&rndf77_1.nseed) * 4) + 6;
	dpifsl_1.lbdd2 = (int)(ranart_(&rndf77_1.nseed) * 4) + 6;
	dpifsm_1.xmdd1 = 1.232f;
	dpifsm_1.xmdd2 = 1.232f;
	if (*srt > dpifsm_1.xmdd1 + dpifsm_1.xmdd2)
	{
		r__1 = dpifsm_1.xmdd1 + dpifsm_1.xmdd2;
		r__2 = dpifsm_1.xmdd1 - dpifsm_1.xmdd2;
		pfinal = sqrt((s - r__1 * r__1) * (s - r__2 * r__2)) / 2.f / *srt;
		if (para8_1.idxsec == 1)
		{
			dpisig_1.sdmdd = fs * pfinal / pinitial * 3.f / 16.f * (xnnfactor * 16.f);
		}
		else if (para8_1.idxsec == 2 || para8_1.idxsec == 4)
		{
			r__1 = dpifsm_1.xmdd1 + dpifsm_1.xmdd2, r__2 = leadng_1.em1 +
														   dpi_1.em2;
			threshold = dmax(r__1, r__2);
			r__1 = *srt - threshold + 2.012f;
			snew = r__1 * r__1;
			if (para8_1.idxsec == 2)
			{
				dpisig_1.sdmdd = fnndpi_(&snew) * pfinal / pinitial * 3.f /
								 16.f * (xnnfactor * 16.f);
			}
			else if (para8_1.idxsec == 4)
			{
				dpisig_1.sdmdd = fnndpi_(&snew) * pfinal / pinitial / 6.f;
			}
		}
		else if (para8_1.idxsec == 3)
		{
			dpisig_1.sdmdd = fs * pfinal / pinitial / 6.f;
		}
	}

	dpifsl_1.lbds1 = (int)(ranart_(&rndf77_1.nseed) * 4) + 6;
	dpifsl_1.lbds2 = (int)(ranart_(&rndf77_1.nseed) * 2) + 10;
	dpifsm_1.xmds1 = 1.232f;
	dpifsm_1.xmds2 = 1.44f;
	if (*srt > dpifsm_1.xmds1 + dpifsm_1.xmds2)
	{
		r__1 = dpifsm_1.xmds1 + dpifsm_1.xmds2;
		r__2 = dpifsm_1.xmds1 - dpifsm_1.xmds2;
		pfinal = sqrt((s - r__1 * r__1) * (s - r__2 * r__2)) / 2.f / *srt;
		if (para8_1.idxsec == 1)
		{
			dpisig_1.sdmds = fs * pfinal / pinitial * 3.f / 16.f * (xnnfactor * 8.f);
		}
		else if (para8_1.idxsec == 2 || para8_1.idxsec == 4)
		{
			r__1 = dpifsm_1.xmds1 + dpifsm_1.xmds2, r__2 = leadng_1.em1 +
														   dpi_1.em2;
			threshold = dmax(r__1, r__2);
			r__1 = *srt - threshold + 2.012f;
			snew = r__1 * r__1;
			if (para8_1.idxsec == 2)
			{
				dpisig_1.sdmds = fnndpi_(&snew) * pfinal / pinitial * 3.f /
								 16.f * (xnnfactor * 8.f);
			}
			else if (para8_1.idxsec == 4)
			{
				dpisig_1.sdmds = fnndpi_(&snew) * pfinal / pinitial / 6.f;
			}
		}
		else if (para8_1.idxsec == 3)
		{
			dpisig_1.sdmds = fs * pfinal / pinitial / 6.f;
		}
	}

	dpifsl_1.lbdp1 = (int)(ranart_(&rndf77_1.nseed) * 4) + 6;
	dpifsl_1.lbdp2 = (int)(ranart_(&rndf77_1.nseed) * 2) + 12;
	dpifsm_1.xmdp1 = 1.232f;
	dpifsm_1.xmdp2 = 1.535f;
	if (*srt > dpifsm_1.xmdp1 + dpifsm_1.xmdp2)
	{
		r__1 = dpifsm_1.xmdp1 + dpifsm_1.xmdp2;
		r__2 = dpifsm_1.xmdp1 - dpifsm_1.xmdp2;
		pfinal = sqrt((s - r__1 * r__1) * (s - r__2 * r__2)) / 2.f / *srt;
		if (para8_1.idxsec == 1)
		{
			dpisig_1.sdmdp = fs * pfinal / pinitial * 3.f / 16.f * (xnnfactor * 8.f);
		}
		else if (para8_1.idxsec == 2 || para8_1.idxsec == 4)
		{
			r__1 = dpifsm_1.xmdp1 + dpifsm_1.xmdp2, r__2 = leadng_1.em1 +
														   dpi_1.em2;
			threshold = dmax(r__1, r__2);
			r__1 = *srt - threshold + 2.012f;
			snew = r__1 * r__1;
			if (para8_1.idxsec == 2)
			{
				dpisig_1.sdmdp = fnndpi_(&snew) * pfinal / pinitial * 3.f /
								 16.f * (xnnfactor * 8.f);
			}
			else if (para8_1.idxsec == 4)
			{
				dpisig_1.sdmdp = fnndpi_(&snew) * pfinal / pinitial / 6.f;
			}
		}
		else if (para8_1.idxsec == 3)
		{
			dpisig_1.sdmdp = fs * pfinal / pinitial / 6.f;
		}
	}

	dpifsl_1.lbss1 = (int)(ranart_(&rndf77_1.nseed) * 2) + 10;
	dpifsl_1.lbss2 = (int)(ranart_(&rndf77_1.nseed) * 2) + 10;
	dpifsm_1.xmss1 = 1.44f;
	dpifsm_1.xmss2 = 1.44f;
	if (*srt > dpifsm_1.xmss1 + dpifsm_1.xmss2)
	{
		r__1 = dpifsm_1.xmss1 + dpifsm_1.xmss2;
		r__2 = dpifsm_1.xmss1 - dpifsm_1.xmss2;
		pfinal = sqrt((s - r__1 * r__1) * (s - r__2 * r__2)) / 2.f / *srt;
		if (para8_1.idxsec == 1)
		{
			dpisig_1.sdmss = fs * pfinal / pinitial * 3.f / 16.f * xnnfactor;
		}
		else if (para8_1.idxsec == 2 || para8_1.idxsec == 4)
		{
			r__1 = dpifsm_1.xmss1 + dpifsm_1.xmss2, r__2 = leadng_1.em1 +
														   dpi_1.em2;
			threshold = dmax(r__1, r__2);
			r__1 = *srt - threshold + 2.012f;
			snew = r__1 * r__1;
			if (para8_1.idxsec == 2)
			{
				dpisig_1.sdmss = fnndpi_(&snew) * pfinal / pinitial * 3.f /
								 16.f * xnnfactor;
			}
			else if (para8_1.idxsec == 4)
			{
				dpisig_1.sdmss = fnndpi_(&snew) * pfinal / pinitial / 6.f;
			}
		}
		else if (para8_1.idxsec == 3)
		{
			dpisig_1.sdmns = fs * pfinal / pinitial / 6.f;
		}
	}

	dpifsl_1.lbsp1 = (int)(ranart_(&rndf77_1.nseed) * 2) + 10;
	dpifsl_1.lbsp2 = (int)(ranart_(&rndf77_1.nseed) * 2) + 12;
	dpifsm_1.xmsp1 = 1.44f;
	dpifsm_1.xmsp2 = 1.535f;
	if (*srt > dpifsm_1.xmsp1 + dpifsm_1.xmsp2)
	{
		r__1 = dpifsm_1.xmsp1 + dpifsm_1.xmsp2;
		r__2 = dpifsm_1.xmsp1 - dpifsm_1.xmsp2;
		pfinal = sqrt((s - r__1 * r__1) * (s - r__2 * r__2)) / 2.f / *srt;
		if (para8_1.idxsec == 1)
		{
			dpisig_1.sdmsp = fs * pfinal / pinitial * 3.f / 16.f * (xnnfactor * 2.f);
		}
		else if (para8_1.idxsec == 2 || para8_1.idxsec == 4)
		{
			r__1 = dpifsm_1.xmsp1 + dpifsm_1.xmsp2, r__2 = leadng_1.em1 +
														   dpi_1.em2;
			threshold = dmax(r__1, r__2);
			r__1 = *srt - threshold + 2.012f;
			snew = r__1 * r__1;
			if (para8_1.idxsec == 2)
			{
				dpisig_1.sdmsp = fnndpi_(&snew) * pfinal / pinitial * 3.f /
								 16.f * (xnnfactor * 2.f);
			}
			else if (para8_1.idxsec == 4)
			{
				dpisig_1.sdmsp = fnndpi_(&snew) * pfinal / pinitial / 6.f;
			}
		}
		else if (para8_1.idxsec == 3)
		{
			dpisig_1.sdmsp = fs * pfinal / pinitial / 6.f;
		}
	}

	dpifsl_1.lbpp1 = (int)(ranart_(&rndf77_1.nseed) * 2) + 12;
	dpifsl_1.lbpp2 = (int)(ranart_(&rndf77_1.nseed) * 2) + 12;
	dpifsm_1.xmpp1 = 1.535f;
	dpifsm_1.xmpp2 = 1.535f;
	if (*srt > dpifsm_1.xmpp1 + dpifsm_1.xmpp2)
	{
		r__1 = dpifsm_1.xmpp1 + dpifsm_1.xmpp2;
		r__2 = dpifsm_1.xmpp1 - dpifsm_1.xmpp2;
		pfinal = sqrt((s - r__1 * r__1) * (s - r__2 * r__2)) / 2.f / *srt;
		if (para8_1.idxsec == 1)
		{
			dpisig_1.sdmpp = fs * pfinal / pinitial * 3.f / 16.f * xnnfactor;
		}
		else if (para8_1.idxsec == 2 || para8_1.idxsec == 4)
		{
			r__1 = dpifsm_1.xmpp1 + dpifsm_1.xmpp2, r__2 = leadng_1.em1 +
														   dpi_1.em2;
			threshold = dmax(r__1, r__2);
			r__1 = *srt - threshold + 2.012f;
			snew = r__1 * r__1;
			if (para8_1.idxsec == 2)
			{
				dpisig_1.sdmpp = fnndpi_(&snew) * pfinal / pinitial * 3.f /
								 16.f * xnnfactor;
			}
			else if (para8_1.idxsec == 4)
			{
				dpisig_1.sdmpp = fnndpi_(&snew) * pfinal / pinitial / 6.f;
			}
		}
		else if (para8_1.idxsec == 3)
		{
			dpisig_1.sdmpp = fs * pfinal / pinitial / 6.f;
		}
	}

	*sdm = dpisig_1.sdmel + dpisig_1.sdmnn + dpisig_1.sdmnd + dpisig_1.sdmns + dpisig_1.sdmnp + dpisig_1.sdmdd + dpisig_1.sdmds +
		   dpisig_1.sdmdp + dpisig_1.sdmss + dpisig_1.sdmsp + dpisig_1.sdmpp;
	if (*ianti == 1)
	{
		dpifsl_1.lbnn1 = -dpifsl_1.lbnn1;
		dpifsl_1.lbnn2 = -dpifsl_1.lbnn2;
		dpifsl_1.lbnd1 = -dpifsl_1.lbnd1;
		dpifsl_1.lbnd2 = -dpifsl_1.lbnd2;
		dpifsl_1.lbns1 = -dpifsl_1.lbns1;
		dpifsl_1.lbns2 = -dpifsl_1.lbns2;
		dpifsl_1.lbnp1 = -dpifsl_1.lbnp1;
		dpifsl_1.lbnp2 = -dpifsl_1.lbnp2;
		dpifsl_1.lbdd1 = -dpifsl_1.lbdd1;
		dpifsl_1.lbdd2 = -dpifsl_1.lbdd2;
		dpifsl_1.lbds1 = -dpifsl_1.lbds1;
		dpifsl_1.lbds2 = -dpifsl_1.lbds2;
		dpifsl_1.lbdp1 = -dpifsl_1.lbdp1;
		dpifsl_1.lbdp2 = -dpifsl_1.lbdp2;
		dpifsl_1.lbss1 = -dpifsl_1.lbss1;
		dpifsl_1.lbss2 = -dpifsl_1.lbss2;
		dpifsl_1.lbsp1 = -dpifsl_1.lbsp1;
		dpifsl_1.lbsp2 = -dpifsl_1.lbsp2;
		dpifsl_1.lbpp1 = -dpifsl_1.lbpp1;
		dpifsl_1.lbpp2 = -dpifsl_1.lbpp2;
	}
	return 0;
}

int crdmbb_(float *px, float *py, float *pz, float *srt, int *i1, int *i2, int *iblock, int *ntag, float *sig, int *nt, int *ianti)
{
	float r__1, r__2, r__3, r__4;

	int s_wsle(cilist *), do_lio(int *, int *, char *, ftnlen),
		e_wsle(void);
	double sqrt(double);
	int s_stop(char *, ftnlen);

	static float s;
	extern int dmelangle_(float *, float *, float *, float *);
	static float x1;
	static int idm, lbm;
	static float pxn, pyn, pzn;
	static int lbb1, lbb2;
	static float e1cm, e2cm, xmb1, pt1d, pt2d, pt3d, xmb2, edcm, pt1i1, pt2i1,
		pt3i1, pt1i2, pt2i2, pt3i2;
	static int ideut;
	static float p1beta, p2beta, pdbeta, pfinal;
	extern double ranart_(int *);
	static float transf;
	extern int rotate_(float *, float *, float *, float *, float *, float *), dmangle_(float *, float *, float *, int *, int *,
																				 float *, int *);

	static cilist io___2347 = {0, 91, 0, 0, 0};
	static cilist io___2348 = {0, 91, 0, 0, 0};
	static cilist io___2359 = {0, 91, 0, 0, 0};
	static cilist io___2360 = {0, 91, 0, 0, 0};
	static cilist io___2365 = {0, 91, 0, 0, 0};
	static cilist io___2366 = {0, 91, 0, 0, 0};
	static cilist io___2367 = {0, 6, 0, 0, 0};

	*iblock = 0;
	*ntag = 0;
	leadng_1.em1 = cc_1.e[*i1 - 1];
	dpi_1.em2 = cc_1.e[*i2 - 1];
	r__1 = *srt;
	s = r__1 * r__1;
	if (*sig <= 0.f)
	{
		return 0;
	}

	if (abs(leadng_1.lb1) == 42)
	{
		ideut = *i1;
		lbm = dpi_1.lb2;
		idm = *i2;
	}
	else
	{
		ideut = *i2;
		lbm = leadng_1.lb1;
		idm = *i1;
	}
	if ((para8_1.idpert == 1 || para8_1.idpert == 2) && dpert_1.dpertp[ideut - 1] != 1.f)
	{
		x1 = ranart_(&rndf77_1.nseed);
		if (x1 <= dpisig_1.sdmel / *sig)
		{
			if (*ianti == 0)
			{
				s_wsle(&io___2347);
				do_lio(&c__9, &c__1, "  d+", (ftnlen)4);
				do_lio(&c__3, &c__1, (char *)&lbm, (ftnlen)sizeof(int));
				do_lio(&c__9, &c__1, " (pert d M elastic) @nt=", (ftnlen)24);
				do_lio(&c__3, &c__1, (char *)&(*nt), (ftnlen)sizeof(int));
				do_lio(&c__9, &c__1, " @prob=", (ftnlen)7);
				do_lio(&c__4, &c__1, (char *)&dpert_1.dpertp[ideut - 1], (ftnlen)sizeof(float));
				e_wsle();
			}
			else
			{
				s_wsle(&io___2348);
				do_lio(&c__9, &c__1, "  d+", (ftnlen)4);
				do_lio(&c__3, &c__1, (char *)&lbm, (ftnlen)sizeof(int));
				do_lio(&c__9, &c__1, " (pert dbar M elastic) @nt=", (ftnlen)27);
				do_lio(&c__3, &c__1, (char *)&(*nt), (ftnlen)sizeof(int));
				do_lio(&c__9, &c__1, " @prob=", (ftnlen)7);
				do_lio(&c__4, &c__1, (char *)&dpert_1.dpertp[ideut - 1], (ftnlen)sizeof(float));
				e_wsle();
			}
			r__1 = leadng_1.em1 + dpi_1.em2;
			r__2 = leadng_1.em1 - dpi_1.em2;
			pfinal = sqrt((s - r__1 * r__1) * (s - r__2 * r__2)) / 2.f / *srt;
			dmelangle_(&pxn, &pyn, &pzn, &pfinal);
			rotate_(px, py, pz, &pxn, &pyn, &pzn);
			r__1 = cc_1.e[ideut - 1];
			r__2 = pxn;
			r__3 = pyn;
			r__4 = pzn;
			edcm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
			pdbeta = pxn * bg_1.betax + pyn * bg_1.betay + pzn * bg_1.betaz;
			transf = bg_1.gamma * (bg_1.gamma * pdbeta / (bg_1.gamma + 1.f) +
								   edcm);
			pt1d = bg_1.betax * transf + pxn;
			pt2d = bg_1.betay * transf + pyn;
			pt3d = bg_1.betaz * transf + pzn;
			bb_1.p[ideut * 3 - 3] = pt1d;
			bb_1.p[ideut * 3 - 2] = pt2d;
			bb_1.p[ideut * 3 - 1] = pt3d;
			*iblock = 504;
			leadng_1.px1 = bb_1.p[*i1 * 3 - 3];
			leadng_1.py1 = bb_1.p[*i1 * 3 - 2];
			leadng_1.pz1 = bb_1.p[*i1 * 3 - 1];
			ee_1.id[*i1 - 1] = 2;
			ee_1.id[*i2 - 1] = 2;
			aa_1.r__[ideut * 3 - 3] = aa_1.r__[idm * 3 - 3];
			aa_1.r__[ideut * 3 - 2] = aa_1.r__[idm * 3 - 2];
			aa_1.r__[ideut * 3 - 1] = aa_1.r__[idm * 3 - 1];
		}
		else
		{
			if (*ianti == 0)
			{
				s_wsle(&io___2359);
				do_lio(&c__9, &c__1, "  d+", (ftnlen)4);
				do_lio(&c__3, &c__1, (char *)&lbm, (ftnlen)sizeof(int));
				do_lio(&c__9, &c__1, " ->BB (pert d destrn) @nt=", (ftnlen)26);
				do_lio(&c__3, &c__1, (char *)&(*nt), (ftnlen)sizeof(int));
				do_lio(&c__9, &c__1, " @prob=", (ftnlen)7);
				do_lio(&c__4, &c__1, (char *)&dpert_1.dpertp[ideut - 1], (ftnlen)sizeof(float));
				e_wsle();
			}
			else
			{
				s_wsle(&io___2360);
				do_lio(&c__9, &c__1, "  d+", (ftnlen)4);
				do_lio(&c__3, &c__1, (char *)&lbm, (ftnlen)sizeof(int));
				do_lio(&c__9, &c__1, " ->BB (pert dbar destrn) @nt=", (ftnlen)29);
				do_lio(&c__3, &c__1, (char *)&(*nt), (ftnlen)sizeof(int));
				do_lio(&c__9, &c__1, " @prob=", (ftnlen)7);
				do_lio(&c__4, &c__1, (char *)&dpert_1.dpertp[ideut - 1], (ftnlen)sizeof(float));
				e_wsle();
			}
			cc_1.e[ideut - 1] = 0.f;
			*iblock = 502;
		}
		return 0;
	}

	*iblock = 502;
	x1 = ranart_(&rndf77_1.nseed);
	if (x1 <= dpisig_1.sdmnn / *sig)
	{
		lbb1 = dpifsl_1.lbnn1;
		lbb2 = dpifsl_1.lbnn2;
		xmb1 = dpifsm_1.xmnn1;
		xmb2 = dpifsm_1.xmnn2;
	}
	else if (x1 <= (dpisig_1.sdmnn + dpisig_1.sdmnd) / *sig)
	{
		lbb1 = dpifsl_1.lbnd1;
		lbb2 = dpifsl_1.lbnd2;
		xmb1 = dpifsm_1.xmnd1;
		xmb2 = dpifsm_1.xmnd2;
	}
	else if (x1 <= (dpisig_1.sdmnn + dpisig_1.sdmnd + dpisig_1.sdmns) / *sig)
	{
		lbb1 = dpifsl_1.lbns1;
		lbb2 = dpifsl_1.lbns2;
		xmb1 = dpifsm_1.xmns1;
		xmb2 = dpifsm_1.xmns2;
	}
	else if (x1 <= (dpisig_1.sdmnn + dpisig_1.sdmnd + dpisig_1.sdmns +
					dpisig_1.sdmnp) /
					   *sig)
	{
		lbb1 = dpifsl_1.lbnp1;
		lbb2 = dpifsl_1.lbnp2;
		xmb1 = dpifsm_1.xmnp1;
		xmb2 = dpifsm_1.xmnp2;
	}
	else if (x1 <= (dpisig_1.sdmnn + dpisig_1.sdmnd + dpisig_1.sdmns +
					dpisig_1.sdmnp + dpisig_1.sdmdd) /
					   *sig)
	{
		lbb1 = dpifsl_1.lbdd1;
		lbb2 = dpifsl_1.lbdd2;
		xmb1 = dpifsm_1.xmdd1;
		xmb2 = dpifsm_1.xmdd2;
	}
	else if (x1 <= (dpisig_1.sdmnn + dpisig_1.sdmnd + dpisig_1.sdmns +
					dpisig_1.sdmnp + dpisig_1.sdmdd + dpisig_1.sdmds) /
					   *sig)
	{
		lbb1 = dpifsl_1.lbds1;
		lbb2 = dpifsl_1.lbds2;
		xmb1 = dpifsm_1.xmds1;
		xmb2 = dpifsm_1.xmds2;
	}
	else if (x1 <= (dpisig_1.sdmnn + dpisig_1.sdmnd + dpisig_1.sdmns +
					dpisig_1.sdmnp + dpisig_1.sdmdd + dpisig_1.sdmds + dpisig_1.sdmdp) /
					   *sig)
	{
		lbb1 = dpifsl_1.lbdp1;
		lbb2 = dpifsl_1.lbdp2;
		xmb1 = dpifsm_1.xmdp1;
		xmb2 = dpifsm_1.xmdp2;
	}
	else if (x1 <= (dpisig_1.sdmnn + dpisig_1.sdmnd + dpisig_1.sdmns +
					dpisig_1.sdmnp + dpisig_1.sdmdd + dpisig_1.sdmds + dpisig_1.sdmdp + dpisig_1.sdmss) /
					   *sig)
	{
		lbb1 = dpifsl_1.lbss1;
		lbb2 = dpifsl_1.lbss2;
		xmb1 = dpifsm_1.xmss1;
		xmb2 = dpifsm_1.xmss2;
	}
	else if (x1 <= (dpisig_1.sdmnn + dpisig_1.sdmnd + dpisig_1.sdmns +
					dpisig_1.sdmnp + dpisig_1.sdmdd + dpisig_1.sdmds + dpisig_1.sdmdp + dpisig_1.sdmss + dpisig_1.sdmsp) /
					   *sig)
	{
		lbb1 = dpifsl_1.lbsp1;
		lbb2 = dpifsl_1.lbsp2;
		xmb1 = dpifsm_1.xmsp1;
		xmb2 = dpifsm_1.xmsp2;
	}
	else if (x1 <= (dpisig_1.sdmnn + dpisig_1.sdmnd + dpisig_1.sdmns +
					dpisig_1.sdmnp + dpisig_1.sdmdd + dpisig_1.sdmds + dpisig_1.sdmdp + dpisig_1.sdmss + dpisig_1.sdmsp + dpisig_1.sdmpp) /
					   *sig)
	{
		lbb1 = dpifsl_1.lbpp1;
		lbb2 = dpifsl_1.lbpp2;
		xmb1 = dpifsm_1.xmpp1;
		xmb2 = dpifsm_1.xmpp2;
	}
	else
	{
		lbb1 = leadng_1.lb1;
		lbb2 = dpi_1.lb2;
		xmb1 = leadng_1.em1;
		xmb2 = dpi_1.em2;
		*iblock = 504;
	}
	ee_1.lb[*i1 - 1] = lbb1;
	cc_1.e[*i1 - 1] = xmb1;
	ee_1.lb[*i2 - 1] = lbb2;
	cc_1.e[*i2 - 1] = xmb2;
	leadng_1.lb1 = ee_1.lb[*i1 - 1];
	dpi_1.lb2 = ee_1.lb[*i2 - 1];
	r__1 = xmb1 + xmb2;
	r__2 = xmb1 - xmb2;
	pfinal = sqrt((s - r__1 * r__1) * (s - r__2 * r__2)) / 2.f / *srt;

	if (*iblock == 502)
	{
		dmangle_(&pxn, &pyn, &pzn, nt, ianti, &pfinal, &lbm);
	}
	else if (*iblock == 504)
	{
		if (*ianti == 0)
		{
			s_wsle(&io___2365);
			do_lio(&c__9, &c__1, " d+", (ftnlen)3);
			do_lio(&c__3, &c__1, (char *)&lbm, (ftnlen)sizeof(int));
			do_lio(&c__9, &c__1, " (regular d M elastic) @evt#", (ftnlen)28);
			do_lio(&c__3, &c__1, (char *)&arevt_1.iaevt, (ftnlen)sizeof(int));
			do_lio(&c__9, &c__1, " @nt=", (ftnlen)5);
			do_lio(&c__3, &c__1, (char *)&(*nt), (ftnlen)sizeof(int));
			do_lio(&c__9, &c__1, " lb1,2=", (ftnlen)7);
			do_lio(&c__3, &c__1, (char *)&leadng_1.lb1, (ftnlen)sizeof(int));
			do_lio(&c__3, &c__1, (char *)&dpi_1.lb2, (ftnlen)sizeof(int));
			e_wsle();
		}
		else
		{
			s_wsle(&io___2366);
			do_lio(&c__9, &c__1, " d+", (ftnlen)3);
			do_lio(&c__3, &c__1, (char *)&lbm, (ftnlen)sizeof(int));
			do_lio(&c__9, &c__1, " (regular dbar M elastic) @evt#", (ftnlen)31);
			do_lio(&c__3, &c__1, (char *)&arevt_1.iaevt, (ftnlen)sizeof(int));
			do_lio(&c__9, &c__1, " @nt=", (ftnlen)5);
			do_lio(&c__3, &c__1, (char *)&(*nt), (ftnlen)sizeof(int));
			do_lio(&c__9, &c__1, " lb1,2=", (ftnlen)7);
			do_lio(&c__3, &c__1, (char *)&leadng_1.lb1, (ftnlen)sizeof(int));
			do_lio(&c__3, &c__1, (char *)&dpi_1.lb2, (ftnlen)sizeof(int));
			e_wsle();
		}
		dmelangle_(&pxn, &pyn, &pzn, &pfinal);
	}
	else
	{
		s_wsle(&io___2367);
		do_lio(&c__9, &c__1, "Wrong iblock number in crdmbb()", (ftnlen)31);
		e_wsle();
		s_stop("", (ftnlen)0);
	}
	rotate_(px, py, pz, &pxn, &pyn, &pzn);
	r__1 = cc_1.e[*i1 - 1];
	r__2 = pxn;
	r__3 = pyn;
	r__4 = pzn;
	e1cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	p1beta = pxn * bg_1.betax + pyn * bg_1.betay + pzn * bg_1.betaz;
	transf = bg_1.gamma * (bg_1.gamma * p1beta / (bg_1.gamma + 1.f) + e1cm);
	pt1i1 = bg_1.betax * transf + pxn;
	pt2i1 = bg_1.betay * transf + pyn;
	pt3i1 = bg_1.betaz * transf + pzn;

	bb_1.p[*i1 * 3 - 3] = pt1i1;
	bb_1.p[*i1 * 3 - 2] = pt2i1;
	bb_1.p[*i1 * 3 - 1] = pt3i1;
	r__1 = cc_1.e[*i2 - 1];
	r__2 = pxn;
	r__3 = pyn;
	r__4 = pzn;
	e2cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	p2beta = -pxn * bg_1.betax - pyn * bg_1.betay - pzn * bg_1.betaz;
	transf = bg_1.gamma * (bg_1.gamma * p2beta / (bg_1.gamma + 1.f) + e2cm);
	pt1i2 = bg_1.betax * transf - pxn;
	pt2i2 = bg_1.betay * transf - pyn;
	pt3i2 = bg_1.betaz * transf - pzn;

	bb_1.p[*i2 * 3 - 3] = pt1i2;
	bb_1.p[*i2 * 3 - 2] = pt2i2;
	bb_1.p[*i2 * 3 - 1] = pt3i2;

	leadng_1.px1 = bb_1.p[*i1 * 3 - 3];
	leadng_1.py1 = bb_1.p[*i1 * 3 - 2];
	leadng_1.pz1 = bb_1.p[*i1 * 3 - 1];
	leadng_1.em1 = cc_1.e[*i1 - 1];
	dpi_1.em2 = cc_1.e[*i2 - 1];
	ee_1.id[*i1 - 1] = 2;
	ee_1.id[*i2 - 1] = 2;
	return 0;
}

int dmangle_(float *pxn, float *pyn, float *pzn, int *nt,
			 int *ianti, float *pfinal, int *lbm)
{
	float r__1;

	double sqrt(double), cos(double), sin(double);
	int s_wsle(cilist *), do_lio(int *, int *, char *, ftnlen),
		e_wsle(void);

	static float c1, s1, t1, ct1, st1;
	extern double ranart_(int *);

	static cilist io___2383 = {0, 91, 0, 0, 0};
	static cilist io___2384 = {0, 91, 0, 0, 0};

	c1 = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
	t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
	r__1 = c1;
	s1 = sqrt(1.f - r__1 * r__1);
	ct1 = cos(t1);
	st1 = sin(t1);
	*pzn = *pfinal * c1;
	*pxn = *pfinal * s1 * ct1;
	*pyn = *pfinal * s1 * st1;
	if (*ianti == 0)
	{
		s_wsle(&io___2383);
		do_lio(&c__9, &c__1, " d+", (ftnlen)3);
		do_lio(&c__3, &c__1, (char *)&(*lbm), (ftnlen)sizeof(int));
		do_lio(&c__9, &c__1, " ->BB (regular d destrn) @evt#", (ftnlen)30);
		do_lio(&c__3, &c__1, (char *)&arevt_1.iaevt, (ftnlen)sizeof(int));
		do_lio(&c__9, &c__1, " @nt=", (ftnlen)5);
		do_lio(&c__3, &c__1, (char *)&(*nt), (ftnlen)sizeof(int));
		do_lio(&c__9, &c__1, " lb1,2=", (ftnlen)7);
		do_lio(&c__3, &c__1, (char *)&leadng_1.lb1, (ftnlen)sizeof(int));
		do_lio(&c__3, &c__1, (char *)&dpi_1.lb2, (ftnlen)sizeof(int));
		e_wsle();
	}
	else
	{
		s_wsle(&io___2384);
		do_lio(&c__9, &c__1, " d+", (ftnlen)3);
		do_lio(&c__3, &c__1, (char *)&(*lbm), (ftnlen)sizeof(int));
		do_lio(&c__9, &c__1, " ->BB (regular dbar destrn) @evt#", (ftnlen)33);
		do_lio(&c__3, &c__1, (char *)&arevt_1.iaevt, (ftnlen)sizeof(int));
		do_lio(&c__9, &c__1, " @nt=", (ftnlen)5);
		do_lio(&c__3, &c__1, (char *)&(*nt), (ftnlen)sizeof(int));
		do_lio(&c__9, &c__1, " lb1,2=", (ftnlen)7);
		do_lio(&c__3, &c__1, (char *)&leadng_1.lb1, (ftnlen)sizeof(int));
		do_lio(&c__3, &c__1, (char *)&dpi_1.lb2, (ftnlen)sizeof(int));
		e_wsle();
	}

	return 0;
}

int dmelangle_(float *pxn, float *pyn, float *pzn, float *pfinal)
{
	float r__1;

	double sqrt(double), cos(double), sin(double);

	static float c1, s1, t1, ct1, st1;
	extern double ranart_(int *);

	c1 = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
	t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
	r__1 = c1;
	s1 = sqrt(1.f - r__1 * r__1);
	ct1 = cos(t1);
	st1 = sin(t1);
	*pzn = *pfinal * c1;
	*pxn = *pfinal * s1 * ct1;
	*pyn = *pfinal * s1 * st1;
	return 0;
}

int sdbelastic_(float *srt, float *sdb)
{
	float r__1;

	static float s, threshold, snew;
	extern double fdbel_(float *);
	static float sdbel;

	*sdb = 0.f;
	sdbel = 0.f;
	return 0;

	if (*srt <= leadng_1.em1 + dpi_1.em2)
	{
		return 0;
	}
	r__1 = *srt;
	s = r__1 * r__1;
	if (para8_1.idxsec == 1 || para8_1.idxsec == 3)
	{
		sdbel = fdbel_(&s);
	}
	else if (para8_1.idxsec == 2 || para8_1.idxsec == 4)
	{
		threshold = leadng_1.em1 + dpi_1.em2;
		r__1 = *srt - threshold + 2.012f;
		snew = r__1 * r__1;
		sdbel = fdbel_(&snew);
	}
	*sdb = sdbel;
	return 0;
}

int crdbel_(float *px, float *py, float *pz, float *srt, int *i1, int *i2, int *iblock, int *ntag, float *sig, int *nt, int *ianti)
{
	float r__1, r__2, r__3, r__4;

	int s_wsle(cilist *), do_lio(int *, int *, char *, ftnlen),
		e_wsle(void);
	double sqrt(double);

	static float s;
	extern int dbelangle_(float *, float *, float *, float *);
	static int idb, lbb;
	static float pxn, pyn, pzn, e1cm, e2cm, pt1d, pt2d, pt3d, edcm, pt1i1,
		pt2i1, pt3i1, pt1i2, pt2i2, pt3i2;
	static int ideut;
	static float p1beta, p2beta, pdbeta, pfinal, transf;
	extern int rotate_(float *, float *, float *, float *, float *, float *);

	static cilist io___2398 = {0, 91, 0, 0, 0};
	static cilist io___2399 = {0, 91, 0, 0, 0};
	static cilist io___2410 = {0, 91, 0, 0, 0};
	static cilist io___2411 = {0, 91, 0, 0, 0};

	*iblock = 0;
	*ntag = 0;
	leadng_1.em1 = cc_1.e[*i1 - 1];
	dpi_1.em2 = cc_1.e[*i2 - 1];
	r__1 = *srt;
	s = r__1 * r__1;
	if (*sig <= 0.f)
	{
		return 0;
	}
	*iblock = 503;

	if (abs(leadng_1.lb1) == 42)
	{
		ideut = *i1;
		lbb = dpi_1.lb2;
		idb = *i2;
	}
	else
	{
		ideut = *i2;
		lbb = leadng_1.lb1;
		idb = *i1;
	}
	if ((para8_1.idpert == 1 || para8_1.idpert == 2) && dpert_1.dpertp[ideut - 1] != 1.f)
	{
		if (*ianti == 0)
		{
			s_wsle(&io___2398);
			do_lio(&c__9, &c__1, "  d+", (ftnlen)4);
			do_lio(&c__3, &c__1, (char *)&lbb, (ftnlen)sizeof(int));
			do_lio(&c__9, &c__1, " (pert d B elastic) @nt=", (ftnlen)24);
			do_lio(&c__3, &c__1, (char *)&(*nt), (ftnlen)sizeof(int));
			do_lio(&c__9, &c__1, " @prob=", (ftnlen)7);
			do_lio(&c__4, &c__1, (char *)&dpert_1.dpertp[ideut - 1], (ftnlen)sizeof(float));
			do_lio(&c__4, &c__1, (char *)&bb_1.p[idb * 3 - 3], (ftnlen)sizeof(float));
			do_lio(&c__4, &c__1, (char *)&bb_1.p[idb * 3 - 2], (ftnlen)sizeof(float));
			do_lio(&c__4, &c__1, (char *)&bb_1.p[ideut * 3 - 3], (ftnlen)sizeof(float));
			do_lio(&c__4, &c__1, (char *)&bb_1.p[ideut * 3 - 2], (ftnlen)sizeof(float));
			e_wsle();
		}
		else
		{
			s_wsle(&io___2399);
			do_lio(&c__9, &c__1, "  d+", (ftnlen)4);
			do_lio(&c__3, &c__1, (char *)&lbb, (ftnlen)sizeof(int));
			do_lio(&c__9, &c__1, " (pert dbar Bbar elastic) @nt=", (ftnlen)30);
			do_lio(&c__3, &c__1, (char *)&(*nt), (ftnlen)sizeof(int));
			do_lio(&c__9, &c__1, " @prob=", (ftnlen)7);
			do_lio(&c__4, &c__1, (char *)&dpert_1.dpertp[ideut - 1], (ftnlen)sizeof(float));
			do_lio(&c__4, &c__1, (char *)&bb_1.p[idb * 3 - 3], (ftnlen)sizeof(float));
			do_lio(&c__4, &c__1, (char *)&bb_1.p[idb * 3 - 2], (ftnlen)sizeof(float));
			do_lio(&c__4, &c__1, (char *)&bb_1.p[ideut * 3 - 3], (ftnlen)sizeof(float));
			do_lio(&c__4, &c__1, (char *)&bb_1.p[ideut * 3 - 2], (ftnlen)sizeof(float));
			e_wsle();
		}
		r__1 = leadng_1.em1 + dpi_1.em2;
		r__2 = leadng_1.em1 - dpi_1.em2;
		pfinal = sqrt((s - r__1 * r__1) * (s - r__2 * r__2)) / 2.f / *srt;
		dbelangle_(&pxn, &pyn, &pzn, &pfinal);
		rotate_(px, py, pz, &pxn, &pyn, &pzn);
		r__1 = cc_1.e[ideut - 1];
		r__2 = pxn;
		r__3 = pyn;
		r__4 = pzn;
		edcm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
		pdbeta = pxn * bg_1.betax + pyn * bg_1.betay + pzn * bg_1.betaz;
		transf = bg_1.gamma * (bg_1.gamma * pdbeta / (bg_1.gamma + 1.f) +
							   edcm);
		pt1d = bg_1.betax * transf + pxn;
		pt2d = bg_1.betay * transf + pyn;
		pt3d = bg_1.betaz * transf + pzn;
		bb_1.p[ideut * 3 - 3] = pt1d;
		bb_1.p[ideut * 3 - 2] = pt2d;
		bb_1.p[ideut * 3 - 1] = pt3d;
		leadng_1.px1 = bb_1.p[*i1 * 3 - 3];
		leadng_1.py1 = bb_1.p[*i1 * 3 - 2];
		leadng_1.pz1 = bb_1.p[*i1 * 3 - 1];
		ee_1.id[*i1 - 1] = 2;
		ee_1.id[*i2 - 1] = 2;
		aa_1.r__[ideut * 3 - 3] = aa_1.r__[idb * 3 - 3];
		aa_1.r__[ideut * 3 - 2] = aa_1.r__[idb * 3 - 2];
		aa_1.r__[ideut * 3 - 1] = aa_1.r__[idb * 3 - 1];
		return 0;
	}

	if (*ianti == 0)
	{
		s_wsle(&io___2410);
		do_lio(&c__9, &c__1, " d+", (ftnlen)3);
		do_lio(&c__3, &c__1, (char *)&lbb, (ftnlen)sizeof(int));
		do_lio(&c__9, &c__1, " (regular d B elastic) @evt#", (ftnlen)28);
		do_lio(&c__3, &c__1, (char *)&arevt_1.iaevt, (ftnlen)sizeof(int));
		do_lio(&c__9, &c__1, " @nt=", (ftnlen)5);
		do_lio(&c__3, &c__1, (char *)&(*nt), (ftnlen)sizeof(int));
		do_lio(&c__9, &c__1, " lb1,2=", (ftnlen)7);
		do_lio(&c__3, &c__1, (char *)&leadng_1.lb1, (ftnlen)sizeof(int));
		do_lio(&c__3, &c__1, (char *)&dpi_1.lb2, (ftnlen)sizeof(int));
		e_wsle();
	}
	else
	{
		s_wsle(&io___2411);
		do_lio(&c__9, &c__1, " d+", (ftnlen)3);
		do_lio(&c__3, &c__1, (char *)&lbb, (ftnlen)sizeof(int));
		do_lio(&c__9, &c__1, " (regular dbar Bbar elastic) @evt#", (ftnlen)34);
		do_lio(&c__3, &c__1, (char *)&arevt_1.iaevt, (ftnlen)sizeof(int));
		do_lio(&c__9, &c__1, " @nt=", (ftnlen)5);
		do_lio(&c__3, &c__1, (char *)&(*nt), (ftnlen)sizeof(int));
		do_lio(&c__9, &c__1, " lb1,2=", (ftnlen)7);
		do_lio(&c__3, &c__1, (char *)&leadng_1.lb1, (ftnlen)sizeof(int));
		do_lio(&c__3, &c__1, (char *)&dpi_1.lb2, (ftnlen)sizeof(int));
		e_wsle();
	}
	r__1 = leadng_1.em1 + dpi_1.em2;
	r__2 = leadng_1.em1 - dpi_1.em2;
	pfinal = sqrt((s - r__1 * r__1) * (s - r__2 * r__2)) / 2.f / *srt;
	dbelangle_(&pxn, &pyn, &pzn, &pfinal);
	rotate_(px, py, pz, &pxn, &pyn, &pzn);
	r__1 = cc_1.e[*i1 - 1];
	r__2 = pxn;
	r__3 = pyn;
	r__4 = pzn;
	e1cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	p1beta = pxn * bg_1.betax + pyn * bg_1.betay + pzn * bg_1.betaz;
	transf = bg_1.gamma * (bg_1.gamma * p1beta / (bg_1.gamma + 1.f) + e1cm);
	pt1i1 = bg_1.betax * transf + pxn;
	pt2i1 = bg_1.betay * transf + pyn;
	pt3i1 = bg_1.betaz * transf + pzn;

	bb_1.p[*i1 * 3 - 3] = pt1i1;
	bb_1.p[*i1 * 3 - 2] = pt2i1;
	bb_1.p[*i1 * 3 - 1] = pt3i1;
	r__1 = cc_1.e[*i2 - 1];
	r__2 = pxn;
	r__3 = pyn;
	r__4 = pzn;
	e2cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	p2beta = -pxn * bg_1.betax - pyn * bg_1.betay - pzn * bg_1.betaz;
	transf = bg_1.gamma * (bg_1.gamma * p2beta / (bg_1.gamma + 1.f) + e2cm);
	pt1i2 = bg_1.betax * transf - pxn;
	pt2i2 = bg_1.betay * transf - pyn;
	pt3i2 = bg_1.betaz * transf - pzn;

	bb_1.p[*i2 * 3 - 3] = pt1i2;
	bb_1.p[*i2 * 3 - 2] = pt2i2;
	bb_1.p[*i2 * 3 - 1] = pt3i2;

	leadng_1.px1 = bb_1.p[*i1 * 3 - 3];
	leadng_1.py1 = bb_1.p[*i1 * 3 - 2];
	leadng_1.pz1 = bb_1.p[*i1 * 3 - 1];
	leadng_1.em1 = cc_1.e[*i1 - 1];
	dpi_1.em2 = cc_1.e[*i2 - 1];
	ee_1.id[*i1 - 1] = 2;
	ee_1.id[*i2 - 1] = 2;
	return 0;
}

double fnndpi_(float *s)
{
	float ret_val, r__1, r__2, r__3;

	double exp(double);

	if (*s <= 4.0481439999999997f)
	{
		ret_val = 0.f;
	}
	else
	{
		r__1 = *s - 4.65f;
		r__2 = *s - 4.65f;
		r__3 = *s - 6.f;
		ret_val = exp(-(r__1 * r__1) / .1f) * 26.f + exp(-(r__2 * r__2) / 2.f) * 4.f + exp(-(r__3 * r__3) / 10.f) * .28f;
	}
	return ret_val;
}

int dbelangle_(float *pxn, float *pyn, float *pzn, float *pfinal)
{
	float r__1;

	double sqrt(double), cos(double), sin(double);

	static float c1, s1, t1, ct1, st1;
	extern double ranart_(int *);

	c1 = 1.f - ranart_(&rndf77_1.nseed) * 2.f;
	t1 = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
	r__1 = c1;
	s1 = sqrt(1.f - r__1 * r__1);
	ct1 = cos(t1);
	st1 = sin(t1);
	*pzn = *pfinal * c1;
	*pxn = *pfinal * s1 * ct1;
	*pyn = *pfinal * s1 * st1;
	return 0;
}

double fdpiel_(float *s)
{
	float ret_val, r__1, r__2;

	double exp(double);

	if (*s <= 4.0481439999999997f)
	{
		ret_val = 0.f;
	}
	else
	{
		r__1 = *s - 4.67f;
		r__2 = *s - 6.25f;
		ret_val = exp(-(r__1 * r__1) / .15f) * 63.f + exp(-(r__2 * r__2) /
														  .3f) *
														  15.f;
	}
	return ret_val;
}

double fdbel_(float *s)
{
	float ret_val, r__1, r__2;

	double exp(double);

	if (*s <= 4.0481439999999997f)
	{
		ret_val = 0.f;
	}
	else
	{
		r__1 = *s - 7.93f;
		r__2 = *s - 7.93f;
		ret_val = exp(-(r__1 * r__1) / .003f) * 2500.f + exp(-(r__2 * r__2) / .1f) * 300.f + 10.f;
	}
	return ret_val;
}
