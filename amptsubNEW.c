#include "f2c.h"
#include <math.h>
struct
{
	float dx, dy, dz, dpx, dpy, dpz;
} gg_;

#define gg_1 gg_

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
	int masspr, massta, iseed, iavoid;
	float dt;
} input1_;

#define input1_1 input1_

struct
{
	int ilab, manyb, ntmax, icoll, insys, ipot, mode, imomen, nfreq,
		icflow, icrho, icou, kpoten, kmul;
} input2_;

#define input2_1 input2_

struct
{
	float plab, elab, zeropt, b0, bi, bm, dencut, cycbox;
} input3_;

#define input3_1 input3_

struct
{
	int iperts;
} imulst_;

#define imulst_1 imulst_

struct
{
	double dpcoal, drcoal, ecritl;
} coal_;

#define coal_1 coal_

struct
{
	int nevent, isoft, isflag, izpc;
} anim_;

#define anim_1 anim_

struct
{
	int ioscar;
} para7_;

#define para7_1 para7_

struct arprnt_1_
{
	float arpar1[100];
	int iapar2[50];
	float arint1[100];
	int iaint2[50];
};

#define arprnt_1 (*(struct arprnt_1_ *)&arprnt_)

struct
{
	int itypar[150001];
	float gxar[150001], gyar[150001], gzar[150001], ftar[150001], pxar[150001],
		pyar[150001], pzar[150001], pear[150001], xmar[150001];
} arprc_;

#define arprc_1 arprc_

struct
{
	double smearp, smearh;
} smearz_;

#define smearz_1 smearz_

struct
{
	int nattzp;
} nzpc_;

#define nzpc_1 nzpc_

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
	int nseed;
} rndf77_;

#define rndf77_1 rndf77_

struct
{
	int idpert, npertd, idxsec;
} para8_;

#define para8_1 para8_

struct
{
	int multi1[1];
} arerc1_;

#define arerc1_1 arerc1_

struct
{
	int ityp1[150001];
	float gx1[150001], gy1[150001], gz1[150001], ft1[150001], px1[150001], py1[150001], pz1[150001], ee1[150001], xm1[150001];
} arprc1_;

#define arprc1_1 arprc1_

struct
{
	float tfdcy[150001], tfdpi[150001], tft[150001];
} tdecay_;

#define tdecay_1 tdecay_

struct
{
	float dpertt[150001], dpertp[150001], dplast[150001],
		dpdcy[150001], dpdpi[150001], dpt[150001], dpp1[150001], dppion[150001];
} dpert_;

#define dpert_1 dpert_

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
	int id[150001], lb[150001];
} ee_;

#define ee_1 ee_

struct
{
	float betax, betay, betaz, gamma;
} bg_;

#define bg_1 bg_

struct
{
	int nnn;
} nn_;

#define nn_1 nn_

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
	float dy1ntb[50], dy1ntp[50], dy1hm[50], dy1kp[50], dy1km[50], dy1k0s[50],
		dy1la[50], dy1lb[50], dy1phi[50], dm1pip[50], dm1pim[50], dmt1pr[50], dmt1pb[50], dmt1kp[50], dm1km[50], dm1k0s[50], dmt1la[50],
		dmt1lb[50], dy1msn[50], dy1pip[50], dy1pim[50], dy1pi0[50], dy1pr[50], dy1pb[50], dy1neg[50], dy1ch[50], de1neg[50], de1ch[50];
} arana1_;

#define arana1_1 arana1_

struct
{
	float dy2ntb[50], dy2ntp[50], dy2hm[50], dy2kp[50], dy2km[50], dy2k0s[50],
		dy2la[50], dy2lb[50], dy2phi[50], dm2pip[50], dm2pim[50], dmt2pr[50], dmt2pb[50], dmt2kp[50], dm2km[50], dm2k0s[50], dmt2la[50],
		dmt2lb[50], dy2msn[50], dy2pip[50], dy2pim[50], dy2pi0[50], dy2pr[50], dy2pb[50], dy2neg[50], dy2ch[50], de2neg[50], de2ch[50];
} arana2_;

#define arana2_1 arana2_

struct
{
	int mul;
} para1_;

#define para1_1 para1_

struct
{
	float yp[900], yt[900];
} hjcrdn_;

#define hjcrdn_1 hjcrdn_

struct
{
	int npj[300], kfpj[150000];
	float pjpx[150000], pjpy[150000], pjpz[150000], pjpe[150000], pjpm[150000];
	int ntj[300], kftj[150000];
	float pjtx[150000], pjty[150000], pjtz[150000], pjte[150000], pjtm[150000];
} hjjet1_;

#define hjjet1_1 hjjet1_

struct
{
	int nsg, njsg[150001], iasg[450003], k1sg[15000100], k2sg[15000100];
	float pxsg[15000100], pysg[15000100], pzsg[15000100], pesg[15000100], pmsg[15000100];
} hjjet2_;

#define hjjet2_1 hjjet2_

struct
{
	double gx0[400001], gy0[400001], gz0[400001], ft0[400001], px0[400001], py0[400001], pz0[400001], e0[400001], xmass0[400001];
	int ityp0[400001];
} prec1_;

#define prec1_1 prec1_

struct
{
	int iaevt, iarun;
} arevt_;

#define arevt_1 arevt_

struct
{
	int iout;
} arout_;

#define arout_1 arout_

struct
{
	double gx5[400001], gy5[400001], gz5[400001], ft5[400001], px5[400001], py5[400001], pz5[400001], e5[400001], xmass5[400001];
	int ityp5[400001];
} prec2_;

#define prec2_1 prec2_

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
	double pxsgs[450003], pysgs[450003], pzsgs[450003], pesgs[450003], pmsgs[450003], gxsgs[450003], gysgs[450003], gzsgs[450003], ftsgs[450003];
	int k1sgs[450003], k2sgs[450003], njsgs[150001];
} soft_;

#define soft_1 soft_

struct
{
	double v2i, eti, xmulti, v2mi, s2mi, xmmult, v2bi, s2bi, xbmult;
} iflow_;

#define iflow_1 iflow_

struct
{
	float v2f, etf, xmultf, v2fpi, xmulpi;
} fflow_;

#define fflow_1 fflow_

struct
{
	int np[150001];
} strg_;

#define strg_1 strg_

struct
{
	double gxfrz[400001], gyfrz[400001], gzfrz[400001], ftfrz[400001],
		pxfrz[400001], pyfrz[400001], pzfrz[400001], efrz[400001], xmfrz[400001], tfrz[302];
	int ifrz[400001], idfrz[400001], itlast;
} frzprc_;

#define frzprc_1 frzprc_

struct
{
	float e_1[100];
	int e_2[50];
	float e_3[100];
	int e_4[50];
} arprnt_ = {1.19f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
			 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
			 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
			 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
			 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
			 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
			 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
			 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 3, 0, 0, 0,
			 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			 0, 0, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
			 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
			 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
			 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
			 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
			 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
			 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f,
			 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0, 0, 0, 0, 0,
			 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			 0};

static int c__9 = 9;
static int c__1 = 1;
static int c__3 = 3;
static int c__4 = 4;
static int c_b96 = 150001;
static double c_b112 = 3.5;
static double c_b114 = 1.86;
static double c_b116 = 2.74;
static double c_b117 = -1.8;
static double c_b119 = -.9;
static double c_b121 = -1.3;
static double c_b122 = -2.3;
static double c_b126 = -2.6;
static double c_b129 = -1.83;
static double c_b130 = -2.09;
static int c__0 = 0;
int artset_(void)
{
	float r__1;
	olist o__1;
	double sqrt(double);
	int f_open(olist *);
	static int iplab;
	coal_1.ecritl = 1.;
	input1_1.massta = 1;
	zz_1.zta = 1;
	input1_1.masspr = 1;
	zz_1.zpr = 1;
	input3_1.plab = 14.6f;
	iplab = 2;
	if (iplab == 2)
	{

		r__1 = input3_1.plab;
		input3_1.elab = sqrt(r__1 * r__1 + .88040689000000005f) - .9383f;
	}
	else
	{
		input3_1.elab = input3_1.plab;
	}
	input3_1.elab *= 1e3f;

	input3_1.zeropt = 0.f;

	input1_1.iseed = 700721;

	imulst_1.iperts = 0;

	input2_1.manyb = 1;
	input3_1.b0 = 1.f;
	input3_1.bi = 0.f;
	input3_1.bm = 0.f;

	input2_1.icoll = -1;

	run_1.num = 1;

	input2_1.insys = 1;

	input2_1.ipot = 3;

	input2_1.mode = 0;
	if (input2_1.icoll == -1)
	{
		input2_1.ipot = 0;
	}

	gg_1.dx = 2.73f;
	gg_1.dy = 2.73f;
	gg_1.dz = 2.73f;

	gg_1.dpx = .6f;
	gg_1.dpy = .6f;
	gg_1.dpz = .6f;

	input1_1.iavoid = 1;

	input2_1.imomen = 1;

	if (input2_1.icoll == -1)
	{
		input2_1.imomen = 3;
	}

	input2_1.nfreq = 10;

	input2_1.icflow = 0;

	input2_1.icrho = 0;

	input2_1.icou = 0;

	input2_1.kpoten = 0;
	input2_1.kmul = 1;

	input3_1.dencut = 15.f;

	input3_1.cycbox = 0.f;

	if ((anim_1.isoft == 4 || anim_1.isoft == 5) && para7_1.ioscar == 2)
	{
		o__1.oerr = 0;
		o__1.ounit = 92;
		o__1.ofnmlen = 25;
		o__1.ofnm = "ana/initial_parton_sm.dat";
		o__1.orl = 0;
		o__1.osta = "UNKNOWN";
		o__1.oacc = 0;
		o__1.ofm = 0;
		o__1.oblnk = 0;
		f_open(&o__1);
	}
	return 0;
}

int arini_(void)
{

	int s_wsle(cilist *), do_lio(int *, int *, char *, int),
		e_wsle(void);
	int s_stop(char *, int);

	static int iflg;
	extern int arini1_(void), artord_(void);

	static cilist io___3 = {0, 6, 0, 0, 0};

	iflg = arprnt_1.iapar2[0];
	switch (iflg)
	{
	case 1:
		return 0;
	case 2:
		return 0;
	case 3:
		arini1_();
	}

	s_wsle(&io___3);
	do_lio(&c__9, &c__1, "IAPAR2(1) must be 1, 2, or 3", (int)28);
	e_wsle();
	s_stop("", (int)0);
	artord_();
	return 0;
}

int arini1_(void)
{

	int i__1;
	float r__1;
	olist o__1;

	int f_open(olist *), s_wsle(cilist *), do_lio(int *, int *, char *, int), e_wsle(void);
	double log(double), cosh(double), sinh(double);
	static int i__, np;
	static float vx, vy, rap, tau0, taui;
	double ranart_(int *);
	static float zsmear;
	static cilist io___7 = {0, 6, 0, 0, 0};
	static cilist io___8 = {0, 6, 0, 0, 0};
	static cilist io___9 = {0, 6, 0, 0, 0};
	static cilist io___10 = {0, 6, 0, 0, 0};
	static cilist io___11 = {0, 6, 0, 0, 0};
	static cilist io___15 = {0, 6, 0, 0, 0};
	static cilist io___16 = {0, 6, 0, 0, 0};
	static cilist io___17 = {0, 6, 0, 0, 0};
	static cilist io___18 = {0, 6, 0, 0, 0};
	static cilist io___19 = {0, 6, 0, 0, 0};

	if (para8_1.idpert == 1 || para8_1.idpert == 2)
	{
		o__1.oerr = 0;
		o__1.ounit = 90;
		o__1.ofnmlen = 17;
		o__1.ofnm = "ana/ampt_pert.dat";
		o__1.orl = 0;
		o__1.osta = "UNKNOWN";
		o__1.oacc = 0;
		o__1.ofm = 0;
		o__1.oblnk = 0;
		f_open(&o__1);
	}

	tau0 = arprnt_1.arpar1[0];
	np = arprnt_1.iaint2[0];

	if (anim_1.isoft == 3 || anim_1.isoft == 4 || anim_1.isoft == 5)
	{
		if (np <= nzpc_1.nattzp)
		{
			return 0;
		}
		i__1 = np;
		for (i__ = nzpc_1.nattzp + 1; i__ <= i__1; ++i__)
		{
			if ((r__1 = arprc_1.pzar[i__ - 1], dabs(r__1)) >= arprc_1.pear[i__ - 1])
			{
				s_wsle(&io___7);
				do_lio(&c__9, &c__1, " IN ARINI1", (int)10);
				e_wsle();
				s_wsle(&io___8);
				do_lio(&c__9, &c__1, "ABS(PZ) .GE. EE for particle ", (int)29);
				do_lio(&c__3, &c__1, (char *)&i__, (int)sizeof(int));
				e_wsle();
				s_wsle(&io___9);
				do_lio(&c__9, &c__1, " FLAV = ", (int)8);
				do_lio(&c__3, &c__1, (char *)&arprc_1.itypar[i__ - 1], (int)sizeof(int));
				do_lio(&c__9, &c__1, " PX = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&arprc_1.pxar[i__ - 1], (int)sizeof(float));
				do_lio(&c__9, &c__1, " PY = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&arprc_1.pyar[i__ - 1], (int)sizeof(float));
				e_wsle();
				s_wsle(&io___10);
				do_lio(&c__9, &c__1, " PZ = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&arprc_1.pzar[i__ - 1], (int)sizeof(float));
				do_lio(&c__9, &c__1, " EE = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&arprc_1.pear[i__ - 1], (int)sizeof(float));
				e_wsle();
				s_wsle(&io___11);
				do_lio(&c__9, &c__1, " XM = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&arprc_1.xmar[i__ - 1], (int)sizeof(float));
				e_wsle();
				rap = 1e6f;
                vx = arprc_1.pxar[i__ - 1] / arprc_1.pear[i__ - 1];
                vy = arprc_1.pyar[i__ - 1] / arprc_1.pear[i__ - 1];
                arprc_1.ftar[i__ - 1] = tau0 * cosh(rap);
                arprc_1.gxar[i__ - 1] += vx * arprc_1.ftar[i__ - 1];
                arprc_1.gyar[i__ - 1] += vy * arprc_1.ftar[i__ - 1];
                arprc_1.gzar[i__ - 1] = tau0 * sinh(rap);
				
			}
			rap = log((arprc_1.pear[i__ - 1] + arprc_1.pzar[i__ - 1]) / (arprc_1.pear[i__ - 1] - arprc_1.pzar[i__ - 1])) * .5f;
		}
		return 0;
	}

	i__1 = np;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		if ((r__1 = arprc_1.pzar[i__ - 1], dabs(r__1)) >= arprc_1.pear[i__ -1])
		{
			s_wsle(&io___15);
			do_lio(&c__9, &c__1, " IN ARINI1", (int)10);
			e_wsle();
			s_wsle(&io___16);
			do_lio(&c__9, &c__1, "ABS(PZ) .GE. EE for particle ", (int)29);
			do_lio(&c__3, &c__1, (char *)&i__, (int)sizeof(int));
			e_wsle();
			s_wsle(&io___17);
			do_lio(&c__9, &c__1, " FLAV = ", (int)8);
			do_lio(&c__3, &c__1, (char *)&arprc_1.itypar[i__ - 1], (int)sizeof(int));
			do_lio(&c__9, &c__1, " PX = ", (int)6);
			do_lio(&c__4, &c__1, (char *)&arprc_1.pxar[i__ - 1], (int)sizeof(float));
			do_lio(&c__9, &c__1, " PY = ", (int)6);
			do_lio(&c__4, &c__1, (char *)&arprc_1.pyar[i__ - 1], (int)sizeof(float));
			e_wsle();
			s_wsle(&io___18);
			do_lio(&c__9, &c__1, " PZ = ", (int)6);
			do_lio(&c__4, &c__1, (char *)&arprc_1.pzar[i__ - 1], (int)sizeof(float));
			do_lio(&c__9, &c__1, " EE = ", (int)6);
			do_lio(&c__4, &c__1, (char *)&arprc_1.pear[i__ - 1], (int)sizeof(float));
			e_wsle();
			s_wsle(&io___19);
			do_lio(&c__9, &c__1, " XM = ", (int)6);
			do_lio(&c__4, &c__1, (char *)&arprc_1.xmar[i__ - 1], (int)sizeof(float));
			e_wsle();
			rap = 1e6f;
		}
		rap = log((arprc_1.pear[i__ - 1] + arprc_1.pzar[i__ - 1]) / (arprc_1.pear[i__ - 1] - arprc_1.pzar[i__ - 1])) * .5f;
		vx = arprc_1.pxar[i__ - 1] / arprc_1.pear[i__ - 1];
		vy = arprc_1.pyar[i__ - 1] / arprc_1.pear[i__ - 1];
		taui = arprc_1.ftar[i__ - 1] + tau0;
		arprc_1.ftar[i__ - 1] = taui * cosh(rap);
		arprc_1.gxar[i__ - 1] += vx * tau0 * cosh(rap);
		arprc_1.gyar[i__ - 1] += vy * tau0 * cosh(rap);
		arprc_1.gzar[i__ - 1] = taui * sinh(rap);
		zsmear = (float)smearz_1.smearh * (ranart_(&rndf77_1.nseed) * 2.f -1.f);
		arprc_1.gzar[i__ - 1] += zsmear;
		if (arprc_1.pxar[i__ - 1] == 0.f && arprc_1.pyar[i__ - 1] == 0.f &&
			arprc_1.pear[i__ - 1] * 2 / hparnt_1.hint1[0] > .99f && (arprc_1.itypar[i__ - 1] == 2112 || arprc_1.itypar[i__ - 1] == 2212))
		{

			taui = 1e-20f;
			arprc_1.ftar[i__ - 1] = taui * cosh(rap);
			arprc_1.gzar[i__ - 1] = taui * sinh(rap) + zsmear;
		}
	}
	return 0;
}

int artord_(void)
{

	int i__1;

	static int i__, np;
	static float ee0[150001], ft0[150001], gx0[150001], gy0[150001], gz0[150001], xm0[150001], px0[150001], py0[150001], pz0[150001];
	static int npar, indx[150001], ityp0[150001];
	extern int arindx_(int *, int *, float *, int *);

	npar = 0;
	np = arprnt_1.iaint2[0];
	i__1 = np;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		ityp0[i__ - 1] = arprc_1.itypar[i__ - 1];
		gx0[i__ - 1] = arprc_1.gxar[i__ - 1];
		gy0[i__ - 1] = arprc_1.gyar[i__ - 1];
		gz0[i__ - 1] = arprc_1.gzar[i__ - 1];
		ft0[i__ - 1] = arprc_1.ftar[i__ - 1];
		px0[i__ - 1] = arprc_1.pxar[i__ - 1];
		py0[i__ - 1] = arprc_1.pyar[i__ - 1];
		pz0[i__ - 1] = arprc_1.pzar[i__ - 1];
		ee0[i__ - 1] = arprc_1.pear[i__ - 1];
		xm0[i__ - 1] = arprc_1.xmar[i__ - 1];
	}
	arindx_(&c_b96, &np, ft0, indx);
	i__1 = np;
	for (i__ = 1; i__ <= i__1; ++i__)
	{

		++npar;

		arprc_1.itypar[npar - 1] = ityp0[indx[i__ - 1] - 1];
		arprc_1.gxar[npar - 1] = gx0[indx[i__ - 1] - 1];
		arprc_1.gyar[npar - 1] = gy0[indx[i__ - 1] - 1];
		arprc_1.gzar[npar - 1] = gz0[indx[i__ - 1] - 1];
		arprc_1.ftar[npar - 1] = ft0[indx[i__ - 1] - 1];
		arprc_1.pxar[npar - 1] = px0[indx[i__ - 1] - 1];
		arprc_1.pyar[npar - 1] = py0[indx[i__ - 1] - 1];
		arprc_1.pzar[npar - 1] = pz0[indx[i__ - 1] - 1];
		arprc_1.pear[npar - 1] = ee0[indx[i__ - 1] - 1];
		arprc_1.xmar[npar - 1] = xm0[indx[i__ - 1] - 1];
	}
	arprnt_1.iaint2[0] = npar;

	return 0;
}

int arini2_(int *k)
{

	int i__1;

	static int i__, ip, irun;

	arerc1_1.multi1[*k - 1] = arprnt_1.iaint2[0];
	i__1 = arerc1_1.multi1[*k - 1];
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		arprc1_1.ityp1[i__ + *k * 150001 - 150002] = arprc_1.itypar[i__ - 1];
		arprc1_1.gx1[i__ + *k * 150001 - 150002] = arprc_1.gxar[i__ - 1];
		arprc1_1.gy1[i__ + *k * 150001 - 150002] = arprc_1.gyar[i__ - 1];
		arprc1_1.gz1[i__ + *k * 150001 - 150002] = arprc_1.gzar[i__ - 1];
		arprc1_1.ft1[i__ + *k * 150001 - 150002] = arprc_1.ftar[i__ - 1];
		arprc1_1.px1[i__ + *k * 150001 - 150002] = arprc_1.pxar[i__ - 1];
		arprc1_1.py1[i__ + *k * 150001 - 150002] = arprc_1.pyar[i__ - 1];
		arprc1_1.pz1[i__ + *k * 150001 - 150002] = arprc_1.pzar[i__ - 1];
		arprc1_1.ee1[i__ + *k * 150001 - 150002] = arprc_1.pear[i__ - 1];
		arprc1_1.xm1[i__ + *k * 150001 - 150002] = arprc_1.xmar[i__ - 1];

		dpert_1.dpp1[i__ + *k * 150001 - 150002] = 1.f;
	}

	for (ip = 1; ip <= 150001; ++ip)
	{
		tdecay_1.tfdcy[ip - 1] = input2_1.ntmax * input1_1.dt;
		tdecay_1.tft[ip - 1] = input2_1.ntmax * input1_1.dt;
	}

	for (irun = 1; irun <= 1; ++irun)
	{
		for (ip = 1; ip <= 150001; ++ip)
		{
			tdecay_1.tfdpi[ip + irun * 150001 - 150002] = input2_1.ntmax *input1_1.dt;
		}
	}
	return 0;
}

int iarflv_(int *ipdg)
{

	int ret_val;

	static float r__;
	double ranart_(int *);

	if (*ipdg == -1114)
	{
		ret_val = -6;
		return ret_val;
	}

	if (*ipdg == -2114)
	{
		ret_val = -7;
		return ret_val;
	}

	if (*ipdg == -2214)
	{
		ret_val = -8;
		return ret_val;
	}

	if (*ipdg == -2224)
	{
		ret_val = -9;
		return ret_val;
	}

	if (*ipdg == -2212)
	{
		ret_val = -1;
		return ret_val;
	}

	if (*ipdg == -2112)
	{
		ret_val = -2;
		return ret_val;
	}

	if (*ipdg == 221)
	{
		ret_val = 0;
		return ret_val;
	}

	if (*ipdg == 2212)
	{
		ret_val = 1;
		return ret_val;
	}

	if (*ipdg == 2112)
	{
		ret_val = 2;
		return ret_val;
	}

	if (*ipdg == -211)
	{
		ret_val = 3;
		return ret_val;
	}

	if (*ipdg == 111)
	{
		ret_val = 4;
		return ret_val;
	}

	if (*ipdg == 211)
	{
		ret_val = 5;
		return ret_val;
	}

	if (*ipdg == 1114)
	{
		ret_val = 6;
		return ret_val;
	}

	if (*ipdg == 2114)
	{
		ret_val = 7;
		return ret_val;
	}

	if (*ipdg == 2214)
	{
		ret_val = 8;
		return ret_val;
	}

	if (*ipdg == 2224)
	{
		ret_val = 9;
		return ret_val;
	}

	if (*ipdg == 3122)
	{
		ret_val = 14;
		return ret_val;
	}

	if (*ipdg == -3122)
	{
		ret_val = -14;
		return ret_val;
	}

	if (*ipdg == 3112)
	{
		ret_val = 15;
		return ret_val;
	}

	if (*ipdg == -3112)
	{
		ret_val = -15;
		return ret_val;
	}

	if (*ipdg == 3212)
	{
		ret_val = 16;
		return ret_val;
	}

	if (*ipdg == -3212)
	{
		ret_val = -16;
		return ret_val;
	}

	if (*ipdg == 3222)
	{
		ret_val = 17;
		return ret_val;
	}

	if (*ipdg == -3222)
	{
		ret_val = -17;
		return ret_val;
	}

	if (*ipdg == -321)
	{
		ret_val = 21;
		return ret_val;
	}

	if (*ipdg == 321)
	{
		ret_val = 23;
		return ret_val;
	}

	if (*ipdg == 311)
	{
		ret_val = 23;
		return ret_val;
	}

	if (*ipdg == -311)
	{
		ret_val = 21;
		return ret_val;
	}

	if (*ipdg == 310 || *ipdg == 130)
	{
		r__ = ranart_(&rndf77_1.nseed);
		if (r__ > .5f)
		{
			ret_val = 23;
		}
		else
		{
			ret_val = 21;
		}
		return ret_val;
	}

	if (*ipdg == -213)
	{
		ret_val = 25;
		return ret_val;
	}

	if (*ipdg == 113)
	{
		ret_val = 26;
		return ret_val;
	}

	if (*ipdg == 213)
	{
		ret_val = 27;
		return ret_val;
	}

	if (*ipdg == 223)
	{
		ret_val = 28;
		return ret_val;
	}

	if (*ipdg == 333)
	{
		ret_val = 29;
		return ret_val;
	}

	if (*ipdg == 323)
	{
		ret_val = 30;
		return ret_val;
	}

	if (*ipdg == -323)
	{
		ret_val = -30;
		return ret_val;
	}

	if (*ipdg == 313)
	{
		ret_val = 30;
		return ret_val;
	}

	if (*ipdg == -313)
	{
		ret_val = -30;
		return ret_val;
	}

	if (*ipdg == 331)
	{
		ret_val = 31;
		return ret_val;
	}

	if (*ipdg == 3312)
	{
		ret_val = 40;
		return ret_val;
	}

	if (*ipdg == -3312)
	{
		ret_val = -40;
		return ret_val;
	}

	if (*ipdg == 3322)
	{
		ret_val = 41;
		return ret_val;
	}

	if (*ipdg == -3322)
	{
		ret_val = -41;
		return ret_val;
	}

	if (*ipdg == 3334)
	{
		ret_val = 45;
		return ret_val;
	}

	if (*ipdg == -3334)
	{
		ret_val = -45;
		return ret_val;
	}

	if (*ipdg == 6666)
	{
		ret_val = 44;
		return ret_val;
	}

	ret_val = *ipdg + 10000;
	return ret_val;
}

int invflv_(int *iart)
{

	int ret_val;

	double ranart_(int *);

	if (*iart == -6)
	{
		ret_val = -1114;
		return ret_val;
	}

	if (*iart == -7)
	{
		ret_val = -2114;
		return ret_val;
	}

	if (*iart == -8)
	{
		ret_val = -2214;
		return ret_val;
	}

	if (*iart == -9)
	{
		ret_val = -2224;
		return ret_val;
	}

	if (*iart == -1)
	{
		ret_val = -2212;
		return ret_val;
	}

	if (*iart == -2)
	{
		ret_val = -2112;
		return ret_val;
	}

	if (*iart == 0)
	{
		ret_val = 221;
		return ret_val;
	}

	if (*iart == 1)
	{
		ret_val = 2212;
		return ret_val;
	}

	if (*iart == 2)
	{
		ret_val = 2112;
		return ret_val;
	}

	if (*iart == 3)
	{
		ret_val = -211;
		return ret_val;
	}

	if (*iart == 4)
	{
		ret_val = 111;
		return ret_val;
	}

	if (*iart == 5)
	{
		ret_val = 211;
		return ret_val;
	}

	if (*iart == 6)
	{
		ret_val = 1114;
		return ret_val;
	}

	if (*iart == 7)
	{
		ret_val = 2114;
		return ret_val;
	}

	if (*iart == 8)
	{
		ret_val = 2214;
		return ret_val;
	}

	if (*iart == 9)
	{
		ret_val = 2224;
		return ret_val;
	}

	if (*iart == 14)
	{
		ret_val = 3122;
		return ret_val;
	}

	if (*iart == -14)
	{
		ret_val = -3122;
		return ret_val;
	}

	if (*iart == 15)
	{
		ret_val = 3112;
		return ret_val;
	}

	if (*iart == -15)
	{
		ret_val = -3112;
		return ret_val;
	}

	if (*iart == 16)
	{
		ret_val = 3212;
		return ret_val;
	}

	if (*iart == -16)
	{
		ret_val = -3212;
		return ret_val;
	}

	if (*iart == 17)
	{
		ret_val = 3222;
		return ret_val;
	}

	if (*iart == -17)
	{
		ret_val = -3222;
		return ret_val;
	}

	if (*iart == 21)
	{

		ret_val = -321;

		return ret_val;
	}

	if (*iart == 23)
	{

		ret_val = 321;

		return ret_val;
	}

	if (*iart == 22)
	{
		ret_val = 130;
		return ret_val;
	}

	if (*iart == 24)
	{
		ret_val = 310;
		return ret_val;
	}

	if (*iart == 25)
	{
		ret_val = -213;
		return ret_val;
	}

	if (*iart == 26)
	{
		ret_val = 113;
		return ret_val;
	}

	if (*iart == 27)
	{
		ret_val = 213;
		return ret_val;
	}

	if (*iart == 28)
	{
		ret_val = 223;
		return ret_val;
	}

	if (*iart == 29)
	{
		ret_val = 333;
		return ret_val;
	}

	if (*iart == 30)
	{
		ret_val = 323;
		if (ranart_(&rndf77_1.nseed) > .5f)
		{
			ret_val = 313;
		}
		return ret_val;
	}

	if (*iart == -30)
	{
		ret_val = -323;
		if (ranart_(&rndf77_1.nseed) > .5f)
		{
			ret_val = -313;
		}
		return ret_val;
	}

	if (*iart == 31)
	{
		ret_val = 331;
		return ret_val;
	}

	if (*iart == 32)
	{
		ret_val = 777;
		return ret_val;
	}

	if (*iart == 40)
	{
		ret_val = 3312;
		return ret_val;
	}

	if (*iart == -40)
	{
		ret_val = -3312;
		return ret_val;
	}

	if (*iart == 41)
	{
		ret_val = 3322;
		return ret_val;
	}

	if (*iart == -41)
	{
		ret_val = -3322;
		return ret_val;
	}

	if (*iart == 45)
	{
		ret_val = 3334;
		return ret_val;
	}

	if (*iart == -45)
	{
		ret_val = -3334;
		return ret_val;
	}

	if (*iart == 44)
	{
		ret_val = 6666;
		return ret_val;
	}

	if (*iart == 42)
	{
		ret_val = 42;
		return ret_val;
	}
	else if (*iart == -42)
	{
		ret_val = -42;
		return ret_val;
	}

	ret_val = *iart - 10000;
	return ret_val;
}

int ardata_(void)
{
	return 0;
}

int arindx_(int *n, int *m, float *arrin, int *indx)
{

	int i__1;
	static int i__, j, l;
	static float q;
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

int newka_(int *icase, int *irun, int *iseed,float *dt, int *nt, int *ictrl, int *i1, int *i2, float *srt, float *pcx, float *pcy, float *pcz, int *iblock)
{

	float r__1, r__2, r__3, r__4, r__5, r__6, r__7, r__8;

	double sqrt(double);

	static float ec;
	static int ic, ik, il, im;
	static float ds;
	static int in, lb1, lb2;
	static float em1, em2;
	static int ik0, ik1, ik2, ik3, im3, im4;
	static float sig, dsr, brel;
	extern int npik_(int *, int *, float *, int *,int *, int *, int *, float *, float *, float *, float *,int *, float *, int *);
	static logical lb1bn, lb2bn, lb1mn, lb2mn;
	static float px1cm, py1cm, pz1cm;
	static logical lb1bn0, lb1bn1, lb2bn1, lb2bn0, lb1mn0, lb2mn0, lb1mn1,lb2mn1, lb1mn2, lb2mn2;
	static int nchrg;
	extern double akpel_(float *);
	static float brsig, bmass, pkaon;
	extern double aknel_(float *);
	static float brsgm;
	static int ipipi;
	extern int kaonn_(float *, float *, int *, int *,float *, int *, int *, int *, int *, int *,float *, float *, float *, float *, int *);
	static float sgsum, sigma0;
	extern double pinsg0_(float *), x2kaon_(float *);
	static float sgsum1, sgsum3, sigela;
	extern double akplam_(float *), aknlam_(float *);
	extern int distce_(int *, int *, float *, float *,float *, float *, float *, int *, float *, float *, float *);
	extern double aknsgm_(float *), akpsgm_(float *);
	static int ielstc;
	extern int nnkaon_(int *, int *, int *,int *, int *, int *, float *, float *, float *, float *,int *);
	double ranart_(int *);
	static float ratiok, sigsgm;
	static int inpion;
	extern int pihypn_(int *, int *, int *, float *, int *, int *, int *, int *, float *, float *,float *, float *, int *, int *);

	*icase = -1;

	nchrg = -100;

	*ictrl = 1;
	lb1 = ee_1.lb[*i1 - 1];
	lb2 = ee_1.lb[*i2 - 1];
	em1 = cc_1.e[*i1 - 1];
	em2 = cc_1.e[*i2 - 1];
	lb1bn = lb1 == 1 || lb1 == 2 || lb1 > 5 && lb1 <= 13;
	lb2bn = lb2 == 1 || lb2 == 2 || lb2 > 5 && lb2 <= 13;
	lb1bn0 = lb1 == 2 || lb1 == 7 || lb1 == 10 || lb1 == 12;
	lb2bn0 = lb2 == 2 || lb2 == 7 || lb2 == 10 || lb2 == 12;
	lb1bn1 = lb1 == 1 || lb1 == 8 || lb1 == 11 || lb1 == 13;
	lb2bn1 = lb2 == 1 || lb2 == 8 || lb2 == 11 || lb2 == 13;
	lb1mn = em1 < .2f || lb1 == 0 || lb1 >= 25 && lb1 <= 29;
	lb2mn = em2 < .2f || lb2 == 0 || lb2 >= 25 && lb2 <= 29;
	lb1mn0 = lb1 == 0 || lb1 == 4 || lb1 == 26 || lb1 == 28 || lb1 == 29;
	lb2mn0 = lb2 == 0 || lb2 == 4 || lb2 == 26 || lb2 == 28 || lb2 == 29;
	lb1mn1 = lb1 == 5 || lb1 == 27;
	lb2mn1 = lb2 == 5 || lb2 == 27;
	lb1mn2 = lb1 == 3 || lb1 == 25;
	lb2mn2 = lb2 == 3 || lb2 == 25;

	if (lb1bn && lb2bn)
	{

		*icase = 1;

		sig = 40.f;
		if (lb1 == 9 && lb2 == 9)
		{
			nchrg = 4;
		}
		if (lb1bn1 && lb2 == 9 || lb2bn1 && lb1 == 9)
		{
			nchrg = 3;
		}
		if (lb1bn0 && lb2 == 9 || lb2bn0 && lb1 == 9 || lb1bn1 && lb2bn1)
		{
			nchrg = 2;
		}
		if (lb1bn1 && lb2bn0 || lb1 == 6 && lb2 == 9 || lb2bn1 && lb1bn0 ||
			lb2 == 6 && lb1 == 9)
		{
			nchrg = 1;
		}
		if (lb1bn0 && lb2bn0 || lb1bn1 && lb2 == 6 || lb2bn1 && lb1 == 6)
		{
			nchrg = 0;
		}
		if (lb1bn0 && lb2 == 6 || lb2bn0 && lb1 == 6)
		{
			nchrg = -1;
		}
		if (lb1 == 6 && lb2 == 6)
		{
			nchrg = -2;
		}

		if (nchrg >= -1 && nchrg <= 2)
		{

			brsig = x2kaon_(srt);
		}
		else
		{
			brsig = 0.f;
		}

		brsig *= 2.f;
	}

	if (lb1bn && lb2mn || lb2bn && lb1mn)
	{

		*icase = 2;
		sig = 20.f;
		sigma0 = pinsg0_(srt);
		brsig = 0.f;
		if (lb1bn1 && lb2mn0 || lb2bn1 && lb1mn0 || lb1bn0 && lb2mn1 ||
			lb2bn0 && lb1mn1 || lb1 == 9 && lb2mn2 || lb2 == 9 && lb1mn2)
		{
			nchrg = 1;

			if (lb1bn1 || lb2bn1)
			{
				brsig = sigma0 * .5f;
			}
			if (lb1bn0 || lb2bn0)
			{
				brsig = sigma0 * 2.f;
			}
		}
		if (lb1bn0 && lb2mn0 || lb2bn0 && lb1mn0 || lb1bn1 && lb2mn2 ||
			lb2bn1 && lb1mn2 || lb1 == 6 && lb2mn1 || lb2 == 6 && lb1mn1)
		{
			nchrg = 0;
			if (lb1bn1 || lb2bn1)
			{

				brsig = sigma0 * 3.f;

				ratiok = .66666666666666663f;
			}
			if (lb1bn0 || lb2bn0)
			{
				brsig = sigma0 * 2.5f;

				ratiok = .2f;
			}
		}
		if (lb1bn0 && lb2mn2 || lb2bn0 && lb1mn2 || lb1 == 6 && lb2mn0 || lb2 == 6 && lb1mn0)
		{
			nchrg = -1;
			if (lb1bn0 || lb2bn0)
			{
				brsig = sigma0;
			}
		}

		if (lb1 == 6 && lb2mn2 || lb2 == 6 && lb1mn2)
		{
			nchrg = -2;
		}

		if (lb1bn1 && lb2mn1 || lb2bn1 && lb1mn1 || lb1 == 9 && lb2mn0 || lb2 == 9 && lb1mn0)
		{
			nchrg = 2;
		}

		if (nchrg >= -2 && nchrg <= 2)
		{
			brsig = sigma0 * 3.f;
		}
	}

	if (lb1bn && (lb2 == 21 || lb2 == -30) || lb2bn && (lb1 == 21 || lb1 ==
																		 -30))
	{

		bmass = .938f;
		if (*srt <= bmass + .498f)
		{

			pkaon = 0.f;
		}
		else
		{

			r__2 = *srt;

			r__3 = bmass;

			r__1 = (r__2 * r__2 - (r__3 * r__3 + .248004f)) / 2.f / bmass;
			pkaon = sqrt(r__1 * r__1 - .248004f);
		}
		sig = 0.f;
		if (lb1 == 1 || lb2 == 1 || lb1 == 8 || lb2 == 8 || lb1 == 11 || lb2 == 11 || lb1 == 13 || lb2 == 13)
		{

			nchrg = 0;
			sigela = akpel_(&pkaon);
			sigsgm = akpsgm_(&pkaon) * 3.f;
			sig = sigela + sigsgm + akplam_(&pkaon);
		}
		if (lb1 == 2 || lb2 == 2 || lb1 == 7 || lb2 == 7 || lb1 == 10 || lb2 == 10 || lb1 == 12 || lb2 == 12)
		{

			nchrg = -1;
			sigela = aknel_(&pkaon);
			sigsgm = aknsgm_(&pkaon) * 2.f;
			sig = sigela + sigsgm + aknlam_(&pkaon);
		}
		if (lb1 == 6 || lb2 == 6)
		{

			nchrg = -2;
			sigela = aknel_(&pkaon);
			sigsgm = aknsgm_(&pkaon);
			sig = sigela + sigsgm;
		}
		if (lb1 == 9 || lb2 == 9)
		{

			nchrg = 1;
			sigela = akpel_(&pkaon);
			sigsgm = akpsgm_(&pkaon) * 2.f;
			sig = sigela + sigsgm + akplam_(&pkaon);
		}

		sigela = (akpel_(&pkaon) + aknel_(&pkaon)) * .5f;
		sigsgm = akpsgm_(&pkaon) * 1.5f + aknsgm_(&pkaon);
		sig = sigela + sigsgm + akplam_(&pkaon);

		if (sig > 1e-7f)
		{

			*icase = 3;
			brel = sigela / sig;
			brsgm = sigsgm / sig;

			brsig = sig;
		}
	}

	if (lb1 >= 14 && lb1 <= 17 && (lb2 >= 3 && lb2 <= 5) || lb2 >= 14 && lb2 <= 17 && (lb1 >= 3 && lb1 <= 5))
	{

		nchrg = -100;
		if (lb1 == 15 && (lb2 == 3 || lb2 == 25) || lb2 == 15 && (lb1 == 3 ||
																  lb1 == 25))
		{
			nchrg = -2;

			bmass = 1.232f;
		}
		if (lb1 == 15 && lb2mn0 || lb2 == 15 && lb1mn0 || (lb1 == 14 || lb1 == 16) && (lb2 == 3 || lb2 == 25) || (lb2 == 14 || lb2 == 16) && (lb1 == 3 || lb1 == 25))
		{
			nchrg = -1;

			bmass = .938f;
		}
		if (lb1 == 15 && (lb2 == 5 || lb2 == 27) || lb2 == 15 && (lb1 == 5 || lb1 == 27) || lb1 == 17 && (lb2 == 3 || lb2 == 25) || lb2 == 17 && (lb1 == 3 || lb1 == 25) || (lb1 == 14 || lb1 == 16) && lb2mn0 || (lb2 == 14 || lb2 == 16) && lb1mn0)
		{
			nchrg = 0;

			bmass = .938f;
		}
		if (lb1 == 17 && lb2mn0 || lb2 == 17 && lb1mn0 || (lb1 == 14 || lb1 == 16) && (lb2 == 5 || lb2 == 27) || (lb2 == 14 || lb2 == 16) && (lb1 == 5 || lb1 == 27))
		{
			nchrg = 1;

			bmass = 1.232f;
		}
		sig = 0.f;
		if (nchrg != -100 && *srt > bmass + .498f)
		{

			*icase = 4;

			r__2 = *srt;

			r__1 = (r__2 * r__2 - 1.1278479999999997f) / 2.f / .938f;
			pkaon = sqrt(r__1 * r__1 - .248004f);

			if (lb1 == 14 || lb2 == 14)
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

			r__1 = *srt;

			r__2 = bmass + .498f;

			r__3 = *srt;

			r__4 = .498f - bmass;

			r__5 = *srt;

			r__6 = em1 + em2;

			r__7 = *srt;

			r__8 = em1 - em2;
			sig = (r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4) /
				  (r__5 * r__5 - r__6 * r__6) / (r__7 * r__7 - r__8 * r__8) * sigma0;

			if (nchrg == -2 || nchrg == 2)
			{
				sig *= 2.f;
			}

			if (lb1 == 14 || lb2 == 14)
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

			brsig = sig;
			if (sig < 1e-7f)
			{
				sig = 1e-7f;
			}
		}

		*icase = 4;
		brsig = sig;

		sigela = 10.f;
		sig += sigela;
		brel = sigela / sig;
	}

	if (*icase == -1)
	{
		*ictrl = -1;
		return 0;
	}
	px1cm = *pcx;
	py1cm = *pcy;
	pz1cm = *pcz;
	ds = sqrt(sig / 31.4f);
	dsr = ds + .1f;

	r__1 = em1 + em2 + .02f;
	ec = r__1 * r__1;

	distce_(i1, i2, &dsr, &ds, dt, &ec, srt, &ic, &px1cm, &py1cm, &pz1cm);
	if (ic == -1)
	{

		*ictrl = -1;

		return 0;
	}

	ik = 0;
	ik0 = 0;
	ik1 = 0;
	ik2 = 0;
	ik3 = 0;
	il = 0;
	im = 0;
	im3 = 0;
	im4 = 0;
	in = 0;
	inpion = 0;
	ipipi = 0;
	sgsum = 0.f;
	sgsum1 = 0.f;
	sgsum3 = 0.f;
	if (*icase == 1)
	{
		++ik;
		if (*srt > 2.8639f)
		{
			++ik0;
			if (em1 < 1.f && em2 < 1.f)
			{
				++ik1;
				sgsum1 += brsig;
			}
			if (em1 > 1.f && em2 > 1.f)
			{
				++ik3;
				sgsum3 += brsig;
			}
			if (em1 > 1.f && em2 < 1.f)
			{
				++ik2;
			}
			if (em1 < 1.f && em2 > 1.f)
			{
				++ik2;
			}
			sgsum += brsig;
		}
	}
	if (*icase == 2)
	{
		++inpion;
	}
	if (*icase == 5)
	{
		++ipipi;
	}

	if (ranart_(&rndf77_1.nseed) > brsig / sig)
	{

		*ictrl = -1;
		return 0;
	}
	++il;

	if (*icase == 1)
	{
		++in;

		nnkaon_(irun, iseed, ictrl, i1, i2, iblock, srt, pcx, pcy, pcz, &nchrg);
	}
	if (*icase == 2)
	{
		++im;

		npik_(irun, iseed, dt, nt, ictrl, i1, i2, srt, pcx, pcy, pcz, &nchrg,
			  &ratiok, iblock);
	}

	if (*icase == 3)
	{
		++im3;

		kaonn_(&brel, &brsgm, irun, iseed, dt, nt, ictrl, i1, i2, iblock, srt,
			   pcx, pcy, pcz, &nchrg);
	}

	if (*icase == 4)
	{
		++im4;

		if (ranart_(&rndf77_1.nseed) < brel)
		{
			ielstc = 1;
		}
		else
		{
			ielstc = 0;
		}

		pihypn_(&ielstc, irun, iseed, dt, nt, ictrl, i1, i2, srt, pcx, pcy,
				pcz, &nchrg, iblock);
	}

	return 0;
}

double x2kaon_(float *srt)
{

	float ret_val, r__1, r__2;
	double d__1;

	double sqrt(double), log(double), pow_dd(double *, double *);

	static float x, f1, f2, f3, smin, sigma1, sigma2, sigma3;

	smin = 2.8639f;
	ret_val = 1e-7f;
	if (*srt < smin)
	{
		return ret_val;
	}
	sigma1 = 2.8f;
	sigma2 = 7.7f;
	sigma3 = 3.9f;

	r__1 = *srt;

	r__2 = smin;
	x = r__1 * r__1 / (r__2 * r__2) + 1e-7f;
	f1 = (1.f / sqrt(x) + 1.f) * log(x) - (1.f - 1.f / sqrt(x)) * 4.f;
	f2 = 1.f - 1.f / sqrt(x) * (log(sqrt(x)) + 1.f);

	r__1 = x;
	d__1 = (double)((x - 1.f) / (r__1 * r__1));
	f3 = pow_dd(&d__1, &c_b112);

	r__1 = 1.f - 1.f / x;
	ret_val = r__1 * (r__1 * r__1) * (sigma1 * f1 + sigma2 * f2) + sigma3 *
																	   f3;
	return ret_val;
}

double pinsg0_(float *srt)
{

	float ret_val, r__1, r__2;
	double d__1;

	double pow_dd(double *, double *);

	static float srt0, ratio;

	srt0 = 1.9339999999999999f;
	if (*srt < srt0)
	{
		ret_val = 0.f;
		return ret_val;
	}

	r__1 = srt0;

	r__2 = *srt;
	ratio = r__1 * r__1 / (r__2 * r__2);
	d__1 = (double)(1.f - ratio);

	r__1 = ratio;
	ret_val = pow_dd(&d__1, &c_b114) * 1.121f * (r__1 * r__1);
	return ret_val;
}

double aknel_(float *pkaon)
{

	float ret_val;
	double d__1;

	double pow_dd(double *, double *);

	static float sigma1;

	if (*pkaon < .5f || *pkaon >= 4.f)
	{
		sigma1 = 0.f;
	}
	if (*pkaon >= .5f && *pkaon < 1.f)
	{
		d__1 = (double)(*pkaon);
		sigma1 = pow_dd(&d__1, &c_b116) * 20.f;
	}
	if (*pkaon >= 1.f && *pkaon < 4.f)
	{
		d__1 = (double)(*pkaon);
		sigma1 = pow_dd(&d__1, &c_b117) * 20.f;
	}
	ret_val = sigma1;
	return ret_val;
}

double akpel_(float *pkaon)
{

	float ret_val;
	double d__1;

	double pow_dd(double *, double *);

	static float sigma2;

	if (*pkaon < .25f || *pkaon >= 4.f)
	{
		sigma2 = 0.f;
	}
	if (*pkaon >= .25f && *pkaon < 4.f)
	{
		d__1 = (double)(*pkaon);
		sigma2 = pow_dd(&d__1, &c_b119) * 13.f;
	}
	ret_val = sigma2;
	return ret_val;
}

double aknsgm_(float *pkaon)
{

	float ret_val;
	double d__1;

	double pow_dd(double *, double *);

	static float sigma2;

	if (*pkaon < .5f || *pkaon >= 6.f)
	{
		sigma2 = 0.f;
	}
	if (*pkaon >= .5f && *pkaon < 1.f)
	{
		d__1 = (double)(*pkaon);
		sigma2 = pow_dd(&d__1, &c_b121) * 1.2f;
	}
	if (*pkaon >= 1.f && *pkaon < 6.f)
	{
		d__1 = (double)(*pkaon);
		sigma2 = pow_dd(&d__1, &c_b122) * 1.2f;
	}
	ret_val = sigma2;
	return ret_val;
}

double akpsgm_(float *pkaon)
{

	float ret_val;
	double d__1;

	double pow_dd(double *, double *);

	static float sigma1;

	if (*pkaon < .2f || *pkaon >= 1.5f)
	{
		sigma1 = 0.f;
	}
	if (*pkaon >= .2f && *pkaon < 1.5f)
	{
		d__1 = (double)(*pkaon);
		sigma1 = pow_dd(&d__1, &c_b117) * .6f;
	}
	ret_val = sigma1;
	return ret_val;
}

double akplam_(float *pkaon)
{

	float ret_val, r__1;
	double d__1;

	double pow_dd(double *, double *);

	static float p, sigma;

	p = *pkaon;
	if (*pkaon < .2f || *pkaon >= 10.f)
	{
		sigma = 0.f;
	}
	if (*pkaon >= .2f && *pkaon < .9f)
	{

		r__1 = p;
		sigma = r__1 * r__1 * 50.f - p * 67.f + 24.f;
	}
	if (*pkaon >= .9f && *pkaon < 10.f)
	{
		d__1 = (double)(*pkaon);
		sigma = pow_dd(&d__1, &c_b126) * 3.f;
	}
	ret_val = sigma;
	return ret_val;
}

double aknlam_(float *pkaon)
{

	float ret_val;

	extern double akplam_(float *);

	ret_val = akplam_(pkaon);
	return ret_val;
}

double aknpsg_(float *pkaon)
{

	float ret_val;
	double d__1;

	double pow_dd(double *, double *);

	static float sigma1;

	if (*pkaon <= .345f)
	{
		d__1 = (double)(*pkaon);
		sigma1 = pow_dd(&d__1, &c_b129) * .624f;
	}
	else
	{
		d__1 = (double)(*pkaon);
		sigma1 = pow_dd(&d__1, &c_b130) * .7f;
	}
	ret_val = sigma1;
	return ret_val;
}

int nnkaon_(int *irun, int *iseed, int *ictrl,
			int *i1, int *i2, int *iblock, float *srt, float *pcx, float *pcy, float *pcz, int *nchrg)
{

	float r__1, r__2, r__3, r__4;

	double sqrt(double);

	static int n;
	static float px[4], py[4], pz[4];
	static int lb1, lb2;
	static float dm3, dm4, e1cm, e2cm, eti1, eti2, pt1i1, pt2i1, pt3i1, pt1i2,
		pt2i2, pt3i2;
	static int iflag;
	static float betak, epcmk, p1beta, p2beta, betaak, epcmak;
	extern int fstate_(int *, float *, float *, float *,
					   float *, float *, float *, int *);
	static float transf;
	extern int rotate_(float *, float *, float *, float *, float *, float *);
	static float pxrota, pyrota, pzrota;

	dm3 = .938f;
	dm4 = .938f;

	n = 0;

	if (*nchrg <= -1 || *nchrg >= 3)
	{
		dm3 = 1.232f;
	}
	if (*nchrg == -2 || *nchrg == 4)
	{
		dm4 = 1.232f;
	}

	*iblock = 0;
	fstate_(iseed, srt, &dm3, &dm4, px, py, pz, &iflag);
	if (iflag < 0)
	{

		*ictrl = -1;
		++n;
		return 0;
	}
	*iblock = 12;

	pxrota = px[0];
	pyrota = py[0];
	pzrota = pz[0];

	rotate_(pcx, pcy, pcz, &pxrota, &pyrota, &pzrota);
	px[0] = pxrota;
	py[0] = pyrota;
	pz[0] = pzrota;

	pxrota = px[1];
	pyrota = py[1];
	pzrota = pz[1];

	rotate_(pcx, pcy, pcz, &pxrota, &pyrota, &pzrota);
	px[1] = pxrota;
	py[1] = pyrota;
	pz[1] = pzrota;

	pxrota = px[2];
	pyrota = py[2];
	pzrota = pz[2];

	rotate_(pcx, pcy, pcz, &pxrota, &pyrota, &pzrota);
	px[2] = pxrota;
	py[2] = pyrota;
	pz[2] = pzrota;

	pxrota = px[3];
	pyrota = py[3];
	pzrota = pz[3];

	rotate_(pcx, pcy, pcz, &pxrota, &pyrota, &pzrota);
	px[3] = pxrota;
	py[3] = pyrota;
	pz[3] = pzrota;
	nn_1.nnn += 2;

	pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 23;
	if (*nchrg == -1 || *nchrg == -2)
	{
	}

	pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = .498f;

	pd_1.lpion[nn_1.nnn - 1 + *irun * 150001 - 150002] = 21;

	pc_1.epion[nn_1.nnn - 1 + *irun * 150001 - 150002] = .498f;

	r__1 = dm3;

	r__2 = px[0];

	r__3 = py[0];

	r__4 = pz[0];
	e1cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	p1beta = px[0] * bg_1.betax + py[0] * bg_1.betay + pz[0] * bg_1.betaz;
	transf = bg_1.gamma * (bg_1.gamma * p1beta / (bg_1.gamma + 1) + e1cm);
	pt1i1 = bg_1.betax * transf + px[0];
	pt2i1 = bg_1.betay * transf + py[0];
	pt3i1 = bg_1.betaz * transf + pz[0];
	eti1 = dm3;

	lb1 = 2;
	if (*nchrg >= -2 && *nchrg <= 1)
	{
		lb1 = 2;
	}

	if (*nchrg == -2 || *nchrg == -1)
	{
		lb1 = 6;
	}

	if (*nchrg == 1 || *nchrg == 2)
	{
		lb1 = 1;
	}
	if (*nchrg == 3 || *nchrg == 4)
	{
		lb1 = 9;
	}

	r__1 = dm4;

	r__2 = px[1];

	r__3 = py[1];

	r__4 = pz[1];
	e2cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	p2beta = px[1] * bg_1.betax + py[1] * bg_1.betay + pz[1] * bg_1.betaz;
	transf = bg_1.gamma * (bg_1.gamma * p2beta / (bg_1.gamma + 1) + e2cm);
	pt1i2 = bg_1.betax * transf + px[1];
	pt2i2 = bg_1.betay * transf + py[1];
	pt3i2 = bg_1.betaz * transf + pz[1];
	eti2 = dm4;

	lb2 = 2;

	if (*nchrg >= -1 || *nchrg <= 1)
	{
		lb2 = 2;
	}
	if (*nchrg == 2 || *nchrg == 3)
	{
		lb2 = 1;
	}
	if (*nchrg == 4)
	{
		lb2 = 9;
	}
	if (*nchrg == -2)
	{
		lb2 = 6;
	}

	bb_1.p[*i1 * 3 - 3] = pt1i1;
	bb_1.p[*i1 * 3 - 2] = pt2i1;
	bb_1.p[*i1 * 3 - 1] = pt3i1;
	cc_1.e[*i1 - 1] = eti1;
	ee_1.lb[*i1 - 1] = lb1;
	bb_1.p[*i2 * 3 - 3] = pt1i2;
	bb_1.p[*i2 * 3 - 2] = pt2i2;
	bb_1.p[*i2 * 3 - 1] = pt3i2;
	cc_1.e[*i2 - 1] = eti2;
	ee_1.lb[*i2 - 1] = lb2;

	r__1 = pc_1.epion[nn_1.nnn - 1 + *irun * 150001 - 150002];

	r__2 = px[2];

	r__3 = py[2];

	r__4 = pz[2];
	epcmk = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	betak = px[2] * bg_1.betax + py[2] * bg_1.betay + pz[2] * bg_1.betaz;
	transf = bg_1.gamma * (bg_1.gamma * betak / (bg_1.gamma + 1.f) + epcmk);
	pb_1.ppion[(nn_1.nnn - 1 + *irun * 150001) * 3 - 450006] = bg_1.betax *
																   transf +
															   px[2];
	pb_1.ppion[(nn_1.nnn - 1 + *irun * 150001) * 3 - 450005] = bg_1.betay *
																   transf +
															   py[2];
	pb_1.ppion[(nn_1.nnn - 1 + *irun * 150001) * 3 - 450004] = bg_1.betaz *
																   transf +
															   pz[2];
	pa_1.rpion[(nn_1.nnn - 1 + *irun * 150001) * 3 - 450006] = aa_1.r__[*i1 *
																			3 -
																		3];
	pa_1.rpion[(nn_1.nnn - 1 + *irun * 150001) * 3 - 450005] = aa_1.r__[*i1 *
																			3 -
																		2];
	pa_1.rpion[(nn_1.nnn - 1 + *irun * 150001) * 3 - 450004] = aa_1.r__[*i1 *
																			3 -
																		1];

	dpert_1.dppion[nn_1.nnn - 1 + *irun * 150001 - 150002] = dpert_1.dpertp[*i1 - 1] * dpert_1.dpertp[*i2 - 1];

	r__1 = pc_1.epion[nn_1.nnn + *irun * 150001 - 150002];

	r__2 = px[3];

	r__3 = py[3];

	r__4 = pz[3];
	epcmak = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	betaak = px[3] * bg_1.betax + py[3] * bg_1.betay + pz[3] * bg_1.betaz;
	transf = bg_1.gamma * (bg_1.gamma * betaak / (bg_1.gamma + 1.f) + epcmak);
	pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = bg_1.betax *
															   transf +
														   px[3];
	pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = bg_1.betay *
															   transf +
														   py[3];
	pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = bg_1.betaz *
															   transf +
														   pz[3];
	pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = aa_1.r__[*i2 * 3 -
																	3];
	pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = aa_1.r__[*i2 * 3 -
																	2];
	pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = aa_1.r__[*i2 * 3 -
																	1];

	dpert_1.dppion[nn_1.nnn + *irun * 150001 - 150002] = dpert_1.dpertp[*i1 -
																		1] *
														 dpert_1.dpertp[*i2 - 1];
	return 0;
}

int lorntz_(int *ilo, float *b, float *pi, float *pj)
{

	double sqrt(double);

	static int i__;
	static float bb, ga, gam, pib, pjb, deno3;

	--pj;
	--pi;
	--b;

	bb = b[1] * b[1] + b[2] * b[2] + b[3] * b[3];
	deno3 = sqrt(1.f - bb);
	if (deno3 == 0.f)
	{
		deno3 = 1e-10f;
	}
	gam = 1.f / deno3;
	ga = gam * gam / (gam + 1.f);
	if (*ilo == 1)
	{
		pib = pi[1] * b[1] + pi[2] * b[2] + pi[3] * b[3];
		pjb = pj[1] * b[1] + pj[2] * b[2] + pj[3] * b[3];
	for (i__ = 1; i__ <= 3; ++i__)
	{
		pi[i__] += b[i__] * (ga * pib + gam * pi[4]);
		pj[i__] += b[i__] * (ga * pjb + gam * pj[4]);
	}
	pi[4] = gam * (pi[4] + pib);
	pj[4] = gam * (pj[4] + pjb);
	return 0;
	}

	pib = pi[1] * b[1] + pi[2] * b[2] + pi[3] * b[3];
	pjb = pj[1] * b[1] + pj[2] * b[2] + pj[3] * b[3];

	for (i__ = 1; i__ <= 3; ++i__)
	{
		pi[i__] += b[i__] * (ga * pib - gam * pi[4]);
		pj[i__] += b[i__] * (ga * pjb - gam * pj[4]);
	}
	pi[4] = gam * (pi[4] - pib);
	pj[4] = gam * (pj[4] - pjb);
	return 0;	
}

int fstate_(int *iseed, float *srt, float *dm3, float *dm4,float *px, float *py, float *pz, int *iflag)
{

	float r__1, r__2, r__3, r__4;

	double sqrt(double);
	int s_wsle(cilist *), do_lio(int *, int *, char *, int),e_wsle(void);
	double log(double), cos(double), sin(double);
	static float v1, v2, pe[4], bbb, fac, aka, ekm, ekp, pio, ptp, rsq, ptp2,rmt3, rmt4, gama, beta, ptkm, ptkp, pzkm, pzkp, ekmax, pkmax,guass, pzcms, xstar, ptkmi2, ptkpl2;
	double ranart_(int *);
	static float resten;
	static int icount;
	static float restms, restpz;
	static cilist io___156 = {0, 1, 0, 0, 0};

	--pz;
	--py;
	--px;

	*iflag = -1;

	pio = 3.1415926f;
	aka = .498f;

	icount = 0;
	ekmax = (*srt - *dm3 - *dm4) / 2.f;
	if (ekmax <= aka)
	{
		return 0;
	}

	r__1 = ekmax;

	r__2 = aka;
	pkmax = sqrt(r__1 * r__1 - r__2 * r__2);
	if (*dm3 <= 0.f || *dm4 <= 0.f)
	{
		s_wsle(&io___156);
		do_lio(&c__9, &c__1, "error: minus mass!!!", (int)20);
		e_wsle();
		return 0;
	}

L50:
	++icount;
	if (icount > 10)
	{
		return 0;
	}
	ptkmi2 = log(ranart_(&rndf77_1.nseed)) * -.24125452352231608f;
	ptkm = sqrt(ptkmi2);
	v1 = ranart_(&rndf77_1.nseed);
	v2 = ranart_(&rndf77_1.nseed);
	r__1 = v1;
	r__2 = v2;
	rsq = r__1 * r__1 + r__2 * r__2;
	if (rsq >= 1.f || rsq <= 0.f)
	{
		v1 = ranart_(&rndf77_1.nseed);
		v2 = ranart_(&rndf77_1.nseed);
		r__1 = v1;
		r__2 = v2;
		rsq = r__1 * r__1 + r__2 * r__2;
	}
	fac = sqrt(log(rsq) * -2.f / rsq);
	guass = v1 * fac;
	if (guass >= 5.f)
	{
		v1 = ranart_(&rndf77_1.nseed);
		v2 = ranart_(&rndf77_1.nseed);
		r__1 = v1;
		r__2 = v2;
		rsq = r__1 * r__1 + r__2 * r__2;
	}
	xstar = guass / 5.f;
	pzkm = pkmax * xstar;

	r__1 = aka;

	r__2 = pzkm;

	r__3 = ptkm;
	ekm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
	if (ranart_(&rndf77_1.nseed) > aka / ekm){
		++icount;
	if (icount > 10)
	{
		return 0;
	}
	ptkmi2 = log(ranart_(&rndf77_1.nseed)) * -.24125452352231608f;
	ptkm = sqrt(ptkmi2);
	}
	bbb = ranart_(&rndf77_1.nseed);
	px[3] = ptkm * cos(pio * 2.f * bbb);
	py[3] = ptkm * sin(pio * 2.f * bbb);
	if (ranart_(&rndf77_1.nseed) > .5f)
	{
		pzkm *= -1.f;
	}
	pz[3] = pzkm;
	pe[2] = ekm;
	v1 = ranart_(&rndf77_1.nseed);
	v2 = ranart_(&rndf77_1.nseed);
	r__1 = v1;
	r__2 = v2;
	rsq = r__1 * r__1 + r__2 * r__2;
	if (rsq >= 1.f || rsq <= 0.f)
	{
        v1 = ranart_(&rndf77_1.nseed);
        v2 = ranart_(&rndf77_1.nseed);
        r__1 = v1;
        r__2 = v2;
        rsq = r__1 * r__1 + r__2 * r__2;
	}
	fac = sqrt(log(rsq) * -2.f / rsq);
	guass = v1 * fac;
	if (guass >= 3.25f)
	{
        v1 = ranart_(&rndf77_1.nseed);
        v2 = ranart_(&rndf77_1.nseed);
        r__1 = v1;
        r__2 = v2;
        rsq = r__1 * r__1 + r__2 * r__2;
	}
	xstar = guass / 3.25f;
	pzkp = pkmax * xstar;
	r__1 = aka;
	r__2 = pzkp;
	r__3 = ptkp;
	ekp = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
	if (ranart_(&rndf77_1.nseed) > aka / ekp)
	{
        ptkpl2 = log(ranart_(&rndf77_1.nseed)) * -.27173913043478259f;
        ptkp = sqrt(ptkpl2);

	}
	bbb = ranart_(&rndf77_1.nseed);
	px[4] = ptkp * cos(pio * 2.f * bbb);
	py[4] = ptkp * sin(pio * 2.f * bbb);
	if (ranart_(&rndf77_1.nseed) > .5f)
	{
		pzkp *= -1.f;
	}
	pz[4] = pzkp;
	pe[3] = ekp;
	resten = *srt - pe[2] - pe[3];
	restpz = -pz[3] - pz[4];

	if (resten <= dabs(restpz))
	{
		++icount;
	if (icount > 10)
	{
		return 0;
	}
	ptkmi2 = log(ranart_(&rndf77_1.nseed)) * -.24125452352231608f;
	ptkm = sqrt(ptkmi2);
	}

	r__1 = resten;

	r__2 = restpz;
	restms = sqrt(r__1 * r__1 - r__2 * r__2);

	if (restms < *dm3 + *dm4)
	{
		++icount;
	if (icount > 10)
	{
		return 0;
	}
	ptkmi2 = log(ranart_(&rndf77_1.nseed)) * -.24125452352231608f;
	ptkm = sqrt(ptkmi2);
	}
	ptp2 = log(ranart_(&rndf77_1.nseed)) * -.3623188405797102f;
	ptp = sqrt(ptp2);
	bbb = ranart_(&rndf77_1.nseed);
	px[2] = ptp * cos(pio * 2.f * bbb);
	py[2] = ptp * sin(pio * 2.f * bbb);
	px[1] = (px[4] + px[3] + px[2]) * -1.f;
	py[1] = (py[4] + py[3] + py[2]) * -1.f;
	r__1 = *dm3;
	r__2 = px[1];
	r__3 = py[1];
	rmt3 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
	r__1 = *dm4;
	r__2 = px[2];
	r__3 = py[2];
	rmt4 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
	if (restms < rmt3 + rmt4)
	{
		++icount;
	if (icount > 10)
	{
		return 0;
	}
	ptkmi2 = log(ranart_(&rndf77_1.nseed)) * -.24125452352231608f;
	ptkm = sqrt(ptkmi2);
	}

	r__1 = restms;
	r__2 = rmt3 + rmt4;
	r__3 = restms;
	r__4 = rmt3 - rmt4;
	pzcms = sqrt((r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4)) /2.f / restms;
	if (ranart_(&rndf77_1.nseed) > .5f)
	{
		pz[1] = pzcms;
		pz[2] = -pzcms;
	}
	else
	{
		pz[1] = -pzcms;
		pz[2] = pzcms;
	}
	beta = restpz / resten;
	r__1 = beta;
	gama = 1.f / sqrt(1.f - r__1 * r__1);
	r__1 = rmt3;
	r__2 = pz[1];
	pz[1] = pz[1] * gama + beta * gama * sqrt(r__1 * r__1 + r__2 * r__2);
	r__1 = rmt4;
	r__2 = pz[2];
	pz[2] = pz[2] * gama + beta * gama * sqrt(r__1 * r__1 + r__2 * r__2);
	r__1 = rmt3;
	r__2 = pz[1];
	pe[0] = sqrt(r__1 * r__1 + r__2 * r__2);
	r__1 = rmt4;
	r__2 = pz[2];
	pe[1] = sqrt(r__1 * r__1 + r__2 * r__2);
	*iflag = 1;
	return 0;
}

int npik_(int *irun, int *iseed, float *dt, int *nt, int *ictrl, int *i1, int *i2, float *srt, float *pcx,
		  float *pcy, float *pcz, int *nchrg, float *ratiok, int *iblock)
{

	float r__1, r__2, r__3, r__4;

	double sqrt(double), cos(double), sin(double);

	static int i__, k, k1, k2;
	static float p1[4], p2[4], p3[4], bb[3], pk;
	static int lb1, lb2;
	static float fai, eip;
	static int ilo;
	static float css, sss, e1cm, e2cm, eti1, eti2, pt1i1, pt2i1, pt3i1, pt1i2,
		pt2i2, pt3i2, rmnp, pznp, px1cm, py1cm, pz1cm, betak, epcmk,
		pkmax, p1beta, p2beta;
	double ranart_(int *);
	static float transf;
	extern int rotate_(float *, float *, float *, float *, float *, float *);
	static float pxrota, pyrota, pzrota;
	extern int lorntz_(int *, float *, float *, float *);

	px1cm = *pcx;
	py1cm = *pcy;
	pz1cm = *pcz;
	*ictrl = 1;
	lb1 = ee_1.lb[*i1 - 1];
	lb2 = ee_1.lb[*i2 - 1];
	k1 = *i1;
	k2 = *i2;

	if (lb2 == 1 || lb2 == 2 || lb2 >= 6 && lb2 <= 13)
	{
		k1 = *i2;
		k2 = *i1;
	}

	ee_1.lb[k1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 2) + 1;
	ee_1.lb[k2 - 1] = 23;

	r__1 = *srt;

	r__2 = *srt;
	pkmax = sqrt((r__1 * r__1 - 3.7403559999999998f) * (r__2 * r__2 -
														.87984399999999985f)) /
			2.f / *srt;
	pk = ranart_(&rndf77_1.nseed) * pkmax;

	css = 1.f - ranart_(&rndf77_1.nseed) * 2.f;

	r__1 = css;
	sss = sqrt(1.f - r__1 * r__1);
	fai = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
	p3[0] = pk * sss * cos(fai);
	p3[1] = pk * sss * sin(fai);
	p3[2] = pk * css;

	r__1 = pk;
	eip = *srt - sqrt(r__1 * r__1 + .248004f);

	r__1 = eip;

	r__2 = pk;
	rmnp = sqrt(r__1 * r__1 - r__2 * r__2);
	for (i__ = 1; i__ <= 3; ++i__)
	{
		bb[i__ - 1] = p3[i__ - 1] * -1.f / eip;
	}

	r__1 = rmnp;

	r__2 = rmnp;
	pznp = sqrt((r__1 * r__1 - 2.0620959999999999f) * (r__2 * r__2 -
													   .19359999999999997f)) /
		   2.f / rmnp;

	css = 1.f - ranart_(&rndf77_1.nseed) * 2.f;

	r__1 = css;
	sss = sqrt(1.f - r__1 * r__1);
	fai = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
	p1[0] = pznp * sss * cos(fai);
	p1[1] = pznp * sss * sin(fai);
	p1[2] = pznp * css;

	r__1 = pznp;
	p1[3] = sqrt(r__1 * r__1 + .87984399999999985f);

	r__1 = pznp;
	p2[3] = sqrt(r__1 * r__1 + .248004f);
	for (i__ = 1; i__ <= 3; ++i__)
	{
		p2[i__ - 1] = p1[i__ - 1] * -1.f;
	}

	ilo = 1;

	lorntz_(&ilo, bb, p1, p2);

	pxrota = p1[0];
	pyrota = p1[1];
	pzrota = p1[2];

	rotate_(pcx, pcy, pcz, &pxrota, &pyrota, &pzrota);
	p1[0] = pxrota;
	p1[1] = pyrota;
	p1[2] = pzrota;

	pxrota = p2[0];
	pyrota = p2[1];
	pzrota = p2[2];

	rotate_(pcx, pcy, pcz, &pxrota, &pyrota, &pzrota);
	p2[0] = pxrota;
	p2[1] = pyrota;
	p2[2] = pzrota;

	pxrota = p3[0];
	pyrota = p3[1];
	pzrota = p3[2];

	rotate_(pcx, pcy, pcz, &pxrota, &pyrota, &pzrota);
	p3[0] = pxrota;
	p3[1] = pyrota;
	p3[2] = pzrota;
	++nn_1.nnn;

	pd_1.lpion[nn_1.nnn + *irun * 150001 - 150002] = 21;

	pc_1.epion[nn_1.nnn + *irun * 150001 - 150002] = .498f;

	r__1 = p1[0];

	r__2 = p1[1];

	r__3 = p1[2];
	e1cm = sqrt(r__1 * r__1 + .87984399999999985f + r__2 * r__2 + r__3 * r__3);
	p1beta = p1[0] * bg_1.betax + p1[1] * bg_1.betay + p1[2] * bg_1.betaz;
	transf = bg_1.gamma * (bg_1.gamma * p1beta / (bg_1.gamma + 1) + e1cm);
	pt1i1 = bg_1.betax * transf + p1[0];
	pt2i1 = bg_1.betay * transf + p1[1];
	pt3i1 = bg_1.betaz * transf + p1[2];
	eti1 = .938f;
	lb1 = ee_1.lb[k1 - 1];

	r__1 = p3[0];

	r__2 = p3[1];

	r__3 = p3[2];
	e2cm = sqrt(r__1 * r__1 + .248004f + r__2 * r__2 + r__3 * r__3);
	p2beta = p3[0] * bg_1.betax + p3[1] * bg_1.betay + p3[2] * bg_1.betaz;
	transf = bg_1.gamma * (bg_1.gamma * p2beta / (bg_1.gamma + 1) + e2cm);
	pt1i2 = bg_1.betax * transf + p3[0];
	pt2i2 = bg_1.betay * transf + p3[1];
	pt3i2 = bg_1.betaz * transf + p3[2];
	eti2 = .498f;
	lb2 = ee_1.lb[k2 - 1];

	bb_1.p[k1 * 3 - 3] = pt1i1;
	bb_1.p[k1 * 3 - 2] = pt2i1;
	bb_1.p[k1 * 3 - 1] = pt3i1;
	cc_1.e[k1 - 1] = eti1;
	ee_1.lb[k1 - 1] = lb1;
	bb_1.p[k2 * 3 - 3] = pt1i2;
	bb_1.p[k2 * 3 - 2] = pt2i2;
	bb_1.p[k2 * 3 - 1] = pt3i2;
	cc_1.e[k2 - 1] = eti2;
	ee_1.lb[k2 - 1] = lb2;

	*iblock = 101;

	r__1 = pc_1.epion[nn_1.nnn + *irun * 150001 - 150002];

	r__2 = p2[0];

	r__3 = p2[1];

	r__4 = p2[2];
	epcmk = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	betak = p2[0] * bg_1.betax + p2[1] * bg_1.betay + p2[2] * bg_1.betaz;
	transf = bg_1.gamma * (bg_1.gamma * betak / (bg_1.gamma + 1.f) + epcmk);
	pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = bg_1.betax *
															   transf +
														   p2[0];
	pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = bg_1.betay *
															   transf +
														   p2[1];
	pb_1.ppion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = bg_1.betaz *
															   transf +
														   p2[2];

	dpert_1.dppion[nn_1.nnn + *irun * 150001 - 150002] = dpert_1.dpertp[*i1 -
																		1] *
														 dpert_1.dpertp[*i2 - 1];

	k = *i2;
	if (ee_1.lb[*i1 - 1] == 1 || ee_1.lb[*i1 - 1] == 2)
	{
		k = *i1;
	}
	pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450006] = aa_1.r__[k * 3 - 3];
	pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450005] = aa_1.r__[k * 3 - 2];
	pa_1.rpion[(nn_1.nnn + *irun * 150001) * 3 - 450004] = aa_1.r__[k * 3 - 1];
	return 0;
}

int pihypn_(int *ielstc, int *irun, int *iseed,
			float *dt, int *nt, int *ictrl, int *i1, int *i2, float *srt, float *pcx, float *pcy, float *pcz, int *nchrg, int *iblock)
{

	float r__1, r__2, r__3, r__4;

	double sqrt(double), cos(double), sin(double);

	static int i__, k1, k2;
	static float p1[4], p2[4], pk;
	static int lb1, lb2;
	static float dm3, dm4, fai, css, sss, e1cm, e2cm, eti1, eti2, pt1i1, pt2i1,
		pt3i1, pt1i2, pt2i2, pt3i2, px1cm, py1cm, pz1cm, pkmax, p1beta,
		p2beta;
	double ranart_(int *);
	static float transf;
	extern int rotate_(float *, float *, float *, float *, float *, float *);
	static float pxrota, pyrota, pzrota;

	px1cm = *pcx;
	py1cm = *pcy;
	pz1cm = *pcz;
	*ictrl = 1;

	if (*ielstc == 1)
	{

		k1 = *i1;
		k2 = *i2;
		dm3 = cc_1.e[k1 - 1];
		dm4 = cc_1.e[k2 - 1];
		*iblock = 10;
	}
	else
	{
		*iblock = 12;

		k1 = *i1;
		k2 = *i2;

		if (ee_1.lb[*i1 - 1] < 14 || ee_1.lb[*i1 - 1] > 17)
		{
			k1 = *i2;
			k2 = *i1;
		}

		ee_1.lb[k1 - 1] = (int)(ranart_(&rndf77_1.nseed) * 2) + 1;
		if (*nchrg == -2)
		{
			ee_1.lb[k1 - 1] = 6;
		}

		if (*nchrg == 2)
		{
			ee_1.lb[k1 - 1] = 9;
		}

		ee_1.lb[k2 - 1] = 21;
		dm3 = .938f;
		if (*nchrg == -2 || *nchrg == 1)
		{
			dm3 = 1.232f;
		}
		dm4 = .498f;
	}

	r__1 = *srt;

	r__2 = dm3 + dm4;

	r__3 = *srt;

	r__4 = dm3 - dm4;
	pkmax = sqrt((r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4)) /
			2.f / *srt;
	pk = pkmax;

	css = 1.f - ranart_(&rndf77_1.nseed) * 2.f;

	r__1 = css;
	sss = sqrt(1.f - r__1 * r__1);
	fai = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
	p1[0] = pk * sss * cos(fai);
	p1[1] = pk * sss * sin(fai);
	p1[2] = pk * css;
	for (i__ = 1; i__ <= 3; ++i__)
	{
		p2[i__ - 1] = p1[i__ - 1] * -1.f;
	}

	pxrota = p1[0];
	pyrota = p1[1];
	pzrota = p1[2];

	rotate_(pcx, pcy, pcz, &pxrota, &pyrota, &pzrota);
	p1[0] = pxrota;
	p1[1] = pyrota;
	p1[2] = pzrota;

	pxrota = p2[0];
	pyrota = p2[1];
	pzrota = p2[2];

	rotate_(pcx, pcy, pcz, &pxrota, &pyrota, &pzrota);
	p2[0] = pxrota;
	p2[1] = pyrota;
	p2[2] = pzrota;

	r__1 = dm3;

	r__2 = p1[0];

	r__3 = p1[1];

	r__4 = p1[2];
	e1cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	p1beta = p1[0] * bg_1.betax + p1[1] * bg_1.betay + p1[2] * bg_1.betaz;
	transf = bg_1.gamma * (bg_1.gamma * p1beta / (bg_1.gamma + 1) + e1cm);
	pt1i1 = bg_1.betax * transf + p1[0];
	pt2i1 = bg_1.betay * transf + p1[1];
	pt3i1 = bg_1.betaz * transf + p1[2];
	eti1 = dm3;
	lb1 = ee_1.lb[k1 - 1];

	r__1 = dm4;

	r__2 = p2[0];

	r__3 = p2[1];

	r__4 = p2[2];
	e2cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	p2beta = p2[0] * bg_1.betax + p2[1] * bg_1.betay + p2[2] * bg_1.betaz;
	transf = bg_1.gamma * (bg_1.gamma * p2beta / (bg_1.gamma + 1) + e2cm);
	pt1i2 = bg_1.betax * transf + p2[0];
	pt2i2 = bg_1.betay * transf + p2[1];
	pt3i2 = bg_1.betaz * transf + p2[2];

	eti2 = dm4;
	lb2 = ee_1.lb[k2 - 1];

	bb_1.p[k1 * 3 - 3] = pt1i1;
	bb_1.p[k1 * 3 - 2] = pt2i1;
	bb_1.p[k1 * 3 - 1] = pt3i1;
	cc_1.e[k1 - 1] = eti1;
	ee_1.lb[k1 - 1] = lb1;
	bb_1.p[k2 * 3 - 3] = pt1i2;
	bb_1.p[k2 * 3 - 2] = pt2i2;
	bb_1.p[k2 * 3 - 1] = pt3i2;
	cc_1.e[k2 - 1] = eti2;
	ee_1.lb[k2 - 1] = lb2;

	return 0;
}

int kaonn_(float *brel, float *brsgm, int *irun, int *iseed, float *dt, int *nt, int *ictrl, int *i1, int *i2, int *iblock, float *srt, float *pcx, float *pcy, float *pcz,
		   int *nchrg)
{

	float r__1, r__2, r__3, r__4;

	double sqrt(double), cos(double), sin(double);

	static int i__, k1, k2;
	static float p1[4], p2[4], pk;
	static int lb1, lb2;
	static float em1, em2, eee, fai, css, rrr, sss, e1cm, e2cm, eti1, eti2,
		pt1i1, pt2i1, pt3i1, pt1i2, pt2i2, pt3i2, px1cm, py1cm, pz1cm,
		pkmax, p1beta, p2beta;
	double ranart_(int *);
	static float transf;
	extern int rotate_(float *, float *, float *, float *, float *, float *);
	static float pxrota, pyrota, pzrota;

	px1cm = *pcx;
	py1cm = *pcy;
	pz1cm = *pcz;
	*ictrl = 1;

	k1 = *i1;
	k2 = *i2;

	if (cc_1.e[*i1 - 1] < .5f && cc_1.e[*i1 - 1] > .01f)
	{
		k1 = *i2;
		k2 = *i1;
	}

	eee = cc_1.e[k2 - 1];

	rrr = ranart_(&rndf77_1.nseed);
	if (rrr < *brel)
	{

		lb1 = ee_1.lb[k1 - 1];
		lb2 = ee_1.lb[k2 - 1];
		em1 = cc_1.e[k1 - 1];
		em2 = cc_1.e[k2 - 1];
		*iblock = 10;
	}
	else
	{
		*iblock = 12;
		if (rrr < *brel + *brsgm)
		{

			em1 = 1.1974f;
			em2 = .138f;

			lb1 = (int)(ranart_(&rndf77_1.nseed) * 3) + 15;
			lb2 = (int)(ranart_(&rndf77_1.nseed) * 3) + 3;
		}
		else
		{

			em1 = 1.1157f;
			em2 = .138f;

			lb1 = 14;

			lb2 = (int)(ranart_(&rndf77_1.nseed) * 3) + 3;
		}
	}
	ee_1.lb[k1 - 1] = lb1;
	ee_1.lb[k2 - 1] = lb2;

	r__1 = *srt;

	r__2 = em1 + em2;

	r__3 = *srt;

	r__4 = em1 - em2;
	pkmax = sqrt((r__1 * r__1 - r__2 * r__2) * (r__3 * r__3 - r__4 * r__4)) /
			2.f / *srt;
	pk = pkmax;

	css = 1.f - ranart_(&rndf77_1.nseed) * 2.f;

	r__1 = css;
	sss = sqrt(1.f - r__1 * r__1);
	fai = ranart_(&rndf77_1.nseed) * 6.2831852000000001f;
	p1[0] = pk * sss * cos(fai);
	p1[1] = pk * sss * sin(fai);
	p1[2] = pk * css;
	for (i__ = 1; i__ <= 3; ++i__)
	{
		p2[i__ - 1] = p1[i__ - 1] * -1.f;
	}

	pxrota = p1[0];
	pyrota = p1[1];
	pzrota = p1[2];

	rotate_(pcx, pcy, pcz, &pxrota, &pyrota, &pzrota);
	p1[0] = pxrota;
	p1[1] = pyrota;
	p1[2] = pzrota;

	pxrota = p2[0];
	pyrota = p2[1];
	pzrota = p2[2];

	rotate_(pcx, pcy, pcz, &pxrota, &pyrota, &pzrota);
	p2[0] = pxrota;
	p2[1] = pyrota;
	p2[2] = pzrota;

	r__1 = em1;

	r__2 = p1[0];

	r__3 = p1[1];

	r__4 = p1[2];
	e1cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	p1beta = p1[0] * bg_1.betax + p1[1] * bg_1.betay + p1[2] * bg_1.betaz;
	transf = bg_1.gamma * (bg_1.gamma * p1beta / (bg_1.gamma + 1) + e1cm);
	pt1i1 = bg_1.betax * transf + p1[0];
	pt2i1 = bg_1.betay * transf + p1[1];
	pt3i1 = bg_1.betaz * transf + p1[2];
	eti1 = em1;

	r__1 = em2;

	r__2 = p2[0];

	r__3 = p2[1];

	r__4 = p2[2];
	e2cm = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	p2beta = p2[0] * bg_1.betax + p2[1] * bg_1.betay + p2[2] * bg_1.betaz;
	transf = bg_1.gamma * (bg_1.gamma * p2beta / (bg_1.gamma + 1) + e2cm);
	pt1i2 = bg_1.betax * transf + p2[0];
	pt2i2 = bg_1.betay * transf + p2[1];
	pt3i2 = bg_1.betaz * transf + p2[2];
	eti2 = em2;

	bb_1.p[k1 * 3 - 3] = pt1i1;
	bb_1.p[k1 * 3 - 2] = pt2i1;
	bb_1.p[k1 * 3 - 1] = pt3i1;
	cc_1.e[k1 - 1] = eti1;
	bb_1.p[k2 * 3 - 3] = pt1i2;
	bb_1.p[k2 * 3 - 2] = pt2i2;
	bb_1.p[k2 * 3 - 1] = pt3i2;
	cc_1.e[k2 - 1] = eti2;

	return 0;
}

int artan1_(void)
{

	int i__1, i__2;
	float r__1, r__2, r__3;

	double sqrt(double), log(double);
	int s_wsle(cilist *), do_lio(int *, int *, char *, int),e_wsle(void);
	static int i__, j;
	static float y, ee;
	static int iy;
	static float xm, px, py, pz, eta;
	static int imt;
	static float xmt;
	static int ieta;
	static float dxmt;
	static int ityp;
	static float ptot;
	extern int luchge_(int *);
	static cilist io___302 = {0, 6, 0, 0, 0};
	static cilist io___303 = {0, 6, 0, 0, 0};
	static cilist io___304 = {0, 6, 0, 0, 0};
	static cilist io___305 = {0, 6, 0, 0, 0};
	static cilist io___306 = {0, 6, 0, 0, 0};

	i__1 = run_1.num;
	for (j = 1; j <= i__1; ++j)
	{
		i__2 = arerc1_1.multi1[j - 1];
		for (i__ = 1; i__ <= i__2; ++i__)
		{
			ityp = arprc1_1.ityp1[i__ + j * 150001 - 150002];
			px = arprc1_1.px1[i__ + j * 150001 - 150002];
			py = arprc1_1.py1[i__ + j * 150001 - 150002];
			pz = arprc1_1.pz1[i__ + j * 150001 - 150002];
			ee = arprc1_1.ee1[i__ + j * 150001 - 150002];
			xm = arprc1_1.xm1[i__ + j * 150001 - 150002];

			if (xm < .01f)
			{
				continue;
			}

			r__1 = px;
			r__2 = py;
			r__3 = pz;
			ptot = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
			eta = log((ptot + pz + 1e-5f) / (ptot - pz + 1e-5f)) * .5f;
			r__1 = px;
			r__2 = py;
			r__3 = xm;
			xmt = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
			if (dabs(pz) >= ee)
			{
				s_wsle(&io___302);
				do_lio(&c__9, &c__1, "IN ARTAN1", (int)9);
				e_wsle();
				s_wsle(&io___303);
				do_lio(&c__9, &c__1, "PARTICLE ", (int)9);
				do_lio(&c__3, &c__1, (char *)&i__, (int)sizeof(int));
				do_lio(&c__9, &c__1, " RUN ", (int)5);
				do_lio(&c__3, &c__1, (char *)&j, (int)sizeof(int));
				do_lio(&c__9, &c__1, "PREC ERR", (int)8);
				e_wsle();
				s_wsle(&io___304);
				do_lio(&c__9, &c__1, " FLAV = ", (int)8);
				do_lio(&c__3, &c__1, (char *)&ityp, (int)sizeof(int));
				do_lio(&c__9, &c__1, " PX = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&px, (int)sizeof(float));
				do_lio(&c__9, &c__1, " PY = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&py, (int)sizeof(float));
				e_wsle();
				s_wsle(&io___305);
				do_lio(&c__9, &c__1, " PZ = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&pz, (int)sizeof(float));
				do_lio(&c__9, &c__1, " EE = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&ee, (int)sizeof(float));
				e_wsle();
				s_wsle(&io___306);
				do_lio(&c__9, &c__1, " XM = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&xm, (int)sizeof(float));
				e_wsle();
			}
			dxmt = xmt - xm;
			y = log((ee + pz) / (ee - pz)) * .5f;

			if (dabs(y) >= 10.f)
			{
				if (y < -1.f || y > 1.f)
			{
                continue;
			}
			if (dxmt >= 2.5f || dxmt == 0.f)
			{
                continue;
			}
			imt = (int)(dxmt / .05f) + 1;
			if (ityp == 211)
			{
				arana1_1.dm1pip[imt - 1] += 1.f / xmt;
			}
			if (ityp == -211)
			{
				arana1_1.dm1pim[imt - 1] += 1.f / xmt;
			}
			if (ityp == 2212)
			{
				arana1_1.dmt1pr[imt - 1] += 1.f / xmt;
			}
			if (ityp == -2212)
			{
				arana1_1.dmt1pb[imt - 1] += 1.f / xmt;
			}
			if (ityp == 321)
			{
				arana1_1.dmt1kp[imt - 1] += 1.f / xmt;
			}
			if (ityp == -321)
			{
				arana1_1.dm1km[imt - 1] += 1.f / xmt;
			}

			if (ityp == 130)
			{
				arana1_1.dm1k0s[imt - 1] += 1.f / xmt;
			}
			if (ityp == 3122)
			{
				arana1_1.dmt1la[imt - 1] += 1.f / xmt;
			}
			if (ityp == -3122)
			{
				arana1_1.dmt1lb[imt - 1] += 1.f / xmt;
                }
			}
			if (dabs(eta) >= 10.f)
			{
				if (y < -1.f || y > 1.f)
			{
                continue;
			}
			if (dxmt >= 2.5f || dxmt == 0.f)
			{
                continue;
			}
			imt = (int)(dxmt / .05f) + 1;
			if (ityp == 211)
			{
				arana1_1.dm1pip[imt - 1] += 1.f / xmt;
			}
			if (ityp == -211)
			{
				arana1_1.dm1pim[imt - 1] += 1.f / xmt;
			}
			if (ityp == 2212)
			{
				arana1_1.dmt1pr[imt - 1] += 1.f / xmt;
			}
			if (ityp == -2212)
			{
				arana1_1.dmt1pb[imt - 1] += 1.f / xmt;
			}
			if (ityp == 321)
			{
				arana1_1.dmt1kp[imt - 1] += 1.f / xmt;
			}
			if (ityp == -321)
			{
				arana1_1.dm1km[imt - 1] += 1.f / xmt;
			}

			if (ityp == 130)
			{
				arana1_1.dm1k0s[imt - 1] += 1.f / xmt;
			}
			if (ityp == 3122)
			{
				arana1_1.dmt1la[imt - 1] += 1.f / xmt;
			}
			if (ityp == -3122)
			{
				arana1_1.dmt1lb[imt - 1] += 1.f / xmt;
                }
            }
			iy = (int)((y + 10.f) / .4f) + 1;
			ieta = (int)((eta + 10.f) / .4f) + 1;
			if (ityp < -1000)
			{
				arana1_1.dy1ntb[iy - 1] += -1.f;
			}
			if (ityp > 1000)
			{
				arana1_1.dy1ntb[iy - 1] += 1.f;
			}
			if (ityp == -2212)
			{
				arana1_1.dy1ntp[iy - 1] += -1.f;
			}
			if (ityp == 2212)
			{
				arana1_1.dy1ntp[iy - 1] += 1.f;
			}

			if (ityp == -2112)
			{
				arana1_1.dy1hm[iy - 1] += 1.f;
			}

			if (luchge_(&ityp) != 0)
			{
				arana1_1.dy1ch[iy - 1] += 1.f;
				arana1_1.de1ch[ieta - 1] += 1.f;
				if (luchge_(&ityp) < 0)
				{
					arana1_1.dy1neg[iy - 1] += 1.f;
					arana1_1.de1neg[ieta - 1] += 1.f;
				}
			}

			if (ityp >= 100 && ityp < 1000 || ityp > -1000 && ityp <= -100)
			{
				arana1_1.dy1msn[iy - 1] += 1.f;
			}
			if (ityp == 211)
			{
				arana1_1.dy1pip[iy - 1] += 1.f;
			}
			if (ityp == -211)
			{
				arana1_1.dy1pim[iy - 1] += 1.f;
			}
			if (ityp == 111)
			{
				arana1_1.dy1pi0[iy - 1] += 1.f;
			}
			if (ityp == 2212)
			{
				arana1_1.dy1pr[iy - 1] += 1.f;
			}
			if (ityp == -2212)
			{
				arana1_1.dy1pb[iy - 1] += 1.f;
			}

			if (ityp == 321)
			{
				arana1_1.dy1kp[iy - 1] += 1.f;
			}
			if (ityp == -321)
			{
				arana1_1.dy1km[iy - 1] += 1.f;
			}

			if (ityp == 130)
			{
				arana1_1.dy1k0s[iy - 1] += 1.f;
			}
			if (ityp == 3122)
			{
				arana1_1.dy1la[iy - 1] += 1.f;
			}
			if (ityp == -3122)
			{
				arana1_1.dy1lb[iy - 1] += 1.f;
			}
			if (ityp == 333)
			{
				arana1_1.dy1phi[iy - 1] += 1.f;
			}
		}
	}
	return 0;
}

int artan2_(void)
{

	int i__1, i__2;
	float r__1, r__2, r__3;
	double sqrt(double), log(double);
	int s_wsle(cilist *), do_lio(int *, int *, char *, int),e_wsle(void);
	static int i__, j;
	static float y, ee;
	static int iy;
	static float xm, px, py, pz, eta;
	static int imt;
	static float xmt;
	static int ieta;
	static float dxmt;
	static int ityp;
	static float ptot;
	extern int luchge_(int *);
	static cilist io___323 = {0, 6, 0, 0, 0};
	static cilist io___324 = {0, 6, 0, 0, 0};
	static cilist io___325 = {0, 6, 0, 0, 0};
	static cilist io___326 = {0, 6, 0, 0, 0};
	static cilist io___327 = {0, 6, 0, 0, 0};
	i__1 = run_1.num;
	for (j = 1; j <= i__1; ++j)
	{
		i__2 = arerc1_1.multi1[j - 1];
		for (i__ = 1; i__ <= i__2; ++i__)
		{
			ityp = arprc1_1.ityp1[i__ + j * 150001 - 150002];
			px = arprc1_1.px1[i__ + j * 150001 - 150002];
			py = arprc1_1.py1[i__ + j * 150001 - 150002];
			pz = arprc1_1.pz1[i__ + j * 150001 - 150002];
			ee = arprc1_1.ee1[i__ + j * 150001 - 150002];
			xm = arprc1_1.xm1[i__ + j * 150001 - 150002];
			r__1 = px;
			r__2 = py;
			r__3 = xm;
			xmt = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
			if (xm < .01f)
			{
				continue;
			}

			r__1 = px;
			r__2 = py;
			r__3 = pz;
			ptot = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
			eta = log((ptot + pz + 1e-5f) / (ptot - pz + 1e-5f)) * .5f;
			if (dabs(pz) >= ee)
			{
				s_wsle(&io___323);
				do_lio(&c__9, &c__1, "IN ARTAN2", (int)9);
				e_wsle();
				s_wsle(&io___324);
				do_lio(&c__9, &c__1, "PARTICLE ", (int)9);
				do_lio(&c__3, &c__1, (char *)&i__, (int)sizeof(int));
				do_lio(&c__9, &c__1, " RUN ", (int)5);
				do_lio(&c__3, &c__1, (char *)&j, (int)sizeof(int));
				do_lio(&c__9, &c__1, "PREC ERR", (int)8);
				e_wsle();
				s_wsle(&io___325);
				do_lio(&c__9, &c__1, " FLAV = ", (int)8);
				do_lio(&c__3, &c__1, (char *)&ityp, (int)sizeof(int));
				do_lio(&c__9, &c__1, " PX = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&px, (int)sizeof(float));
				do_lio(&c__9, &c__1, " PY = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&py, (int)sizeof(float));
				e_wsle();
				s_wsle(&io___326);
				do_lio(&c__9, &c__1, " PZ = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&pz, (int)sizeof(float));
				do_lio(&c__9, &c__1, " EE = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&ee, (int)sizeof(float));
				e_wsle();
				s_wsle(&io___327);
				do_lio(&c__9, &c__1, " XM = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&xm, (int)sizeof(float));
				e_wsle();
			}
			dxmt = xmt - xm;
			y = log((ee + pz) / (ee - pz)) * .5f;

			if (dabs(y) >= 10.f)
			{
                if (y < -1.f || y > 1.f)
			{
				continue;
			}
			if (dxmt >= 2.5f || dxmt == 0.f)
			{
				continue;
			}
			imt = (int)(dxmt / .05f) + 1;
			if (ityp == 211)
			{
				arana2_1.dm2pip[imt - 1] += 1.f / xmt;
			}
			if (ityp == -211)
			{
				arana2_1.dm2pim[imt - 1] += 1.f / xmt;
			}
			if (ityp == 2212)
			{
				arana2_1.dmt2pr[imt - 1] += 1.f / xmt;
			}
			if (ityp == -2212)
			{
				arana2_1.dmt2pb[imt - 1] += 1.f / xmt;
			}
			if (ityp == 321)
			{
				arana2_1.dmt2kp[imt - 1] += 1.f / xmt;
			}
			if (ityp == -321)
			{
				arana2_1.dm2km[imt - 1] += 1.f / xmt;
			}

			if (ityp == 130)
			{
				arana2_1.dm2k0s[imt - 1] += 1.f / xmt;
			}
			if (ityp == 3122)
			{
				arana2_1.dmt2la[imt - 1] += 1.f / xmt;
			}
			if (ityp == -3122)
			{
				arana2_1.dmt2lb[imt - 1] += 1.f / xmt;
                }
			}

			if (dabs(eta) >= 10.f)
			{
				if (y < -1.f || y > 1.f)
			{
				continue;
			}
			if (dxmt >= 2.5f || dxmt == 0.f)
			{
				continue;
			}
			imt = (int)(dxmt / .05f) + 1;
			if (ityp == 211)
			{
				arana2_1.dm2pip[imt - 1] += 1.f / xmt;
			}
			if (ityp == -211)
			{
				arana2_1.dm2pim[imt - 1] += 1.f / xmt;
			}
			if (ityp == 2212)
			{
				arana2_1.dmt2pr[imt - 1] += 1.f / xmt;
			}
			if (ityp == -2212)
			{
				arana2_1.dmt2pb[imt - 1] += 1.f / xmt;
			}
			if (ityp == 321)
			{
				arana2_1.dmt2kp[imt - 1] += 1.f / xmt;
			}
			if (ityp == -321)
			{
				arana2_1.dm2km[imt - 1] += 1.f / xmt;
			}

			if (ityp == 130)
			{
				arana2_1.dm2k0s[imt - 1] += 1.f / xmt;
			}
			if (ityp == 3122)
			{
				arana2_1.dmt2la[imt - 1] += 1.f / xmt;
			}
			if (ityp == -3122)
			{
				arana2_1.dmt2lb[imt - 1] += 1.f / xmt;
			}
			}
			iy = (int)((y + 10.f) / .4f) + 1;
			ieta = (int)((eta + 10.f) / .4f) + 1;
			if (ityp < -1000)
			{
				arana2_1.dy2ntb[iy - 1] += -1.f;
			}
			if (ityp > 1000)
			{
				arana2_1.dy2ntb[iy - 1] += 1.f;
			}
			if (ityp == -2212)
			{
				arana2_1.dy2ntp[iy - 1] += -1.f;
			}
			if (ityp == 2212)
			{
				arana2_1.dy2ntp[iy - 1] += 1.f;
			}
			if (ityp == -2112)
			{
				arana2_1.dy2hm[iy - 1] += 1.f;
			}
			if (luchge_(&ityp) != 0)
			{
				arana2_1.dy2ch[iy - 1] += 1.f;
				arana2_1.de2ch[ieta - 1] += 1.f;
				if (luchge_(&ityp) < 0)
				{
					arana2_1.dy2neg[iy - 1] += 1.f;
					arana2_1.de2neg[ieta - 1] += 1.f;
				}
			}

			if (ityp >= 100 && ityp < 1000 || ityp > -1000 && ityp <= -100)
			{
				arana2_1.dy2msn[iy - 1] += 1.f;
			}
			if (ityp == 211)
			{
				arana2_1.dy2pip[iy - 1] += 1.f;
			}
			if (ityp == -211)
			{
				arana2_1.dy2pim[iy - 1] += 1.f;
			}
			if (ityp == 111)
			{
				arana2_1.dy2pi0[iy - 1] += 1.f;
			}
			if (ityp == 2212)
			{
				arana2_1.dy2pr[iy - 1] += 1.f;
			}
			if (ityp == -2212)
			{
				arana2_1.dy2pb[iy - 1] += 1.f;
			}

			if (ityp == 321)
			{
				arana2_1.dy2kp[iy - 1] += 1.f;
			}
			if (ityp == -321)
			{
				arana2_1.dy2km[iy - 1] += 1.f;
			}

			if (ityp == 130)
			{
				arana2_1.dy2k0s[iy - 1] += 1.f;
			}
			if (ityp == 3122)
			{
				arana2_1.dy2la[iy - 1] += 1.f;
			}
			if (ityp == -3122)
			{
				arana2_1.dy2lb[iy - 1] += 1.f;
			}
			if (ityp == 333)
			{
				arana2_1.dy2phi[iy - 1] += 1.f;
			}
		}
	}
	return 0;
}

int artout_(int *nevnt)
{

	static char fmt_333[] = "(2(f12.5,1x))";
	float r__1, r__2;
	olist o__1;
	int f_open(olist *), s_wsfe(cilist *), do_fio(int *, char *, int), e_wsfe(void), s_wsle(cilist *), do_lio(int *, int *, char *, int), e_wsle(void);
	static int i__;
	static float ymid, scale1, scale2;
	static cilist io___337 = {0, 30, 0, fmt_333, 0};
	static cilist io___338 = {0, 31, 0, fmt_333, 0};
	static cilist io___339 = {0, 32, 0, fmt_333, 0};
	static cilist io___340 = {0, 37, 0, fmt_333, 0};
	static cilist io___341 = {0, 38, 0, fmt_333, 0};
	static cilist io___342 = {0, 39, 0, fmt_333, 0};
	static cilist io___343 = {0, 40, 0, fmt_333, 0};
	static cilist io___344 = {0, 41, 0, fmt_333, 0};
	static cilist io___345 = {0, 42, 0, fmt_333, 0};
	static cilist io___346 = {0, 33, 0, fmt_333, 0};
	static cilist io___347 = {0, 34, 0, fmt_333, 0};
	static cilist io___348 = {0, 35, 0, fmt_333, 0};
	static cilist io___349 = {0, 36, 0, fmt_333, 0};
	static cilist io___350 = {0, 43, 0, fmt_333, 0};
	static cilist io___351 = {0, 44, 0, fmt_333, 0};
	static cilist io___352 = {0, 45, 0, fmt_333, 0};
	static cilist io___353 = {0, 46, 0, fmt_333, 0};
	static cilist io___354 = {0, 47, 0, fmt_333, 0};
	static cilist io___355 = {0, 48, 0, fmt_333, 0};
	static cilist io___356 = {0, 50, 0, fmt_333, 0};
	static cilist io___357 = {0, 51, 0, fmt_333, 0};
	static cilist io___358 = {0, 52, 0, fmt_333, 0};
	static cilist io___359 = {0, 53, 0, fmt_333, 0};
	static cilist io___360 = {0, 54, 0, fmt_333, 0};
	static cilist io___361 = {0, 55, 0, fmt_333, 0};
	static cilist io___362 = {0, 56, 0, fmt_333, 0};
	static cilist io___363 = {0, 57, 0, fmt_333, 0};
	static cilist io___364 = {0, 58, 0, fmt_333, 0};
	static cilist io___365 = {0, 0, 0, 0, 0};
	static cilist io___366 = {0, 0, 0, 0, 0};
	static cilist io___367 = {0, 30, 0, fmt_333, 0};
	static cilist io___368 = {0, 31, 0, fmt_333, 0};
	static cilist io___369 = {0, 32, 0, fmt_333, 0};
	static cilist io___370 = {0, 37, 0, fmt_333, 0};
	static cilist io___371 = {0, 38, 0, fmt_333, 0};
	static cilist io___372 = {0, 39, 0, fmt_333, 0};
	static cilist io___373 = {0, 40, 0, fmt_333, 0};
	static cilist io___374 = {0, 41, 0, fmt_333, 0};
	static cilist io___375 = {0, 42, 0, fmt_333, 0};
	static cilist io___376 = {0, 33, 0, fmt_333, 0};
	static cilist io___377 = {0, 34, 0, fmt_333, 0};
	static cilist io___378 = {0, 35, 0, fmt_333, 0};
	static cilist io___379 = {0, 36, 0, fmt_333, 0};
	static cilist io___380 = {0, 43, 0, fmt_333, 0};
	static cilist io___381 = {0, 44, 0, fmt_333, 0};
	static cilist io___382 = {0, 45, 0, fmt_333, 0};
	static cilist io___383 = {0, 46, 0, fmt_333, 0};
	static cilist io___384 = {0, 47, 0, fmt_333, 0};
	static cilist io___385 = {0, 48, 0, fmt_333, 0};
	static cilist io___386 = {0, 50, 0, fmt_333, 0};
	static cilist io___387 = {0, 51, 0, fmt_333, 0};
	static cilist io___388 = {0, 52, 0, fmt_333, 0};
	static cilist io___389 = {0, 53, 0, fmt_333, 0};
	static cilist io___390 = {0, 54, 0, fmt_333, 0};
	static cilist io___391 = {0, 55, 0, fmt_333, 0};
	static cilist io___392 = {0, 56, 0, fmt_333, 0};
	static cilist io___393 = {0, 57, 0, fmt_333, 0};
	static cilist io___394 = {0, 58, 0, fmt_333, 0};

	o__1.oerr = 0;
	o__1.ounit = 30;
	o__1.ofnmlen = 17;
	o__1.ofnm = "ana/dndy_netb.dat";
	o__1.orl = 0;
	o__1.osta = "UNKNOWN";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	o__1.oerr = 0;
	o__1.ounit = 31;
	o__1.ofnmlen = 17;
	o__1.ofnm = "ana/dndy_netp.dat";
	o__1.orl = 0;
	o__1.osta = "UNKNOWN";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	o__1.oerr = 0;
	o__1.ounit = 32;
	o__1.ofnmlen = 15;
	o__1.ofnm = "ana/dndy_nb.dat";
	o__1.orl = 0;
	o__1.osta = "UNKNOWN";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	o__1.oerr = 0;
	o__1.ounit = 33;
	o__1.ofnmlen = 16;
	o__1.ofnm = "ana/dndy_neg.dat";
	o__1.orl = 0;
	o__1.osta = "UNKNOWN";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	o__1.oerr = 0;
	o__1.ounit = 34;
	o__1.ofnmlen = 15;
	o__1.ofnm = "ana/dndy_ch.dat";
	o__1.orl = 0;
	o__1.osta = "UNKNOWN";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	o__1.oerr = 0;
	o__1.ounit = 35;
	o__1.ofnmlen = 16;
	o__1.ofnm = "ana/dnde_neg.dat";
	o__1.orl = 0;
	o__1.osta = "UNKNOWN";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	o__1.oerr = 0;
	o__1.ounit = 36;
	o__1.ofnmlen = 15;
	o__1.ofnm = "ana/dnde_ch.dat";
	o__1.orl = 0;
	o__1.osta = "UNKNOWN";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	o__1.oerr = 0;
	o__1.ounit = 37;
	o__1.ofnmlen = 15;
	o__1.ofnm = "ana/dndy_kp.dat";
	o__1.orl = 0;
	o__1.osta = "UNKNOWN";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	o__1.oerr = 0;
	o__1.ounit = 38;
	o__1.ofnmlen = 15;
	o__1.ofnm = "ana/dndy_km.dat";
	o__1.orl = 0;
	o__1.osta = "UNKNOWN";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);

	o__1.oerr = 0;
	o__1.ounit = 39;
	o__1.ofnmlen = 16;
	o__1.ofnm = "ana/dndy_k0l.dat";
	o__1.orl = 0;
	o__1.osta = "UNKNOWN";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	o__1.oerr = 0;
	o__1.ounit = 40;
	o__1.ofnmlen = 15;
	o__1.ofnm = "ana/dndy_la.dat";
	o__1.orl = 0;
	o__1.osta = "UNKNOWN";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	o__1.oerr = 0;
	o__1.ounit = 41;
	o__1.ofnmlen = 15;
	o__1.ofnm = "ana/dndy_lb.dat";
	o__1.orl = 0;
	o__1.osta = "UNKNOWN";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	o__1.oerr = 0;
	o__1.ounit = 42;
	o__1.ofnmlen = 16;
	o__1.ofnm = "ana/dndy_phi.dat";
	o__1.orl = 0;
	o__1.osta = "UNKNOWN";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);

	o__1.oerr = 0;
	o__1.ounit = 43;
	o__1.ofnmlen = 18;
	o__1.ofnm = "ana/dndy_meson.dat";
	o__1.orl = 0;
	o__1.osta = "UNKNOWN";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	o__1.oerr = 0;
	o__1.ounit = 44;
	o__1.ofnmlen = 16;
	o__1.ofnm = "ana/dndy_pip.dat";
	o__1.orl = 0;
	o__1.osta = "UNKNOWN";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	o__1.oerr = 0;
	o__1.ounit = 45;
	o__1.ofnmlen = 16;
	o__1.ofnm = "ana/dndy_pim.dat";
	o__1.orl = 0;
	o__1.osta = "UNKNOWN";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	o__1.oerr = 0;
	o__1.ounit = 46;
	o__1.ofnmlen = 16;
	o__1.ofnm = "ana/dndy_pi0.dat";
	o__1.orl = 0;
	o__1.osta = "UNKNOWN";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	o__1.oerr = 0;
	o__1.ounit = 47;
	o__1.ofnmlen = 15;
	o__1.ofnm = "ana/dndy_pr.dat";
	o__1.orl = 0;
	o__1.osta = "UNKNOWN";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	o__1.oerr = 0;
	o__1.ounit = 48;
	o__1.ofnmlen = 15;
	o__1.ofnm = "ana/dndy_pb.dat";
	o__1.orl = 0;
	o__1.osta = "UNKNOWN";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);

	o__1.oerr = 0;
	o__1.ounit = 50;
	o__1.ofnmlen = 19;
	o__1.ofnm = "ana/dndmtdy_pip.dat";
	o__1.orl = 0;
	o__1.osta = "UNKNOWN";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	o__1.oerr = 0;
	o__1.ounit = 51;
	o__1.ofnmlen = 23;
	o__1.ofnm = "ana/dndmtdy_0_1_pim.dat";
	o__1.orl = 0;
	o__1.osta = "UNKNOWN";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	o__1.oerr = 0;
	o__1.ounit = 52;
	o__1.ofnmlen = 18;
	o__1.ofnm = "ana/dndmtdy_pr.dat";
	o__1.orl = 0;
	o__1.osta = "UNKNOWN";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	o__1.oerr = 0;
	o__1.ounit = 53;
	o__1.ofnmlen = 18;
	o__1.ofnm = "ana/dndmtdy_pb.dat";
	o__1.orl = 0;
	o__1.osta = "UNKNOWN";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	o__1.oerr = 0;
	o__1.ounit = 54;
	o__1.ofnmlen = 18;
	o__1.ofnm = "ana/dndmtdy_kp.dat";
	o__1.orl = 0;
	o__1.osta = "UNKNOWN";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	o__1.oerr = 0;
	o__1.ounit = 55;
	o__1.ofnmlen = 22;
	o__1.ofnm = "ana/dndmtdy_0_5_km.dat";
	o__1.orl = 0;
	o__1.osta = "UNKNOWN";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	o__1.oerr = 0;
	o__1.ounit = 56;
	o__1.ofnmlen = 19;
	o__1.ofnm = "ana/dndmtdy_k0s.dat";
	o__1.orl = 0;
	o__1.osta = "UNKNOWN";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	o__1.oerr = 0;
	o__1.ounit = 57;
	o__1.ofnmlen = 18;
	o__1.ofnm = "ana/dndmtdy_la.dat";
	o__1.orl = 0;
	o__1.osta = "UNKNOWN";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
	o__1.oerr = 0;
	o__1.ounit = 58;
	o__1.ofnmlen = 18;
	o__1.ofnm = "ana/dndmtdy_lb.dat";
	o__1.orl = 0;
	o__1.osta = "UNKNOWN";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);

	scale1 = 1.f / (float)(*nevnt * run_1.num) / .4f;
	scale2 = 1.f / (float)(*nevnt * run_1.num) / .05f / 2.f;

	for (i__ = 1; i__ <= 50; ++i__)
	{
		ymid = (i__ - .5f) * .4f - 10.f;
		s_wsfe(&io___337);
		do_fio(&c__1, (char *)&ymid, (int)sizeof(float));
		r__1 = scale1 * arana1_1.dy1ntb[i__ - 1];
		do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
		e_wsfe();
		s_wsfe(&io___338);
		do_fio(&c__1, (char *)&ymid, (int)sizeof(float));
		r__1 = scale1 * arana1_1.dy1ntp[i__ - 1];
		do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
		e_wsfe();
		s_wsfe(&io___339);
		do_fio(&c__1, (char *)&ymid, (int)sizeof(float));
		r__1 = scale1 * arana1_1.dy1hm[i__ - 1];
		do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
		e_wsfe();
		s_wsfe(&io___340);
		do_fio(&c__1, (char *)&ymid, (int)sizeof(float));
		r__1 = scale1 * arana1_1.dy1kp[i__ - 1];
		do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
		e_wsfe();
		s_wsfe(&io___341);
		do_fio(&c__1, (char *)&ymid, (int)sizeof(float));
		r__1 = scale1 * arana1_1.dy1km[i__ - 1];
		do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
		e_wsfe();
		s_wsfe(&io___342);
		do_fio(&c__1, (char *)&ymid, (int)sizeof(float));
		r__1 = scale1 * arana1_1.dy1k0s[i__ - 1];
		do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
		e_wsfe();
		s_wsfe(&io___343);
		do_fio(&c__1, (char *)&ymid, (int)sizeof(float));
		r__1 = scale1 * arana1_1.dy1la[i__ - 1];
		do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
		e_wsfe();
		s_wsfe(&io___344);
		do_fio(&c__1, (char *)&ymid, (int)sizeof(float));
		r__1 = scale1 * arana1_1.dy1lb[i__ - 1];
		do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
		e_wsfe();
		s_wsfe(&io___345);
		do_fio(&c__1, (char *)&ymid, (int)sizeof(float));
		r__1 = scale1 * arana1_1.dy1phi[i__ - 1];
		do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
		e_wsfe();
		s_wsfe(&io___346);
		do_fio(&c__1, (char *)&ymid, (int)sizeof(float));
		r__1 = scale1 * arana1_1.dy1neg[i__ - 1];
		do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
		e_wsfe();
		s_wsfe(&io___347);
		do_fio(&c__1, (char *)&ymid, (int)sizeof(float));
		r__1 = scale1 * arana1_1.dy1ch[i__ - 1];
		do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
		e_wsfe();
		s_wsfe(&io___348);
		do_fio(&c__1, (char *)&ymid, (int)sizeof(float));
		r__1 = scale1 * arana1_1.de1neg[i__ - 1];
		do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
		e_wsfe();
		s_wsfe(&io___349);
		do_fio(&c__1, (char *)&ymid, (int)sizeof(float));
		r__1 = scale1 * arana1_1.de1ch[i__ - 1];
		do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
		e_wsfe();
		s_wsfe(&io___350);
		do_fio(&c__1, (char *)&ymid, (int)sizeof(float));
		r__1 = scale1 * arana1_1.dy1msn[i__ - 1];
		do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
		e_wsfe();
		s_wsfe(&io___351);
		do_fio(&c__1, (char *)&ymid, (int)sizeof(float));
		r__1 = scale1 * arana1_1.dy1pip[i__ - 1];
		do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
		e_wsfe();
		s_wsfe(&io___352);
		do_fio(&c__1, (char *)&ymid, (int)sizeof(float));
		r__1 = scale1 * arana1_1.dy1pim[i__ - 1];
		do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
		e_wsfe();
		s_wsfe(&io___353);
		do_fio(&c__1, (char *)&ymid, (int)sizeof(float));
		r__1 = scale1 * arana1_1.dy1pi0[i__ - 1];
		do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
		e_wsfe();
		s_wsfe(&io___354);
		do_fio(&c__1, (char *)&ymid, (int)sizeof(float));
		r__1 = scale1 * arana1_1.dy1pr[i__ - 1];
		do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
		e_wsfe();
		s_wsfe(&io___355);
		do_fio(&c__1, (char *)&ymid, (int)sizeof(float));
		r__1 = scale1 * arana1_1.dy1pb[i__ - 1];
		do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
		e_wsfe();
		if (arana1_1.dm1pip[i__ - 1] != 0.f)
		{
			s_wsfe(&io___356);
			r__1 = (i__ - .5f) * .05f;
			do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
			r__2 = scale2 * arana1_1.dm1pip[i__ - 1];
			do_fio(&c__1, (char *)&r__2, (int)sizeof(float));
			e_wsfe();
		}
		if (arana1_1.dm1pim[i__ - 1] != 0.f)
		{
			s_wsfe(&io___357);
			r__1 = (i__ - .5f) * .05f;
			do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
			r__2 = scale2 * .1f * arana1_1.dm1pim[i__ - 1];
			do_fio(&c__1, (char *)&r__2, (int)sizeof(float));
			e_wsfe();
		}
		if (arana1_1.dmt1pr[i__ - 1] != 0.f)
		{
			s_wsfe(&io___358);
			r__1 = (i__ - .5f) * .05f;
			do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
			r__2 = scale2 * arana1_1.dmt1pr[i__ - 1];
			do_fio(&c__1, (char *)&r__2, (int)sizeof(float));
			e_wsfe();
		}
		if (arana1_1.dmt1pb[i__ - 1] != 0.f)
		{
			s_wsfe(&io___359);
			r__1 = (i__ - .5f) * .05f;
			do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
			r__2 = scale2 * arana1_1.dmt1pb[i__ - 1];
			do_fio(&c__1, (char *)&r__2, (int)sizeof(float));
			e_wsfe();
		}
		if (arana1_1.dmt1kp[i__ - 1] != 0.f)
		{
			s_wsfe(&io___360);
			r__1 = (i__ - .5f) * .05f;
			do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
			r__2 = scale2 * arana1_1.dmt1kp[i__ - 1];
			do_fio(&c__1, (char *)&r__2, (int)sizeof(float));
			e_wsfe();
		}
		if (arana1_1.dm1km[i__ - 1] != 0.f)
		{
			s_wsfe(&io___361);
			r__1 = (i__ - .5f) * .05f;
			do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
			r__2 = scale2 * .5f * arana1_1.dm1km[i__ - 1];
			do_fio(&c__1, (char *)&r__2, (int)sizeof(float));
			e_wsfe();
		}
		if (arana1_1.dm1k0s[i__ - 1] != 0.f)
		{
			s_wsfe(&io___362);
			r__1 = (i__ - .5f) * .05f;
			do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
			r__2 = scale2 * arana1_1.dm1k0s[i__ - 1];
			do_fio(&c__1, (char *)&r__2, (int)sizeof(float));
			e_wsfe();
		}
		if (arana1_1.dmt1la[i__ - 1] != 0.f)
		{
			s_wsfe(&io___363);
			r__1 = (i__ - .5f) * .05f;
			do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
			r__2 = scale2 * arana1_1.dmt1la[i__ - 1];
			do_fio(&c__1, (char *)&r__2, (int)sizeof(float));
			e_wsfe();
		}
		if (arana1_1.dmt1lb[i__ - 1] != 0.f)
		{
			s_wsfe(&io___364);
			r__1 = (i__ - .5f) * .05f;
			do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
			r__2 = scale2 * arana1_1.dmt1lb[i__ - 1];
			do_fio(&c__1, (char *)&r__2, (int)sizeof(float));
			e_wsfe();
		}
	}

	for (i__ = 30; i__ <= 48; ++i__)
	{
		io___365.ciunit = i__;
		s_wsle(&io___365);
		do_lio(&c__9, &c__1, "after hadron evolution", (int)22);
		e_wsle();
	}
	for (i__ = 50; i__ <= 58; ++i__)
	{
		io___366.ciunit = i__;
		s_wsle(&io___366);
		do_lio(&c__9, &c__1, "after hadron evolution", (int)22);
		e_wsle();
	}
	for (i__ = 1; i__ <= 50; ++i__)
	{
		ymid = (i__ - .5f) * .4f - 10.f;
		s_wsfe(&io___367);
		do_fio(&c__1, (char *)&ymid, (int)sizeof(float));
		r__1 = scale1 * arana2_1.dy2ntb[i__ - 1];
		do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
		e_wsfe();
		s_wsfe(&io___368);
		do_fio(&c__1, (char *)&ymid, (int)sizeof(float));
		r__1 = scale1 * arana2_1.dy2ntp[i__ - 1];
		do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
		e_wsfe();
		s_wsfe(&io___369);
		do_fio(&c__1, (char *)&ymid, (int)sizeof(float));
		r__1 = scale1 * arana2_1.dy2hm[i__ - 1];
		do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
		e_wsfe();
		s_wsfe(&io___370);
		do_fio(&c__1, (char *)&ymid, (int)sizeof(float));
		r__1 = scale1 * arana2_1.dy2kp[i__ - 1];
		do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
		e_wsfe();
		s_wsfe(&io___371);
		do_fio(&c__1, (char *)&ymid, (int)sizeof(float));
		r__1 = scale1 * arana2_1.dy2km[i__ - 1];
		do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
		e_wsfe();
		s_wsfe(&io___372);
		do_fio(&c__1, (char *)&ymid, (int)sizeof(float));
		r__1 = scale1 * arana2_1.dy2k0s[i__ - 1];
		do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
		e_wsfe();
		s_wsfe(&io___373);
		do_fio(&c__1, (char *)&ymid, (int)sizeof(float));
		r__1 = scale1 * arana2_1.dy2la[i__ - 1];
		do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
		e_wsfe();
		s_wsfe(&io___374);
		do_fio(&c__1, (char *)&ymid, (int)sizeof(float));
		r__1 = scale1 * arana2_1.dy2lb[i__ - 1];
		do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
		e_wsfe();
		s_wsfe(&io___375);
		do_fio(&c__1, (char *)&ymid, (int)sizeof(float));
		r__1 = scale1 * arana2_1.dy2phi[i__ - 1];
		do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
		e_wsfe();
		s_wsfe(&io___376);
		do_fio(&c__1, (char *)&ymid, (int)sizeof(float));
		r__1 = scale1 * arana2_1.dy2neg[i__ - 1];
		do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
		e_wsfe();
		s_wsfe(&io___377);
		do_fio(&c__1, (char *)&ymid, (int)sizeof(float));
		r__1 = scale1 * arana2_1.dy2ch[i__ - 1];
		do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
		e_wsfe();
		s_wsfe(&io___378);
		do_fio(&c__1, (char *)&ymid, (int)sizeof(float));
		r__1 = scale1 * arana2_1.de2neg[i__ - 1];
		do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
		e_wsfe();
		s_wsfe(&io___379);
		do_fio(&c__1, (char *)&ymid, (int)sizeof(float));
		r__1 = scale1 * arana2_1.de2ch[i__ - 1];
		do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
		e_wsfe();
		s_wsfe(&io___380);
		do_fio(&c__1, (char *)&ymid, (int)sizeof(float));
		r__1 = scale1 * arana2_1.dy2msn[i__ - 1];
		do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
		e_wsfe();
		s_wsfe(&io___381);
		do_fio(&c__1, (char *)&ymid, (int)sizeof(float));
		r__1 = scale1 * arana2_1.dy2pip[i__ - 1];
		do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
		e_wsfe();
		s_wsfe(&io___382);
		do_fio(&c__1, (char *)&ymid, (int)sizeof(float));
		r__1 = scale1 * arana2_1.dy2pim[i__ - 1];
		do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
		e_wsfe();
		s_wsfe(&io___383);
		do_fio(&c__1, (char *)&ymid, (int)sizeof(float));
		r__1 = scale1 * arana2_1.dy2pi0[i__ - 1];
		do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
		e_wsfe();
		s_wsfe(&io___384);
		do_fio(&c__1, (char *)&ymid, (int)sizeof(float));
		r__1 = scale1 * arana2_1.dy2pr[i__ - 1];
		do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
		e_wsfe();
		s_wsfe(&io___385);
		do_fio(&c__1, (char *)&ymid, (int)sizeof(float));
		r__1 = scale1 * arana2_1.dy2pb[i__ - 1];
		do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
		e_wsfe();

		if (arana2_1.dm2pip[i__ - 1] != 0.f)
		{
			s_wsfe(&io___386);
			r__1 = (i__ - .5f) * .05f;
			do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
			r__2 = scale2 * arana2_1.dm2pip[i__ - 1];
			do_fio(&c__1, (char *)&r__2, (int)sizeof(float));
			e_wsfe();
		}
		if (arana2_1.dm2pim[i__ - 1] != 0.f)
		{
			s_wsfe(&io___387);
			r__1 = (i__ - .5f) * .05f;
			do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
			r__2 = scale2 * .1f * arana2_1.dm2pim[i__ - 1];
			do_fio(&c__1, (char *)&r__2, (int)sizeof(float));
			e_wsfe();
		}
		if (arana2_1.dmt2pr[i__ - 1] != 0.f)
		{
			s_wsfe(&io___388);
			r__1 = (i__ - .5f) * .05f;
			do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
			r__2 = scale2 * arana2_1.dmt2pr[i__ - 1];
			do_fio(&c__1, (char *)&r__2, (int)sizeof(float));
			e_wsfe();
		}
		if (arana2_1.dmt2pb[i__ - 1] != 0.f)
		{
			s_wsfe(&io___389);
			r__1 = (i__ - .5f) * .05f;
			do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
			r__2 = scale2 * arana2_1.dmt2pb[i__ - 1];
			do_fio(&c__1, (char *)&r__2, (int)sizeof(float));
			e_wsfe();
		}
		if (arana2_1.dmt2kp[i__ - 1] != 0.f)
		{
			s_wsfe(&io___390);
			r__1 = (i__ - .5f) * .05f;
			do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
			r__2 = scale2 * arana2_1.dmt2kp[i__ - 1];
			do_fio(&c__1, (char *)&r__2, (int)sizeof(float));
			e_wsfe();
		}
		if (arana2_1.dm2km[i__ - 1] != 0.f)
		{
			s_wsfe(&io___391);
			r__1 = (i__ - .5f) * .05f;
			do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
			r__2 = scale2 * .5f * arana2_1.dm2km[i__ - 1];
			do_fio(&c__1, (char *)&r__2, (int)sizeof(float));
			e_wsfe();
		}
		if (arana2_1.dm2k0s[i__ - 1] != 0.f)
		{
			s_wsfe(&io___392);
			r__1 = (i__ - .5f) * .05f;
			do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
			r__2 = scale2 * arana2_1.dm2k0s[i__ - 1];
			do_fio(&c__1, (char *)&r__2, (int)sizeof(float));
			e_wsfe();
		}
		if (arana2_1.dmt2la[i__ - 1] != 0.f)
		{
			s_wsfe(&io___393);
			r__1 = (i__ - .5f) * .05f;
			do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
			r__2 = scale2 * arana2_1.dmt2la[i__ - 1];
			do_fio(&c__1, (char *)&r__2, (int)sizeof(float));
			e_wsfe();
		}
		if (arana2_1.dmt2lb[i__ - 1] != 0.f)
		{
			s_wsfe(&io___394);
			r__1 = (i__ - .5f) * .05f;
			do_fio(&c__1, (char *)&r__1, (int)sizeof(float));
			r__2 = scale2 * arana2_1.dmt2lb[i__ - 1];
			do_fio(&c__1, (char *)&r__2, (int)sizeof(float));
			e_wsfe();
		}
	}
	return 0;
}

int hjana1_(void)
{

	static int iw = 0;

	int i__1, i__2;
	float r__1, r__2, r__3;

	int s_wsle(cilist *), do_lio(int *, int *, char *, int),
		e_wsle(void);
	double log(double), sqrt(double);

	static int i__, j;
	static float y1, y2, pe;
	static int ir;
	static float pm;
	static int iy;
	static float px, py, pz, yr, rap;
	static int imt;
	static float xmt, dyg1[50], dyp1[50];
	static int nisg;
	static float dxmt;
	static int ityp;
	static float dyg1c[50], deyg1[50], dmyg1[200], deyp1[50], seyg1[50], dmyp1[200], smyg1[200], seyp1[50], snyg1[50], smyp1[200], snyp1[50];
	static int nsubg, nisgs, nsubp, isevt, isrun;
	static float deyg1c[50], dmyg1c[50], dnrin1[50], seyg1c[50], dnrpj1[50],
		dnrtg1[50], smyg1c[50], snyg1c[50], snrin1[50], dnrtt1[50],
		snrpj1[50], snrtg1[50], snrtt1[50];
	static int nsubgs, nsubps;

	static cilist io___438 = {0, 6, 0, 0, 0};
	static cilist io___439 = {0, 6, 0, 0, 0};
	static cilist io___440 = {0, 6, 0, 0, 0};
	static cilist io___441 = {0, 6, 0, 0, 0};
	static cilist io___447 = {0, 6, 0, 0, 0};
	static cilist io___448 = {0, 6, 0, 0, 0};
	static cilist io___449 = {0, 6, 0, 0, 0};
	static cilist io___450 = {0, 6, 0, 0, 0};
	static cilist io___451 = {0, 6, 0, 0, 0};
	static cilist io___452 = {0, 6, 0, 0, 0};
	static cilist io___453 = {0, 6, 0, 0, 0};
	static cilist io___454 = {0, 6, 0, 0, 0};
	static cilist io___459 = {0, 6, 0, 0, 0};
	static cilist io___460 = {0, 6, 0, 0, 0};
	static cilist io___461 = {0, 6, 0, 0, 0};
	static cilist io___462 = {0, 6, 0, 0, 0};
	static cilist io___463 = {0, 6, 0, 0, 0};
	static cilist io___464 = {0, 6, 0, 0, 0};
	static cilist io___465 = {0, 6, 0, 0, 0};
	static cilist io___466 = {0, 6, 0, 0, 0};
	static cilist io___467 = {0, 6, 0, 0, 0};
	static cilist io___468 = {0, 6, 0, 0, 0};

	nsubg = 0;
	nisgs = 0;
	nsubgs = 0;
	nsubps = 0;
	if (isevt == arevt_1.iaevt && isrun == arevt_1.iarun)
	{
		for (i__ = 1; i__ <= 200; ++i__)
		{
			dmyp1[i__ - 1] = smyp1[i__ - 1];
			dmyg1[i__ - 1] = smyg1[i__ - 1];
		}
		for (i__ = 1; i__ <= 50; ++i__)
		{
			dyp1[i__ - 1] = snyp1[i__ - 1];
			deyp1[i__ - 1] = seyp1[i__ - 1];
			dyg1[i__ - 1] = snyg1[i__ - 1];
			deyg1[i__ - 1] = seyg1[i__ - 1];
			dnrpj1[i__ - 1] = snrpj1[i__ - 1];
			dnrtg1[i__ - 1] = snrtg1[i__ - 1];
			dnrin1[i__ - 1] = snrin1[i__ - 1];
			dnrtt1[i__ - 1] = snrtt1[i__ - 1];
			dyg1c[i__ - 1] = snyg1c[i__ - 1];
			dmyg1c[i__ - 1] = smyg1c[i__ - 1];
			deyg1c[i__ - 1] = seyg1c[i__ - 1];
		}
		nsubp = nsubps;
		nsubg = nsubgs;
		nisg = nisgs;
	}
	else
	{
		for (i__ = 1; i__ <= 200; ++i__)
		{
			smyp1[i__ - 1] = dmyp1[i__ - 1];
			smyg1[i__ - 1] = dmyg1[i__ - 1];
		}
		for (i__ = 1; i__ <= 50; ++i__)
		{
			snyp1[i__ - 1] = dyp1[i__ - 1];
			seyp1[i__ - 1] = deyp1[i__ - 1];
			snyg1[i__ - 1] = dyg1[i__ - 1];
			seyg1[i__ - 1] = deyg1[i__ - 1];
			snrpj1[i__ - 1] = dnrpj1[i__ - 1];
			snrtg1[i__ - 1] = dnrtg1[i__ - 1];
			snrin1[i__ - 1] = dnrin1[i__ - 1];
			snrtt1[i__ - 1] = dnrtt1[i__ - 1];
			snyg1c[i__ - 1] = dyg1c[i__ - 1];
			smyg1c[i__ - 1] = dmyg1c[i__ - 1];
			seyg1c[i__ - 1] = deyg1c[i__ - 1];
		}
		nsubps = nsubp;
		nsubgs = nsubg;
		nisgs = nisg;
		isevt = arevt_1.iaevt;
		isrun = arevt_1.iarun;
		++iw;
	}

	i__1 = hparnt_1.ihnt2[0];
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		i__2 = hjjet1_1.npj[i__ - 1];
		for (j = 1; j <= i__2; ++j)
		{
			ityp = hjjet1_1.kfpj[i__ + j * 300 - 301];
			px = hjjet1_1.pjpx[i__ + j * 300 - 301];
			py = hjjet1_1.pjpy[i__ + j * 300 - 301];
			pz = hjjet1_1.pjpz[i__ + j * 300 - 301];
			pe = hjjet1_1.pjpe[i__ + j * 300 - 301];
			pm = hjjet1_1.pjpm[i__ + j * 300 - 301];
			if (dabs(pz) >= pe)
			{
				s_wsle(&io___438);
				do_lio(&c__9, &c__1, " IN HJANA1, PROJ STR ", (int)21);
				do_lio(&c__3, &c__1, (char *)&i__, (int)sizeof(int));
				do_lio(&c__9, &c__1, " PART ", (int)6);
				do_lio(&c__3, &c__1, (char *)&j, (int)sizeof(int));
				e_wsle();
				s_wsle(&io___439);
				do_lio(&c__9, &c__1, " FLAV = ", (int)8);
				do_lio(&c__3, &c__1, (char *)&ityp, (int)sizeof(int));
				do_lio(&c__9, &c__1, " PX = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&px, (int)sizeof(float));
				do_lio(&c__9, &c__1, " PY = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&py, (int)sizeof(float));
				e_wsle();
				s_wsle(&io___440);
				do_lio(&c__9, &c__1, " PZ = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&pz, (int)sizeof(float));
				do_lio(&c__9, &c__1, " EE = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&pe, (int)sizeof(float));
				e_wsle();
				s_wsle(&io___441);
				do_lio(&c__9, &c__1, " XM = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&pm, (int)sizeof(float));
				e_wsle();
			}
			rap = log((pe + pz) / (pe - pz)) * .5f;

			r__1 = px;

			r__2 = py;

			r__3 = pm;
			xmt = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
			dxmt = xmt - pm;
			iy = (int)(dabs(rap) / .2f) + 1;
			if (iy > 50)
			{
                imt = (int)(dxmt / .05f) + 1;
                if (rap > 1.f || rap <= -1.f)
			{
				continue;
			}
			if (imt > 200)
			{
				continue;
			}
			dmyp1[imt - 1] += 1.f / xmt;
			if (ityp == 21)
			{
				dmyg1[imt - 1] += 1.f / xmt;
			}
			}
			dyp1[iy - 1] += 1.f;
			deyp1[iy - 1] += xmt;
			if (ityp == 21)
			{
				dyg1[iy - 1] += 1.f;
				deyg1[iy - 1] += xmt;
			}
		}
	}
	i__1 = hparnt_1.ihnt2[2];
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		i__2 = hjjet1_1.ntj[i__ - 1];
		for (j = 1; j <= i__2; ++j)
		{
			ityp = hjjet1_1.kftj[i__ + j * 300 - 301];
			px = hjjet1_1.pjtx[i__ + j * 300 - 301];
			py = hjjet1_1.pjty[i__ + j * 300 - 301];
			pz = hjjet1_1.pjtz[i__ + j * 300 - 301];
			pe = hjjet1_1.pjte[i__ + j * 300 - 301];
			pm = hjjet1_1.pjtm[i__ + j * 300 - 301];
			if (dabs(pz) >= pe)
			{
				s_wsle(&io___447);
				do_lio(&c__9, &c__1, " IN HJANA1, TARG STR ", (int)21);
				do_lio(&c__3, &c__1, (char *)&i__, (int)sizeof(int));
				do_lio(&c__9, &c__1, " PART ", (int)6);
				do_lio(&c__3, &c__1, (char *)&j, (int)sizeof(int));
				e_wsle();
				s_wsle(&io___448);
				do_lio(&c__9, &c__1, " FLAV = ", (int)8);
				do_lio(&c__3, &c__1, (char *)&ityp, (int)sizeof(int));
				do_lio(&c__9, &c__1, " PX = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&px, (int)sizeof(float));
				do_lio(&c__9, &c__1, " PY = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&py, (int)sizeof(float));
				e_wsle();
				s_wsle(&io___449);
				do_lio(&c__9, &c__1, " PZ = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&pz, (int)sizeof(float));
				do_lio(&c__9, &c__1, " EE = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&pe, (int)sizeof(float));
				e_wsle();
				s_wsle(&io___450);
				do_lio(&c__9, &c__1, " XM = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&pm, (int)sizeof(float));
				e_wsle();
			}
			rap = log((pe + pz) / (pe - pz)) * .5f;

			r__1 = px;

			r__2 = py;

			r__3 = pm;
			xmt = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
			dxmt = xmt - pm;
			iy = (int)(dabs(rap) / .2f) + 1;
			if (iy > 50)
			{
				if (rap > 1.f || rap <= -1.f)
			{
				continue;
			}
			imt = (int)(dxmt / .05f) + 1;
			if (imt > 200)
			{
				continue;
			}
			dmyp1[imt - 1] += 1.f / xmt;
			if (ityp == 21)
			{
				dmyg1[imt - 1] += 1.f / xmt;
			}
			}
			dyp1[iy - 1] += 1.f;
			deyp1[iy - 1] += xmt;
			if (ityp == 21)
			{
				dyg1[iy - 1] += 1.f;
				deyg1[iy - 1] += xmt;
			}
		}
	}
	i__1 = hjjet2_1.nsg;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		i__2 = hjjet2_1.njsg[i__ - 1];
		for (j = 1; j <= i__2; ++j)
		{
			ityp = hjjet2_1.k2sg[i__ + j * 150001 - 150002];
			px = hjjet2_1.pxsg[i__ + j * 150001 - 150002];
			py = hjjet2_1.pysg[i__ + j * 150001 - 150002];
			pz = hjjet2_1.pzsg[i__ + j * 150001 - 150002];
			pe = hjjet2_1.pesg[i__ + j * 150001 - 150002];
			pm = hjjet2_1.pmsg[i__ + j * 150001 - 150002];
			if (dabs(pz) >= pe)
			{
				s_wsle(&io___451);
				do_lio(&c__9, &c__1, " IN HJANA1, INDP STR ", (int)21);
				do_lio(&c__3, &c__1, (char *)&i__, (int)sizeof(int));
				do_lio(&c__9, &c__1, " PART ", (int)6);
				do_lio(&c__3, &c__1, (char *)&j, (int)sizeof(int));
				e_wsle();
				s_wsle(&io___452);
				do_lio(&c__9, &c__1, " FLAV = ", (int)8);
				do_lio(&c__3, &c__1, (char *)&ityp, (int)sizeof(int));
				do_lio(&c__9, &c__1, " PX = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&px, (int)sizeof(float));
				do_lio(&c__9, &c__1, " PY = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&py, (int)sizeof(float));
				e_wsle();
				s_wsle(&io___453);
				do_lio(&c__9, &c__1, " PZ = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&pz, (int)sizeof(float));
				do_lio(&c__9, &c__1, " EE = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&pe, (int)sizeof(float));
				e_wsle();
				s_wsle(&io___454);
				do_lio(&c__9, &c__1, " XM = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&pm, (int)sizeof(float));
				e_wsle();
			}
			rap = log((pe + pz) / (pe - pz)) * .5f;
			r__1 = px;
			r__2 = py;
			r__3 = pm;
			xmt = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
			dxmt = xmt - pm;
			iy = (int)(dabs(rap) / .2f) + 1;
			if (iy > 50)
			{
				if (rap > 1.f || rap <= -1.f)
			{
				continue;
			}
			imt = (int)(dxmt / .05f) + 1;
			if (imt > 200)
			{
				continue;
			}
			dmyp1[imt - 1] += 1.f / xmt;
			if (ityp == 21)
			{
				dmyg1[imt - 1] += 1.f / xmt;
			}
			}
			dyp1[iy - 1] += 1.f;
			deyp1[iy - 1] += xmt;
			if (ityp == 21)
			{
				dyg1[iy - 1] += 1.f;
				deyg1[iy - 1] += xmt;
			}
		}
	}
	i__1 = hparnt_1.ihnt2[0];
	for (i__ = 1; i__ <= i__1; ++i__)
	{

		r__1 = hjcrdn_1.yp[i__ * 3 - 3];

		r__2 = hjcrdn_1.yp[i__ * 3 - 2];
		yr = sqrt(r__1 * r__1 + r__2 * r__2);
		ir = (int)(yr / .2f) + 1;

		if (ir > 50 || ir < 1)
		{
			continue;
		}
		dnrpj1[ir - 1] += 1.f;
		dnrtt1[ir - 1] += 1.f;
	}
	i__1 = hparnt_1.ihnt2[2];
	for (i__ = 1; i__ <= i__1; ++i__)
	{

		r__1 = hjcrdn_1.yt[i__ * 3 - 3];

		r__2 = hjcrdn_1.yt[i__ * 3 - 2];
		yr = sqrt(r__1 * r__1 + r__2 * r__2);
		ir = (int)(yr / .2f) + 1;
		if (ir > 50 || ir < 1)
		{
			continue;
		}
		dnrtg1[ir - 1] += 1.f;
		dnrtt1[ir - 1] += 1.f;
	}
	i__1 = hjjet2_1.nsg;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		y1 = (hjcrdn_1.yp[hjjet2_1.iasg[i__ - 1] * 3 - 3] + hjcrdn_1.yt[hjjet2_1.iasg[i__ + 150000] * 3 - 3]) * .5f;
		y2 = (hjcrdn_1.yp[hjjet2_1.iasg[i__ - 1] * 3 - 2] + hjcrdn_1.yt[hjjet2_1.iasg[i__ + 150000] * 3 - 2]) * .5f;

		r__1 = y1;

		r__2 = y2;
		yr = sqrt(r__1 * r__1 + r__2 * r__2);
		ir = (int)(yr / .2f) + 1;
		if (ir > 50 || ir < 1)
		{
			continue;
		}
		dnrin1[ir - 1] += 1.f;
		dnrtt1[ir - 1] += 1.f;
	}
	i__1 = para1_1.mul;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		ityp = prec1_1.ityp0[i__ - 1];
		px = (float)prec1_1.px0[i__ - 1];
		py = (float)prec1_1.py0[i__ - 1];
		pz = (float)prec1_1.pz0[i__ - 1];
		pe = (float)prec1_1.e0[i__ - 1];
		pm = (float)prec1_1.xmass0[i__ - 1];
		if (dabs(pz) >= pe)
		{
			s_wsle(&io___459);
			do_lio(&c__9, &c__1, " IN HJANA1, GLUON ", (int)18);
			do_lio(&c__3, &c__1, (char *)&i__, (int)sizeof(int));
			e_wsle();
			s_wsle(&io___460);
			do_lio(&c__9, &c__1, " FLAV = ", (int)8);
			do_lio(&c__3, &c__1, (char *)&ityp, (int)sizeof(int));
			do_lio(&c__9, &c__1, " PX = ", (int)6);
			do_lio(&c__4, &c__1, (char *)&px, (int)sizeof(float));
			do_lio(&c__9, &c__1, " PY = ", (int)6);
			do_lio(&c__4, &c__1, (char *)&py, (int)sizeof(float));
			e_wsle();
			s_wsle(&io___461);
			do_lio(&c__9, &c__1, " PZ = ", (int)6);
			do_lio(&c__4, &c__1, (char *)&pz, (int)sizeof(float));
			do_lio(&c__9, &c__1, " EE = ", (int)6);
			do_lio(&c__4, &c__1, (char *)&pe, (int)sizeof(float));
			e_wsle();
			s_wsle(&io___462);
			do_lio(&c__9, &c__1, " XM = ", (int)6);
			do_lio(&c__4, &c__1, (char *)&pm, (int)sizeof(float));
			e_wsle();
		}
		rap = log((pe + pz) / (pe - pz)) * .5f;
		r__1 = px;
		r__2 = py;
		r__3 = pm;
		xmt = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
		dxmt = xmt - pm;
		iy = (int)(dabs(rap) / .2f) + 1;
		if (iy > 50)
		{
			if (rap > 1.f || rap <= -1.f)
		{
			continue;
		}
		imt = (int)(dxmt / .05f) + 1;
		if (imt > 50)
		{
			continue;
		}
        dmyg1c[imt - 1] += 1.f / xmt;
		}
		dyg1c[iy - 1] += 1.f;
		deyg1c[iy - 1] += xmt;
	}

	i__1 = hparnt_1.ihnt2[0];
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		i__2 = hjjet1_1.npj[i__ - 1];
		for (j = 1; j <= i__2; ++j)
		{
			++nsubp;
			if (hjjet1_1.kfpj[i__ + j * 300 - 301] == 21)
			{
				++nsubg;
			}
		}
	}
	i__1 = hparnt_1.ihnt2[2];
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		i__2 = hjjet1_1.ntj[i__ - 1];
		for (j = 1; j <= i__2; ++j)
		{
			++nsubp;
			if (hjjet1_1.kftj[i__ + j * 300 - 301] == 21)
			{
				++nsubg;
			}
		}
	}
	i__1 = hjjet2_1.nsg;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		i__2 = hjjet2_1.njsg[i__ - 1];
		for (j = 1; j <= i__2; ++j)
		{
			++nsubp;
			if (hjjet2_1.k2sg[i__ + j * 150001 - 150002] == 21)
			{
				++nsubg;
			}
		}
	}
	nisg += hjjet2_1.nsg;
	if (arout_1.iout == 1)
	{

		s_wsle(&io___463);
		do_lio(&c__9, &c__1, " in HJANA1 ", (int)11);
		e_wsle();
		s_wsle(&io___464);
		do_lio(&c__9, &c__1, " total number of partons = ", (int)27);
		i__1 = nsubp / iw;
		do_lio(&c__3, &c__1, (char *)&i__1, (int)sizeof(int));
		e_wsle();
		s_wsle(&io___465);
		do_lio(&c__9, &c__1, " total number of gluons = ", (int)26);
		i__1 = nsubg / iw;
		do_lio(&c__3, &c__1, (char *)&i__1, (int)sizeof(int));
		e_wsle();
		s_wsle(&io___466);
		do_lio(&c__9, &c__1, " number of projectile strings = ", (int)32);
		do_lio(&c__3, &c__1, (char *)&hparnt_1.ihnt2[0], (int)sizeof(int));
		e_wsle();
		s_wsle(&io___467);
		do_lio(&c__9, &c__1, " number of target strings = ", (int)28);
		do_lio(&c__3, &c__1, (char *)&hparnt_1.ihnt2[2], (int)sizeof(int));
		e_wsle();
		s_wsle(&io___468);
		do_lio(&c__9, &c__1, " number of independent strings = ", (int)33);
		i__1 = nisg / iw;
		do_lio(&c__3, &c__1, (char *)&i__1, (int)sizeof(int));
		e_wsle();
	}

	return 0;
}

int hjan1a_(void)
{

	static int iw = 0;

	int i__1;
	double d__1, d__2;

	double sqrt(double);

	static int i__, it, igx, igy;
	static float dtg1a[50], stg1a[50];
	static int isevt, isrun;
	extern int hjan1b_(void);
	static float dgxg1a[50], dgyg1a[50], sgxg1a[50], sgyg1a[50];

	if (isevt == arevt_1.iaevt && isrun == arevt_1.iarun)
	{
		for (i__ = 1; i__ <= 50; ++i__)
		{
			dgxg1a[i__ - 1] = sgxg1a[i__ - 1];
			dgyg1a[i__ - 1] = sgyg1a[i__ - 1];
			dtg1a[i__ - 1] = stg1a[i__ - 1];
		}
	}
	else
	{
		for (i__ = 1; i__ <= 50; ++i__)
		{
			sgxg1a[i__ - 1] = dgxg1a[i__ - 1];
			sgyg1a[i__ - 1] = dgyg1a[i__ - 1];
			stg1a[i__ - 1] = dtg1a[i__ - 1];
		}
		isevt = arevt_1.iaevt;
		isrun = arevt_1.iarun;
		++iw;
	}
	i__1 = para1_1.mul;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		igx = (int)((d__1 = prec2_1.gx5[i__ - 1], (float)abs(d__1)) /.2f) +1;
		if (igx > 50 || igx < 1)
		{
            igy = (int)((d__1 = prec2_1.gy5[i__ - 1], (float)abs(d__1)) /.2f) +1;
            if (igy > 50 || igy < 1)
            {
                d__1 = prec2_1.ft5[i__ - 1];
                d__2 = prec2_1.gz5[i__ - 1];
                it = (int)((float)sqrt(d__1 * d__1 - d__2 * d__2) / .2f) + 1;
                if (it > 50 || it < 1)
                {
                    continue;
                }
                    dtg1a[it - 1] += 1.f;
            }
                    dgyg1a[igy - 1] += 1.f;	
        }
                    dgxg1a[igx - 1] += 1.f;
    }
	hjan1b_();
	return 0;
}

int hjan1b_(void)
{

	static int iw = 0;
	int i__1;
	float r__1, r__2;
	double d__1, d__2;
	double sqrt(double);
	static int i__, j, k;
	static float r0;
	static int ir, it;
	static float gx0, gy0, tau7, dtg1b[50], stg1b[50];
	static int isevt, isrun;
	static float dnrg1b[50], snrg1b[50];

	if (isevt == arevt_1.iaevt && isrun == arevt_1.iarun)
	{
		for (i__ = 1; i__ <= 50; ++i__)
		{
			dnrg1b[i__ - 1] = snrg1b[i__ - 1];
			dtg1b[i__ - 1] = stg1b[i__ - 1];
		}
	}
	else
	{
		for (i__ = 1; i__ <= 50; ++i__)
		{
			snrg1b[i__ - 1] = dnrg1b[i__ - 1];
			stg1b[i__ - 1] = dtg1b[i__ - 1];
		}
		isevt = arevt_1.iaevt;
		isrun = arevt_1.iarun;
		++iw;
	}

	i__1 = para1_1.mul;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		j = ilist8_1.lstrg1[i__ - 1];
		if (j <= srec1_1.nsp)
		{
			k = j;
			gx0 = hjcrdn_1.yp[j * 3 - 3];
			gy0 = hjcrdn_1.yp[j * 3 - 2];
		}
		else if (j <= srec1_1.nsp + srec1_1.nst)
		{
			k = j - srec1_1.nsp;
			gx0 = hjcrdn_1.yt[k * 3 - 3];
			gy0 = hjcrdn_1.yt[k * 3 - 2];
		}
		else
		{
			k = j - srec1_1.nsp - srec1_1.nst;
			gx0 = (hjcrdn_1.yp[hjjet2_1.iasg[k - 1] * 3 - 3] + hjcrdn_1.yt[hjjet2_1.iasg[k + 150000] * 3 - 3]) * .5f;
			gy0 = (hjcrdn_1.yp[hjjet2_1.iasg[k - 1] * 3 - 2] + hjcrdn_1.yt[hjjet2_1.iasg[k + 150000] * 3 - 2]) * .5f;
		}

		r__1 = (float)prec2_1.gx5[i__ - 1] - gx0;

		r__2 = (float)prec2_1.gy5[i__ - 1] - gy0;
		r0 = sqrt(r__1 * r__1 + r__2 * r__2);
		ir = (int)(r0 / .2f) + 1;
		if (ir > 50 || ir < 1)
		{
			d__1 = prec2_1.ft5[i__ - 1];
            d__2 = prec2_1.gz5[i__ - 1];
            tau7 = sqrt((float)(d__1 * d__1 - d__2 * d__2));
            it = (int)(tau7 / .2f) + 1;
            if (it > 50 || it < 1)
            {
                continue
            }
            dtg1b[it - 1] += 1.f;
		}
		dnrg1b[ir - 1] += 1.f;	
	}

	return 0;
}

int hjana2_(void)
{
	static int iw = 0;
	int i__1, i__2;
	float r__1, r__2, r__3;
	double d__1, d__2;
	int s_wsle(cilist *), do_lio(int *, int *, char *, int),e_wsle(void);
	double log(double), sqrt(double);
	static int i__, j;
	static float pe;
	static int nj, ir;
	static float pm;
	static int it, iy;
	static float px, py, pz, yr, rap;
	static int imt;
	static float xmt, dyg2[50], dyp2[50];
	static int nisg;
	static float dxmt;
	static int ityp;
	static float dyg2c[50], deyg2[50], dtin2[50], dmyg2[200], deyp2[50], dtpj2[50], dttg2[50], seyg2[50], dmyp2[200], stin2[50], smyg2[200],seyp2[50], snyg2[50], stpj2[50], sttg2[50], smyp2[200], snyp2[50];
	static int nsubg, nisgs, nsubp, isevt, isrun;
	extern int hjan2a_(void), hjan2b_(void);
	static float deyg2c[50], dmyg2c[50], dnrin2[50], seyg2c[50], dnrpj2[50],dnrtg2[50], smyg2c[50], snyg2c[50], snrin2[50], dnrtt2[50],snrpj2[50], snrtg2[50], dttot2[50], snrtt2[50], sttot2[50];
	static int nsubgs, nsubps;
	static cilist io___549 = {0, 6, 0, 0, 0};
	static cilist io___550 = {0, 6, 0, 0, 0};
	static cilist io___551 = {0, 6, 0, 0, 0};
	static cilist io___552 = {0, 6, 0, 0, 0};
	static cilist io___558 = {0, 6, 0, 0, 0};
	static cilist io___559 = {0, 6, 0, 0, 0};
	static cilist io___560 = {0, 6, 0, 0, 0};
	static cilist io___561 = {0, 6, 0, 0, 0};
	static cilist io___563 = {0, 6, 0, 0, 0};
	static cilist io___564 = {0, 6, 0, 0, 0};
	static cilist io___565 = {0, 6, 0, 0, 0};
	static cilist io___566 = {0, 6, 0, 0, 0};
	static cilist io___570 = {0, 6, 0, 0, 0};
	static cilist io___571 = {0, 6, 0, 0, 0};
	static cilist io___572 = {0, 6, 0, 0, 0};
	static cilist io___573 = {0, 6, 0, 0, 0};
	static cilist io___574 = {0, 6, 0, 0, 0};
	static cilist io___575 = {0, 6, 0, 0, 0};
	static cilist io___576 = {0, 6, 0, 0, 0};
	static cilist io___577 = {0, 6, 0, 0, 0};
	static cilist io___578 = {0, 6, 0, 0, 0};
	static cilist io___579 = {0, 6, 0, 0, 0};

	if (isevt == arevt_1.iaevt && isrun == arevt_1.iarun)
	{
		for (i__ = 1; i__ <= 200; ++i__)
		{
			dmyp2[i__ - 1] = smyp2[i__ - 1];
			dmyg2[i__ - 1] = smyg2[i__ - 1];
		}
		for (i__ = 1; i__ <= 50; ++i__)
		{
			dyp2[i__ - 1] = snyp2[i__ - 1];
			deyp2[i__ - 1] = seyp2[i__ - 1];
			dyg2[i__ - 1] = snyg2[i__ - 1];
			deyg2[i__ - 1] = seyg2[i__ - 1];
			dnrpj2[i__ - 1] = snrpj2[i__ - 1];
			dnrtg2[i__ - 1] = snrtg2[i__ - 1];
			dnrin2[i__ - 1] = snrin2[i__ - 1];
			dnrtt2[i__ - 1] = snrtt2[i__ - 1];
			dtpj2[i__ - 1] = stpj2[i__ - 1];
			dttg2[i__ - 1] = sttg2[i__ - 1];
			dtin2[i__ - 1] = stin2[i__ - 1];
			dttot2[i__ - 1] = sttot2[i__ - 1];
			dyg2c[i__ - 1] = snyg2c[i__ - 1];
			dmyg2c[i__ - 1] = smyg2c[i__ - 1];
			deyg2c[i__ - 1] = seyg2c[i__ - 1];
		}
		nsubp = nsubps;
		nsubg = nsubgs;
		nisg = nisgs;
	}
	else
	{
		for (i__ = 1; i__ <= 200; ++i__)
		{
			smyp2[i__ - 1] = dmyp2[i__ - 1];
			smyg2[i__ - 1] = dmyg2[i__ - 1];
		}
		for (i__ = 1; i__ <= 50; ++i__)
		{
			snyp2[i__ - 1] = dyp2[i__ - 1];
			seyp2[i__ - 1] = deyp2[i__ - 1];
			snyg2[i__ - 1] = dyg2[i__ - 1];
			seyg2[i__ - 1] = deyg2[i__ - 1];
			snrpj2[i__ - 1] = dnrpj2[i__ - 1];
			snrtg2[i__ - 1] = dnrtg2[i__ - 1];
			snrin2[i__ - 1] = dnrin2[i__ - 1];
			snrtt2[i__ - 1] = dnrtt2[i__ - 1];
			stpj2[i__ - 1] = dtpj2[i__ - 1];
			sttg2[i__ - 1] = dttg2[i__ - 1];
			stin2[i__ - 1] = dtin2[i__ - 1];
			sttot2[i__ - 1] = dttot2[i__ - 1];
			snyg2c[i__ - 1] = dyg2c[i__ - 1];
			smyg2c[i__ - 1] = dmyg2c[i__ - 1];
			seyg2c[i__ - 1] = deyg2c[i__ - 1];
		}
		nsubps = nsubp;
		nsubgs = nsubg;
		nisgs = nisg;
		isevt = arevt_1.iaevt;
		isrun = arevt_1.iarun;
		++iw;
	}

	if (anim_1.isoft == 3 || anim_1.isoft == 4 || anim_1.isoft == 5)
	{
		goto L510;
	}

	i__1 = hparnt_1.ihnt2[0];
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		i__2 = hjjet1_1.npj[i__ - 1];
		for (j = 1; j <= i__2; ++j)
		{
			ityp = hjjet1_1.kfpj[i__ + j * 300 - 301];
			px = hjjet1_1.pjpx[i__ + j * 300 - 301];
			py = hjjet1_1.pjpy[i__ + j * 300 - 301];
			pz = hjjet1_1.pjpz[i__ + j * 300 - 301];
			pe = hjjet1_1.pjpe[i__ + j * 300 - 301];
			pm = hjjet1_1.pjpm[i__ + j * 300 - 301];

			if (dabs(pz) >= pe)
			{
				s_wsle(&io___549);
				do_lio(&c__9, &c__1, " IN HJANA2, PROJ STR ", (int)21);
				do_lio(&c__3, &c__1, (char *)&i__, (int)sizeof(int));
				do_lio(&c__9, &c__1, " PART ", (int)6);
				do_lio(&c__3, &c__1, (char *)&j, (int)sizeof(int));
				e_wsle();
				s_wsle(&io___550);
				do_lio(&c__9, &c__1, " FLAV = ", (int)8);
				do_lio(&c__3, &c__1, (char *)&ityp, (int)sizeof(int));
				do_lio(&c__9, &c__1, " PX = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&px, (int)sizeof(float));
				do_lio(&c__9, &c__1, " PY = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&py, (int)sizeof(float));
				e_wsle();
				s_wsle(&io___551);
				do_lio(&c__9, &c__1, " PZ = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&pz, (int)sizeof(float));
				do_lio(&c__9, &c__1, " EE = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&pe, (int)sizeof(float));
				e_wsle();
				s_wsle(&io___552);
				do_lio(&c__9, &c__1, " XM = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&pm, (int)sizeof(float));
				e_wsle();
			}
			rap = log((pe + pz) / (pe - pz)) * .5f;
			r__1 = px;
			r__2 = py;
			r__3 = pm;
			xmt = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
			dxmt = xmt - pm;
			iy = (int)(dabs(rap) / .2f) + 1;
			if (iy > 50)
			{
				if (rap > 1.f || rap <= -1.f)
			{
				continue;
			}
            imt = (int)(dxmt / .05f) + 1;
			if (imt > 200)
			{
				continue;
			}
			dmyp2[imt - 1] += 1.f / xmt;
			if (ityp == 21)
            {
				dmyg2[imt - 1] += 1.f / xmt;
                }
            }
			dyp2[iy - 1] += 1.f;
			deyp2[iy - 1] += xmt;
			if (ityp == 21)
			{
				dyg2[iy - 1] += 1.f;
				deyg2[iy - 1] += xmt;
			}
			
		}
	}
	i__1 = hparnt_1.ihnt2[2];
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		i__2 = hjjet1_1.ntj[i__ - 1];
		for (j = 1; j <= i__2; ++j)
		{
			ityp = hjjet1_1.kftj[i__ + j * 300 - 301];
			px = hjjet1_1.pjtx[i__ + j * 300 - 301];
			py = hjjet1_1.pjty[i__ + j * 300 - 301];
			pz = hjjet1_1.pjtz[i__ + j * 300 - 301];
			pe = hjjet1_1.pjte[i__ + j * 300 - 301];
			pm = hjjet1_1.pjtm[i__ + j * 300 - 301];

			if (dabs(pz) >= pe)
			{
				s_wsle(&io___558);
				do_lio(&c__9, &c__1, " IN HJANA2, TARG STR ", (int)21);
				do_lio(&c__3, &c__1, (char *)&i__, (int)sizeof(int));
				do_lio(&c__9, &c__1, " PART ", (int)6);
				do_lio(&c__3, &c__1, (char *)&j, (int)sizeof(int));
				e_wsle();
				s_wsle(&io___559);
				do_lio(&c__9, &c__1, " FLAV = ", (int)8);
				do_lio(&c__3, &c__1, (char *)&ityp, (int)sizeof(int));
				do_lio(&c__9, &c__1, " PX = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&px, (int)sizeof(float));
				do_lio(&c__9, &c__1, " PY = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&py, (int)sizeof(float));
				e_wsle();
				s_wsle(&io___560);
				do_lio(&c__9, &c__1, " PZ = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&pz, (int)sizeof(float));
				do_lio(&c__9, &c__1, " EE = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&pe, (int)sizeof(float));
				e_wsle();
				s_wsle(&io___561);
				do_lio(&c__9, &c__1, " XM = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&pm, (int)sizeof(float));
				e_wsle();
			}
			rap = log((pe + pz) / (pe - pz)) * .5f;
			r__1 = px;
			r__2 = py;
			r__3 = pm;
			xmt = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
			dxmt = xmt - pm;
			iy = (int)(dabs(rap) / .2f) + 1;
			if (iy > 50)
			{
				if (rap > 1.f || rap <= -1.f)
			{
				continue;
			}
			imt = (int)(dxmt / .05f) + 1;
			if (imt > 200)
			{
				continue;
			}
			dmyp2[imt - 1] += 1.f / xmt;
			if (ityp == 21)
            {
				dmyg2[imt - 1] += 1.f / xmt;
			}
			}
			dyp2[iy - 1] += 1.f;
			deyp2[iy - 1] += xmt;
			if (ityp == 21)
			{
				dyg2[iy - 1] += 1.f;
				deyg2[iy - 1] += xmt;
			}		
		}
	}

L510:
	i__1 = hjjet2_1.nsg;
	for (i__ = 1; i__ <= i__1; ++i__)
	{

		nj = hjjet2_1.njsg[i__ - 1];
		if (anim_1.isoft == 3 || anim_1.isoft == 4 || anim_1.isoft == 5)
		{
			nj = soft_1.njsgs[i__ - 1];
		}
		i__2 = nj;
		for (j = 1; j <= i__2; ++j)
		{

			ityp = hjjet2_1.k2sg[i__ + j * 150001 - 150002];
			px = hjjet2_1.pxsg[i__ + j * 150001 - 150002];
			py = hjjet2_1.pysg[i__ + j * 150001 - 150002];
			pz = hjjet2_1.pzsg[i__ + j * 150001 - 150002];
			pe = hjjet2_1.pesg[i__ + j * 150001 - 150002];
			pm = hjjet2_1.pmsg[i__ + j * 150001 - 150002];

			if (anim_1.isoft == 3 || anim_1.isoft == 4 || anim_1.isoft == 5)
			{
				ityp = soft_1.k2sgs[i__ + j * 150001 - 150002];
				px = (float)soft_1.pxsgs[i__ + j * 150001 - 150002];
				py = (float)soft_1.pysgs[i__ + j * 150001 - 150002];
				pz = (float)soft_1.pzsgs[i__ + j * 150001 - 150002];
				pe = (float)soft_1.pesgs[i__ + j * 150001 - 150002];
				pm = (float)soft_1.pmsgs[i__ + j * 150001 - 150002];
			}

			if (dabs(pz) >= pe)
			{
				s_wsle(&io___563);
				do_lio(&c__9, &c__1, " IN HJANA2, INDP STR ", (int)21);
				do_lio(&c__3, &c__1, (char *)&i__, (int)sizeof(int));
				do_lio(&c__9, &c__1, " PART ", (int)6);
				do_lio(&c__3, &c__1, (char *)&j, (int)sizeof(int));
				e_wsle();
				s_wsle(&io___564);
				do_lio(&c__9, &c__1, " FLAV = ", (int)8);
				do_lio(&c__3, &c__1, (char *)&ityp, (int)sizeof(int));
				do_lio(&c__9, &c__1, " PX = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&px, (int)sizeof(float));
				do_lio(&c__9, &c__1, " PY = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&py, (int)sizeof(float));
				e_wsle();
				s_wsle(&io___565);
				do_lio(&c__9, &c__1, " PZ = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&pz, (int)sizeof(float));
				do_lio(&c__9, &c__1, " EE = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&pe, (int)sizeof(float));
				e_wsle();
				s_wsle(&io___566);
				do_lio(&c__9, &c__1, " XM = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&pm, (int)sizeof(float));
				e_wsle();
			}
			rap = log((pe + pz) / (pe - pz)) * .5f;
			r__1 = px;
			r__2 = py;
			r__3 = pm;
			xmt = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
			dxmt = xmt - pm;
			iy = (int)(dabs(rap) / .2f) + 1;
			if (iy > 50)
			{
				if (rap > 1.f || rap <= -1.f)
			{
				continue;
			}
			imt = (int)(dxmt / .05f) + 1;
			if (imt > 200)
			{
				continue;
			}
			dmyp2[imt - 1] += 1.f / xmt;
			if (ityp == 21)
			{
				dmyg2[imt - 1] += 1.f / xmt;
			}
			}
			dyp2[iy - 1] += 1.f;
			deyp2[iy - 1] += xmt;
			if (ityp == 21)
			{
				dyg2[iy - 1] += 1.f;
				deyg2[iy - 1] += xmt;
			}
		}
	}

	if (anim_1.isoft == 3 || anim_1.isoft == 4 || anim_1.isoft == 5)
	{
		goto L520;
	}
	i__1 = hparnt_1.ihnt2[0];
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		j = i__;
		d__1 = srec2_1.zt1[j - 1];
		d__2 = srec2_1.zt2[j - 1];
		yr = sqrt((float)(d__1 * d__1 + d__2 * d__2));
		ir = (int)(yr / .2f) + 1;
		if (ir > 50 || ir < 1)
		{
            it = (int)((float)srec2_1.ataui[j - 1] / .2f) + 1;
            if (it > 50 || it < 1)
            {
                continue;
                }
                dtpj2[it - 1] += 1.f;
                dttot2[it - 1] += 1.f;
                }
                dnrpj2[ir - 1] += 1.f;
                dnrtt2[ir - 1] += 1.f;
	}
	i__1 = hparnt_1.ihnt2[2];
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		j = i__ + hparnt_1.ihnt2[0];
		d__1 = srec2_1.zt1[j - 1];
		d__2 = srec2_1.zt2[j - 1];
		yr = sqrt((float)(d__1 * d__1 + d__2 * d__2));
		ir = (int)(yr / .2f) + 1;
		if (ir > 50 || ir < 1)
		{
            it = (int)((float)srec2_1.ataui[j - 1] / .2f) + 1;
            if (it > 50 || it < 1)
            {
                continue;
                }
                dttg2[it - 1] += 1.f;
                dttot2[it - 1] += 1.f;
                }
                dnrtg2[ir - 1] += 1.f;
                dnrtt2[ir - 1] += 1.f;
	}

L520:
	i__1 = hjjet2_1.nsg;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		j = i__ + hparnt_1.ihnt2[0] + hparnt_1.ihnt2[2];
		if (anim_1.isoft == 3 || anim_1.isoft == 4 || anim_1.isoft == 5)
		{
			j = i__;
		}
		d__1 = srec2_1.zt1[j - 1];
		d__2 = srec2_1.zt2[j - 1];
		yr = sqrt((float)(d__1 * d__1 + d__2 * d__2));
		ir = (int)(yr / .2f) + 1;
		if (ir > 50 || ir < 1)
		{
            it = (int)((float)srec2_1.ataui[j - 1] / .2f) + 1;
            if (it > 50 || it < 1)
            {
                continue;
                }
                dtin2[it - 1] += 1.f;
                dttot2[it - 1] += 1.f;
                }
                dnrin2[ir - 1] += 1.f;
                dnrtt2[ir - 1] += 1.f;
	}
	i__1 = para1_1.mul;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		ityp = prec2_1.ityp5[i__ - 1];
		px = (float)prec2_1.px5[i__ - 1];
		py = (float)prec2_1.py5[i__ - 1];
		pz = (float)prec2_1.pz5[i__ - 1];
		pe = (float)prec2_1.e5[i__ - 1];
		pm = (float)prec2_1.xmass5[i__ - 1];

		if (dabs(pz) >= pe)
		{
			s_wsle(&io___570);
			do_lio(&c__9, &c__1, " IN HJANA2, GLUON ", (int)18);
			do_lio(&c__3, &c__1, (char *)&i__, (int)sizeof(int));
			e_wsle();
			s_wsle(&io___571);
			do_lio(&c__9, &c__1, " FLAV = ", (int)8);
			do_lio(&c__3, &c__1, (char *)&ityp, (int)sizeof(int));
			do_lio(&c__9, &c__1, " PX = ", (int)6);
			do_lio(&c__4, &c__1, (char *)&px, (int)sizeof(float));
			do_lio(&c__9, &c__1, " PY = ", (int)6);
			do_lio(&c__4, &c__1, (char *)&py, (int)sizeof(float));
			e_wsle();
			s_wsle(&io___572);
			do_lio(&c__9, &c__1, " PZ = ", (int)6);
			do_lio(&c__4, &c__1, (char *)&pz, (int)sizeof(float));
			do_lio(&c__9, &c__1, " EE = ", (int)6);
			do_lio(&c__4, &c__1, (char *)&pe, (int)sizeof(float));
			e_wsle();
			s_wsle(&io___573);
			do_lio(&c__9, &c__1, " XM = ", (int)6);
			do_lio(&c__4, &c__1, (char *)&pm, (int)sizeof(float));
			e_wsle();
		}
		rap = log((pe + pz) / (pe - pz)) * .5f;
		r__1 = px;
		r__2 = py;
		r__3 = pm;
		xmt = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
		dxmt = xmt - pm;
		iy = (int)(dabs(rap) / .2f) + 1;
		if (iy > 50)
		{
            if (rap > 1.f || rap <= -1.f)
            {
                continue;
            }
                imt = (int)(dxmt / .05f) + 1;
                if (imt > 50)
                {
                    continue;
                }
            dmyg2c[imt - 1] += 1.f / xmt;
        }
        dyg2c[iy - 1] += 1.f;
        deyg2c[iy - 1] += xmt;
	}

	if (anim_1.isoft == 3 || anim_1.isoft == 4 || anim_1.isoft == 5)
	{
		goto L530;
	}

	i__1 = hparnt_1.ihnt2[0];
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		i__2 = hjjet1_1.npj[i__ - 1];
		for (j = 1; j <= i__2; ++j)
		{
			++nsubp;
			if (hjjet1_1.kfpj[i__ + j * 300 - 301] == 21)
			{
				++nsubg;
			}
		}
	}
	i__1 = hparnt_1.ihnt2[2];
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		i__2 = hjjet1_1.ntj[i__ - 1];
		for (j = 1; j <= i__2; ++j)
		{
			++nsubp;
			if (hjjet1_1.kftj[i__ + j * 300 - 301] == 21)
			{
				++nsubg;
			}
		}
	}

L530:
	i__1 = hjjet2_1.nsg;
	for (i__ = 1; i__ <= i__1; ++i__)
	{

		nj = hjjet2_1.njsg[i__ - 1];
		if (anim_1.isoft == 3 || anim_1.isoft == 4 || anim_1.isoft == 5)
		{
			nj = soft_1.njsgs[i__ - 1];
		}
		i__2 = nj;
		for (j = 1; j <= i__2; ++j)
		{

			++nsubp;

			if (anim_1.isoft == 3 || anim_1.isoft == 4 || anim_1.isoft == 5)
			{
				if (soft_1.k2sgs[i__ + j * 150001 - 150002] == 21)
				{
					++nsubg;
				}
			}
			else
			{
				if (hjjet2_1.k2sg[i__ + j * 150001 - 150002] == 21)
				{
					++nsubg;
				}
			}
		}
	}

	nisg += hjjet2_1.nsg;
	if (arout_1.iout == 1)
	{

		s_wsle(&io___574);
		do_lio(&c__9, &c__1, " in HJANA2 ", (int)11);
		e_wsle();
		s_wsle(&io___575);
		do_lio(&c__9, &c__1, " total number of partons = ", (int)27);
		i__1 = nsubp / iw;
		do_lio(&c__3, &c__1, (char *)&i__1, (int)sizeof(int));
		e_wsle();
		s_wsle(&io___576);
		do_lio(&c__9, &c__1, " total number of gluons = ", (int)26);
		i__1 = nsubg / iw;
		do_lio(&c__3, &c__1, (char *)&i__1, (int)sizeof(int));
		e_wsle();
		s_wsle(&io___577);
		do_lio(&c__9, &c__1, " number of projectile strings = ", (int)32);
		do_lio(&c__3, &c__1, (char *)&hparnt_1.ihnt2[0], (int)sizeof(int));
		e_wsle();
		s_wsle(&io___578);
		do_lio(&c__9, &c__1, " number of target strings = ", (int)28);
		do_lio(&c__3, &c__1, (char *)&hparnt_1.ihnt2[2], (int)sizeof(int));
		e_wsle();
		s_wsle(&io___579);
		do_lio(&c__9, &c__1, " number of independent strings = ", (int)33);
		i__1 = nisg / iw;
		do_lio(&c__3, &c__1, (char *)&i__1, (int)sizeof(int));
		e_wsle();
	}
	hjan2a_();
	hjan2b_();
	return 0;
}

int hjan2a_(void)
{

	static int iw = 0;
	int i__1, i__2;
	float r__1;
	double d__1, d__2;
	double sqrt(double);
	static int i__, j, it, igx, igy;
	static float dtg2a[50], dtp2a[50], stg2a[50], stp2a[50];
	static int isevt, isrun;
	static float dgxg2a[50], dgyg2a[50], dgxp2a[50], dgyp2a[50], sgxg2a[50],sgyg2a[50], sgxp2a[50], sgyp2a[50];

	if (isevt == arevt_1.iaevt && isrun == arevt_1.iarun)
	{
		for (i__ = 1; i__ <= 50; ++i__)
		{
			dgxp2a[i__ - 1] = sgxp2a[i__ - 1];
			dgyp2a[i__ - 1] = sgyp2a[i__ - 1];
			dtp2a[i__ - 1] = stp2a[i__ - 1];
			dgxg2a[i__ - 1] = sgxg2a[i__ - 1];
			dgyg2a[i__ - 1] = sgyg2a[i__ - 1];
			dtg2a[i__ - 1] = stg2a[i__ - 1];
		}
	}
	else
	{
		for (i__ = 1; i__ <= 50; ++i__)
		{
			sgxp2a[i__ - 1] = dgxp2a[i__ - 1];
			sgyp2a[i__ - 1] = dgyp2a[i__ - 1];
			stp2a[i__ - 1] = dtp2a[i__ - 1];
			sgxg2a[i__ - 1] = dgxg2a[i__ - 1];
			sgyg2a[i__ - 1] = dgyg2a[i__ - 1];
			stg2a[i__ - 1] = dtg2a[i__ - 1];
		}
		isevt = arevt_1.iaevt;
		isrun = arevt_1.iarun;
		++iw;
	}

	i__1 = hparnt_1.ihnt2[0];
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		i__2 = hjjet1_1.npj[i__ - 1];
		for (j = 1; j <= i__2; ++j)
		{
			if (hjjet1_1.kfpj[i__ + j * 300 - 301] != 21)
			{
				igx = (int)((r__1 = hjcrdn_1.yp[i__ * 3 - 3], dabs(r__1)) / .2f) + 1;
				if (igx > 50 || igx < 1)
				{
					igy = (int)((r__1 = hjcrdn_1.yp[i__ * 3 - 2], dabs(r__1)) / .2f) + 1;
                    if (igy > 50 || igy < 1)
                    {
                        it = 1;
                        dtp2a[it - 1] += 1.f;
                    }
				dgyp2a[igy - 1] += 1.f;
				}
				dgxp2a[igx - 1] += 1.f;
			}
		}
	}
	i__1 = hparnt_1.ihnt2[2];
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		i__2 = hjjet1_1.ntj[i__ - 1];
		for (j = 1; j <= i__2; ++j)
		{
			if (hjjet1_1.kftj[i__ + j * 300 - 301] != 21)
			{
				igx = (int)((r__1 = hjcrdn_1.yt[i__ * 3 - 3], dabs(r__1)) / .2f) + 1;
				if (igx > 50 || igx < 1)
				{
                    igy = (int)((r__1 = hjcrdn_1.yt[i__ * 3 - 2], dabs(r__1)) / .2f) + 1;
                    if (igy > 50 || igy < 1)
                    {
                        it = 1;
                        dtp2a[it - 1] += 1.f;
				}
				dgyp2a[igy - 1] += 1.f;
				}
				dgxp2a[igx - 1] += 1.f;
			}
		}
	}
	i__1 = hjjet2_1.nsg;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		i__2 = hjjet2_1.njsg[i__ - 1];
		for (j = 1; j <= i__2; ++j)
		{
			if (hjjet2_1.k2sg[i__ + j * 150001 - 150002] != 21)
			{
				igx = (int)((r__1 = (hjcrdn_1.yp[hjjet2_1.iasg[i__ - 1] *3 -3] +hjcrdn_1.yt[hjjet2_1.iasg[i__ + 150000] * 3 - 3]) *.5f,dabs(r__1)) /.2f) +1;
				if (igx > 50 || igx < 1)
				{
					igy = (int)((r__1 = (hjcrdn_1.yp[hjjet2_1.iasg[i__ - 1] *3 -2] +hjcrdn_1.yt[hjjet2_1.iasg[i__ + 150000] * 3 - 2]) *.5f,dabs(r__1)) /.2f) +1;
                    if (igy > 50 || igy < 1)
				{
                    it = 1;
                    dtp2a[it - 1] += 1.f;
				}
				dgyp2a[igy - 1] += 1.f;
				}
				dgxp2a[igx - 1] += 1.f;	
			}
		}
	}
	i__1 = para1_1.mul;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		igx = (int)((r__1 = (float)prec2_1.gx5[i__ - 1], dabs(r__1)) /.2f) +1;
		if (igx > 50 || igx < 1)
		{
			igy = (int)((r__1 = (float)prec2_1.gy5[i__ - 1], dabs(r__1)) /.2f) +1;
		if (igy > 50 || igy < 1)
		{
            d__1 = prec2_1.ft5[i__ - 1];
            d__2 = prec2_1.gz5[i__ - 1];
            it = (int)(sqrt((float)(d__1 * d__1 - d__2 * d__2)) / .2f) + 1;
            if (it > 50 || it < 1)
            {
			continue;
            }
            dtg2a[it - 1] += 1.f;
            dtp2a[it - 1] += 1.f;
		}
		dgyg2a[igy - 1] += 1.f;
		dgyp2a[igy - 1] += 1.f;		
		}
		dgxg2a[igx - 1] += 1.f;
		dgxp2a[igx - 1] += 1.f;
	}

	return 0;
}

int hjan2b_(void)
{

	static int iw = 0;
	int i__1;
	float r__1, r__2;
	double d__1, d__2;
	double sqrt(double);
	static int i__, j;
	static float r0;
	static int ir, it;
	static float gx0, gy0, tau7, dtau, dtg2b[50], stg2b[50];
	static int isevt, isrun;
	static float dnrg2b[50], snrg2b[50];

	if (isevt == arevt_1.iaevt && isrun == arevt_1.iarun)
	{
		for (i__ = 1; i__ <= 50; ++i__)
		{
			dnrg2b[i__ - 1] = snrg2b[i__ - 1];
			dtg2b[i__ - 1] = stg2b[i__ - 1];
		}
	}
	else
	{
		for (i__ = 1; i__ <= 50; ++i__)
		{
			snrg2b[i__ - 1] = dnrg2b[i__ - 1];
			stg2b[i__ - 1] = dtg2b[i__ - 1];
		}
		isevt = arevt_1.iaevt;
		isrun = arevt_1.iarun;
		++iw;
	}

	i__1 = para1_1.mul;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		j = ilist8_1.lstrg1[i__ - 1];
		gx0 = (float)srec2_1.zt1[j - 1];
		gy0 = (float)srec2_1.zt2[j - 1];

		r__1 = (float)prec2_1.gx5[i__ - 1] - gx0;

		r__2 = (float)prec2_1.gy5[i__ - 1] - gy0;
		r0 = sqrt(r__1 * r__1 + r__2 * r__2);
		ir = (int)(r0 / .2f) + 1;
		if (ir > 50 || ir < 1)
		{
            d__1 = prec2_1.ft5[i__ - 1];
            d__2 = prec2_1.gz5[i__ - 1];
            tau7 = sqrt((float)(d__1 * d__1 - d__2 * d__2));
            dtau = tau7 - (float)srec2_1.ataui[j - 1];
            it = (int)(dtau / .2f) + 1;
            if (it > 25 || it < -24)
            {
                continue;
                }
                dtg2b[it + 24] += 1.f;
                }
                dnrg2b[ir - 1] += 1.f;
	}

	return 0;
}

int hjana3_(void)
{
	static int iw = 0;
	int i__1, i__2;
	float r__1, r__2, r__3;
	double sqrt(double);
	int s_wsle(cilist *), do_lio(int *, int *, char *, int),e_wsle(void);
	double log(double);
	static int i__, j;
	static float y, ee;
	static int iy;
	static float xm, px, py, pz;
	static int imt;
	static float xmt, dxmt;
	static int ityp;
	static float deyh3[50], dmyh3[50], dndyh3[50];
	static cilist io___626 = {0, 6, 0, 0, 0};
	static cilist io___627 = {0, 6, 0, 0, 0};
	static cilist io___628 = {0, 6, 0, 0, 0};
	static cilist io___629 = {0, 6, 0, 0, 0};
	static cilist io___630 = {0, 6, 0, 0, 0};
	++iw;
	i__1 = run_1.num;
	for (j = 1; j <= i__1; ++j)
	{
		i__2 = arerc1_1.multi1[j - 1];
		for (i__ = 1; i__ <= i__2; ++i__)
		{
			ityp = arprc1_1.ityp1[i__ + j * 150001 - 150002];
			if (ityp > -100 && ityp < 100)
			{
				continue;
			}
			px = arprc1_1.px1[i__ + j * 150001 - 150002];
			py = arprc1_1.py1[i__ + j * 150001 - 150002];
			pz = arprc1_1.pz1[i__ + j * 150001 - 150002];
			ee = arprc1_1.ee1[i__ + j * 150001 - 150002];
			xm = arprc1_1.xm1[i__ + j * 150001 - 150002];
			r__1 = px;
			r__2 = py;
			r__3 = xm;
			xmt = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
			if (dabs(pz) >= ee)
			{
				s_wsle(&io___626);
				do_lio(&c__9, &c__1, "IN HJANA3", (int)9);
				e_wsle();
				s_wsle(&io___627);
				do_lio(&c__9, &c__1, " PARTICLE ", (int)10);
				do_lio(&c__3, &c__1, (char *)&i__, (int)sizeof(int));
				do_lio(&c__9, &c__1, " RUN ", (int)5);
				do_lio(&c__3, &c__1, (char *)&j, (int)sizeof(int));
				do_lio(&c__9, &c__1, "PREC ERR", (int)8);
				e_wsle();
				s_wsle(&io___628);
				do_lio(&c__9, &c__1, " FLAV = ", (int)8);
				do_lio(&c__3, &c__1, (char *)&ityp, (int)sizeof(int));
				do_lio(&c__9, &c__1, " PX = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&px, (int)sizeof(float));
				do_lio(&c__9, &c__1, " PY = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&py, (int)sizeof(float));
				e_wsle();
				s_wsle(&io___629);
				do_lio(&c__9, &c__1, " PZ = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&pz, (int)sizeof(float));
				do_lio(&c__9, &c__1, " EE = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&ee, (int)sizeof(float));
				e_wsle();
				s_wsle(&io___630);
				do_lio(&c__9, &c__1, " XM = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&xm, (int)sizeof(float));
				e_wsle();
			}
			dxmt = xmt - xm;
			y = log((ee + pz) / (ee - pz)) * .5f;

			iy = (int)((y + 10.f) / .2f) + 1;
			if (iy > 50)
			{
				if (y < -1.f || y >= 1.f)
			{
				continue;
			}
			imt = (int)(dxmt / .05f) + 1;
			if (imt > 50)
			{
				continue;
			}
			dmyh3[imt - 1] += 1.f / xmt;
			}
			dndyh3[iy - 1] += 1.f;
			deyh3[iy - 1] += xmt;
		}
	}

	return 0;
}

int hjana4_(void)
{
	static int iw = 0;
	int i__1, i__2;
	float r__1, r__2, r__3;
	double sqrt(double);
	int s_wsle(cilist *), do_lio(int *, int *, char *, int),e_wsle(void);
	double log(double);
	static int i__, j;
	static float y, ee;
	static int iy;
	static float xm, px, py, pz;
	static int imt;
	static float xmt, dxmt;
	static int ityp;
	static float deyh4[50], dmyh4[50], dndyh4[50];
	static cilist io___648 = {0, 6, 0, 0, 0};
	static cilist io___649 = {0, 6, 0, 0, 0};
	static cilist io___650 = {0, 6, 0, 0, 0};
	static cilist io___651 = {0, 6, 0, 0, 0};
	static cilist io___652 = {0, 6, 0, 0, 0};
	++iw;
	i__1 = run_1.num;
	for (j = 1; j <= i__1; ++j)
	{
		i__2 = arerc1_1.multi1[j - 1];
		for (i__ = 1; i__ <= i__2; ++i__)
		{
			ityp = arprc1_1.ityp1[i__ + j * 150001 - 150002];
			if (ityp > -100 && ityp < 100)
			{
				continue;
			}
			px = arprc1_1.px1[i__ + j * 150001 - 150002];
			py = arprc1_1.py1[i__ + j * 150001 - 150002];
			pz = arprc1_1.pz1[i__ + j * 150001 - 150002];
			ee = arprc1_1.ee1[i__ + j * 150001 - 150002];
			xm = arprc1_1.xm1[i__ + j * 150001 - 150002];
			r__1 = px;
			r__2 = py;
			r__3 = xm;
			xmt = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
			if (dabs(pz) >= ee)
			{
				s_wsle(&io___648);
				do_lio(&c__9, &c__1, "IN HJANA4", (int)9);
				e_wsle();
				s_wsle(&io___649);
				do_lio(&c__9, &c__1, " PARTICLE ", (int)10);
				do_lio(&c__3, &c__1, (char *)&i__, (int)sizeof(int));
				do_lio(&c__9, &c__1, " RUN ", (int)5);
				do_lio(&c__3, &c__1, (char *)&j, (int)sizeof(int));
				do_lio(&c__9, &c__1, "PREC ERR", (int)8);
				e_wsle();
				s_wsle(&io___650);
				do_lio(&c__9, &c__1, " FLAV = ", (int)8);
				do_lio(&c__3, &c__1, (char *)&ityp, (int)sizeof(int));
				do_lio(&c__9, &c__1, " PX = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&px, (int)sizeof(float));
				do_lio(&c__9, &c__1, " PY = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&py, (int)sizeof(float));
				e_wsle();
				s_wsle(&io___651);
				do_lio(&c__9, &c__1, " PZ = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&pz, (int)sizeof(float));
				do_lio(&c__9, &c__1, " EE = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&ee, (int)sizeof(float));
				e_wsle();
				s_wsle(&io___652);
				do_lio(&c__9, &c__1, " XM = ", (int)6);
				do_lio(&c__4, &c__1, (char *)&xm, (int)sizeof(float));
				e_wsle();
			}
			dxmt = xmt - xm;
			y = log((ee + pz) / (ee - pz)) * .5f;

			iy = (int)((y + 10.f) / .2f) + 1;
			if (iy > 50)
			{
				if (y < -1.f || y >= 1.f)
			{
				continue;
			}
			imt = (int)(dxmt / .05f) + 1;
			if (imt > 50)
			{
				continue;
			}
			dmyh4[imt - 1] += 1.f / xmt;
			}
			dndyh4[iy - 1] += 1.f;
			deyh4[iy - 1] += xmt;			
		}
	}

	return 0;
}

int zpstrg_(void)
{

	int i__1;
	double d__1, d__2;
	double sqrt(double);
	static int i__, j;
	static float bb;
	static double tau7;
	static int nstr;
	static double shift;
	static int istrg;
	if (anim_1.isoft == 5)
	{
		i__1 = para1_1.mul;
		for (i__ = 1; i__ <= i__1; ++i__)
		{
			prec2_1.ityp5[i__ - 1] = frzprc_1.idfrz[i__ - 1];
			prec2_1.gx5[i__ - 1] = frzprc_1.gxfrz[i__ - 1];
			prec2_1.gy5[i__ - 1] = frzprc_1.gyfrz[i__ - 1];
			prec2_1.gz5[i__ - 1] = frzprc_1.gzfrz[i__ - 1];
			prec2_1.ft5[i__ - 1] = frzprc_1.ftfrz[i__ - 1];
			prec2_1.px5[i__ - 1] = frzprc_1.pxfrz[i__ - 1];
			prec2_1.py5[i__ - 1] = frzprc_1.pyfrz[i__ - 1];
			prec2_1.pz5[i__ - 1] = frzprc_1.pzfrz[i__ - 1];
			prec2_1.e5[i__ - 1] = frzprc_1.efrz[i__ - 1];
			prec2_1.xmass5[i__ - 1] = frzprc_1.xmfrz[i__ - 1];
		}
	}

	for (i__ = 1; i__ <= 150001; ++i__)
	{
		srec2_1.ataui[i__ - 1] = 0.;
		srec2_1.zt1[i__ - 1] = 0.;
		srec2_1.zt2[i__ - 1] = 0.;

		srec2_1.zt3[i__ - 1] = 0.;
		strg_1.np[i__ - 1] = 0;
	}
	i__1 = para1_1.mul;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		istrg = ilist8_1.lstrg1[i__ - 1];

		d__1 = prec2_1.ft5[i__ - 1];

		d__2 = prec2_1.gz5[i__ - 1];
		tau7 = sqrt(d__1 * d__1 - d__2 * d__2);
		srec2_1.ataui[istrg - 1] += tau7;
		srec2_1.zt1[istrg - 1] += prec2_1.gx5[i__ - 1];
		srec2_1.zt2[istrg - 1] += prec2_1.gy5[i__ - 1];
		srec2_1.zt3[istrg - 1] += prec2_1.gz5[i__ - 1];
		++strg_1.np[istrg - 1];
	}
	nstr = srec1_1.nsp + srec1_1.nst + srec1_1.nsi;

	if (anim_1.isoft == 3 || anim_1.isoft == 4 || anim_1.isoft == 5)
	{
		i__1 = nstr;
		for (i__ = 1; i__ <= i__1; ++i__)
		{
			if (strg_1.np[i__ - 1] != 0)
			{
				srec2_1.ataui[i__ - 1] /= strg_1.np[i__ - 1];
				srec2_1.zt1[i__ - 1] /= strg_1.np[i__ - 1];
				srec2_1.zt2[i__ - 1] /= strg_1.np[i__ - 1];
				srec2_1.zt3[i__ - 1] /= strg_1.np[i__ - 1];
			}
		}
		return 0;
	}

	i__1 = nstr;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		if (strg_1.np[i__ - 1] != 0)
		{
			srec2_1.ataui[i__ - 1] /= strg_1.np[i__ - 1];
			srec2_1.zt1[i__ - 1] /= strg_1.np[i__ - 1];
			srec2_1.zt2[i__ - 1] /= strg_1.np[i__ - 1];
			srec2_1.zt3[i__ - 1] /= strg_1.np[i__ - 1];
		}
		else
		{
			if (i__ <= srec1_1.nsp)
			{
				j = i__;
				srec2_1.zt1[i__ - 1] = (double)hjcrdn_1.yp[j * 3 - 3];
				srec2_1.zt2[i__ - 1] = (double)hjcrdn_1.yp[j * 3 - 2];
				srec2_1.zt3[i__ - 1] = 0.;
			}
			else if (i__ > srec1_1.nsp && i__ <= srec1_1.nsp + srec1_1.nst)
			{
				j = i__ - srec1_1.nsp;
				srec2_1.zt1[i__ - 1] = (double)hjcrdn_1.yt[j * 3 - 3];
				srec2_1.zt2[i__ - 1] = (double)hjcrdn_1.yt[j * 3 - 2];
				srec2_1.zt3[i__ - 1] = 0.;
			}
			else
			{
				j = i__ - srec1_1.nsp - srec1_1.nst;
				srec2_1.zt1[i__ - 1] = (double)(hjcrdn_1.yp[hjjet2_1.iasg[j - 1] * 3 - 3] + hjcrdn_1.yt[hjjet2_1.iasg[j + 150000] * 3 - 3]) * .5;
				srec2_1.zt2[i__ - 1] = (double)(hjcrdn_1.yp[hjjet2_1.iasg[j - 1] * 3 - 2] + hjcrdn_1.yt[hjjet2_1.iasg[j + 150000] * 3 - 2]) * .5;
				srec2_1.zt3[i__ - 1] = 0.;
			}
		}
	}

	bb = hparnt_1.hint1[18];
	i__1 = nstr;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		if (strg_1.np[i__ - 1] != 0)
		{
			shift = 0.;
		}
		else
		{
			shift = (double)bb * .5;
		}
		if (i__ <= srec1_1.nsp)
		{
			srec2_1.zt1[i__ - 1] += shift;
		}
		else if (i__ > srec1_1.nsp && i__ <= srec1_1.nsp + srec1_1.nst)
		{
			srec2_1.zt1[i__ - 1] -= shift;
		}
	}

	return 0;
}

double ranart_(int *nseed)
{

	float ret_val;
	extern double rand(int *);
	ret_val = rand(&c__0);
	return ret_val;
}
