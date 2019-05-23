#include "f2c.h"
struct
{
	int ioscar;
} para7_;

#define para7_1 para7_

struct
{
	int lblast[150001];
	float xlast[600004], plast[600004];
	int nlast;
} hbt_;

#define hbt_1 hbt_

struct
{
	int masspr, massta, iseed, iavoid;
	float dt;
} input1_;

#define input1_1 input1_

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
	int itimeh;
	float bimp;
} lastt_;

#define lastt_1 lastt_

struct
{
	float tfdcy[150001], tfdpi[150001], tft[150001];
} tdecay_;

#define tdecay_1 tdecay_

struct
{
	int iaevt, iarun;
} arevt_;

#define arevt_1 arevt_

struct
{
	float efrm;
	int npart1, npart2;
} snn_;

#define snn_1 snn_

struct
{
	int nelt, ninthj, nelp, ninp;
} hjglbr_;

#define hjglbr_1 hjglbr_

struct
{
	float ftsv[150001], ftsvt[150001];
} ftmax_;

#define ftmax_1 ftmax_

struct
{
	float dpertt[150001], dpertp[150001], dplast[150001],
		dpdcy[150001], dpdpi[150001], dpt[150001], dpp1[150001], dppion[150001];
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
	int idpert, npertd, idxsec;
} para8_;

#define para8_1 para8_

struct
{
	double enenew, pxnew, pynew, pznew;
} lor_;

#define lor_1 lor_

struct
{
	float ptwo[10];
} decom_;

#define decom_1 decom_

struct
{
	int nseed;
} rndf77_;

#define rndf77_1 rndf77_

struct
{
	int katt[600004];
	float patt[600004];
} hmain2_;

#define hmain2_1 hmain2_

struct
{
	int natt;
	float eatt;
	int jatt, nt, np, n0, n01, n10, n11;
} hmain1_;

#define hmain1_1 hmain1_

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
	int lstrg0[400001], lpart0[400001];
} ilist7_;

#define ilist7_1 ilist7_

struct
{
	int itypar[150001];
	float gxar[150001], gyar[150001], gzar[150001], ftar[150001], pxar[150001],
		pyar[150001], pzar[150001], pear[150001], xmar[150001];
} arprc_;

#define arprc_1 arprc_

struct
{
	int nnozpc, itypn[4001];
	float gxn[4001], gyn[4001], gzn[4001], ftn[4001], pxn[4001], pyn[4001],
		pzn[4001], een[4001], xmn[4001];
} noprec_;

#define noprec_1 noprec_

struct
{
	double pxsgs[450003], pysgs[450003], pzsgs[450003], pesgs[450003], pmsgs[450003], gxsgs[450003], gysgs[450003], gzsgs[450003], ftsgs[450003];
	int k1sgs[450003], k2sgs[450003], njsgs[150001];
} soft_;

#define soft_1 soft_

struct
{
	int nevent, isoft, isflag, izpc;
} anim_;

#define anim_1 anim_

struct
{
	double vxp0[400001], vyp0[400001], vzp0[400001];
} precpa_;

#define precpa_1 precpa_

struct
{
	double gxp[3], gyp[3], gzp[3], ftp[3], pxp[3], pyp[3], pzp[3], pep[3],
		pmp[3];
} loclco_;

#define loclco_1 loclco_

struct
{
	int nsg, njsg[150001], iasg[450003], k1sg[15000100], k2sg[15000100];
	float pxsg[15000100], pysg[15000100], pzsg[15000100], pesg[15000100], pmsg[15000100];
} hjjet2_;

#define hjjet2_1 hjjet2_

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
	double gxp0[3], gyp0[3], gzp0[3], ft0fom;
} prtn23_;

#define prtn23_1 prtn23_

struct
{
	int nattzp;
} nzpc_;

#define nzpc_1 nzpc_

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
	double dpcoal, drcoal, ecritl;
} coal_;

#define coal_1 coal_

struct
{
	int iap, izp, iat, izt;
} oscar1_;

#define oscar1_1 oscar1_

struct
{
	char frame[8], amptvn[25];
} oscar2_;

#define oscar2_1 oscar2_

struct
{
	int ilab, manyb, ntmax, icoll, insys, ipot, mode, imomen, nfreq,
		icflow, icrho, icou, kpoten, kmul;
} input2_;

#define input2_1 input2_

struct
{
	int n, k[45000];
	float p[45000], v[45000];
} lujets_;

#define lujets_1 lujets_

struct
{
	int kchg[1500];
	float pmas[2000], parf[2000], vckm[16];
} ludat2_;

#define ludat2_1 ludat2_

struct
{
	int mdcy[1500], mdme[4000];
	float brat[2000];
	int kfdp[10000];
} ludat3_;

#define ludat3_1 ludat3_

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
	int nsav, iksdcy;
} resdcy_;

#define resdcy_1 resdcy_

struct
{
	int lb1;
	float px1, py1, pz1, em1, e1, xfnl, yfnl, zfnl, tfnl, px1n, py1n, pz1n,
		dp1n;
} leadng_;

#define leadng_1 leadng_

struct
{
	double gx5[400001], gy5[400001], gz5[400001], ft5[400001], px5[400001], py5[400001], pz5[400001], e5[400001], xmass5[400001];
	int ityp5[400001];
} prec2_;

#define prec2_1 prec2_

struct
{
	double gxfrz[400001], gyfrz[400001], gzfrz[400001], ftfrz[400001],
		pxfrz[400001], pyfrz[400001], pzfrz[400001], efrz[400001], xmfrz[400001], tfrz[302];
	int ifrz[400001], idfrz[400001], itlast;
} frzprc_;

#define frzprc_1 frzprc_

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
	double ct[400001], ot[400001], tlarge;
} ilist5_;

#define ilist5_1 ilist5_

static int c__1 = 1;
static int c__9 = 9;
static int c__5 = 5;
static int c__3 = 3;
static int c_n1 = -1;
static int c_b205 = 150001;
static int c__2 = 2;
static int c__4 = 4;

int hbtout_(int *nnew, int *nt, int *ntmax)
{

	static char fmt_190[] = "(3(i7),f10.4,5x,6(i4))";
	static char fmt_200[] = "(i6,2(1x,f8.3),1x,f10.3,1x,f6.3,4(1x,f8.2))";
	static char fmt_250[] = "(i5,2(1x,f8.3),1x,f10.3,2(1x,f7.1),1x,f8.2,1x,f7.2,1x,e9.3)";
	static char fmt_201[] = "(i6,2(1x,f8.3),1x,f10.3,1x,f6.3,4(1x,e8.2))";
	static char fmt_251[] = "(i5,2(1x,f8.3),1x,f10.3,4(1x,e8.2),1x,e9.3)";
	int i__1, i__2;
	float r__1, r__2, r__3, r__4, r__5, r__6;
	double sqrt(double);
	int s_wsfe(cilist *), do_fio(int *, char *, ftnlen), e_wsfe(void),s_wsle(cilist *), do_lio(int *, int *, char *, ftnlen),e_wsle(void);

	static int i__, ii;
	static float dr;
	static int ip, ip2;
	static float ene, xnew[3];
	static int newkp[150001];
	static float deltat;
	extern int hoscar_(void);
	static int iplast, lastkp[150001], ndpert;
	extern int invflv_(int *);
	static cilist io___13 = {0, 16, 0, fmt_190, 0};
	static cilist io___14 = {0, 90, 0, fmt_190, 0};
	static cilist io___15 = {0, 16, 0, fmt_200, 0};
	static cilist io___16 = {0, 90, 0, fmt_250, 0};
	static cilist io___17 = {0, 99, 0, 0, 0};
	static cilist io___18 = {0, 16, 0, fmt_200, 0};
	static cilist io___19 = {0, 90, 0, fmt_250, 0};
	static cilist io___20 = {0, 99, 0, 0, 0};
	static cilist io___21 = {0, 16, 0, fmt_201, 0};
	static cilist io___22 = {0, 90, 0, fmt_251, 0};
	static cilist io___23 = {0, 99, 0, 0, 0};

	i__1 = max(hbt_1.nlast, *nnew);
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		lastkp[i__ - 1] = 0;
	}
	i__1 = *nnew;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		newkp[i__ - 1] = 0;
	}

	i__1 = *nnew;
	for (ip = 1; ip <= i__1; ++ip)
	{
		i__2 = hbt_1.nlast;
		for (iplast = 1; iplast <= i__2; ++iplast)
		{
			if (bb_1.p[ip * 3 - 3] == hbt_1.plast[(iplast << 2) - 4] &&
				bb_1.p[ip * 3 - 2] == hbt_1.plast[(iplast << 2) - 3] &&
				bb_1.p[ip * 3 - 1] == hbt_1.plast[(iplast << 2) - 2] &&
				cc_1.e[ip - 1] == hbt_1.plast[(iplast << 2) - 1] &&
				ee_1.lb[ip - 1] == hbt_1.lblast[iplast - 1] &&
				dpert_1.dpertp[ip - 1] == dpert_1.dplast[iplast - 1] &&
				lastkp[iplast - 1] == 0)
			{

				deltat = *nt * input1_1.dt - hbt_1.xlast[(iplast << 2) - 1];

				r__1 = hbt_1.plast[(iplast << 2) - 4];

				r__2 = hbt_1.plast[(iplast << 2) - 3];

				r__3 = hbt_1.plast[(iplast << 2) - 2];

				r__4 = hbt_1.plast[(iplast << 2) - 1];
				ene = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);

				for (ii = 1; ii <= 3; ++ii)
				{
					xnew[ii - 1] = hbt_1.xlast[ii + (iplast << 2) - 5] +hbt_1.plast[ii + (iplast << 2) - 5] / ene *deltat;
				}

				r__1 = aa_1.r__[ip * 3 - 3] - xnew[0];

				r__2 = aa_1.r__[ip * 3 - 2] - xnew[1];

				r__3 = aa_1.r__[ip * 3 - 1] - xnew[2];
				dr = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);

				if (dr <= .01f)
				{
					lastkp[iplast - 1] = 1;
					newkp[ip - 1] = 1;

					goto L100;
				}
			}
		}
	L100:;
	}

	i__1 = *nnew;
	for (ip = 1; ip <= i__1; ++ip)
	{
		if (newkp[ip - 1] == 0)
		{
			i__2 = *nnew;
			for (iplast = 1; iplast <= i__2; ++iplast)
			{
				if (lastkp[iplast - 1] == 0)
				{

					hbt_1.xlast[(iplast << 2) - 4] = aa_1.r__[ip * 3 - 3];
					hbt_1.xlast[(iplast << 2) - 3] = aa_1.r__[ip * 3 - 2];
					hbt_1.xlast[(iplast << 2) - 2] = aa_1.r__[ip * 3 - 1];
					hbt_1.xlast[(iplast << 2) - 1] = *nt * input1_1.dt;

					if (*nt == *ntmax)
					{

						if (tdecay_1.tfdcy[ip - 1] > *ntmax * input1_1.dt +
														 .001f)
						{
							hbt_1.xlast[(iplast << 2) - 1] = tdecay_1.tfdcy[ip - 1];
						}
						else if (ftmax_1.ftsv[ip - 1] > (*ntmax - 1) *
															input1_1.dt)
						{
							hbt_1.xlast[(iplast << 2) - 1] = ftmax_1.ftsv[ip - 1];
						}
					}
					hbt_1.plast[(iplast << 2) - 4] = bb_1.p[ip * 3 - 3];
					hbt_1.plast[(iplast << 2) - 3] = bb_1.p[ip * 3 - 2];
					hbt_1.plast[(iplast << 2) - 2] = bb_1.p[ip * 3 - 1];
					hbt_1.plast[(iplast << 2) - 1] = cc_1.e[ip - 1];
					hbt_1.lblast[iplast - 1] = ee_1.lb[ip - 1];
					lastkp[iplast - 1] = 1;

					dpert_1.dplast[iplast - 1] = dpert_1.dpertp[ip - 1];
					goto L150;
				}
			}
		}
	L150:;
	}

	if (*nnew < hbt_1.nlast)
	{
		i__1 = hbt_1.nlast;
		for (iplast = 1; iplast <= i__1; ++iplast)
		{
			if (lastkp[iplast - 1] == 0)
			{
				i__2 = hbt_1.nlast;
				for (ip2 = iplast + 1; ip2 <= i__2; ++ip2)
				{
					if (lastkp[ip2 - 1] == 1)
					{
						hbt_1.xlast[(iplast << 2) - 4] = hbt_1.xlast[(ip2 << 2) - 4];
						hbt_1.xlast[(iplast << 2) - 3] = hbt_1.xlast[(ip2 << 2) - 3];
						hbt_1.xlast[(iplast << 2) - 2] = hbt_1.xlast[(ip2 << 2) - 2];
						hbt_1.xlast[(iplast << 2) - 1] = hbt_1.xlast[(ip2 << 2) - 1];
						hbt_1.plast[(iplast << 2) - 4] = hbt_1.plast[(ip2 << 2) - 4];
						hbt_1.plast[(iplast << 2) - 3] = hbt_1.plast[(ip2 << 2) - 3];
						hbt_1.plast[(iplast << 2) - 2] = hbt_1.plast[(ip2 << 2) - 2];
						hbt_1.plast[(iplast << 2) - 1] = hbt_1.plast[(ip2 << 2) - 1];
						hbt_1.lblast[iplast - 1] = hbt_1.lblast[ip2 - 1];
						lastkp[iplast - 1] = 1;

						dpert_1.dplast[iplast - 1] = dpert_1.dplast[ip2 - 1];
						goto L170;
					}
				}
			}
		L170:;
		}
	}
	hbt_1.nlast = *nnew;

	if (*nt == *ntmax)
	{

		ndpert = 0;
		i__1 = hbt_1.nlast;
		for (ip = 1; ip <= i__1; ++ip)
		{
			if (dpert_1.dplast[ip - 1] > .99999f && dpert_1.dplast[ip - 1] <
														1.00001f)
			{
			}
			else
			{
				++ndpert;
			}
		}

		s_wsfe(&io___13);
		do_fio(&c__1, (char *)&arevt_1.iaevt, (ftnlen)sizeof(int));
		do_fio(&c__1, (char *)&arevt_1.iarun, (ftnlen)sizeof(int));
		i__1 = hbt_1.nlast - ndpert;
		do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(int));
		do_fio(&c__1, (char *)&lastt_1.bimp, (ftnlen)sizeof(float));
		do_fio(&c__1, (char *)&snn_1.npart1, (ftnlen)sizeof(int));
		do_fio(&c__1, (char *)&snn_1.npart2, (ftnlen)sizeof(int));
		do_fio(&c__1, (char *)&hjglbr_1.nelp, (ftnlen)sizeof(int));
		do_fio(&c__1, (char *)&hjglbr_1.ninp, (ftnlen)sizeof(int));
		do_fio(&c__1, (char *)&hjglbr_1.nelt, (ftnlen)sizeof(int));
		do_fio(&c__1, (char *)&hjglbr_1.ninthj, (ftnlen)sizeof(int));
		e_wsfe();

		if (para8_1.idpert == 1 || para8_1.idpert == 2)
		{
			s_wsfe(&io___14);
			do_fio(&c__1, (char *)&arevt_1.iaevt, (ftnlen)sizeof(int));
			do_fio(&c__1, (char *)&arevt_1.iarun, (ftnlen)sizeof(int));
			do_fio(&c__1, (char *)&ndpert, (ftnlen)sizeof(int));
			do_fio(&c__1, (char *)&lastt_1.bimp, (ftnlen)sizeof(float));
			do_fio(&c__1, (char *)&snn_1.npart1, (ftnlen)sizeof(int));
			do_fio(&c__1, (char *)&snn_1.npart2, (ftnlen)sizeof(int));
			do_fio(&c__1, (char *)&hjglbr_1.nelp, (ftnlen)sizeof(int));
			do_fio(&c__1, (char *)&hjglbr_1.ninp, (ftnlen)sizeof(int));
			do_fio(&c__1, (char *)&hjglbr_1.nelt, (ftnlen)sizeof(int));
			do_fio(&c__1, (char *)&hjglbr_1.ninthj, (ftnlen)sizeof(int));
			e_wsfe();
		}
		i__1 = hbt_1.nlast;
		for (ip = 1; ip <= i__1; ++ip)
		{

			r__1 = hbt_1.plast[(ip << 2) - 2];

			r__2 = hbt_1.plast[(ip << 2) - 1];
			if (hbt_1.plast[(ip << 2) - 4] == 0.f && hbt_1.plast[(ip << 2) - 3] == 0.f && sqrt(r__1 * r__1 + r__2 * r__2) * 2 / hparnt_1.hint1[0] > .99f && (hbt_1.lblast[ip - 1] == 1 || hbt_1.lblast[ip - 1] == 2))
			{

				if (dpert_1.dplast[ip - 1] > .99999f && dpert_1.dplast[ip - 1] < 1.00001f)
				{
					s_wsfe(&io___15);
					i__2 = invflv_(&hbt_1.lblast[ip - 1]);
					do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(int));
					do_fio(&c__1, (char *)&hbt_1.plast[(ip << 2) - 4], (ftnlen)sizeof(float));
					do_fio(&c__1, (char *)&hbt_1.plast[(ip << 2) - 3], (ftnlen)sizeof(float));
					do_fio(&c__1, (char *)&hbt_1.plast[(ip << 2) - 2], (ftnlen)sizeof(float));
					do_fio(&c__1, (char *)&hbt_1.plast[(ip << 2) - 1], (ftnlen)sizeof(float));
					do_fio(&c__1, (char *)&hbt_1.xlast[(ip << 2) - 4], (ftnlen)sizeof(float));
					do_fio(&c__1, (char *)&hbt_1.xlast[(ip << 2) - 3], (ftnlen)sizeof(float));
					do_fio(&c__1, (char *)&hbt_1.xlast[(ip << 2) - 2], (ftnlen)sizeof(float));

					r__2 = hbt_1.plast[(ip << 2) - 2];

					r__3 = hbt_1.plast[(ip << 2) - 1];
					r__1 = sqrt(r__2 * r__2 + r__3 * r__3) * 1e-20f /
						   hbt_1.plast[(ip << 2) - 1];
					do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(float));
					e_wsfe();
				}
				else
				{
					if (para8_1.idpert == 1 || para8_1.idpert == 2)
					{
						s_wsfe(&io___16);
						i__2 = invflv_(&hbt_1.lblast[ip - 1]);
						do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(int));
						do_fio(&c__1, (char *)&hbt_1.plast[(ip << 2) - 4], (ftnlen)sizeof(float));
						do_fio(&c__1, (char *)&hbt_1.plast[(ip << 2) - 3], (ftnlen)sizeof(float));
						do_fio(&c__1, (char *)&hbt_1.plast[(ip << 2) - 2], (ftnlen)sizeof(float));
						do_fio(&c__1, (char *)&hbt_1.xlast[(ip << 2) - 4], (ftnlen)sizeof(float));
						do_fio(&c__1, (char *)&hbt_1.xlast[(ip << 2) - 3], (ftnlen)sizeof(float));
						do_fio(&c__1, (char *)&hbt_1.xlast[(ip << 2) - 2], (ftnlen)sizeof(float));

						r__2 = hbt_1.plast[(ip << 2) - 2];

						r__3 = hbt_1.plast[(ip << 2) - 1];
						r__1 = sqrt(r__2 * r__2 + r__3 * r__3) * 1e-20f /
							   hbt_1.plast[(ip << 2) - 1];
						do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(float));
						do_fio(&c__1, (char *)&dpert_1.dplast[ip - 1], (ftnlen)sizeof(float));
						e_wsfe();
					}
					else
					{
						s_wsle(&io___17);
						do_lio(&c__9, &c__1, "Unexpected perturbative partic"
											 "les",
							   (ftnlen)33);
						e_wsle();
					}
				}
			}
			else
			{

				r__5 = (r__1 = hbt_1.xlast[(ip << 2) - 4], dabs(r__1)), r__6 = (r__2 = hbt_1.xlast[(ip << 2) - 3], dabs(r__2)),
				r__5 = max(r__5, r__6), r__6 = (r__3 = hbt_1.xlast[(ip << 2) - 2], dabs(r__3)), r__5 = max(r__5, r__6), r__6 = (r__4 = hbt_1.xlast[(ip << 2) - 1], dabs(r__4));
				if (dmax(r__5, r__6) < 9999.f)
				{
					if (dpert_1.dplast[ip - 1] > .99999f && dpert_1.dplast[ip - 1] < 1.00001f)
					{
						s_wsfe(&io___18);
						i__2 = invflv_(&hbt_1.lblast[ip - 1]);
						do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(int));
						do_fio(&c__1, (char *)&hbt_1.plast[(ip << 2) - 4], (ftnlen)sizeof(float));
						do_fio(&c__1, (char *)&hbt_1.plast[(ip << 2) - 3], (ftnlen)sizeof(float));
						do_fio(&c__1, (char *)&hbt_1.plast[(ip << 2) - 2], (ftnlen)sizeof(float));
						do_fio(&c__1, (char *)&hbt_1.plast[(ip << 2) - 1], (ftnlen)sizeof(float));
						do_fio(&c__1, (char *)&hbt_1.xlast[(ip << 2) - 4], (ftnlen)sizeof(float));
						do_fio(&c__1, (char *)&hbt_1.xlast[(ip << 2) - 3], (ftnlen)sizeof(float));
						do_fio(&c__1, (char *)&hbt_1.xlast[(ip << 2) - 2], (ftnlen)sizeof(float));
						do_fio(&c__1, (char *)&hbt_1.xlast[(ip << 2) - 1], (ftnlen)sizeof(float));
						e_wsfe();
					}
					else
					{
						if (para8_1.idpert == 1 || para8_1.idpert == 2)
						{
							s_wsfe(&io___19);
							i__2 = invflv_(&hbt_1.lblast[ip - 1]);
							do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(int));
							do_fio(&c__1, (char *)&hbt_1.plast[(ip << 2) - 4],
								   (ftnlen)sizeof(float));
							do_fio(&c__1, (char *)&hbt_1.plast[(ip << 2) - 3],
								   (ftnlen)sizeof(float));
							do_fio(&c__1, (char *)&hbt_1.plast[(ip << 2) - 2],
								   (ftnlen)sizeof(float));
							do_fio(&c__1, (char *)&hbt_1.xlast[(ip << 2) - 4],
								   (ftnlen)sizeof(float));
							do_fio(&c__1, (char *)&hbt_1.xlast[(ip << 2) - 3],
								   (ftnlen)sizeof(float));
							do_fio(&c__1, (char *)&hbt_1.xlast[(ip << 2) - 2],
								   (ftnlen)sizeof(float));
							do_fio(&c__1, (char *)&hbt_1.xlast[(ip << 2) - 1],
								   (ftnlen)sizeof(float));
							do_fio(&c__1, (char *)&dpert_1.dplast[ip - 1], (ftnlen)sizeof(float));
							e_wsfe();
						}
						else
						{
							s_wsle(&io___20);
							do_lio(&c__9, &c__1, "Unexpected perturbative pa"
												 "rticles",
								   (ftnlen)33);
							e_wsle();
						}
					}
				}
				else
				{

					if (dpert_1.dplast[ip - 1] > .99999f && dpert_1.dplast[ip - 1] < 1.00001f)
					{
						s_wsfe(&io___21);
						i__2 = invflv_(&hbt_1.lblast[ip - 1]);
						do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(int));
						do_fio(&c__1, (char *)&hbt_1.plast[(ip << 2) - 4], (ftnlen)sizeof(float));
						do_fio(&c__1, (char *)&hbt_1.plast[(ip << 2) - 3], (ftnlen)sizeof(float));
						do_fio(&c__1, (char *)&hbt_1.plast[(ip << 2) - 2], (ftnlen)sizeof(float));
						do_fio(&c__1, (char *)&hbt_1.plast[(ip << 2) - 1], (ftnlen)sizeof(float));
						do_fio(&c__1, (char *)&hbt_1.xlast[(ip << 2) - 4], (ftnlen)sizeof(float));
						do_fio(&c__1, (char *)&hbt_1.xlast[(ip << 2) - 3], (ftnlen)sizeof(float));
						do_fio(&c__1, (char *)&hbt_1.xlast[(ip << 2) - 2], (ftnlen)sizeof(float));
						do_fio(&c__1, (char *)&hbt_1.xlast[(ip << 2) - 1], (ftnlen)sizeof(float));
						e_wsfe();
					}
					else
					{
						if (para8_1.idpert == 1 || para8_1.idpert == 2)
						{
							s_wsfe(&io___22);
							i__2 = invflv_(&hbt_1.lblast[ip - 1]);
							do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(int));
							do_fio(&c__1, (char *)&hbt_1.plast[(ip << 2) - 4],
								   (ftnlen)sizeof(float));
							do_fio(&c__1, (char *)&hbt_1.plast[(ip << 2) - 3],
								   (ftnlen)sizeof(float));
							do_fio(&c__1, (char *)&hbt_1.plast[(ip << 2) - 2],
								   (ftnlen)sizeof(float));
							do_fio(&c__1, (char *)&hbt_1.xlast[(ip << 2) - 4],
								   (ftnlen)sizeof(float));
							do_fio(&c__1, (char *)&hbt_1.xlast[(ip << 2) - 3],
								   (ftnlen)sizeof(float));
							do_fio(&c__1, (char *)&hbt_1.xlast[(ip << 2) - 2],
								   (ftnlen)sizeof(float));
							do_fio(&c__1, (char *)&hbt_1.xlast[(ip << 2) - 1],
								   (ftnlen)sizeof(float));
							do_fio(&c__1, (char *)&dpert_1.dplast[ip - 1], (ftnlen)sizeof(float));
							e_wsfe();
						}
						else
						{
							s_wsle(&io___23);
							do_lio(&c__9, &c__1, "Unexpected perturbative pa"
												 "rticles",
								   (ftnlen)33);
							e_wsle();
						}
					}
				}
			}
		}
		if (para7_1.ioscar == 1)
		{
			hoscar_();
		}
	}

	return 0;
}

int decomp_(float *px0, float *py0, float *pz0, float *xm0)
{

	double d__1, d__2, d__3, d__4, d__5, d__6;

	double sqrt(double), cos(double), sin(double);
	int s_wsle(cilist *), do_lio(int *, int *, char *, ftnlen),
		e_wsle(void);

	static double ds, de0, de1, de2, gam, dpx, dpy, dpz, dcth, dpcm, dphi,
		dbex, dbey, dbez, beta2;
	extern double ranart_(int *);
	extern int lorenz_(double *, double *,
					   double *, double *, double *, double *,
					   double *);

	static cilist io___39 = {0, 6, 0, 0, 0};

	dcth = (double)ranart_(&rndf77_1.nseed) * 2. - 1.;
	dphi = (double)(ranart_(&rndf77_1.nseed) * hparnt_1.hipr1[39]) * 2.;

	d__1 = (double)(*xm0);
	ds = d__1 * d__1;

	d__1 = (double)(decom_1.ptwo[8] + decom_1.ptwo[9]);

	d__2 = (double)(decom_1.ptwo[8] - decom_1.ptwo[9]);
	dpcm = sqrt((ds - d__1 * d__1) * (ds - d__2 * d__2) / ds / 4.);
	dpz = dpcm * dcth;

	d__1 = dcth;
	dpx = dpcm * sqrt(1. - d__1 * d__1) * cos(dphi);

	d__1 = dcth;
	dpy = dpcm * sqrt(1. - d__1 * d__1) * sin(dphi);

	d__1 = (double)decom_1.ptwo[8];

	d__2 = dpcm;
	de1 = sqrt(d__1 * d__1 + d__2 * d__2);

	d__1 = (double)decom_1.ptwo[9];

	d__2 = dpcm;
	de2 = sqrt(d__1 * d__1 + d__2 * d__2);

	d__1 = (double)(*px0);

	d__2 = (double)(*py0);

	d__3 = (double)(*pz0);

	d__4 = (double)(*xm0);
	de0 = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3 + d__4 * d__4);
	dbex = (double)(*px0) / de0;
	dbey = (double)(*py0) / de0;
	dbez = (double)(*pz0) / de0;

	d__1 = dbex;

	d__2 = dbey;

	d__3 = dbez;
	beta2 = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
	gam = 1. / sqrt(1. - beta2);
	if (beta2 >= .9999999999999)
	{
		s_wsle(&io___39);
		do_lio(&c__9, &c__1, "1", (ftnlen)1);
		do_lio(&c__5, &c__1, (char *)&dbex, (ftnlen)sizeof(double));
		do_lio(&c__5, &c__1, (char *)&dbey, (ftnlen)sizeof(double));
		do_lio(&c__5, &c__1, (char *)&dbez, (ftnlen)sizeof(double));
		do_lio(&c__5, &c__1, (char *)&beta2, (ftnlen)sizeof(double));
		do_lio(&c__5, &c__1, (char *)&gam, (ftnlen)sizeof(double));
		e_wsle();
	}

	d__1 = -dbex;
	d__2 = -dbey;
	d__3 = -dbez;
	lorenz_(&de1, &dpx, &dpy, &dpz, &d__1, &d__2, &d__3);
	decom_1.ptwo[0] = (float)lor_1.pxnew;
	decom_1.ptwo[2] = (float)lor_1.pynew;
	decom_1.ptwo[4] = (float)lor_1.pznew;
	decom_1.ptwo[6] = (float)lor_1.enenew;

	d__1 = -dpx;
	d__2 = -dpy;
	d__3 = -dpz;
	d__4 = -dbex;
	d__5 = -dbey;
	d__6 = -dbez;
	lorenz_(&de2, &d__1, &d__2, &d__3, &d__4, &d__5, &d__6);
	decom_1.ptwo[1] = (float)lor_1.pxnew;
	decom_1.ptwo[3] = (float)lor_1.pynew;
	decom_1.ptwo[5] = (float)lor_1.pznew;
	decom_1.ptwo[7] = (float)lor_1.enenew;

	return 0;
}

int htop_(void)
{

	static char fmt_200[] = "(i6,2(1x,f8.3),1x,f10.3,1x,f6.3,4(1x,f8.2))";
	static char fmt_201[] = "(i6,2(1x,f8.3),1x,f10.3,1x,f6.3,4(1x,e8.2))";

	int i__1, i__2;
	float r__1, r__2, r__3;
	double d__1, d__2, d__3, d__4, d__5, d__6;

	int s_wsle(cilist *), do_lio(int *, int *, char *, ftnlen),
		e_wsle(void), i_sign(int *, int *), pow_ii(int *, int *), s_wsfe(cilist *), do_fio(int *, char *, ftnlen),
		e_wsfe(void);

	static int nsmbbbar, nsmmeson, i__, i1, i2, i3, i4, id, it[4], ipar,
		npar;
	static float xmdq, rnum;
	static int idabs;
	static float ftime, ptwox, ptwoy, ptwoz;
	extern int decomp_(float *, float *, float *, float *);
	static int ipamax;
	extern double ranart_(int *);
	static int inozpc;
	extern double ulmass_(int *);

	static cilist io___47 = {0, 92, 0, 0, 0};
	static cilist io___48 = {0, 92, 0, 0, 0};
	static cilist io___49 = {0, 92, 0, 0, 0};
	static cilist io___50 = {0, 92, 0, 0, 0};
	static cilist io___51 = {0, 92, 0, 0, 0};
	static cilist io___61 = {0, 92, 0, fmt_200, 0};
	static cilist io___62 = {0, 92, 0, fmt_201, 0};
	static cilist io___67 = {0, 92, 0, fmt_200, 0};
	static cilist io___68 = {0, 92, 0, fmt_201, 0};
	static cilist io___69 = {0, 20, 0, 0, 0};
	static cilist io___70 = {0, 20, 0, 0, 0};

	npar = 0;
	noprec_1.nnozpc = 0;

	if ((anim_1.isoft == 4 || anim_1.isoft == 5) && para7_1.ioscar == 2)
	{
		nsmbbbar = 0;
		nsmmeson = 0;
		i__1 = hmain1_1.natt;
		for (i__ = 1; i__ <= i__1; ++i__)
		{
			id = arprc_1.itypar[i__ - 1];
			idabs = abs(id);
			i2 = idabs / 10 % 10;
			if (arprc_1.pxar[i__ - 1] == 0.f && arprc_1.pyar[i__ - 1] == 0.f && arprc_1.pear[i__ - 1] >= hparnt_1.hint1[0] / 2 * .99f && (id == 2112 || id == 2212))
			{
			}
			else if (idabs > 1000 && i2 != 0)
			{

				++nsmbbbar;
			}
			else if (idabs > 100 && idabs < 1000 || idabs > 10000)
			{

				++nsmmeson;
			}
		}
		s_wsle(&io___47);
		i__1 = nsmbbbar * 3 + (nsmmeson << 1);
		do_lio(&c__3, &c__1, (char *)&i__1, (ftnlen)sizeof(int));
		e_wsle();
		s_wsle(&io___48);
		do_lio(&c__9, &c__1, " is the total # of initial partons after strin"
							 "g  melting",
			   (ftnlen)56);
		e_wsle();
		s_wsle(&io___49);
		do_lio(&c__9, &c__1, "String melting converts ", (ftnlen)24);
		do_lio(&c__3, &c__1, (char *)&nsmbbbar, (ftnlen)sizeof(int));
		do_lio(&c__9, &c__1, " baryons &", (ftnlen)10);
		do_lio(&c__3, &c__1, (char *)&nsmmeson, (ftnlen)sizeof(int));
		do_lio(&c__9, &c__1, " mesons", (ftnlen)7);
		e_wsle();
		s_wsle(&io___50);
		do_lio(&c__9, &c__1, "Total # of initial particles= ", (ftnlen)30);
		do_lio(&c__3, &c__1, (char *)&hmain1_1.natt, (ftnlen)sizeof(int));
		e_wsle();
		s_wsle(&io___51);
		do_lio(&c__9, &c__1, "Total # of initial particles (gamma,e,muon,..."
							 ")  not entering ZPC= ",
			   (ftnlen)67);
		i__1 = hmain1_1.natt - nsmbbbar - nsmmeson;
		do_lio(&c__3, &c__1, (char *)&i__1, (ftnlen)sizeof(int));
		e_wsle();
	}

	i__1 = hmain1_1.natt;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		id = arprc_1.itypar[i__ - 1];
		idabs = abs(id);
		i4 = idabs / 1000 % 10;
		i3 = idabs / 100 % 10;
		i2 = idabs / 10 % 10;
		i1 = idabs % 10;
		rnum = ranart_(&rndf77_1.nseed);

		r__1 = arprc_1.pxar[i__ - 1];

		r__2 = arprc_1.pyar[i__ - 1];

		r__3 = arprc_1.xmar[i__ - 1];
		ftime = arprc_1.pear[i__ - 1] * .197f / (r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
		inozpc = 0;
		it[0] = 0;
		it[1] = 0;
		it[2] = 0;
		it[3] = 0;

		if (arprc_1.pxar[i__ - 1] == 0.f && arprc_1.pyar[i__ - 1] == 0.f &&
			arprc_1.pear[i__ - 1] >= hparnt_1.hint1[0] / 2 * .99f && (id == 2112 || id == 2212))
		{

			inozpc = 1;
		}
		else if (idabs > 1000 && i2 != 0)
		{

			if ((i4 == 1 || i4 == 2) && i4 == i3 || i4 == 3 && i3 == 3)
			{
				if (i1 == 2)
				{
					if (rnum <= .5f)
					{
						it[0] = i4;
						it[1] = i3 * 1000 + i2 * 100 + 1;
					}
					else if (rnum <= .66666666666666663f)
					{
						it[0] = i4;
						it[1] = i3 * 1000 + i2 * 100 + 3;
					}
					else
					{
						it[0] = i2;
						it[1] = i4 * 1000 + i3 * 100 + 3;
					}
				}
				else if (i1 == 4)
				{
					if (rnum <= .66666666666666663f)
					{
						it[0] = i4;
						it[1] = i3 * 1000 + i2 * 100 + 3;
					}
					else
					{
						it[0] = i2;
						it[1] = i4 * 1000 + i3 * 100 + 3;
					}
				}
			}
			else if (i4 == 1 || i4 == 2)
			{
				if (i1 == 2)
				{
					if (rnum <= .5f)
					{
						it[0] = i2;
						it[1] = i4 * 1000 + i3 * 100 + 1;
					}
					else if (rnum <= .66666666666666663f)
					{
						it[0] = i2;
						it[1] = i4 * 1000 + i3 * 100 + 3;
					}
					else
					{
						it[0] = i4;
						it[1] = i3 * 1000 + i2 * 100 + 3;
					}
				}
				else if (i1 == 4)
				{
					if (rnum <= .66666666666666663f)
					{
						it[0] = i2;
						it[1] = i4 * 1000 + i3 * 100 + 3;
					}
					else
					{
						it[0] = i4;
						it[1] = i3 * 1000 + i2 * 100 + 3;
					}
				}
			}
			else if (i4 >= 3)
			{
				it[0] = i4;
				if (i3 < i2)
				{
					it[1] = i2 * 1000 + i3 * 100 + 1;
				}
				else
				{
					it[1] = i3 * 1000 + i2 * 100 + 3;
				}
			}

			if (id < 0)
			{
				it[0] = -it[0];
				it[1] = -it[1];
			}

			if (anim_1.isoft == 4 || anim_1.isoft == 5)
			{
				it[2] = it[1] / 1000 % 10;
				it[3] = it[1] / 100 % 10;
			}
		}
		else if (idabs > 100 && idabs < 1000 || idabs > 10000)
		{

			if (i3 == i2)
			{
				if (i3 == 1 || i3 == 2)
				{
					if (rnum <= .5f)
					{
						it[0] = 1;
						it[1] = -1;
					}
					else
					{
						it[0] = 2;
						it[1] = -2;
					}
				}
				else
				{
					it[0] = i3;
					it[1] = -i3;
				}
			}
			else
			{
				if (i_sign(&c__1, &id) * pow_ii(&c_n1, &i3) == 1)
				{
					it[0] = i3;
					it[1] = -i2;
				}
				else
				{
					it[0] = i2;
					it[1] = -i3;
				}
			}
		}
		else
		{

			inozpc = 1;
		}

		if (inozpc == 1)
		{
			soft_1.njsgs[i__ - 1] = 0;
			++noprec_1.nnozpc;
			noprec_1.itypn[noprec_1.nnozpc - 1] = arprc_1.itypar[i__ - 1];
			noprec_1.pxn[noprec_1.nnozpc - 1] = arprc_1.pxar[i__ - 1];
			noprec_1.pyn[noprec_1.nnozpc - 1] = arprc_1.pyar[i__ - 1];
			noprec_1.pzn[noprec_1.nnozpc - 1] = arprc_1.pzar[i__ - 1];
			noprec_1.een[noprec_1.nnozpc - 1] = arprc_1.pear[i__ - 1];
			noprec_1.xmn[noprec_1.nnozpc - 1] = arprc_1.xmar[i__ - 1];
			noprec_1.gxn[noprec_1.nnozpc - 1] = arprc_1.gxar[i__ - 1];
			noprec_1.gyn[noprec_1.nnozpc - 1] = arprc_1.gyar[i__ - 1];
			noprec_1.gzn[noprec_1.nnozpc - 1] = arprc_1.gzar[i__ - 1];
			noprec_1.ftn[noprec_1.nnozpc - 1] = arprc_1.ftar[i__ - 1];
		}
		else
		{
			soft_1.njsgs[i__ - 1] = 2;
			decom_1.ptwo[8] = ulmass_(it);
			decom_1.ptwo[9] = ulmass_(&it[1]);
			decomp_(&hmain2_1.patt[i__ - 1], &hmain2_1.patt[i__ + 150000], &hmain2_1.patt[i__ + 300001], &arprc_1.xmar[i__ - 1]);
			ipamax = 2;
			if ((anim_1.isoft == 4 || anim_1.isoft == 5) && abs(it[1]) > 1000)
			{
				ipamax = 1;
			}
			i__2 = ipamax;
			for (ipar = 1; ipar <= i__2; ++ipar)
			{
				++npar;
				prec1_1.ityp0[npar - 1] = it[ipar - 1];
				prec1_1.px0[npar - 1] = (double)decom_1.ptwo[ipar - 1];
				prec1_1.py0[npar - 1] = (double)decom_1.ptwo[ipar + 1];
				prec1_1.pz0[npar - 1] = (double)decom_1.ptwo[ipar + 3];
				prec1_1.e0[npar - 1] = (double)decom_1.ptwo[ipar + 5];
				prec1_1.xmass0[npar - 1] = (double)decom_1.ptwo[ipar + 7];
				prec1_1.gx0[npar - 1] = (double)arprc_1.gxar[i__ - 1];
				prec1_1.gy0[npar - 1] = (double)arprc_1.gyar[i__ - 1];
				prec1_1.gz0[npar - 1] = (double)arprc_1.gzar[i__ - 1];
				prec1_1.ft0[npar - 1] = (double)ftime;
				ilist7_1.lstrg0[npar - 1] = i__;
				ilist7_1.lpart0[npar - 1] = ipar;
				precpa_1.vxp0[npar - 1] = (double)(hmain2_1.patt[i__ - 1] / hmain2_1.patt[i__ + 450002]);
				precpa_1.vyp0[npar - 1] = (double)(hmain2_1.patt[i__ +
																	 150000] /
													   hmain2_1.patt[i__ + 450002]);
				precpa_1.vzp0[npar - 1] = (double)(hmain2_1.patt[i__ +
																	 300001] /
													   hmain2_1.patt[i__ + 450002]);

				if ((anim_1.isoft == 4 || anim_1.isoft == 5) &&
					para7_1.ioscar == 2)
				{

					d__5 = (d__1 = prec1_1.gx0[npar - 1], abs(d__1)), d__6 = (d__2 = prec1_1.gy0[npar - 1], abs(d__2)), d__5 = max(d__5, d__6), d__6 = (d__3 = prec1_1.gz0[npar - 1], abs(d__3)), d__5 = max(d__5, d__6), d__6 = (d__4 = prec1_1.ft0[npar - 1], abs(d__4));
					if (max(d__5, d__6) < 9999.)
					{
						s_wsfe(&io___61);
						do_fio(&c__1, (char *)&prec1_1.ityp0[npar - 1], (ftnlen)sizeof(int));
						do_fio(&c__1, (char *)&prec1_1.px0[npar - 1], (ftnlen)sizeof(double));
						do_fio(&c__1, (char *)&prec1_1.py0[npar - 1], (ftnlen)sizeof(double));
						do_fio(&c__1, (char *)&prec1_1.pz0[npar - 1], (ftnlen)sizeof(double));
						do_fio(&c__1, (char *)&prec1_1.xmass0[npar - 1], (ftnlen)sizeof(double));
						do_fio(&c__1, (char *)&prec1_1.gx0[npar - 1], (ftnlen)sizeof(double));
						do_fio(&c__1, (char *)&prec1_1.gy0[npar - 1], (ftnlen)sizeof(double));
						do_fio(&c__1, (char *)&prec1_1.gz0[npar - 1], (ftnlen)sizeof(double));
						do_fio(&c__1, (char *)&prec1_1.ft0[npar - 1], (ftnlen)sizeof(double));
						e_wsfe();
					}
					else
					{
						s_wsfe(&io___62);
						do_fio(&c__1, (char *)&prec1_1.ityp0[npar - 1], (ftnlen)sizeof(int));
						do_fio(&c__1, (char *)&prec1_1.px0[npar - 1], (ftnlen)sizeof(double));
						do_fio(&c__1, (char *)&prec1_1.py0[npar - 1], (ftnlen)sizeof(double));
						do_fio(&c__1, (char *)&prec1_1.pz0[npar - 1], (ftnlen)sizeof(double));
						do_fio(&c__1, (char *)&prec1_1.xmass0[npar - 1], (ftnlen)sizeof(double));
						do_fio(&c__1, (char *)&prec1_1.gx0[npar - 1], (ftnlen)sizeof(double));
						do_fio(&c__1, (char *)&prec1_1.gy0[npar - 1], (ftnlen)sizeof(double));
						do_fio(&c__1, (char *)&prec1_1.gz0[npar - 1], (ftnlen)sizeof(double));
						do_fio(&c__1, (char *)&prec1_1.ft0[npar - 1], (ftnlen)sizeof(double));
						e_wsfe();
					}
				}
			}

			if ((anim_1.isoft == 4 || anim_1.isoft == 5) && abs(it[1]) > 1000)
			{
				soft_1.njsgs[i__ - 1] = 3;
				xmdq = decom_1.ptwo[9];
				decom_1.ptwo[8] = ulmass_(&it[2]);
				decom_1.ptwo[9] = ulmass_(&it[3]);

				ptwox = decom_1.ptwo[1];
				ptwoy = decom_1.ptwo[3];
				ptwoz = decom_1.ptwo[5];
				decomp_(&ptwox, &ptwoy, &ptwoz, &xmdq);

				for (ipar = 1; ipar <= 2; ++ipar)
				{
					++npar;
					prec1_1.ityp0[npar - 1] = it[ipar + 1];
					prec1_1.px0[npar - 1] = (double)decom_1.ptwo[ipar -
																	 1];
					prec1_1.py0[npar - 1] = (double)decom_1.ptwo[ipar +
																	 1];
					prec1_1.pz0[npar - 1] = (double)decom_1.ptwo[ipar +
																	 3];
					prec1_1.e0[npar - 1] = (double)decom_1.ptwo[ipar + 5];
					prec1_1.xmass0[npar - 1] = (double)decom_1.ptwo[ipar + 7];
					prec1_1.gx0[npar - 1] = (double)arprc_1.gxar[i__ - 1];
					prec1_1.gy0[npar - 1] = (double)arprc_1.gyar[i__ - 1];
					prec1_1.gz0[npar - 1] = (double)arprc_1.gzar[i__ - 1];
					prec1_1.ft0[npar - 1] = (double)ftime;
					ilist7_1.lstrg0[npar - 1] = i__;
					ilist7_1.lpart0[npar - 1] = ipar + 1;
					precpa_1.vxp0[npar - 1] = (double)(hmain2_1.patt[i__ - 1] / hmain2_1.patt[i__ + 450002]);
					precpa_1.vyp0[npar - 1] = (double)(hmain2_1.patt[i__ + 150000] / hmain2_1.patt[i__ + 450002]);
					precpa_1.vzp0[npar - 1] = (double)(hmain2_1.patt[i__ + 300001] / hmain2_1.patt[i__ + 450002]);

					if ((anim_1.isoft == 4 || anim_1.isoft == 5) &&
						para7_1.ioscar == 2)
					{

						d__5 = (d__1 = prec1_1.gx0[npar - 1], abs(d__1)),
						d__6 = (d__2 = prec1_1.gy0[npar - 1], abs(
																  d__2)),
						d__5 = max(d__5, d__6), d__6 = (d__3 = prec1_1.gz0[npar - 1], abs(d__3)), d__5 = max(d__5, d__6), d__6 = (d__4 = prec1_1.ft0[npar - 1], abs(d__4));
						if (max(d__5, d__6) < 9999.)
						{
							s_wsfe(&io___67);
							do_fio(&c__1, (char *)&prec1_1.ityp0[npar - 1], (ftnlen)sizeof(int));
							do_fio(&c__1, (char *)&prec1_1.px0[npar - 1], (ftnlen)sizeof(double));
							do_fio(&c__1, (char *)&prec1_1.py0[npar - 1], (ftnlen)sizeof(double));
							do_fio(&c__1, (char *)&prec1_1.pz0[npar - 1], (ftnlen)sizeof(double));
							do_fio(&c__1, (char *)&prec1_1.xmass0[npar - 1], (ftnlen)sizeof(double));
							do_fio(&c__1, (char *)&prec1_1.gx0[npar - 1], (ftnlen)sizeof(double));
							do_fio(&c__1, (char *)&prec1_1.gy0[npar - 1], (ftnlen)sizeof(double));
							do_fio(&c__1, (char *)&prec1_1.gz0[npar - 1], (ftnlen)sizeof(double));
							do_fio(&c__1, (char *)&prec1_1.ft0[npar - 1], (ftnlen)sizeof(double));
							e_wsfe();
						}
						else
						{
							s_wsfe(&io___68);
							do_fio(&c__1, (char *)&prec1_1.ityp0[npar - 1], (ftnlen)sizeof(int));
							do_fio(&c__1, (char *)&prec1_1.px0[npar - 1], (ftnlen)sizeof(double));
							do_fio(&c__1, (char *)&prec1_1.py0[npar - 1], (ftnlen)sizeof(double));
							do_fio(&c__1, (char *)&prec1_1.pz0[npar - 1], (ftnlen)sizeof(double));
							do_fio(&c__1, (char *)&prec1_1.xmass0[npar - 1], (ftnlen)sizeof(double));
							do_fio(&c__1, (char *)&prec1_1.gx0[npar - 1], (ftnlen)sizeof(double));
							do_fio(&c__1, (char *)&prec1_1.gy0[npar - 1], (ftnlen)sizeof(double));
							do_fio(&c__1, (char *)&prec1_1.gz0[npar - 1], (ftnlen)sizeof(double));
							do_fio(&c__1, (char *)&prec1_1.ft0[npar - 1], (ftnlen)sizeof(double));
							e_wsfe();
						}
					}
				}
			}
		}
	}
	para1_1.mul = npar;

	if ((anim_1.isoft == 4 || anim_1.isoft == 5) && para7_1.ioscar == 2)
	{
		if (hmain1_1.natt - nsmbbbar - nsmmeson != noprec_1.nnozpc)
		{
			s_wsle(&io___69);
			do_lio(&c__9, &c__1, "Problem with the total # of initial partic"
								 "les (gamma,e,muon,...) not entering ZPC",
				   (ftnlen)81);
			e_wsle();
		}
		if (nsmbbbar * 3 + (nsmmeson << 1) != npar)
		{
			s_wsle(&io___70);
			do_lio(&c__9, &c__1, "Problem with the total # of initial parton"
								 "s   after string melting",
				   (ftnlen)66);
			e_wsle();
		}
	}

	return 0;
}

int ptoh_(void)
{

	int i__1, i__2, i__3;
	float r__1, r__2, r__3, r__4;
	double d__1, d__2, d__3, d__4;

	double sqrt(double);
	int i_sign(int *, int *), pow_ii(int *, int *),
		s_wsle(cilist *), do_lio(int *, int *, char *, ftnlen),
		e_wsle(void);

	static int i__;
	static double e1, e2, e3;
	static int k1, k2, k3, kf, ki, kj, kk, ix, kf1, kf2;
	static double px1, py1, pz1, px2, py2, pz2, px3, py3, pz3, gam;
	static int ibs;
	static double bex, bey, bez;
	static int kdq, isg, npi0;
	static float ppi0, tau0;
	static int kmax;
	static double etot;
	static int indx[150001], kmin;
	static double beta2;
	static int k1abs, k2abs, k3abs;
	static float prho0;
	static int ndiag[150001], npich, iuudd, nuudd, ktemp, inatt, mstj24;
	static double ftavg0;
	extern int index1_(int *, int *, double *,
					   int *);
	static double gxavg0, gyavg0, gzavg0;
	extern int coales_(void);
	static double xmdiag[150001], drlocl;
	extern int locldr_(int *, double *);
	static int nrhoch;
	extern double ranart_(int *);
	static int ipartn, idqspn;
	static double xmpair;
	static int imspin;
	extern double ulmass_(int *);
	extern int lorenz_(double *, double *,
					   double *, double *, double *, double *,
					   double *);

	static cilist io___131 = {0, 6, 0, 0, 0};

	coales_();

	mstj24 = ludat1_1.mstj[23];
	ludat1_1.mstj[23] = 0;
	nuudd = 0;
	npich = 0;
	nrhoch = 0;
	ppi0 = 1.f;
	prho0 = 0.f;

	i__1 = hjjet2_1.nsg;
	for (isg = 1; isg <= i__1; ++isg)
	{
		if (soft_1.njsgs[isg - 1] != 0)
		{
			++hmain1_1.natt;
			k1 = soft_1.k2sgs[isg - 1];
			k1abs = abs(k1);
			px1 = soft_1.pxsgs[isg - 1];
			py1 = soft_1.pysgs[isg - 1];
			pz1 = soft_1.pzsgs[isg - 1];
			k2 = soft_1.k2sgs[isg + 150000];
			k2abs = abs(k2);
			px2 = soft_1.pxsgs[isg + 150000];
			py2 = soft_1.pysgs[isg + 150000];
			pz2 = soft_1.pzsgs[isg + 150000];

			e1 = soft_1.pesgs[isg - 1];
			e2 = soft_1.pesgs[isg + 150000];

			d__1 = e1 + e2;

			d__2 = px1 + px2;

			d__3 = py1 + py2;

			d__4 = pz1 + pz2;
			xmpair = sqrt(d__1 * d__1 - d__2 * d__2 - d__3 * d__3 - d__4 * d__4);
			ibs = 2;
			imspin = 0;
			if (k1 == -k2 && abs(k1) <= 2 && soft_1.njsgs[isg - 1] == 2)
			{
				++nuudd;
				xmdiag[nuudd - 1] = xmpair;
				ndiag[nuudd - 1] = hmain1_1.natt;
			}
			k3 = 0;
			if ((anim_1.isoft == 4 || anim_1.isoft == 5) && soft_1.njsgs[isg - 1] == 3)
			{
				k3 = soft_1.k2sgs[isg + 300001];
				k3abs = abs(k3);
				px3 = soft_1.pxsgs[isg + 300001];
				py3 = soft_1.pysgs[isg + 300001];
				pz3 = soft_1.pzsgs[isg + 300001];
				e3 = soft_1.pesgs[isg + 300001];

				d__1 = e1 + e2 + e3;

				d__2 = px1 + px2 + px3;

				d__3 = py1 + py2 + py3;

				d__4 = pz1 + pz2 + pz3;
				xmpair = sqrt(d__1 * d__1 - d__2 * d__2 - d__3 * d__3 - d__4 * d__4);
			}

			if (anim_1.isoft == 3 && (k1abs > 1000 || k2abs > 1000))
			{
				if (k1abs > 1000)
				{
					kdq = k1abs;
					kk = k2abs;
				}
				else
				{
					kdq = k2abs;
					kk = k1abs;
				}
				ki = kdq / 1000 % 10;
				kj = kdq / 100 % 10;
				if (kdq % 10 == 1)
				{
					idqspn = 0;
				}
				else
				{
					idqspn = 1;
				}

				if (kk > ki)
				{
					ktemp = kk;
					kk = kj;
					kj = ki;
					ki = ktemp;
				}
				else if (kk > kj)
				{
					ktemp = kk;
					kk = kj;
					kj = ktemp;
				}

				if (ki != kj && ki != kk && kj != kk)
				{
					if (idqspn == 0)
					{
						kf = ki * 1000 + kk * 100 + kj * 10 + ibs;
					}
					else
					{
						kf = ki * 1000 + kj * 100 + kk * 10 + ibs;
					}
				}
				else if (ki == kj && ki == kk)
				{

					kf = ki * 1000 + kj * 100 + kk * 10 + 4;
				}
				else
				{
					kf = ki * 1000 + kj * 100 + kk * 10 + ibs;
				}

				if (kf == 2112 || kf == 2212)
				{
					i__2 = kf + 2;
					if ((r__1 = (float)xmpair - ulmass_(&kf), dabs(r__1)) > (r__2 = (float)xmpair - ulmass_(&i__2), dabs(r__2)))
					{
						kf += 2;
					}
				}
				if (k1 < 0)
				{
					kf = -kf;
				}
			}
			else if ((anim_1.isoft == 4 || anim_1.isoft == 5) &&
					 soft_1.njsgs[isg - 1] == 3)
			{
				if (k1abs > k2abs)
				{
					ki = k1abs;
					kk = k2abs;
				}
				else
				{
					ki = k2abs;
					kk = k1abs;
				}
				if (k3abs > ki)
				{
					kj = ki;
					ki = k3abs;
				}
				else if (k3abs < kk)
				{
					kj = kk;
					kk = k3abs;
				}
				else
				{
					kj = k3abs;
				}

				if (ki == kj && ki == kk)
				{

					ibs = 4;
					kf = ki * 1000 + kj * 100 + kk * 10 + ibs;
				}
				else if (ki != kj && ki != kk && kj != kk)
				{

					ibs = 2;
					kf1 = ki * 1000 + kj * 100 + kk * 10 + ibs;
					kf2 = ki * 1000 + kk * 100 + kj * 10 + ibs;
					kf = kf1;
					if ((r__1 = (float)xmpair - ulmass_(&kf1), dabs(r__1)) > (r__2 = (float)xmpair - ulmass_(&kf2), dabs(r__2)))
					{
						kf = kf2;
					}
				}
				else
				{
					ibs = 2;
					kf = ki * 1000 + kj * 100 + kk * 10 + ibs;

					if (kf == 2112 || kf == 2212)
					{
						i__2 = kf + 2;
						if ((r__1 = (float)xmpair - ulmass_(&kf), dabs(r__1)) > (r__2 = (float)xmpair - ulmass_(&i__2),
																				dabs(r__2)))
						{
							kf += 2;
						}
					}
				}
				if (k1 < 0)
				{
					kf = -kf;
				}
			}
			else
			{
				if (k1abs == k2abs)
				{
					if (k1abs <= 2)
					{

						kf = 0;
					}
					else if (k1abs <= 3)
					{

						kf = 333;
					}
					else
					{
						kf = k1abs * 100 + k1abs * 10 + (imspin << 1) + 1;
					}
				}
				else
				{
					if (k1abs > k2abs)
					{
						kmax = k1abs;
						kmin = k2abs;
					}
					else if (k1abs < k2abs)
					{
						kmax = k2abs;
						kmin = k1abs;
					}
					i__2 = k1 + k2;
					kf = (kmax * 100 + kmin * 10 + (imspin << 1) + 1) *
						 i_sign(&c__1, &i__2) * pow_ii(&c_n1, &kmax);

					if (abs(kf) % 10 == 1)
					{
						i__2 = abs(kf);
						i__3 = abs(kf) + 2;
						if ((r__1 = (float)xmpair - ulmass_(&i__2), dabs(r__1)) > (r__2 = (float)xmpair - ulmass_(&i__3),
																				  dabs(r__2)))
						{
							kf = (abs(kf) + 2) * i_sign(&c__1, &kf);
						}
					}
				}
			}
			arprc_1.itypar[hmain1_1.natt - 1] = kf;
			hmain2_1.katt[hmain1_1.natt - 1] = kf;
			if (abs(kf) == 211)
			{
				++npich;
			}
			else if (abs(kf) == 213)
			{
				++nrhoch;
			}
		}
	}

	if (nuudd != 0)
	{
		ppi0 = (float)(npich / 2) / (float)nuudd;
		prho0 = (float)(nrhoch / 2) / (float)nuudd;
	}

	npi0 = 0;
	i__1 = hjjet2_1.nsg;
	for (isg = 1; isg <= i__1; ++isg)
	{
		if (soft_1.k2sgs[isg - 1] == -soft_1.k2sgs[isg + 150000] && (i__2 = soft_1.k2sgs[isg - 1], abs(i__2)) <= 2 && soft_1.njsgs[isg - 1] == 2)
		{
			if (ranart_(&rndf77_1.nseed) <= ppi0)
			{
				++npi0;
			}
		}
	}

	if (nuudd > 1)
	{
		index1_(&c_b205, &nuudd, xmdiag, indx);
	}
	else
	{
		indx[0] = 1;
	}

	i__1 = nuudd;
	for (ix = 1; ix <= i__1; ++ix)
	{
		iuudd = indx[ix - 1];
		inatt = ndiag[iuudd - 1];
		if (ix <= npi0)
		{
			kf = 111;
		}
		else if (ranart_(&rndf77_1.nseed) <= prho0 / (1 - ppi0 + 1e-5f))
		{
			kf = 113;
		}
		else
		{

			if (ranart_(&rndf77_1.nseed) <= .5f)
			{
				kf = 221;
			}
			else
			{
				kf = 223;
			}
		}
		arprc_1.itypar[inatt - 1] = kf;
		hmain2_1.katt[inatt - 1] = kf;
	}

	inatt = 0;
	i__1 = hjjet2_1.nsg;
	for (isg = 1; isg <= i__1; ++isg)
	{
		if (soft_1.njsgs[isg - 1] != 0)
		{
			++inatt;
			k1 = soft_1.k2sgs[isg - 1];
			k1abs = abs(k1);
			px1 = soft_1.pxsgs[isg - 1];
			py1 = soft_1.pysgs[isg - 1];
			pz1 = soft_1.pzsgs[isg - 1];
			k2 = soft_1.k2sgs[isg + 150000];
			k2abs = abs(k2);
			px2 = soft_1.pxsgs[isg + 150000];
			py2 = soft_1.pysgs[isg + 150000];
			pz2 = soft_1.pzsgs[isg + 150000];
			e1 = soft_1.pesgs[isg - 1];
			e2 = soft_1.pesgs[isg + 150000];

			if (soft_1.njsgs[isg - 1] == 2)
			{
				arprc_1.pxar[inatt - 1] = (float)(px1 + px2);
				arprc_1.pyar[inatt - 1] = (float)(py1 + py2);
				arprc_1.pzar[inatt - 1] = (float)(pz1 + pz2);
				hmain2_1.patt[inatt - 1] = arprc_1.pxar[inatt - 1];
				hmain2_1.patt[inatt + 150000] = arprc_1.pyar[inatt - 1];
				hmain2_1.patt[inatt + 300001] = arprc_1.pzar[inatt - 1];
				etot = e1 + e2;
			}
			else if ((anim_1.isoft == 4 || anim_1.isoft == 5) &&
					 soft_1.njsgs[isg - 1] == 3)
			{
				px3 = soft_1.pxsgs[isg + 300001];
				py3 = soft_1.pysgs[isg + 300001];
				pz3 = soft_1.pzsgs[isg + 300001];
				e3 = soft_1.pesgs[isg + 300001];
				arprc_1.pxar[inatt - 1] = (float)(px1 + px2 + px3);
				arprc_1.pyar[inatt - 1] = (float)(py1 + py2 + py3);
				arprc_1.pzar[inatt - 1] = (float)(pz1 + pz2 + pz3);
				hmain2_1.patt[inatt - 1] = arprc_1.pxar[inatt - 1];
				hmain2_1.patt[inatt + 150000] = arprc_1.pyar[inatt - 1];
				hmain2_1.patt[inatt + 300001] = arprc_1.pzar[inatt - 1];
				etot = e1 + e2 + e3;
			}
			arprc_1.xmar[inatt - 1] = ulmass_(&arprc_1.itypar[inatt - 1]);

			r__1 = arprc_1.pxar[inatt - 1];

			r__2 = arprc_1.pyar[inatt - 1];

			r__3 = arprc_1.pzar[inatt - 1];

			r__4 = arprc_1.xmar[inatt - 1];
			arprc_1.pear[inatt - 1] = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
			hmain2_1.patt[inatt + 450002] = arprc_1.pear[inatt - 1];
			hmain1_1.eatt += arprc_1.pear[inatt - 1];
			ipartn = soft_1.njsgs[isg - 1];
			i__2 = ipartn;
			for (i__ = 1; i__ <= i__2; ++i__)
			{
				loclco_1.ftp[i__ - 1] = soft_1.ftsgs[isg + i__ * 150001 -
													 150002];
				loclco_1.gxp[i__ - 1] = soft_1.gxsgs[isg + i__ * 150001 -
													 150002];
				loclco_1.gyp[i__ - 1] = soft_1.gysgs[isg + i__ * 150001 -
													 150002];
				loclco_1.gzp[i__ - 1] = soft_1.gzsgs[isg + i__ * 150001 -
													 150002];
				loclco_1.pxp[i__ - 1] = soft_1.pxsgs[isg + i__ * 150001 -
													 150002];
				loclco_1.pyp[i__ - 1] = soft_1.pysgs[isg + i__ * 150001 -
													 150002];
				loclco_1.pzp[i__ - 1] = soft_1.pzsgs[isg + i__ * 150001 -
													 150002];
				loclco_1.pmp[i__ - 1] = soft_1.pmsgs[isg + i__ * 150001 -
													 150002];
				loclco_1.pep[i__ - 1] = soft_1.pesgs[isg + i__ * 150001 -
													 150002];
			}
			locldr_(&ipartn, &drlocl);

			tau0 = arprnt_1.arpar1[0];
			ftavg0 = prtn23_1.ft0fom + (double)tau0;
			gxavg0 = 0.;
			gyavg0 = 0.;
			gzavg0 = 0.;
			i__2 = ipartn;
			for (i__ = 1; i__ <= i__2; ++i__)
			{
				gxavg0 += prtn23_1.gxp0[i__ - 1] / ipartn;
				gyavg0 += prtn23_1.gyp0[i__ - 1] / ipartn;
				gzavg0 += prtn23_1.gzp0[i__ - 1] / ipartn;
			}
			bex = (double)arprc_1.pxar[inatt - 1] / etot;
			bey = (double)arprc_1.pyar[inatt - 1] / etot;
			bez = (double)arprc_1.pzar[inatt - 1] / etot;

			d__1 = bex;

			d__2 = bey;

			d__3 = bez;
			beta2 = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
			gam = 1. / sqrt(1. - beta2);
			if (beta2 >= .9999999999999)
			{
				s_wsle(&io___131);
				do_lio(&c__9, &c__1, "2", (ftnlen)1);
				do_lio(&c__5, &c__1, (char *)&bex, (ftnlen)sizeof(double));
				do_lio(&c__5, &c__1, (char *)&bey, (ftnlen)sizeof(double));
				do_lio(&c__5, &c__1, (char *)&bez, (ftnlen)sizeof(double));
				do_lio(&c__5, &c__1, (char *)&beta2, (ftnlen)sizeof(double));
				do_lio(&c__5, &c__1, (char *)&gam, (ftnlen)sizeof(double));
				e_wsle();
			}

			d__1 = -bex;
			d__2 = -bey;
			d__3 = -bez;
			lorenz_(&ftavg0, &gxavg0, &gyavg0, &gzavg0, &d__1, &d__2, &d__3);
			arprc_1.gxar[inatt - 1] = (float)lor_1.pxnew;
			arprc_1.gyar[inatt - 1] = (float)lor_1.pynew;
			arprc_1.gzar[inatt - 1] = (float)lor_1.pznew;
			arprc_1.ftar[inatt - 1] = (float)lor_1.enenew;
		}
	}

	nzpc_1.nattzp = hmain1_1.natt;
	ludat1_1.mstj[23] = mstj24;

	return 0;
}

int coales_(void)
{

	int i__1, i__2, i__3;

	int s_wsle(cilist *), do_lio(int *, int *, char *, ftnlen),
		e_wsle(void);
	int s_stop(char *, ftnlen);
	double sqrt(double);

	static int j, ip;
	static double dp0, dp1[2], dr0, dr1[2];
	static int ipi, isg, jsg, ipmin, ipmax, iover[150001];
	extern int exchge_(int *, int *, int *,
					   int *);
	static double dplocl;
	extern int locldr_(int *, double *);
	static double drlocl;
	static int ibaryn;

	static cilist io___134 = {0, 6, 0, 0, 0};

	i__1 = hjjet2_1.nsg;
	for (isg = 1; isg <= i__1; ++isg)
	{
		iover[isg - 1] = 0;
	}

	i__1 = hjjet2_1.nsg;
	for (isg = 1; isg <= i__1; ++isg)
	{
		if (soft_1.njsgs[isg - 1] != 2 || iover[isg - 1] == 1)
		{
			goto L150;
		}

		if (soft_1.k2sgs[isg - 1] < 0)
		{
			s_wsle(&io___134);
			do_lio(&c__9, &c__1, "Antiquark appears in quark loop; stop", (ftnlen)37);
			e_wsle();
			s_stop("", (ftnlen)0);
		}

		for (j = 1; j <= 2; ++j)
		{
			loclco_1.ftp[j - 1] = soft_1.ftsgs[isg + j * 150001 - 150002];
			loclco_1.gxp[j - 1] = soft_1.gxsgs[isg + j * 150001 - 150002];
			loclco_1.gyp[j - 1] = soft_1.gysgs[isg + j * 150001 - 150002];
			loclco_1.gzp[j - 1] = soft_1.gzsgs[isg + j * 150001 - 150002];
			loclco_1.pxp[j - 1] = soft_1.pxsgs[isg + j * 150001 - 150002];
			loclco_1.pyp[j - 1] = soft_1.pysgs[isg + j * 150001 - 150002];
			loclco_1.pzp[j - 1] = soft_1.pzsgs[isg + j * 150001 - 150002];
			loclco_1.pmp[j - 1] = soft_1.pmsgs[isg + j * 150001 - 150002];
			loclco_1.pep[j - 1] = soft_1.pesgs[isg + j * 150001 - 150002];
		}
		locldr_(&c__2, &drlocl);
		dr0 = drlocl;

		dp0 = sqrt((loclco_1.pep[0] * loclco_1.pep[1] - loclco_1.pxp[0] * loclco_1.pxp[1] - loclco_1.pyp[0] * loclco_1.pyp[1] -
					loclco_1.pzp[0] * loclco_1.pzp[1] - loclco_1.pmp[0] * loclco_1.pmp[1]) *
				   2);

		i__2 = hjjet2_1.nsg;
		for (jsg = 1; jsg <= i__2; ++jsg)
		{

			if (jsg == isg || iover[jsg - 1] == 1)
			{
				goto L120;
			}
			if (soft_1.njsgs[jsg - 1] == 2)
			{
				ipmin = 2;
				ipmax = 2;
			}
			else if (soft_1.njsgs[jsg - 1] == 3 && soft_1.k2sgs[jsg - 1] <
													   0)
			{
				ipmin = 1;
				ipmax = 3;
			}
			else
			{
				goto L120;
			}
			i__3 = ipmax;
			for (ip = ipmin; ip <= i__3; ++ip)
			{
				dplocl = sqrt((loclco_1.pep[0] * soft_1.pesgs[jsg + ip * 150001 - 150002] - loclco_1.pxp[0] * soft_1.pxsgs[jsg + ip * 150001 - 150002] - loclco_1.pyp[0] * soft_1.pysgs[jsg + ip * 150001 - 150002] -
							   loclco_1.pzp[0] * soft_1.pzsgs[jsg + ip * 150001 -
															  150002] -
							   loclco_1.pmp[0] * soft_1.pmsgs[jsg + ip * 150001 - 150002]) *
							  2);

				if (dplocl > coal_1.dpcoal)
				{
					goto L120;
				}
				loclco_1.ftp[1] = soft_1.ftsgs[jsg + ip * 150001 - 150002];
				loclco_1.gxp[1] = soft_1.gxsgs[jsg + ip * 150001 - 150002];
				loclco_1.gyp[1] = soft_1.gysgs[jsg + ip * 150001 - 150002];
				loclco_1.gzp[1] = soft_1.gzsgs[jsg + ip * 150001 - 150002];
				loclco_1.pxp[1] = soft_1.pxsgs[jsg + ip * 150001 - 150002];
				loclco_1.pyp[1] = soft_1.pysgs[jsg + ip * 150001 - 150002];
				loclco_1.pzp[1] = soft_1.pzsgs[jsg + ip * 150001 - 150002];
				loclco_1.pmp[1] = soft_1.pmsgs[jsg + ip * 150001 - 150002];
				loclco_1.pep[1] = soft_1.pesgs[jsg + ip * 150001 - 150002];
				locldr_(&c__2, &drlocl);

				if (drlocl > coal_1.drcoal)
				{
					goto L120;
				}

				if (dp0 > coal_1.dpcoal || dr0 > coal_1.drcoal || drlocl < dr0)
				{
					dp0 = dplocl;
					dr0 = drlocl;
					exchge_(&isg, &c__2, &jsg, &ip);
				}
			}
		L120:;
		}
		if (dp0 <= coal_1.dpcoal && dr0 <= coal_1.drcoal)
		{
			iover[isg - 1] = 1;
		}
	L150:;
	}

	i__1 = hjjet2_1.nsg;
	for (isg = 1; isg <= i__1; ++isg)
	{
		if (soft_1.njsgs[isg - 1] != 2 || iover[isg - 1] == 1)
		{
			goto L250;
		}

		for (j = 1; j <= 2; ++j)
		{
			loclco_1.ftp[j - 1] = soft_1.ftsgs[isg + j * 150001 - 150002];
			loclco_1.gxp[j - 1] = soft_1.gxsgs[isg + j * 150001 - 150002];
			loclco_1.gyp[j - 1] = soft_1.gysgs[isg + j * 150001 - 150002];
			loclco_1.gzp[j - 1] = soft_1.gzsgs[isg + j * 150001 - 150002];
			loclco_1.pxp[j - 1] = soft_1.pxsgs[isg + j * 150001 - 150002];
			loclco_1.pyp[j - 1] = soft_1.pysgs[isg + j * 150001 - 150002];
			loclco_1.pzp[j - 1] = soft_1.pzsgs[isg + j * 150001 - 150002];
			loclco_1.pmp[j - 1] = soft_1.pmsgs[isg + j * 150001 - 150002];
			loclco_1.pep[j - 1] = soft_1.pesgs[isg + j * 150001 - 150002];
		}
		locldr_(&c__2, &drlocl);
		dr0 = drlocl;
		dp0 = sqrt((loclco_1.pep[0] * loclco_1.pep[1] - loclco_1.pxp[0] * loclco_1.pxp[1] - loclco_1.pyp[0] * loclco_1.pyp[1] -
					loclco_1.pzp[0] * loclco_1.pzp[1] - loclco_1.pmp[0] * loclco_1.pmp[1]) *
				   2);

		i__2 = hjjet2_1.nsg;
		for (jsg = 1; jsg <= i__2; ++jsg)
		{
			if (jsg == isg || iover[jsg - 1] == 1)
			{
				goto L220;
			}
			if (soft_1.njsgs[jsg - 1] == 2)
			{
				ipmin = 1;
				ipmax = 1;
			}
			else if (soft_1.njsgs[jsg - 1] == 3 && soft_1.k2sgs[jsg - 1] >
													   0)
			{
				ipmin = 1;
				ipmax = 3;
			}
			else
			{
				goto L220;
			}
			i__3 = ipmax;
			for (ip = ipmin; ip <= i__3; ++ip)
			{
				dplocl = sqrt((loclco_1.pep[1] * soft_1.pesgs[jsg + ip * 150001 - 150002] - loclco_1.pxp[1] * soft_1.pxsgs[jsg + ip * 150001 - 150002] - loclco_1.pyp[1] * soft_1.pysgs[jsg + ip * 150001 - 150002] -
							   loclco_1.pzp[1] * soft_1.pzsgs[jsg + ip * 150001 -
															  150002] -
							   loclco_1.pmp[1] * soft_1.pmsgs[jsg + ip * 150001 - 150002]) *
							  2);

				if (dplocl > coal_1.dpcoal)
				{
					goto L220;
				}
				loclco_1.ftp[0] = soft_1.ftsgs[jsg + ip * 150001 - 150002];
				loclco_1.gxp[0] = soft_1.gxsgs[jsg + ip * 150001 - 150002];
				loclco_1.gyp[0] = soft_1.gysgs[jsg + ip * 150001 - 150002];
				loclco_1.gzp[0] = soft_1.gzsgs[jsg + ip * 150001 - 150002];
				loclco_1.pxp[0] = soft_1.pxsgs[jsg + ip * 150001 - 150002];
				loclco_1.pyp[0] = soft_1.pysgs[jsg + ip * 150001 - 150002];
				loclco_1.pzp[0] = soft_1.pzsgs[jsg + ip * 150001 - 150002];
				loclco_1.pmp[0] = soft_1.pmsgs[jsg + ip * 150001 - 150002];
				loclco_1.pep[0] = soft_1.pesgs[jsg + ip * 150001 - 150002];
				locldr_(&c__2, &drlocl);

				if (drlocl > coal_1.drcoal)
				{
					goto L220;
				}

				if (dp0 > coal_1.dpcoal || dr0 > coal_1.drcoal || drlocl < dr0)
				{
					dp0 = dplocl;
					dr0 = drlocl;
					exchge_(&isg, &c__1, &jsg, &ip);
				}
			}
		L220:;
		}
		if (dp0 <= coal_1.dpcoal && dr0 <= coal_1.drcoal)
		{
			iover[isg - 1] = 1;
		}
	L250:;
	}

	i__1 = hjjet2_1.nsg;
	for (isg = 1; isg <= i__1; ++isg)
	{
		if (soft_1.njsgs[isg - 1] != 3 || iover[isg - 1] == 1)
		{
			goto L350;
		}
		ibaryn = soft_1.k2sgs[isg - 1];

		for (j = 1; j <= 2; ++j)
		{
			loclco_1.ftp[j - 1] = soft_1.ftsgs[isg + j * 150001 - 150002];
			loclco_1.gxp[j - 1] = soft_1.gxsgs[isg + j * 150001 - 150002];
			loclco_1.gyp[j - 1] = soft_1.gysgs[isg + j * 150001 - 150002];
			loclco_1.gzp[j - 1] = soft_1.gzsgs[isg + j * 150001 - 150002];
			loclco_1.pxp[j - 1] = soft_1.pxsgs[isg + j * 150001 - 150002];
			loclco_1.pyp[j - 1] = soft_1.pysgs[isg + j * 150001 - 150002];
			loclco_1.pzp[j - 1] = soft_1.pzsgs[isg + j * 150001 - 150002];
			loclco_1.pmp[j - 1] = soft_1.pmsgs[isg + j * 150001 - 150002];
			loclco_1.pep[j - 1] = soft_1.pesgs[isg + j * 150001 - 150002];
		}
		locldr_(&c__2, &drlocl);
		dr1[0] = drlocl;
		dp1[0] = sqrt((loclco_1.pep[0] * loclco_1.pep[1] - loclco_1.pxp[0] * loclco_1.pxp[1] - loclco_1.pyp[0] * loclco_1.pyp[1] -
					   loclco_1.pzp[0] * loclco_1.pzp[1] - loclco_1.pmp[0] * loclco_1.pmp[1]) *
					  2);

		loclco_1.ftp[1] = soft_1.ftsgs[isg + 300001];
		loclco_1.gxp[1] = soft_1.gxsgs[isg + 300001];
		loclco_1.gyp[1] = soft_1.gysgs[isg + 300001];
		loclco_1.gzp[1] = soft_1.gzsgs[isg + 300001];
		loclco_1.pxp[1] = soft_1.pxsgs[isg + 300001];
		loclco_1.pyp[1] = soft_1.pysgs[isg + 300001];
		loclco_1.pzp[1] = soft_1.pzsgs[isg + 300001];
		loclco_1.pmp[1] = soft_1.pmsgs[isg + 300001];
		loclco_1.pep[1] = soft_1.pesgs[isg + 300001];
		locldr_(&c__2, &drlocl);
		dr1[1] = drlocl;
		dp1[1] = sqrt((loclco_1.pep[0] * loclco_1.pep[1] - loclco_1.pxp[0] * loclco_1.pxp[1] - loclco_1.pyp[0] * loclco_1.pyp[1] -
					   loclco_1.pzp[0] * loclco_1.pzp[1] - loclco_1.pmp[0] * loclco_1.pmp[1]) *
					  2);

		i__2 = hjjet2_1.nsg;
		for (jsg = 1; jsg <= i__2; ++jsg)
		{
			if (jsg == isg || iover[jsg - 1] == 1)
			{
				goto L320;
			}
			if (soft_1.njsgs[jsg - 1] == 2)
			{
				if (ibaryn > 0)
				{
					ipmin = 1;
				}
				else
				{
					ipmin = 2;
				}
				ipmax = ipmin;
			}
			else if (soft_1.njsgs[jsg - 1] == 3 && ibaryn * soft_1.k2sgs[jsg - 1] > 0)
			{
				ipmin = 1;
				ipmax = 3;
			}
			else
			{
				goto L320;
			}
			i__3 = ipmax;
			for (ip = ipmin; ip <= i__3; ++ip)
			{
				dplocl = sqrt((loclco_1.pep[0] * soft_1.pesgs[jsg + ip * 150001 - 150002] - loclco_1.pxp[0] * soft_1.pxsgs[jsg + ip * 150001 - 150002] - loclco_1.pyp[0] * soft_1.pysgs[jsg + ip * 150001 - 150002] -
							   loclco_1.pzp[0] * soft_1.pzsgs[jsg + ip * 150001 -
															  150002] -
							   loclco_1.pmp[0] * soft_1.pmsgs[jsg + ip * 150001 - 150002]) *
							  2);

				if (dplocl > coal_1.dpcoal)
				{
					goto L320;
				}
				loclco_1.ftp[1] = soft_1.ftsgs[jsg + ip * 150001 - 150002];
				loclco_1.gxp[1] = soft_1.gxsgs[jsg + ip * 150001 - 150002];
				loclco_1.gyp[1] = soft_1.gysgs[jsg + ip * 150001 - 150002];
				loclco_1.gzp[1] = soft_1.gzsgs[jsg + ip * 150001 - 150002];
				loclco_1.pxp[1] = soft_1.pxsgs[jsg + ip * 150001 - 150002];
				loclco_1.pyp[1] = soft_1.pysgs[jsg + ip * 150001 - 150002];
				loclco_1.pzp[1] = soft_1.pzsgs[jsg + ip * 150001 - 150002];
				loclco_1.pmp[1] = soft_1.pmsgs[jsg + ip * 150001 - 150002];
				loclco_1.pep[1] = soft_1.pesgs[jsg + ip * 150001 - 150002];
				locldr_(&c__2, &drlocl);

				if (drlocl > coal_1.drcoal)
				{
					goto L320;
				}

				ipi = 0;
				if (dp1[0] > coal_1.dpcoal || dr1[0] > coal_1.drcoal)
				{
					ipi = 2;
					if ((dp1[1] > coal_1.dpcoal || dr1[1] > coal_1.drcoal) &&
						dr1[1] > dr1[0])
					{
						ipi = 3;
					}
				}
				else if (dp1[1] > coal_1.dpcoal || dr1[1] > coal_1.drcoal)
				{
					ipi = 3;
				}
				else if (dr1[0] < dr1[1])
				{
					if (drlocl < dr1[1])
					{
						ipi = 3;
					}
				}
				else if (dr1[1] <= dr1[0])
				{
					if (drlocl < dr1[0])
					{
						ipi = 2;
					}
				}
				if (ipi != 0)
				{
					dp1[ipi - 2] = dplocl;
					dr1[ipi - 2] = drlocl;
					exchge_(&isg, &ipi, &jsg, &ip);
				}
			}
		L320:;
		}
		if (dp1[0] <= coal_1.dpcoal && dr1[0] <= coal_1.drcoal && dp1[1] <= coal_1.dpcoal && dr1[1] <= coal_1.drcoal)
		{
			iover[isg - 1] = 1;
		}
	L350:;
	}

	return 0;
}

int exchge_(int *isg, int *ipi, int *jsg,
			int *ipj)
{
	static int k1, k2;
	static double pe, ft, pm, gx, gy, gz, px, py, pz;

	k1 = soft_1.k1sgs[*isg + *ipi * 150001 - 150002];
	k2 = soft_1.k2sgs[*isg + *ipi * 150001 - 150002];
	px = soft_1.pxsgs[*isg + *ipi * 150001 - 150002];
	py = soft_1.pysgs[*isg + *ipi * 150001 - 150002];
	pz = soft_1.pzsgs[*isg + *ipi * 150001 - 150002];
	pe = soft_1.pesgs[*isg + *ipi * 150001 - 150002];
	pm = soft_1.pmsgs[*isg + *ipi * 150001 - 150002];
	gx = soft_1.gxsgs[*isg + *ipi * 150001 - 150002];
	gy = soft_1.gysgs[*isg + *ipi * 150001 - 150002];
	gz = soft_1.gzsgs[*isg + *ipi * 150001 - 150002];
	ft = soft_1.ftsgs[*isg + *ipi * 150001 - 150002];
	soft_1.k1sgs[*isg + *ipi * 150001 - 150002] = soft_1.k1sgs[*jsg + *ipj * 150001 - 150002];
	soft_1.k2sgs[*isg + *ipi * 150001 - 150002] = soft_1.k2sgs[*jsg + *ipj * 150001 - 150002];
	soft_1.pxsgs[*isg + *ipi * 150001 - 150002] = soft_1.pxsgs[*jsg + *ipj * 150001 - 150002];
	soft_1.pysgs[*isg + *ipi * 150001 - 150002] = soft_1.pysgs[*jsg + *ipj * 150001 - 150002];
	soft_1.pzsgs[*isg + *ipi * 150001 - 150002] = soft_1.pzsgs[*jsg + *ipj * 150001 - 150002];
	soft_1.pesgs[*isg + *ipi * 150001 - 150002] = soft_1.pesgs[*jsg + *ipj * 150001 - 150002];
	soft_1.pmsgs[*isg + *ipi * 150001 - 150002] = soft_1.pmsgs[*jsg + *ipj * 150001 - 150002];
	soft_1.gxsgs[*isg + *ipi * 150001 - 150002] = soft_1.gxsgs[*jsg + *ipj * 150001 - 150002];
	soft_1.gysgs[*isg + *ipi * 150001 - 150002] = soft_1.gysgs[*jsg + *ipj * 150001 - 150002];
	soft_1.gzsgs[*isg + *ipi * 150001 - 150002] = soft_1.gzsgs[*jsg + *ipj * 150001 - 150002];
	soft_1.ftsgs[*isg + *ipi * 150001 - 150002] = soft_1.ftsgs[*jsg + *ipj * 150001 - 150002];
	soft_1.k1sgs[*jsg + *ipj * 150001 - 150002] = k1;
	soft_1.k2sgs[*jsg + *ipj * 150001 - 150002] = k2;
	soft_1.pxsgs[*jsg + *ipj * 150001 - 150002] = px;
	soft_1.pysgs[*jsg + *ipj * 150001 - 150002] = py;
	soft_1.pzsgs[*jsg + *ipj * 150001 - 150002] = pz;
	soft_1.pesgs[*jsg + *ipj * 150001 - 150002] = pe;
	soft_1.pmsgs[*jsg + *ipj * 150001 - 150002] = pm;
	soft_1.gxsgs[*jsg + *ipj * 150001 - 150002] = gx;
	soft_1.gysgs[*jsg + *ipj * 150001 - 150002] = gy;
	soft_1.gzsgs[*jsg + *ipj * 150001 - 150002] = gz;
	soft_1.ftsgs[*jsg + *ipj * 150001 - 150002] = ft;

	return 0;
}

int locldr_(int *icall, double *drlocl)
{

	int i__1, i__2;
	double d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8, d__9;

	double sqrt(double);
	int s_wsle(cilist *), do_lio(int *, int *, char *, ftnlen),
		e_wsle(void);

	static int j;
	static double dt0, gam, bex, bey, bez, pep0[3], ftp0[3], pxp0[3],
		pyp0[3], pzp0[3];
	static int imin, imax;
	static double etot, beta2;
	static int ilate, istep, iearly;
	extern int lorenz_(double *, double *,
					   double *, double *, double *, double *,
					   double *);

	static cilist io___166 = {0, 6, 0, 0, 0};
	static cilist io___167 = {0, 6, 0, 0, 0};
	static cilist io___168 = {0, 6, 0, 0, 0};
	static cilist io___177 = {0, 6, 0, 0, 0};

	if (*icall == 2)
	{
		etot = loclco_1.pep[0] + loclco_1.pep[1];
		bex = (loclco_1.pxp[0] + loclco_1.pxp[1]) / etot;
		bey = (loclco_1.pyp[0] + loclco_1.pyp[1]) / etot;
		bez = (loclco_1.pzp[0] + loclco_1.pzp[1]) / etot;

		for (j = 1; j <= 2; ++j)
		{

			d__1 = bex;

			d__2 = bey;

			d__3 = bez;
			beta2 = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
			gam = 1. / sqrt(1. - beta2);
			if (beta2 >= .9999999999999)
			{
				s_wsle(&io___166);
				do_lio(&c__9, &c__1, "4", (ftnlen)1);
				do_lio(&c__5, &c__1, (char *)&loclco_1.pxp[0], (ftnlen)sizeof(double));
				do_lio(&c__5, &c__1, (char *)&loclco_1.pxp[1], (ftnlen)sizeof(double));
				do_lio(&c__5, &c__1, (char *)&loclco_1.pyp[0], (ftnlen)sizeof(double));
				do_lio(&c__5, &c__1, (char *)&loclco_1.pyp[1], (ftnlen)sizeof(double));
				do_lio(&c__5, &c__1, (char *)&loclco_1.pzp[0], (ftnlen)sizeof(double));
				do_lio(&c__5, &c__1, (char *)&loclco_1.pzp[1], (ftnlen)sizeof(double));
				do_lio(&c__5, &c__1, (char *)&loclco_1.pep[0], (ftnlen)sizeof(double));
				do_lio(&c__5, &c__1, (char *)&loclco_1.pep[1], (ftnlen)sizeof(double));
				do_lio(&c__5, &c__1, (char *)&loclco_1.pmp[0], (ftnlen)sizeof(double));
				do_lio(&c__5, &c__1, (char *)&loclco_1.pmp[1], (ftnlen)sizeof(double));

				d__2 = loclco_1.pxp[0];

				d__3 = loclco_1.pyp[0];

				d__4 = loclco_1.pzp[0];

				d__5 = loclco_1.pmp[0];
				d__1 = sqrt(d__2 * d__2 + d__3 * d__3 + d__4 * d__4 + d__5 * d__5) / loclco_1.pep[0];
				do_lio(&c__5, &c__1, (char *)&d__1, (ftnlen)sizeof(double));

				d__7 = loclco_1.pxp[0];

				d__8 = loclco_1.pyp[0];

				d__9 = loclco_1.pzp[0];
				d__6 = sqrt(d__7 * d__7 + d__8 * d__8 + d__9 * d__9) /
					   loclco_1.pep[0];
				do_lio(&c__5, &c__1, (char *)&d__6, (ftnlen)sizeof(double));
				e_wsle();
				s_wsle(&io___167);
				do_lio(&c__9, &c__1, "4a", (ftnlen)2);
				d__1 = loclco_1.pxp[0] + loclco_1.pxp[1];
				do_lio(&c__5, &c__1, (char *)&d__1, (ftnlen)sizeof(double));
				d__2 = loclco_1.pyp[0] + loclco_1.pyp[1];
				do_lio(&c__5, &c__1, (char *)&d__2, (ftnlen)sizeof(double));
				d__3 = loclco_1.pzp[0] + loclco_1.pzp[1];
				do_lio(&c__5, &c__1, (char *)&d__3, (ftnlen)sizeof(double));
				do_lio(&c__5, &c__1, (char *)&etot, (ftnlen)sizeof(double));
				e_wsle();
				s_wsle(&io___168);
				do_lio(&c__9, &c__1, "4b", (ftnlen)2);
				do_lio(&c__5, &c__1, (char *)&bex, (ftnlen)sizeof(double));
				do_lio(&c__5, &c__1, (char *)&bey, (ftnlen)sizeof(double));
				do_lio(&c__5, &c__1, (char *)&bez, (ftnlen)sizeof(double));
				do_lio(&c__5, &c__1, (char *)&beta2, (ftnlen)sizeof(double));
				do_lio(&c__5, &c__1, (char *)&gam, (ftnlen)sizeof(double));
				e_wsle();
			}

			lorenz_(&loclco_1.ftp[j - 1], &loclco_1.gxp[j - 1], &loclco_1.gyp[j - 1], &loclco_1.gzp[j - 1], &bex, &bey, &bez);
			prtn23_1.gxp0[j - 1] = lor_1.pxnew;
			prtn23_1.gyp0[j - 1] = lor_1.pynew;
			prtn23_1.gzp0[j - 1] = lor_1.pznew;
			ftp0[j - 1] = lor_1.enenew;
			lorenz_(&loclco_1.pep[j - 1], &loclco_1.pxp[j - 1], &loclco_1.pyp[j - 1], &loclco_1.pzp[j - 1], &bex, &bey, &bez);
			pxp0[j - 1] = lor_1.pxnew;
			pyp0[j - 1] = lor_1.pynew;
			pzp0[j - 1] = lor_1.pznew;
			pep0[j - 1] = lor_1.enenew;
		}

		if (ftp0[0] >= ftp0[1])
		{
			ilate = 1;
			iearly = 2;
		}
		else
		{
			ilate = 2;
			iearly = 1;
		}
		prtn23_1.ft0fom = ftp0[ilate - 1];

		dt0 = ftp0[ilate - 1] - ftp0[iearly - 1];
		prtn23_1.gxp0[iearly - 1] += pxp0[iearly - 1] / pep0[iearly - 1] *
									 dt0;
		prtn23_1.gyp0[iearly - 1] += pyp0[iearly - 1] / pep0[iearly - 1] *
									 dt0;
		prtn23_1.gzp0[iearly - 1] += pzp0[iearly - 1] / pep0[iearly - 1] *
									 dt0;

		d__1 = prtn23_1.gxp0[ilate - 1] - prtn23_1.gxp0[iearly - 1];

		d__2 = prtn23_1.gyp0[ilate - 1] - prtn23_1.gyp0[iearly - 1];

		d__3 = prtn23_1.gzp0[ilate - 1] - prtn23_1.gzp0[iearly - 1];
		*drlocl = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
	}
	else if (*icall == 3)
	{
		etot = loclco_1.pep[0] + loclco_1.pep[1] + loclco_1.pep[2];
		bex = (loclco_1.pxp[0] + loclco_1.pxp[1] + loclco_1.pxp[2]) / etot;
		bey = (loclco_1.pyp[0] + loclco_1.pyp[1] + loclco_1.pyp[2]) / etot;
		bez = (loclco_1.pzp[0] + loclco_1.pzp[1] + loclco_1.pzp[2]) / etot;

		d__1 = bex;

		d__2 = bey;

		d__3 = bez;
		beta2 = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
		gam = 1. / sqrt(1. - beta2);
		if (beta2 >= .9999999999999)
		{
			s_wsle(&io___177);
			do_lio(&c__9, &c__1, "5", (ftnlen)1);
			do_lio(&c__5, &c__1, (char *)&bex, (ftnlen)sizeof(double));
			do_lio(&c__5, &c__1, (char *)&bey, (ftnlen)sizeof(double));
			do_lio(&c__5, &c__1, (char *)&bez, (ftnlen)sizeof(double));
			do_lio(&c__5, &c__1, (char *)&beta2, (ftnlen)sizeof(double));
			do_lio(&c__5, &c__1, (char *)&gam, (ftnlen)sizeof(double));
			e_wsle();
		}

		for (j = 1; j <= 3; ++j)
		{
			lorenz_(&loclco_1.ftp[j - 1], &loclco_1.gxp[j - 1], &loclco_1.gyp[j - 1], &loclco_1.gzp[j - 1], &bex, &bey, &bez);
			prtn23_1.gxp0[j - 1] = lor_1.pxnew;
			prtn23_1.gyp0[j - 1] = lor_1.pynew;
			prtn23_1.gzp0[j - 1] = lor_1.pznew;
			ftp0[j - 1] = lor_1.enenew;
			lorenz_(&loclco_1.pep[j - 1], &loclco_1.pxp[j - 1], &loclco_1.pyp[j - 1], &loclco_1.pzp[j - 1], &bex, &bey, &bez);
			pxp0[j - 1] = lor_1.pxnew;
			pyp0[j - 1] = lor_1.pynew;
			pzp0[j - 1] = lor_1.pznew;
			pep0[j - 1] = lor_1.enenew;
		}

		if (ftp0[0] > ftp0[1])
		{
			ilate = 1;
			if (ftp0[2] > ftp0[0])
			{
				ilate = 3;
			}
		}
		else
		{
			ilate = 2;
			if (ftp0[2] >= ftp0[1])
			{
				ilate = 3;
			}
		}
		prtn23_1.ft0fom = ftp0[ilate - 1];

		if (ilate == 1)
		{
			imin = 2;
			imax = 3;
			istep = 1;
		}
		else if (ilate == 2)
		{
			imin = 1;
			imax = 3;
			istep = 2;
		}
		else if (ilate == 3)
		{
			imin = 1;
			imax = 2;
			istep = 1;
		}

		i__1 = imax;
		i__2 = istep;
		for (iearly = imin; i__2 < 0 ? iearly >= i__1 : iearly <= i__1;
			 iearly += i__2)
		{
			dt0 = ftp0[ilate - 1] - ftp0[iearly - 1];
			prtn23_1.gxp0[iearly - 1] += pxp0[iearly - 1] / pep0[iearly - 1] *
										 dt0;
			prtn23_1.gyp0[iearly - 1] += pyp0[iearly - 1] / pep0[iearly - 1] *
										 dt0;
			prtn23_1.gzp0[iearly - 1] += pzp0[iearly - 1] / pep0[iearly - 1] *
										 dt0;
		}
	}

	return 0;
}

int hoscar_(void)
{

	static int nff = 0;

	static char fmt_101[] = "(a10)";
	static char fmt_111[] = "(a12)";
	static char fmt_102[] = "(a4,1x,a20,1x,\002(\002,i3,\002,\002,i3,\002)+"
							"(\002,i3,\002,\002,i3,\002)\002,2x,a4,2x,e10.4,2x,i8)";
	static char fmt_112[] = "(\002# Center-of-mass energy/nucleon-pair is"
							"\002,f12.3,\002GeV\002)";
	static char fmt_103[] = "(i10,2x,i10,2x,f8.3,2x,f8.3)";
	static char fmt_104[] = "(i10,2x,i10,2x,9(e12.6,2x))";

	int i__1, i__2;
	float r__1, r__2, r__3, r__4;

	int s_wsfe(cilist *), do_fio(int *, char *, ftnlen), e_wsfe(void);
	int s_copy(char *, char *, ftnlen, ftnlen);
	int s_cmp(char *, char *, ftnlen, ftnlen);
	double sqrt(double);

	static int i__;
	static float ene, phi, xmp, xmt;
	static char code[8];
	static float ebeam;
	static char reffra[8];
	static int ievent;
	extern int invflv_(int *);
	static int ntestp;

	static cilist io___182 = {0, 19, 0, fmt_101, 0};
	static cilist io___183 = {0, 19, 0, fmt_111, 0};
	static cilist io___190 = {0, 19, 0, fmt_102, 0};
	static cilist io___193 = {0, 19, 0, fmt_112, 0};
	static cilist io___194 = {0, 19, 0, fmt_103, 0};
	static cilist io___197 = {0, 19, 0, fmt_104, 0};

	if (nff == 0)
	{
		s_wsfe(&io___182);
		do_fio(&c__1, "OSCAR1997A", (ftnlen)10);
		e_wsfe();
		s_wsfe(&io___183);
		do_fio(&c__1, "final_id_p_x", (ftnlen)12);
		e_wsfe();
		s_copy(code, "AMPT", (ftnlen)8, (ftnlen)4);
		if (s_cmp(oscar2_1.frame, "CMS", (ftnlen)8, (ftnlen)3) == 0)
		{
			s_copy(reffra, "nncm", (ftnlen)8, (ftnlen)4);
			xmp = (oscar1_1.izp * .93828f + (oscar1_1.iap - oscar1_1.izp) *
												.939457f) /
				  oscar1_1.iap;
			xmt = (oscar1_1.izt * .93828f + (oscar1_1.iat - oscar1_1.izt) *
												.939457f) /
				  oscar1_1.iat;

			r__1 = snn_1.efrm;

			r__2 = xmp;

			r__3 = xmt;
			ebeam = (r__1 * r__1 - r__2 * r__2 - r__3 * r__3) / 2.f / xmt;
		}
		else if (s_cmp(oscar2_1.frame, "LAB", (ftnlen)8, (ftnlen)3) == 0)
		{
			s_copy(reffra, "lab", (ftnlen)8, (ftnlen)3);
			ebeam = snn_1.efrm;
		}
		else
		{
			s_copy(reffra, "unknown", (ftnlen)8, (ftnlen)7);
			ebeam = 0.f;
		}
		ntestp = 1;
		s_wsfe(&io___190);
		do_fio(&c__1, code, (ftnlen)8);
		do_fio(&c__1, oscar2_1.amptvn, (ftnlen)25);
		do_fio(&c__1, (char *)&oscar1_1.iap, (ftnlen)sizeof(int));
		do_fio(&c__1, (char *)&oscar1_1.izp, (ftnlen)sizeof(int));
		do_fio(&c__1, (char *)&oscar1_1.iat, (ftnlen)sizeof(int));
		do_fio(&c__1, (char *)&oscar1_1.izt, (ftnlen)sizeof(int));
		do_fio(&c__1, reffra, (ftnlen)8);
		do_fio(&c__1, (char *)&ebeam, (ftnlen)sizeof(float));
		do_fio(&c__1, (char *)&ntestp, (ftnlen)sizeof(int));
		e_wsfe();
		nff = 1;
		ievent = 1;
		phi = 0.f;
		if (s_cmp(oscar2_1.frame, "CMS", (ftnlen)8, (ftnlen)3) == 0)
		{
			s_wsfe(&io___193);
			do_fio(&c__1, (char *)&snn_1.efrm, (ftnlen)sizeof(float));
			e_wsfe();
		}
	}

	s_wsfe(&io___194);
	do_fio(&c__1, (char *)&ievent, (ftnlen)sizeof(int));
	do_fio(&c__1, (char *)&hbt_1.nlast, (ftnlen)sizeof(int));
	do_fio(&c__1, (char *)&lastt_1.bimp, (ftnlen)sizeof(float));
	do_fio(&c__1, (char *)&phi, (ftnlen)sizeof(float));
	e_wsfe();

	i__1 = hbt_1.nlast;
	for (i__ = 1; i__ <= i__1; ++i__)
	{

		r__1 = hbt_1.plast[(i__ << 2) - 4];

		r__2 = hbt_1.plast[(i__ << 2) - 3];

		r__3 = hbt_1.plast[(i__ << 2) - 2];

		r__4 = hbt_1.plast[(i__ << 2) - 1];
		ene = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
		s_wsfe(&io___197);
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(int));
		i__2 = invflv_(&hbt_1.lblast[i__ - 1]);
		do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(int));
		do_fio(&c__1, (char *)&hbt_1.plast[(i__ << 2) - 4], (ftnlen)sizeof(float));
		do_fio(&c__1, (char *)&hbt_1.plast[(i__ << 2) - 3], (ftnlen)sizeof(float));
		do_fio(&c__1, (char *)&hbt_1.plast[(i__ << 2) - 2], (ftnlen)sizeof(float));
		do_fio(&c__1, (char *)&ene, (ftnlen)sizeof(float));
		do_fio(&c__1, (char *)&hbt_1.plast[(i__ << 2) - 1], (ftnlen)sizeof(float));
		do_fio(&c__1, (char *)&hbt_1.xlast[(i__ << 2) - 4], (ftnlen)sizeof(float));
		do_fio(&c__1, (char *)&hbt_1.xlast[(i__ << 2) - 3], (ftnlen)sizeof(float));
		do_fio(&c__1, (char *)&hbt_1.xlast[(i__ << 2) - 2], (ftnlen)sizeof(float));
		do_fio(&c__1, (char *)&hbt_1.xlast[(i__ << 2) - 1], (ftnlen)sizeof(float));
		e_wsfe();
	}
	++ievent;

	return 0;
}

int getnp_(void)
{

	int i__1;
	float r__1, r__2;

	double sqrt(double);

	static int i__, nspec1, nspec2;
	static float epsiln, pztarg, pzproj;

	if (hmain1_1.natt == 0)
	{
		snn_1.npart1 = 0;
		snn_1.npart2 = 0;
		return 0;
	}

	r__1 = hparnt_1.hint1[5];

	r__2 = hparnt_1.hint1[7];
	pzproj = sqrt(r__1 * r__1 - r__2 * r__2);

	r__1 = hparnt_1.hint1[6];

	r__2 = hparnt_1.hint1[8];
	pztarg = sqrt(r__1 * r__1 - r__2 * r__2);
	epsiln = .01f;

	nspec1 = 0;
	nspec2 = 0;
	i__1 = hmain1_1.natt;
	for (i__ = 1; i__ <= i__1; ++i__)
	{
		if ((hmain2_1.katt[i__ - 1] == 2112 || hmain2_1.katt[i__ - 1] == 2212) && hmain2_1.patt[i__ - 1] == 0.f && hmain2_1.patt[i__ + 150000] == 0.f)
		{

			r__1 = 0.f, r__2 = pzproj - epsiln;
			if (hmain2_1.patt[i__ + 300001] > dmax(r__1, r__2))
			{
				++nspec1;
			}
			else if (hmain2_1.patt[i__ + 300001] < -pztarg + epsiln)
			{
				++nspec2;
			}
		}
	}
	snn_1.npart1 = hparnt_1.ihnt2[0] - nspec1;
	snn_1.npart2 = hparnt_1.ihnt2[2] - nspec2;
	return 0;
}

int resdec_(int *i1, int *nt, int *nnn, float *wid, int *idecay)
{

	int i__1;
	float r__1, r__2, r__3, r__4;

	int s_wsle(cilist *), do_lio(int *, int *, char *, ftnlen),
		e_wsle(void);
	double sqrt(double), log(double);
	int s_stop(char *, ftnlen);

	static int kf, ip;
	static float em1a, tau0;
	static int idau;
	static float enet;
	static int irun;
	static float esave;
	static int kdaut, ndaut, ksave;
	static float dpdecp;
	static int lbdaut;
	extern int iarflv_(int *);
	extern int ludecy_(int *);
	static float taudcy;
	extern double ranart_(int *);
	static float xmsave;
	extern int invflv_(int *);
	static float pxsave, pysave, pzsave;
	extern int lulist_(int *);

	static cilist io___208 = {0, 6, 0, 0, 0};
	static cilist io___213 = {0, 10, 0, 0, 0};
	static cilist io___214 = {0, 89, 0, 0, 0};

	irun = *idecay;
	if (leadng_1.lb1 == 0 || leadng_1.lb1 == 25 || leadng_1.lb1 == 26 ||
		leadng_1.lb1 == 27 || leadng_1.lb1 == 28 || leadng_1.lb1 == 29 ||
		abs(leadng_1.lb1) == 30 || leadng_1.lb1 == 24 || abs(leadng_1.lb1) >= 6 && abs(leadng_1.lb1) <= 9 || abs(leadng_1.lb1) == 16)
	{
		kf = invflv_(&leadng_1.lb1);
	}
	else
	{
		return 0;
	}

	ip = 1;

	lujets_1.n = 1;
	lujets_1.k[ip - 1] = 1;
	lujets_1.k[ip + 17999] = 0;
	lujets_1.k[ip + 26999] = 0;
	lujets_1.k[ip + 35999] = 0;

	lujets_1.k[ip + 8999] = kf;
	lujets_1.p[ip - 1] = leadng_1.px1;
	lujets_1.p[ip + 8999] = leadng_1.py1;
	lujets_1.p[ip + 17999] = leadng_1.pz1;
	em1a = leadng_1.em1;

	if ((leadng_1.lb1 == 0 || leadng_1.lb1 == 28) && leadng_1.em1 <
														 .43500000000000005f)
	{
		leadng_1.em1 = .43500000000000005f;
	}
	else if (leadng_1.lb1 >= 25 && leadng_1.lb1 <= 27 && leadng_1.em1 < .30000000000000004f)
	{
		leadng_1.em1 = .30000000000000004f;
	}
	else if (abs(leadng_1.lb1) == 30 && leadng_1.em1 < .65800000000000003f)
	{
		leadng_1.em1 = .65800000000000003f;
	}
	else if (abs(leadng_1.lb1) >= 6 && abs(leadng_1.lb1) <= 9 &&
			 leadng_1.em1 < 1.1000000000000001f)
	{
		leadng_1.em1 = 1.1000000000000001f;
	}
	if (leadng_1.em1 >= em1a + .01f)
	{
		s_wsle(&io___208);
		do_lio(&c__9, &c__1, "Mass increase in resdec():", (ftnlen)26);
		do_lio(&c__3, &c__1, (char *)&(*nt), (ftnlen)sizeof(int));
		r__1 = leadng_1.em1 - em1a;
		do_lio(&c__4, &c__1, (char *)&r__1, (ftnlen)sizeof(float));
		do_lio(&c__3, &c__1, (char *)&leadng_1.lb1, (ftnlen)sizeof(int));
		e_wsle();
	}

	r__1 = leadng_1.em1;

	r__2 = leadng_1.px1;

	r__3 = leadng_1.py1;

	r__4 = leadng_1.pz1;
	leadng_1.e1 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3 + r__4 * r__4);
	lujets_1.p[ip + 26999] = leadng_1.e1;
	lujets_1.p[ip + 35999] = leadng_1.em1;

	dpdecp = dpert_1.dpertp[*i1 - 1];
	ludecy_(&ip);

	if (*nt == input2_1.ntmax)
	{
		tau0 = .19733f / *wid;
		taudcy = tau0 * -1.f * log(1.f - ranart_(&rndf77_1.nseed));
		ndaut = lujets_1.n - resdcy_1.nsav;
		if (ndaut <= 1)
		{
			s_wsle(&io___213);
			do_lio(&c__9, &c__1, "note: ndaut(<1)=", (ftnlen)16);
			do_lio(&c__3, &c__1, (char *)&ndaut, (ftnlen)sizeof(int));
			e_wsle();
			s_wsle(&io___214);
			do_lio(&c__9, &c__1, "note: ndaut(<1)=", (ftnlen)16);
			do_lio(&c__3, &c__1, (char *)&ndaut, (ftnlen)sizeof(int));
			e_wsle();
			lulist_(&c__2);
			s_stop("", (ftnlen)0);
		}

		taudcy = taudcy * leadng_1.e1 / leadng_1.em1;
		leadng_1.tfnl += taudcy;
		leadng_1.xfnl += leadng_1.px1 / leadng_1.e1 * taudcy;
		leadng_1.yfnl += leadng_1.py1 / leadng_1.e1 * taudcy;
		leadng_1.zfnl += leadng_1.pz1 / leadng_1.e1 * taudcy;

		if (lujets_1.n >= resdcy_1.nsav + 2)
		{
			i__1 = lujets_1.n;
			for (idau = resdcy_1.nsav + 2; idau <= i__1; ++idau)
			{
				kdaut = lujets_1.k[idau + 8999];
				if (kdaut == 221 || kdaut == 113 || kdaut == 213 || kdaut == -213 || kdaut == 310)
				{

					ksave = kdaut;
					pxsave = lujets_1.p[idau - 1];
					pysave = lujets_1.p[idau + 8999];
					pzsave = lujets_1.p[idau + 17999];
					esave = lujets_1.p[idau + 26999];
					xmsave = lujets_1.p[idau + 35999];
					lujets_1.k[idau + 8999] = lujets_1.k[resdcy_1.nsav + 9000];
					lujets_1.p[idau - 1] = lujets_1.p[resdcy_1.nsav];
					lujets_1.p[idau + 8999] = lujets_1.p[resdcy_1.nsav + 9000];
					lujets_1.p[idau + 17999] = lujets_1.p[resdcy_1.nsav +
														  18000];
					lujets_1.p[idau + 26999] = lujets_1.p[resdcy_1.nsav +
														  27000];
					lujets_1.p[idau + 35999] = lujets_1.p[resdcy_1.nsav +
														  36000];
					lujets_1.k[resdcy_1.nsav + 9000] = ksave;
					lujets_1.p[resdcy_1.nsav] = pxsave;
					lujets_1.p[resdcy_1.nsav + 9000] = pysave;
					lujets_1.p[resdcy_1.nsav + 18000] = pzsave;
					lujets_1.p[resdcy_1.nsav + 27000] = esave;
					lujets_1.p[resdcy_1.nsav + 36000] = xmsave;

					goto L111;
				}
			}
		}
	L111:

		enet = 0.f;
		i__1 = lujets_1.n;
		for (idau = resdcy_1.nsav + 1; idau <= i__1; ++idau)
		{
			enet += lujets_1.p[idau + 26999];
		}
	}

	i__1 = lujets_1.n;
	for (idau = resdcy_1.nsav + 1; idau <= i__1; ++idau)
	{
		kdaut = lujets_1.k[idau + 8999];
		lbdaut = iarflv_(&kdaut);

		if (*nt == input2_1.ntmax && (kdaut == 130 || kdaut == 310 || abs(kdaut) == 311))
		{
			if (kdaut == 130)
			{
				lbdaut = 22;
			}
			else if (kdaut == 310)
			{
				lbdaut = 24;
			}
			else if (abs(kdaut) == 311)
			{
				if (ranart_(&rndf77_1.nseed) < .5f)
				{
					lbdaut = 22;
				}
				else
				{
					lbdaut = 24;
				}
			}
		}

		if (idau == resdcy_1.nsav + 1)
		{
			ee_1.lb[*i1 - 1] = lbdaut;
			cc_1.e[*i1 - 1] = lujets_1.p[idau + 35999];
			leadng_1.px1n = lujets_1.p[idau - 1];
			leadng_1.py1n = lujets_1.p[idau + 8999];
			leadng_1.pz1n = lujets_1.p[idau + 17999];

			leadng_1.dp1n = dpdecp;
		}
		else
		{
			++(*nnn);
			pd_1.lpion[*nnn + irun * 150001 - 150002] = lbdaut;
			pc_1.epion[*nnn + irun * 150001 - 150002] = lujets_1.p[idau +
																   35999];
			pb_1.ppion[(*nnn + irun * 150001) * 3 - 450006] = lujets_1.p[idau - 1];
			pb_1.ppion[(*nnn + irun * 150001) * 3 - 450005] = lujets_1.p[idau + 8999];
			pb_1.ppion[(*nnn + irun * 150001) * 3 - 450004] = lujets_1.p[idau + 17999];
			pa_1.rpion[(*nnn + irun * 150001) * 3 - 450006] = leadng_1.xfnl;
			pa_1.rpion[(*nnn + irun * 150001) * 3 - 450005] = leadng_1.yfnl;
			pa_1.rpion[(*nnn + irun * 150001) * 3 - 450004] = leadng_1.zfnl;
			tdecay_1.tfdpi[*nnn + irun * 150001 - 150002] = leadng_1.tfnl;

			dpert_1.dppion[*nnn + irun * 150001 - 150002] = dpdecp;
		}
	}

	return 0;
}

int inidcy_(void)
{

	lujets_1.n = 1;
	resdcy_1.nsav = lujets_1.n;
	return 0;
}

int local_(double *t)
{

	int i__1, i__2;
	double d__1, d__2, d__3, d__4;

	int s_wsle(cilist *), do_lio(int *, int *, char *, ftnlen),
		e_wsle(void);
	int s_stop(char *, ftnlen);
	double sqrt(double), cosh(double);

	static double x0, y0;
	static int ip, it;
	static double drt, eta0, rap0, detdy;
	static int itest;
	static double xtest, ytest, etcrit, ettest;

	static cilist io___226 = {0, 1, 0, 0, 0};

	for (it = 1; it <= 301; ++it)
	{
		if (*t >= frzprc_1.tfrz[it - 1] && *t < frzprc_1.tfrz[it])
		{
			if (it == frzprc_1.itlast)
			{
				return 0;
			}
			else
			{
				frzprc_1.itlast = it;
				goto L50;
			}
		}
	}
	s_wsle(&io___226);
	do_lio(&c__9, &c__1, "local time out of range in LOCAL, stop", (ftnlen)38);
	do_lio(&c__5, &c__1, (char *)&(*t), (ftnlen)sizeof(double));
	do_lio(&c__3, &c__1, (char *)&it, (ftnlen)sizeof(int));
	e_wsle();
	s_stop("", (ftnlen)0);
L50:

	i__1 = para1_1.mul;
	for (ip = 1; ip <= i__1; ++ip)
	{

		if (frzprc_1.ifrz[ip - 1] == 1)
		{
			goto L200;
		}
		if (it == 301)
		{

			etcrit = 1e6;
			goto L150;
		}
		else
		{

			etcrit = coal_1.ecritl * 2. / 3.;
		}

		if (*t < prec2_1.ft5[ip - 1])
		{
			goto L200;
		}
		rap0 = prec5_1.rap[ip - 1];
		eta0 = prec5_1.eta[ip - 1];
		x0 = prec2_1.gx5[ip - 1] + prec4_1.vx[ip - 1] * (*t - prec2_1.ft5[ip - 1]);
		y0 = prec2_1.gy5[ip - 1] + prec4_1.vy[ip - 1] * (*t - prec2_1.ft5[ip - 1]);
		detdy = 0.;
		i__2 = para1_1.mul;
		for (itest = 1; itest <= i__2; ++itest)
		{

			if (itest == ip || *t < prec2_1.ft5[itest - 1])
			{
				goto L100;
			}
			ettest = prec5_1.eta[itest - 1];
			xtest = prec2_1.gx5[itest - 1] + prec4_1.vx[itest - 1] * (*t -
																	  prec2_1.ft5[itest - 1]);
			ytest = prec2_1.gy5[itest - 1] + prec4_1.vy[itest - 1] * (*t -
																	  prec2_1.ft5[itest - 1]);

			d__1 = xtest - x0;

			d__2 = ytest - y0;
			drt = sqrt(d__1 * d__1 + d__2 * d__2);

			if ((d__1 = ettest - eta0, abs(d__1)) <= 1. && drt <= 1.)
			{

				d__2 = prec2_1.px5[itest - 1];

				d__3 = prec2_1.py5[itest - 1];

				d__4 = prec2_1.xmass5[itest - 1];
				detdy += sqrt(d__2 * d__2 + d__3 * d__3 + d__4 * d__4) * .5;
			}
		L100:;
		}

		d__1 = cosh(eta0);
		detdy = detdy * (d__1 * d__1) / (*t * 3.1416 * 1. * cosh(rap0));

	L150:
		if (detdy <= etcrit)
		{
			frzprc_1.ifrz[ip - 1] = 1;
			frzprc_1.idfrz[ip - 1] = prec2_1.ityp5[ip - 1];
			frzprc_1.pxfrz[ip - 1] = prec2_1.px5[ip - 1];
			frzprc_1.pyfrz[ip - 1] = prec2_1.py5[ip - 1];
			frzprc_1.pzfrz[ip - 1] = prec2_1.pz5[ip - 1];
			frzprc_1.efrz[ip - 1] = prec2_1.e5[ip - 1];
			frzprc_1.xmfrz[ip - 1] = prec2_1.xmass5[ip - 1];
			if (*t > prec2_1.ft5[ip - 1])
			{
				frzprc_1.gxfrz[ip - 1] = x0;
				frzprc_1.gyfrz[ip - 1] = y0;
				frzprc_1.gzfrz[ip - 1] = prec2_1.gz5[ip - 1] + prec4_1.vz[ip - 1] * (*t - prec2_1.ft5[ip - 1]);
				frzprc_1.ftfrz[ip - 1] = *t;
			}
			else
			{

				frzprc_1.gxfrz[ip - 1] = prec2_1.gx5[ip - 1];
				frzprc_1.gyfrz[ip - 1] = prec2_1.gy5[ip - 1];
				frzprc_1.gzfrz[ip - 1] = prec2_1.gz5[ip - 1];
				frzprc_1.ftfrz[ip - 1] = prec2_1.ft5[ip - 1];
			}
		}
	L200:;
	}

	return 0;
}

int inifrz_(void)
{
	static int it;
	static double step1, step2, step3, step4;

	step1 = .1;
	step2 = 1.;
	step3 = 10.;
	step4 = 100.;

	for (it = 1; it <= 101; ++it)
	{
		frzprc_1.tfrz[it - 1] = (double)(it - 1) * step1 + 0.;
	}
	for (it = 102; it <= 191; ++it)
	{
		frzprc_1.tfrz[it - 1] = (double)(it - 101) * step2 + 10.;
	}
	for (it = 192; it <= 281; ++it)
	{
		frzprc_1.tfrz[it - 1] = (double)(it - 191) * step3 + 100.;
	}
	for (it = 282; it <= 301; ++it)
	{
		frzprc_1.tfrz[it - 1] = (double)(it - 281) * step4 + 1e3;
	}
	frzprc_1.tfrz[301] = ilist5_1.tlarge;

	return 0;
}
