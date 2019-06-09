#include "f2c.h"
#include <math.h>
struct {
    int natt;
    double eatt;
    int jatt, nt, np, n0, n01, n10, n11;
} hmain1_;

#define hmain1_1 hmain1_

struct {
    double hipr1[100];
    int ihpr2[50];
    double hint1[100];
    int ihnt2[50];
} hparnt_;

#define hparnt_1 hparnt_

struct {
    int mstu[200];
    double paru[200];
    int mstj[200];
    double parj[200];
} ludat1_;

#define ludat1_1 ludat1_

struct {
    double arpar1[100];
    int iapar2[50];
    double arint1[100];
    int iaint2[50];
} arprnt_;

#define arprnt_1 arprnt_

struct {
    int iout;
} arout_;

#define arout_1 arout_

struct {
    int iaevt, iarun;
} arevt_;

#define arevt_1 arevt_

struct {
    float smearp, smearh;
} smearz_;

#define smearz_1 smearz_

struct {
    int nseed;
} rndf77_;

#define rndf77_1 rndf77_

struct {
    int nevent, isoft, isflag, izpc;
} anim_;

#define anim_1 anim_

struct {
    float dpcoal, drcoal, ecritl;
} coal_;

#define coal_1 coal_

struct {
    double efrm;
    int npart1, npart2;
} snn_;

#define snn_1 snn_

struct {
    float xmp, xmu, alpha, rscut2, cutof2;
} para2_;

#define para2_1 para2_

struct {
    int ioscar;
} para7_;

#define para7_1 para7_

struct {
    int idpert, npertd, idxsec;
} para8_;

#define para8_1 para8_

struct {
    int iseedp;
} rndm3_;

#define rndm3_1 rndm3_

struct {
    int num;
} run_;

#define run_1 run_

struct {
    int masspr, massta, iseed, iavoid;
    double dt;
} input1_;

#define input1_1 input1_

struct {
    int ilab, manyb, ntmax, icoll, insys, ipot, mode, imomen, nfreq, 
	    icflow, icrho, icou, kpoten, kmul;
} input2_;

#define input2_1 input2_

struct {
    int iap, izp, iat, izt;
} oscar1_;

#define oscar1_1 oscar1_

struct {
    char frame[8], amptvn[25];
} oscar2_;

#define oscar2_1 oscar2_

struct {
    int nsav, iksdcy;
} resdcy_;

#define resdcy_1 resdcy_

struct {
    int iphidcy;
} phidcy_;

#define phidcy_1 phidcy_



static int c__4 = 4;
static int c__1 = 1;
static int c__3 = 3;
static int c__5 = 5;
static int c__9 = 9;
static float c_b92 = .3333;


 int MAIN__(void)
{
   
    static char fmt_111[] = "(a8)";
    static char fmt_50[] = "(\002 \002/11x,\002#############################"
	    "#####################\002/1x,10x,\002#      YaY (A Multi-Phase Transport) model      #\002/1x,10x,\002#               Version"
	    " \002,a20,\002     #\002/1x,10x,\002#                10/01/2008 "
	    "                     #\002/1x,10x,\002##########################"
	    "########################\002/1x,10x,\002 \002)";

   
    int i__1, i__2;
    float d__1;
    olist o__1;
    cllist cl__1;
    int f_open(olist *), s_rsle(cilist *), do_lio(int *, int *, char *, int), e_rsle(void), s_rsfe(cilist *), do_fio(int *, char *, int), e_rsfe(void), f_clos(cllist *);
    int s_copy(char *, char *, int, int);
    int s_wsfe(cilist *), e_wsfe(void), s_wsle(cilist *), e_wsle(void);
    float pow_dd(float *, float *);
    int s_stop(char *, int);
    static int j, k;
    static double bmin, bmax;
    static char targ[8], proj[8];
    static int ipop, iamax;
    extern int arini_(void), getnp_(void),artmn_(void);
    extern float rand(int *);     
    static int imiss, nevnt;
    extern int hjana3_(void), hjana4_(void), arini2_(int *), artan1_(void), artan2_(void);
    static int ihjsed;
    extern int hijing_(char *, double *, double *, int);
    static int nseedr;
    extern int hijset_(double *, char *, char *, char *,int *, int *, int *, int *, int, int,int), inizpc_(void), artset_(void), artout_(int *);
    static cilist io___1 = { 0, 24, 0, 0, 0 };
    static cilist io___2 = { 0, 24, 0, fmt_111, 0 };
    static cilist io___3 = { 0, 24, 0, fmt_111, 0 };
    static cilist io___5 = { 0, 24, 0, fmt_111, 0 };
    static cilist io___7 = { 0, 24, 0, 0, 0 };
    static cilist io___8 = { 0, 24, 0, 0, 0 };
    static cilist io___9 = { 0, 24, 0, 0, 0 };
    static cilist io___10 = { 0, 24, 0, 0, 0 };
    static cilist io___11 = { 0, 24, 0, 0, 0 };
    static cilist io___13 = { 0, 24, 0, 0, 0 };
    static cilist io___15 = { 0, 24, 0, 0, 0 };
    static cilist io___17 = { 0, 24, 0, 0, 0 };
    static cilist io___18 = { 0, 24, 0, 0, 0 };
    static cilist io___19 = { 0, 24, 0, 0, 0 };
    static cilist io___20 = { 0, 24, 0, 0, 0 };
    static cilist io___21 = { 0, 24, 0, 0, 0 };
    static cilist io___22 = { 0, 24, 0, 0, 0 };
    static cilist io___24 = { 0, 24, 0, 0, 0 };
    static cilist io___25 = { 0, 24, 0, 0, 0 };
    static cilist io___26 = { 0, 24, 0, 0, 0 };
    static cilist io___27 = { 0, 24, 0, 0, 0 };
    static cilist io___28 = { 0, 24, 0, 0, 0 };
    static cilist io___29 = { 0, 24, 0, 0, 0 };
    static cilist io___30 = { 0, 24, 0, 0, 0 };
    static cilist io___31 = { 0, 24, 0, 0, 0 };
    static cilist io___32 = { 0, 24, 0, 0, 0 };
    static cilist io___33 = { 0, 24, 0, 0, 0 };
    static cilist io___34 = { 0, 24, 0, 0, 0 };
    static cilist io___36 = { 0, 24, 0, 0, 0 };
    static cilist io___37 = { 0, 24, 0, 0, 0 };
    static cilist io___38 = { 0, 24, 0, 0, 0 };
    static cilist io___39 = { 0, 24, 0, 0, 0 };
    static cilist io___40 = { 0, 24, 0, 0, 0 };
    static cilist io___41 = { 0, 6, 0, fmt_50, 0 };
    static cilist io___42 = { 0, 12, 0, fmt_50, 0 };
    static cilist io___43 = { 0, 6, 0, 0, 0 };
    static cilist io___44 = { 0, 5, 0, 0, 0 };
    static cilist io___46 = { 0, 6, 0, 0, 0 };
    static cilist io___47 = { 0, 12, 0, 0, 0 };
    static cilist io___51 = { 0, 6, 0, 0, 0 };
    static cilist io___53 = { 0, 6, 0, 0, 0 };
    static cilist io___54 = { 0, 6, 0, 0, 0 };
    o__1.oerr = 0;
    o__1.ounit = 24;
    o__1.ofnmlen = 10;
    o__1.ofnm = "input.ampt";
    o__1.orl = 0;
    o__1.osta = "UNKNOWN";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    o__1.oerr = 0;
    o__1.ounit = 12;
    o__1.ofnmlen = 11;
    o__1.ofnm = "ana/version";
    o__1.orl = 0;
    o__1.osta = "UNKNOWN";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    s_rsle(&io___1);
    do_lio(&c__4, &c__1, (char *)&snn_1.efrm, (int)sizeof(double));
    e_rsle();

    s_rsfe(&io___2);
    do_fio(&c__1, oscar2_1.frame, (int)8);
    e_rsfe();
    s_rsfe(&io___3);
    do_fio(&c__1, proj, (int)8);
    e_rsfe();
    s_rsfe(&io___5);
    do_fio(&c__1, targ, (int)8);
    e_rsfe();
    s_rsle(&io___7);
    do_lio(&c__3, &c__1, (char *)&oscar1_1.iap, (int)sizeof(int));
    e_rsle();
    s_rsle(&io___8);
    do_lio(&c__3, &c__1, (char *)&oscar1_1.izp, (int)sizeof(int));
    e_rsle();
    s_rsle(&io___9);
    do_lio(&c__3, &c__1, (char *)&oscar1_1.iat, (int)sizeof(int));
    e_rsle();
    s_rsle(&io___10);
    do_lio(&c__3, &c__1, (char *)&oscar1_1.izt, (int)sizeof(int));
    e_rsle();
    s_rsle(&io___11);
    do_lio(&c__3, &c__1, (char *)&nevnt, (int)sizeof(int));
    e_rsle();
    s_rsle(&io___13);
    do_lio(&c__4, &c__1, (char *)&bmin, (int)sizeof(double));
    e_rsle();
    s_rsle(&io___15);
    do_lio(&c__4, &c__1, (char *)&bmax, (int)sizeof(double));
    e_rsle();

    s_rsle(&io___17);
    do_lio(&c__3, &c__1, (char *)&anim_1.isoft, (int)sizeof(int));
    e_rsle();

    s_rsle(&io___18);
    do_lio(&c__3, &c__1, (char *)&input2_1.ntmax, (int)sizeof(int));
    e_rsle();
    s_rsle(&io___19);
    do_lio(&c__4, &c__1, (char *)&input1_1.dt, (int)sizeof(double));
    e_rsle();

    s_rsle(&io___20);
    do_lio(&c__4, &c__1, (char *)&ludat1_1.parj[40], (int)sizeof(double));
    e_rsle();
    s_rsle(&io___21);
    do_lio(&c__4, &c__1, (char *)&ludat1_1.parj[41], (int)sizeof(double));
    e_rsle();


    s_rsle(&io___22);
    do_lio(&c__3, &c__1, (char *)&ipop, (int)sizeof(int));
    e_rsle();
    if (ipop == 1) {
	hparnt_1.ihpr2[10] = 3;
    }

    s_rsle(&io___24);
    do_lio(&c__4, &c__1, (char *)&ludat1_1.parj[4], (int)sizeof(double));
    e_rsle();

    s_rsle(&io___25);
    do_lio(&c__3, &c__1, (char *)&hparnt_1.ihpr2[5], (int)sizeof(int));
    e_rsle();

    s_rsle(&io___26);
    do_lio(&c__3, &c__1, (char *)&hparnt_1.ihpr2[3], (int)sizeof(int));
    e_rsle();

    s_rsle(&io___27);
    do_lio(&c__4, &c__1, (char *)&hparnt_1.hipr1[13], (int)sizeof(double));
    e_rsle();

    s_rsle(&io___28);
    do_lio(&c__4, &c__1, (char *)&hparnt_1.hipr1[7], (int)sizeof(double));
    e_rsle();

    s_rsle(&io___29);
    do_lio(&c__5, &c__1, (char *)&para2_1.xmu, (int)sizeof(float));
    e_rsle();
    s_rsle(&io___30);
    do_lio(&c__3, &c__1, (char *)&anim_1.izpc, (int)sizeof(int));
    e_rsle();
    s_rsle(&io___31);
    do_lio(&c__5, &c__1, (char *)&para2_1.alpha, (int)sizeof(float));
    e_rsle();

    s_rsle(&io___32);
    do_lio(&c__5, &c__1, (char *)&coal_1.dpcoal, (int)sizeof(float));
    e_rsle();
    s_rsle(&io___33);
    do_lio(&c__5, &c__1, (char *)&coal_1.drcoal, (int)sizeof(float));
    e_rsle();

    s_rsle(&io___34);
    do_lio(&c__3, &c__1, (char *)&ihjsed, (int)sizeof(int));
    e_rsle();

    s_rsle(&io___36);
    do_lio(&c__3, &c__1, (char *)&rndf77_1.nseed, (int)sizeof(int));
    e_rsle();
    s_rsle(&io___37);
    do_lio(&c__3, &c__1, (char *)&rndm3_1.iseedp, (int)sizeof(int));
    e_rsle();
    s_rsle(&io___38);
    do_lio(&c__3, &c__1, (char *)&resdcy_1.iksdcy, (int)sizeof(int));
    e_rsle();
    s_rsle(&io___39);
    do_lio(&c__3, &c__1, (char *)&phidcy_1.iphidcy, (int)sizeof(int));
    e_rsle();

    s_rsle(&io___40);
    do_lio(&c__3, &c__1, (char *)&para7_1.ioscar, (int)sizeof(int));
    e_rsle();

    para8_1.idpert = 0;
    para8_1.npertd = 1;
    para8_1.idxsec = 1;

    cl__1.cerr = 0;
    cl__1.cunit = 24;
    cl__1.csta = 0;
    f_clos(&cl__1);

    if (anim_1.isoft == 1) {
	s_copy(oscar2_1.amptvn, "1.21 (Default)", (int)25, (int)14);
    } else if (anim_1.isoft == 4) {
	s_copy(oscar2_1.amptvn, "2.21 (StringMelting)", (int)25, (int)
		20);
    } else {
	s_copy(oscar2_1.amptvn, "Test-Only", (int)25, (int)9);
    }
    s_wsfe(&io___41);
    do_fio(&c__1, oscar2_1.amptvn, (int)25);
    e_wsfe();
    s_wsfe(&io___42);
    do_fio(&c__1, oscar2_1.amptvn, (int)25);
    e_wsfe();

    if (ihjsed == 11) {
	s_wsle(&io___43);
	do_lio(&c__9, &c__1, "# Read in NSEED in HIJING at run time (e.g. 20"
		"030819):", (int)54);
	e_wsle();
    }
    s_rsle(&io___44);
    do_lio(&c__3, &c__1, (char *)&nseedr, (int)sizeof(int));
    e_rsle();
    if (ihjsed == 11) {
	rndf77_1.nseed = nseedr;
    }

    if (rndf77_1.nseed % 2 == 0) {
	++rndf77_1.nseed;
    }
    if (ihjsed == 11) {
	s_wsle(&io___46);
	do_lio(&c__9, &c__1, "#   read in: ", (int)13);
	do_lio(&c__3, &c__1, (char *)&rndf77_1.nseed, (int)sizeof(int))
		;
	e_wsle();
	s_wsle(&io___47);
	do_lio(&c__9, &c__1, "# Read in NSEED in HIJING at run time:", (
		int)38);
	do_lio(&c__3, &c__1, (char *)&rndf77_1.nseed, (int)sizeof(int))
		;
	e_wsle();
    }
    cl__1.cerr = 0;
    cl__1.cunit = 12;
    cl__1.csta = 0;
    f_clos(&cl__1);

    rand(&rndf77_1.nseed);


    hparnt_1.ihpr2[9] = 1;

    arprnt_1.arpar1[0] = .7f;



    smearz_1.smearp = 0.;
    iamax = max(oscar1_1.iap,oscar1_1.iat);
    d__1 = (float) iamax;
    smearz_1.smearh = pow_dd(&d__1, &c_b92) * 1.2 / ((float) snn_1.efrm /
	     2 / .938);
    anim_1.nevent = nevnt;


    o__1.oerr = 0;
    o__1.ounit = 16;
    o__1.ofnmlen = 12;
    o__1.ofnm = "ana/ampt.dat";
    o__1.orl = 0;
    o__1.osta = "UNKNOWN";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    o__1.oerr = 0;
    o__1.ounit = 14;
    o__1.ofnmlen = 11;
    o__1.ofnm = "ana/zpc.dat";
    o__1.orl = 0;
    o__1.osta = "UNKNOWN";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);



    hijset_(&snn_1.efrm, oscar2_1.frame, proj, targ, &oscar1_1.iap, &
	    oscar1_1.izp, &oscar1_1.iat, &oscar1_1.izt, (int)8, (int)8, 
	    (int)8);
    artset_();
    inizpc_();

    i__1 = nevnt;
    for (j = 1; j <= i__1; ++j) {
	arevt_1.iaevt = j;
	i__2 = run_1.num;
	for (k = 1; k <= i__2; ++k) {
	    arevt_1.iarun = k;
	    if (arevt_1.iaevt == nevnt && arevt_1.iarun == run_1.num) {
		arout_1.iout = 1;
	    }
	    s_wsle(&io___51);
	    do_lio(&c__9, &c__1, " EVENT ", (int)7);
	    do_lio(&c__3, &c__1, (char *)&j, (int)sizeof(int));
	    do_lio(&c__9, &c__1, ", RUN ", (int)6);
	    do_lio(&c__3, &c__1, (char *)&k, (int)sizeof(int));
	    e_wsle();
	    imiss = 0;
L100:
	    hijing_(oscar2_1.frame, &bmin, &bmax, (int)8);
	    arprnt_1.iaint2[0] = hmain1_1.natt;

	    getnp_();

	    if (hparnt_1.ihpr2[19] == 0) {
		goto L2000;
	    }



	    if (hmain1_1.natt == 0) {
		++imiss;
		if (imiss <= 20) {
		    s_wsle(&io___53);
		    do_lio(&c__9, &c__1, "repeated event: natt=0,j,imiss=", (
			    int)31);
		    do_lio(&c__3, &c__1, (char *)&j, (int)sizeof(int));
		    do_lio(&c__3, &c__1, (char *)&imiss, (int)sizeof(
			    int));
		    e_wsle();
		    goto L100;
		} else {
		    s_wsle(&io___54);
		    do_lio(&c__9, &c__1, "missed event: natt=0,j=", (int)
			    23);
		    do_lio(&c__3, &c__1, (char *)&j, (int)sizeof(int));
		    e_wsle();
		    goto L2000;
		}
	    }

	    arini_();
	    arini2_(&k);

	}

	artan1_();
	hjana3_();
	artmn_();
	hjana4_();
	artan2_();
L2000:
	;
    }

    artout_(&nevnt);

    s_stop("", (int)0);
    return 0;
}

 int ampt_ () { MAIN__ (); return 0; }
