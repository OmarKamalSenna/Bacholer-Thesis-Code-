#include "f2c.h"

/* Common Block Declarations */

struct {
    integer natt;
    real eatt;
    integer jatt, nt, np, n0, n01, n10, n11;
} hmain1_;

#define hmain1_1 hmain1_

struct {
    real hipr1[100];
    integer ihpr2[50];
    real hint1[100];
    integer ihnt2[50];
} hparnt_;

#define hparnt_1 hparnt_

struct {
    integer mstu[200];
    real paru[200];
    integer mstj[200];
    real parj[200];
} ludat1_;

#define ludat1_1 ludat1_

struct {
    real arpar1[100];
    integer iapar2[50];
    real arint1[100];
    integer iaint2[50];
} arprnt_;

#define arprnt_1 arprnt_

struct {
    integer iout;
} arout_;

#define arout_1 arout_

struct {
    integer iaevt, iarun;
} arevt_;

#define arevt_1 arevt_

struct {
    doublereal smearp, smearh;
} smearz_;

#define smearz_1 smearz_

struct {
    integer nseed;
} rndf77_;

#define rndf77_1 rndf77_

struct {
    integer nevent, isoft, isflag, izpc;
} anim_;

#define anim_1 anim_

struct {
    doublereal dpcoal, drcoal, ecritl;
} coal_;

#define coal_1 coal_

struct {
    real efrm;
    integer npart1, npart2;
} snn_;

#define snn_1 snn_

struct {
    doublereal xmp, xmu, alpha, rscut2, cutof2;
} para2_;

#define para2_1 para2_

struct {
    integer ioscar;
} para7_;

#define para7_1 para7_

struct {
    integer idpert, npertd, idxsec;
} para8_;

#define para8_1 para8_

struct {
    integer iseedp;
} rndm3_;

#define rndm3_1 rndm3_

struct {
    integer num;
} run_;

#define run_1 run_

struct {
    integer masspr, massta, iseed, iavoid;
    real dt;
} input1_;

#define input1_1 input1_

struct {
    integer ilab, manyb, ntmax, icoll, insys, ipot, mode, imomen, nfreq, 
	    icflow, icrho, icou, kpoten, kmul;
} input2_;

#define input2_1 input2_

struct {
    integer iap, izp, iat, izt;
} oscar1_;

#define oscar1_1 oscar1_

struct {
    char frame[8], amptvn[25];
} oscar2_;

#define oscar2_1 oscar2_

struct {
    integer nsav, iksdcy;
} resdcy_;

#define resdcy_1 resdcy_

struct {
    integer iphidcy;
} phidcy_;

#define phidcy_1 phidcy_

/* Table of constant values */

static integer c__4 = 4;
static integer c__1 = 1;
static integer c__3 = 3;
static integer c__5 = 5;
static integer c__9 = 9;
static doublereal c_b92 = .3333;

/* .....driver program for A Multi-Phase Transport model */
/* Main program */ int MAIN__(void)
{
    /* Format strings */
    static char fmt_111[] = "(a8)";
    static char fmt_50[] = "(\002 \002/11x,\002#############################"
	    "#####################\002/1x,10x,\002#      AMPT (A Multi-Phase "
	    "Transport) model      #\002/1x,10x,\002#               Version"
	    " \002,a20,\002     #\002/1x,10x,\002#                10/01/2008 "
	    "                     #\002/1x,10x,\002##########################"
	    "########################\002/1x,10x,\002 \002)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer f_open(olist *), s_rsle(cilist *), do_lio(integer *, integer *, 
	    char *, ftnlen), e_rsle(void), s_rsfe(cilist *), do_fio(integer *,
	     char *, ftnlen), e_rsfe(void), f_clos(cllist *);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfe(cilist *), e_wsfe(void), s_wsle(cilist *), e_wsle(void);
    double pow_dd(doublereal *, doublereal *);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer j, k;
    static real bmin, bmax;
    static char targ[8], proj[8];
    static integer ipop, iamax;
    extern /* Subroutine */ int arini_(void), srand_(integer *), getnp_(void),
	     artmn_(void);
    static integer imiss, nevnt;
    extern /* Subroutine */ int hjana3_(void), hjana4_(void), arini2_(integer 
	    *), artan1_(void), artan2_(void);
    static integer ihjsed;
    extern /* Subroutine */ int hijing_(char *, real *, real *, ftnlen);
    static integer nseedr;
    extern /* Subroutine */ int hijset_(real *, char *, char *, char *, 
	    integer *, integer *, integer *, integer *, ftnlen, ftnlen, 
	    ftnlen), inizpc_(void), artset_(void), artout_(integer *);

    /* Fortran I/O blocks */
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



/*     parton coalescence radii in case of string melting: */
/*     initialization value for parton cascade: */
/*     initialization value for hadron cascade: */
/* **************** */
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
    do_lio(&c__4, &c__1, (char *)&snn_1.efrm, (ftnlen)sizeof(real));
    e_rsle();
/*     format-read characters (for ALPHA compilers): */
    s_rsfe(&io___2);
    do_fio(&c__1, oscar2_1.frame, (ftnlen)8);
    e_rsfe();
    s_rsfe(&io___3);
    do_fio(&c__1, proj, (ftnlen)8);
    e_rsfe();
    s_rsfe(&io___5);
    do_fio(&c__1, targ, (ftnlen)8);
    e_rsfe();
    s_rsle(&io___7);
    do_lio(&c__3, &c__1, (char *)&oscar1_1.iap, (ftnlen)sizeof(integer));
    e_rsle();
    s_rsle(&io___8);
    do_lio(&c__3, &c__1, (char *)&oscar1_1.izp, (ftnlen)sizeof(integer));
    e_rsle();
    s_rsle(&io___9);
    do_lio(&c__3, &c__1, (char *)&oscar1_1.iat, (ftnlen)sizeof(integer));
    e_rsle();
    s_rsle(&io___10);
    do_lio(&c__3, &c__1, (char *)&oscar1_1.izt, (ftnlen)sizeof(integer));
    e_rsle();
    s_rsle(&io___11);
    do_lio(&c__3, &c__1, (char *)&nevnt, (ftnlen)sizeof(integer));
    e_rsle();
    s_rsle(&io___13);
    do_lio(&c__4, &c__1, (char *)&bmin, (ftnlen)sizeof(real));
    e_rsle();
    s_rsle(&io___15);
    do_lio(&c__4, &c__1, (char *)&bmax, (ftnlen)sizeof(real));
    e_rsle();
/*     flag to select default AMPT or string melting: */
    s_rsle(&io___17);
    do_lio(&c__3, &c__1, (char *)&anim_1.isoft, (ftnlen)sizeof(integer));
    e_rsle();
/*     read initialization value for hadron cascade: */
    s_rsle(&io___18);
    do_lio(&c__3, &c__1, (char *)&input2_1.ntmax, (ftnlen)sizeof(integer));
    e_rsle();
    s_rsle(&io___19);
    do_lio(&c__4, &c__1, (char *)&input1_1.dt, (ftnlen)sizeof(real));
    e_rsle();
/*     parj(41) and (42) are a and b parameters in Lund string fragmentation: */
    s_rsle(&io___20);
    do_lio(&c__4, &c__1, (char *)&ludat1_1.parj[40], (ftnlen)sizeof(real));
    e_rsle();
    s_rsle(&io___21);
    do_lio(&c__4, &c__1, (char *)&ludat1_1.parj[41], (ftnlen)sizeof(real));
    e_rsle();
/*     IHPR2(11)=3 (or 2) allows the popcorn mechanism in PYTHIA and */
/*     increase the net-baryon stopping in rapidity (value HIJING is 1): */
    s_rsle(&io___22);
    do_lio(&c__3, &c__1, (char *)&ipop, (ftnlen)sizeof(integer));
    e_rsle();
    if (ipop == 1) {
	hparnt_1.ihpr2[10] = 3;
    }
/*     PARJ(5) controls the fraction of BMBbar vs BBbar in popcorn: */
    s_rsle(&io___24);
    do_lio(&c__4, &c__1, (char *)&ludat1_1.parj[4], (ftnlen)sizeof(real));
    e_rsle();
/*     shadowing flag in HIJING: */
    s_rsle(&io___25);
    do_lio(&c__3, &c__1, (char *)&hparnt_1.ihpr2[5], (ftnlen)sizeof(integer));
    e_rsle();
/*     quenching flag in HIJING: */
    s_rsle(&io___26);
    do_lio(&c__3, &c__1, (char *)&hparnt_1.ihpr2[3], (ftnlen)sizeof(integer));
    e_rsle();
/*     quenching rate when quenching flag is on (=1.0 GeV/fm): */
    s_rsle(&io___27);
    do_lio(&c__4, &c__1, (char *)&hparnt_1.hipr1[13], (ftnlen)sizeof(real));
    e_rsle();
/*     Minimum pt of hard or semihard scatterings in HIJING: D=2.0 GeV. */
    s_rsle(&io___28);
    do_lio(&c__4, &c__1, (char *)&hparnt_1.hipr1[7], (ftnlen)sizeof(real));
    e_rsle();
/*     read initialization value for parton cascade: */
    s_rsle(&io___29);
    do_lio(&c__5, &c__1, (char *)&para2_1.xmu, (ftnlen)sizeof(doublereal));
    e_rsle();
    s_rsle(&io___30);
    do_lio(&c__3, &c__1, (char *)&anim_1.izpc, (ftnlen)sizeof(integer));
    e_rsle();
    s_rsle(&io___31);
    do_lio(&c__5, &c__1, (char *)&para2_1.alpha, (ftnlen)sizeof(doublereal));
    e_rsle();
/*     quark coalescence radii in momentum and space for string melting: */
    s_rsle(&io___32);
    do_lio(&c__5, &c__1, (char *)&coal_1.dpcoal, (ftnlen)sizeof(doublereal));
    e_rsle();
    s_rsle(&io___33);
    do_lio(&c__5, &c__1, (char *)&coal_1.drcoal, (ftnlen)sizeof(doublereal));
    e_rsle();
/*     flag: read in HIJING random # seed at runtime(1) or from input.ampt(D=0): */
    s_rsle(&io___34);
    do_lio(&c__3, &c__1, (char *)&ihjsed, (ftnlen)sizeof(integer));
    e_rsle();
/*     2 seeds for random number generators in HIJING/hadron cascade and ZPC: */
    s_rsle(&io___36);
    do_lio(&c__3, &c__1, (char *)&rndf77_1.nseed, (ftnlen)sizeof(integer));
    e_rsle();
    s_rsle(&io___37);
    do_lio(&c__3, &c__1, (char *)&rndm3_1.iseedp, (ftnlen)sizeof(integer));
    e_rsle();
    s_rsle(&io___38);
    do_lio(&c__3, &c__1, (char *)&resdcy_1.iksdcy, (ftnlen)sizeof(integer));
    e_rsle();
    s_rsle(&io___39);
    do_lio(&c__3, &c__1, (char *)&phidcy_1.iphidcy, (ftnlen)sizeof(integer));
    e_rsle();
/*     flag for OSCAR output for final partons and hadrons: */
    s_rsle(&io___40);
    do_lio(&c__3, &c__1, (char *)&para7_1.ioscar, (ftnlen)sizeof(integer));
    e_rsle();
/* lin-5/2008 disabled flags for perturbative treatment of deuterons: */
    para8_1.idpert = 0;
    para8_1.npertd = 1;
    para8_1.idxsec = 1;

    cl__1.cerr = 0;
    cl__1.cunit = 24;
    cl__1.csta = 0;
    f_clos(&cl__1);

    if (anim_1.isoft == 1) {
	s_copy(oscar2_1.amptvn, "1.21 (Default)", (ftnlen)25, (ftnlen)14);
    } else if (anim_1.isoft == 4) {
	s_copy(oscar2_1.amptvn, "2.21 (StringMelting)", (ftnlen)25, (ftnlen)
		20);
    } else {
	s_copy(oscar2_1.amptvn, "Test-Only", (ftnlen)25, (ftnlen)9);
    }
    s_wsfe(&io___41);
    do_fio(&c__1, oscar2_1.amptvn, (ftnlen)25);
    e_wsfe();
    s_wsfe(&io___42);
    do_fio(&c__1, oscar2_1.amptvn, (ftnlen)25);
    e_wsfe();
/*     when ihjsed=11: use environment variable at run time for HIJING nseed: */
    if (ihjsed == 11) {
	s_wsle(&io___43);
	do_lio(&c__9, &c__1, "# Read in NSEED in HIJING at run time (e.g. 20"
		"030819):", (ftnlen)54);
	e_wsle();
    }
    s_rsle(&io___44);
    do_lio(&c__3, &c__1, (char *)&nseedr, (ftnlen)sizeof(integer));
    e_rsle();
    if (ihjsed == 11) {
	rndf77_1.nseed = nseedr;
    }
/*     an odd number is needed for the random number generator: */
    if (rndf77_1.nseed % 2 == 0) {
	++rndf77_1.nseed;
    }
    if (ihjsed == 11) {
	s_wsle(&io___46);
	do_lio(&c__9, &c__1, "#   read in: ", (ftnlen)13);
	do_lio(&c__3, &c__1, (char *)&rndf77_1.nseed, (ftnlen)sizeof(integer))
		;
	e_wsle();
	s_wsle(&io___47);
	do_lio(&c__9, &c__1, "# Read in NSEED in HIJING at run time:", (
		ftnlen)38);
	do_lio(&c__3, &c__1, (char *)&rndf77_1.nseed, (ftnlen)sizeof(integer))
		;
	e_wsle();
    }
    cl__1.cerr = 0;
    cl__1.cunit = 12;
    cl__1.csta = 0;
    f_clos(&cl__1);
/*     9/26/03 random number generator for f77 compiler: */
    srand_(&rndf77_1.nseed);

/* .....turn on warning messages */
    hparnt_1.ihpr2[9] = 1;
/*     string formation time: */
    arprnt_1.arpar1[0] = .7f;
/*     smearp is the smearing halfwidth on parton z0, */
/*     set to 0 for now to avoid overflow in eta. */
/*     smearh is the smearing halfwidth on string production point z0. */
    smearz_1.smearp = 0.;
    iamax = max(oscar1_1.iap,oscar1_1.iat);
    d__1 = (doublereal) iamax;
    smearz_1.smearh = pow_dd(&d__1, &c_b92) * 1.2 / ((doublereal) snn_1.efrm /
	     2 / .938);
    anim_1.nevent = nevnt;

/*     AMPT momentum and space info at freezeout: */
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
/* test off for resonance (phi, K*) studies: */
/*      OPEN (17, FILE = 'ana/res-gain.dat', STATUS = 'UNKNOWN') */
/*      OPEN (18, FILE = 'ana/res-loss.dat', STATUS = 'UNKNOWN') */
    hijset_(&snn_1.efrm, oscar2_1.frame, proj, targ, &oscar1_1.iap, &
	    oscar1_1.izp, &oscar1_1.iat, &oscar1_1.izt, (ftnlen)8, (ftnlen)8, 
	    (ftnlen)8);
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
	    do_lio(&c__9, &c__1, " EVENT ", (ftnlen)7);
	    do_lio(&c__3, &c__1, (char *)&j, (ftnlen)sizeof(integer));
	    do_lio(&c__9, &c__1, ", RUN ", (ftnlen)6);
	    do_lio(&c__3, &c__1, (char *)&k, (ftnlen)sizeof(integer));
	    e_wsle();
	    imiss = 0;
L100:
	    hijing_(oscar2_1.frame, &bmin, &bmax, (ftnlen)8);
	    arprnt_1.iaint2[0] = hmain1_1.natt;
/*     evaluate Npart (from primary NN collisions) for both proj and targ: */
	    getnp_();
/*     switch for final parton fragmentation: */
	    if (hparnt_1.ihpr2[19] == 0) {
		goto L2000;
	    }
/*     In the unlikely case of no interaction (even after loop of 20 in HIJING), */
/*     still repeat the event to get an interaction */
/*     (this may have an additional "trigger" effect): */
	    if (hmain1_1.natt == 0) {
		++imiss;
		if (imiss <= 20) {
		    s_wsle(&io___53);
		    do_lio(&c__9, &c__1, "repeated event: natt=0,j,imiss=", (
			    ftnlen)31);
		    do_lio(&c__3, &c__1, (char *)&j, (ftnlen)sizeof(integer));
		    do_lio(&c__3, &c__1, (char *)&imiss, (ftnlen)sizeof(
			    integer));
		    e_wsle();
		    goto L100;
		} else {
		    s_wsle(&io___54);
		    do_lio(&c__9, &c__1, "missed event: natt=0,j=", (ftnlen)
			    23);
		    do_lio(&c__3, &c__1, (char *)&j, (ftnlen)sizeof(integer));
		    e_wsle();
		    goto L2000;
		}
	    }
/* .....ART initialization and run */
	    arini_();
	    arini2_(&k);
/* L1000: */
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

    s_stop("", (ftnlen)0);
    return 0;
} /* MAIN__ */

/* Main program alias */ int ampt_ () { MAIN__ (); return 0; }
