#include "fdmdpd2.h"

/* $Id: surfpot.c 365 2013-07-24 15:35:45Z fuhrmans $ */

#define EXTEND


struct surfpot
{
        double Vpos[TYPE_MAX];
        double rhoC[TYPE_MAX];
        double sigma[TYPE_MAX];
        double C93[TYPE_MAX];
        double A3;
        double B3;
        double C3;
        double A9;
        double B9;
        double C9;
        double F_R0[TYPE_MAX];
        double V_R0[TYPE_MAX];
};

struct surfbump
{
        double Vpos[TYPE_MAX];
        double rhoC[TYPE_MAX];
        double sigma[TYPE_MAX];
        int    reponly[TYPE_MAX];
        double Pos[2];
        double Rbump;
        double H;
        double C93[TYPE_MAX];
        double A3;
        double B3;
        double C3;
        double A9;
        double B9;
        double C9;
        double F_R0[TYPE_MAX];
        double V_R0[TYPE_MAX];
};

struct surfcuff
{
        double rhoC[TYPE_MAX];
        double sigma;
        int    reponly[TYPE_MAX];
        double center[3];
        int    dir;
        double rcuff;
        double halflength;
        //double FtL;                       //Force times effective length - read in from file
        double F;                         //used for computation, calculated from FtL, halflength and sigma - not anymore!!!
        double m;                         //mass for the changing cylinder radius
        //double eff_len;                   //norm for better comparison of cuffs of different length
        //double Ctot[TYPE_MAX];
        //double Aatt;
        //double Batt;
        //double Catt;
        //double Arep;
        //double Brep;
        //double Crep;
        double C93[TYPE_MAX];
        double A3;
        double B3;
        double C3;
        double A9;
        double B9;
        double C9;
        double F_R0[TYPE_MAX];
        double V_R0[TYPE_MAX];
#ifdef EXTEND
	double halflength_times_r2;
#endif
        FILE*  fh;
	int    block;
};

        static const double Rc = 3.; //Cutoff set to 3xSigma
        static const double R1 = 2.; //R1 set to 2xSigma
        static const double R0 = 1.11; //R0 set to 1.11xSigma (truncates too high forces...)

static void print_surfforce(struct surfpot *surfpot);
static void print_surfpot(struct surfpot *surfpot);
static double force(const double d, const int t, const struct surfpot *surfpot) __attribute__((const));
static double pot(const double d, const int t, const struct surfpot *surfpot) __attribute__((const));
static double force_b(const double d, const int t, const struct surfbump *surfbump) __attribute__((const));
static double force_c(const double d, const int t, const struct surfcuff *surfcuff) __attribute__((const));

/* MF THESE FUNCTIONS ARE NOT NEEDED SINCE THE LINE IS ALWAYS PARALLEL TO EITHER X, Y OR Z *****************************
void vsub(VEC a, VEC b, VEC c)
{
        int i;
        for (i=0 ; i<3 ; i++)
                c[i]=a[i]-b[i];
}

double vlen(VEC a)
{
        double l = 0.;
        int i;
        for (i=0 ; i<3 ; i++)
                l+=a[i]*a[i];
        return sqrt(l);
}

void vprod(VEC a, VEC b, VEC c)
{
        c[0]=a[1]*b[2]-a[2]*b[1];
        c[1]=a[2]*b[0]-a[0]*b[2];
        c[2]=a[0]*b[1]-a[1]*b[0];
}

double PLdist(VEC x1, VEC x2, VEC p) //calculates distance between point p and the line defined by x1 and x2
{
        VEC p_x0, p_x1, x1_x0;
        vsub(p,x0,p_x0);
        vsub(p,x1,p_x1);
        vsub(x1,x0,x1_x0);
        vprod(p_x0,p_x1,xprod);
        return vlen(xprod)/vlen(x1_x0);
}
***********************************************************************************************************************/

void surfpot_print_header(FILE *FH, struct surfpot *surfpot)
{
        if (surfpot == NULL) return;

        int i;

        fprintf(FH, "# surfpot ");
        fprintf(FH, "Vpos= ");
        for (i=0 ; i<TYPE_MAX ; i++) {
                fprintf(FH, "%lg ",surfpot->Vpos[i]);
        }
        fprintf(FH, "VrhoC= ");
        for (i=0 ; i<TYPE_MAX ; i++) {
                fprintf(FH, "%lg ",surfpot->rhoC[i]);
        }
        fprintf(FH, "Vsigma= ");
        for (i=0 ; i<TYPE_MAX ; i++) {
                fprintf(FH, "%lg ",surfpot->sigma[i]);
        }
        fprintf(FH, "\n");

}

void surfbump_print_header(FILE *FH, struct surfbump *surfbump)
{
        if (surfbump == NULL) return;

        int i;

        fprintf(FH, "# surfbump ");
        fprintf(FH, "Vpos= ");
        for (i=0 ; i<TYPE_MAX ; i++) {
                fprintf(FH, "%lg ",surfbump->Vpos[i]);
        }
        fprintf(FH, "VrhoC= ");
        for (i=0 ; i<TYPE_MAX ; i++) {
                fprintf(FH, "%lg ",surfbump->rhoC[i]);
        }
        fprintf(FH, "Vsigma= ");
        for (i=0 ; i<TYPE_MAX ; i++) {
                fprintf(FH, "%lg ",surfbump->sigma[i]);
        }
        fprintf(FH, "reponly= ");
        for (i=0 ; i<TYPE_MAX ; i++) {
                fprintf(FH, "%d ",surfbump->reponly[i]);
        }
        fprintf(FH, "BPos= ");
        for (i=0 ; i<2 ; i++) {
                fprintf(FH, "%lg ",surfbump->Pos[i]);
        }
        fprintf(FH, "BRad= ");
        fprintf(FH, "%lg ",surfbump->Rbump);
        fprintf(FH, "BH= ");
        fprintf(FH, "%lg ",surfbump->H);
        fprintf(FH, "\n");
}

void surfcuff_print_header(FILE *FH, struct surfcuff *surfcuff)
{
        if (surfcuff == NULL) return;

        int i;

        fprintf(FH, "# surfcuff ");
        fprintf(FH, "org= ");
        for (i=0 ; i<3 ; i++) {
                fprintf(FH, "%lg ",surfcuff->center[i]);
        }
        fprintf(FH, "dir= ");
        fprintf(FH, "%d ",surfcuff->dir);
        fprintf(FH, "rad= ");
        fprintf(FH, "%lg ",surfcuff->rcuff);
        fprintf(FH, "halflength= ");
        fprintf(FH, "%lg ",surfcuff->halflength);
        fprintf(FH, "F= ");
        fprintf(FH, "%lg ",surfcuff->F);
        fprintf(FH, "m= ");
        fprintf(FH, "%lg ",surfcuff->m);
        fprintf(FH, "C= ");
        for (i=0 ; i<TYPE_MAX ; i++) {
                //fprintf(FH, "%lg ",surfcuff->Ctot[i]);
                fprintf(FH, "%lg ",surfcuff->rhoC[i]);
        }
        fprintf(FH, "sigma= ");
        fprintf(FH, "%lg ",surfcuff->sigma);
        fprintf(FH, "reponly= ");
        for (i=0 ; i<TYPE_MAX ; i++) {
                fprintf(FH, "%d ",surfcuff->reponly[i]);
        }
        fprintf(FH, "\n");
}

void surfpot_read_header(FILE *FH, struct beads *b)
{
        int i;

        b->surfpot = malloc(sizeof(*b->surfpot));
        if (b->surfpot == NULL) novm("surfpot");

        if (ismaster) {
        	printf("Using surface potential\n");
	}

	fscanf(FH, "# surfpot Vpos= ");
	for (i = 0; i < ARRAY_SIZE(b->surfpot->Vpos); i++) {
		if (fscanf(FH, "%lg ", &b->surfpot->Vpos[i]) != 1)
			fatal(EINVAL, "surfpot Vpos!");
	}
	fscanf(FH, "VrhoC= ");
	for (i = 0; i < ARRAY_SIZE(b->surfpot->rhoC); i++) {
		if (fscanf(FH, "%lg ", &b->surfpot->rhoC[i]) != 1)
			fatal(EINVAL, "surfpot rhoC");
	}
	fscanf(FH, "Vsigma= ");
	for (i = 0; i < ARRAY_SIZE(b->surfpot->sigma); i++) {
		if (fscanf(FH, "%lg ", &b->surfpot->sigma[i]) != 1)
			fatal(EINVAL, "surfpot sigma");
	}


        b->surfpot->A3 = - 3*( (3+4)*Rc - (3+1)*R1 ) / pow(Rc,3+2) / pow(Rc-R1,2);
        b->surfpot->A9 = - 9*( (9+4)*Rc - (9+1)*R1 ) / pow(Rc,9+2) / pow(Rc-R1,2);
        b->surfpot->B3 =   3*( (3+3)*Rc - (3+1)*R1 ) / pow(Rc,3+2) / pow(Rc-R1,3);
        b->surfpot->B9 =   9*( (9+3)*Rc - (9+1)*R1 ) / pow(Rc,9+2) / pow(Rc-R1,3);
        b->surfpot->C9 = 1./CUBE(CUBE(Rc)) - b->surfpot->A9/3.*CUBE(Rc-R1)-b->surfpot->B9/4.*SQR(SQR(Rc-R1));
        b->surfpot->C3 = 1./CUBE(Rc) - b->surfpot->A3/3.*CUBE(Rc-R1)-b->surfpot->B3/4.*SQR(SQR(Rc-R1));
        for (i = 0 ; i < TYPE_MAX ; i++) {
                b->surfpot->C93[i] = b->surfpot->rhoC[i]*M_PI/6.;
                b->surfpot->F_R0[i] = b->surfpot->C93[i]*(9./pow(R0,10)-3./pow(R0,4));
                b->surfpot->V_R0[i] = b->surfpot->C93[i]*(pow(R0,-9) - b->surfpot->C9 - pow(R0,-3) + b->surfpot->C3);
        }


        /* loesch mich */
        //char buf[LINE_MAX];
        //fgets(buf, sizeof(buf), FH);
        //fprintf(stderr, "##%s##\n", buf);

        print_surfforce(b->surfpot);
        print_surfpot(b->surfpot);
}

void surfbump_read_header(FILE *FH, struct beads *b)
{
        int i;

        b->surfbump = malloc(sizeof(*b->surfbump));
        if (b->surfbump == NULL) novm("surfbump");

        if (ismaster) {
        	printf("Using surface potential with bump\n");
	}

	fscanf(FH, "# surfbump Vpos= ");
	for (i = 0; i < ARRAY_SIZE(b->surfbump->Vpos); i++) {
		if (fscanf(FH, "%lg ", &b->surfbump->Vpos[i]) != 1)
			fatal(EINVAL, "surfbump Vpos!");
	}
	fscanf(FH, "VrhoC= ");
	for (i = 0; i < ARRAY_SIZE(b->surfbump->rhoC); i++) {
		if (fscanf(FH, "%lg ", &b->surfbump->rhoC[i]) != 1)
			fatal(EINVAL, "surfbump rhoC");
	}
	fscanf(FH, "Vsigma= ");
	for (i = 0; i < ARRAY_SIZE(b->surfbump->sigma); i++) {
		if (fscanf(FH, "%lg ", &b->surfbump->sigma[i]) != 1)
			fatal(EINVAL, "surfbump sigma");
	}
	fscanf(FH, "reponly= ");
	for (i = 0; i < ARRAY_SIZE(b->surfbump->reponly); i++) {
		if (fscanf(FH, "%d ", &b->surfbump->reponly[i]) != 1)
			fatal(EINVAL, "surfbump reponly");
	}
	fscanf(FH, "BPos= ");
	for (i = 0; i < ARRAY_SIZE(b->surfbump->Pos); i++) {
		if (fscanf(FH, "%lg ", &b->surfbump->Pos[i]) != 1)
			fatal(EINVAL, "surfbump Pos");
	}
	fscanf(FH, "BRad= ");
	if (fscanf(FH, "%lg ", &b->surfbump->Rbump) != 1)
	        fatal(EINVAL, "surfbump BRad");
	fscanf(FH, "BH= ");
	if (fscanf(FH, "%lg ", &b->surfbump->H) != 1)
	        fatal(EINVAL, "surfbump BH");

        for (i = 0 ; i < TYPE_MAX ; i++) {
                if (b->surfbump->H < b->surfbump->sigma[i]*Rc) {
                        fatal(EINVAL, "BH must be >= sigma*Rc for all types");
                }
        }

        b->surfbump->A3 = - 3*( (3+4)*Rc - (3+1)*R1 ) / pow(Rc,3+2) / pow(Rc-R1,2);
        b->surfbump->A9 = - 9*( (9+4)*Rc - (9+1)*R1 ) / pow(Rc,9+2) / pow(Rc-R1,2);
        b->surfbump->B3 =   3*( (3+3)*Rc - (3+1)*R1 ) / pow(Rc,3+2) / pow(Rc-R1,3);
        b->surfbump->B9 =   9*( (9+3)*Rc - (9+1)*R1 ) / pow(Rc,9+2) / pow(Rc-R1,3);
        b->surfbump->C9 = 1./CUBE(CUBE(Rc)) - b->surfbump->A9/3.*CUBE(Rc-R1)-b->surfbump->B9/4.*SQR(SQR(Rc-R1));
        b->surfbump->C3 = 1./CUBE(Rc) - b->surfbump->A3/3.*CUBE(Rc-R1)-b->surfbump->B3/4.*SQR(SQR(Rc-R1));
        for (i = 0 ; i < TYPE_MAX ; i++) {
                b->surfbump->C93[i] = b->surfbump->rhoC[i]*M_PI/6.;
                b->surfbump->F_R0[i] = b->surfbump->C93[i]*(9./pow(R0,10)-3./pow(R0,4));
                b->surfbump->V_R0[i] = b->surfbump->C93[i]*(pow(R0,-9)-b->surfbump->C9-pow(R0,-3)+b->surfbump->C3);
        }
        print_surfbump_forces(b,b->surfbump);
}

void surfcuff_read_header(FILE *FH, struct beads *b)
{
        int i;

        b->surfcuff = malloc(sizeof(*b->surfcuff));
        if (b->surfcuff == NULL) novm("surfcuff");

        if (ismaster) {
        	printf("Using cuff potential\n");
	}

	fscanf(FH, "# surfcuff org= ");
	for (i = 0; i < ARRAY_SIZE(b->surfcuff->center); i++) {
		if (fscanf(FH, "%lg ", &b->surfcuff->center[i]) != 1)
			fatal(EINVAL, "surfcuff org!");
	}
	fscanf(FH, "dir= ");
	if (fscanf(FH, "%d ", &b->surfcuff->dir) != 1)
	        fatal(EINVAL, "surfcuff dir");
	fscanf(FH, "rad= ");
	if (fscanf(FH, "%lg ", &b->surfcuff->rcuff) != 1)
	        fatal(EINVAL, "surfcuff rad");
	fscanf(FH, "halflength= ");
	if (fscanf(FH, "%lg ", &b->surfcuff->halflength) != 1)
	        fatal(EINVAL, "surfcuff length");
	fscanf(FH, "F= ");
	if (fscanf(FH, "%lg ", &b->surfcuff->F) != 1)
	        fatal(EINVAL, "surfcuff F");
	fscanf(FH, "m= ");
	if (fscanf(FH, "%lg ", &b->surfcuff->m) != 1)
	        fatal(EINVAL, "surfcuff m");
	fscanf(FH, "C= ");
	//for (i = 0; i < ARRAY_SIZE(b->surfcuff->Ctot); i++) {
	//	if (fscanf(FH, "%lg ", &b->surfcuff->Ctot[i]) != 1)
	for (i = 0; i < ARRAY_SIZE(b->surfcuff->rhoC); i++) {
		if (fscanf(FH, "%lg ", &b->surfcuff->rhoC[i]) != 1)
			fatal(EINVAL, "surfcuff C");
	}
	fscanf(FH, "sigma= ");
	if (fscanf(FH, "%lg ", &b->surfcuff->sigma) != 1)
		fatal(EINVAL, "surfcuff sigma");
	fscanf(FH, "reponly= ");
	for (i = 0; i < ARRAY_SIZE(b->surfcuff->reponly); i++) {
		if (fscanf(FH, "%d ", &b->surfcuff->reponly[i]) != 1)
			fatal(EINVAL, "surfcuff reponly");
	}

        if (b->surfcuff->rcuff < 0. || b->surfcuff->halflength < 0.)
                fatal(EINVAL, "rad and length must be >= 0");
        if (b->surfcuff->rcuff < b->surfcuff->sigma*Rc)
                fatal(EINVAL, "potentials may not overlap in center of cuff");
        if (b->surfcuff->dir > 2 || b->surfcuff->dir < 0)
                fatal(EINVAL, "dir must be 0, 1 or 2");
        if (b->surfcuff->m <= 0)
                fatal(EINVAL, "mass must be positive");

        //surfcuff->eff_len = 2. * (surfcuff->halflength + .25*M_PI*b->surfcuff->sigma*Rc);
        //surfcuff->F = surfcuff->FtL / (2. * (surfcuff->halflength + .25*M_PI*b->surfcuff->sigma*Rc)); //F is normalized according to the length of the cuff
                                                                                                      //the effective length used for normalization is
                                                                                                      //2. * (surfcuff->halflength + .25*M_PI*b->surfcuff->sigma*Rc)
                                                                                                      //i.e. the integrated component perpendicular to the mantle of
                                                                                                      //the force unit vector
//#define ATT (6)
//#define REP (12)
//
//        b->surfcuff->Aatt = - ATT*( (ATT+4)*Rc - (ATT+1)*R1 ) / pow(Rc,ATT+2) / pow(Rc-R1,2);
//        b->surfcuff->Arep = - REP*( (REP+4)*Rc - (REP+1)*R1 ) / pow(Rc,REP+2) / pow(Rc-R1,2);
//        b->surfcuff->Batt =   ATT*( (ATT+3)*Rc - (ATT+1)*R1 ) / pow(Rc,ATT+2) / pow(Rc-R1,3);
//        b->surfcuff->Brep =   REP*( (REP+3)*Rc - (REP+1)*R1 ) / pow(Rc,REP+2) / pow(Rc-R1,3);
//        b->surfcuff->Crep = 1./pow(Rc,REP) - b->surfcuff->Arep/3.*CUBE(Rc-R1)-b->surfcuff->Brep/4.*SQR(SQR(Rc-R1));
//        b->surfcuff->Catt = 1./pow(Rc,ATT) - b->surfcuff->Aatt/3.*CUBE(Rc-R1)-b->surfcuff->Batt/4.*SQR(SQR(Rc-R1));
//        for (i = 0 ; i < TYPE_MAX ; i++) {
//                b->surfcuff->Ctot[i] = b->surfcuff->C[i]*M_PI/6.;
//                b->surfcuff->F_R0[i] = b->surfcuff->Ctot[i]*(REP/pow(R0,REP+1)-ATT/pow(R0,ATT+1));
//                b->surfcuff->V_R0[i] = b->surfcuff->Ctot[i]*(pow(R0,-REP) - b->surfcuff->Crep - pow(R0,-ATT) + b->surfcuff->Catt);
//        }
        b->surfcuff->A3 = - 3*( (3+4)*Rc - (3+1)*R1 ) / pow(Rc,3+2) / pow(Rc-R1,2);
        b->surfcuff->A9 = - 9*( (9+4)*Rc - (9+1)*R1 ) / pow(Rc,9+2) / pow(Rc-R1,2);
        b->surfcuff->B3 =   3*( (3+3)*Rc - (3+1)*R1 ) / pow(Rc,3+2) / pow(Rc-R1,3);
        b->surfcuff->B9 =   9*( (9+3)*Rc - (9+1)*R1 ) / pow(Rc,9+2) / pow(Rc-R1,3);
        b->surfcuff->C9 = 1./CUBE(CUBE(Rc)) - b->surfcuff->A9/3.*CUBE(Rc-R1)-b->surfcuff->B9/4.*SQR(SQR(Rc-R1));
        b->surfcuff->C3 = 1./CUBE(Rc) - b->surfcuff->A3/3.*CUBE(Rc-R1)-b->surfcuff->B3/4.*SQR(SQR(Rc-R1));
        for (i = 0 ; i < TYPE_MAX ; i++) {
                b->surfcuff->C93[i] = b->surfcuff->rhoC[i]*M_PI/6.;
                b->surfcuff->F_R0[i] = b->surfcuff->C93[i]*(9./pow(R0,10)-3./pow(R0,4));
                b->surfcuff->V_R0[i] = b->surfcuff->C93[i]*(pow(R0,-9) - b->surfcuff->C9 - pow(R0,-3) + b->surfcuff->C3);
        }
        //print_surfcuff_forces(b,b->surfcuff);
        if (ismaster) {
	        char fn[PATH_MAX];
                sprintf(fn, "cuff%09lu.dat", b->step);
                b->surfcuff->fh = fopen(fn,"w");
                if (b->surfcuff->fh == NULL)
                        fatal(EIO, "Couldn't open surfcuff.dat for writing");
                fprintf(b->surfcuff->fh,"# using a surfcuff with F=%lg, m=%lg and halflength=%lg\n",
                                b->surfcuff->F,b->surfcuff->m,b->surfcuff->halflength);
                fprintf(b->surfcuff->fh,"# time, rcuff, perpendicular force on mantle, circumference at zero transition of force\n");
        } 
#ifdef EXTEND
	b->surfcuff->halflength_times_r2 = b->surfcuff->halflength * SQR(b->surfcuff->rcuff);
#endif

	surfcuff_print_forces(b,b->surfcuff);
}

static double force(const double d, const int t, const struct surfpot *surfpot)
{
        if (d >= Rc) {
                return 0.;
                //Vpot also 0
        }
        else if (d > R1) {
                double dmr1_2 = SQR(d-R1);
                double dmr1_3 = CUBE(d-R1);
                //double dmr1_4 = SQR(dmr1_2);
                //double d3  = CUBE(d);
                double d4  = SQR(SQR(d));
                //double d9  = SQR(d4) * d;
                double d10 = SQR(d4) * SQR(d);
                //Vpot = surfpot->C93[t] * ( 1./d9 - surfpot->A9/3.*dmr1_3 - surfpot->B9*dmr1_4 - surfpot->C9
                //                         - 1./d3 + surfpot->A3/3.*dmr1_3 + surfpot->B3*dmr1_4 + surfpot->C3);
                return surfpot->C93[t] * ( 9./d10 + surfpot->A9*dmr1_2 + surfpot->B9*dmr1_3 
                                         - 3./ d4 - surfpot->A3*dmr1_2 - surfpot->B3*dmr1_3 );
        }
        else if (d > R0) {
                //double d3  = CUBE(d);
                double d4  = SQR(SQR(d));
                //double d9  = SQR(d4) * d;
                double d10 = SQR(d4) * SQR(d);
                //Vpot = surfpot->C93[t] * ( 1./d9 - surfpot->C9 - 1./d3 + surfpot->C3 );
                return surfpot->C93[t]*(9./d10-3./d4);
        }
        else { //d is below r0xSigma now
                //F=const=F(r0)
                //Vpot=V_R0[t];
                return surfpot->F_R0[t];
        }
}

//PROBLEM: IF I USE VARIABLE EXPONENTS, I WILL HAVE TO USE POW TO CALCULATE THE POWERS OF THE DISTANCES, WHICH MAY BE SLOW...
//SO - DECIDE ON FIXED EXPONENTS AFTER ALL???
static double force_c(const double d, const int t, const struct surfcuff *surfcuff) //USES VARIABLE EXPONENTS - MF MAYBE UPDATE ALL FORCE CALLS TO THIS FUNCTIONS...
{
        if (d >= Rc) {
                return 0.;
                //Vpot also 0
        }
        else if (d > R1) {
                double dmr1_2 = SQR(d-R1);
                double dmr1_3 = CUBE(d-R1);
                //double dmr1_4 = SQR(dmr1_2);
                //double d3  = CUBE(d);
                double d4  = SQR(SQR(d));
                //double d9  = SQR(d4) * d;
                double d10 = SQR(d4) * SQR(d);
                //Vpot = surfcuff->C93[t] * ( 1./d9 - surfcuff->A9/3.*dmr1_3 - surfcuff->B9*dmr1_4 - surfcuff->C9
                //                         - 1./d3 + surfcuff->A3/3.*dmr1_3 + surfcuff->B3*dmr1_4 + surfcuff->C3);
                //return surfcuff->Ctot[t] * ( REP / d10 + surfcuff->Arep*dmr1_2 + surfcuff->Brep*dmr1_3 
                //                           - ATT /  d4 - surfcuff->Aatt*dmr1_2 - surfcuff->Batt*dmr1_3 );
                return surfcuff->C93[t] * ( 9./d10 + surfcuff->A9*dmr1_2 + surfcuff->B9*dmr1_3 
                                         - 3./ d4 - surfcuff->A3*dmr1_2 - surfcuff->B3*dmr1_3 );
        }
        else if (d > R0) {
                //double d3  = CUBE(d);
                double d4  = SQR(SQR(d));
                //double d9  = SQR(d4) * d;
                double d10 = SQR(d4) * SQR(d);
                //Vpot = surfcuff->C93[t] * ( 1./d9 - surfcuff->C9 - 1./d3 + surfcuff->C3 );
                //return surfcuff->Ctot[t]*(REP/d10-ATT/d4);
                return surfcuff->C93[t]*(9./d10-3./d4);
        }
        else { //d is below r0xSigma now
                //F=const=F(r0)
                //Vpot=V_R0[t];
                return surfcuff->F_R0[t];
        }
}

static double pot(const double d, const int t, const struct surfpot *surfpot)
{
        if (d >= Rc) {
                return 0.;
                //Vpot also 0
        }
        else if (d > R1) {
                double dmr1_2 = SQR(d-R1);
                double dmr1_3 = CUBE(d-R1);
                double dmr1_4 = SQR(dmr1_2);
                double d3  = CUBE(d);
                double d4  = SQR(SQR(d));
                double d9  = SQR(d4) * d;
                //double d10 = SQR(d4) * SQR(d);
                return surfpot->C93[t] * ( 1./d9 - surfpot->A9/3.*dmr1_3 - surfpot->B9/4.*dmr1_4 - surfpot->C9
                                         - 1./d3 + surfpot->A3/3.*dmr1_3 + surfpot->B3/4.*dmr1_4 + surfpot->C3);
                //return surfpot->C93[t] * ( 9./d10 + surfpot->A9*dmr1_2 + surfpot->B9*dmr1_3 
                //                         - 3./ d4 - surfpot->A3*dmr1_2 - surfpot->B3*dmr1_3 );
        }
        else if (d > R0) {
                double d3  = CUBE(d);
                double d4  = SQR(SQR(d));
                double d9  = SQR(d4) * d;
                //double d10 = SQR(d4) * SQR(d);
                return surfpot->C93[t] * ( 1./d9 - surfpot->C9 - 1./d3 + surfpot->C3 );
                //return surfpot->C93[t]*(9./d10-3./d4);
        }
        else { //d is below r0xSigma now
                //F=const=F(r0)
                return surfpot->V_R0[t] + surfpot->F_R0[t] * (R0 - d);
                //return surfpot->F_R0[t];
        }
}

static double force_b(const double d, const int t, const struct surfbump *surfbump)
{
        //ATTENTION: Expressions for Vpot are missing C9 and C3 and have to be corrected, if uncommented
        if (d >= Rc) {
                return 0.;
                //Vpot also 0
        }
        else if (d > R1) {
                double dmr1_2 = SQR(d-R1);
                double dmr1_3 = CUBE(d-R1);
                //double dmr1_4 = SQR(dmr1_2);
                //double d3  = CUBE(d);
                double d4  = SQR(SQR(d));
                //double d9  = SQR(d4) * d;
                double d10 = SQR(d4) * SQR(d);
                //Vpot = surfbump->C93[t] * ( 1./d9 - surfbump->A9/3.*dmr1_3 - surfbump->B9*dmr1_4
                //                          - 1./d3 + surfbump->A3/3.*dmr1_3 + surfbump->B3*dmr1_4 );
                return surfbump->C93[t] * ( 9./d10 + surfbump->A9*dmr1_2 + surfbump->B9*dmr1_3 
                                          - 3./ d4 - surfbump->A3*dmr1_2 - surfbump->B3*dmr1_3 );
        }
        else if (d > R0) {
                //double d3  = CUBE(d);
                double d4  = SQR(SQR(d));
                //double d9  = SQR(d4) * d;
                double d10 = SQR(d4) * SQR(d);
                //Vpot = surfbump->C93[t] * ( 1./d9 - 1./d3 );
                return surfbump->C93[t]*(9./d10-3./d4);
        }
        else { //d is below r0xSigma now
                //Vpot=V_R0[t];
                return surfbump->F_R0[t];
        }
}

void surfpot_calc_forces(struct beads *b, struct surfpot *surfpot)
{
        if (surfpot == NULL) return;

        int i;
        VEC *x = b->x_intra;
        VEC *f = b->f_intra;
        PASSPORT *p = b->passport + col_index * b->groupsize;

        for (i = 0; i < b->groupsize; i++) {
                if (!GET_EXISTS(p[i])) continue; // in order to skip ghost particles
                int t = GET_TYPE(p[i]);
                double d = (x[i][NORMAL]-surfpot->Vpos[t])/surfpot->sigma[t];

                f[i][NORMAL] += force(d, t, surfpot);
                //fprintf(stderr,"f+=%lg d/sigma=%lg\n", force(d, t, surfpot), d);
        }
}


void surfbump_calc_forces(struct beads *b, struct surfbump *surfbump)
{
        if (surfbump == NULL) return;


        int i,idim;
	const VEC lh = {.5 * b->l[0], .5 * b->l[1], .5 * b->l[2]};
        const int dim[3] = {TANG1,TANG2,NORMAL};
        VEC XZdir;
        double base_curv_ref[TYPE_MAX][3];
        double top_curv_ref[TYPE_MAX][3];
        for (i = 0 ; i < TYPE_MAX ; i++) {
                base_curv_ref[i][0] = surfbump->Rbump + surfbump->sigma[i]*Rc;
                base_curv_ref[i][1] = 0.;
                base_curv_ref[i][2] =  surfbump->sigma[i]*Rc;
                top_curv_ref[i][0] = surfbump->Rbump;
                top_curv_ref[i][1] = 0.;
                top_curv_ref[i][2] = surfbump->H;
        }
        VEC XYdir;
        double d_xy, d_xy2, vecl, vecl2;

        VEC *x = b->x_intra;
        VEC *f = b->f_intra;
        PASSPORT *p = b->passport + col_index * b->groupsize;

        for (i = 0; i < b->groupsize; i++) {
                if (!GET_EXISTS(p[i])) continue; // in order to skip ghost particles
                int t = GET_TYPE(p[i]);
                if (surfbump->rhoC[t] == 0) continue;
                for (idim = 0, d_xy2 = 0. ; idim < 2 ; idim++) {
                        XYdir[idim] = x[i][dim[idim]] - b->surfbump->Pos[idim];
                        XYdir[idim] += (XYdir[idim] < -lh[dim[idim]]) ? b->l[dim[idim]]: 0.; // PBC
                        XYdir[idim] -= (XYdir[idim] >  lh[dim[idim]]) ? b->l[dim[idim]]: 0.; // minimum distance MF
                        d_xy2 += SQR(XYdir[idim]);
                }
                d_xy = sqrt(d_xy2);
                //XYdir[2] = x[i][dim[2]] - surfbump->Vpos[t];




                if (d_xy >= surfbump->sigma[t]*Rc + surfbump->Rbump) {
                        // particle is not within horizontal range of the bump
                        double d = (x[i][NORMAL]-surfbump->Vpos[t])/surfbump->sigma[t];
                        double frc = force_b(d, t, surfbump);
                        if (!surfbump->reponly[t] || frc>0) {
                                f[i][NORMAL] += frc;
                        } 
//#define STRANGETOP
#ifndef STRANGETOP                        
                } else if (d_xy <= surfbump->Rbump) {
                        // particle is above flat top of bump
                        double d = (x[i][NORMAL] - surfbump->Vpos[t] - surfbump->H)/surfbump->sigma[t];
                        double frc = force_b(d, t, surfbump);
                        if (!surfbump->reponly[t] || frc>0) {
                                f[i][NORMAL] += frc;
                        } 
#else
                } else if (d_xy <= surfbump->Rbump) {
                        // particle is above flat top of bump
                        // ALTERED VERSION!!!!! STRANGETOP IS HERE!!!
                        double d = (x[i][NORMAL] - surfbump->Vpos[t] - surfbump->H)/surfbump->sigma[t];
                        double frc = force_b(d, t, surfbump);
                        if (/*!surfbump->reponly[t] ||*/ frc>0) {
                                f[i][NORMAL] += frc;
                        }
#endif
                } else {
                        if (x[i][NORMAL] - surfbump->Vpos[t] < surfbump->sigma[t]*Rc) {
                                // particle in curved region at base of bump
                                // we first deal with the particle's projection into the T1,N plane (e.g. XZ)
                                // for this we can use the radial distance we already know and x[N]
                                VEC XZcoord = {d_xy, 0., x[i][NORMAL]-surfbump->Vpos[t]};
                                for (idim = 0, vecl2 = 0. ; idim < 3 ; idim++) {
                                        XZdir[idim] = base_curv_ref[t][idim] - XZcoord[idim];
                                        vecl2 += SQR(XZdir[idim]);
                                }
                                vecl = sqrt(vecl2);
                                for (idim = 0 ; idim < 3 ; idim++) {
                                        XZdir[idim] /= vecl;
                                } // now we have the planar direction (unit vector pointing away from surface) and the distance from the reference
                                //double d = surfbump->sigma[t]*Rc - vecl;
                                double d = Rc - vecl / surfbump->sigma[t];
                                double frc = force_b(d, t, surfbump);
                                if (!surfbump->reponly[t] || frc>0) {
                                        VEC dir;
                                        for (idim = 0 ; idim < 2 ; idim++) { //scale XYdir so that it has a length corresponding to the X component of XZdir
                                                dir[idim] = XYdir[idim] / d_xy * XZdir[0];
                                        }
                                        dir[2] = XZdir[2];
                                        for (idim = 0 ; idim < 3 ; idim++) {
                                                f[i][dim[idim]] += frc * dir[idim];
                                        }
                                }
                        } else if (x[i][NORMAL]-surfbump->Vpos[t] <= b->surfbump->H) {
                                // particle next to vertical side of bump
                                // XYdir has already been calculated while measuring d_xy
                                double d = (d_xy - surfbump->Rbump) / surfbump->sigma[t];
                                double frc = force_b(d, t, surfbump);
                                if (!surfbump->reponly[t] || frc>0) {
                                        for (idim = 0 ; idim < 2 ; idim++) {
                                                XYdir[idim] /= d_xy;
                                                f[i][dim[idim]] += frc * XYdir[idim]; 
                                        }
                                }
                        } else if (x[i][NORMAL]-surfbump->Vpos[t] <  b->surfbump->H + surfbump->sigma[t]*Rc) {
                                // particle next to round tip of bump (but not above the flat top)
                                // we first deal with the particle's projection into the T1,N plane (e.g. XZ)
                                // for this we can use the radial distance we already know and x[N]
                                VEC XZcoord = {d_xy, 0., x[i][NORMAL]-surfbump->Vpos[t]};
                                for (idim = 0, vecl2 = 0. ; idim < 3 ; idim++) {
                                        XZdir[idim] = XZcoord[idim] - top_curv_ref[t][idim];
                                        vecl2 += SQR(XZdir[idim]);
                                }
                                vecl = sqrt(vecl2);
                                for (idim = 0 ; idim < 3 ; idim++) {
                                        XZdir[idim] /= vecl;
                                } // now we have the planar direction (unit vector pointing away from surface) and the distance from the reference
                                //double d = vecl - surfbump->sigma[t]*Rc;
                                double d = vecl / surfbump->sigma[t] ;
                                double frc = force_b(d, t, surfbump);
                                if (!surfbump->reponly[t] || frc>0) {
                                        VEC dir;
                                        for (idim = 0 ; idim < 2 ; idim++) { //scale XYdir so that it has a length corresponding to the X component of XZdir
                                                dir[idim] = XYdir[idim] / d_xy * XZdir[0];
                                        }
                                        dir[2] = XZdir[2];
#ifndef STRANGETOP                        
                                        for (idim = 0 ; idim < 3 ; idim++) {
                                                f[i][dim[idim]] += frc * dir[idim];
                                        }
#else
                                        if (frc>0) {
                                                for (idim = 0 ; idim < 3 ; idim++) {
                                                        f[i][dim[idim]] += frc * dir[idim];
                                                }
                                        // STRANGETOP - attractive forces in horizontal direction are ignored in order to stretch the membrane
                                        // one might consider amplifying repulsive horizontal forces
                                        // attractive forces in normal direction are applied normally (i.e. also if attractive)
                                        } else {
                                                f[i][NORMAL] += frc * dir[NORMAL];
                                        }
#endif
                                }
                        }
                }
                //fprintf(stderr,"f+=%lg d/sigma=%lg\n", force(d, t, surfbump), d);
        }
}

void surfcuff_calc_forces(struct beads *b, struct surfcuff *surfcuff)
{
        if (surfcuff == NULL) return;

        int i,idim,t;
	const VEC lh = {.5 * b->l[0], .5 * b->l[1], .5 * b->l[2]};
        //const int dim[3] = {TANG1,TANG2,NORMAL};

        VEC tolinedir; //direction normal to the line running through the cuff's center
        double d_toline, d_toline2, d_online, frc, tmp, dist, f_on_mantle=0., sign;
        /****** list of distances ***************************************************************************************
         *      d_toline: initially length of tolinedir                                                                 *
         *                later it will become the distance to the cuff's mantle if the cuff was long enough            *
         *      dist: distance to the reference point for the force in this region                                      *
         *      d_online: initially the displacement along the direction of the cuff relative to the cuff's center,     *
         *                later it will be changed to be relative to one of the cuff's end                              *
         ****************************************************************************************************************/

        //int offline_ref_dir = (surfcuff->dir == 0) ? 1 : 0;

        //for (idim = 0 , i = 0 ; idim < 3 ; idim++) {
        //        if (idim != surfcuff->dir) {
        //                offline_dims[i] = idim;
        //                i++;
        //        }
        //}


        //double offline_ref = surfcuff->center[offline_ref_dir] + surfcuff->rcuff;
        //double offline_ref = surfcuff->center[offline_dims[0]] + surfcuff->rcuff; //either offline_dim would work...
        //oops, I guess we don't need this and it is wrong - should be only rcuff!

        //double online_ref = surfcuff->center[surfcuff->dir] + surfcuff->halflength; //probably not needed?

        VEC *x = b->x_intra;
        VEC *f = b->f_intra;
        PASSPORT *p = b->passport + col_index * b->groupsize;

        for (i = 0; i < b->groupsize; i++) {
                if (!GET_EXISTS(p[i])) continue; // in order to skip ghost particles
                t = GET_TYPE(p[i]);
                if (surfcuff->rhoC[t] == 0) continue;
                for (idim = 0, d_toline2=0. ; idim < 3 ; idim++) {
                        tmp = x[i][idim] - surfcuff->center[idim];
                        tmp += (tmp < -lh[idim]) ? b->l[idim]: 0.; // PBC
                        tmp -= (tmp >  lh[idim]) ? b->l[idim]: 0.; // minimum distance MF
                        if (idim == surfcuff->dir) {
                                d_online = tmp;
                                tolinedir[idim] = 0.;
                        } else {
                                tolinedir[idim] = -tmp;
                                d_toline2 += SQR(tmp);
                        }
                }
                d_toline = sqrt(d_toline2);
		for (idim = 0 ; idim < 3 ; idim++) { //d_toline should be direction to mantle, not to center line
			tolinedir[idim] *= (surfcuff->rcuff-d_toline)/d_toline;
		}
		//fprintf(stderr,"d_online:%lg d_toline:%lg tolinedir:%lg %lg %lg\n",d_online,d_toline,tolinedir[0],tolinedir[1],tolinedir[2]);
                sign = 1.; //standard assumption is that the particle is inside or next to the cuff
                if (d_toline > surfcuff->rcuff)
                        sign = -1.; //bead is somehow outside of cuff resulting into a different force direction...

                if (fabs(d_online) >= surfcuff->halflength + Rc*surfcuff->sigma) {
                //particle is before or behind cuff and out of range
                        continue;
                } else if (fabs(d_online) <= surfcuff->halflength) {
                //particle is inside cuff
			//fprintf(stderr,"particle is inside cuff\n"); //MF_DEBUG
                        d_toline = surfcuff->rcuff - d_toline;
                        dist = fabs(d_toline);
                } else {
                //particle is before or behind cuff but potentially within range
                        if (d_online < 0.) {    //before
                                d_online += surfcuff->halflength;
				//fprintf(stderr,"particle is before cuff\n"); //MF_DEBUG
                        } else {                //behind
                                d_online -= surfcuff->halflength;
				//fprintf(stderr,"particle is behind cuff\n"); //MF_DEBUG
                        }
                        tmp = surfcuff->rcuff - d_toline; //tmp is now the distance between particle and cylinder's mantle (if the cylinder were long enough...)
                        
                        for (idim = 0, d_toline2=0. ; idim < 3 ; idim++) {
                                if (idim == surfcuff->dir) {
                                        tolinedir[idim] = d_online;
                                }
                                d_toline2+=SQR(tolinedir[idim]); //update and overwrite d_toline
                        }
                        //tolinedir is now the direction of the force
                        d_toline = tmp;
                        dist = sqrt(d_toline2);
                }
                //we now have the distance and the direction of the force, so: apply force!
                //(and calculate force on the mantle while we do it...)
                frc = force_c(dist,t,surfcuff);
                if (!surfcuff->reponly[t] || frc>0.) {
			//fprintf(stderr,"applying force of %lg\n",frc); //MF_DEBUG
                        for (idim = 0 ; idim < 3 ; idim++) {
                                if (idim == surfcuff->dir) {
                                        f[i][idim] += frc*tolinedir[idim]/dist; //note that we have to normalize the direction to unit length
#ifdef EXTEND
                                        f_on_mantle -=frc*tolinedir[idim]/dist;
#endif
                                } else {
                                        f[i][idim] += sign*frc*tolinedir[idim]/dist; //note that we have to normalize the direction to unit length
                                        f_on_mantle -= sign*frc*tolinedir[idim]/dist;
                                }
                        }
                }
        }
	//sum up f_on_mantle contributions from different processors
	MPI_Allreduce(MPI_IN_PLACE, &f_on_mantle, 1, MPI_DOUBLE, MPI_SUM, comm_grid);
        //UPDATE RADIUS
        //we take into account that a small cylinder radius entails a small mantle area and thus a weaker force
        //for this we substract d0 = V0/F0+R0 from the cylinder radius - this is the distance at which the potential goes through zero
        double circ = 2.*M_PI*(surfcuff->rcuff-(surfcuff->V_R0[1]/surfcuff->F_R0[1]+R0)); //using headgroups particles as reference
	////BUT FOR NOW WE IGNORE THIS!!! AND USE THE ACTUAL RADIUS AS RADIUS INSTEAD...
        //double circ = 2.*M_PI*(surfcuff->rcuff/*-(surfcuff->V_R0[1]/surfcuff->F_R0[1]+R0)*/); //using headgroups particles as reference
        surfcuff->rcuff -= (surfcuff->F * circ - f_on_mantle) / surfcuff->m * b->dt ;
        if (surfcuff->rcuff < surfcuff->sigma*Rc) {
                surfcuff->rcuff = surfcuff->sigma*Rc;
                fprintf(stderr,"surfcuff: radius has become smaller than interaction range - setting radius equal to interaction range\n");
        }
#ifdef EXTEND
	surfcuff->halflength = surfcuff->halflength_times_r2 / (SQR(surfcuff->rcuff));
#endif
	//OUTPUT
        if(ismaster)
                fprintf(surfcuff->fh,"%lg %lg %lg %lg\n",b->time,surfcuff->rcuff,f_on_mantle,circ);
                //fprintf(stdout,"%lg %lg %lg %lg\n",b->time,surfcuff->rcuff,f_on_mantle,circ); //MF_DEBUG
}

void surfcuff_print_forces(struct beads *b, struct surfcuff *surfcuff)
{
        if (surfcuff == NULL) return;
	if (!ismaster) return; //this only needs to be done on one processor
	if (surfcuff->dir != 0) {
		fprintf(stderr, "Won't create plot of surfcuff potential, because the direction isn't x.\n");
		return;
	}

        FILE *surf_out = fopen("surfcuff.gp", "w");
        if (surf_out == NULL) fatal(EIO, "Couldn't open surfcuff.gp");
        fprintf(surf_out,"set style line 1 lt 1 lw 2 lc rgb \"black\"\n");
        fprintf(surf_out,"set style line 2 lt 1 lw 2 lc rgb \"red\"\n");
        fprintf(surf_out,"set style line 3 lt 1 lw 2 lc rgb \"blue\"\n");
        fprintf(surf_out,"set style line 4 lt 1 lw 2 lc rgb \"green\"\n");
        fprintf(surf_out,"set style line 5 lt 1 lw 2 lc rgb \"orange\"\n");
        fprintf(surf_out,"set style arrow 1 ls 1\n");
        fprintf(surf_out,"set style arrow 2 ls 2\n");
        fprintf(surf_out,"set style arrow 3 ls 3\n");
        fprintf(surf_out,"set style arrow 4 ls 4\n");
        fprintf(surf_out,"set style arrow 5 ls 5\n");

        int idim,t;
	const VEC lh = {.5 * b->l[0], .5 * b->l[1], .5 * b->l[2]};
        //const int dim[3] = {TANG1,TANG2,NORMAL};

        VEC tolinedir; //direction normal to the line running through the cuff's center
        double d_toline, d_toline2, d_online, frc, tmp, dist, sign;
        /****** list of distances ***************************************************************************************
         *      d_toline: initially length of tolinedir                                                                 *
         *                later it will become the distance to the cuff's mantle if the cuff was long enough            *
         *      dist: distance to the reference point for the force in this region                                      *
         *      d_online: initially the displacement along the direction of the cuff relative to the cuff's center,     *
         *                later it will be changed to be relative to one of the cuff's end                              *
         ****************************************************************************************************************/

        //int offline_ref_dir = (surfcuff->dir == 0) ? 1 : 0;

        //for (idim = 0 , i = 0 ; idim < 3 ; idim++) {
        //        if (idim != surfcuff->dir) {
        //                offline_dims[i] = idim;
        //                i++;
        //        }
        //}


        //double offline_ref = surfcuff->center[offline_ref_dir] + surfcuff->rcuff;
        //double offline_ref = surfcuff->center[offline_dims[0]] + surfcuff->rcuff; //either offline_dim would work...
        //oops, I guess we don't need this and it is wrong - should be only rcuff!

        //double online_ref = surfcuff->center[surfcuff->dir] + surfcuff->halflength; //probably not needed?


	//MF_TODO loop should be over some points, not particles
	double xxx=0, phi, theta=0.1, rad;
	int color,sgn;
	double coord[3], dr_xy[2];
	double tip[3];
	//MF in order to follow the cuff potential at a equipotential distance, we have to vary the x coordinate (xxx) when within the cuff
	//MF and the angle at the edge (theta) when next to the cuff
	//MF the while loop below does that, and the if ... else ... decides, along which parameter to travel
	while (/*xxx <= surfcuff->halflength &&*/ theta <= .5*M_PI) {
	  for (rad = .1*surfcuff->sigma ; rad <= surfcuff->rcuff+surfcuff->sigma*Rc ; rad+=.5*surfcuff->sigma) {
            for (phi = 0 ; phi < 2.*M_PI ; phi+=99999.5) { //angle around cuff axis
	     for (sgn = 1 ; sgn >= -1 ; sgn-=2) {
	      fprintf(stderr,"xxx:%lg rad:%lg phi:%lg sgn:%d theta:%lg\n",xxx,rad,phi,sgn,theta);
	      if (xxx < surfcuff->halflength) {
	        coord[0] = surfcuff->center[0] + sgn*xxx;
	        //xxx += .2;
		coord[1] = surfcuff->center[1] + sin(phi)*(-rad+surfcuff->rcuff);
		coord[2] = surfcuff->center[2] + cos(phi)*(-rad+surfcuff->rcuff);
	      } else { //angle along edge
		coord[1] = surfcuff->center[1] + sin(phi)*(surfcuff->rcuff);
		coord[2] = surfcuff->center[2] + cos(phi)*(surfcuff->rcuff);
		dr_xy[0] = sgn*rad*cos(theta); //create direction vector in xy plane
		dr_xy[1] = rad*sin(theta);
		coord[1] -= dr_xy[1] * sin(phi); //rotate direction vector around cuff axis
		coord[2] -= dr_xy[1] * cos(phi);
		coord[0] = dr_xy[0];
		coord[0] += surfcuff->center[0] + sgn * surfcuff->halflength; //add center of cuff
		//theta += .2;
	      }
	      //we have the coordinate now...
		//fprintf(stderr,"xxx:%lg rad:%lg phi:%lg sgn:%d theta:%lg coord:%lg %lg %lg\n",xxx,rad,phi,sgn,theta,coord[0],coord[1],coord[2]);
                t = 1;
                if (surfcuff->rhoC[t] == 0) continue;
		
                for (idim = 0, d_toline2=0. ; idim < 3 ; idim++) {
                        tmp = coord[idim] - surfcuff->center[idim];
                        tmp += (tmp < -lh[idim]) ? b->l[idim]: 0.; // PBC
                        tmp -= (tmp >  lh[idim]) ? b->l[idim]: 0.; // minimum distance MF
                        if (idim == surfcuff->dir) {
                                d_online = tmp;
                                tolinedir[idim] = 0.;
                        } else {
                                tolinedir[idim] = -tmp;
                                d_toline2 += SQR(tmp);
                        }
                }
                d_toline = sqrt(d_toline2);
		for (idim = 0 ; idim < 3 ; idim++) { //d_toline should be direction to mantle, not to center line
			tolinedir[idim] *= (surfcuff->rcuff-d_toline)/d_toline;
		}
		//fprintf(stderr,"d_online:%lg d_toline:%lg tolinedir:%lg %lg %lg\n",d_online,d_toline,tolinedir[0],tolinedir[1],tolinedir[2]);
                sign = 1.; //standard assumption is that the particle is inside or next to the cuff
                if (d_toline > surfcuff->rcuff)
                        sign = -1.; //bead is somehow outside of cuff resulting into a different force direction...

                if (fabs(d_online) >= surfcuff->halflength + Rc*surfcuff->sigma) {
                //particle is before or behind cuff and out of range
                	fprintf(stderr,"particle is before or behind cuff and out of range\n"); //MF_DEBUG
                        continue;
                } else if (fabs(d_online) <= surfcuff->halflength) {
                //particle is inside cuff
			fprintf(stderr,"particle is inside cuff\n"); //MF_DEBUG
                        d_toline = surfcuff->rcuff - d_toline;
			color = 1;
                        dist = fabs(d_toline);
                } else {
                //particle is before or behind cuff but potentially within range
                        if (d_online < 0.) {    //before
                                d_online += surfcuff->halflength;
				color = 2;
				fprintf(stderr,"particle is before cuff\n"); //MF_DEBUG
                        } else {                //behind
                                d_online -= surfcuff->halflength;
				color = 3;
				fprintf(stderr,"particle is behind cuff\n"); //MF_DEBUG
                        }
                        tmp = surfcuff->rcuff - d_toline; //tmp is now the distance between particle and cylinder's mantle (if the cylinder were long enough...)
                        
                        for (idim = 0, d_toline2=0. ; idim < 3 ; idim++) {
                                if (idim == surfcuff->dir) {
                                        tolinedir[idim] = d_online;
                                //} else {
                                //        tolinedir[idim] *= tmp; //give vector a length corresponding to above mentioned distance
                                }
                                d_toline2+=SQR(tolinedir[idim]); //update and overwrite d_toline
                        }
                        //tolinedir is now the direction of the force
                        d_toline = tmp;
                        dist = sqrt(d_toline2);
                }
                //we now have the distance and the direction of the force, so: apply force!
                //(and calculate force on the mantle while we do it...)
                frc = force_c(dist,t,surfcuff);
                if (!surfcuff->reponly[t] || frc>0.) {
                        for (idim = 0 ; idim < 3 ; idim++) {
				if (idim == surfcuff->dir) {
                                	tip[idim] = coord[idim] + frc*tolinedir[idim]/dist; //note that we have to normalize the direction to unit length
				} else {
                                	tip[idim] = coord[idim] + sign*frc*tolinedir[idim]/dist; //note that we have to normalize the direction to unit length
				}
                        }
                }
                if (frc!=0.) {
                        fprintf(surf_out,"set arrow from %lg, %lg, %lg to %lg, %lg, %lg as %d\n",coord[0],coord[1],coord[2],tip[0],tip[1],tip[2],color);
		}
	      }
            }     
	  }
	  if (xxx < surfcuff->halflength) {
	    xxx += .2;
	    fprintf(stderr,"increasing xxx\n");
	  } else {
	    theta += .2;
	    fprintf(stderr,"increasing theta\n");
	  }
        }
	//finish gnuplot file by plotting an arrow along the axis of the cylinder
	color = 4;
        fprintf(surf_out,"set arrow from %lg, %lg, %lg to %lg, %lg, %lg as %d\n",
	surfcuff->center[0]-surfcuff->halflength,surfcuff->center[1],surfcuff->center[2],
	surfcuff->center[0]+surfcuff->halflength,surfcuff->center[1],surfcuff->center[2],color);
	fprintf(surf_out,"set xrange [0:%lg]\n",b->l[0]);
	fprintf(surf_out,"set yrange [0:%lg]\n",b->l[1]);
	fprintf(surf_out,"set zrange [0:%lg]\n",b->l[2]);
	fprintf(surf_out,"splot 0\n");

        fclose(surf_out);
}

void print_surfbump_forces(struct beads *b, struct surfbump *surfbump)
{

        int i,idim;
	const VEC lh = {.5 * b->l[0], .5 * b->l[1], .5 * b->l[2]};
        const int dim[3] = {TANG1,TANG2,NORMAL};
        VEC XZdir;
        double base_curv_ref[TYPE_MAX][3];
        double top_curv_ref[TYPE_MAX][3];
        for (i = 0 ; i < TYPE_MAX ; i++) {
                base_curv_ref[i][0] = surfbump->Rbump + surfbump->sigma[i]*Rc;
                base_curv_ref[i][1] = 0.;
                base_curv_ref[i][2] =  surfbump->sigma[i]*Rc;
                top_curv_ref[i][0] = surfbump->Rbump;
                top_curv_ref[i][1] = 0.;
                top_curv_ref[i][2] = surfbump->H;
        }
        VEC XYdir;
        double d_xy, d_xy2, vecl, vecl2;

        VEC coord,fvec;
        FILE *surf_out = fopen("surfbump.gp", "w");
        if (surf_out == NULL) fatal(EIO, "Couldn't open surfbump.gp");
        fprintf(surf_out,"set style line 1 lt 1 lw 2 lc rgb \"black\"\n");
        fprintf(surf_out,"set style line 2 lt 1 lw 2 lc rgb \"red\"\n");
        fprintf(surf_out,"set style line 3 lt 1 lw 2 lc rgb \"blue\"\n");
        fprintf(surf_out,"set style line 4 lt 1 lw 2 lc rgb \"green\"\n");
        fprintf(surf_out,"set style line 5 lt 1 lw 2 lc rgb \"orange\"\n");
        fprintf(surf_out,"set style arrow 1 ls 1\n");
        fprintf(surf_out,"set style arrow 2 ls 2\n");
        fprintf(surf_out,"set style arrow 3 ls 3\n");
        fprintf(surf_out,"set style arrow 4 ls 4\n");
        fprintf(surf_out,"set style arrow 5 ls 5\n");

        for (coord[dim[0]]= surfbump->Pos[0]-surfbump->Rbump-surfbump->sigma[1]*Rc*2. ; coord[0]< surfbump->Pos[0]+surfbump->Rbump+surfbump->sigma[1]*Rc*2.; coord[0]+=0.2) {
//          for (coord[dim[1]]= surfbump->Pos[1]-surfbump->Rbump-surfbump->sigma[1]*Rc*2. ; coord[1]< surfbump->Pos[1]+surfbump->Rbump+surfbump->sigma[1]*Rc*2.; coord[1]+=2.0) {
            for (coord[dim[2]]= surfbump->Vpos[1]-3 ; coord[2] < surfbump->Vpos[1]+surfbump->H+surfbump->sigma[1]*Rc*2. ; coord[2]+=0.2 ) {

                coord[dim[1]] = surfbump->Pos[1];
                //coord[dim[2]] = surfbump->Vpos[1];

                double frc = 0.;
                int color = 1;
                for (idim=0;idim<3;idim++) {
                //        coord[idim]=drand48()*b->l[idim];
                        fvec[idim]=0.;
                }
                int t = 1;
                if (surfbump->rhoC[t] == 0) continue;
                for (idim = 0, d_xy2 = 0. ; idim < 2 ; idim++) {
                        XYdir[idim] = coord[dim[idim]] - b->surfbump->Pos[idim];
                        XYdir[idim] += (XYdir[idim] < -lh[dim[idim]]) ? b->l[dim[idim]]: 0.; // PBC
                        XYdir[idim] -= (XYdir[idim] >  lh[dim[idim]]) ? b->l[dim[idim]]: 0.; // minimum distance MF
                        d_xy2 += SQR(XYdir[idim]);
                }
                d_xy = sqrt(d_xy2);
                //XYdir[2] = coord[dim[2]] - surfbump->Vpos[t];




                if (d_xy >= surfbump->sigma[t]*Rc + surfbump->Rbump) {
                        // particle is not within horizontal range of the bump
                        color = 1;
                        double d = (coord[NORMAL]-surfbump->Vpos[t])/surfbump->sigma[t];
                        frc = force_b(d, t, surfbump);
                        if (!surfbump->reponly[t] || frc>0) {
                                fvec[NORMAL] += frc;
                        } 
                } else if (d_xy <= surfbump->Rbump) {
                        // particle is above flat top of bump
                        color = 2;
                        double d = (coord[NORMAL] - surfbump->Vpos[t] - surfbump->H)/surfbump->sigma[t];
                        frc = force_b(d, t, surfbump);
                        if (!surfbump->reponly[t] || frc>0) {
                                fvec[NORMAL] += frc;
                        } 
                } else {
                        if (coord[NORMAL] - surfbump->Vpos[t] < surfbump->sigma[t]*Rc) {
                                // particle in curved region at base of bump
                                // we first deal with the particle's projection into the T1,N plane (e.g. XZ)
                                // for this we can use the radial distance we already know and x[N]
                                color = 3;
                                VEC XZcoord = {d_xy, 0., coord[NORMAL]-surfbump->Vpos[t]};
                                for (idim = 0, vecl2 = 0. ; idim < 3 ; idim++) {
                                        XZdir[idim] = base_curv_ref[t][idim] - XZcoord[idim];
                                        vecl2 += SQR(XZdir[idim]);
                                }
                                vecl = sqrt(vecl2);
                                for (idim = 0 ; idim < 3 ; idim++) {
                                        XZdir[idim] /= vecl;
                                } // now we have the planar direction (unit vector pointing away from surface) and the distance from the reference
                                //double d = surfbump->sigma[t]*Rc - vecl;
                                double d = Rc - vecl / surfbump->sigma[t]; // MFnew CHECK!!!
                                frc = force_b(d, t, surfbump);
                                if (!surfbump->reponly[t] || frc>0) {
                                        VEC dir;
                                        for (idim = 0 ; idim < 2 ; idim++) { //scale XYdir so that it has a length corresponding to the X component of XZdir
                                                dir[idim] = XYdir[idim] / d_xy * XZdir[0];
                                        }
                                        dir[2] = XZdir[2];
                                        for (idim = 0 ; idim < 3 ; idim++) {
                                                fvec[dim[idim]] += frc * dir[idim];
                                        }
                                }
                        } else if (coord[NORMAL]-surfbump->Vpos[t] <= b->surfbump->H) {
                                // particle next to vertical side of bump
                                // XYdir has already been calculated while measuring d_xy
                                color = 4;
                                double d = (d_xy - surfbump->Rbump) / surfbump->sigma[t];
                                frc = force_b(d, t, surfbump);
                                if (!surfbump->reponly[t] || frc>0) {
                                        for (idim = 0 ; idim < 2 ; idim++) {
                                                XYdir[idim] /= d_xy;
                                                fvec[dim[idim]] += frc * XYdir[idim]; 
                                        }
                                }
                        } else if (coord[NORMAL]-surfbump->Vpos[t] <  b->surfbump->H + surfbump->sigma[t]*Rc) {
                                // particle next to round tip of bump (but not above the flat top)
                                // we first deal with the particle's projection into the T1,N plane (e.g. XZ)
                                // for this we can use the radial distance we already know and x[N]
                                color = 5;
                                VEC XZcoord = {d_xy, 0., coord[NORMAL]-surfbump->Vpos[t]};
                                for (idim = 0, vecl2 = 0. ; idim < 3 ; idim++) {
                                        XZdir[idim] = XZcoord[idim] - top_curv_ref[t][idim];
                                        vecl2 += SQR(XZdir[idim]);
                                }
                                vecl = sqrt(vecl2);
                                for (idim = 0 ; idim < 3 ; idim++) {
                                        XZdir[idim] /= vecl;
                                } // now we have the planar direction (unit vector pointing away from surface) and the distance from the reference
                                //double d = vecl - surfbump->sigma[t]*Rc;
                                double d = vecl / surfbump->sigma[t] ; //MFnew CHECK!!!
                                frc = force_b(d, t, surfbump);
                                //fprintf(stderr,"d= %lg frc= %lg\n",d,frc);
                                if (!surfbump->reponly[t] || frc>0) {
                                        VEC dir;
                                        for (idim = 0 ; idim < 2 ; idim++) { //scale XYdir so that it has a length corresponding to the X component of XZdir
                                                dir[idim] = XYdir[idim] / d_xy * XZdir[0];
                                        }
                                        dir[2] = XZdir[2];
                                        for (idim = 0 ; idim < 3 ; idim++) {
                                                fvec[dim[idim]] += frc * dir[idim];
                                        }
                                }
                        }
                }
                if (frc!=0.) {
                        fprintf(surf_out,"set arrow from %lg, %lg, %lg to %lg, %lg, %lg as %d\n",coord[0],coord[1],coord[2],coord[0]+fvec[0],coord[1]+fvec[1],coord[2]+fvec[2],color);
                        //i++;
                }
//            }
          }
        }
        //fprintf(surf_out,"splot [0:%lg] [0:%lg] [0:%lg] %lg\n",b->l[dim[0]],b->l[dim[1]],b->l[dim[2]],surfbump->Vpos[1]);
        fprintf(surf_out,"splot [%lg:%lg] [%lg:%lg] [%lg:%lg] %lg\n",
        surfbump->Pos[0]-surfbump->Rbump-surfbump->sigma[1]*Rc*2., surfbump->Pos[0]+surfbump->Rbump+surfbump->sigma[1]*Rc*2.,
        surfbump->Pos[1]-surfbump->Rbump-surfbump->sigma[1]*Rc*2., surfbump->Pos[1]+surfbump->Rbump+surfbump->sigma[1]*Rc*2.,
        surfbump->Vpos[1]-3, surfbump->Vpos[1]+surfbump->H+surfbump->sigma[1]*Rc*2.,surfbump->Vpos[1]);
        fclose(surf_out);
}



static void print_surfforce(struct surfpot *surfpot)
{
        FILE *surf_out = fopen("surfforce.dat", "w");
        if (surf_out == NULL) fatal(EIO, "Couldn't open surfforce.dat");

        int t = 1;
        double z;

        for (z = 0; z < (Rc+1)*surfpot->sigma[t]; z+=surfpot->sigma[t]/100) {
                //fprintf (surf_out,"%lg %lg\n",z+surfpot->Vpos[t],surfpot->C93[t]*force(z/surfpot->sigma[t],t,surfpot));
                fprintf (surf_out,"%lg %lg\n",z,surfpot->C93[t]*force(z/surfpot->sigma[t],t,surfpot));
        }
        fclose(surf_out);
}

static void print_surfpot(struct surfpot *surfpot)
{
        FILE *surf_out = fopen("surfpot.dat", "w");
        if (surf_out == NULL) fatal(EIO, "Couldn't open surfpot.dat");

        int t = 1;
        double z;

        for (z = 0; z < (Rc+1)*surfpot->sigma[t]; z+=surfpot->sigma[t]/100) {
                //fprintf (surf_out,"%lg %lg\n",z+surfpot->Vpos[t],surfpot->C93[t]*pot(z/surfpot->sigma[t],t,surfpot));
                fprintf (surf_out,"%lg %lg\n",z,surfpot->C93[t]*pot(z/surfpot->sigma[t],t,surfpot));
        }
        fclose(surf_out);
}

void surfpot_free(struct surfpot *surfpot)
{
        if (surfpot == NULL) return;
        free(surfpot);
}

void surfbump_free(struct surfbump *surfbump)
{
        if (surfbump == NULL) return;
        free(surfbump);
}

void surfcuff_free(struct surfcuff *surfcuff)
{
        if (surfcuff == NULL) return;
        if(ismaster)
                fclose(surfcuff->fh);
        free(surfcuff);
}

void surfcuff_reset_beads(struct beads *b, struct surfcuff *surfcuff) {
	//Should be called in advance, so we can have random coordinates in the starting configuration!
        if (surfcuff == NULL) return;
        int i, j, k;
        int idx = 0;

        /* b->local_n is the number of molecules which belong to this
         * MPI process */ 
        for (i = 0; i < b->local_n; i++) {
                /* b->local_n_b[i] is the block index of molecule i */
                int thisblock = b->local_n_b[i];


                /* the total number of molecules in this block in all
                 * processes can be found in b->n[thisblock]. */

                int N = b->N[thisblock];
                if (thisblock == surfcuff->block) { //we found the surfcuff block
                	/* b->N[] is an array with b->blocks entries, specifying
                	 * the number of beads per molecule in this block */
                        for (j = idx; j < idx + N; j++) {
                                for (k = 0; k < 3; k++) {
                                        /* the unfolded, i.e., absolute,
                                         * coordinates of this bead are located
                                         * in b->x_intra[][]. */
                                        //double x = b->x_intra[idx][k];
                                        double x = b->x_intra[j][k]; //MF correction

                                        /* b->f_intra contains the forces on
                                         * the beads belonging to the local
                                         * process. Each bead exists exactly
                                         * in one of all processes in
                                         * f_intra[][]. */
                                        //b->f_intra[idx][k] += .1 * x;
                                        b->f_intra[j][k] += .1 * x; //MF correction


					//MF I GUESS I SHOULD UPDATE THE POSITIONS IN b->xv,
					//   THE POSITIONS IN b->x_intra are calculated in set_xf_intra according to b->xv
					//   AND WOULD OTHERWISE OVERWRITE OUR RESET POSITIONS AGAIN...
					//   also set velocities to 0 for good measure

					//   think if we need to do the periodic image tracking...
                                }
                        }
                }
                idx += N;
        }
        return 0.;

}
