#include "fdmdpd2.h"

struct usp
{
        int blocks[2];
        double k;
        double x0s[4]; //MF_WARNING only works for restraining in plane at the moment
        double coms[6];
        //double comsMI2d[4];
        double alpha;
        /* missing are Energy, displacement */
        //int newround;
        double tt[2][9]; //Traegheitstensor
        FILE* fh;
};


//        static const double Rc = 3.; //Cutoff set to 3xSigma
//        static const double R1 = 2.; //R1 set to 2xSigma
//        static const double R0 = 1.11; //R0 set to 1.11xSigma (truncates too high forces...)


int eigval(double A[9], double w[3]);


void usp_print_header(FILE *FH, struct usp *usp)
{
        if (usp == NULL) return;

        int i;

        fprintf(FH, "# usp ");
        fprintf(FH, "blocks= ");
        for (i=0 ; i<2 ; i++) {
                fprintf(FH, "%d ",usp->blocks[i]);
        }
        fprintf(FH, "xa0= ");
        for (i=0 ; i<2 ; i++) {
                fprintf(FH, "%lg ",usp->x0s[i]);
        }
        fprintf(FH, "xb0= ");
        for (i=0 ; i<2 ; i++) {
                fprintf(FH, "%lg ",usp->x0s[2+i]);
        }
        fprintf(FH, "k= %lg",usp->k);
        fprintf(FH, "\n");

}


void usp_read_header(FILE *FH, struct beads *b)
{
        int i;

        b->usp = malloc(sizeof(*b->usp));
        if (b->usp == NULL) novm("usp");

        printf("Restraining two blocks with umbrella potentials (usp).\n");

	fscanf(FH, "# usp blocks= ");
	for (i = 0; i < 2 ; i++) {
		if (fscanf(FH, "%d ", &b->usp->blocks[i]) != 1)
			fatal(EINVAL, "usp blocks!");
                if (b->usp->blocks[i] >= b->blocks)
                        fatal(EINVAL, "usp - block doesn't exist!");
	}
	fscanf(FH, "xa0= ");
	for (i = 0; i < 2 ; i++) {
		if (fscanf(FH, "%lg ", &b->usp->x0s[i]) != 1)
			fatal(EINVAL, "usp xa0!");
	}
	fscanf(FH, "xb0= ");
	for (i = 0; i < 2 ; i++) {
		if (fscanf(FH, "%lg ", &b->usp->x0s[2+i]) != 1)
			fatal(EINVAL, "usp xb0!");
	}
	fscanf(FH, "k= ");
        if (fscanf(FH, "%lg ", &b->usp->k) != 1)
			fatal(EINVAL, "usp k!");


        if (ismaster) {
                //MF_TODO use consecutive filenames by adding the stepnumber/startingtime into the name
	        char fn[PATH_MAX];
                sprintf(fn, "usp%09lu.dat", b->step); // MF BUSY HERE
                b->usp->fh = fopen(fn,"w");
                if (b->usp->fh == NULL)
                        fatal(EIO, "Couldn't open usp.dat for writing");
                fprintf(b->usp->fh,"# blocks %d and %d were restrained to %lg %lg and %lg %lg in the tangential plane with k %lg\n",
                                b->usp->blocks[0],b->usp->blocks[1],b->usp->x0s[0],b->usp->x0s[1],b->usp->x0s[2],b->usp->x0s[3],b->usp->k);
                fprintf(b->usp->fh,"# time, T1/2/N pos of block 1, orientation of block 1, T1/2/N pos of block 2, orientation of block 2\n");
        }
}


void usp_get_coms(struct beads *b, struct usp *usp)
{
        int i, j, k;
        int idx = 0;

        /* b->local_n is the number of molecules which belong to this
         * MPI process */ 
        const int dim[3] = {TANG1,TANG2,NORMAL};
        for (j=0 ; j<6 ; j++)
                usp->coms[j]=0.;
	for (i = 0; i < b->local_n; i++) { // GET COMS
                //FOR EACH ATOM CHECK IF THE BLOCK IS SUBJECT TO THE UMBRELLA POTENTIAL
                int thisblock = b->local_n_b[i];
                /* b->N[] is an array with b->blocks entries, specifying
                 * the number of beads per molecule in this block */
		int N = b->N[thisblock];
                if (thisblock==usp->blocks[0] || thisblock==usp->blocks[1]) {
                        int blck=1; //IDENTIFY BLOCK AS 2nd OR 1st OF THE SPECIFIED BLOCKS IN THE HEADER
                        if (thisblock==usp->blocks[0])
                                blck=0;

                        //fprintf(stderr,"found block %d with number %d\n",blck,usp->blocks[blck]);


                        /* the total number of molecules in this block in all
                         * processes can be found in b->n[thisblock]. */

		        for (j = idx; j < idx + N; j++) {
		        	for (k = 0; k < 3; k++) {
                                        /* the unfolded, i.e., absolute,
                                         * coordinates of this bead are located
                                         * in b->x_intra[][]. */

                                        //THE COORDINATES ARE ADDED INTO THE COMS
                                        usp->coms[3*blck+k] += b->x_intra[j][dim[k]];
                                        //fprintf(stderr,"com[%d]=%lg, x[i]=%lg, dim[i]=%d\n",k,usp->coms[3*blck+k],b->x_intra[j][dim[k]],dim[k]);

		        	}
		        }
                }
		idx += N;
	}
	MPI_Allreduce(MPI_IN_PLACE, usp->coms, 6, MPI_DOUBLE, MPI_SUM, comm_grid);
        for (i=0 ; i<2 ; i++) //blocks //NORMALIZE COMS
                for (j=0 ; j<3 ; j++) {//dims
                        usp->coms[3*i+j] /= b->n[usp->blocks[i]] * b->N[usp->blocks[i]];
                        //fprintf(stderr,"block:%d dim:%d com:%lg based on n:%d times N:%d\n",i,j,usp->coms[3*i+j],b->n[usp->blocks[i]],b->N[usp->blocks[i]]);
                }
}

void usp_apply_forces(struct beads *b, struct usp *usp)
{
        int i, j, k;
        int idx = 0;
        //double fdir[4];
        //double f[2];
        //double dist[2];
        double f[4];
        double comsMI2d[4];

        const int dim[3] = {TANG1,TANG2,NORMAL};

        //CALCULATE FORCES
        //calculate distance. for this, the minimum image convention has to be applied.
        for (i=0 ; i<2 ; i++) {//dimension, 2d
                for (j=0 ; j<2 ; j++) {//restrained block
                        comsMI2d[2*j+i] = usp->coms[3*j+i];
                        //for (k=-1 ; k<=1 ; k+=2) {//sign
                                k = (comsMI2d[2*j+i] - usp->x0s[2*j+i] > 0) ? -1 : 1;
                                //fprintf(stderr,"distance calculation i:%d j:%d k:%d\n",i,j,k);
                                //fprintf(stderr,"comsMI2d:%lg x0s:%lg k*b->l:%lg\n",comsMI2d[2*j+i],usp->x0s[2*j+i],k*b->l[dim[i]]);
                                while ( fabs(comsMI2d[2*j+i] + k*b->l[dim[i]] - usp->x0s[2*j+i]) < fabs(comsMI2d[2*j+i] - usp->x0s[2*j+i]) ) {
                                        comsMI2d[2*j+i] += k*b->l[dim[i]];
                                }
                        //}
                }
        }
        //calculate force
        for (j=0 ; j<2 ; j++) { //restrained block
                for (i=0 ; i<2 ; i++) { //dimension, 2d
                        f[2*j+i] = -usp->k * (comsMI2d[2*j+i] - usp->x0s[2*j+i]);
                }
        }


        //APPLY FORCES.
	for (i = 0; i < b->local_n; i++) {
                //FOR EACH ATOM CHECK IF THE BLOCK IS SUBJECT TO THE UMBRELLA POTENTIAL
                int thisblock = b->local_n_b[i];
		int N = b->N[thisblock];
                if (thisblock==usp->blocks[0] || thisblock==usp->blocks[1]) {
                        int blck=1; //IDENTIFY BLOCK AS 2nd OR 1st OF THE SPECIFIED BLOCKS IN THE HEADER
                        if (thisblock==usp->blocks[0])
                                blck=0;


		        for (j = idx; j < idx + N; j++) {
		        	for (k = 0; k < 2; k++) { //dimension, 2d
                                        //apply forces
					b->f_intra[j][k] += f[2*blck+k];

		        	}
		        }
                }
		idx += N;
	}
}

void usp_eval_print(struct beads *b, struct usp *usp)
{
        if (usp == NULL) return;
        int i, j, k;
        int idx = 0;
        double xrel[3];

        const int dim[3] = {TANG1,TANG2,NORMAL};

        for (j=0 ; j<2 ; j++)
                for (k=0 ; k<9 ; k++)
                        usp->tt[j][k]=0.;
        //DETERMINE ORIENTATION
	for (i = 0; i < b->local_n; i++) {
                //FOR EACH ATOM CHECK IF THE BLOCK IS SUBJECT TO THE UMBRELLA POTENTIAL
                int thisblock = b->local_n_b[i];
		int N = b->N[thisblock];
                if (thisblock==usp->blocks[0] || thisblock==usp->blocks[1]) {
                        int blck=1; //IDENTIFY BLOCK AS 2nd OR 1st OF THE SPECIFIED BLOCKS IN THE HEADER
                        if (thisblock==usp->blocks[0])
                                blck=0;


		        for (j = idx; j < idx + N; j++) {
                                for (k=0 ; k<3 ; k++) { //dimension, 3d
                                        xrel[k] = b->x_intra[j][dim[k]] - usp->coms[3*blck+k];  // Transformation to relative coordinates
                                        //fprintf(stderr,"block:%d bead:%d xrel[%d]=%lg\n",blck,j,k,xrel[k]);
                                }

                                usp->tt[blck][0*3+0] += ((xrel[1]*xrel[1]) + (xrel[2]*xrel[2]));   // Traegheitstensor
                                usp->tt[blck][0*3+1] -= (xrel[0]*xrel[1]);
                                usp->tt[blck][0*3+2] -= (xrel[0]*xrel[2]);

                                usp->tt[blck][1*3+1] += ((xrel[0]*xrel[0]) + (xrel[2]*xrel[2]));
                                usp->tt[blck][1*3+2] -= (xrel[1]*xrel[2]);

                                usp->tt[blck][2*3+2] += ((xrel[0]*xrel[0]) + (xrel[1]*xrel[1]));

		        }
                }
		idx += N;
	}
        if(ismaster)
                fprintf(usp->fh,"%lg",b->time);

        for (k=0 ; k<2 ; k++) { //blocks
	        MPI_Allreduce(MPI_IN_PLACE, usp->tt[k], 9, MPI_DOUBLE, MPI_SUM, comm_grid);

                usp->tt[k][1*3+0]  = usp->tt[k][0*3+1]; //symmetry
                usp->tt[k][2*3+0]  = usp->tt[k][0*3+2];
                usp->tt[k][2*3+1]  = usp->tt[k][1*3+2]; //tt[k] is now complete!

                //for (i=0 ; i<3 ; i++) {
                //        for (j=0 ; j<3 ; j++) {
                //                fprintf(stderr,"block:%d tt[%d][%d]:%lg\n",k,i,j,usp->tt[k][3*i+j]);
                //        }
                //}


                // Diagonalisieren des Tr??heitstensors durch L??ung einer kubischen Gleichung
                // kubische cardanische Gleichung:
                // x^3 + ra*x^2 + sa*x + ta = 0

                /*
                double ra,sa,ta,p,q,rho,phi;
                ra = -(usp->tt[k][0*3+0] + usp->tt[k][1*3+1] + usp->tt[k][2*3+2]);
                sa = (usp->tt[k][0*3+0] * (usp->tt[k][1*3+1] + usp->tt[k][2*3+2])
                    - usp->tt[k][1*3+0]*usp->tt[k][0*3+1] - usp->tt[k][2*3+0]*usp->tt[k][0*3+2]
                    + usp->tt[k][1*3+1]*usp->tt[k][2*3+2] - usp->tt[k][2*3+1]*usp->tt[k][1*3+2]);
                ta = -(usp->tt[k][1*3+0]*(usp->tt[k][2*3+1]*usp->tt[k][0*3+2] - usp->tt[k][0*3+1]*usp->tt[k][2*3+2]) 
                     + usp->tt[k][2*3+0]*(usp->tt[k][0*3+1]*usp->tt[k][1*3+2] - usp->tt[k][1*3+1]*usp->tt[k][0*3+2])
                     + usp->tt[k][0*3+0]*(usp->tt[k][1*3+1]*usp->tt[k][2*3+2] - usp->tt[k][2*3+1]*usp->tt[k][1*3+2]) );

                // reduzierte Form der kubischen Gleichung y^3 + p*y + q = 0
                // rho = sqrt(-p^3/27.)  cos phi = -q/(2*rho)

                p = (3.*sa - ra*ra)/3.;
                q = ((2.*ra*ra*ra)/27.) - (ra*sa/3.) + ta;

                rho = sqrt(-(p*p*p)/27.);
                phi = acos(-q/(2.*rho));

                //fprintf(stderr,"ra:%lg sa:%lg ta:%lg p:%lg q:%lg rho:%lg phi:%lg\n",ra,sa,ta,p,q,rho,phi);

                double y[3];
                double lambda[3];
                double maxlambda[3];
                double evec[3];
                y[0] = 2.*pow(rho,1./3.)*cos(phi/3.);
                y[1] = 2.*pow(rho,1./3.)*cos(phi/3. + 2.*M_PI/3.);
                y[2] = 2.*pow(rho,1./3.)*cos(phi/3. + 4.*M_PI/3.);
                */

                double evec[3];
                double lambda2[3];
                double maxlambda[3];

                eigval(usp->tt[k],lambda2); // get Eigenvalues

                maxlambda[0] = MAX(lambda2[0], MAX(lambda2[1],lambda2[2])); // groesstes Tragheitsmoment

                if (maxlambda[0] == lambda2[0])
                {
                    maxlambda[1] = MAX(lambda2[1],lambda2[2]);
                    maxlambda[2] = MIN(lambda2[1],lambda2[2]); 
                }
                else if (maxlambda[0] == lambda2[1])
                {
                    maxlambda[1] = MAX(lambda2[0],lambda2[2]);
                    maxlambda[2] = MIN(lambda2[0],lambda2[2]); 
                }
                else
                {
                    maxlambda[1] = MAX(lambda2[0],lambda2[1]);
                    maxlambda[2] = MIN(lambda2[0],lambda2[1]);
                }

                //fprintf(stderr,"block %d has eigenvalues in descending order %lg %lg %lg\n",k,maxlambda[0],maxlambda[1],maxlambda[2]);
                //for (i=0;i<3;i++) {
                //        for (j=0;j<3;j++) {
                //                fprintf(stderr,"%lg ",usp->tt[k][i*3+j]);
                //        }
                //        fprintf(stderr,"\n");
                //}



                //fprintf(stderr,"y:%lg %lg %lg\n",y[0],y[1],y[2]);

                /*
                for(j=0; j<3; j++)
                        lambda[j] = y[j] - ra/3.; // Eigenwerte (Haupttraegheitsmomente)

                //fprintf(stderr,"lambda:%lg %lg %lg\n",lambda[0],lambda[1],lambda[2]);

                maxlambda[0] = MAX(lambda[0], MAX(lambda[1],lambda[2])); // groesstes Tragheitsmoment

                    if (maxlambda[0] == lambda[0])
                    {
                  	maxlambda[1] = MAX(lambda[1],lambda[2]);
                  	maxlambda[2] = MIN(lambda[1],lambda[2]); 
                    }
                    else if (maxlambda[0] == lambda[1])
                    {
                  	maxlambda[1] = MAX(lambda[0],lambda[2]);
                  	maxlambda[2] = MIN(lambda[0],lambda[2]); 
                    }
                    else
                    {
                  	maxlambda[1] = MAX(lambda[0],lambda[1]);
                  	maxlambda[2] = MIN(lambda[0],lambda[1]);
                    }

                */

                //fprintf(stderr,"maxlambda:%lg %lg %lg\n",maxlambda[0],maxlambda[1],maxlambda[2]);

                // Bestimmung des Eigenvektors zum gr??ten Eigenwert, e.g. oblates Objekt
                //evec[2] = 1.;
                //evec[1] = ((usp->tt[k][0*3+0]-maxlambda[0])*usp->tt[k][1*3+2]*evec[2] - usp->tt[k][1*3+0]*usp->tt[k][0*3+2]*evec[2]);
                //evec[1] = evec[1] / (usp->tt[k][1*3+0]*usp->tt[k][0*3+1] + (maxlambda[0]-usp->tt[k][0*3+0])*(usp->tt[k][1*3+1]-maxlambda[0]));
                //evec[0] = (usp->tt[k][0*3+1]*evec[1] + usp->tt[k][0*3+2]*evec[2]) / (maxlambda[0]-usp->tt[k][0*3+0]);

                // Bestimmung des Eigenvektors zum kleinsten Eigenwert, e.g. prolates Objekt
                evec[2] = 1.;
                evec[1] = ((usp->tt[k][0*3+0]-maxlambda[2])*usp->tt[k][1*3+2]*evec[2] - usp->tt[k][1*3+0]*usp->tt[k][0*3+2]*evec[2]);
                evec[1] = evec[1] / (usp->tt[k][1*3+0]*usp->tt[k][0*3+1] + (maxlambda[2]-usp->tt[k][0*3+0])*(usp->tt[k][1*3+1]-maxlambda[2]));
                evec[0] = (usp->tt[k][0*3+1]*evec[1] + usp->tt[k][0*3+2]*evec[2]) / (maxlambda[2]-usp->tt[k][0*3+0]);

                double xlen = sqrt(evec[0]*evec[0]+evec[1]*evec[1]+evec[2]*evec[2]);

                for(j=0; j<3; j++)
                        evec[j] /= xlen;    // normierter Eigenvektor = Richtungsvektor der Haupttr??heitsachse

                //OUTPUT
                if(ismaster)
                        fprintf(usp->fh," %lg %lg %lg %lg %lg %lg",usp->coms[3*k],usp->coms[3*k+1],usp->coms[3*k+2],evec[0],evec[1],evec[2]);

        }
        if(ismaster)
                fprintf(usp->fh,"\n");
}

void usp_wrapper(struct beads *b, struct usp *usp)
{
        if (usp == NULL) return;
        usp_get_coms(b,b->usp);
        usp_apply_forces(b,b->usp);
}

void usp_free(struct usp *usp)
{
        if (usp == NULL) return;
        if(ismaster)
                fclose(usp->fh);
        free(usp);
}


int eigval(double A[9], double w[3])
// Calculates the eigenvalues of a symmetric 3x3 matrix A using Cardano's
// analytical algorithm.
// Only the diagonal and upper triangular parts of A are accessed. The access
// is read-only.
// ----------------------------------------------------------------------------
// Parameters:
//   A: The symmetric input matrix
//   w: Storage buffer for eigenvalues

// Constants
#define M_SQRT3    1.73205080756887729352744634151   // sqrt(3)

{
  double m, c1, c0;
  
  // Determine coefficients of characteristic poynomial. We write
  //       | a   d   f  |
  //  A =  | d*  b   e  |
  //       | f*  e*  c  |
  double de = A[0*3+1] * A[1*3+2];                                    // d * e
  double dd = SQR(A[0*3+1]);                                         // d^2
  double ee = SQR(A[1*3+2]);                                         // e^2
  double ff = SQR(A[0*3+2]);                                         // f^2
  m  = A[0*3+0] + A[1*3+1] + A[2*3+2];
  c1 = (A[0*3+0]*A[1*3+1] + A[0*3+0]*A[2*3+2] + A[1*3+1]*A[2*3+2])        // a*b + a*c + b*c - d^2 - e^2 - f^2
          - (dd + ee + ff);
  c0 = A[2*3+2]*dd + A[0*3+0]*ee + A[1*3+1]*ff - A[0*3+0]*A[1*3+1]*A[2*3+2]
            - 2.0 * A[0*3+2]*de;                                     // c*d^2 + a*e^2 + b*f^2 - a*b*c - 2*f*d*e)

  double p, sqrt_p, q, c, s, phi;
  p = SQR(m) - 3.0*c1;
  q = m*(p - (3.0/2.0)*c1) - (27.0/2.0)*c0;
  sqrt_p = sqrt(fabs(p));

  phi = 27.0 * ( 0.25*SQR(c1)*(p - c1) + c0*(q + 27.0/4.0*c0));
  phi = (1.0/3.0) * atan2(sqrt(fabs(phi)), q);
  
  c = sqrt_p*cos(phi);
  s = (1.0/M_SQRT3)*sqrt_p*sin(phi);

  w[1]  = (1.0/3.0)*(m - c);
  w[2]  = w[1] + s;
  w[0]  = w[1] + c;
  w[1] -= s;

  return 0;
}
