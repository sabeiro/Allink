#include "fdmdpd2.h"

/* functions like these are typically called from forces_calc() in main.c 
 * they have to be called after set_xf_intra() and before add_forces()*/
struct blockforce
{
        int b;
        VEC f;
};

struct multiblockforce
{
        int start_b; // block number from which on the blocks are subjected to a force
        int n; // number of consecutive blocks which are subject to the force - THESE SHOULD BE AT THE END OF THE CONFIGURATION FILE!!!
	// MF TODO BYSY HERE: allocate these vectors in read_header explicitely, using 3*n as length
        VEC *min; // position of potential minimum for each block
        double* coms; // calculated center of mass for each block
	double k; // force constant
};

void blockforce_print_header(FILE *FH, struct blockforce *blockforce)
{
        if (blockforce == NULL) return;

        int i;

        fprintf(FH, "# blockforce b= %d ",blockforce->b);
        fprintf(FH, "f= ");
        for (i=0 ; i<3 ; i++) {
                fprintf(FH, "%lg ",blockforce->f[i]);
        }
        fprintf(FH, "\n");
}

void multiblockforce_print_header(FILE *FH, struct multiblockforce *multiblockforce)
{
        if (multiblockforce == NULL) return;

        int i,j;

        fprintf(FH, "# multiblockforce start_b= %d ",multiblockforce->start_b);
        fprintf(FH, "n= %d ",multiblockforce->n);
        fprintf(FH, "minima= ");
        for (i=0 ; i<multiblockforce->n ; i++) {
	        for (j=0 ; j<3 ; j++) {
        	        fprintf(FH, "%lg ",multiblockforce->min[i][j]);
        	}
	}
        fprintf(FH, "k= %lg ",multiblockforce->k);
        fprintf(FH, "\n");
}

void blockforce_read_header(FILE *FH, struct beads *b)
{
        int i;

        b->blockforce = malloc(sizeof(*b->blockforce));
        if (b->blockforce == NULL) novm("blockforce");


	fscanf(FH, "# blockforce b= ");
	if (fscanf(FH, "%d ", &b->blockforce->b) != 1)
		fatal(EINVAL, "blockforce b!");
	fscanf(FH, "f= ");
	for (i = 0; i < ARRAY_SIZE(b->blockforce->f); i++) {
		if (fscanf(FH, "%lg ", &b->blockforce->f[i]) != 1)
			fatal(EINVAL, "blockforce f");
	}

        printf("Using blockforce! f= ");
	for (i = 0; i < ARRAY_SIZE(b->blockforce->f); i++) 
                printf("%lg ",b->blockforce->f[i]);
        printf("acting on block %d.\n", b->blockforce->b);
}

void multiblockforce_read_header(FILE *FH, struct beads *b)
{
        int i,j;

        b->multiblockforce = malloc(sizeof(*b->multiblockforce));
        if (b->multiblockforce == NULL) novm("multiblockforce");


	fscanf(FH, "# multiblockforce start_b= ");
	if (fscanf(FH, "%d ", &b->multiblockforce->start_b) != 1)
		fatal(EINVAL, "multiblockforce start_b!");
	fscanf(FH, "n= ");
	if (fscanf(FH, "%d ", &b->multiblockforce->n) != 1)
		fatal(EINVAL, "multiblockforce n!");

        b->multiblockforce->coms = calloc(b->multiblockforce->n*3,sizeof(*b->multiblockforce->coms));
        if (b->multiblockforce->coms == NULL) novm("multiblockforce->coms");
        b->multiblockforce->min = calloc(b->multiblockforce->n,sizeof(*b->multiblockforce->min));
        if (b->multiblockforce->coms == NULL) novm("multiblockforce->min");

	fscanf(FH, "minima= ");
	for (i = 0; i < b->multiblockforce->n; i++) {
		for (j = 0; j < 3 ; j++)
			if (fscanf(FH, "%lg ", &b->multiblockforce->min[i][j]) != 1)
				fatal(EINVAL, "multiblockforce min");
	}
	fscanf(FH, "k= ");
        if (fscanf(FH, "%lg ", &b->multiblockforce->k) != 1)
			fatal(EINVAL, "multiblockforce k!");


        printf("Using multiblockforce! k=%lg starting from block %d for %d blocks with the following minima:\n",b->multiblockforce->k,b->multiblockforce->start_b,b->multiblockforce->n);
	for (i = 0; i < ARRAY_SIZE(b->multiblockforce->min); i++) {
		for (j = 0; j < 3 ; j++) {
                	printf("%lg ",b->multiblockforce->min[i][j]);
		}
	}
	printf("\n");
}

double apply_blockforce(struct beads *b, struct blockforce *blockforce)
{
        if (blockforce == NULL) return;
        int i, j, k;
        int idx = 0;

	for (i = 0; i < b->local_n; i++) {
                int thisblock = b->local_n_b[i];

		int N = b->N[thisblock];


		if (thisblock == blockforce->b) {
			for (j = idx; j < idx + N; j++) {
				for (k = 0; k < 3; k++) {
					b->f_intra[j][k] += blockforce->f[k];
				}
			}
		}
		idx += N;
	}
        return 0.;
}

void multiblockforce_get_coms(struct beads *b, struct multiblockforce *multiblockforce) //adapted from usp_get_coms
{
        int i, j, k;
        int idx = 0;

        /* b->local_n is the number of molecules which belong to this
         * MPI process */ 
        //const int dim[3] = {TANG1,TANG2,NORMAL};
        for (i=0 ; i<multiblockforce->n ; i++)
        	for (j=0 ; j<3 ; j++)
                	multiblockforce->coms[i*3+j]=0.;
	for (i = 0; i < b->local_n; i++) { // GET COMS
                //FOR EACH ATOM CHECK IF THE BLOCK IS SUBJECT TO THE UMBRELLA POTENTIAL
                int thisblock = b->local_n_b[i];
                /* b->N[] is an array with b->blocks entries, specifying
                 * the number of beads per molecule in this block */
		int N = b->N[thisblock];
                if (thisblock >= multiblockforce->start_b) {



                        /* the total number of molecules in this block in all
                         * processes can be found in b->n[thisblock]. */

		        for (j = idx; j < idx + N; j++) {
		        	for (k = 0; k < 3; k++) {
                                        /* the unfolded, i.e., absolute,
                                         * coordinates of this bead are located
                                         * in b->x_intra[][]. */

                                        //THE COORDINATES ARE ADDED INTO THE COMS
                                        //usp->coms[3*blck+k] += b->x_intra[j][dim[k]];
                                        multiblockforce->coms[(thisblock-multiblockforce->start_b)*3+k] += b->x_intra[j][k];
                                        //fprintf(stderr,"com[%d]=%lg, x[i]=%lg, dim[i]=%d\n",k,usp->coms[3*blck+k],b->x_intra[j][dim[k]],dim[k]);
		        	}
		        }
                }
		idx += N;
	}
	MPI_Allreduce(MPI_IN_PLACE, multiblockforce->coms, multiblockforce->n*3, MPI_DOUBLE, MPI_SUM, comm_grid);
        for (i=0 ; i<multiblockforce->n ; i++) //blocks //NORMALIZE COMS
                for (j=0 ; j<3 ; j++) {//dims
                        multiblockforce->coms[i*3+j] /= b->n[multiblockforce->start_b+i] * b->N[multiblockforce->start_b+i];
                }
}

void multiblockforce_apply_forces(struct beads *b, struct multiblockforce *multiblockforce) //adapted from usp_apply_forces
{
        int i, j, k, block, dim;
        int idx = 0;
        double f[multiblockforce->n*3];
        double comsMI3d[multiblockforce->n*3]; //center of mass with minimum image correction

        //const int dim[3] = {TANG1,TANG2,NORMAL};

        //CALCULATE FORCES
        //calculate distance. for this, the minimum image convention has to be applied.
	for (block = 0 ; block < multiblockforce->n ; block++) {
		for (dim = 0 ; dim < 3 ; dim++) {
                        comsMI3d[3*block+dim] = multiblockforce->coms[3*block+dim];
			fprintf(stderr,"COM of block %d dim %d before MIC: %lg\n",block,dim,multiblockforce->coms[3*block+dim]);
                        k = (comsMI3d[3*block+dim] - multiblockforce->min[block][dim] > 0) ? -1 : 1;
                        //while ( fabs(comsMI3d[2*j+i] + k*b->l[dim[i]] - usp->x0s[2*j+i]) < fabs(comsMI3d[2*j+i] - usp->x0s[2*j+i]) ) {
                        //        comsMI3d[2*j+i] += k*b->l[dim[i]];
			//}
                        while ( fabs(comsMI3d[3*block+dim] + k*b->l[dim] - multiblockforce->min[block][dim]) < fabs(comsMI3d[3*block+dim] - multiblockforce->min[block][dim]) ) {
                                comsMI3d[3*block+dim] += k*b->l[dim];
                        }
			fprintf(stderr,"COM of block %d dim %d after MIC: %lg\n",block,dim,multiblockforce->coms[3*block+dim]);
		}
	}
        //calculate force
	for (block = 0 ; block < multiblockforce->n ; block++) {
		for (dim = 0 ; dim < 3 ; dim++) {
                        f[3*block+dim] = -multiblockforce->k * (comsMI3d[3*block+dim] - multiblockforce->min[block][dim]);
			fprintf(stderr,"calculated force is %d %d: %lg\n",block,dim,f[3*block+dim]); //MF DEBUG
                }
        }


        //APPLY FORCES.
	for (i = 0; i < b->local_n; i++) {
                //FOR EACH ATOM CHECK IF THE BLOCK IS SUBJECT TO THE UMBRELLA POTENTIAL
                int thisblock = b->local_n_b[i];
		int N = b->N[thisblock];
                if (thisblock >= multiblockforce->start_b) {
		        for (j = idx; j < idx + N; j++) {
		        	for (k = 0; k < 3; k++) { //dimension
                                        //apply forces
					b->f_intra[j][k] += f[3*(thisblock - multiblockforce->start_b)+k];
		        	}
		        }
                }
		idx += N;
	}
}

void multiblockforce_wrapper(struct beads *b, struct multiblockforce *multiblockforce)
{
        if (multiblockforce == NULL) return;
        multiblockforce_get_coms(b,b->multiblockforce);
        multiblockforce_apply_forces(b,b->multiblockforce);
}

void blockforce_free(struct blockforce *blockforce)
{
        if (blockforce == NULL) return;
        free(blockforce);
}

void multiblockforce_free(struct multiblockforce *multiblockforce)
{
        if (multiblockforce == NULL) return;
        free(multiblockforce->min);
        free(multiblockforce->coms);
        free(multiblockforce);
}
