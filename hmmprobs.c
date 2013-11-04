
/* [[file:~/Work/simsec/org/simsec.org::*C][block-35]] */
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <unistd.h>
#include <stdarg.h>
#include <math.h>

#define J 4 /* Number of alleles */
#define NUMPOPS 2
#define K 3 /* Number of ancestry states */

void usage() ;
void error(char *fmt, ...) ;

int main(int argc, char **argv) {
    int T = 0, t, j, k, c, a, a1, a2, pop, n, Q=50, y[Q], ploidy = -1, nstates ;
    double *popfreq, e[Q],//eps, //eps=-1,
        p[NUMPOPS][J],
        Pr_y_given_z[K],
        Pr_a_given_g[J][J][J],
        Pr_g_given_z[K], Pr_y_given_g ;

    while((c = getopt(argc, argv, "p:t:n:")) != -1) {
        switch(c) {
        case 'p': ploidy = atoi(optarg) ; break ;
        case 't': T = atoi(optarg) ; break ;
		/* read in number of reads per site */
		case 'n': Q = atoi(optarg) ; break ;
        case '?': usage() ;
        }
    }

    if(!(T > 0 && (ploidy == 1 || ploidy == 2))) usage() ;
    nstates = (ploidy == 2 ? K : 2) ;

    for(a1 = 0 ; a1 < J ; a1++) {
        for(a2 = 0 ; a2 < J ; a2++) {
            for(a = 0 ; a < J ; a++) {
                n = (a == a1) + (ploidy == 2 && a == a2) ;
                Pr_a_given_g[a][a1][a2] = n/(double)ploidy;
                //Pr_a_given_g[a][a1][a2] = (n/(double)ploidy) * (1 - 3*eps/4) + (1 - n/(double)ploidy) * eps/4 ;
            }
        }
    }

    for(t = 0 ; t < T ; t++) {
        for(j = 0 ; j < Q ; j++)
            if(scanf("%d", y+j) != 1) error("(hmmprob) Error reading allele %d on line %d\n", j, t+1) ;
        for(pop = 0 ; pop < NUMPOPS ; pop++) {
            for(j = 0 ; j < J ; j++)
                if(scanf("%lf", p[pop] + j) != 1)
                    error("(hmmprob) Error reading probability of allele %d in population %d on line %d\n", j, pop, t+1) ;
        }

        //if(scanf("%lf", &eps) != 1) error("(hmmprob) Error reading eps %f on line %d\n", eps, t+1) ;
		/* Read in eps for each read */
        for(j = 0 ; j < Q ; j++)
            if(scanf("%lf", e+j) != 1) error("(hmmprob) Error reading eps of allele %d on line %d\n", j, t+1) ;
			
        for(k = 0 ; k < nstates ; k++) Pr_y_given_z[k] = 0 ;

        for(a1 = 0 ; a1 < J ; a1++) {
            for(a2 = 0 ; a2 < J ; a2++) {
                if( ploidy == 1 && a2 > 0 ) break ;
                /* Treat each read as an independent event and multiply Pr_y_given_g  */
                Pr_y_given_g = 1 ;
                for(a = 0 ; a < Q ; a++) {
                   // Pr_y_given_g *= pow(Pr_a_given_g[a][a1][a2] * (1 - 3*eps/4) + (1 - Pr_a_given_g[a][a1][a2]) * eps/4, y[a]) ;
					if (y[a]!=5)
						Pr_y_given_g *= Pr_a_given_g[y[a]][a1][a2] * (1 - 3*e[a]/4) + (1 - Pr_a_given_g[y[a]][a1][a2]) * e[a]/4 ;
				}
                /* Pr this genotype given ancestry */
                if(ploidy == 2) {
                    Pr_g_given_z[0] = p[0][a1] * p[0][a2] ;
                    Pr_g_given_z[1] = p[0][a1] * p[1][a2] ;
                    Pr_g_given_z[2] = p[1][a1] * p[1][a2] ;
                }
                else {
                    Pr_g_given_z[0] = p[0][a1] ;
                    Pr_g_given_z[1] = p[1][a1] ;
                }
                /* Increment average over genotypes */
                for(k = 0 ; k < nstates ; k++)
                    Pr_y_given_z[k] += Pr_y_given_g * Pr_g_given_z[k] ;
            }
        }
        for(k = 0 ; k < nstates ; k++)
            printf("%e%s", Pr_y_given_z[k], k < (nstates-1) ? "\t" : "\n") ;
    }
    return 0 ;
}
void error(char *fmt, ...) {
    va_list args;

    fflush(stderr);
    
    va_start(args, fmt);
    vfprintf(stderr, fmt, args);
    va_end(args);
    
    fflush(stderr) ;
    exit(2) ;
}

void usage() {
    error("hmm -t numobs -p ploidy < inputfile\n") ;
    //error("hmm -t numobs -e epsilon < inputfile\n") ;
}
/* block-35 ends here */
