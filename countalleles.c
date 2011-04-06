
/* [[file:~/Work/simsec/org/simsec.org::*Counting%20alleles][block-28]] */
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <unistd.h>
#include <stdarg.h>

#define NUMSTATES 5
enum {NA=-1, A=0, C=1, G=2, T=3, N=4} ;

void usage() ;
void error(char *fmt, ...) ;
int encode(char *seq, int seqlen, int *iseq) ;

int main() {
    int linelen, i, s ;
    char *seq = NULL ;
    size_t maxlinelength=0 ;
    size_t maxiseqlength=0 ;
	 int counts[NUMSTATES] ;
    int *iseq = NULL ;

    while( (linelen = getline(&seq, &maxlinelength, stdin)) > 0 ) {
        if(maxlinelength > maxiseqlength) {
			   maxiseqlength = maxlinelength ;
				iseq = realloc(iseq, maxiseqlength * sizeof(int));
        }
        linelen-- ;
        encode(seq, linelen, iseq) ;
        for(s = 0 ; s < NUMSTATES ; s++) counts[s] = 0 ;
        for(i = 0 ; i < linelen ; i++)
            counts[iseq[i]]++ ;
        for(s = 0 ; s < NUMSTATES ; s++)
            printf("%d%s", counts[s], s == NUMSTATES-1 ? "" : "\t") ;
        printf("\n") ;
    }    
}

int encode(char *seq, int seqlen, int *iseq) {
    int i, nmiss=0 ;
    char c ;
    
    for(i = 0 ; i < seqlen ; i++) {
        c = seq[i] ;
        if(c == 'N') ++nmiss ;
        iseq[i] = 
            c == 'A' ? 0 :
            c == 'C' ? 1 :
            c == 'G' ? 2 :
            c == 'T' ? 3 :
            c == 'N' ? 4 :
            NA ;
        if(iseq[i] == NA)
            error("Invalid base: %c\n", c) ;
    }
    return nmiss ;
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
    error("aslink -c numchars -d maxdiff [-n maxmiss]\n") ;
}
/* block-28 ends here */
