#include "libecoPCR/ecoPCR.h"
#include <getopt.h>
#include <stdlib.h>

#define VERSION "0.1"

/* ----------------------------------------------- */
/* printout verbose mode                           */                                           
/* ----------------------------------------------- */
static void printTaxon(ecotx_t *taxon){
	printf("# taxid : %d | rank : %d | name : %s \n\n",taxon->taxid, taxon->rank, taxon->name);
}

/* ----------------------------------------------- */
/* printout help                                   */                                           
/* ----------------------------------------------- */
#define PP fprintf(stdout, 

static void PrintHelp()
{
        PP      "\n------------------------------------------\n");
        PP      " ecoisundertaxon Version %s\n", VERSION);
        PP      "------------------------------------------\n");
        PP      " synopsis : searching relationship in taxonomy\n");
        PP      " usage: ecoisundertaxon [options] database\n");
        PP      "------------------------------------------\n");
        PP      " options:\n");
        PP      " -1    : [FIRST] taxomic id of the hypothetical son\n\n");
        PP      " -2    : [SECOND] taxonomic id of the hypothetical parent\n\n");
        PP      " -h    : [H]elp - print <this> help\n\n");
        PP      " -v    : [V]erbose mode. Display taxonomic information for both\n");
        PP      "       : taxonomic id.\n\n");
        PP      "------------------------------------------\n");
        PP      " database : to match the expected format, the database\n");
        PP      " has to be formated first by the ecoPCRFormat.py program located.\n");
        PP      " in the tools directory. Type the radical only, leaving out the extension\n");
        PP      "------------------------------------------\n\n");
        PP		" https://www.grenoble.prabi.fr/trac/ecoPCR/wiki");
        PP      "------------------------------------------\n\n");        
}

#undef PP

/* ----------------------------------------------- */
/* printout usage and exit                         */
/* ----------------------------------------------- */

#define PP fprintf(stderr, 

static void ExitUsage(stat)
        int stat;
{
        PP      "usage: ecoisundertaxon [-1 taxid] [-2 taxid] [-v] [-h] datafile\n");
        PP      "type \"ecoisundertaxon -h\" for help\n");

        if (stat)
            exit(stat);
}

#undef  PP

/* ----------------------------------------------- */
/* MAIN						                       */
/* ----------------------------------------------- */

int main(int argc, char **argv){
	int32_t     	carg 		= 0;
	int32_t			taxid_1		= 0;
	int32_t			taxid_2		= 0;
	int32_t			verbose		= 0;
	int32_t			errflag		= 0;
	ecotaxonomy_t 	*taxonomy	= NULL;
	ecotx_t			*son		= NULL;
	ecotx_t			*parent		= NULL;
		
	
	while ((carg = getopt(argc, argv, "1:2:vh")) != -1) {
	    switch (carg) {                              
	        case '1':
	           sscanf(optarg,"%d",&taxid_1);
	           break;
	
	        case '2':     
	           sscanf(optarg,"%d",&taxid_2);
	           break;
	           
	        case 'v':     
	           verbose = 1;
	           break;
	           
	        case 'h':     
	           PrintHelp();
           	   exit(0);
	           break;
	           
	        case '?':
            	errflag++;
	     }
	}
	
	if ((argc -= optind) != 1)
    	errflag++;
    	
    if (errflag)
    	ExitUsage(errflag);
  		
	taxonomy = read_taxonomy(argv[optind],0);
	
	son = eco_findtaxonbytaxid(taxonomy, taxid_1);
	
	if (verbose){
		parent = eco_findtaxonbytaxid(taxonomy, taxid_2);
		printTaxon(son);
		printTaxon(parent);
	}
	
	if (eco_isundertaxon(son, taxid_2))
		printf("# taxid_1 (%d) is son of taxid_2 (%d)\n",taxid_1, taxid_2);
	else
		printf("# taxid_1 (%d) is NOT son of taxid_2 (%d)\n",taxid_1, taxid_2);
	
	return 0;
}