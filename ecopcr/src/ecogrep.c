#include "libecoPCR/ecoPCR.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <getopt.h>
#include <stdlib.h>
#include <sys/stat.h>


#define VERSION "0.1"

void getLineContent(char *stream, ecoseq_t *seq, ecoseq_t *oligoseq_1, ecoseq_t *oligoseq_2){
	
	int		i;
	char	*buffer;

 	for(	i=0, buffer = strtok(stream,"|");
			buffer != NULL;
			i++, buffer = strtok(NULL,"|"))
	{
			switch (i) {
				case 0:
					seq->AC = strdup(buffer);
					break;
				case 2:
					sscanf(buffer,"%d",&seq->taxid);
					break;
				case 13:
					oligoseq_1->SQ = strdup(buffer);
					oligoseq_1->SQ_length = strlen(buffer);
					break;
				case 16:
					oligoseq_2->SQ = strdup(buffer);
					oligoseq_2->SQ_length = strlen(buffer);
					break;
				case 20:
					seq->SQ = strdup(buffer);
					seq->SQ_length = strlen(buffer);
					break;
				default:
					break;
			}
	}
}


void freememory(char **tab, int32_t num){
	int32_t i;
	for (i=0;i<num-1;i++){
		ECOFREE(tab[i],"Error in freememory function");	
	}
}

/**
 * Check if pattern match a string using mamberall libapat function
 * @param	line		array containing sequence information
 * @param	pattern		array containing patterns to test on the sequence
 * @param	numpattern	number of pattern in pattern array
 * @param	error_max	error rate allowed by the user
 * 
 * @return	int			1 if a pattern match, else 0
 **/
int	ispatternmatching(ecoseq_t *seq, PatternPtr pattern){
	if (pattern != NULL)
	{
		SeqPtr 	apatseq	= NULL;
		apatseq=ecoseq2apatseq(seq,apatseq,0);
		return ManberAll(apatseq,pattern,0,0,apatseq->seqlen) > 0;
	}
	else return 0;
}

/* ----------------------------------------------- */
/* printout help                                   */                                           
/* ----------------------------------------------- */
#define PP fprintf(stdout, 

static void PrintHelp()
{
        PP      "\n------------------------------------------\n");
        PP      " ecogrep Version %s\n", VERSION);
        PP      "------------------------------------------\n");
        PP      " synopsis : filtering ecoPCR result based on\n");
        PP      " taxonomic id filter and regular expression pattern\n");
        PP      " usage: ecogrep [options] filename\n");
        PP      "------------------------------------------\n");
        PP      " options:\n");
        PP      " -1    : [FIRST] : compare the given pattern with direct strand oligonucleotide\n\n");
        PP      " -2    : [SECOND] : compare the given pattern with reverse strand oligonucleotide\n\n");
        PP      " -d    : [D]atabase containing taxonomic information\n\n");
        PP      " -e    : [E]rrors : max error allowed in pattern match (option-1, -2 and -p) (0 by default)\n\n"); 
        PP      " -p    : [P]attern : oligonucleotide pattern\n\n"); 
        PP      " -h    : [H]elp : print <this> help\n\n");
        PP      " -i    : [I]gnore subtree under given taxonomic id\n\n");
        PP      " -r    : [R]estrict search to subtree under given taxomic id\n\n");
        PP      " -v    : in[V]ert the sense of matching, to select non-matching lines.\n\n");
        PP      " argument:\n");
        PP      " ecoPCR ouput file name\n");
        PP      "------------------------------------------\n\n");
        PP		" http://www.grenoble.prabi.fr/trac/ecoPCR/\n");
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
        PP      "usage: ecogrep [-d database] [-p pattern] [-i taxid] [-r taxid] [-v] [-h] <file name>\n");
        PP      "type \"ecogrep -h\" for help\n");

        if (stat)
            exit(stat);
}

#undef  PP

/* ----------------------------------------------- */
/* MAIN						                       */
/* ----------------------------------------------- */

#define LINE_BUFF_SIZE 10000

int main(int argc, char **argv){
	int32_t     	carg 			= 0;
	int32_t			r				= 0;			// number of restricted taxid
	int32_t			i				= 0;			// number of ignored taxid
	int32_t			v				= 0;			// stores if -v mode is active
	int32_t			k				= 0;			// file counter
	int32_t			errflag			= 0;			
	int32_t       	error_max		= 0;			// stores the error rate allowed by the user
	int32_t			matchingresult 	= 0;			// stores number of matching result
	
	ecotaxonomy_t	*taxonomy;						// stores the taxonomy
	
	ecoseq_t		*seq				= NULL;		// stores sequence info
	ecoseq_t		*oligoseq_1			= NULL;		// stores the oligo_1 info
	ecoseq_t		*oligoseq_2			= NULL;		// stores the oligo_2 info
	
	char			*database			= NULL;		// stores the database path (for taxonomy)
	
	char			*p					= NULL;		// pattern for sequence 
	PatternPtr		pattern				= NULL;		// stores the build pattern for sequence 
	char			*o1					= NULL;		// pattern for direct strand oligo
	PatternPtr		oligo_1				= NULL;		// stores the build pattern	for direct strand oligo
	char			*o2					= NULL;		// pattern for reverse strand oligo
	PatternPtr		oligo_2				= NULL;		// stores the build pattern for reverse strand oligo
		
	int32_t		  	*restricted_taxid 	= NULL;		// stores the restricted taxid
	int32_t       	*ignored_taxid		= NULL;		// stores the ignored taxid

	FILE			*file				= NULL;		// stores the data stream, stdin by default
	char			*stream				= ECOMALLOC(sizeof(char *)*LINE_BUFF_SIZE,"error stream buffer allocation"); 
	char			*orig				= ECOMALLOC(sizeof(char *)*LINE_BUFF_SIZE,"error orig buffer allocation"); 

	int 			is_ignored 			= 0;		
	int				is_included			= 0;
	int				is_matching			= 0;
	int 			match_o1			= 0;
	int 			match_o2			= 0;
	int				good				= 0;
	
	seq = new_ecoseq();
	oligoseq_1 = new_ecoseq();
	oligoseq_2 = new_ecoseq();
	
	/**
	 * Parse commande line options
	 **/
	while ((carg = getopt(argc, argv, "1:2:p:d:i:r:e:vh")) != -1) {

	    switch (carg) {                              
	        case '1':
	           o1 = ECOMALLOC(strlen(optarg)+1,
	           						"Error on o1 allocation");
	           strcpy(o1,optarg);
	           break;

	        case '2':
	           o2 = ECOMALLOC(strlen(optarg)+1,
	           						"Error on o2 allocation");
	           strcpy(o2,optarg);
	           break;
	           	        
	        case 'd':
	           database = ECOMALLOC(strlen(optarg)+1,
	           						"Error on datafile allocation");
	           strcpy(database,optarg);
	           break;
	
	        case 'i':   
	           ignored_taxid = ECOREALLOC(	ignored_taxid,
	           								sizeof(int32_t)*(i+1),
	           								"Error on ignored_taxid reallocation");
	           sscanf(optarg,"%d",&ignored_taxid[i]);
	           i++;
	           break;
	        
	        case 'r': 
	           restricted_taxid = ECOREALLOC(	restricted_taxid, 
	           									sizeof(int32_t)*(r+1),
	           									"Error on restricted_taxid reallocation");  
	           sscanf(optarg,"%d",&restricted_taxid[r]);
	           r++;
	           break;
	           	           
	        case 'v':     
	           v = 1;
	           break;
	           
	        case 'h':     
	           PrintHelp();
           	   exit(0);
	           break;
	           
	        case 'e':	
	           sscanf(optarg,"%d",&error_max);
    	       break;

	        case 'p':	
	           p = ECOMALLOC(strlen(optarg)+1,
	           						"Error on pattern allocation");		   
	           strcpy(p,optarg);
	           break;
	           
	        case '?':
            	errflag++;
	     }
	}

	/**
	 * Check sequence pattern length and build it in PatternPtr format
	 **/
	if(p)
	{
		if (strlen(p) > 32){
			printf("# Sorry, ecogrep doesn't handle pattern longer than 32 characters.\
					\n# Please check it out : %s\n",p);
			exit(EXIT_FAILURE);
		}
		else if ( (pattern = buildPattern(p,error_max)) == NULL)
			exit(EXIT_FAILURE);
	}
	
			
			
	/**
	 * Check o1 pattern length and build it in PatternPtr format
	 **/
	if(o1)
	{
		if (strlen(o1) > 32){
			printf("# Sorry, ecogrep doesn't handle pattern longer than 32 characters.\
					\n# Please check it out : %s\n",o1);
			exit(EXIT_FAILURE);
		}
		else if ( (oligo_1 = buildPattern(o1,error_max)) == NULL)
			exit(EXIT_FAILURE);
	}

	/**
	 * Check o2 pattern length and build it in PatternPtr format
	 **/
	if(o2)
	{
		if (strlen(o2) > 32){
			printf("# Sorry, ecogrep doesn't handle pattern longer than 32 characters.\
					\n# Please check it out : %s\n",o2);
			exit(EXIT_FAILURE);
		}
		else if ( (oligo_2 = buildPattern(o2,error_max)) == NULL)
			exit(EXIT_FAILURE);
	}

	/**
	 * try to get the database name from environment variable
	 * if no database name specified in the -d option
	 **/
	if (database == NULL)
	{
		database = getenv("ECOPCRDB");
		if (database == NULL)
			errflag++;
	}
	
	/**
	 * check at leat one processing is asked
	 * either patterns or taxid filters
	 **/			
	if ( !p && !o1 && !o2 && restricted_taxid == NULL && ignored_taxid == NULL )
	{
		errflag++;
	}		
	if (errflag)
    	ExitUsage(errflag);

	/**
	 * Get the taxonomy back
	 **/
	taxonomy = read_taxonomy(database,0);
	 
	/**
	 * Parse the stream
	 **/
	for (k=0 ; argc >= optind ; optind++, k++){

		matchingresult = 0;

		if ( (file = fopen(argv[optind], "r")) == NULL)
		{
			if (isatty(fileno(stdin)) == 0)
			{
		 		file = stdin;
		 		printf("# Processing standard input...\n");
			}
	 		else
	 			break;	
		}
		else
			printf("# Processing %s...\n",argv[optind]);
		
		while( fgets(stream, LINE_BUFF_SIZE, file) != NULL ){
		 		
		 	if (stream[0]!= '#')
		 	{
		 		
		 		stream[LINE_BUFF_SIZE-1]=0;
		 		
		 		strcpy(orig,stream);
		 		
				getLineContent(stream,seq,oligoseq_1,oligoseq_2);	

				/* -----------------------------------------------*/
				/* is ignored if at least one option -i 		  */
				/*  AND											  */
				/* if current sequence is son of taxid 		 	  */
				/* -----------------------------------------------*/ 	
			    is_ignored = ( (i > 0) && (eco_is_taxid_included(	taxonomy, 
			                            							ignored_taxid, 
			                            							i, 
			                            							seq->taxid))
			                 );				

				/* -----------------------------------------------*/
				/* is included if no -r option					  */
				/*  OR											  */
				/*  if current sequence is son of taxid			  */
				/* -----------------------------------------------*/ 				
				is_included = ( (r == 0) || (eco_is_taxid_included(	taxonomy, 
			                            									restricted_taxid, 
			                            	 								r, 
			                            									seq->taxid))
			                  );

				/* ----------------------------------------------------------- */
				/* match if no pattern or if pattern match current sequence    */
				/* ----------------------------------------------------------- */ 
			    is_matching = ( !p || (ispatternmatching(seq,pattern)));
			    
			    /* ---------------------------------------------------------------------------- */
				/* match if no direct oligo pattern or if pattern match current direct oligo    */
				/* ---------------------------------------------------------------------------- */ 
			    match_o1 = (!o1 || (ispatternmatching(oligoseq_1,oligo_1)));
			    
			    /* ------------------------------------------------------------------------------- */
				/* match if no revesrse oligo pattern or if pattern match current reverse oligo    */
				/* ------------------------------------------------------------------------------- */ 
			    match_o2 = (!o2 || (ispatternmatching(oligoseq_2,oligo_2)));
			        				
			    good = (is_included  && 	is_matching && 	!is_ignored && match_o1 && match_o2);
			    
			    if (v)
			       good=!good;
			       			
				if 	( good )
				{
					printf("%s",orig);
	     			matchingresult++;
				}
		 	}
		}
		if ( file != stdin )
    		fclose(file);
    		
    	printf("# %d matching result(s)\n#\n",matchingresult);
	}
	 	
	/**
	 * clean and free before leaving
	 **/
	ECOFREE(orig,"Error in free orig");
	ECOFREE(stream,"Error in free stream");
	ECOFREE(ignored_taxid,"Error in free stream");
	ECOFREE(restricted_taxid,"Error in free stream");
	
	return 0;
}
