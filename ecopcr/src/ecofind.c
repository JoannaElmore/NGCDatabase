#include "libecoPCR/ecoPCR.h"
#include <regex.h>
#include <string.h>
#include <stdlib.h>
#include <getopt.h>
#include <stdio.h>
#define VERSION "0.1"

/**
 * display the result
 **/

void displayPath(ecotx_t *taxon, ecotaxonomy_t *taxonomy){

	if (taxon != taxon->parent){
		displayPath(taxon->parent,taxonomy);
		printf(";");
	}
	if (rank_index("no rank",taxonomy->ranks) != taxon->rank)
		printf("%s:", taxonomy->ranks->label[taxon->rank]);
	printf("%s", taxon->name);
}


static void printresult(ecotx_t *taxon,econame_t* name,ecotaxonomy_t *taxonomy, int32_t pathDisplay)
{
	char* rankname;
	char* classname;
	char* matchedname=taxon->name;
	
	classname="scientific name";
	if (name)
	{
		classname=name->classname;
		matchedname=name->name;
	}
	
	rankname= taxonomy->ranks->label[taxon->rank];	
	
	printf("%10d \t| %15s \t|\t %-50s \t|\t %15s \t|\t %s",
			taxon->taxid,
			rankname,
			matchedname,
			classname,
			taxon->name);
	if (pathDisplay) {
	    printf("\t|\t");
	    displayPath(taxon, taxonomy);
	}
    printf("\n");
}

/**
 * display header before printing any result
 **/
static void printheader(int32_t pathDisplay)
{
	printf("# %12s \t| %15s \t|\t %-50s \t|\t %-15s \t|\t %s%s\n#\n",
			"taxonomy id",
			"taxonomy rank",
			"name",
			"class name",
			"scientific name",
			pathDisplay ? "\t|\t path":"");
}


/**
 * display son's list for given taxon
 **/
static void get_son(ecotaxonomy_t *taxonomy, ecotx_t *taxon, int32_t *count, char *rankname, int32_t pathDisplay)
{
	int32_t		i;
	ecotx_t 	*current_taxon;
	
	for (	i = 0, current_taxon = taxonomy->taxons->taxon;
			i < taxonomy->taxons->count; 
			i++, current_taxon++)
	{	

		if (taxon != current_taxon && taxon->taxid == current_taxon->parent->taxid)
		{
			if (rankname == NULL || !strcmp(rankname,taxonomy->ranks->label[current_taxon->rank]))
			{
				printresult(current_taxon, NULL, taxonomy, pathDisplay);
				(*count)++;	
			}
			get_son(taxonomy,current_taxon,count,rankname, pathDisplay);
		}
	}
}



/**
 * display list of rank filter option (-l option)
 **/
static void listfilteroptions(ecorankidx_t *ranks)
{
	int32_t		i;
	
	printf("#\n");
	
	for (	i=0;
			i < ranks->count;
			i++)
	{
		printf("#  %s \n",ranks->label[i]);
	}
	
	printf("#\n");
}


/* ---------------------------------------- */
/* get back on given taxid taxonomic parent	*/
/* and display it							*/
/* ---------------------------------------- */
void gettaxidparents(int32_t taxid, ecotaxonomy_t *taxonomy, char *rankname, int32_t pathDisplay)
{
	ecotx_t *next_parent;
	int32_t	c = 0;
	
	next_parent = eco_findtaxonbytaxid(taxonomy, taxid);		
	
	printheader(pathDisplay);
	
	printresult(next_parent, NULL,taxonomy, pathDisplay);
		
	while ( strcmp(next_parent->name, "root") )
	{
		next_parent = next_parent->parent;
		if (rankname == NULL || !strcmp(rankname,taxonomy->ranks->label[next_parent->rank]))
		{
			printresult(next_parent, NULL,taxonomy, pathDisplay);
			c++;
		}
	}
	
	printf("#  %d parent(s) found\n#\n",c);
}


/**
 * printout usage and exit                     
 **/
#define PP fprintf(stderr, 

static void ExitUsage(stat)
        int stat;
{
        PP      "usage: ecofind [-d database] [-h] [-l] [-P] [-r taxonomic rank] [-p taxid] [-s taxid] <taxon name pattern> ... \n");
        PP      "type \"ecofind -h\" for help\n");
        if (stat)
            exit(stat);
}

#undef  PP

/**
 * printout help
 **/
#define PP fprintf(stdout, 

static void PrintHelp()
{
        PP      "------------------------------------------\n");
        PP      " ecofind Version %s\n", VERSION);
        PP      "------------------------------------------\n");
        PP      "synopsis : searching for taxonomic and rank and\n");
        PP		"           taxonomy id for given regular  expression patterns\n\n");
        PP      "usage: ecofind [options] <patterns>\n");
        PP      "------------------------------------------\n");
        PP      "options:\n");
        PP      "-a : [A]ll enable the search on all alternative names and not only scientific names.\n\n");        
        PP      "-d : [D]atabase containing the taxonomy.\n");
        PP		"     To match the expected format, the database\n");
        PP		"     has to be formated first by the ecoPCRFormat.py\n");
        PP		"     program located in the tools directory.\n");
        PP		"     Write the database radical without any extension.\n\n");
        PP      "-h : [H]elp - print <this> help\n\n");
        PP      "-l : [L]ist all taxonomic rank available for -r option\n\n");
        PP      "-P : [P]ath : add a column containing the full path for each displayed taxon\n\n");
        PP      "-p : [P]arents : specifiying this option displays all parental tree's information for the given taxid.\n\n");
        PP      "-r : [R]estrict to given taxonomic rank\n\n");      
        PP      "-s : [S]ons: specifiying this option displays all subtree's information for the given taxid.\n\n");
        PP      "-P : Display taxonomic [P]ath as suplementary column in output\n\n");
        PP      "arguments:\n");        
        PP      "<taxon> name pattern bearing regular expressions\n\n");
        PP      "------------------------------------------\n");
        PP		" http://www.grenoble.prabi.fr/trac/ecoPCR/\n");
        PP      "------------------------------------------\n\n");
}

/* ----------------------------------------------- */

#define PATTERN_NUMBER 	10
#define PATTERN_LENGHT 	40
#define	RESULT_LENGTH	100

int main(int argc, char **argv) 
{
	int32_t     	carg		= 0;
	int32_t			nummatch 	= 0;
	int32_t 		k,j 		= 0;
	int32_t     	errflag 	= 0;
	int32_t			tax_count	= 0;
	int32_t         alternative = 0;
	char			*prefix		= NULL;
	ecotaxonomy_t 	*taxonomy;	
	econame_t		*name;
	int32_t         name_count;
	
	int				re_error;
	int				re_match;
	regex_t			re_preg;
	
	int32_t 		uptree 		= 0;
	int32_t 		subtree 	= 0;
	char			*rankname	= NULL;
	int32_t			rankfilter 	= 1;
	int32_t			list 		= 0;
	int32_t			path 		= 0;
	
	ecotx_t 		*subtree_parent;
	int32_t			count_son = 0;


	while ((carg = getopt(argc, argv, "had:p:s:r:lP")) != -1) {
	    switch (carg) {                        
	       case 's':               	/* path to the database     */
	          sscanf(optarg,"%d",&subtree);
	          break;

	       case 'r':               	/* rank filter     */
	          rankname = ECOMALLOC(strlen(optarg)+1,"allocation rankname");
	          strcpy(rankname,optarg);
	          rankfilter = 0;
	          break;

	       case 'd':               	/* path to the database     */
	          prefix = ECOMALLOC(strlen(optarg)+1,"allocation prefix");
	          strcpy(prefix,optarg);
	          break;
	          
	       case 'l':               	/* list rank filter options */
	          list = 1;
	          break;

	       case 'P':               	/* Path output option */
	          path=1;
	          break;
	
	       case 'a':               	/* allow alternative names  */
	          alternative = 1;
	          break;
	
	       case 'h':               	/* display help	           	*/      
	          PrintHelp();
	          exit(0);
	          break;
	          
	       case 'p':               	/* taxid		           	*/      
	          sscanf(optarg,"%d",&uptree);
	          break;	          
	           
	       case '?':               	/* bad option   	        */
	          errflag++;
	     }
	}
	
	if ((argc - optind) < 1)
        errflag++;
	
	if (prefix == NULL)
	{
		prefix = getenv("ECOPCRDB");
		if (prefix == NULL)
			errflag++;
	}
		
	if (errflag && !uptree && !rankname && !subtree && !list)
		ExitUsage(errflag);
		
	/**
	 * load taxonomy using libecoPCR functions
	 **/
	printf("# \n#  opening %s database\n",prefix);
	
	taxonomy = read_taxonomy(prefix,1);
	tax_count = taxonomy->taxons->count;
	name_count = taxonomy->names->count;
		
	
	/* ----------------------------------------	*/
	/* list -r option possibility		    	*/
	/* ----------------------------------------	*/		
	if (list)
	{
		listfilteroptions(taxonomy->ranks);
		return 0;
	}
		
	/* ----------------------------------------	*/
	/* display taxid parent if -t option    	*/
	/* specified in command line				*/
	/* ----------------------------------------	*/	
	if (uptree)
	{
		gettaxidparents(uptree,taxonomy,rankname, path);
		return 0;	
	}
	
	/* ----------------------------------------	*/
	/* display taxid sons if -s option    		*/
	/* specified in command line				*/
	/* ----------------------------------------	*/		
	if (subtree)
	{
		printheader(path);
		subtree_parent = eco_findtaxonbytaxid(taxonomy,subtree);
		printresult(subtree_parent, NULL,taxonomy, path);
		get_son(taxonomy, subtree_parent,&count_son,rankname, path);
		printf("#  %d son(s) found\n#\n",count_son);
		return 0;
	}
	
	printf("#  %d taxons\n", tax_count);
		
	/**
	 * parse taxonomy
	 **/
	for (k=optind;k<argc;k++)
	{
		printf("#\n#  searching for '%s' pattern\n",argv[k]); 
		
		re_error = regcomp (&re_preg, argv[k], REG_NOSUB | REG_EXTENDED | REG_ICASE);
		if (re_error)
		{
			fprintf(stderr,"#  misformed pattern '%s'\n",argv[k]);
			exit(1);
		}
		
		nummatch=0;
		
	    printheader(path);
		
		for (j=0,name=taxonomy->names->names;
		     j < name_count;
		     name++,j++)
		{  			
   			
   			if(rankname)
   				rankfilter = !(strcmp(rankname,taxonomy->ranks->label[name->taxon->rank]));

   			re_match = regexec (&re_preg, name->name, 0, NULL, 0);
   			
   			if (!re_match && (alternative || name->is_scientificname) && rankfilter)
   			{
   				printresult(name->taxon,name,taxonomy, path);
   				nummatch++;
   			}
   			
		} 		
			
		printf("#  %d records found \n",nummatch);
	   	regfree(&re_preg);
	}		
			
	return 0;
}


