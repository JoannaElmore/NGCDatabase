#include "libecoPCR/ecoPCR.h"
#include "libthermo/nnparams.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <getopt.h>


#define VERSION "1.0.1"


/* ----------------------------------------------- */
/* printout help                                   */                                           
/* ----------------------------------------------- */
#define PP fprintf(stdout, 

static void PrintHelp()
{
        PP      "------------------------------------------\n");
        PP      " ecoPCR Version %s\n", VERSION);
        PP      "------------------------------------------\n");
        PP      "synopsis : searching for sequence and taxonomy hybriding with given primers\n");
        PP      "usage: ecoPCR [options] <nucleotidic patterns>\n");
        PP      "------------------------------------------\n");
        PP      "options:\n");
        PP      "-a    : Salt concentration in M for Tm computation (default 0.05 M)\n\n");
        PP      "-c    : Consider that the database sequences are [c]ircular\n\n");
        PP      "-d    : [D]atabase : to match the expected format, the database\n");
        PP      "        has to be formatted first by the ecoPCRFormat.py program located.\n");
        PP      "        in the tools directory.\n");
        PP      "        ecoPCRFormat.py creates three file types :\n");
        PP      "            .sdx : contains the sequences\n");
        PP      "            .tdx : contains information concerning the taxonomy\n");
        PP      "            .rdx : contains the taxonomy rank\n\n");
        PP      "        ecoPCR needs all the file type. As a result, you have to write the\n");
        PP      "        database radical without any extension. For example /ecoPCRDB/gbmam\n\n");
        PP      "-D    : Keeps the specified number of nucleotides on each side of the in silico \n");
        PP      "        amplified sequences (including the amplified DNA fragment plus the two target \n");
        PP      "        sequences of the primers).\n\n");
        PP      "-e    : [E]rror : max errors allowed by oligonucleotide (0 by default)\n\n");
        PP      "-h    : [H]elp - print <this> help\n\n");
        PP      "-i    : [I]gnore the given taxonomy id.\n");
        PP      "        Taxonomy id are available using the ecofind program.\n");
        PP      "        see its help typing ecofind -h for more information.\n\n");        
        PP      "-k    : [K]ingdom mode : set the kingdom mode\n");
        PP      "        super kingdom mode by default.\n\n");
        PP      "-l    : minimum [L]ength : define the minimum amplication length. \n\n");
        PP      "-L    : maximum [L]ength : define the maximum amplicationlength. \n\n");
        PP      "-m    : Salt correction method for Tm computation (SANTALUCIA : 1\n");
		PP      "        or OWCZARZY:2, default=1)\n\n");
        PP      "-r    : [R]estricts the search to the given taxonomic id.\n");
        PP      "        Taxonomy id are available using the ecofind program.\n");
        PP      "        see its help typing ecofind -h for more information.\n\n");
        PP      "\n");
        PP      "------------------------------------------\n");
        PP      "first argument : oligonucleotide for direct strand\n\n");
        PP      "second argument : oligonucleotide for reverse strand\n\n");
        PP      "------------------------------------------\n");
        PP      "Table result description : \n");
        PP      "column 1 : accession number\n");
        PP      "column 2 : sequence length\n");
        PP      "column 3 : taxonomic id\n");
        PP      "column 4 : rank\n");
        PP      "column 5 : species taxonomic id\n");
        PP      "column 6 : scientific name\n");
        PP      "column 7 : genus taxonomic id\n");
        PP      "column 8 : genus name\n");
        PP      "column 9 : family taxonomic id\n");
        PP      "column 10 : family name\n");
        PP      "column 11 : super kingdom taxonomic id\n");
        PP      "column 12 : super kingdom name\n");
        PP      "column 13 : strand (direct or reverse)\n");
        PP      "column 14 : first oligonucleotide\n");
        PP      "column 15 : number of errors for the first strand\n");
        PP      "column 16 : Tm for hybridization of primer 1 at this site\n");
        PP      "column 17 : second oligonucleotide\n");
        PP      "column 18 : number of errors for the second strand\n");
        PP      "column 19 : Tm for hybridization of primer 1 at this site\n");
        PP      "column 20 : amplification length\n");
        PP      "column 21 : sequence\n");
        PP      "column 22 : definition\n");
        PP      "------------------------------------------\n");
        PP		" https://git.metabarcoding.org/obitools/ecopcr/wikis/home\n");
        PP      "------------------------------------------\n\n");        
        PP      "\n");

}

#undef PP

/* ----------------------------------------------- */
/* printout usage and exit                         */
/* ----------------------------------------------- */

#define PP fprintf(stderr, 

static void ExitUsage(stat)
        int stat;
{
        PP      "usage: ecoPCR [-d database] [-l value] [-L value] [-e value] [-r taxid] [-i taxid] [-k] oligo1 oligo2\n");
        PP      "type \"ecoPCR -h\" for help\n");

        if (stat)
            exit(stat);
}

#undef  PP

void printRepeat(ecoseq_t *seq,
				 char* primer1, char* primer2,
				 PNNParams tparm,
                 PatternPtr o1, PatternPtr o2,
                 char strand, 
                 char kingdom,
                 int32_t pos1, int32_t pos2,
                 int32_t err1, int32_t err2,
                 ecotaxonomy_t *taxonomy,
                 int32_t delta)
{
	char     *AC;
	int32_t  seqlength;
	int32_t	 taxid;
	int32_t  species_taxid;
	int32_t  genus_taxid;
	int32_t  family_taxid;
	int32_t  superkingdom_taxid;
	char     *rank;
	char     *scientificName;
	char     *genus_name;
	char     *family_name;
	char     *superkingdom_name;
	
	ecotx_t  *taxon;
	ecotx_t  *main_taxon;
	
	char     oligo1[MAX_PAT_LEN+1];
	char     oligo2[MAX_PAT_LEN+1];
	
	int32_t  error1;
	int32_t  error2;
	int32_t  ldelta,rdelta;
	
	char     *amplifia = NULL;
	int32_t  amplength;
	double   tm1,tm2;
	double   tm=0;
	
	int32_t i;

	AC = seq->AC;
	seqlength = seq->SQ_length;

	
	main_taxon    = &taxonomy->taxons->taxon[seq->taxid];
	taxid         = main_taxon->taxid;
	scientificName= main_taxon->name;				
	rank          = taxonomy->ranks->label[main_taxon->rank];
	taxon         = eco_getspecies(main_taxon,taxonomy);
	if (taxon)
		{
			species_taxid = taxon->taxid;
			scientificName= taxon->name;				
		}
	else 
		species_taxid = -1;
		
	taxon         = eco_getgenus((taxon) ? taxon:main_taxon,taxonomy);
	if (taxon)
		{
			genus_taxid = taxon->taxid;
			genus_name= taxon->name;				
		}
	else 
		{
			genus_taxid = -1;
			genus_name  = "###";
		}
	
	taxon         = eco_getfamily((taxon) ? taxon:main_taxon,taxonomy);
	if (taxon)
		{
			family_taxid = taxon->taxid;
			family_name= taxon->name;				
		}
	else 
		{
			family_taxid = -1;
			family_name  = "###";
		}
	
	if (kingdom)
		taxon         = eco_getkingdom((taxon) ? taxon:main_taxon,taxonomy);
	else
		taxon         = eco_getsuperkingdom((taxon) ? taxon:main_taxon,taxonomy);
		
	if (taxon)
		{
			superkingdom_taxid = taxon->taxid;
			superkingdom_name= taxon->name;				
		}
	else 
		{
			superkingdom_taxid = -1;
			superkingdom_name  = "###";
		}
	


	ldelta=(pos1 <= delta)?pos1:delta;



	/*rdelta=((pos2+delta)>=seqlength)?seqlength-pos2-1:delta;        */
	rdelta=((pos2+delta)>=seqlength)?seqlength-pos2:delta;

	amplifia = getSubSequence(seq->SQ,pos1-ldelta,pos2+rdelta);
	amplength= strlen(amplifia)-rdelta-ldelta;
	
	if (strand=='R')
	{

	    ecoComplementSequence(amplifia);

		strncpy(oligo1,amplifia + rdelta ,o2->patlen);

		oligo1[o2->patlen]=0;
		error1=err2;

		strncpy(oligo2, amplifia + rdelta + amplength - o1->patlen,o1->patlen);
		oligo2[o1->patlen]=0;
		error2=err1;
		
		if (delta==0)
			amplifia+=o2->patlen;
		else
		{
			delta=ldelta;
			ldelta=rdelta+o2->patlen;
			rdelta=delta+o1->patlen;
		}
	}
	else /* strand == 'D' */
	{
		strncpy(oligo1,amplifia+ldelta,o1->patlen);
		oligo1[o1->patlen]=0;
		error1=err1;
		
		strncpy(oligo2,amplifia + ldelta + amplength - o2->patlen,o2->patlen);
		oligo2[o2->patlen]=0;
		error2=err2;
		
		if (delta==0)
			amplifia+=o1->patlen;
		else
		{
			ldelta+=o1->patlen;
			rdelta+=o2->patlen;
		}
		
	}
	
	
	ecoComplementSequence(oligo2);
	if(delta==0)
		amplifia[amplength - o2->patlen - o1->patlen]=0;
	else
	{
		delta=ldelta+rdelta+amplength-o1->patlen-o2->patlen;
		for (i=0;i<ldelta;i++)
			amplifia[i]|=32;
		for (i=1;i<=rdelta;i++)
			amplifia[delta-i]|=32;

		amplifia[delta]=0;

	}
	
	tm1=nparam_CalcTwoTM(tparm,oligo1,primer1,o1->patlen) - 273.15;
	tm2=nparam_CalcTwoTM(tparm,oligo2,primer2,o2->patlen) - 273.15;
	tm = (tm1 < tm2) ? tm1:tm2;
	printf("%-15s | %9d | %8d | %-20s | %8d | %-30s | %8d | %-30s | %8d | %-30s | %8d | %-30s | %c | %-32s | %2d | %5.2f | %-32s | %2d | %5.2f | %5d | %s | %s\n",
			AC,
			seqlength,
			taxid,
			rank,
			species_taxid,
			scientificName,
			genus_taxid,
			genus_name,
			family_taxid,
			family_name,
			superkingdom_taxid,
			superkingdom_name,
			strand,
			oligo1,
			error1,
			tm1,
			oligo2,
			error2,
			tm2,
			amplength - o1->patlen - o2->patlen,
			amplifia,
			seq->DE
		  );

}

int main(int argc, char **argv) 
{
	ecoseq_t      *seq;
	ecotaxonomy_t *taxonomy;
	char          *scname;
	char          head[11];
	char          tail[11];
	
	int           carg;
	
	char          *oligo1=NULL;
	char          *oligo2=NULL;
	
	PatternPtr    o1;
	PatternPtr    o2;
	PatternPtr    o1c;
	PatternPtr 	  o2c;
	
	int32_t       delta=0;
	int32_t       lmin=0;
	int32_t       lmax=0;
	int32_t       error_max=0;
	int32_t       errflag=0;
	char          kingdom_mode=0;
	
	char          *prefix = NULL;
	
	int32_t       checkedSequence = 0;
	int32_t       positiveSequence= 0;
	int32_t       amplifiatCount  = 0;
	
	int32_t       o1Hits;
	int32_t       o2Hits;
	int32_t       o1cHits;
	int32_t       o2cHits;
	
	int32_t		  begin;
	int32_t       length;
	
	SeqPtr        apatseq=NULL;
	StackiPtr     stktmp;
	
	int32_t       i;
	int32_t       j;
	int32_t       posi;
	int32_t		  posj;
	int32_t       erri;
	int32_t		  errj;
	
	int32_t		  *restricted_taxid = NULL;
	int32_t       *ignored_taxid	= NULL;
	int32_t		  r=0;
	int32_t		  g=0;	
	int32_t		  circular=0;
	
	int32_t		  saltmethod=SALT_METHOD_SANTALUCIA;
	double		  salt=0.05;
	CNNParams     tparm;

    while ((carg = getopt(argc, argv, "hcd:l:L:e:i:r:km:a:tD:")) != -1) {
    	
     switch (carg) {
                                /* -------------------- */
        case 'd':               /* database name        */
                                /* -------------------- */
           prefix = ECOMALLOC(strlen(optarg)+1,
                              "Error on prefix allocation");
           strcpy(prefix,optarg);
           break;

                                /* -------------------- */
        case 'h':               /* help                 */
                                /* -------------------- */
           PrintHelp();
           exit(0);
           break;

           	   	   	   	   	    /* ------------------------- */
        case 'D':               /* min amplification lenght  */
        						/* ------------------------- */
        	sscanf(optarg,"%d",&delta);
        	break;

                                /* ------------------------- */
        case 'l':               /* min amplification lenght  */
                                /* ------------------------- */
           sscanf(optarg,"%d",&lmin);
           break;

                                /* -------------------------- */
        case 'L':               /* max amplification lenght   */
                                /* -------------------------- */
          sscanf(optarg,"%d",&lmax);
          break;

                                /* -------------------- */
        case 'e':               /* error max            */
                                /* -------------------- */
          sscanf(optarg,"%d",&error_max);
          break;
          		
          						/* -------------------- */
        case 'k':               /* set the kingdom mode */
          kingdom_mode = 1;		/* -------------------- */
          break;
         
         						/* ------------------------------------------ */
        case 'r':               /* stores the restricting search taxonomic id */
        						/* ------------------------------------------ */
          restricted_taxid = ECOREALLOC(restricted_taxid,sizeof(int32_t)*(r+1),
          									"Error on restricted_taxid reallocation");
          sscanf(optarg,"%d",&restricted_taxid[r]);
          r++;		
          break;
          
         						/* --------------------------------- */
        case 'i':               /* stores the taxonomic id to ignore */
        						/* --------------------------------- */
		  ignored_taxid = ECOREALLOC(ignored_taxid,sizeof(int32_t)*(g+1),
          									"Error on excluded_taxid reallocation");
          sscanf(optarg,"%d",&ignored_taxid[g]);
          g++;
          break;
          
                                /* -------------------- */
        case 'c':               /* stores the taxonomic id to ignore */
        						/* --------------------------------- */
		  circular = 1;
          break;
          

					/* --------------------------------- */
		case 'm':               /* set salt method                   */
					/* --------------------------------- */
		sscanf(optarg,"%d",&(saltmethod));
		break;

					/* --------------------------------- */
		case 'a':               /* set salt 	         */
					/* --------------------------------- */
		sscanf(optarg,"%lf",&(salt));
		break;

		case '?':               /* bad option           */
                                /* -------------------- */
            errflag++;
   		}

	}
	
	/**
	 * check the path to the database is given as last argument
	 */
	if ((argc -= optind) == 2)
	{
		
		oligo1 = ECOMALLOC(strlen(argv[optind])+1,
                              "Error on oligo1 allocation");
		strcpy(oligo1,argv[optind]);
		optind++;
		oligo2 = ECOMALLOC(strlen(argv[optind])+1,
                              "Error on oligo1 allocation");
		strcpy(oligo2,argv[optind]);
		
		if (circular)
		{
			circular = strlen(oligo1);
			if (strlen(oligo2)>(size_t)circular)
				circular = strlen(oligo2);
		}
	}
	else
		errflag++;
	
	if (prefix == NULL)
	{
		prefix = getenv("ECOPCRDB");
		if (prefix == NULL)
			errflag++;
	}
	
	nparam_InitParams(&tparm,DEF_CONC_PRIMERS,
				 	  DEF_CONC_PRIMERS,
					  salt,
					  saltmethod);
	                  
    if (!oligo1 || !oligo2)
    		errflag++;
	
	if (errflag)
		ExitUsage(errflag);
		
	o1 = buildPattern(oligo1,error_max);
	o2 = buildPattern(oligo2,error_max);
	
	o1c = complementPattern(o1);
	o2c = complementPattern(o2);
	
	printf("#@ecopcr-v2\n");
	printf("#\n");
	printf("# ecoPCR version %s\n",VERSION);
	printf("# direct  strand oligo1 : %-32s ; oligo2c : %32s\n", o1->cpat,o2c->cpat);
	printf("# reverse strand oligo2 : %-32s ; oligo1c : %32s\n", o2->cpat,o1c->cpat);
	printf("# max error count by oligonucleotide : %d\n",error_max);
	
	double tm,tm1,tm2;

	tm1=nparam_CalcSelfTM(&tparm,o1->cpat,o1->patlen) - 273.15;
	tm2=nparam_CalcSelfTM(&tparm,o2->cpat,o2->patlen) - 273.15;
	tm = (tm1 < tm2) ? tm1:tm2;

	printf("# optimal Tm for primers 1 : %5.2f\n",tm1);
	printf("# optimal Tm for primers 2 : %5.2f\n",tm2);

	printf("# database : %s\n",prefix);
	if (lmin && lmax)
		printf("# amplifiat length between [%d,%d] bp\n",lmin,lmax);
	else if (lmin)
		printf("# amplifiat length larger than %d bp\n",lmin);
	else if (lmax)
		printf("# amplifiat length smaller than %d bp\n",lmax);
	if (kingdom_mode)
		printf("# output in kingdom mode\n");
	else
		printf("# output in superkingdom mode\n");
	if (circular)
		printf("# DB sequences are considered as circular\n");
	else
		printf("# DB sequences are considered as linear\n");
	printf("#\n");

	taxonomy = read_taxonomy(prefix,0);

	seq = ecoseq_iterator(prefix);
		
	checkedSequence = 0;
	positiveSequence= 0;
	amplifiatCount  = 0;
	
	while(seq)
	{	
		checkedSequence++;
		/**
		* check if current sequence should be included
		**/
		if ( (r == 0) || 
		     (eco_is_taxid_included(taxonomy, 
		                            restricted_taxid, 
		                            r, 
		                            taxonomy->taxons->taxon[seq->taxid].taxid)
		     ) 
		   )
		     if ((g == 0) || 
		         !(eco_is_taxid_included(taxonomy, 
		                                 ignored_taxid, 
		                                 g, 
		                                 taxonomy->taxons->taxon[seq->taxid].taxid)
		          )
		        )
		     {
		
				//scname = taxonomy->taxons->taxon[seq->taxid].name;
				//strncpy(head,seq->SQ,10);
				//head[10]=0;
				//strncpy(tail,seq->SQ+seq->SQ_length-10,10);
				//tail[10]=0;
		
				apatseq=ecoseq2apatseq(seq,apatseq,circular);
				
				o1Hits = ManberAll(apatseq,o1,0,0,apatseq->seqlen+apatseq->circular);
				o2cHits= 0;
				
				if (o1Hits)
				{
					stktmp = apatseq->hitpos[0];
					begin = stktmp->val[0] + o1->patlen;
					
					if (lmax)
						length= stktmp->val[stktmp->top-1] + o1->patlen - begin + lmax + o2->patlen;
					else
						length= apatseq->seqlen - begin;
					
					if (circular)
					{
						begin = 0;
						length=apatseq->seqlen+circular;
					}	
					o2cHits = ManberAll(apatseq,o2c,1,begin,length);
		
					if (o2cHits)
						for (i=0; i < o1Hits;i++)
						{
							posi = apatseq->hitpos[0]->val[i];
														
							if (posi < apatseq->seqlen)
							{
								erri = apatseq->hiterr[0]->val[i];
								for (j=0; j < o2cHits; j++)
								{
									posj  =apatseq->hitpos[1]->val[j];
									
									if (posj < apatseq->seqlen)
									{
										posj+=o2c->patlen;
										// printf("coucou %d %d %d\n",posi,posj,apatseq->seqlen);
										errj = apatseq->hiterr[1]->val[j];
										length = 0;
										if (posj > posi)
											length = posj - posi - o1->patlen - o2->patlen;
										if (posj < posi)
											length = posj + apatseq->seqlen - posi - o1->patlen - o2->patlen;
										if ((length>0) &&	// For when primers touch or overlap
											(!lmin || (length >= lmin)) &&
											(!lmax || (length <= lmax)))
										{
										    printRepeat(seq,oligo1,oligo2,&tparm,o1,o2c,'D',kingdom_mode,posi,posj,erri,errj,taxonomy,delta);
											//printf("%s\tD\t%s...%s (%d)\t%d\t%d\t%d\t%d\t%s\n",seq->AC,head,tail,seq->SQ_length,o1Hits,o2cHits,posi,posj,scname);
										}
									}
								}
							}
						}
				}
					
				o2Hits = ManberAll(apatseq,o2,2,0,apatseq->seqlen);
				o1cHits= 0;
				if (o2Hits)
				{
					stktmp = apatseq->hitpos[2];
					begin = stktmp->val[0] + o2->patlen;
					
					if (lmax)
						length= stktmp->val[stktmp->top-1] + o2->patlen - begin + lmax + o1->patlen;
					else
						length= apatseq->seqlen - begin;
						
					if (circular)
					{
						begin = 0;
						length=apatseq->seqlen+circular;
					}	

					o1cHits = ManberAll(apatseq,o1c,3,begin,length);
					
					if (o1cHits)
						for (i=0; i < o2Hits;i++)
						{
							posi = apatseq->hitpos[2]->val[i];

							if (posi < apatseq->seqlen)
							{
								erri = apatseq->hiterr[2]->val[i];
								for (j=0; j < o1cHits; j++)
								{
									posj=apatseq->hitpos[3]->val[j];
									if (posj < apatseq->seqlen)
									{
										posj+=o1c->patlen;
										errj=apatseq->hiterr[3]->val[j];

										length = 0;
										if (posj > posi)
											length = posj - posi + 1  - o2->patlen - o1->patlen; /* - o1->patlen : deleted by <EC> (prior to the OBITools3) */
										if (posj < posi)
											length = posj + apatseq->seqlen - posi - o1->patlen - o2->patlen;
										if ((length>0) &&	// For when primers touch or overlap
											(!lmin || (length >= lmin)) &&
											(!lmax || (length <= lmax)))
										{
										    printRepeat(seq,oligo1,oligo2,&tparm,o2,o1c,'R',kingdom_mode,posi,posj,erri,errj,taxonomy,delta);
											//printf("%s\tR\t%s...%s (%d)\t%d\t%d\t%d\t%d\t%s\n",seq->AC,head,tail,seq->SQ_length,o2Hits,o1cHits,posi,posj,scname);
										}
									}
								}
							}
						}
				}	
				
		} /* End of taxonomic selection */
		
		delete_ecoseq(seq);
		
		seq = ecoseq_iterator(NULL);
	}
	
	ECOFREE(restricted_taxid, "Error: could not free restricted_taxid\n");
	ECOFREE(ignored_taxid, "Error: could not free excluded_taxid\n");
		
	return 0;
}
