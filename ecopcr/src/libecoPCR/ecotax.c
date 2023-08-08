#include "ecoPCR.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

static ecotx_t *readnext_ecotaxon(FILE *f,ecotx_t *taxon);

 /** 
 * Open the taxonomy database 
 * @param	pointer to the database (.tdx file) 
 * @return	a ecotxidx_t structure 
 */
ecotxidx_t     *read_taxonomyidx(const char *filename,const char *filename2)
{
	int32_t      count;
	int32_t      count2;
	FILE         *f;
	FILE         *f2;
	ecotxidx_t *index;
	int32_t      i;
	
	f  = open_ecorecorddb(filename,&count,1);
	f2 = open_ecorecorddb(filename2,&count2,0);


	
	index = (ecotxidx_t*) ECOMALLOC(sizeof(ecotxidx_t) + sizeof(ecotx_t) * (count+count2-1),
	                                  "Allocate taxonomy");
	 
	index->count=count+count2;

	fprintf(stderr,"Reading %d taxa...\n",count);
	for (i=0; i < count; i++){
		readnext_ecotaxon(f,&(index->taxon[i]));
		index->taxon[i].parent=index->taxon + (int32_t)index->taxon[i].parent;
	}					

	if (count2>0)
		fprintf(stderr,"Reading %d local taxa...\n",count2);
	else
		fprintf(stderr,"No local taxon\n");

	for (i=0; i < count2; i++){
		readnext_ecotaxon(f2,&(index->taxon[count+i]));
		index->taxon[count+i].parent=index->taxon + (int32_t)index->taxon[count+i].parent;
	}

	return index;
}


int32_t delete_taxonomy(ecotxidx_t *index)
{
	int32_t i;
	
	if (index)
	{
		for (i=0; i< index->count; i++)
			if (index->taxon[i].name)
				ECOFREE(index->taxon[i].name,"Free scientific name");
				
		ECOFREE(index,"Free Taxonomy");
		
		return 0;
	}
	
	return 1;
}



int32_t delete_taxon(ecotx_t *taxon)
{
	if (taxon)
	{
		if (taxon->name)
			ECOFREE(taxon->name,"Free scientific name");
			
		ECOFREE(taxon,"Free Taxon");
		
		return 0;
	}
		
	return 1;
}


/**
 * Read the database for a given taxon a save the data 
 * into the taxon structure(if any found)
 * @param	*f	pointer to FILE type returned by fopen
 * @param	*taxon	pointer to the structure
 * 
 * @return	a ecotx_t structure if any taxon found else NULL 
 */
ecotx_t *readnext_ecotaxon(FILE *f,ecotx_t *taxon)
{
	
	ecotxformat_t *raw;
	int32_t  rs;
	
	raw = read_ecorecord(f,&rs);
	
	if (!raw)
		return NULL;

	if (is_big_endian())
	{
		raw->namelength = swap_int32_t(raw->namelength);
		raw->parent     = swap_int32_t(raw->parent);
		raw->rank       = swap_int32_t(raw->rank);
		raw->taxid      = swap_int32_t(raw->taxid); 
	}
	
	taxon->parent = (ecotx_t*)raw->parent;
	taxon->taxid  = raw->taxid;
	taxon->rank   = raw->rank;
	
	taxon->name   = ECOMALLOC((raw->namelength+1) * sizeof(char),
	                          "Allocate taxon scientific name");
	                          
	strncpy(taxon->name,raw->name,raw->namelength);
		
	return taxon;
}


ecotaxonomy_t    *read_taxonomy(const char *prefix,int32_t readAlternativeName)
{
	ecotaxonomy_t *tax;
	char          *filename;
	char          *filename2;
	int           buffsize;
	
	tax = ECOMALLOC(sizeof(ecotaxonomy_t),
	                "Allocate taxonomy structure");
	
	buffsize = strlen(prefix)+10;
	
	filename = ECOMALLOC(buffsize,
	                     "Allocate filename");
	filename2= ECOMALLOC(buffsize,
	                     "Allocate filename");
	
	snprintf(filename,buffsize,"%s.rdx",prefix);
	
	tax->ranks = read_rankidx(filename);
	
	snprintf(filename,buffsize,"%s.tdx",prefix);
	snprintf(filename2,buffsize,"%s.ldx",prefix);
	
	tax->taxons = read_taxonomyidx(filename,filename2);
	
	if (readAlternativeName)
	{
  	   snprintf(filename,buffsize,"%s.ndx",prefix);
	   tax->names=read_nameidx(filename,tax);
	}
	else
	   tax->names=NULL;
	return tax;
	
}



int32_t delete_ecotaxonomy(ecotaxonomy_t *taxonomy)
{
	if (taxonomy)
	{
		if (taxonomy->ranks)
			ECOFREE(taxonomy->ranks,"Free rank index");
			
		if (taxonomy->taxons)
			ECOFREE(taxonomy->taxons,"Free taxon index");
			
		ECOFREE(taxonomy,"Free taxonomy structure");
		
		return 0;
	}
	
	return 1;
}

ecotx_t *eco_findtaxonatrank(ecotx_t *taxon,
                                int32_t rankidx)
{
	ecotx_t *current_taxon;
	ecotx_t *next_taxon;
	
	current_taxon = taxon;
	next_taxon    = current_taxon->parent;
	
	while ((current_taxon!=next_taxon) &&  // I' am the root node
		   (current_taxon->rank!=rankidx))
		   {
		   	current_taxon = next_taxon;
		   	next_taxon    = current_taxon->parent;
		   }
		   
	if (current_taxon->rank==rankidx)
		return current_taxon;
	else
		return NULL;
}

/**
 * Get back information concerning a taxon from a taxonomic id
 * @param 	*taxonomy 	the taxonomy database
 * @param	taxid		the taxonomic id 
 * 
 * @result	a ecotx_t structure containing the taxonimic information  
 **/
ecotx_t *eco_findtaxonbytaxid(ecotaxonomy_t *taxonomy, 
							  int32_t taxid)
{
	ecotx_t    *current_taxon;
	int32_t     taxoncount;
	int32_t     i;
	
	taxoncount=taxonomy->taxons->count;
	
	for (current_taxon=taxonomy->taxons->taxon,
	     i=0;
	     i < taxoncount;
	     i++,
	     current_taxon++){
	     if (current_taxon->taxid==taxid){
	     	return current_taxon;
	     }
	 }
	
	return (ecotx_t*)NULL;	
}

/**
 * Find out if taxon is son of other taxon (identified by its taxid)
 * @param	*taxon son 		taxon
 * @param	parent_taxid 	taxonomic id of the other taxon
 * 
 * @return 	1 is the other taxid math a parent taxid, else 0
 **/
int eco_isundertaxon(ecotx_t *taxon, 
						int other_taxid)
{
	ecotx_t *next_parent;
	
	next_parent = taxon->parent;	
	
	while ( (other_taxid != next_parent->taxid) && 
			(strcmp(next_parent->name, "root")) )
	{
		next_parent = next_parent->parent;
	}
		
	if (other_taxid == next_parent->taxid)
		return 1;
	else
		return 0;
}

ecotx_t *eco_getspecies(ecotx_t *taxon,
						ecotaxonomy_t *taxonomy)
{
	static ecotaxonomy_t *tax=NULL;
	static int32_t		 rankindex=-1;
	
	if (taxonomy && tax!=taxonomy)
	{
		rankindex = rank_index("species",taxonomy->ranks);
		tax=taxonomy;
	}
		
	if (!tax || rankindex < 0)
		ECOERROR(ECO_ASSERT_ERROR,"No taxonomy defined");
		
	return eco_findtaxonatrank(taxon,rankindex);
}

ecotx_t *eco_getgenus(ecotx_t *taxon,
						ecotaxonomy_t *taxonomy)
{
	static ecotaxonomy_t *tax=NULL;
	static int32_t		 rankindex=-1;
	
	if (taxonomy && tax!=taxonomy)
	{
		rankindex = rank_index("genus",taxonomy->ranks);
		tax=taxonomy;
	}
		
	if (!tax || rankindex < 0)
		ECOERROR(ECO_ASSERT_ERROR,"No taxonomy defined");
		
	return eco_findtaxonatrank(taxon,rankindex);
}


ecotx_t *eco_getfamily(ecotx_t *taxon,
						ecotaxonomy_t *taxonomy)
{
	static ecotaxonomy_t *tax=NULL;
	static int32_t		 rankindex=-1;
	
	if (taxonomy && tax!=taxonomy)
	{
		rankindex = rank_index("family",taxonomy->ranks);
		tax=taxonomy;
	}
		
	if (!tax || rankindex < 0)
		ECOERROR(ECO_ASSERT_ERROR,"No taxonomy defined");
		
	return eco_findtaxonatrank(taxon,rankindex);
}

ecotx_t *eco_getkingdom(ecotx_t *taxon,
						ecotaxonomy_t *taxonomy)
{
	static ecotaxonomy_t *tax=NULL;
	static int32_t		 rankindex=-1;
	
	if (taxonomy && tax!=taxonomy)
	{
		rankindex = rank_index("kingdom",taxonomy->ranks);
		tax=taxonomy;
	}
		
	if (!tax || rankindex < 0)
		ECOERROR(ECO_ASSERT_ERROR,"No taxonomy defined");
		
	return eco_findtaxonatrank(taxon,rankindex);
}

ecotx_t *eco_getsuperkingdom(ecotx_t *taxon,
						ecotaxonomy_t *taxonomy)
{
	static ecotaxonomy_t *tax=NULL;
	static int32_t		 rankindex=-1;
	
	if (taxonomy && tax!=taxonomy)
	{
		rankindex = rank_index("superkingdom",taxonomy->ranks);
		tax=taxonomy;
	}
		
	if (!tax || rankindex < 0)
		ECOERROR(ECO_ASSERT_ERROR,"No taxonomy defined");
		
	return eco_findtaxonatrank(taxon,rankindex);
}
