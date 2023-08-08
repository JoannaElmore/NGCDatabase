#!/usr/bin/env python

import struct
import sys
import os
import gzip


#####
#
# Generic file function
#
#####

class Filter(object):
    
    
    def __init__(self,path):
        self._path = path
        self._taxonFile =  "%s.tdx" % self._path
        self._ranksFile =  "%s.rdx" % self._path
        self._namesFile =  "%s.ndx" % self._path
        self._taxonomy, self._index, self._ranks, self._name = self.__readNodeTable()


    def __universalOpen(self,file):
        if isinstance(file,str):
            if file[-3:] == '.gz':
                rep = gzip.open(file)
            else:
                rep = open(file)
        else:
            rep = file
        return rep

    def __universalTell(self,file):
        if isinstance(file, gzip.GzipFile):
            file=file.myfileobj
        return file.tell()
    
    def __fileSize(self,file):
        if isinstance(file, gzip.GzipFile):
            file=file.myfileobj
        pos = file.tell()
        file.seek(0,2)
        length = file.tell()
        file.seek(pos,0)
        return length
    
    def __progressBar(self,pos,max,reset=False,delta=[]):
        if reset:
            del delta[:]
        if not delta:
            delta.append(time.time())
            delta.append(time.time())
    
        delta[1]=time.time()
        elapsed = delta[1]-delta[0]
        percent = float(pos)/max * 100
        remain = time.strftime('%H:%M:%S',time.gmtime(elapsed / percent * (100-percent)))
        bar = '#' * int(percent/2)
        bar+= '|/-\\-'[pos % 5]
        bar+= ' ' * (50 - int(percent/2))
        sys.stderr.write('\r%5.1f %% |%s] remain : %s' %(percent,bar,remain))




    #####
    #
    # Iterator functions
    #
    #####
    
    
    
    def __ecoRecordIterator(self,file):
        file = self.__universalOpen(file)
        (recordCount,) = struct.unpack('> I',file.read(4))
    
        for i in xrange(recordCount):
            (recordSize,)=struct.unpack('>I',file.read(4))
            record = file.read(recordSize)
            yield record
    
               
    def __ecoNameIterator(self):
        for record in self.__ecoRecordIterator(self._namesFile):
            lrecord = len(record)
            lnames  = lrecord - 16
            (isScientificName,namelength,classLength,indextaxid,names)=struct.unpack('> I I I I %ds' % lnames, record)
            name=names[:namelength]
            classname=names[namelength:]
            yield (name,classname,indextaxid)
    
    
    def __ecoTaxonomicIterator(self):
        for record in self.__ecoRecordIterator(self._taxonFile):
            lrecord = len(record)
            lnames  = lrecord - 16
            (taxid,rankid,parentidx,nameLength,name)=struct.unpack('> I I I I %ds' % lnames, record)
            yield  (taxid,rankid,parentidx,name)
    
    
    def __ecoSequenceIterator(self,file):
        for record in self.__ecoRecordIterator(file):
            lrecord = len(record)
            lnames  = lrecord - (4*4+20)
            (taxid,seqid,deflength,seqlength,cptseqlength,string)=struct.unpack('> I 20s I I I %ds' % lnames, record)
            de = string[:deflength]
            seq = gzip.zlib.decompress(string[deflength:])
            yield  (taxid,seqid,deflength,seqlength,cptseqlength,de,seq)
    
            
    def __ecoRankIterator(self):
        for record in self.__ecoRecordIterator(self._ranksFile):
            yield  record
    
    
    #####
    #
    # Indexes
    #
    #####
    
    def __ecoNameIndex(self):
        indexName = [x for x in self.__ecoNameIterator()]
        return indexName

    def __ecoRankIndex(self):
        rank = [r for r in self.__ecoRankIterator()]
        return rank

    def __ecoTaxonomyIndex(self):
        taxonomy = []
        index = {}
        i = 0;
        for x in self.__ecoTaxonomicIterator():
            taxonomy.append(x)
            index[x[0]] = i 
            i = i + 1
        return taxonomy, index

    def __readNodeTable(self):
        taxonomy, index = self.__ecoTaxonomyIndex()
        ranks = self.__ecoRankIndex()
        name = self.__ecoNameIndex()
        return taxonomy,index,ranks,name


    def findTaxonByTaxid(self,taxid):
        return self._taxonomy[self._index[taxid]]



    #####
    #
    # PUBLIC METHODS
    #
    #####


    def subTreeIterator(self, taxid):
        "return subtree for given taxonomic id "
        idx = self._index[taxid]
        yield self._taxonomy[idx]
        for t in self._taxonomy:
            if t[2] == idx:
                for subt in self.subTreeIterator(t[0]):
                    yield subt
    
    
    def parentalTreeIterator(self, taxid):
        """
           return parental tree for given taxonomic id starting from
           first ancester to the root.
        """
        taxon=self.findTaxonByTaxid(taxid)
        while taxon[2]!= 0: 
            yield taxon
            taxon = self._taxonomy[taxon[2]]
        yield self._taxonomy[0]
    
    
    def ecoPCRResultIterator(self, file):
        "iteration on ecoPCR result file"
        file = self.__universalOpen(file)
        data = ColumnFile(file,
                          sep='|',
                          types=(str,int,int,
                                 str,int,str,
                                 int,str,int,
                                 str,int,str,
                                 str,str,int,
                                 str,int,int,
                                 str,str),skip='#')
        
        for ac, sq_len, taxid,\
            rank, sp_taxid, species,\
            ge_taxid, genus, fa_taxid,\
            family, sk_taxid, s_kgdom,\
            strand, oligo_1, error_1,\
            oligo_2, error_2, amp_len,\
            sq_des, definition in data:
            
            yield {'ac':ac, 'sq_len':sq_len, 'taxid':taxid,
                   'rank':rank, 'sp_taxid':sp_taxid, 'species':species,
                   'ge_taxid':ge_taxid, 'genus':genus, 'fa_taxid':fa_taxid,
                   'family':family, 'sk_taxid':sk_taxid, 's_kgdom':s_kgdom,
                   'strand':strand, 'oligo_1':oligo_1, 'error_1':error_1,
                   'oligo_2':oligo_2, 'error_2':error_2, 'amp_len':amp_len,
                   'sq_des':sq_des, 'definition':definition}
    
    def rankFilter(self,rankid,filter):
        return self._ranks[rankid] == filter


    def lastCommonTaxon(self,taxid_1, taxid_2): 
        t1 = [x[0] for x in self.parentalTreeIterator(taxid_1)]
        t2 = [x[0] for x in self.parentalTreeIterator(taxid_2)]
        t1.reverse()
        t2.reverse()
        count = t1 < t2 and len(t1) or len(t2)
        for i in range(count):
            if t1[i] != t2[i]:
               return t1[i-1]
    
    
    

class ColumnFile(object):
    
    def __init__(self,stream,sep=None,strip=True,types=None,skip=None):
        if isinstance(stream,str):
            self._stream = open(stream)
        elif hasattr(stream,'next'):
            self._stream = stream
        else:
            raise ValueError,'stream must be string or an iterator'
        self._delimiter=sep
        self._strip=strip
        if types:
            self._types=[x for x in types]
            for i in xrange(len(self._types)):
                if self._types[i] is bool:
                    self._types[i]=ColumnFile.str2bool
        else:
            self._types=None
        self._skip = skip
            
    def str2bool(x):
        return bool(eval(x.strip()[0].upper(),{'T':True,'V':True,'F':False}))
                    
    str2bool = staticmethod(str2bool)
            
        
    def __iter__(self):
        return self
    
    def next(self):
        ligne = self._stream.next()
        while ligne[0] == self._skip:
            ligne = self._stream.next()
        data = ligne.split(self._delimiter)
        if self._strip or self._types:
            data = [x.strip() for x in data]
        if self._types:
            it = self.endLessIterator(self._types)
            data = [x[1](x[0]) for x in ((y,it.next()) for y in data)]
        return data
    
    def endLessIterator(self,endedlist):
        for x in endedlist:
            yield x
        while(1):
            yield endedlist[-1]


class Table(list):
    
    def __init__(self, headers, types):
        list.__init__(self)
        self.headers = headers
        self.types = types
        self.lines = []
        
    def printTable(self):
        for h in self.headers:
            print "\t%s\t|" % h,
        print "\n"
        for l in self.lines:
            for c in l:
                print "\t%s\t|" % c
            print "\n"
            
    def getColumn(self,n):
        print "\t%s\n" % self.header[n]
        for i in range(len(self.lines)):
            print "\t%s\n" % i[n]
        
        
        

