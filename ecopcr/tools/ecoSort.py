#!/usr/bin/env python

import struct
import sys
import os
import gzip
import re
import string

from reportlab.graphics.shapes import *
from reportlab.graphics.charts.barcharts import VerticalBarChart
from reportlab.graphics.charts.piecharts import Pie
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.lib import colors
from reportlab.platypus import *


#####
#
# Generic file function
#
#####

class Filter(object):
    """
       This object provides different iterators and method :
        * findTaxonByTaxid
        * subTreeIterator
        * parentalTreeIterator
        * ecoPCRResultIterator
        * rankFilter
        * lastCommonTaxon
        
        see each method for more informations
    """
    
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


    #####
    #
    # PUBLIC METHODS
    #
    #####

    def findTaxonByTaxid(self,taxid):
        """
            Returns a list containing [taxid,rankid,parent_index,nameLength,name]
            It takes one argument : a taxonomic id
        """
        return self._taxonomy[self._index[taxid]]


    def subTreeIterator(self, taxid):
        """
            Returns subtree for given taxid from first child
            to last child. It takes one argument : a taxonomic id
        """
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
        """
           iteration on ecoPCR result file"
           It returns a dictionnary
        """
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
        """
            boolean telling whether rankid match filter
            takes 2 arguments : 
            1- rankid
            2- filter (i.e genus)
        """
        return self._ranks[rankid] == filter


    def lastCommonTaxon(self,taxid_1, taxid_2):
        """
            returns a last common parent for two given taxon.
            It starts from the root and goes down the tree
            until their parents diverge.
            It takes 2 arguments :
            1- taxid 1
            2- taxid 2
        """
        t1 = [x[0] for x in self.parentalTreeIterator(taxid_1)]
        t2 = [x[0] for x in self.parentalTreeIterator(taxid_2)]
        t1.reverse()
        t2.reverse()
        count = t1 < t2 and len(t1) or len(t2)
        for i in range(count):
            if t1[i] != t2[i]:
               return t1[i-1]
        return t1[count-1]
       
    
    
    

class ColumnFile(object):
    """
        cut an ecoPCR file into a list
    """
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
        self._oligo = {}
            
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
        
    def getOligo(self,line):
        if line[2:8] == 'direct':
            r = re.compile('(?<=direct  strand oligo1 : )[A-Z]+(?= *)')
            self._oligo['o1'] = r.findall(line)
        if line[2:9] == 'reverse':
            r = re.compile('(?<=reverse strand oligo2 : )[A-Z]+(?= *)')
            self._oligo['o2'] = r.findall(line) 
        return None
            
            


###########
#
# DATA STRUCTURE
#
###########


class ecoTable(list):
    """
        Data object inheriting from list
    """
    def __init__(self, headers, types):
        list.__init__(self)
        self.headers = headers
        self.types = types


    def __setitem__ (self,key,value):
        """
            method overloaded to check data types
        """
        for e in range(len(value)):
            value[e] = self.types[e](value[e])
        list.__setitem__(self,key,value)
        
    def __getitem__(self,index):
        newtable = ecoTable(self.headers,self.types)
        if isinstance(index,slice):
            newtable.extend(list.__getitem__(self,index))
        else:
            newtable.append(list.__getitem__(self,index))
        
        return newtable
    
    def getColumns(self,columnList):
        newhead = [self.headers[x] for x in columnList]
        newtype = [self.types[x] for x in columnList]
        newtable = ecoTable(newhead,newtype)
        for line in self:
            newtable.append([line[x] for x in columnList])
        
        return newtable


###########
#
# PARSE FUNCTIONS
#
###########

def _parseOligoResult(filter,file,strand):
    seq = {}
    
    if strand == 'direct':
        key = 'oligo_1'
    elif strand == 'reverse':
        key = 'oligo_2'
        
    for s in filter.ecoPCRResultIterator(file):
        o = s[key]
        taxid = s['taxid']
        if not seq.has_key(o):
            seq[o] = [1,taxid]
        else:
            seq[o][0] = seq[o][0] + 1
            seq[o][1] = filter.lastCommonTaxon(seq[o][1],taxid)
    return seq

    
def _parseTaxonomyResult(table):
    tax = {}
    for l in table:
        taxid = l[2]
        scName = l[3]
        count = l[1]
        if not tax.has_key(taxid):
            tax[taxid] = [1,scName,count]
        else:
            tax[taxid][0] = tax[taxid][0] + 1
            tax[taxid][2] = tax[taxid][2] + count
    return tax


def _sortTable(e1,e2):
    e1 = e1[1]
    e2 = e2[1]
    if e1 < e2:
        return 1
    if e1 > e2:
        return -1
    return 0 


def _countOligoMismatch(o1,o2):
    """
        define mismatch between two oligonucleotids
        return number of mismatch
    """ 
    mmatch = 0
    if len(o1) < len(o2):
        mmatch = int(len(o2) - len(o1))
    for i in range(len(o1)):
        try:
            if o1[i] != o2[i]:
                mmatch = mmatch + 1
        except:
            mmatch = mmatch + 1    

    return mmatch

###########
#
# TOOLS FUNCTIONS
#
###########

def customSort(table,x,y):
    """
        
    """
    x = x-1
    y = y-1
    h = (table.headers[x],table.headers[y])
    t = (table.types[x],table.types[y])
    cTable = ecoTable(h,t)
    
    tmp = {}
    
    for l in table:
        if tmp.has_key(l[x]):
            tmp[l[x]] = tmp[l[x]] + l[y]
        else:
            tmp[l[x]] = l[y]
    
    for k,v in tmp.items():
        cTable.append([k,v])
    
    return cTable
    

def countColumnOccurrence(table,x):
    x = x-1
    h = (table.headers[x],"count")
    t = (table.types[x],int)
    cTable = Table(h,t)
    
    tmp = {}
    
    for l in table:
        if tmp.has_key(l[x]):
            tmp[l[x]] = tmp[l[x]] + 1
        else:
            tmp[l[x]] = 1
    
    for k,v in tmp.items():
        cTable.append([k,v])
    
    return cTable


def buildSpecificityTable(table):
    header =  ("mismatch","taxon","count")
    type = (int,str,int)
    speTable = ecoTable(header,type)
    
    tmp = {}
    for l in table:
        if not tmp.has_key(l[5]):
            tmp[l[5]] = {}
        if not tmp[l[5]].has_key(l[3]):
            tmp[l[5]][l[3]] = l[1]
        else:
            tmp[l[5]][l[3]] = tmp[l[5]][l[3]] + l[1]
    
    for mismatch in tmp.items():
        for taxon,count in mismatch[1].items():
            speTable.append([mismatch[0],taxon,count])     

    return speTable
   

def buildOligoTable(table, file, filter, oligoRef, strand='direct'):
    """
        Fills and sorts a Table object with ecoPCR result file
        
        Takes 4 arguments
        1- Table object
        2- ecoPCR result file path
        3- Filter Object
        4- the oligo used in ecoPCR
        5- the oligo type : direct or reverse
        
    """
    seq = _parseOligoResult(filter, file, strand)
    
    i = 0
    for oligo, info in seq.items():
        table.append(0)
        count, lctTaxid = info[0], info[1]
        scName = filter.findTaxonByTaxid(info[1])[3]
        rank = filter._ranks[filter.findTaxonByTaxid(info[1])[1]]
        mismatch = _countOligoMismatch(oligoRef,oligo)
        table[i]=[oligo,count,lctTaxid,scName,rank,mismatch]
        i = i + 1
        
    table.sort(_sortTable)

 
def buildTaxonomicTable(table):
    """
        Fill a Table object with a taxonomic synthesis
    """
    taxHeaders = ("scName","numOfOligo","numOfAmpl","taxid")
    taxTypes = (str, int, int, int)
    taxTable = ecoTable(taxHeaders, taxTypes)
    
    tax = _parseTaxonomyResult(table) 
  
    i = 0
    for taxid, info in tax.items():
        taxTable.append(0)
        numOfOligo, scName, numOfAmpl = info[0], info[1], info[2]
        taxTable[i]=[scName,numOfOligo,numOfAmpl,taxid]
        i = i + 1 

    taxTable.sort(_sortTable)

    return taxTable

def _parseSequenceResult(filter,file,id):
    sequences = {}
    idIndex = {}
    
    if id == 'family':
        key = 'fa_taxid'
    elif id == 'genus':
        key = 'ge_taxid'
    else:
        key = 'taxid'
    
    for s in filter.ecoPCRResultIterator(file):
        seq = s['sq_des']
        id = s[key]
        if not idIndex.has_key(id):
            idIndex[id] = []
        if not sequences.has_key(seq):
            sequences[seq] = [id]
        else:
            sequences[seq].append(id)
    return sequences, idIndex


def _sameValuesInList(array):
    for i in range(len(array)-1):
        if array[i] != array[i+1]:
            return False
    return True


def _sortSequences(file,filter):
    
    sequences, idIndex = _parseSequenceResult(filter,file,'species')
    
    for s,id in sequences.items():
        if len(id) == 1 or _sameValuesInList(id):
            idIndex[id[0]].append(1)
        else:
            for e in id:
                idIndex[e].append(0)
    
   
    for id,values in idIndex.items():
        idIndex[id] = float(values.count(1)) / float(len(values)) * 100

            
    identified = {}
    non_identified = {}
    ambiguous = {}
     
    return sequences, idIndex

def getIntraSpeciesDiversity(table,file,filter):
    
    intraDiv = {} 
    
    seq, idIndex = _sortSequences(file,filter)
    
    for id,percent in idIndex.items():
        if percent == 100:
            intraDiv[id] = [0,[]]
            for seq,idList in sequences.items():
                if id in idList:
                    intraDiv[id][0] = intraDiv[id][0] + 1
                    intraDiv[id][1].append(seq)
                    
    for id, values in intraDiv.items():
        table.append(id,values[0],values[1])
                    
    

###########
#
# OUTPUT FUNCTIONS
#
###########

def printTable(table):
    """
        Displays the content a of Table object
        Take 1 arguments
        1- Table object
    """

    format = ("%20s | " * len(table.headers))[:-3]
    print format % tuple([str(e) for e in table.headers ]) +"\n" + "-"*23*len(table.headers)
    for l in table:
        print format % tuple([str(e) for e in l ])
    print "# %d results" % len(table)
       
        
def saveAsCSV(table,path):
    """
        Creates a csv file from a Table object
        Takes 2 arguments
        1- Table object
        2- path of the file-to-be    
    """
    file = open(path,'w')
    file.write(','.join([str(e) for e in table.headers ]) + "\n")
    for l in table:
        file.write(','.join([str(e) for e in l ]) + "\n")
    file.close()


def grepTable(table,col,pattern):
    """
        Filters a Table object with regular expression
        Takes 3 arguments :
        1- Table object
        2- number of column to match with
        3- regular expression pattern
        
        Returns a Table object
    """
    col = col -1
    p = re.compile(pattern, re.IGNORECASE)
    out = ecoTable(table.headers,table.types)
    for l in table:
        if p.search(l[col]):
            out.append(l)
    return out


###########
#
# GRAPH FUNCTIONS
#
###########

class EcoGraph(object):
    
    def __init__(self):
        self._styles = getSampleStyleSheet()
        
        self._element = []
        self._element.append(self._box(Paragraph("EcoPCR report",  self._styles['Title'])))
        self._element.append(Spacer(0, 0.5 * cm))
    
    def _box(self,flt, center=True):
        box_style = [('BOX', (0, 0), (-1, -1), 0.5, colors.lightgrey)]
        if center:
            box_style += [('ALIGN', (0, 0), (-1, -1), 'CENTER')]
        return Table([[flt]], style=box_style)

    def _addChart(self,chart,title):
        drawing = Drawing(300, 250)
        drawing.add(chart)
        self._element.append(self._box(Paragraph(title, self._styles['Normal'])))
        self._element.append(self._box(drawing))
        self._element.append(Spacer(0, 0.5 * cm))
    
    def _formatData(self,table):
        data, label = [],[]
        for i in range(len(table)):
            label.append(table[i][0])
            data.append(table[i][1])
        return data, label
     
    def makePie(self, table, title):
        data, label = self._formatData(table)
        pie = Pie()
        pie.x = 100
        pie.y = 100
        pie.data = data
        pie.labels = label
        self._addChart(pie, title)
    
    def makeHistogram(self, table, title):    
        data, label = self._formatData(table)
        data = [tuple(data)]
        
        histo = VerticalBarChart() 
        histo.x = 10 
        histo.y = 70 
        histo.height = 150 
        histo.width = 300
        histo.bars.strokeWidth = 1
        histo.barSpacing = 1
        histo.barLabels.dy = 5
        histo.barLabelFormat = '%d'
        histo.barLabels.fontSize = 9 - (len(data[0])/10)
        histo.data = data
        
        histo.categoryAxis.labels.boxAnchor = 'e'
        histo.categoryAxis.labels.textAnchor = 'start'        
        histo.categoryAxis.labels.dx = -40
        histo.categoryAxis.labels.dy = -50
        histo.categoryAxis.labels.angle = 45
        histo.categoryAxis.labels.width = 10
        histo.categoryAxis.labels.height = 4
        histo.categoryAxis.categoryNames = label
        histo.categoryAxis.strokeWidth = 1
        histo.categoryAxis.labels.fontSize = 8
                
        histo.valueAxis.valueMin = min(data[0])*0.7 
        histo.valueAxis.valueMax = max(data[0])
        step = (max(data[0]) - min(data[0])) / 10
        histo.valueAxis.valueStep = step > 1 and step or 1
        
        self._addChart(histo, title)
        
    def makeReport(self,path):
        doc = SimpleDocTemplate(path)       
        doc.build(self._element)

    
###################### 

 
def init():
    file = "/Users/bessiere/Documents/workspace/ecoPCR/src/toto.tmp"
    oligo = {'o1':'ATGTTTAAAA','o2':'ATGGGGGTATTG'}
     
    filter = Filter("/ecoPCRDB/gbmam")
    
    headers = ("oligo", "Num", "LCT Taxid", "Sc Name", "Rank", "distance")
    types = (str, int, int, str, str, int)
    
    o1Table = ecoTable(headers, types)
    o2Table = ecoTable(headers, types)
    
    buildOligoTable(o1Table, file, filter, oligo['o1'], 'direct')
    buildOligoTable(o2Table, file, filter, oligo['o2'], 'reverse')

 
    taxTable = buildTaxonomicTable(o1Table)
    speTable = buildSpecificityTable(o1Table)
       
    return o1Table, o2Table, taxTable



def start():
    file = "/Users/bessiere/Documents/workspace/ecoPCR/src/toto.tmp"
    filter = Filter("/ecoPCRDB/gbmam")
       
    speHeaders = ("taxid","num of seq","list of seq")
    speTypes = (int,int,list) 
    speTable = ecoTable(speHeaders,speTypes)
        
    getIntraSpeciesDiversity(speTable, file, filter)



