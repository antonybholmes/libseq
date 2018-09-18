#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 30 10:23:09 2018

@author: antony
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 14 12:24:32 2018

@author: antony
"""

import sys
import numpy as np
import libbam
import libdna
import struct
import os

#SAMTOOLS='/ifs/scratch/cancer/Lab_RDF/abh2138/tools/samtools-1.8/bin/samtools'
SAMTOOLS = '/ifs/scratch/cancer/Lab_RDF/abh2138/tools/samtools-0.1.19/samtools'

CHRS = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20,', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM']
BIN_WIDTHS = {10, 100, 1000, 10000, 100000, 1000000, 10000000}
BIN_WIDTH = 100

MAGIC_NUMBER_OFFSET_BYTES = 0
BIN_SIZE_OFFSET_BYTES = MAGIC_NUMBER_OFFSET_BYTES + 4
BIN_WIDTH_OFFSET_BYTES = BIN_SIZE_OFFSET_BYTES + 1
N_BINS_OFFSET_BYTES = BIN_WIDTH_OFFSET_BYTES + 4
BINS_OFFSET_BYTES = N_BINS_OFFSET_BYTES + 4

class BinCountWriter(object):
    def __init__(self, bam, genome, samtools=SAMTOOLS):
        self.__bam = bam
        self.__dir = os.path.dirname(os.path.abspath(bam))
        self.__genome = genome
        self.__samtools = samtools
        # cache counts
        self.__read_map = np.zeros(280000000, dtype=int)
        
    def _reset(self):
        self.__read_map.fill(0)
    
    def _write(self, chr):
        if "_" in chr:
            # only encode official chr
            return
        
        max_i = np.max(np.where(self.__read_map > 0))
        max_bin = max_i // BIN_WIDTH
        bins = max_bin + 1
        
        block_map = np.zeros(bins, dtype=int)
        
        i = 0
        
        for b in range(0, bins):
            c = int(round(np.floor(np.mean(self.__read_map[i:(i + BIN_WIDTH)]))))
            
            block_map[b] = c
            
            i += BIN_WIDTH
            
        bin_size_bits = 8
        maxc = 0
        
        for c in block_map:
            if c > 255:
                bin_size_bits = 16
            
            if c > 65535:
                bin_size_bits = 32
                
            if c > maxc:
                maxc = c
        
        print('Blocks', block_map.size)
        print('Block size', bin_size_bits, maxc)
        
        out = os.path.join(self.__dir, '{}.{}.{}bit.bc'.format(chr, self.__genome, bin_size_bits))
        
        print('Writing to {}...'.format(out))
        
        f = open(out, 'wb')
        f.write(struct.pack('I', 42))
        # Write the bin size in bytes, either 1, 2, or 4
        f.write(struct.pack('B', bin_size_bits // 8))
        f.write(struct.pack('I', BIN_WIDTH))
        # so we know how many bytes are used to represent a count
        #f.write(struct.pack('B', size_in_bytes))
        #f.write(struct.pack('I', max_i))
        #f.write(bytes(read_map))
        
        f.write(struct.pack('I', len(block_map)))
        
        if bin_size_bits == 8:
            for c in block_map:
                f.write(struct.pack('B', c))
        elif bin_size_bits == 16:
            for c in block_map:
                f.write(struct.pack('H', c))
        else:
            for c in block_map:
                f.write(struct.pack('I', c))
                
        f.close()
    
    def write(self, chr):
        sam = libbam.SamReader(self.__bam, samtools=self.__samtools)
        
        # reset the counts
        self._reset()
        
        c = 0
        
        started = False
        
        for read in sam:
            if read.chr == chr:
                started = True
            else:
                if started:
                    # We can stop as we are no longer on this chr
                    break
                else:
                    # until we find the chr, skip processing the read
                    continue
            
            #for i in range(read.pos, read.pos + read.length):
            #    read_map[i] += 1 #.add(read.qname)
            
            self.__read_map[read.pos:(read.pos + read.length)] += 1
            
            if c % 100000 == 0:
                print('Processed', str(c), 'reads...')
          
            c += 1
            
        if c == 0:
            return
        
        self._write(chr)
    
    def write_all(self):
        sam = libbam.SamReader(self.__bam, samtools=self.__samtools)
        
        chr = ''
        c = 0

        for read in sam:
            if read.chr != chr:
                if c > 0:
                    self._write(chr)
                
                self._reset()
                c = 0
                chr = read.chr
                
                print('Processing', chr, '...')
            
            self.__read_map[read.pos:(read.pos + read.length)] += 1
            
            if c % 100000 == 0:
                print('Processed', str(c), 'reads...')
          
            c += 1
            
        # Process what is remaining
        if c > 0:
            self._write(chr)
            

class BinCountReader(object):
    def __init__(self, dir, genome):
        self.__dir = dir
        self.__genome = genome
        self.__file_map = {}
        
    def _get_file(self, chr):
        # Need to search for exact chr within file otherwise chr1 can map
        # to chr10 etc.
        s = chr + '.'
        
        if chr not in self.__file_map:
            for file in os.listdir(self.__dir):
                if s in file and self.__genome in file:
                    self.__file_map[chr] = os.path.join(self.__dir, file)
                    break
                
        return self.__file_map[chr]
    
    @staticmethod
    def _get_magic_num(file):
        f = open(file, 'rb')
        s = f.read(4)
        f.close()
        
        # return 42
        return struct.unpack('I', s)[0] #return int.from_bytes(s, byteorder='little', signed=False)
    
    def get_magic_num(self, chr):
        """
        Return the magic check number 42 for determining if the file is
        encoded correctly.
        
        Parameters
        ----------
        chr : str
            Chromosome of interest.
        
        Return
        ------
        int
            Should return 42. If not, the file is corrupt or not being
            decoded correctly.
        """
        
        return BinCountReader._get_magic_num(self._get_file(chr))
    
    @staticmethod
    def _get_bin_size(file):
        """
        Returns the bin size in bytes
        
        Parameters
        ----------
        file : str
            Path to bin count file.
        
        Returns
        -------
        int
            Bin size in bytes, either 1, 2, or 4.
        """
        
        f = open(file, 'rb')
        f.seek(BIN_SIZE_OFFSET_BYTES)
        s = f.read(1)
        f.close()
        
        # return in bytes
        return struct.unpack('B', s)[0] #return int.from_bytes(s, byteorder='little', signed=False) // 8
    
    def get_bin_size(self, chr):
        """
        Returns the bin size in bytes
        
        Parameters
        ----------
        chr : str
            Chromosome of interest.
        
        Returns
        -------
        int
            Bin size in bytes, either 1, 2, or 4.
        """
        
        return BinCountReader._get_bin_size(self._get_file(chr))
    
    @staticmethod
    def _get_bin_width(file):
        f = open(file, 'rb')
        f.seek(BIN_WIDTH_OFFSET_BYTES)
        s = f.read(4)
        f.close()
        
        # return in bytes
        return struct.unpack('I', s)[0] #return int.from_bytes(s, byteorder='little', signed=False) // 8
    
    def get_bin_width(self, chr):
        return BinCountReader._get_bin_width(self._get_file(chr))
    
    @staticmethod
    def _get_bin_count(file):
        f = open(file, 'rb')
        f.seek(N_BINS_OFFSET_BYTES)
        s = f.read(4)
        f.close()
        
        # return in bytes
        return struct.unpack('I', s)[0] #int.from_bytes(s, byteorder='little', signed=False)
    
    def get_bin_count(self, chr):
        return BinCountReader._get_bin_count(self._get_file(chr))
    
    @staticmethod
    def _get_counts(file, loc, bin_size, bin_width=BIN_WIDTH):
        if bin_width not in BIN_WIDTHS:
            return []
        
        sb = loc.start // bin_width
        eb = loc.end // bin_width
        n = max(1, eb - sb)
        
        sa = sb * bin_size
        sn = n * bin_size
        
        f = open(file, 'rb')
        # start address of seek
        
        #print('sa', BINS_OFFSET_BYTES, sa)
        
        f.seek(BINS_OFFSET_BYTES + sa)
        # read block of bytes to reduce file io
        d = f.read(sn)
        f.close()
        
        ret = np.zeros(n, dtype=int)
        
        if bin_size == 4:
            d_type = 'I'
        elif bin_size == 2:
            d_type = 'H'
        else:
            d_type = 'B'
            
        di = 0
        dn = len(d)
            
        for i in range(0, n):
            if di >= dn:
                break
            
            ret[i] = struct.unpack(d_type, d[di:(di + bin_size)])[0]
            di += bin_size
                
        return ret
    
    def get_counts(self, loc, bin_width=BIN_WIDTH):
        """
        Returns the mean counts of a binned region spanning the genomic
        coordinates.
        
        Parameters
        ----------
        loc : str or libdna.Loc
            Either a genomic coordinate of the form 'chr1:1-2' or a
            libdna.Loc object
        bin : int, optional
            The size of the desired bins within the genomic coordinate.
            If bins are joined together, the average of each bin rounded to
            the nearest int is returned.
        
        Return
        ------
        list
            List of counts in each bin (ints).
        """
        
        
        
        loc = libdna.parse_loc(loc)
        
        if loc is None:
            return []
        
        file = self._get_file(loc.chr)
        bin_size = self.get_bin_size(loc.chr)
        
        #print(file)
        #print(self.get_magic_num(loc.chr))
        #print(bin_size)
        #print(self.get_bin_width(loc.chr))
        #print(self.get_bin_count(loc.chr))
        
        d = BinCountReader._get_counts(file, loc, bin_size)
        
        if len(d) == 0:
            return []
                
        if bin_width != BIN_WIDTH:
            sb = loc.start // bin_width
            eb = loc.end // bin_width
            n = max(1, eb - sb)

            d2 = np.zeros(n, dtype=int)
            
            i = 0
                
            if bin_width > BIN_WIDTH:
                f = bin_width // BIN_WIDTH
                
                # How many sets of f bins we can scan from d, If length of d
                # is not an exact multiple of bin, n will be smaller than
                n = d.size // f
            
                for i2 in range(0, n):
                    d2[i2] = int(np.round(np.mean(d[i:(i + f)])))
                    
                    i += f
            else:
                # if we select a window less than BIN_WIDTH, use the BIN_WIDTH
                # window data to pad the smaller bins, thus for example if
                # BIN_WIDTH = 100 and we choose a bin size of 10, if
                # d_BIN_WIDTH[0] = 1 then d_bin[0:10] = 1 since we just 
                # duplicate the values for the 10 bins that represent 1
                # large bin. In some sense, this functionality is of little
                # practical use and is discouraged.
                f = BIN_WIDTH // bin_width
                
                for i2 in range(0, n):
                    d2[i2] = d[i]
                    
                    if (i2 + 1) % f == 0:
                        i += 1
             
            d = d2
        
        return d
