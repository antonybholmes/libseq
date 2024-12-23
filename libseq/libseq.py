#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 30 10:23:09 2018

@author: antony
"""

import collections
import math
import numpy as np
import libbam
import libdna
import struct
import os


# SAMTOOLS='/ifs/scratch/cancer/Lab_RDF/abh2138/tools/samtools-1.8/bin/samtools'
# SAMTOOLS = '/ifs/scratch/cancer/Lab_RDF/abh2138/tools/samtools-0.1.19/samtools'

CHRS = [
    "chr1",
    "chr2",
    "chr3",
    "chr4",
    "chr5",
    "chr6",
    "chr7",
    "chr8",
    "chr9",
    "chr10",
    "chr11",
    "chr12",
    "chr13",
    "chr14",
    "chr15",
    "chr16",
    "chr17",
    "chr18",
    "chr19",
    "chr20,",
    "chr21",
    "chr22",
    "chrX",
    "chrY",
    "chrM",
]

BIN_WIDTHS = {100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000}

POWER = {
    100: 2,
    1000: 3,
    10000: 4,
    100000: 5,
    1000000: 6,
    10000000: 7,
    100000000: 8,
    1000000000: 9,
}

MIN_BIN_WIDTH = 100

MAGIC_NUMBER_OFFSET_BYTES = 0
BIN_SIZE_OFFSET_BYTES = MAGIC_NUMBER_OFFSET_BYTES + 4
BIN_WIDTH_OFFSET_BYTES = BIN_SIZE_OFFSET_BYTES + 4
N_BINS_OFFSET_BYTES = BIN_WIDTH_OFFSET_BYTES + 4
BINS_OFFSET_BYTES = N_BINS_OFFSET_BYTES + 4

NO_DATA = np.zeros(0, dtype=int)


class BinCountWriter:
    def __init__(
        self,
        sample,
        bam,
        genome,
        bin_width=100,
        platform="ChIP-seq",
        stat="mean",
        mode="round2",
        outdir="trackbin",
    ):
        self._sample = sample
        self.bam = bam
        self._dir = os.path.dirname(os.path.abspath(bam))

        self._genome = genome
        self._platform = platform
        # self.samtools = samtools
        #self._power = POWER[bin_width]
        self._bin_width = bin_width
        self._stat = stat
        self._mode = mode
        # cache counts
        self._read_map = np.zeros(280000000, dtype=int)
        self._bin_map = collections.defaultdict(int)
        self.sum_c = 0
        self._outdir = outdir  # os.path.join(outdir, genome, self._sample)

    def _reset(self):
        self._read_map.fill(0)
        self._bin_map.clear()

    def _write(self, chr):
        if "_" in chr:
            # only encode official chr
            return

        max_i = np.max(np.where(self._read_map > 0))
        max_bin = math.floor(max_i / self._bin_width)
        bins = max_bin + 1

        block_map = np.zeros(bins, dtype=int)

        i = 0

        for b in range(0, bins):
            if self._stat == "count":
                c = self._bin_map[b]
            elif self._stat == "max":
                c = np.max(self._read_map[i : (i + self._bin_width)])
            else:
                # mean reads per bin
                c = int(
                    round(np.floor(np.mean(self._read_map[i : (i + self._bin_width)])))
                )

            block_map[b] = c

            self.sum_c += c

            i += self._bin_width

        # calculate how many bytes are needed to store each
        # block. We can dynamically vary this to save space
        bin_size_bytes = 1
        maxc = 0

        for c in block_map:
            if c > 255:
                bin_size_bytes = max(bin_size_bytes, 2)

            if c > 65535:
                bin_size_bytes = max(bin_size_bytes, 4)

            if c > maxc:
                maxc = c

        print("Blocks", block_map.size)
        print("Block size", bin_size_bytes, maxc)

        out = os.path.join(
            self._outdir,
            f"{chr}_bw{self._bin_width}_c{self._stat}_{self._genome}.trackbin",
        )

        print(f"Writing to {out}..., block size {bin_size_bytes}")

        with open(out, "wb") as f:
            f.write(struct.pack("=I", 42))
            # Write the bin size in bytes, either 1, 2, or 4
            f.write(struct.pack("=I", bin_size_bytes))
            f.write(struct.pack("=I", self._bin_width))
            # so we know how many bytes are used to represent a count
            # f.write(struct.pack('B', size_in_bytes))
            # f.write(struct.pack('I', max_i))
            # f.write(bytes(read_map))

            f.write(struct.pack("=I", len(block_map)))

            if bin_size_bytes == 1:
                for c in block_map:
                    f.write(struct.pack("=B", c))
            elif bin_size_bytes == 2:
                for c in block_map:
                    f.write(struct.pack("=H", c))
            else:
                # use a full 32bit int
                for c in block_map:
                    f.write(struct.pack("=I", c))

    def _write_sql(self, chr: str):
        if "_" in chr:
            # only encode official chr
            return

        dir = os.path.join(self._outdir, f"bin{self._bin_width}")
        os.makedirs(dir, exist_ok=True)

        max_i = np.max(np.where(self._read_map > 0))
        max_bin = math.floor(max_i / self._bin_width)
        bins = max_bin + 1

        block_map = np.zeros(bins, dtype=int)

        i = 0

        print("writing sql", chr, self._mode, self._stat, self._bin_width)

        for b in range(0, bins):
            if self._stat == "count":
                c = self._bin_map[b]
            elif self._stat == "max":
                c = np.max(self._read_map[i : (i + self._bin_width)])
            else:
                # mean reads per bin
                c = int(round(np.mean(self._read_map[i : (i + self._bin_width)])))

            if self._mode == "round2":
                # round to nearest multiple of 2 so that we reduce
                # bin variation to make smaller bins
                c = int(np.round(c / 2)) * 2

            block_map[b] = c

            self.sum_c += c

            i += self._bin_width

        out = os.path.join(
            dir,
            f"{chr}_bin{self._bin_width}_{self._genome}.sql",
        )

        print(f"Writing to {out}...")

        with open(out, "w") as f:
            print("BEGIN TRANSACTION;", file=f)
            print(
                f"INSERT INTO track (genome, platform, name, chr, bin_width, stat_mode) VALUES ('{self._genome}', '{self._platform}', '{self._sample}', '{chr}', {self._bin_width}, '{self._stat}');",
                file=f,
            )
            print("COMMIT;", file=f)

            print("BEGIN TRANSACTION;", file=f)

            current_count = block_map[0]
            start = 0

            for i, c in enumerate(block_map):
                if c != current_count:
                    if current_count > 0:
                        print(
                            f"INSERT INTO bins (start, end, reads) VALUES ({start}, {i}, {current_count});",
                            file=f,
                        )

                    current_count = c
                    start = i
                    running = False

            # the last bin has an end 1 past the last index of the bins
            print(
                f"INSERT INTO bins (start, end, reads) VALUES ({start}, {len(block_map)}, {current_count});",
                file=f,
            )

            print("COMMIT;", file=f)

    def _write_count(self, reads: int):
        with open(
            os.path.join(self._outdir, f"reads_{self._genome}.txt"),
            "w",
        ) as f:
            print(str(reads), file=f)

    def _write_track_sql(self, reads: int):
        id = f"{self._genome}:{self._platform}:{self._sample}"

        with open(
            os.path.join(self._outdir, "track.sql"),
            "w",
        ) as f:
            print("BEGIN TRANSACTION;", file=f)
            print(
                f"INSERT INTO track (public_id, genome, platform, name, reads, stat_mode) VALUES ('{id}', '{self._genome}', '{self._platform}', '{self._sample}', {reads}, '{self._stat}');",
                file=f,
            )
            print("COMMIT;", file=f)

        # f = open(os.path.join(self.dir, 'bc.reads.{}.{}.txt'.format(self.genome, self.mode)), 'w')
        # f.write(str(reads))
        # f.close()

        # f = open(os.path.join(self.dir, 'counts.{}.{}bw.{}.bc'.format(self.genome, self.power, self.mode)), 'wb')
        # f.write(struct.pack('>I', self.sum_c))
        # f.close()

        # f = open(os.path.join(self.dir, 'bc.counts.{}.{}bw.{}.txt'.format(self.genome, self.power, self.mode)), 'w')
        # f.write(str(self.sum_c))
        # f.close()

    def write(self, chr):
        sam = libbam.BamReader(self.bam)

        # reset the counts
        self._reset()

        self.sum_c = 0

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

            # for i in range(read.pos, read.pos + read.length):
            #    read_map[i] += 1 #.add(read.qname)

            self._read_map[read.pos : (read.pos + read.length)] += 1

            if c % 100000 == 0:
                print("Processed", str(c), "reads...")

            c += 1

        if c == 0:
            return

        self._write(chr)

    def write_all_chr(self):
        sam = libbam.BamReader(self.bam)

        chr = ""
        c = 0

        self.sum_c = 0

        reads = 0

        for read in sam:
            # assume bam is sorted, then
            # when we switch to a new chr, reset
            # counts and begin aggregating againcd
            if read.chr != chr:
                if c > 0:
                    self._write(chr)

                self._reset()
                c = 0
                chr = read.chr

                print("Processing", chr, "...")

            if self._stat == "count":
                sb = math.floor(read.pos / self._bin_width)
                eb = math.floor((read.pos + read.length - 1) / self._bin_width)

                # unique reads in each bin
                for b in range(sb, eb + 1):
                    self._bin_map[b] += 1

            # count reads at a 1bp resolution, which we can aggregate
            # later
            self._read_map[read.pos : (read.pos + read.length)] += 1

            if c % 100000 == 0:
                print("Processed", str(c), "reads...")

            reads += 1
            c += 1

        # Process what is remaining
        if c > 0:
            self._write(chr)

        self._write_count(reads)

    def write_all_chr_sql(self, paired=False):
        # ensure out dir exists
        os.makedirs(self._outdir, exist_ok=True)

        reader = libbam.BamReader(self.bam, paired=paired)

        chrs = reader.chrs()

        chr = ""
        c = 0

        self.sum_c = 0

        reads = 0

        for chr in chrs:
            self._reset()
            c = 0
            print("Processing sql", chr, "...")

            for read in reader.reads(chr):
                if self._stat == "count":
                    sb = math.floor(read.pos / self._bin_width)
                    eb = math.floor((read.pos + read.length - 1) / self._bin_width)

                    # unique reads in each bin
                    for b in range(sb, eb + 1):
                        self._bin_map[b] += 1

                # count reads at a 1bp resolution, which we can aggregate
                # later
                self._read_map[read.pos : (read.pos + read.length)] += 1

                if c % 100000 == 0:
                    print("Processed", str(c), "reads...")

                reads += 1
                c += 1

            self._write_sql(chr)

        self._write_track_sql(reads)


class BinCountReader:
    def __init__(self, dir, genome, mode="max"):
        self._dir = dir
        self.genome = genome
        self.mode = mode

        self.file_map = collections.defaultdict(lambda: collections.defaultdict(str))
        self.bin_size_map = collections.defaultdict(
            lambda: collections.defaultdict(int)
        )
        self._reads = -1

    def _get_file(self, chr, power):
        if power not in self.file_map[chr]:
            # Need to search for exact chr within file otherwise chr1 can map
            # to chr10 etc.
            s = chr + "."
            p = "{}bw".format(power)

            for file in os.listdir(self._dir):
                if (
                    s in file
                    and self.genome in file
                    and p in file
                    and self.mode in file
                ):
                    self.file_map[chr][power] = file  # os.path.join(self.dir, file)
                    break

        return self.file_map[chr][power]

    def _read_data(self, file, n, seek=0):
        f = open(os.path.join(self._dir, file), "rb")
        f.seek(seek)
        data = f.read(n)
        f.close()
        return data

    def _get_magic_num(self, file):
        data = self.read(file, 4)

        # return 42
        return struct.unpack("=I", data)[
            0
        ]  # return int.from_bytes(s, byteorder='little', signed=False)

    def get_magic_num(self, chr, power):
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

        return self._get_magic_num(self._get_file(chr, power))

    def _get_bin_size(self, file):
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

        s = self._read_data(file, 1, seek=BIN_SIZE_OFFSET_BYTES)

        #        f = open(file, 'rb')
        #        f.seek(BIN_SIZE_OFFSET_BYTES)
        #        s = f.read(1)
        #        f.close()

        # return in bytes
        return struct.unpack("=B", s)[0]

    def get_bin_size(self, chr, power):
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
        if power not in self.bin_size_map[chr]:
            self.bin_size_map[chr][power] = self._get_bin_size(
                self._get_file(chr, power)
            )

        return self.bin_size_map[chr][power]

    @property
    def reads(self):
        """
        Return the number of reads in the library.

        Returns
        -------
        int
            The number of reads in the library.
        """

        if self._reads == -1:
            # Cache the number of reads
            file = "reads.{}.{}.bc".format(self.genome, self.mode)

            if os.path.exists(file):
                f = open(file, "rb")
                self._reads = struct.unpack("=I", f.read(4))[0]
                f.close()

        return self._reads

    def _get_bin_width(self, file):
        #        f = open(file, 'rb')
        #        f.seek(BIN_WIDTH_OFFSET_BYTES)
        #        s = f.read(4)
        #        f.close()

        s = self._read_data(file, 4, seek=BIN_WIDTH_OFFSET_BYTES)

        # return in bytes
        return struct.unpack("=I", s)[0]

    def get_bin_width(self, chr):
        return self._get_bin_width(self._get_file(chr))

    def _get_bin_count(self, file):
        #        f = open(file, 'rb')
        #        f.seek(N_BINS_OFFSET_BYTES)
        #        s = f.read(4)
        #        f.close()

        s = self._read_data(file, 4, seek=N_BINS_OFFSET_BYTES)

        # return in bytes
        return struct.unpack("=I", s)[0]

    def get_bin_count(self, chr):
        return self._get_bin_count(self._get_file(chr))

    def _get_counts(self, file, loc, bin_width, bin_size):
        sb = math.floor(loc.start / bin_width)
        eb = math.floor(loc.end / bin_width)
        n = max(1, eb - sb)

        sa = sb
        sn = n

        if bin_size > 1:
            sa *= bin_size
            sn *= bin_size

        #        f = open(file, 'rb')

        # start address of seek

        # print('sa', BINS_OFFSET_BYTES, sa)

        #        f.seek(BINS_OFFSET_BYTES + sa)
        #        # read block of bytes to reduce file io
        #        d = f.read(sn)
        #        f.close()

        d = self._read_data(file, sn, seek=BINS_OFFSET_BYTES + sa)

        ret = np.zeros(n, dtype=int)

        if bin_size == 4:
            d_type = "=I"
        elif bin_size == 2:
            d_type = "=H"
        else:
            d_type = "=B"

        di = 0
        dn = len(d)

        for i in range(0, n):
            if di >= dn:
                break

            ret[i] = struct.unpack(d_type, d[di : (di + bin_size)])[0]
            di += bin_size

        return ret

    def get_counts(self, loc, bin_width):
        """
        Returns the counts of a binned region spanning the genomic
        coordinates.

        Parameters
        ----------
        loc : str or libdna.Loc
            Either a genomic coordinate of the form 'chr1:1-2' or a
            libdna.Loc object
        bin_width : int
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
            return NO_DATA

        #        if bin_width != BIN_WIDTH:
        #            # In order to ensure that the averages are consistent as
        #            # we move around, reset the loc so that the start
        #            # corresponds with the bin. To explain why consider this
        #            # The stored bin width is 100. If we choose a bin width
        #            # of 1000 and slide the window around a little, we may
        #            # pick a start point at 700,800 etc. Since the desired
        #            # bin width is 1000, we will alter which bins are used to
        #            # create the average and thus as we slide the average will
        #            # jump about
        #            s = loc.start // bin_width * bin_width
        #            e = loc.end // bin_width * bin_width
        #
        #            loc = libdna.Loc(loc.chr, s, e)

        power = POWER[bin_width]

        file = self._get_file(loc.chr, power)

        if file is None or file == "":
            return NO_DATA

        bin_size = self.get_bin_size(loc.chr, power)

        # print(file)
        # print(self.get_magic_num(loc.chr))
        # print(bin_size)
        # print(self.get_bin_width(loc.chr))
        # print(self.get_bin_count(loc.chr))

        d = self._get_counts(file, loc, bin_width, bin_size)

        if len(d) == 0:
            return NO_DATA

        if bin_width < MIN_BIN_WIDTH:
            # Take averages when bins are not the same size as the
            # reference, e.g. if bin width is 1000, we take the mean
            # of the 10 bins that will fit in that bin

            sb = math.floor(loc.start / bin_width)
            eb = math.floor(loc.end / bin_width)

            i = 0

            # if we select a window less than BIN_WIDTH, use the BIN_WIDTH
            # window data to pad the smaller bins, thus for example if
            # BIN_WIDTH = 100 and we choose a bin size of 10, if
            # d_BIN_WIDTH[0] = 1 then d_bin[0:10] = 1 since we just
            # duplicate the values for the 10 bins that represent 1
            # large bin. In some sense, this functionality is of little
            # practical use and is discouraged.
            f = math.floor(MIN_BIN_WIDTH / bin_width)
            n = max(1, eb - sb)

            d2 = np.zeros(n, dtype=int)

            for i2 in range(0, n):
                d2[i2] = d[i]

                if (i2 + 1) % f == 0:
                    i += 1

            d = d2

        return d


class S3BinCountReader(BinCountReader):
    def __init__(self, bucket, genome, mode="max"):
        super().__init__(bucket, genome, mode=mode)
        self.fs = s3fs.S3FileSystem(anon=True)

    def _get_file(self, chr, power):
        if power not in self.file_map[chr]:
            # Need to search for exact chr within file otherwise chr1 can map
            # to chr10 etc.
            s = chr + "."
            p = "{}bw".format(power)

            for file in self.fs.ls(self._dir):
                if (
                    s in file
                    and self.genome in file
                    and p in file
                    and self.mode in file
                ):
                    self.file_map[chr][power] = file  # os.path.join(self.dir, file)
                    break

        return self.file_map[chr][power]

    def _read_data(self, file, n, seek=0):
        f = self.fs.open(file, "rb")
        f.seek(seek)
        data = f.read(n)
        return data
