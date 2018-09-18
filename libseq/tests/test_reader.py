import libdna
import logging

class TestReader(TestCase):
    def read(self):
        bam = '/ifs/scratch/cancer/Lab_RDF/abh2138/ChIP_seq/data/samples/hg19/rdf/katia/CB4_BCL6_RK040/reads/CB4_BCL6_RK040_hg19.sorted.rmdup.bam' #sys.argv[1]
        chr = 'chr1' #sys.argv[2]

        #bc = BinCountWriter(bam)
        #bc.write_all()
        #bc.create(chr)
    
        br = BinCountReader('/ifs/scratch/cancer/Lab_RDF/abh2138/ChIP_seq/data/samples/hg19/rdf/katia/CB4_BCL6_RK040/reads', 'hg19')
        print(br.get_counts('chr1:100000-190900', bin_width=100))
