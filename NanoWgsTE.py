import os
import pandas as pd
import pysam
from Bio import SeqIO

def fq2fa(fastq):
    cnt_reads = 0
    with open (fastq + '.fasta', 'w') as outFasta:
        for seq in SeqIO.parse(fastq, 'fastq'):
            SeqIO.write(seq, outFasta, 'fasta')
            cnt_reads += 1
    print('Number of reads', cnt_reads)
    return(fastq + '.fasta')
def BLASt_vs_Target(reads, target_fasta, blastpath = "blastn", makeblastdb = 'makeblastdb'):
    outFile = reads + "blast_target.tab"
    os.system('{1} -in {0} -out {0} -dbtype nucl'.format(target_fasta, makeblastdb))
    blastCMD = '{3} -query {0} -db {1} -outfmt 6 -num_threads 150 -word_size 11 -out {2}'.format(reads, target_fasta, outFile, blastpath)
    print(blastCMD)
    os.system(blastCMD)
    return(outFile)

def mapppingBamSort(reads, genome_fasta, mm2 = 'minimap2', samtools_path="samtools", bamtools_path="bamtools"):
    os.system('{3} -ax map-ont -t 100 {0} {1} > {2}.sam'.format(genome_fasta, reads, reads, mm2))
    sam_file = '{0}.sam'.format(reads)
    bam_file = sam_file.rsplit('.', 1)[0] + ".bam"
    sort_bam_file = bam_file.rsplit(r'/', 1)[0] + r"/sorted_" + bam_file.rsplit(r'/', 1)[1]

    ##sam to bam
    os.system('{2} view -Sb {0} > {1}'.format(sam_file, bam_file, samtools_path))

    ##sort bam
    os.system('{2} sort -in {0} -out {1}'.format(bam_file, sort_bam_file, bamtools_path))

    ##index
    os.system('{1} index -in {0}'.format(sort_bam_file, bamtools_path))

    return(sort_bam_file)

def selectTargetReads(blastRes, fastq, minBlastHitLen = 500):
    target_fastq_selected = "{}Selected_target_reads.fastq".format(fastq)
    #select blast positive reads
    reads = pd.read_csv(blastRes, sep = "\t", prefix = "V", header = None)

    target_reads = reads[reads['V3'] > minBlastHitLen]
    ids_target_reads = set(target_reads["V0"])
    print('Number of selected reads:', len(ids_target_reads))

    ## create fastq file with selected reads
    cnt_se = 0
    with open(target_fastq_selected,'w') as outFtarget:
        for seq in SeqIO.parse(fastq, 'fastq'):
            if seq.id in ids_target_reads:
                cnt_se += 1
                SeqIO.write(seq, outFtarget, 'fastq')
    print('Number of  reads in out fastq file:', cnt_se)
    return(target_fastq_selected)

def _getInsertion(align, min_indel=1000):
    bases_before_insertion = 0
    # print(align.reference_name, align.reference_start, align.cigartuples)
    for cigar in align.cigartuples:
        if cigar[0] == 1 and cigar[1] > min_indel:
            return (align.reference_start + bases_before_insertion, cigar[1])
        elif cigar[0] == 0 or cigar[0] == 2:
            bases_before_insertion += cigar[1]
    return False


def getSplittedPositions(bam_file, mapping_quality=40, min_len_clipped=1000):
    ### read alignments and collect split coordinates
    split_positions = {}  # chromosome:split positions
    cnt = 0
    cnt_i = 0
    non_cnt = 0
    for algns in pysam.AlignmentFile(bam_file, 'rb'):
        if algns.cigartuples and algns.reference_name != None:
            if algns.mapping_quality >= mapping_quality and not algns.is_supplementary:
                # add chromosome to the dictionary
                if algns.reference_name not in split_positions:
                    split_positions[algns.reference_name] = []
                ##check end of the reads
                if (algns.cigartuples[-1][0] == 5 or algns.cigartuples[-1][0] == 4) and algns.cigartuples[-1][
                    1] > min_len_clipped:
                    split_positions[algns.reference_name].append(algns.reference_end)
                    cnt_i += 1

                ##check start of the reads
                if (algns.cigartuples[0][0] == 5 or algns.cigartuples[0][0] == 4) and algns.cigartuples[0][
                    1] > min_len_clipped:
                    split_positions[algns.reference_name].append(algns.reference_start)
                    cnt_i += 1

                isInsPresent = _getInsertion(algns)
                if isInsPresent:
                    split_positions[algns.reference_name].append(isInsPresent[0])
                    print("INSERTION:", algns.reference_name, isInsPresent[0], 'Size:', isInsPresent[1])

        else:
            non_cnt += 1

    print('Number of splitted ends in the reads', cnt_i)
    print('Number of aligns without reference', non_cnt)
    return (split_positions)

## merge split positions and count
def mergeSplitCount(split_positions, min_distance = 100):
    merged_split_positions = {}
    cnt = 0
    for chromosome in split_positions:
        merged_split_positions[chromosome] = []
        new_splits = []
        split_positions[chromosome] = sorted(split_positions[chromosome])
        for spl_pos in split_positions[chromosome]:
            if not new_splits:
                new_splits.append([spl_pos,spl_pos,1])
            else:
                if abs(new_splits[-1][0] - spl_pos) <= min_distance:
                    new_splits[-1][2] += 1
                    new_splits[-1][1] = spl_pos
                else:
                    cnt += 1
                    print("{4}\t{0}:{1}..{2}\t{3}".format(chromosome,new_splits[-1][0],new_splits[-1][1], new_splits[-1][-1], cnt))
                    new_splits.append([spl_pos,spl_pos,1])
        merged_split_positions[chromosome] = new_splits
    return(merged_split_positions)

def getBed(merged_split_positions, outFile_name):
    #position as sum merged intervals
    with open(outFile_name + ".bed", 'w') as outFile, open(outFile_name + ".readcount", 'w') as outCNT:
        cnt = 0
        for chromosome in merged_split_positions:
            for spl_pos in merged_split_positions[chromosome]:
                start, end = int(spl_pos[0]), int(spl_pos[1])
                id = "{0}:{1}..{2}".format(chromosome, start, end)
                start, end = start - 50, end + 50
                outFile.write("{0}\t{1}\t{2}\t{3}\n".format(chromosome, start, end, id))
                outCNT.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(chromosome, start, end, id, spl_pos[-1] ))

def main(genome_fasta, fastq, target_fasta, outBed, min_blast_hit=500,
         mapping_quality=10, min_len_clipped=500, blastpath='blastn',
         makeblastdb="makeblastdb", minimap2_path='minimap2', samtools_path='samtools', bamtools_path='bamtools'):
    # 1. fq to fasta
    print('########## 1. Fastq to Fasta conversion #######')
    read_fasta = fq2fa(fastq)

    # 2. blast reads vs vasta
    print('########## 2. BLAST reads vs target #######')

    blastRes = BLASt_vs_Target(read_fasta,
                               target_fasta, blastpath = blastpath, makeblastdb = makeblastdb)

    # 3. get fastq of target reads
    print('########## 3. Selection of the reads with similarity to the target #######')
    selected_fastq_reads = selectTargetReads(blastRes, fastq, minBlastHitLen=min_blast_hit)

    # 4. mapping, bam sorting and indexing
    print('########## 4. Mapping the selected reads to the genome #######')
    sorted_bam = mapppingBamSort(selected_fastq_reads, genome_fasta, mm2=minimap2_path,
                                 samtools_path=samtools_path, bamtools_path=bamtools_path)

    ## 5. find split positions
    print('########## 5. Identification of the genomic regions with insertions by capture of clipped reads #######')

    split_positions = getSplittedPositions(sorted_bam, mapping_quality = mapping_quality, min_len_clipped = min_len_clipped)

    ## 6. find insertion sites supported by the clipped reads
    print('########## 6. Count number of reads supporting insertions #######')

    splitted_coords_merged = mergeSplitCount(split_positions, min_distance=100)

    # 7. get BED file with the coordinates
    print('########## 7. Bedfile and count file generation #######')
    getBed(splitted_coords_merged, outBed)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Script to find regions of insertions using Nanopore reads')
    parser.add_argument('genome_fasta', help='path to target sequence in fasta format')
    parser.add_argument('fastq', help='path to fastq file of reads')
    parser.add_argument('target_fasta', help='path to target sequence in fasta format')
    parser.add_argument('outBed', help='path to the output bed file')
    parser.add_argument('-mbh', '--min_blast_hit', help='minimum length of BLAST hit for reads selection', default=500, type=int)
    parser.add_argument('-q', '--map_q', help='minimum mapping quality', default=40,type=int)
    parser.add_argument('-mlc', '--min_len_clipped', help='minimum length of the clipped part', default=500,type=int)
    parser.add_argument('-bp', '--blastn_path', help='path to BLASTn program', default='blastn')
    parser.add_argument('-mdbp', '--makeblastdb_path', help='path to makeblastdb program', default='makeblastdb')
    parser.add_argument('-samtp', '--samtools_path', help='path to samtools program', default='bamtools')
    parser.add_argument('-bamtp', '--bamtools_path', help='path to bamtools program', default='bamtools')
    parser.add_argument('-mm2', '--minimap2_path', help='path to minimap2 program', default="minimap2")


    args = parser.parse_args()

    print(args.fastq)

    main(args.genome_fasta, args.fastq, args.target_fasta, args.outBed, min_blast_hit = args.min_blast_hit,
         mapping_quality=args.map_q, min_len_clipped = args.min_len_clipped, minimap2_path = args.minimap2_path,
         samtools_path = args.samtools_path, bamtools_path = args.bamtools_path )
