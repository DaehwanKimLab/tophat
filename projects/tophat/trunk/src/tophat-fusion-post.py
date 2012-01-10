#!/usr/bin/env python


"""
tophat-fusion-post.py

Created by Daehwan Kim on 2011-05-05.
Copyright (c) 2011 Daehwan Kim. All rights reserved.
"""

import sys, getopt, warnings
import os, subprocess, errno
import string, re
from datetime import datetime, date, time
import math
from sets import Set


use_message = '''
TopHat maps short sequences from spliced transcripts to whole genomes.

Usage:
    tophat-fusion-post [options] <bowtie_index>
    
Options:
    -v/--version
    -o/--output-dir                <string>    [ default: ./tophatfusion_out ]
    --num_fusion_reads             <int>       [ default: 3                  ]
    --num_fusion_pairs             <int>       [ default: 2                  ]
    --num_fusion_both              <int>       [ default: 5                  ]

    --fusion-read-mismatches       <int>       [ default: 2                  ]
    --fusion-multireads            <int>       [ default: 2                  ]

    -p/--num-threads               <int>       [ default: 1                  ]

    --skip-fusion-kmer
    --skip-filter-fusion
    --skip-blast
    --skip-read-dist
    --skip-html

    --tex-table
'''


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

        
output_dir = "./tophatfusion_out/"
logging_dir = output_dir + "logs/"
tmp_dir = output_dir + "tmp/"


class TopHatFusionParams:
    def __init__(self,
                 keep_tmp = ""):
        self.keep_tmp = keep_tmp

        self.num_fusion_reads = 3
        self.num_fusion_pairs = 2
        self.num_fusion_both = 0

        self.fusion_read_mismatches = 2
        self.fusion_multireads = 2

        self.num_threads = 1
        
        self.skip_fusion_kmer = False
        self.skip_filter_fusion = False
        self.skip_blast = False
        self.skip_read_dist = False
        self.skip_html = False

        self.tex_table = False
            
    def check(self):
        if False:
            die("Error: arg to --num-threads must be greater than 0")

    def parse_options(self, argv):
        try:
            opts, args = getopt.getopt(argv[1:], "hp:", 
                                        ["version",
                                         "help",  
                                         "output-dir=",
                                         "num-fusion-reads=",
                                         "num-fusion-pairs=",
                                         "num-fusion-both=",
                                         "fusion-read-mismatches=",
                                         "fusion-multireads=",
                                         "num-threads=",
                                         "skip-fusion-kmer",
                                         "skip-filter-fusion",
                                         "skip-blast",
                                         "skip-read-dist",
                                         "skip-html",
                                         "tex-table"])
        except getopt.error, msg:
            raise Usage(msg)

        for option, value in opts:
            if option in ("-v", "--version"):
                print "TopHat v%s" % (get_version())
                sys.exit(0)
            if option in ("-h", "--help"):
                raise Usage(use_message)
            if option == "--num-fusion-reads":
                self.num_fusion_reads = int(value)
            if option == "--num-fusion-pairs":
                self.num_fusion_pairs = int(value)
            if option == "--num-fusion-both":
                self.num_fusion_both = int(value)
            if option == "--fusion-read-mismatches":
                self.fusion_read_mismatches = int(value)
            if option == "--fusion-multireads":
                self.fusion_multireads = int(value)
            if option in ("-p", "--num-threads"):
                self.num_threads = int(value)
            if option == "--skip-fusion-kmer":
                self.skip_fusion_kmer = True
            if option == "--skip-filter-fusion":
                self.skip_filter_fusion = True
            if option == "--skip-blast":
                self.skip_blast = True
            if option == "--skip-read-dist":
                self.skip_read_dist = True
            if option == "--skip-html":
                self.skip_html = True
            if option == "--tex-table":
                self.tex_table = True

        if len(args) < 1:
            raise Usage(use_message)
        
        return args

    
# Returns the current time in a nice format
def right_now():
    curr_time = datetime.now()
    return curr_time.strftime("%c")


def reverse_complement(seq):
    result = ""
    for nt in seq:
        base = nt
        if nt == 'A':
            base = 'T'
        elif nt == 'a':
            base = 't'
        elif nt == 'C':
            base = 'G'
        elif nt == 'c':
            base = 'g'
        elif nt == 'G':
            base = 'C'
        elif nt == 'g':
            base = 'c'
        elif nt == 'T':
            base = 'A'
        elif nt == 't':
            base = 'a'

        result = base + result

    return result


# Ensures that the output, logging, and temp directories are present. If not, 
# they are created
def prepare_output_dir():
    if os.path.exists(output_dir):
        pass
    else:        
        os.mkdir(output_dir)
        
    if os.path.exists(logging_dir):
        pass
    else:        
        os.mkdir(logging_dir)
        
    if os.path.exists(tmp_dir):
        pass
    else:
        try:
          os.makedirs(tmp_dir)
        except OSError, o:
          die("\nError creating directory %s (%s)" % (tmp_dir, o))

          
def map_fusion_kmer(bwt_idx_prefix, params):
    def get_fusion_seq():
        seq_dic = {}
        for dir in sorted(os.listdir('.')):
            if string.find(dir, "tophat_") != 0:
                continue

            fusion = dir + "/fusions.out"
            if not os.path.exists(fusion):
                continue
        
            fusion_file = open(fusion, 'r')
            fusion_file.readline()
            for line in fusion_file:
                left_seq, right_seq = line[:-1].split('\t@\t')[2:4]
                left_seq = left_seq.split(' ')[0]
                right_seq = right_seq.split(' ')[1]

                if len(left_seq) < 23 or len(right_seq) < 23:
                    continue

                seq_dic[left_seq[-23:]] = 1
                seq_dic[right_seq[:23]] = 1
            
            fusion_file.close()

        fusion_seq_fa = open(output_dir + "fusion_seq.fa", 'w')
        for seq in seq_dic.keys():
            print >> fusion_seq_fa, ">%s" % seq
            print >> fusion_seq_fa, seq
        
        fusion_seq_fa.close()

    def convert_bowtie():
        bwt_dic = {}
        bwtout_file = open(output_dir + "fusion_seq.bwtout", 'r')
        for line in bwtout_file:
            seq, temp, chr, coord = line[:-1].split('\t')[0:4]
            if seq in bwt_dic:
                bwt_dic[seq].append(chr + ":" + coord)
            else:
                bwt_dic[seq] = [chr + ":" + coord]

        bwtout_file.close()

        kmer_map = open(fusion_kmer_file_name, 'w')
        for seq, chrs in bwt_dic.items():
            print >> kmer_map, "%s\t%s" % (seq, ','.join(chrs))

        kmer_map.close()

    print >> sys.stderr, "[%s] Extracting 23-mer around fusions and mapping them using Bowtie" % right_now()

    fusion_kmer_file_name = output_dir + "fusion_seq.map"

    if not os.path.exists(fusion_kmer_file_name):
        get_fusion_seq()
        cmd = ['bowtie', '-p', '8', '-a', '-n', '3', '-m', '100', bwt_idx_prefix, '-f', '%sfusion_seq.fa' % output_dir]
        subprocess.call(cmd, stdout=open(output_dir + 'fusion_seq.bwtout', 'w'), stderr=open('/dev/null', 'w'))
        convert_bowtie()

    
def filter_fusion(params):
    def filter_fusion_impl(fusion, refGene_list, ensGene_list, ensGtp_dic, seq_chr_dic, fusion_gene_list):
        def gene_exists(gene_list, chr, coord, dir, is_left):
            min = 0
            max = len(gene_list) - 1

            while max - min >= 0:
                mid = (min + max) / 2
                gene = gene_list[mid]
                ref_chr = gene[1]

                if chr != ref_chr:
                    if chr < ref_chr:
                        max = mid - 1
                    else:
                        min = mid + 1
                    continue

                left_coord = gene[2]
                right_coord = gene[3]

                if coord >= left_coord and coord <= right_coord:
                    left_coords = gene[5]
                    right_coords = gene[6]
                    sense = gene[7]

                    belong = False
                    where = "outside"
                    for i in range(len(left_coords)):
                        # gives some relax!
                        relax = 3

                        left = int(left_coords[i]) - 1
                        right = int(right_coords[i]) - 1
                        if coord <= right + relax:
                            if coord < left - relax:
                                where = "intron%d(%d-%d)" % (i, int(right_coords[i-1]), left - 1)
                            else:
                                if ((is_left and dir == "f") or (not is_left and dir == "r")) and abs(coord - right) <= relax:
                                    belong = True

                                if ((is_left and dir == "r") or (not is_left and dir == "f")) and abs(coord - left) <= relax:
                                    belong = True

                                where = "exon%d(%d-%d)" % (i + 1, left, right)
                            break

                    return [gene[0], gene[4], where, belong, sense]
                elif coord < left_coord:
                    max = mid - 1
                else:
                    min = mid + 1

            return ["N/A", "N/A", "N/A", False, "N/A"]

        def how_diff(first, second):
            seq_len = len(first)

            min_value = 10000
            prev, curr = [0 for i in range(seq_len)], [0 for i in range(seq_len)]
            for j in range(seq_len):
                for i in range(seq_len):
                    value = 10000
                    if first[i] == second[j]:
                        match = 0
                    else:
                        match = 1

                    # right
                    if i == 0:
                        value = j * 2 + match
                    elif j > 0:
                        value = prev[i] + 2

                    temp_value = 10000

                    # down
                    if j == 0:
                        temp_value = i * 2 + match
                    elif i > 0:
                        temp_value = curr[i-1] + 2

                    if temp_value < value:
                        value = temp_value

                    # match
                    if i > 0 and j > 0:
                        temp_value = prev[i-1] + match

                    if temp_value < value:
                        value = temp_value

                    curr[i] = value

                    if (i == seq_len - 1 or j == seq_len - 1) and value < min_value:
                        min_value = value

                prev, curr = curr, prev

            return min_value

        kmer_len = len(seq_chr_dic.keys()[0])
        sample_name = fusion.split("/")[0][len("tophat_"):]

        data = os.getcwd().split('/')[-1]

        fusion_file = open(fusion, 'r')
        fusion_file.readline()
        for line in fusion_file:
            info, sim, left_seq_org, right_seq_org, left_dist, right_dist, pair_list = line[:-1].split('\t@\t')[:7]

            info = info.split('\t')
            sim = sim.split(' ')
            left_seq = left_seq_org.replace(' ', '')
            right_seq = right_seq_org.replace(' ', '')

            num_reads = int(info[4])
            num_pair_ends = int(info[5])
            num_pair_ends_fusion = int(info[6])
            num_pair_ends_both = int(num_pair_ends + num_pair_ends_fusion * 0.5)
            num_unsupport_reads = int(info[7])
            left_ext = int(info[8])
            right_ext = int(info[9])
            sym = float(info[10])

            half_len = len(left_seq) / 2
            chr1, chr2 = info[0].split('-')[:2]
            coord1, coord2 = int(info[1]), int(info[2])
            dir = info[3]

            if chr1 == "chrM" or chr2 == "chrM":
                continue

            if string.find(sample_name, "single") != -1:
                single = True
            else:
                single = False

            if string.find(data, "maher") != -1:
                extent = min(10, num_reads)
                if single:
                    if left_ext < 25 + extent * 2 or right_ext < 25 + extent * 2:
                        continue
                else:
                    if left_ext < 14 + extent or right_ext < 14 + extent:
                        continue
            else:
                if left_ext < 16 or right_ext < 16:
                    continue

            both = num_reads + num_pair_ends_both
            all = both

            if num_pair_ends > num_reads * 50:
                continue

            # daehwan: num_pair_ends ==> num_pair_ends_both
            if num_reads < params.num_fusion_reads or \
                    num_pair_ends < params.num_fusion_pairs or\
                    both < params.num_fusion_both:
                continue

            if (chr1 != chr2 and num_unsupport_reads > num_reads) or \
                    (chr1 == chr2 and num_unsupport_reads > all + num_pair_ends + 5):
                continue

            pairs = []
            if num_pair_ends >= 1:
                pairs = pair_list.strip().split()

                left, right = pairs[0].split(':')
                if abs(int(left)) + abs(int(right)) > 2000:
                    continue

                # daehwan
                pairs = pairs[:200]

            # are the sequences around the breakpoint different enough?
            if int(sim[0]) < 8:
                continue

            # is the reads distributed symmetrically?
            if sym >= 22 + max(0, 6 - num_reads):
                continue

            max_intron_len = 100000
            if chr1 == chr2 and dir == "ff":
                coord_dif = coord2 - coord1
                if coord_dif > 0 and coord_dif < max_intron_len:
                    continue

            if not left_seq[half_len-kmer_len:half_len] in seq_chr_dic or not right_seq[half_len:half_len+kmer_len] in seq_chr_dic:
                continue

            left_chrs = seq_chr_dic[left_seq[half_len-kmer_len:half_len]]
            right_chrs = seq_chr_dic[right_seq[half_len:half_len+kmer_len]]

            if chr1 == chr2:
                max_intron_len = min(max_intron_len, abs(coord1 - coord2) * 9 / 10)

            same = False
            for chr_coord in left_chrs:
                chr, coord = chr_coord.split(':')
                coord = int(coord)
                if chr == chr2 and abs(coord - coord2) < max_intron_len:
                    same = True
                    break

            if same:
                continue

            for chr_coord in right_chrs:
                chr, coord = chr_coord.split(':')
                coord = int(coord)
                if chr == chr1 and abs(coord - coord1) < max_intron_len:
                    same = True
                    break

            if same:
                continue

            def find_gene(chr, coord, one_dir, is_left):
                result = []
                for gene_list in [refGene_list, ensGene_list]:
                    result.append(gene_exists(gene_list, chr, coord, one_dir, is_left))

                if result[0][0] == "N/A":
                    return result[1] + result[1][:2]
                else:
                    return result[0] + result[1][:2]

            dir = info[3]
            gene1, gene1_name, gene1_where, gene1_belong, gene1_sense, ens_gene1, ens_gene1_name = find_gene(chr1, coord1, dir[0], True)
            gene2, gene2_name, gene2_where, gene2_belong, gene2_sense, ens_gene2, ens_gene2_name = find_gene(chr2, coord2, dir[1], False)

            if ens_gene1_name == ens_gene2_name or ens_gene1 == ens_gene2 or gene1 == "N/A" or gene2 == "N/A" or \
                    (string.find(gene1, "ENS") == 0 and string.find(gene2, "ENS") == 0) or gene1_name == gene2_name:
                continue

            """
            if gene1_name in ensGtp_dic:
                gene1_name = ensGtp_dic[gene1_name]

            if gene2_name in ensGtp_dic:
                gene2_name = ensGtp_dic[gene2_name]
            """

            left_diff = how_diff(left_seq[half_len - 20:half_len], right_seq[half_len - 20:half_len])
            if left_diff <= 8:
                continue

            right_diff = how_diff(left_seq[half_len:half_len+20], right_seq[half_len:half_len+20])
            if right_diff <= 8:
                continue

            if left_diff + right_diff < 20:
                continue

            left_dist = left_dist.strip().split(' ')
            right_dist = right_dist.strip().split(' ')

            for i in range(len(left_dist)):
                left_dist[i] = '%d' % min(9, int(left_dist[i]))
                right_dist[i] ='%d' % min(9, int(right_dist[i]))

            swap = False
            if (dir == 'ff' and gene1_sense == '-' and gene2_sense == '-') or \
                    (dir == 'fr' and gene1_sense == '-' and gene2_sense == '+') or \
                    (dir == 'rf' and gene1_sense == '+' and gene2_sense == '-'):
                swap = True

            if swap:
                if dir == 'ff':
                    dir = 'rr'

                info[0] = "%s-%s" % (chr2, chr1)
                info[1:3] = info[1:3][::-1]
                info[3] = dir
                info[8:10] = info[8:10][::-1]
                left_seq_org, right_seq_org = reverse_complement(right_seq_org), reverse_complement(left_seq_org)
                left_dist, right_dist = right_dist, left_dist
                gene1_name, gene1_where, gene1_sense, gene2_name, gene2_where, gene2_sense = \
                    gene2_name, gene2_where, gene2_sense, gene1_name, gene1_where, gene1_sense
                for j in range(len(pairs)):
                    pair = pairs[j].split(':')
                    pairs[j] = ':'.join(pair[::-1])

            # check if this is due to trans-splicing.
            if gene1_sense in "+-" and gene2_sense in "+-":
                if ((dir == 'ff' or dir == 'rr') and gene1_sense != gene2_sense) or \
                    (dir == 'fr' and (gene1_sense != '+' or gene2_sense != '-')) or \
                    (dir == 'rf' and (gene1_sense != '-' or gene2_sense != '+')):
                    # print >> sys.stderr, "fusion due to trans-splicing", info
                    # print >> sys.stderr, gene1, gene1_name, gene1_where, gene1_sense
                    # print >> sys.stderr, gene2, gene2_name, gene2_where, gene2_sense
                    None

            fusion_gene = []
            fusion_gene.append(sample_name + ' ' + ' '.join(info[:10]))
            fusion_gene.append(left_seq_org)
            fusion_gene.append(right_seq_org)
            fusion_gene.append("%s %s" % (''.join(left_dist[::-1]), ''.join(right_dist)))

            fusion_gene.append("%s %s %s %s" % (gene1_name, gene1_where, gene2_name, gene2_where))
            fusion_gene.append(" ".join(pairs))

            fusion_gene_list.append(fusion_gene)

        fusion_file.close

    print >> sys.stderr, "[%s] Filtering fusions" % right_now()
    
    seq_chr_dic = {}
    seq_chr_file = open(output_dir + "fusion_seq.map", 'r')
    for line in seq_chr_file:
        seq, chrs = line[:-1].split('\t')
        chrs = chrs.split(',')
        seq_chr_dic[seq] = chrs
    
    seq_chr_file.close()

    re_chr = re.compile(r'^chr(\d+|[XYM])$')
    re_mir = re.compile(r'^(MIR)')

    def read_genes(gene_file_name, offset = 1, id = -4):
        gene_list = []
        gene_file = open(gene_file_name, 'r')
        for line in gene_file:
            line = line[:-1].split('\t')[offset:]
            num_exons = int(line[7])
            left_coords = line[8].split(',')[:num_exons]
            right_coords = line[9].split(',')[:num_exons]
            
            if re_chr.findall(line[1]) and not re_mir.findall(line[id]):
                gene_list.append([line[0], line[1], int(line[3]), int(line[4]), line[id], left_coords, right_coords, line[2]])

        gene_file.close()
        return gene_list

    refGene_list = read_genes("refGene_sorted.txt")
    ensGene_list = read_genes("ensGene.txt")

    ensGtp_dic = {}
    ensGtp_file = open("ensGtp.txt", "r")
    for line in ensGtp_file:
        transcript, gene, protein = line[:-1].split()
        ensGtp_dic[transcript] = protein
        ensGtp_dic[gene] = protein
    ensGtp_file.close()

    fusion_gene_list = []
    for file in sorted(os.listdir(".")):
        if string.find(file, "tophat_") != 0:
            continue

        fusion_file = file + "/fusions.out"
        if not os.path.exists(fusion_file):
            continue

        print >> sys.stderr, "\tProcessing:", fusion_file
        filter_fusion_impl(fusion_file, refGene_list, ensGene_list, ensGtp_dic, seq_chr_dic, fusion_gene_list)

    fusion_out_file = output_dir + "potential_fusion.txt"
    print >> sys.stderr, '\t%d fusions are output in %s' % (len(fusion_gene_list), fusion_out_file)
    output_file = open(fusion_out_file, 'w')
    for fusion_gene in fusion_gene_list:
        for line in fusion_gene:
            print >> output_file, line
    
    output_file.close()


def parallel_work(pids, work):
    child = -1
    for i in range(len(pids)):
        if pids[i] == 0:
            child = i
            break

    while child == -1:
        status = os.waitpid(0, 0)
        for i in range(len(pids)):
            if status[0] == pids[i]:
                child = i
                pids[i] = 0
                break

    child_id = os.fork()
    if child_id == 0:
        work()
        os._exit(os.EX_OK)
    else:
        # print >> sys.stderr, '\t\t>> thread %d: %d' % (child, child_id)
        pids[child] = child_id


def wait_pids(pids):
    for pid in pids:
        if pid > 0:
            os.waitpid(pid, 0)
            
    
def do_blast(params):
    print >> sys.stderr, "[%s] Blasting 50-mers around fusions" % right_now()

    file_name = output_dir + "potential_fusion.txt"
    blast_genomic = output_dir + "blast_genomic"
    blast_nt = output_dir + "blast_nt"

    if not os.path.exists(blast_genomic):
        os.system("mkdir %s" % blast_genomic);

    if not os.path.exists(blast_nt):
        os.system("mkdir %s" % blast_nt)

    pids = [0 for i in range(params.num_threads)]

    count = 0
    line_no = 0
        
    check_list = []
    output_list = []
    output = ""
    file = open(file_name, 'r')
    for line in file:
        if line_no % 6 == 0:
            count += 1
            
        if line_no % 6 == 1:
            left_seq = line[:-1].split(" ")[0]

        if line_no % 6 == 2:
            right_seq = line[:-1].split(" ")[1]

        if line_no % 6 == 4:
            def blast(database, seq, outdir):
                file_name = "%s/%s" % (outdir, seq)
                if os.path.exists(file_name):
                    return

                blast_cmd = "echo %s | blastn -db %s -evalue 1e-10 -word_size 28"
                output = os.popen(blast_cmd % (seq, database)).read()

                if str.find(output, "No hits found") != -1:
                    blast_cmd = "echo %s | blastn -db %s -evalue 1e-5"
                    output = os.popen(blast_cmd % (seq, database)).read()
                
                pos1 = str.find(output, ">ref|")
                pos2 = str.find(output, "Database: ", pos1)

                if pos1 != -1 and pos1 < pos2:
                    output = output[pos1:pos2].rstrip()
                else:
                    output = ""

                file = open(file_name, "w")
                file.write(output)
                file.close()

            seq = left_seq + right_seq

            def work():
                blast("blast_human/human_genomic", left_seq, blast_genomic)
                blast("blast_human/human_genomic", right_seq, blast_genomic)
                blast("blast_human/human_genomic", seq, blast_genomic)
                blast("blast_human/nt", left_seq, blast_nt)
                blast("blast_human/nt", right_seq, blast_nt)
                blast("blast_human/nt", seq, blast_nt)
            
            if not os.path.exists(output_dir + "blast_nt/" + seq ):
                print >> sys.stderr,  "\t%d. %s" % (count, line[:-1])
                if params.num_threads <= 1:
                    work()
                else:
                    parallel_work(pids, work)
                
        line_no += 1

    if params.num_threads > 1:
        wait_pids(pids)

        
def read_dist(params):
    def check_region(sample_name, sam, check_list):
        def output_region(sample_name, fusion, reads):
                chr_ref, pos1_ref, pos2_ref, dir_ref = fusion
                chr1_ref, chr2_ref = chr_ref.split('-')

                if len(reads) <= 0:
                    return

                check_file_name = "%s/%s_%s_%d_%d_%s" % (check_dir, sample_name, chr_ref, pos1_ref, pos2_ref, dir_ref)
                check_file = open(check_file_name, "w")

                def my_cmp(a, b):
                    if a[3] and not b[3]:
                        return -1
                    elif not a[3] and b[3]:
                        return 1

                    if a[3]:
                        if dir_ref[0] == "f":
                            return a[4] - b[4]
                        else:
                            return b[4] - a[4]
                    else:
                        if dir_ref[1] == "f":
                            return a[4] - b[4]
                        else:
                            return b[4] - a[4]

                reads = sorted(reads, cmp=my_cmp)
                begin_pos = reads[0][4]
                dist_to_breakpoint = abs(pos1_ref - begin_pos)
                for read in reads:
                    seq_pos = 0
                    read_id, chr1, chr2, before, left_pos, right_pos, cigars, seq, qual, mismatch, left = read[:-1]
                    if (left and left_readid_dic[read_id] > params.fusion_multireads) or\
                            (not left and right_readid_dic[read_id] > params.fusion_multireads):
                        continue

                    space_len = 0
                    if before:
                        space_len = abs(left_pos - begin_pos)
                    else:
                        space_len = abs(pos1_ref - begin_pos) + 2
                        if dir_ref[1] == "f":
                            space_len += (left_pos - pos2_ref)
                        else:
                            space_len += (pos2_ref - left_pos)

                    qual_str = seq_str = ' ' * space_len

                    saw_fusion = False
                    for cigar in cigars:
                        length = int(cigar[:-1])
                        cigar_op = cigar[-1]

                        if cigar_op in "Mm":
                            seq_str += seq[seq_pos:seq_pos + length]
                            qual_str += qual[seq_pos:seq_pos + length]
                            seq_pos += length

                        if cigar_op in "Nn":
                            seq_str += ('|' * length)
                            qual_str += ('|' * length)

                        if cigar_op in "Dd":
                            seq_str += ("D" * length)
                            qual_str += ("D" * length)

                        if cigar_op in "Ii":
                            seq_pos += length

                        if cigar_op in "F":
                            seq_str += " "
                            qual_str += " "

                    cigar_str = ''.join(cigars)
                    prefix_str = '%s %s %d %d %s' % (chr1, chr2, left_pos, right_pos, cigar_str)
                    prefix_str += ' ' * (60 - len(prefix_str))

                    L = "L"
                    if not left:
                        L = "R"

                    begin, end = max(0, dist_to_breakpoint - within), dist_to_breakpoint + within
                    output = '%s %s\n' % (prefix_str, seq_str[begin:end])
                    check_file.write(output)
                    
                check_file.close()

        
        cigar_re = re.compile('\d+\w')
        within = 300

        reads_list = [[] for i in range(len(check_list))]
        left_readid_dic = {}
        right_readid_dic = {}
        sam_file = open(sam, 'r')
        contigs = {}
        for line in sam_file:
            if line[0] == '@':
                if line[1:3] == 'SQ':
                    line = line[:-1].split('\t')
                    contig = line[1].split(':')[1]
                    if contig not in contigs:
                        contigs[contig] = len(contigs)

                continue

            line = line[:-1].split('\t')
            read_id, chr1, left_pos, cigar_str, seq, qual, mismatch = line[0], line[2], int(line[3]) - 1, line[5], line[9], line[10], line[11]
            mismatch = int(mismatch.split(':')[-1])

            secondary_fusion_alignment = False
            for field in line[10:]:
                if field[:5] == "FR:Z:":
                    if field[5] == '2':
                        secondary_fusion_alignment = True
                        break

                    chr1, left_pos, cigar_str, seq, qual = field[7:].split()
                    left_pos = int(left_pos) - 1

            if secondary_fusion_alignment:
                continue

            flag = int(line[1])
            left = (flag & 64 != 0)
            antisense = (flag & 16 != 0)

            if left:
                if not read_id in left_readid_dic:
                    left_readid_dic[read_id] = 1
                else:
                    left_readid_dic[read_id] += 1
            else:
                if not read_id in right_readid_dic:
                    right_readid_dic[read_id] = 1
                else:
                    right_readid_dic[read_id] += 1

            if '-' in chr1:
                chr1, chr2 = chr1.split('-')
            else:
                chr2 = chr1

            right_pos = left_pos
            cigars = cigar_re.findall(cigar_str)
            dir = "ff"
            saw_fusion = False
            for cigar in cigars:
                length = int(cigar[:-1])
                cigar_op = cigar[-1]

                if cigar_op in "MDN":
                    right_pos += length

                    if saw_fusion:
                        dir = dir[0] + "f"
                    else:
                        dir = "f" + dir[1]

                if cigar_op in "mdn":
                    right_pos -= length

                    if saw_fusion:
                        dir = dir[0] + "r"
                    else:
                        dir = "r" + dir[1]

                if cigar_op in "F":
                    if dir[0] == "f":
                        pos1 = right_pos - 1
                    else:
                        pos1 = right_pos + 1

                    pos2 = right_pos = length - 1
                    saw_fusion = True

                if cigar_op in "iIdD":
                    mismatch -= length

            if not saw_fusion:
                report_idx = -1
                for i in range(len(check_list)):
                    chr_ref, pos1_ref, pos2_ref = check_list[i][:3]
                    chr1_ref, chr2_ref = chr_ref.split('-')

                    # a fusion point can be given like "chr9-chr1"
                    if (contigs[chr1_ref] < contigs[chr1] or (contigs[chr1_ref] == contigs[chr1] and left_pos - pos1_ref > 1000000)) and \
                            (contigs[chr2_ref] < contigs[chr1] or (contigs[chr2_ref] == contigs[chr1] and left_pos - pos2_ref > 1000000)):
                        report_idx = i
                    else:
                        break                    

                if report_idx >= 0:
                    for j in range(report_idx+1):
                        output_region(sample_name, check_list[j], reads_list[j])

                    check_list = check_list[report_idx+1:]
                    reads_list = reads_list[report_idx+1:]

                    if len(check_list) <= 0:
                        break

            fusion_read = saw_fusion
            if mismatch > params.fusion_read_mismatches:
                continue

            for i in range(len(check_list)):
                chr_ref, pos1_ref, pos2_ref, dir_ref = check_list[i]
                chr1_ref, chr2_ref = chr_ref.split('-')

                def reverse_read(fusion_read, dir, antisense, left_pos, right_pos, fusion_left, seq, qual, cigars):
                    if dir == 'ff' or not fusion_read:
                        left_pos, right_pos = right_pos - 1, left_pos - 1
                    elif dir == 'fr':
                        left_pos, right_pos = right_pos + 1, left_pos - 1
                    elif dir == 'rf':
                        left_pos, right_pos = right_pos - 1, left_pos + 1

                    reversed_cigars = []
                    cigars.reverse()
                    for i in range(len(cigars)):
                        cigar = cigars[i]
                        opcode = cigar[-1]

                        if opcode == 'F':
                            reversed_cigars.append("%dF" % fusion_left)                  
                        else:
                            if str.islower(opcode):
                                opcode = str.upper(opcode)
                            else:
                                opcode = str.lower(opcode)

                            cigar = cigar[:-1] + opcode
                            reversed_cigars.append(cigar)

                    if dir == 'ff':
                        dir = 'rr'

                    return dir, antisense == False, left_pos, right_pos, reverse_complement(seq), qual[::-1], reversed_cigars

                if fusion_read and chr1 == chr2_ref and chr2 == chr1_ref and pos1 == pos2_ref and pos2 == pos1_ref:
                    dir, antisense, left_pos, right_pos, seq, qual, cigars = \
                        reverse_read(fusion_read, dir, antisense, left_pos, right_pos, pos1, seq, qual, cigars)

                    chr1, chr2 = chr2, chr1
                    pos1, pos2 = pos2, pos1

                if (fusion_read and chr1 == chr1_ref and chr2 == chr2_ref) or \
                        (not fusion_read and (chr1 == chr1_ref or chr1 == chr2_ref)):
                    if not fusion_read:
                        do_not_continue = False
                        if chr1 == chr1_ref:
                            if abs(left_pos - pos1_ref) <= 10000 or abs(right_pos - pos1_ref) <= 10000:
                                do_not_continue = True

                        if chr1 == chr2_ref:
                            if abs(left_pos - pos2_ref) <= 10000 or abs(right_pos - pos2_ref) <= 10000:
                                do_not_continue = True

                        if not do_not_continue:
                            continue

                    seq_temp = seq
                    qual_temp = qual
                    cigars_temp = cigars
                    left_pos_temp, right_pos_temp = left_pos, right_pos

                    before = True
                    if fusion_read:
                        if dir != dir_ref:
                            continue

                        if pos1 != pos1_ref or pos2 != pos2_ref:
                            continue
                    else:
                        do_not_continue = False
                        if chr1 == chr1_ref:
                            if dir_ref[0] == "f" and not (right_pos - pos1_ref >= 10 or pos1_ref - right_pos >= within):
                                do_not_continue = True
                            if dir_ref[0] == "r" and not (pos1_ref - left_pos >= 10 or left_pos - pos1_ref >= within):
                                do_not_continue = True

                        if chr1 == chr2_ref:
                            if dir_ref[1] == "f" and not (pos2_ref - left_pos >= 10 or left_pos - pos2_ref >= within):
                                do_not_continue = True
                                before = False
                            if dir_ref[1] == "r" and not (right_pos - pos2_ref >= 10 or pos2_ref - right_pos >= within):
                                do_not_continue = True
                                before = False

                        if not do_not_continue:
                            continue                    

                        if (before and dir_ref[0] == "r") or (not before and dir_ref[1] == "r"):
                            dir, antisense, left_pos_temp, right_pos_temp, seq_temp, qual_temp, cigars_temp = \
                                reverse_read(fusion_read, dir, antisense, left_pos_temp, right_pos_temp, 0, seq_temp, qual_temp, cigars_temp)


                    reads_list[i].append([read_id, chr1, chr2, before, left_pos_temp, right_pos_temp, cigars_temp, seq_temp, qual_temp, mismatch, left, line])

        for i in range(len(check_list)):
            output_region(sample_name, check_list[i], reads_list[i])

        sam_file.close()


    print >> sys.stderr, "[%s] Generating read distributions around fusions" % right_now()

    file_name = output_dir + "potential_fusion.txt"
    file = open(file_name, 'r')

    check_dir = output_dir + "check"
    check_list = {}
    line_no = 0
    for line in file:
        if line_no % 6 == 0:
            temp_list = line[:-1].split(' ')
            sample_name = temp_list[0]

            if not sample_name in check_list:
                check_list[sample_name] = []

            if not os.path.exists("%s/%s_%s" % (check_dir, sample_name, "_".join(temp_list[1:5]))):
                check_list[sample_name].append(' '.join(temp_list[1:5]))                
            
        line_no += 1

    file.close()

    if not os.path.exists(check_dir):
        os.system("mkdir %s" % check_dir)

    pids = [0 for i in range(params.num_threads)]

    for sample_name, list in check_list.items():
        sam_file_name = 'tophat_%s/accepted_hits.sam' % sample_name
        if not os.path.exists(sam_file_name):
            bam_file_name = 'tophat_%s/accepted_hits.bam' % sample_name
            if not os.path.exists(bam_file_name):
                continue

            cmd = ['samtools', 'view', '-h', bam_file_name]
            subprocess.call(cmd, stdout=open(sam_file_name, 'w'), stderr=open('/dev/null', 'w'))

        increment = 50
        for i in range((len(list) + increment - 1) / increment):
            temp_list = list[i*increment:(i+1)*increment]
            print >> sys.stderr, '\t%s (%d-%d)' % (sample_name, i*increment + 1, min((i+1)*increment, len(list)))

            check_list = []
            for fusion in temp_list:
                print >> sys.stderr, '\t\t%s' % fusion
                fusion = fusion.split()
                check_list.append([fusion[0], int(fusion[1]), int(fusion[2]), fusion[3]])

            def work():
                if len(check_list) > 0:
                    check_region(sample_name, sam_file_name, check_list)

                for j in range(len(temp_list)):
                    list_j = temp_list[j].split(" ")
                    check_file_name = "%s/%s_%s" % (dir, sample_name, "_".join(list_j))
                    check_file = open(check_file_name, "w")
                    check_file.write(check_output_list[j])
                    check_file.close()

            if params.num_threads <= 1:
                work()
            else:
                parallel_work(pids, work)

    if params.num_threads > 1:
        wait_pids(pids)


                
def generate_html(params):
    def blast_output(database, seq):
        blast_output_filename = "%s/%s" % (database, seq)
        if os.path.exists(blast_output_filename):
            file = open(blast_output_filename, "r")
            return file.read() + "\n"

        return ""

    def read_fusion_list(fusion_list):
        re_find_chromosomes = re.compile(r'Homo sapiens chromosome (\d+|[XY])')
        re_find_identities = re.compile(r'Identities = (\d+)\/\d+ \((\d+)%\)')
        re_find_exon = re.compile(r'exon\d+\((\d+-\d+)\)')

        edgren_list = []
        if os.path.exists("fusion_info.txt"):
            edgren_file = open("fusion_info.txt", 'r')
            for fusion in edgren_file:
                sample_name, gene1_name, coord1, gene2_name, coord2 = fusion[:-1].split()

                chr1, coord1, sense1 = coord1.split(':')
                chr1 = "chr" + chr1
                chr2, coord2, sense2 = coord2.split(':')
                chr2 = "chr" + chr2

                left1, right1 = coord1.split('-')
                left2, right2 = coord2.split('-')

                if sense1 == sense2:
                    dir = "ff"
                elif sense1 == '1':
                    dir = "fr"
                else:
                    dir = "rf"

                edgren_list.append([sample_name, gene1_name, chr1, int(left1), int(right1), gene2_name, chr2, int(left2), int(right2), dir])

            edgren_file.close()

        line_no = 0
        do_not_add = False
        output = {}
        file_name = output_dir + "potential_fusion.txt"
        file = open(file_name, 'r')
        for line in file:
            if line_no % 6 == 0:
                if output and not do_not_add:
                    fusion_list.append(output)

                do_not_add = False
                temp_list = line[:-1].split(' ')
                sample_name = temp_list[0]

                chr = temp_list[1]
                chr1, chr2 = chr.split('-')
                left = int(temp_list[2])
                right = int(temp_list[3])
                dir = temp_list[4]

                output = {}
                output["sample_name"] = sample_name
                output["chr"] = chr
                output["chr1"] = chr1
                output["chr2"] = chr2
                output["left_coord"] = left
                output["right_coord"] = right
                output["dir"] = dir
                output["stats"] = temp_list[5:]

                def found_in_fusion_list(chr1, chr2, left, right, dir):
                    for fusion in edgren_list:
                        ref_chr1, ref_left1, ref_right1, ref_chr2, ref_left2, ref_right2, ref_dir = fusion[2:5] + fusion[6:]

                        if dir == "rr":
                            dir = "ff"


                        if dir != ref_dir:
                            second_event = True
                        else:
                            second_event = False

                        if (chr1 == ref_chr1 and chr2 == ref_chr2 and left >= ref_left1 and left <= ref_right1 and right >= ref_left2 and right <= ref_right2) or\
                                (chr1 == ref_chr2 and chr2 == ref_chr1 and left >= ref_left2 and left <= ref_right2 and right >= ref_left1 and right <= ref_right1):
                            return True, second_event

                    return False, False

                output["found_in_list"], output["second_event"] = found_in_fusion_list(chr1, chr2, left, right, dir)

                # daehwan - I had decided to delete the following checking routine as of Jul, 2011.
                if False:
                # if not output["found_in_list"] and chr1 == chr2:
                    # if ((dir == "ff" and right - left < 1000000 and right - left > 0) or (dir == "rr" and left - right < 1000000 and left - right > 0)):
                    if abs(right - left) < 100000:
                        do_not_add = True

            elif line_no % 6 == 1:
                output["left_seq"] = line[:-1].split(" ")[0]

            elif line_no % 6 == 2:
                output["right_seq"] = line[:-1].split(" ")[1]

                left_seq = output["left_seq"]
                right_seq = output["right_seq"]
                seq = left_seq + right_seq

                temp_output = blast_output(output_dir + "blast_genomic", left_seq)
                left_chromosomes = Set(re_find_chromosomes.findall(temp_output))

                temp_output = blast_output(output_dir + "blast_genomic", right_seq)
                right_chromosomes = Set(re_find_chromosomes.findall(temp_output))

                # daehwan
                if string.find(output["sample_name"], "single") == -1:
                    chr1, chr2 = output["chr"].split("-")
                    if chr1 != chr2 and left_chromosomes.intersection(right_chromosomes):
                        do_not_add = True

                    temp_output = blast_output(output_dir + "blast_genomic", seq) + blast_output(output_dir + "blast_nt", seq)
                    for identity in re_find_identities.findall(temp_output):
                        query = int(identity[0])
                        percent = int(identity[1])
                        if query + percent > 160:
                            do_not_add = True
                            break

            elif line_no % 6 == 3:
                output["depth"] = line[:-1]

            elif line_no % 6 == 4:
                output["gene"] = line[:-1]
                temp = output["gene"].split()
                output["gene1"], output["gene2"] = temp[0], temp[2]

            else:
                pair_str = line.strip()
                if len(pair_str) > 0:
                    pairs = pair_str.split(" ")
                else:
                    pairs = []

                output["pair_coords"] = pairs

                read_output = ""
                read_output_file = output_dir + "check/%s_%s_%d_%d_%s" % (output["sample_name"], output["chr"],
                                                                          output["left_coord"], output["right_coord"], output["dir"])

                if not do_not_add and os.path.exists(read_output_file):
                    read_output = open(read_output_file, "r").read()

                    def unique(list, left, right):
                        color_len = 300
                        result, lcolor, rcolor = [], [0 for i in range(color_len)], [0 for i in range(color_len)]
                        keys = []
                        F_start, F_end = 1000000, 0
                        for item in list:
                            key = item.split(" ")
                            if len(key) < 6 :
                                continue

                            chr1, chr2, pos1, pos2, cigars = key[:5]
                            pos1, pos2 = int(pos1), int(pos2)

                            def color(list, left, right):
                                for i in range(left, right):
                                    if i >= len(list):
                                        break

                                    list[i] += 1                            

                            if 'F' in cigars:
                                dist1 = abs(pos1 - left)
                                dist2 = abs(pos2 - right)

                                color(lcolor, 0, dist1)
                                color(rcolor, 0, dist2)

                            elif abs(pos2 - left) < color_len:
                                color(lcolor, abs(pos2 - left), abs(pos1 - left))

                            elif abs(pos1 - right) < color_len:
                                color(rcolor, abs(pos1 - right), abs(pos2 - right))

                            if "F" in key[4] or len(keys) == 0 or key[:2] != keys[-1][:2] or abs(int(key[2]) - int(keys[-1][2])) >= 3:
                                keys.append(key[:3])
                                result.append(item)

                                if "F" in key[4]:
                                    index = len(result) - 1
                                    if index < F_start:
                                        F_start = index

                                    if F_end < index:
                                        F_end = index

                        if F_end == 0:
                            F_start, F_end = 0, len(result)

                        temp_result = result[max(0, F_start - 50):(F_end+50)]
                        result = []
                        for i in range(len(temp_result)):
                            item = temp_result[i].split(" ")
                            if i == 0:
                                trim = 0
                                for ch in item[40:-1]:
                                    trim += 1

                            item = item[:40] + item[40 + trim:]
                            result.append(" ".join(item))

                        return result, lcolor, rcolor

                    if len(read_output) < 1024 * 1024 * 1024:
                        read_output = read_output.split("\n")
                        read_output, lcolor, rcolor = unique(read_output, output["left_coord"], output["right_coord"])

                        def stat(lcolor, rcolor, dist):
                            lcount, lsum, rcount, rsum = 1, 0, 1, 0
                            lgap, lpass, rgap, rpass = 0, False, 0, False
                            for i in range(dist):
                                if lcolor[i] > 0:
                                    lcount += 1
                                    lsum += lcolor[i]

                                    if lgap > 0:
                                        lpass = True
                                else:
                                    if not lpass:
                                        lgap += 1

                                if rcolor[i] > 0:
                                    rcount += 1
                                    rsum += rcolor[i]

                                    if rgap > 0:
                                        rpass = True
                                else:
                                    if not rpass:
                                        rgap += 1

                            if not lpass:
                                lgap = 0

                            if not rpass:
                                rgap = 0

                            return lcount, lsum / lcount, lgap, rcount, rsum / rcount, rgap

                        temp =  output["gene"].split()
                        gene1_loc, gene2_loc = temp[1], temp[3]

                        def find_exon_len(gene_loc, coord, dir, is_left):
                            exon_len = 1000000
                            for loc in re_find_exon.findall(gene_loc):
                                start, end = loc.split('-')
                                start, end = int(start), int(end)
                                if (is_left and dir == 'f') or (not is_left and dir == 'r'):
                                    exon_len = coord - start + 1
                                else:
                                    exon_len = end - coord + 1
                                return exon_len

                            return exon_len

                        lcount_min, rcount_min, diff_max = 150, 150, 120

                        lcount_exon_min = find_exon_len(gene1_loc, output["left_coord"], output["dir"][0], True)
                        rcount_exon_min = find_exon_len(gene2_loc, output["right_coord"], output["dir"][1], False)

                        lcount_min = min(lcount_min, lcount_exon_min - 20)
                        rcount_min = min(rcount_min, rcount_exon_min - 20)
                        diff_max = min(diff_max, abs(lcount_min - rcount_min) + 20)

                        if lcount_exon_min < 1000 and rcount_exon_min < 1000:
                            diff_max = max(diff_max, abs(lcount_exon_min - rcount_exon_min) + 20)                        

                        lcount, lavg, lgap, rcount, ravg, rgap = stat(lcolor, rcolor, len(lcolor))
                        if lcount <= lcount_min or rcount <= rcount_min or lgap / lcount > 0.1 or rgap / rcount > 0.1:
                            if abs(min(lcount, lcount_exon_min) - min(rcount, rcount_exon_min)) > diff_max or lcount < 60 or rcount < 60:
                                do_not_add = True
                                print >> sys.stderr, lcount_min, rcount_min, diff_max
                                if output["found_in_list"]:
                                    print >> sys.stderr, output
                                    print >> sys.stderr, lcount, rcount, abs(lcount - rcount), lgap, rgap

                        def derivation(color, length, avg):
                            der = 0.0
                            for i in range(length):
                                diff = 1.0 - float(color[i]) / float(max(1, avg))
                                der += diff * diff

                            return math.sqrt(der / max(1, length))

                        output["read_output"] = read_output
                        output["lcount"] = lcount
                        output["lavg"] = lavg
                        output["lgap"] = lgap
                        output["lder"] = derivation(lcolor, lcount_min, lavg)
                        output["rcount"] = rcount
                        output["ravg"] = ravg
                        output["rgap"] = rgap
                        output["rder"] = derivation(rcolor, rcount_min, ravg)
                        print >> sys.stderr, "\t%d" % (line_no / 6 + 1)

                    else:
                        read_output = [r"too many - not shown"]

                if not "read_output" in output:
                    output["read_output"] = ""
                    output["lcount"] = 0
                    output["lavg"] = 0
                    output["lgap"] = 0
                    output["lder"] = 0
                    output["rcount"] = 0
                    output["ravg"] = 0
                    output["rgap"] = 0
                    output["rder"] = 0

            line_no += 1

        if output and not do_not_add:
            fusion_list.append(output)

        file.close()
        
    def cluster_fusion(fusion_list, cluster_list):
        cluster_dist = 500000

        cluster_temp_list = []
        parent = [i for i in range(len(fusion_list))]

        for i in range(len(fusion_list)):
            fusion = fusion_list[i]

            cluster = {}
            cluster["index"] = [i]
            cluster["chr"] = fusion["chr"]
            cluster["left1"] = cluster["left2"] = fusion["left_coord"]
            cluster["right1"] = cluster["right2"] = fusion["right_coord"]
            cluster["dir"] = fusion["dir"]
            cluster["found_in_list"] = fusion["found_in_list"]

            cluster_temp_list.append(cluster)

        def parent_index(parent, i):
            parent_i = parent[i]
            if i == parent_i:
                return i
            else:
                parent[i] = parent[parent_i]
                return parent_index(parent, parent[i])

        for i in range(len(fusion_list) - 1):
            parent_i = parent_index(parent, i)
            for j in range(i+1, len(fusion_list)):
                parent_j = parent_index(parent, j)

                if parent_i == parent_j:
                    continue

                cluster_i = cluster_temp_list[parent_i]
                cluster_j = cluster_temp_list[parent_j]

                if cluster_i["chr"] != cluster_j["chr"] or cluster_i["dir"] != cluster_j["dir"]:
                    continue

                i_left1 = cluster_i["left1"]
                i_left2 = cluster_i["left2"]

                j_left1 = cluster_j["left1"]
                j_left2 = cluster_j["left2"]

                if abs(i_left1 - j_left1) > cluster_dist or abs(i_left2 - j_left2) > cluster_dist or\
                        abs(i_left2 - j_left1) > cluster_dist or abs(i_left1 - j_left2) > cluster_dist:
                    continue

                i_right1 = cluster_i["right1"]
                i_right2 = cluster_i["right2"]

                j_right1 = cluster_j["right1"]
                j_right2 = cluster_j["right2"]

                if abs(i_right1 - j_right1) > cluster_dist or abs(i_right2 - j_right2) > cluster_dist or\
                        abs(i_right2 - j_right1) > cluster_dist or abs(i_right1 - j_right2) > cluster_dist:
                    continue

                cluster_i["left1"] = min(i_left1, j_left1)
                cluster_i["left2"] = max(i_left2, j_left2)

                cluster_i["right1"] = min(i_right1, j_right1)
                cluster_i["right2"] = max(i_right2, j_right2)

                parent[j] = i
                cluster_i["index"].extend(cluster_j["index"])

                if cluster_j["found_in_list"]:
                    cluster_i["found_in_list"] = True

        cluster_temp_list2 = []
        for i in range(len(parent)):
            if i == parent[i]:
                cluster_temp_list2.append(cluster_temp_list[i])

        def num_samples(indices):
            samples = []
            for i in indices:
                sample = fusion_list[i]["sample_name"]
                if not sample in samples:
                    samples.append(sample)

            return len(samples)                    

        def cmp(a, b):
            def pair_count(indices):
                final_score = -1000000.0
                for index in indices:
                    dist = 1000000

                    stats = fusion_list[index]["stats"]
                    num_read = int(stats[0])
                    pair = int(stats[1])
                    pair_fusion = int(stats[2])
                    anti_read = int(stats[3])
                    anti_pair = int(stats[4])
                    lcount, lavg, lgap, lder = fusion_list[index]["lcount"], fusion_list[index]["lavg"], fusion_list[index]["lgap"], fusion_list[index]["lder"] 
                    rcount, ravg, rgap, rder = fusion_list[index]["rcount"], fusion_list[index]["ravg"], fusion_list[index]["rgap"], fusion_list[index]["rder"]
                    pair_coords = fusion_list[index]["pair_coords"]
                    if len(pair_coords) > 0:
                        pair = 0
                        for pair_coord in pair_coords:
                            pair_coord = pair_coord.split(":")
                            temp_dist = abs(int(pair_coord[0])) + abs(int(pair_coord[1]))
                            if temp_dist < dist:
                                dist = temp_dist

                            if temp_dist < 2000:
                                pair += 1

                    anti_read += 0.5
                    if pair == 0:
                        rate = num_read / anti_read
                    else:
                        rate = pair / anti_read

                    max_avg = 300
                    score = lcount + rcount + (min(max_avg, lavg) + min(max_avg, ravg)) - abs(lcount - rcount) \
                        - min(max_avg, abs(lavg - ravg)) - (lgap + rgap) - (lder + rder) * max_avg \
                        + rate - min(dist, 1000)

                    if score > final_score:
                        final_score = score

                    fusion_list[index]["score"] = score
                    fusion_list[index]["test"] = "lcount: %d, lavg: %d, lgap: %d, lder: %f, rcount: %d, ravg: %d, rgap: %f, rder: %f, rate: %f, dist: %d, score: %f" % \
                        (lcount, lavg, lgap, lder, rcount, ravg, rgap, rder, rate, dist, score)

                return int(final_score)

            a_indices = a["index"]
            b_indices = b["index"]

            a_gene, b_gene = 0, 0
            def known_genes(indices):
                num = 0
                for index in indices:
                    temp_num = 0
                    fusion = fusion_list[index]

                    if fusion["gene1"] != "N/A":
                        temp_num += 1

                    if fusion["gene2"] != "N/A":
                        temp_num += 1

                    if temp_num > num:
                        num = temp_num

                return num

            a_gene = known_genes(a_indices)
            b_gene = known_genes(b_indices)

            if a_gene != b_gene:
                return b_gene - a_gene

            a_num_sample, b_num_sample = num_samples(a_indices), num_samples(b_indices)
            a_score = pair_count(a_indices)
            b_score = pair_count(b_indices)

            return b_score - a_score

        cluster_temp_list = sorted(cluster_temp_list2, cmp=cmp)

        max_num_fusions = 500
        for i in range(min(max_num_fusions, len(cluster_temp_list))):
            do_not_add = False
            indices = cluster_temp_list[i]["index"]

            if num_samples(indices) > 5:
                do_not_add = True

            if not do_not_add:
                def cmp(a, b):
                    return int(fusion_list[b]["score"] - fusion_list[a]["score"])

                cluster_temp_list[i]["index"] = sorted(cluster_temp_list[i]["index"], cmp=cmp)
                cluster_list.append(cluster_temp_list[i])

    def generate_html_impl(fusion_list, cluster_list):
        html_file_name = output_dir + 'result.html'
        html_file = open(html_file_name, 'w')

        if params.tex_table:
            tex_table_file_name = output_dir + 'table.tex'
            tex_table_file = open(tex_table_file_name, 'w')
            print >> tex_table_file, r'\documentclass{article}'
            print >> tex_table_file, r'\usepackage{graphicx}'
            print >> tex_table_file, r'\begin{document}'
            print >> tex_table_file, r'\pagestyle{empty}'
            print >> tex_table_file, r"\center{\scalebox{0.7}{"
            print >> tex_table_file, r"\begin{tabular}{| c | c | c | c | c | c | c |}"
            print >> tex_table_file, r"\hline"
            print >> tex_table_file, r"SAMPLE ID & Fusion genes (left-right) & Chromosomes (left-right) & " + \
                r"5$'$ position & 3$'$ position & Spanning reads & Spanning pairs \\"

        print >> html_file, r'<HTML>'
        print >> html_file, r'<HEAD>'
        print >> html_file, r'<TITLE>result</TITLE>'
        print >> html_file, r'<META NAME="description" CONTENT="result">'
        print >> html_file, r'<META NAME="keywords" CONTENT="result">'
        print >> html_file, r'<META NAME="resource-type" CONTENT="document">'
        print >> html_file, r'<META NAME="distribution" CONTENT="global">'
        print >> html_file, r'<META HTTP-EQUIV="Content-Style-Type" CONTENT="text/css">'
        print >> html_file, r'</HEAD>'
        print >> html_file, r'<BODY text="#000000" bgcolor="#FFFFFF">'

        print >> html_file, r'<H1><BR>candidate fusion list</H1>'
        print >> html_file, r'The following tables show fusion candidates where fusions are grouped based on their genomic locations (<a href="#detail">table description</a>)<BR>'
        print >> html_file, r'<a href="#sample">sample information</a><BR>'
        print >> html_file, r'Known fusions are found at <a href="#mcl">Mitelman Database</a> as of Feb 11, 2011<P>'

        def cmp(i, j):
            i_left = fusion_list[i]["left_coord"]
            j_left = fusion_list[j]["left_coord"]
            if i_left == j_left:
                return i - j
            else:
                return i_left - j_left

        indices_list = []
        for c in cluster_list:
            indices_list += c["index"]

        indices_list = sorted(indices_list)
        print >> sys.stderr, "\tnum of fusions: %d" % len(cluster_list)

        num_fusions = 0
        for c in range(len(cluster_list)):
            cluster = cluster_list[c]
            indices = cluster["index"]

            num_fusions += len(indices)

            top_index = indices[0]
            indices = sorted(indices, cmp=cmp)

            chr1, chr2 = cluster["chr"].split('-')

            if params.tex_table:
                print >> tex_table_file, r'\hline'

            found_in_list = ""
            if cluster["found_in_list"]:
                found_in_list = "(confirmed)"

            print >> html_file, r'<P><P><P><BR>'
            print >> html_file, r'%d. %s %s %s' % (c+1, cluster["chr"], cluster["dir"], found_in_list)

            print >> html_file, r'<TABLE CELLPADDING=3 BORDER="1">'
            for i in indices:
                fusion = fusion_list[i]
                sample_name = fusion["sample_name"]
                sample_name = string.split(sample_name, '_')[0]

                output = ""
                stats = fusion["stats"]
                stats = [int(stat) for stat in stats]

                if stats[1] > 0:
                    pair_support_html = r"\htmlref{%d}{pair_%d}" % (stats[1], i)
                else:
                    pair_support_html = r"0"

                if params.tex_table:
                    if sample_name == "lane1":
                        sample_name2 = "thyroid"
                    else:
                        sample_name2 = "testes"
                        
                    print >> tex_table_file, r'%s & %s & %s & %d & %d & %d & %d \\' % \
                        (sample_name2,
                         fusion["gene1"] + "-" + fusion["gene2"],
                         chr1[3:] + "-" + chr2[3:],
                         fusion["left_coord"],
                         fusion["right_coord"],
                         stats[0],
                         stats[1] + stats[2])

                print >> html_file, r'<TR><TD ALIGN="LEFT"><a href="#fusion_%d">%s</a></TD>' % (i, sample_name)
                print >> html_file, r'<TD ALIGN="LEFT">%s</TD>' % fusion["gene1"]
                print >> html_file, r'<TD ALIGN="LEFT">%s</TD>' % chr1
                print >> html_file, r'<TD ALIGN="RIGHT">%d</TD>' % fusion["left_coord"]
                print >> html_file, r'<TD ALIGN="LEFT">%s</TD>' % fusion["gene2"]
                print >> html_file, r'<TD ALIGN="LEFT">%s</TD>' % chr2
                print >> html_file, r'<TD ALIGN="RIGHT">%d</TD>' % fusion["right_coord"]
                print >> html_file, r'<TD ALIGN="RIGHT"><a href="#read_%d">%d</a></TD>' % (i, stats[0])
                print >> html_file, r'<TD ALIGN="RIGHT"><a href="#pair_%d">%d</a></TD>' % (i, stats[1])
                print >> html_file, r'<TD ALIGN="RIGHT">%d</TD>' % stats[2]
                # print >> html_file, r'<TD ALIGN="RIGHT">%.2f</TD></TR>' % fusion["score"]
                
            print >> html_file, r'</TABLE>'

        # sample
        print >> html_file, r'<H1><A NAME="sample"></A><BR>sample list</H1>'
        print >> html_file, r'<TABLE CELLPADDING=3 BORDER="1">'
        print >> html_file, r'<TR><TD ALIGN="LEFT">sample name</TD>'
        print >> html_file, r'<TD ALIGN="RIGHT">fragment length</TD>'
        print >> html_file, r'<TD ALIGN="RIGHT">read length</TD>'
        print >> html_file, r'<TD ALIGN="RIGHT"># of fragments</TD></TR>'

        if os.path.exists("sample_info.txt"):
            sample_file = open("sample_info.txt", "r")
            for line in sample_file:
                line = line[:-1].split("\t")
                print >> html_file, r'<TR><TD ALIGN="LEFT">%s</TD>' % line[0]
                print >> html_file, r'<TD ALIGN="RIGHT">%s</TD>' % line[1]
                print >> html_file, r'<TD ALIGN="RIGHT">%s</TD>' % line[2]
                print >> html_file, r'<TD ALIGN="RIGHT">%s</TD></TR>' % line[3]

            sample_file.close()

        print >> html_file, r'</TABLE>'

        # description
        print >> html_file, r'<H1><A NAME="detail"></A><BR>table description</H1>'
        print >> html_file, r'1. Sample name in which a fusion is identified <BR>'
        print >> html_file, r'2. Gene on the "left" side of the fusion <BR>'
        print >> html_file, r'3. Chromosome ID on the left <BR>'
        print >> html_file, r'4. Coordinates on the left <BR>'
        print >> html_file, r'5. Gene on the "right" side <BR>'
        print >> html_file, r'6. Chromosome ID on the right <BR>'
        print >> html_file, r'7. Coordinates on the right <BR>'
        print >> html_file, r'8. Number of spanning reads <BR>'
        print >> html_file, r'9. Number of spanning mate pairs <BR>'
        print >> html_file, r'10. Number of spanning mate pairs where one end spans a fusion <BR>'
        print >> html_file, r'If you follow the the 9th column, it shows coordinates "number1:number2" where one end is located at a distance of "number1" bases from the left genomic coordinate of a fusion and "number2" is similarly defined'

        for i in indices_list:
            fusion = fusion_list[i]
            sample_name = fusion["sample_name"]
            chr = fusion["chr"]
            left = fusion["left_coord"]
            right = fusion["right_coord"]
            dir = fusion["dir"]

            left_seq = fusion["left_seq"]
            right_seq = fusion["right_seq"]
            seq = left_seq + right_seq

            pair_coords = fusion["pair_coords"]

            left_blast_genomic = blast_output(output_dir + "blast_genomic", left_seq)
            right_blast_genomic = blast_output(output_dir + "blast_genomic", right_seq)
            blast_genomic = blast_output(output_dir + "blast_genomic", seq)

            left_blast_nt = blast_output(output_dir + "blast_nt", left_seq)
            right_blast_nt = blast_output(output_dir + "blast_nt", right_seq)
            blast_nt = blast_output(output_dir + "blast_nt", seq)

            print >> html_file, r'<H1></A><A NAME="fusion_%d"></A><BR>%s %s %d %d %s</H1>' % (i, sample_name, chr, left, right, dir)
            print >> html_file, r'<PRE>'
            print >> html_file, r' '.join(fusion["stats"])
            print >> html_file, r'%s %s' % (left_seq, right_seq)
            print >> html_file, r'%s' % fusion["depth"]
            print >> html_file, r'%s' % fusion["gene"]
            print >> html_file, r'</PRE>'

            print >> html_file, r'<H2>blast search - human genome</A></H2>'
            print >> html_file, r'<H3>left</H3>'
            print >> html_file, r'<PRE>%s</PRE>' % left_blast_genomic
            print >> html_file, r'<H3>right</H3>'
            print >> html_file, r'<PRE>%s</PRE>' % right_blast_genomic

            print >> html_file, r'<H2>blast search - nt</A></H2>'
            print >> html_file, r'<H3>left</H3>'
            print >> html_file, r'<PRE>%s</PRE>' % left_blast_nt
            print >> html_file, r'<H3>right</H3>'
            print >> html_file, r'<PRE>%s</PRE>' % right_blast_nt

            print >> html_file, r'<H2><A NAME="read_%d"></A><BR>reads</H2>' % i
            read_output = fusion["read_output"]
            print >> html_file, r'<PRE>%s</PRE>' % '\n'.join(read_output)

            print >> html_file, r'<H2><A NAME="pair_%d"></A><BR>%d pairs</H2>' % (i, len(pair_coords))
            print >> html_file, r'<PRE>%s</PRE>' % '\n'.join(pair_coords)

        if os.path.exists("mcl"):
            print >> html_file, r'<H1><A NAME="mcl"></A><BR>Mitelman Database of Chromosome Aberrations and Gene Fusions in Cancer</H1>'
            print >> html_file, r'<TABLE CELLPADDING=3 BORDER="1">'
            mcl_file = open("mcl", "r")
            for line in mcl_file:
                line = line[:-1].split("\t")
                if len(line) < 8:
                    continue

                print >> html_file, r'<TR><TD ALIGN="LEFT">%s</TD>' % line[0]
                print >> html_file, r'<TD ALIGN="LEFT">%s</TD>' % line[1]
                print >> html_file, r'<TD ALIGN="RIGHT">%s</TD>' % line[2]
                print >> html_file, r'<TD ALIGN="RIGHT">%s</TD>' % line[3]
                print >> html_file, r'<TD ALIGN="LEFT">%s</TD>' % line[4]
                print >> html_file, r'<TD ALIGN="LEFT">%s</TD>' % line[5]
                print >> html_file, r'<TD ALIGN="LEFT">%s</TD>' % line[6]
                print >> html_file, r'<TD ALIGN="LEFT">%s</TD></TR>' % line[7]

            mcl_file.close()
            print >> html_file, r'</TABLE>'

        print >> html_file, r'</BODY>'
        print >> html_file, r'</HTML>'
        html_file.close()
        
        if params.tex_table:
            print >> tex_table_file, r"\hline"
            print >> tex_table_file, r"\end{tabular}}}"
            print >> tex_table_file, r'\end{document}'
            tex_table_file.close()
            os.system("pdflatex --output-directory=%s %s" % (output_dir, tex_table_file_name))

    print >> sys.stderr, "[%s] Reporting final fusion candidates in html format" % right_now()

    fusion_list = []
    read_fusion_list(fusion_list)

    cluster_list = []
    cluster_fusion(fusion_list, cluster_list)

    generate_html_impl(fusion_list, cluster_list)
    

# Format a DateTime as a pretty string.  
# FIXME: Currently doesn't support days!
def formatTD(td):
  hours = td.seconds // 3600
  minutes = (td.seconds % 3600) // 60
  seconds = td.seconds % 60
  return '%02d:%02d:%02d' % (hours, minutes, seconds) 


# Generate a new temporary filename in the user's tmp directory
def tmp_name():
    tmp_root = tmp_dir
    if os.path.exists(tmp_root):
        pass
    else:        
        os.mkdir(tmp_root)
    return tmp_root + os.tmpnam().split('/')[-1] 


def die(msg=None):
 if msg is not None: 
    print >> sys.stderr, msg
    sys.exit(1)

    
# rough equivalent of the 'which' command to find external programs
# (current script path is tested first, then PATH envvar)
def which(program):
    def is_executable(fpath):
        return os.path.exists(fpath) and os.access(fpath, os.X_OK)
    fpath, fname = os.path.split(program)
    if fpath:
        if is_executable(program):
            return program
    else:
        progpath = os.path.join(bin_dir, program)
        if is_executable(progpath):
           return progpath
        for path in os.environ["PATH"].split(os.pathsep):
           progpath = os.path.join(path, program)
           if is_executable(progpath):
              return progpath
    return None


def die(msg=None):
 if msg is not None: 
    print >> sys.stderr, msg
    sys.exit(1)

    
def prog_path(program):
    progpath=which(program)
    if progpath == None:
        die("Error locating program: "+program)
    return progpath


def get_version():
   return "0.2.0"


def main(argv=None):
    warnings.filterwarnings("ignore", "tmpnam is a potential security risk")
    
    # Initialize default parameter values
    params = TopHatFusionParams()
    try:
        if argv is None:
            argv = sys.argv
            args = params.parse_options(argv)
            params.check()

        bwt_idx_prefix = args[0]

        print >> sys.stderr, "[%s] Beginning TopHat-Fusion post-processing run (v%s)" % (right_now(), get_version())
        print >> sys.stderr, "-----------------------------------------------" 
        
        start_time = datetime.now()
        prepare_output_dir()

        if not params.skip_fusion_kmer:
            map_fusion_kmer(bwt_idx_prefix, params)

        if not params.skip_filter_fusion:
            filter_fusion(params)

        if not params.skip_blast:
            do_blast(params)

        if not params.skip_read_dist:
            read_dist(params)

        if not params.skip_html:
            generate_html(params)
        
        global run_log
        run_log = open(logging_dir + "run.log", "w", 0)
        global run_cmd
        run_cmd = " ".join(argv)
        print >> run_log, run_cmd

        finish_time = datetime.now()
        duration = finish_time - start_time
        print >> sys.stderr,"-----------------------------------------------"
        print >> sys.stderr, "Run complete [%s elapsed]" %  formatTD(duration)
        
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, "\tfor detailed help see http://tophat-fusion.sourceforge.net/manual.html"
        return 2


if __name__ == "__main__":
    sys.exit(main())
