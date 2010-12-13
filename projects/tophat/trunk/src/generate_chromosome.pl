#!/usr/bin/env python

"""
generate_chromosome.pl
"""

import sys
import random
from operator import itemgetter

use_message = '''
Given a chromosome, the script introduces indels and snp.

Usage:
    generate_chromosome.pl input.fastq input.gtf
    , which gives input_var.fastq, input_var.gtf, and variant.list
'''

def driver(fasta_filename, gtf_filename):
    random.seed(0)
    
    def add_var_filename(filename):
        filename_array = filename.split('.')
        return filename_array[0] + "_var." + filename_array[1]

    gtf_out_filename = add_var_filename(gtf_filename)
    fasta_out_filename = add_var_filename(fasta_filename)

    fasta_out_file = open(fasta_out_filename, 'w')
    gtf_out_file = open(gtf_out_filename, 'w')

    fasta_file = open(fasta_filename, 'r')
    print >> fasta_out_file, fasta_file.readline()[:-1]

    def output_fasta(fasta, fata_width, fasta_out_file, end):
        temp_pos = 0
        while temp_pos + fasta_width <= len(fasta):
            print >> fasta_out_file, fasta[temp_pos:temp_pos+fasta_width]
            temp_pos += fasta_width

        if end and temp_pos < len(fasta):
            print >> fasta_out_file, fasta[temp_pos:]
            temp_pos = len(fasta)
            
        return temp_pos

    exon_list = []
    gtf_file = open(gtf_filename, 'r')
    for exon_line in gtf_file:
        exon = exon_line[:-1].split()
        exon[3], exon[4] = int(exon[3]), int(exon[4])
        exon_list.append(exon)
                         
    gtf_file.close

    variant_list_file = open("variant.list", 'w')

    exon_list2 = []
    for exon in sorted(exon_list, key=itemgetter(3)):
        if len(exon_list2) == 0 or exon_list2[-1][4] < exon[3]:
            exon_list2.append(exon)

    exon_list = exon_list2
    
    # 1-based offset
    curr_pos = 1
    fasta = fasta_file.readline()[:-1]
    fasta_width = len(fasta)
    offset_add = 0
    for i in range(len(exon_list)):
        exon = exon_list[i]
        exon_start, exon_end = exon[3] + offset_add, exon[4] + offset_add
        exon_length = exon_end - exon_start + 1

        for fasta_line in fasta_file:
            fasta += fasta_line[:-1]
            if curr_pos + len(fasta) - 1 < exon_start:
                temp_pos = output_fasta(fasta, fasta_width, fasta_out_file, False)
                curr_pos += temp_pos
                fasta = fasta[temp_pos:]
                continue

            if curr_pos + len(fasta) - 1 >= exon_end:
                break

        def get_random_index(p_list):
            p = random.random()
            for i in range(len(p_list)):
                p -= p_list[i]
                if p < 0.0:
                    return i

            return len(p_list) - 1

        def get_random_string(length, avoid_list = []):
            result = ""
            for i in range(length):
                rand_base = "ACGT"[random.randrange(4)]
                while rand_base in avoid_list:
                    rand_base = "ACGT"[random.randrange(4)]
                result += rand_base

            return result

        j = exon_start
        while j <= exon_end:
            if random.randrange(10000) == 0:
                insertion_pos = j - curr_pos
                insertion_len = get_random_index([0.6, 0.3, 0.1]) + 1
                
                fasta = fasta[:insertion_pos] + get_random_string(insertion_len)  + fasta[insertion_pos:]
                offset_add += insertion_len
                exon_end += insertion_len
                exon_length += insertion_len

                print >> variant_list_file, "ins\t%d\t%d" % (j - 1, j + insertion_len)
                
                j += (insertion_len + 1)
                next

            if random.randrange(10000) == 0:
                deletion_pos = j - curr_pos
                deletion_len = get_random_index([0.6, 0.3, 0.1]) + 1
                if j + deletion_len - 1 > exon_end:
                    next
                
                fasta = fasta[:deletion_pos] + fasta[deletion_pos + deletion_len:]
                offset_add -= deletion_len
                exon_end -= deletion_len
                exon_length -= deletion_len

                print >> variant_list_file, "del\t%d\t%d" % (j, deletion_len)
                next

            if random.randrange(1000) == 0:
                snp_pos = j - curr_pos
                fasta = fasta[:snp_pos] + get_random_string(1, fasta[snp_pos])[0] + fasta[snp_pos + 1:]

                print >> variant_list_file, "snp\t%d" % j
                
                j += 1
                next

            j += 1

        exon_list[i][3] = exon_start
        exon_list[i][4] = exon_end

    variant_list_file.close

    for fasta_line in fasta_file:
        fasta += fasta_line[:-1]
    output_fasta(fasta, fasta_width, fasta_out_file, True)
    
    fasta_file.close
    fasta_out_file.close
    gtf_out_file.close

    # test the correctness
    variant_list_file = open('variant.list', 'r')

    variant_list = []
    for variant in variant_list_file:
        variant_list.append(variant[:-1].split())
    variant_list_file.close
    
    fasta_file = open(fasta_filename, 'r')
    fasta_file.readline()
    
    fasta_out_file = open(fasta_out_filename, 'r')
    fasta_out_file.readline()

    fasta, fasta_out = "", ""
    fasta_pos = 1
    max_event_length = 30
    for variant in variant_list:
        while True:
            if fasta_pos + len(fasta_out) - 1 >= int(variant[1]) + max_event_length:
                break

            fasta_out += fasta_out_file.readline()[:-1]
            while len(fasta) < len(fasta_out):
                fasta += fasta_file.readline()[:-1]

            if fasta_pos + len(fasta_out) - 1 >= int(variant[1]):
                continue

            min_len = min(len(fasta), len(fasta_out))
            if fasta[:min_len] == fasta_out[:min_len]:
                fasta_pos += min_len
                fasta = fasta[min_len:]
                fasta_out = fasta_out[min_len:]
            else:
                print >> sys.stderr, fasta_pos, variant
                print >> sys.stderr, len(fasta), fasta
                print >> sys.stderr, len(fasta_out), fasta_out
                print >> sys.stderr, "error - this is not correct"
                exit(1)

        event_pos = int(variant[1]) - fasta_pos
                
        if variant[0] == "ins":
            if fasta[:event_pos+1] != fasta_out[:event_pos+1]:
                print >> sys.stderr, "error before insetion"
                exit(1)
                
            insertion_len = int(variant[2]) - int(variant[1]) - 1            
            fasta = fasta[event_pos + 1:]
            fasta_out = fasta_out[event_pos + 1 + insertion_len:]
            fasta_pos += (event_pos + 1 + insertion_len)
            
        elif variant[0] == "del":
            if fasta[:event_pos] != fasta_out[:event_pos]:
                print >> sys.stderr, "error before del"
                exit(1)

            deletion_len = int(variant[2])
            fasta = fasta[event_pos + deletion_len:]
            fasta_out = fasta_out[event_pos:]
            fasta_pos += event_pos

        elif variant[0] == "snp":
            if fasta[:event_pos] != fasta_out[:event_pos]:
                print >> sys.stderr, "error before snp"
                exit(1)

            fasta = fasta[event_pos+1:]
            fasta_out = fasta_out[event_pos+1:]
            fasta_pos += (event_pos + 1)
        
            
    fasta_out_file.close
    fasta_file.close

    print >> sys.stderr, "successfully generated!, take a look at variant.list"
    

if __name__ == "__main__":
    if len(sys.argv) == 3:
        fasta_filename = sys.argv[-2]
        gtf_filename = sys.argv[-1]
        driver(fasta_filename, gtf_filename)
    else:
        print use_message;
