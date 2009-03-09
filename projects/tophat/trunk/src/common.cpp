/*
 *  common.cpp
 *  TopHat
 *
 *  Created by Cole Trapnell on 11/26/08.
 *  Copyright 2008 Cole Trapnell. All rights reserved.
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>

#include "getopt.h"
#include "common.h"

using namespace std;

int insert_len = 250;
int insert_len_std_dev = 20;
int max_mate_inner_dist = -1; 

uint32_t min_anchor_len = 5;
uint32_t min_intron_length = 70;
uint32_t max_intron_length = 20000;
uint32_t min_exon_length = 100; 
int island_extension = 25;

ReadFormat reads_format = FASTA;

extern void print_usage();

/**
 * Parse an int out of optarg and enforce that it be at least 'lower';
 * if it is less than 'lower', than output the given error message and
 * exit with an error and a usage message.
 */

int parseInt(int lower, const char *errmsg, void (*print_usage)()) {
    long l;
    char *endPtr= NULL;
    l = strtol(optarg, &endPtr, 10);
    if (endPtr != NULL) {
        if (l < lower) {
            cerr << errmsg << endl;
            print_usage();
            exit(1);
        }
        return (int32_t)l;
    }
    cerr << errmsg << endl;
    print_usage();
    exit(1);
    return -1;
}