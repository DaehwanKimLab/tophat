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

#include "common.h"

int insert_len = 250;
int insert_len_std_dev = 20;
int max_mate_inner_dist = -1; 

uint32_t min_anchor_len = 5;
uint32_t min_intron_length = 70;
uint32_t max_intron_length = 20000;

uint32_t min_exon_length = 50;