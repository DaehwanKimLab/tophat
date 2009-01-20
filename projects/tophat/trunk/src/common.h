#ifndef COMMON_H
#define COMMON_H
/*
 *  common.h
 *  TopHat
 *
 *  Created by Cole Trapnell on 11/26/08.
 *  Copyright 2008 Cole Trapnell. All rights reserved.
 *
 */
#include <stdint.h>
#include <cassert>

extern int insert_len;
extern int insert_len_std_dev;
extern int max_mate_inner_dist; 

extern uint32_t min_anchor_len;
extern uint32_t min_intron_length;
extern uint32_t max_intron_length;

#endif
