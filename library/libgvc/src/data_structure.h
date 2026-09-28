#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>

uint16_t decode_16bit_id(unsigned char *payload, size_t bit_idx, 
                         uint8_t word_size, uint16_t mask);

void decode_ids(unsigned char *payload, uint16_t* ids, uint16_t num_ids);
