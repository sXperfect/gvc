#include "data_structure.h"

uint16_t decode_16bit_id(unsigned char *payload, size_t bit_idx, 
                         uint8_t word_size, uint16_t mask)
{
    size_t byte_idx = bit_idx / 8;
    size_t bit_in_byte_idx = bit_idx % 8;
    uint32_t val = payload[byte_idx] << 16 | payload[byte_idx+1] << 8 | payload[byte_idx+2];
    val >>= (24-word_size-bit_in_byte_idx);
    return val & mask;
}

void decode_ids(unsigned char *payload, uint16_t* ids, uint16_t num_ids)
{
    // TODO: word_size is maximum 15
    // TODO: Fix possible out-of-bound
    uint8_t word_size = (uint8_t) log2(num_ids);
    if (word_size > 8 || word_size < 8){

        size_t bit_idx = 0;
        uint16_t mask = (1<<word_size) - 1;
        for (int i = 0; i< num_ids; i++){
            ids[i] = decode_16bit_id(payload, bit_idx, word_size, mask);
            bit_idx += word_size;
        }
    } else if (word_size == 8){
        for (int i = 0; i< num_ids; i++){
            ids[i] = (uint16_t) payload[i];
        }
    }
}


void decode_amax_vec(unsigned char *payload, uint8_t* amax, 
                    uint16_t num_amax, uint8_t word_size)
{
    size_t bit_idx = 0;
    uint16_t mask = (1<<word_size) - 1;
    for (int i = 0; i< num_amax; i++){
        amax[i] = decode_16bit_id(payload, bit_idx, word_size, mask);
        bit_idx += word_size;
    }
}