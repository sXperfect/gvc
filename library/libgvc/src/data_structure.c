#include "data_structure.h"

static uint8_t bits_per_id(uint16_t num_ids)
{
    if (num_ids <= 1) {
        return 0;
    }

    uint16_t max_id = (uint16_t)(num_ids - 1);
    uint8_t bits = 0;
    while (max_id != 0) {
        ++bits;
        max_id >>= 1;
    }
    return bits;
}

uint16_t decode_16bit_id(unsigned char *payload, size_t bit_idx,
                         uint8_t word_size, uint16_t mask)
{
    uint16_t value = 0;
    for (uint8_t i = 0; i < word_size; ++i) {
        const size_t current_bit = bit_idx + i;
        const size_t byte_idx = current_bit / 8;
        const uint8_t bit_in_byte = (uint8_t)(7 - (current_bit % 8));
        value = (uint16_t)((value << 1) |
                          ((payload[byte_idx] >> bit_in_byte) & 1u));
    }
    return (uint16_t)(value & mask);
}

void decode_ids(unsigned char *payload, uint16_t *ids, uint16_t num_ids)
{
    if (num_ids == 0) {
        return;
    }

    const uint8_t word_size = bits_per_id(num_ids);
    if (word_size == 0) {
        ids[0] = 0;
        return;
    }

    const uint16_t mask = (uint16_t)((1u << word_size) - 1u);
    size_t bit_idx = 0;

    for (uint16_t i = 0; i < num_ids; ++i) {
        ids[i] = decode_16bit_id(payload, bit_idx, word_size, mask);
        bit_idx += word_size;
    }
}
