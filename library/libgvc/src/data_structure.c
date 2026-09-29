#include "data_structure.h"

static uint8_t bits_per_id(uint16_t num_ids)
{
    uint16_t value;
    uint8_t bits = 0;

    if (num_ids <= 1) {
        return 0;
    }

    value = (uint16_t)(num_ids - 1);
    while (value != 0) {
        bits++;
        value >>= 1;
    }
    return bits;
}

uint16_t decode_16bit_id(unsigned char *payload, size_t bit_idx,
                         uint8_t word_size, uint16_t mask)
{
    uint16_t value = 0;
    uint8_t bit;

    for (bit = 0; bit < word_size; bit++) {
        size_t current = bit_idx + bit;
        uint8_t source = payload[current / 8];
        uint8_t source_bit = (uint8_t)((source >> (7 - (current % 8))) & 1U);
        value = (uint16_t)((value << 1) | source_bit);
    }
    return (uint16_t)(value & mask);
}

int decode_ids_checked(const unsigned char *payload, size_t payload_len,
                       uint16_t *ids, uint16_t num_ids)
{
    uint8_t word_size;
    size_t required_bytes;
    size_t bit_idx = 0;
    uint16_t mask;
    uint16_t i;
    uint8_t seen[8192] = {0};

    if (num_ids == 0) {
        return payload_len == 0 ? 0 : -5;
    }
    if (ids == NULL) {
        return -1;
    }
    if (num_ids == 1) {
        if (payload_len != 0) {
            return -5;
        }
        ids[0] = 0;
        return 0;
    }

    word_size = bits_per_id(num_ids);
    required_bytes =
        ((size_t)num_ids * (size_t)word_size + 7U) / 8U;
    if (payload == NULL || payload_len < required_bytes) {
        return -2;
    }
    if (payload_len > required_bytes) {
        return -5;
    }

    mask = (uint16_t)((1U << word_size) - 1U);
    for (i = 0; i < num_ids; i++) {
        uint16_t value = decode_16bit_id(
            (unsigned char *)payload, bit_idx, word_size, mask
        );
        uint16_t seen_byte;
        uint8_t seen_mask;

        if (value >= num_ids) {
            return -3;
        }

        seen_byte = (uint16_t)(value >> 3);
        seen_mask = (uint8_t)(1U << (value & 7U));
        if ((seen[seen_byte] & seen_mask) != 0U) {
            return -4;
        }
        seen[seen_byte] |= seen_mask;

        ids[i] = value;
        bit_idx += word_size;
    }

    return 0;
}

void decode_ids(unsigned char *payload, uint16_t *ids, uint16_t num_ids)
{
    uint8_t word_size;
    size_t required_bytes;

    if (num_ids <= 1) {
        (void)decode_ids_checked(payload, 0, ids, num_ids);
        return;
    }

    word_size = bits_per_id(num_ids);
    required_bytes =
        ((size_t)num_ids * (size_t)word_size + 7U) / 8U;
    (void)decode_ids_checked(payload, required_bytes, ids, num_ids);
}


void decode_amax_vec(unsigned char *payload, uint8_t *amax,
                     uint16_t num_amax, uint8_t word_size)
{
    size_t bit_idx = 0;
    uint16_t mask = (uint16_t)((1U << word_size) - 1U);
    uint16_t i;

    for (i = 0; i < num_amax; i++) {
        amax[i] = (uint8_t)decode_16bit_id(
            payload, bit_idx, word_size, mask
        );
        bit_idx += word_size;
    }
}
