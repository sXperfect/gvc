#ifndef GVC_DATA_STRUCTURE_H
#define GVC_DATA_STRUCTURE_H

#include <stddef.h>
#include <stdint.h>

uint16_t decode_16bit_id(unsigned char *payload, size_t bit_idx,
                         uint8_t word_size, uint16_t mask);

/*
 * Returns 0 on success, -1 for an invalid output buffer, -2 for truncation,
 * -3 for an out-of-range id, -4 for duplicate ids, and -5 for trailing bytes.
 */
int decode_ids_checked(const unsigned char *payload, size_t payload_len,
                       uint16_t *ids, uint16_t num_ids);

/* Legacy ABI. New callers should prefer decode_ids_checked(). */
void decode_ids(unsigned char *payload, uint16_t *ids, uint16_t num_ids);

#endif
