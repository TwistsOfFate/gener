unsigned long long get_mask_from_hash_len(int hash_len)
{
    unsigned long long mask = 0;
    if (hash_len > 32) {
        mask = ~0;
    } else {
        for (auto i = 0; i < hash_len * 2; ++i) {
            mask = (mask << 1) + 1;
        }
    }
    return mask;
}
