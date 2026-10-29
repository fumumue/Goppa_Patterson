int  pubkey_systematic(const unsigned short syndrome_matrix[N][K], uint64_t *T_out, size_t *fail_col);
void pubkey_systematic_print(const uint64_t *tt, size_t show_rows, size_t show_words);
int  pubkey_systematic_write(const char *path, const uint64_t *tt);