#ifndef RW_TRIPLET_H
#define RW_TRIPLET_H
/**
 * @file rw_triplet.h
 * @brief The stages of the Random Window (RW) search selected by `rw_opt` (see RW_OPT_T3 etc. in util_io.h and
 *        doc/rw_triplet_degeneracies.tex): the weight-3 degeneracies T_3 of (H, L) (Stage 1), their contraction and
 *        the known generators G (Stage 2), the canonical classes of the minimum-weight codewords (Stage 5), and the
 *        expansion of the collected codewords by T_3 (rw_opt bit 32).
 *
 * A weight-3 degeneracy (a triplet) is a set {u, v, w} of three nonzero columns of H with h_u + h_v + h_w = 0 and
 * l_u + l_v + l_w = 0 (the columns of L): adding it to a codeword does not change its syndromes, e.g., a CNOT hook
 * triangle of a circuit DEM.  The triplets are found exactly by Lemma 1 of the notes: the smallest row of the union of
 * the three columns is in exactly two of them, which have this row as their smallest row, so that a lookup of the
 * Zobrist hash of h_u + h_v (and l_u + l_v) for each pair {u, v} of columns with the same smallest row finds every
 * triplet once.  They give the pair table (u, v) -> w used by the contraction of the candidate codewords (2 -> 1:
 * {u, v} -> {w}; 3 -> 0), by the canonical form of a codeword (2 <-> 2 exchanges across two triplets with a common
 * column), and by the expansion of the collected codewords.
 *
 * author: Leonid Pryadko <leonid.pryadko@ucr.edu>
 */

#include <stdint.h>
#include "util_io.h"

/** @brief The weight-3 degeneracies T_3 of (H, L) and their tables (Stage 1, see rw_t3_build()). */
typedef struct {
  int n;               /**< number of columns of H */
  long long ntrip;     /**< |T_3|, the number of triplets */
  long long probes;    /**< hash probes of the enumeration: the pairs of columns with the same smallest row */
  int complete;        /**< 1: T_3 enumerated (0: skipped, more probes than the budget; then ntrip = 0) */
  double ms;           /**< time of the enumeration in ms */
  int npure;           /**< number of pure-logical columns (h_u = 0, l_u != 0): codewords of weight 1 */
  int *pure;           /**< the pure-logical columns (increasing) */
  uint64_t *pkey;      /**< the pair table (open addressing): key (u << 32) | v of the pair u < v ... */
  int *pval;           /**< ... and w + 1 for the triplet {u, v, w} (0: an empty slot) */
  uint32_t pmask;      /**< number of slots - 1 */
  int *tp;             /**< the triplets of the column u: the pairs tq[2k], tq[2k+1], k = tp[u]..tp[u+1]-1 ... */
  int *tq;             /**< ... the two other columns of each triplet with u */
  int *gp;             /**< the triplet graph G_3: the distinct neighbors gi[gp[u]..gp[u+1]-1] of u */
  int *gi;
} rw_t3_t;

/**
 * @brief Enumerate the weight-3 degeneracies of (H, L) (Lemma 1 of doc/rw_triplet_degeneracies.tex).
 *
 * The columns are those of H (rows of HT, sorted) and the column masks of L (NULL for a classical code: then no
 * triplets are enumerated, as a weight-3 vector with zero syndrome is a codeword).  Zero columns of H are skipped (the
 * pure-logical ones are listed in `pure`).  If the enumeration needs more than `budget` hash probes, it is skipped
 * (`complete` = 0).  The tables are freed (and the counts kept) if there are no triplets.
 *
 * @param HT The transpose of H (n rows: the rows of H in each column of H, sorted).
 * @param m The number of rows of H.
 * @param maskL The column masks of L (NULL: classical code).
 * @param budget The largest number of hash probes.
 * @return The new structure (free with rw_t3_free()).
 */
rw_t3_t *rw_t3_build(const csr_t * const HT, const int m, const colmask_t * const maskL, const long long budget);

/** @brief Free the structure allocated by rw_t3_build() (NULL is ignored); returns NULL. */
rw_t3_t *rw_t3_free(rw_t3_t *t);

/** @brief The third column w of the triplet {u, v, w} in T_3 (-1 if none; u != v). */
static inline int rw_t3_pair(const rw_t3_t * const t, int u, int v) {
  if (u > v) {
    const int x = u;
    u = v;
    v = x;
  }
  const uint64_t key = ((uint64_t)(uint32_t)u << 32) | (uint32_t)v;
  uint32_t s = (uint32_t)((key * 0x9e3779b97f4a7c15ULL) >> 32) & t->pmask;
  while (t->pval[s]) {
    if (t->pkey[s] == key) return t->pval[s] - 1;
    s = (s + 1) & t->pmask;
  }
  return -1;
}

/**
 * @brief Greedy T_3 contraction of a candidate codeword (Stage 2): while two of its columns u, v form a triplet {u, v,
 *        w}, replace them by w (2 -> 1), or remove all three if w is in it too (3 -> 0).  The syndromes of H and L do
 *        not change.
 * @param t The triplets (with ntrip > 0).
 * @param supp The columns of the candidate (any order; changed in place, not sorted).
 * @param cnt Their number.
 * @return The new number of columns.
 */
int rw_t3_contract(const rw_t3_t * const t, int * const supp, int cnt);

/** @brief The known generators G (finG rows) for the RW stages (Stage 2b: shortening, the realized-generator test, and
 *         the half-overlap exchanges of Stage 5), see rw_gdb_build(). */
typedef struct {
  int nrows;            /**< number of nonzero rows of G kept */
  int n;                /**< number of columns */
  int *rp, *ri;         /**< the rows: sorted distinct columns ri[rp[g]..rp[g+1]-1] */
  int *cp, *ci;         /**< the rows of the column u: ci[cp[u]..cp[u+1]-1] */
  int wmax;             /**< the largest row weight w_G */
  int even;             /**< number of rows of even weight */
  unsigned char *wpres; /**< wpres[w] = 1 if a row has weight w (w = 0..wmax) */
  uint64_t *hkey;       /**< the Zobrist hashes of the rows (open addressing) ... */
  int *hval;            /**< ... and the row + 1 (0: empty slot) */
  uint32_t hmask;
} rw_gdb_t;

/**
 * @brief The tables of the known generators: the nonzero rows of G (repeated columns of a row count once), the
 *        column-to-row incidence lists, and the hash of the rows (the Zobrist hash of a set of columns is the XOR of
 *        `zcol` over them).
 * @param G The generator matrix (rows: elements of ker H with zero syndrome of L).
 * @param zcol The random 64-bit values of the columns (n of them), see rw_zobrist_cols().
 * @return The new structure (free with rw_gdb_free()), NULL if G has no nonzero rows.
 */
rw_gdb_t *rw_gdb_build(const csr_t * const G, const uint64_t * const zcol);

/** @brief Free the structure allocated by rw_gdb_build() (NULL is ignored); returns NULL. */
rw_gdb_t *rw_gdb_free(rw_gdb_t *g);

/** @brief The random 64-bit values of n columns (fixed seed) for the Zobrist hashes of sets of columns. */
uint64_t *rw_zobrist_cols(const int n);

/** @brief The Zobrist hash of the set of `cnt` columns `supp` (any order) with the weight and a `salt` (never 0). */
static inline uint64_t rw_set_key(const uint64_t * const zcol, const int * const supp, const int cnt,
                                  const uint64_t salt) {
  uint64_t h = (salt + (uint64_t)cnt) * 0x9e3779b97f4a7c15ULL;
  h ^= h >> 29;
  for (int i = 0; i < cnt; i++) h ^= zcol[supp[i]];
  return h ? h : 1;
}

/** @brief The workspace of a thread for the RW stages (see rw_tws_init()). */
typedef struct {
  int n;               /**< number of columns */
  uint32_t *mark;      /**< column marks of the generation `gen` */
  uint32_t gen;
  int *pos;            /**< the position of a marked column in the current support */
  int *ov;             /**< the overlaps of the current support with the rows of G (gdb) */
  uint32_t *rmark;     /**< rows of G with an overlap in the generation `rgen` (listed in `rows`) */
  uint32_t rgen;
  int *rows;
  int nrows_t;         /**< number of rows listed */
  int *ent;            /**< the canonical form: triples (hub, column, partner), see rw_canon() */
  size_t entcap;       /**< room for this many triples */
  int *tmp;            /**< a support of up to n columns */
  int *sub, *rest;     /**< for rw_cw_reduce() (n columns each) */
  cols_ws_t *irr;      /**< the irreducibility test (csr_cols_rank()) */
} rw_tws_t;

/**
 * @brief The workspace of a thread for the RW stages.
 * @param n The number of columns.
 * @param g The known generators (NULL: none).
 * @param m The number of rows of H (for the irreducibility test).
 * @return The new workspace (free with rw_tws_free()).
 */
rw_tws_t *rw_tws_init(const int n, const rw_gdb_t * const g, const int m);

/** @brief Free the workspace allocated by rw_tws_init() (NULL is ignored); returns NULL. */
rw_tws_t *rw_tws_free(rw_tws_t *ws);

/**
 * @brief Steepest-descent shortening of a candidate codeword by the rows of G (Theorem 5 of the notes): while a row g
 *        meets the support in more than half of its columns, the support is replaced by its sum with the row of the
 *        largest gain 2|c n g| - |g| (the smallest row of the largest gain).
 * @param g The known generators.
 * @param ws The workspace of the thread.
 * @param supp The columns of the candidate (any order, distinct; changed in place, not sorted; the weight only
 *        decreases).
 * @param cnt Their number.
 * @return The new number of columns.
 */
int rw_gdb_shorten(const rw_gdb_t * const g, rw_tws_t * const ws, int * const supp, int cnt);

/**
 * @brief Whether the set of `cnt` columns `supp` (any order) with the Zobrist hash `key` (rw_set_key() with salt 0) is
 *        a row of G (a realized generator, Theorem 5(c)).
 */
int rw_gdb_is_row(const rw_gdb_t * const g, rw_tws_t * const ws, const int * const supp, const int cnt,
                  const uint64_t key);

/**
 * @brief The canonical form of a codeword (Stage 5): with triplets (t->ntrip > 0), the greedy shared-hub 2 <-> 2
 *        canonicalization canon_T3 of Proposition 4 (the lexicographically smallest support c + T1 + T2 for two
 *        triplets T1 = {u_a, v_a, z}, T2 = {u_b, v_b, z} with u_a, u_b in c, v_a, v_b not in c, while it is smaller);
 *        otherwise, with even rows of G, canon_E of Proposition 6 (c + g for an even row g with |c n g| = |g|/2); and
 *        otherwise the support itself.  The weight and the syndromes do not change, so that two codewords of different
 *        logical classes never have the same canonical form.
 * @param t The triplets (NULL or ntrip = 0: none).
 * @param g The known generators (NULL: none).
 * @param ws The workspace of the thread.
 * @param supp The sorted columns of the codeword (replaced by its canonical form, sorted).
 * @param cnt Their number.
 */
void rw_canon(const rw_t3_t * const t, const rw_gdb_t * const g, rw_tws_t * const ws, int * const supp,
              const int cnt);

/**
 * @brief An irreducible part of a reducible codeword (irred=1): while a nonempty proper subset of the columns has zero
 *        syndrome of H, keep the lighter one of the subset and its complement which has a nonzero syndrome of L.
 * @param HT The transpose of H.
 * @param maskL The column masks of L (NULL: a classical code).
 * @param ws The workspace of the thread (with `irr`, `sub`, `rest`).
 * @param arr The sorted columns (replaced by the part, sorted).
 * @param w Their number.
 * @return The weight of the part.
 */
int rw_cw_reduce(const csr_t * const HT, const colmask_t * const maskL, rw_tws_t * const ws, int * const arr, int w);

/** @brief The table of the canonical classes of the minimum weight in the uniform RW steps (Stage 5; under the codeword
 *         mutex): the uniform-step incidence count of each class and the set of the raw codewords, keyed by 64-bit
 *         Zobrist hashes (see rw_set_key()). */
typedef struct {
  uint64_t *key;       /**< the keys of the classes (0: empty slot) */
  uint32_t *inc;       /**< their uniform-step incidence counts */
  size_t cap;          /**< number of slots (a power of 2) */
  size_t n;            /**< number of classes */
  uint64_t *rkey;      /**< the keys of the raw codewords (0: empty slot) */
  size_t rcap, rn;
  int full;            /**< 1: too many classes or raw codewords (RW_CANON_MAX): the table stopped growing */
} rw_canon_tab_t;

/** @brief The largest number of classes (and of raw codewords) in the table of rw_canon_tab_t. */
#define RW_CANON_MAX ((size_t)1 << 24)

/** @brief A new empty table of the canonical classes. */
rw_canon_tab_t *rw_canon_tab_new(void);

/** @brief Free the table (NULL is ignored); returns NULL. */
rw_canon_tab_t *rw_canon_tab_free(rw_canon_tab_t *tab);

/** @brief Remove all classes and raw codewords (for a new minimum weight). */
void rw_canon_tab_clear(rw_canon_tab_t * const tab);

/**
 * @brief The uniform-step incidence count of the class with the key (inserted with count 0 if new, `*is_new` = 1);
 *        NULL if the table is full.
 */
uint32_t *rw_canon_tab_class(rw_canon_tab_t * const tab, const uint64_t key, int * const is_new);

/** @brief Add the key of a raw codeword: 1 if new, 0 if known or the table is full. */
int rw_canon_tab_raw(rw_canon_tab_t * const tab, const uint64_t key);

/**
 * @brief Expand the collected codewords by T_3 (rw_opt bit 32, before outC is written): starting from the codewords of
 *        weight up to min_w + max(dW, 0) in the hash, a breadth-first search adds the equal-weight 2 <-> 2 exchanges
 *        across two triplets with a common column (Proposition 4) and the 1 -> 2 splits c - {z} + {u, v} for the
 *        triplets {u, v, z} with z in c and u, v not in c (weight + 1), up to the weight min_w + max(dW, 0), with the
 *        irreducibility test (irred=1) and maxC (see codeword_add_maybe()).
 * @param p The parameters (the codeword hash).
 * @param t The triplets (ntrip > 0).
 * @param HT The transpose of H.
 * @param ws The workspace (with `irr`).
 * @param limit The largest number of codewords to add.
 * @return The number of codewords added (negative: -1 - added if stopped by `limit`).
 */
long long rw_t3_expand(params_t * const p, const rw_t3_t * const t, const csr_t * const HT, rw_tws_t * const ws,
                       const long long limit);

#endif /* RW_TRIPLET_H */
