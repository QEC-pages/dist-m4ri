/**
 * @file rw_triplet.c
 * @brief The stages of the RW search selected by `rw_opt`: the weight-3 degeneracies T_3 (Stage 1), the contraction
 *        and the known generators (Stage 2), the canonical classes (Stage 5), and the expansion of the collected
 *        codewords by T_3 (rw_opt bit 32); see rw_triplet.h and doc/rw_triplet_degeneracies.tex.
 *
 * author: Leonid Pryadko <leonid.pryadko@ucr.edu>
 */

#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "rw_triplet.h"

/* The splitmix64 output function (a bijective 64-bit hash) and generator */
static inline uint64_t rw_mix64(uint64_t z) {
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  return z ^ (z >> 31);
}

static inline uint64_t rw_splitmix64(uint64_t * const s) { return rw_mix64(*s += 0x9e3779b97f4a7c15ULL); }

static double rw_wall(void) { /* seconds (monotonic clock) */
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

static void *rw_xmalloc(const size_t n) {
  void * const x = malloc(n ? n : 1);
  if (!x) ERROR("memory allocation");
  return x;
}

static void *rw_xcalloc(const size_t n, const size_t s) {
  void * const x = calloc(n ? n : 1, s ? s : 1);
  if (!x) ERROR("memory allocation");
  return x;
}

static int rw_cmp_int(const void *a, const void *b) {
  const int x = *(const int *)a, y = *(const int *)b;
  return (x > y) - (x < y);
}

static int rw_cmp_pair(const void *a, const void *b) { /* pairs of ints, lexicographic */
  const int * const x = a, * const y = b;
  if (x[0] != y[0]) return (x[0] > y[0]) - (x[0] < y[0]);
  return (x[1] > y[1]) - (x[1] < y[1]);
}

/* Insertion sort of a short array */
static inline void rw_isort(int * const a, const int n) {
  for (int i = 1; i < n; i++) {
    const int x = a[i];
    int j = i - 1;
    while (j >= 0 && a[j] > x) {
      a[j + 1] = a[j];
      j--;
    }
    a[j + 1] = x;
  }
}

static inline void rw_sort(int * const a, const int n) {
  if (n <= 32) rw_isort(a, n);
  else qsort(a, (size_t)n, sizeof(int), rw_cmp_int);
}

/* ---------------------------------------------------------------------------------------------------------------- *
 * Stage 1: the weight-3 degeneracies
 * ---------------------------------------------------------------------------------------------------------------- */

/* Is {a, b, c} a triplet: h_a + h_b = h_c (distinct sorted rows ci[cp[j]..cp[j+1]-1]) and l_a + l_b = l_c? */
static int t3_match(const int * const cp, const int * const ci, const colmask_t * const cm, const int a, const int b,
                    const int c) {
  const int kw = cm->kw;
  const word * const la = cm->mask + (size_t)a * kw, * const lb = cm->mask + (size_t)b * kw;
  const word * const lc = cm->mask + (size_t)c * kw;
  for (int w = 0; w < kw; w++)
    if ((la[w] ^ lb[w]) != lc[w]) return 0;
  int ia = cp[a], ib = cp[b], ic = cp[c];
  const int ea = cp[a + 1], eb = cp[b + 1], ec = cp[c + 1];
  while (ia < ea || ib < eb) {
    int r;
    if (ib >= eb) r = ci[ia++];
    else if (ia >= ea) r = ci[ib++];
    else if (ci[ia] < ci[ib]) r = ci[ia++];
    else if (ci[ib] < ci[ia]) r = ci[ib++];
    else { /* a common row cancels */
      ia++;
      ib++;
      continue;
    }
    if (ic >= ec || ci[ic] != r) return 0;
    ic++;
  }
  return ic == ec;
}

static void t3_pair_insert(rw_t3_t * const t, int u, int v, const int w) {
  if (u > v) {
    const int x = u;
    u = v;
    v = x;
  }
  const uint64_t key = ((uint64_t)(uint32_t)u << 32) | (uint32_t)v;
  uint32_t s = (uint32_t)((key * 0x9e3779b97f4a7c15ULL) >> 32) & t->pmask;
  while (t->pval[s]) {
    if (t->pkey[s] == key) { /* (repeated columns w: the smallest one, whatever the order of the rows of H) */
      if (w + 1 < t->pval[s]) t->pval[s] = w + 1;
      return;
    }
    s = (s + 1) & t->pmask;
  }
  t->pkey[s] = key;
  t->pval[s] = w + 1;
}

rw_t3_t *rw_t3_build(const csr_t * const HT, const int m, const colmask_t * const cm, const long long budget) {
  const double t0 = rw_wall();
  const int n = HT->rows;
  rw_t3_t * const t = rw_xcalloc(1, sizeof(rw_t3_t));
  t->n = n;
  t->complete = 1;
  if (!cm || n <= 0 || m <= 0) {
    t->ms = 1e3 * (rw_wall() - t0);
    return t;
  }
  const int kw = cm->kw;
  /* the distinct rows of each column (a repeated entry counts once, as in mzd_from_csr()) */
  int * const cp = rw_xmalloc(((size_t)n + 1) * sizeof(int));
  int * const ci = rw_xmalloc(((size_t)HT->p[n] + 1) * sizeof(int));
  int nz = 0;
  for (int j = 0; j < n; j++) {
    cp[j] = nz;
    for (int k = HT->p[j]; k < HT->p[j + 1]; k++)
      if (nz == cp[j] || ci[nz - 1] != HT->i[k]) ci[nz++] = HT->i[k];
  }
  cp[n] = nz;
  /* the Zobrist hash of each column (random keys of the rows of H and L), and the columns by their smallest row */
  uint64_t st = 0x517cc1b727220a95ULL;
  uint64_t * const zr = rw_xmalloc((size_t)m * sizeof(uint64_t));
  for (int r = 0; r < m; r++) zr[r] = rw_splitmix64(&st);
  uint64_t * const zl = rw_xmalloc((size_t)kw * 64 * sizeof(uint64_t));
  for (int r = 0; r < kw * 64; r++) zl[r] = rw_splitmix64(&st);
  uint64_t * const h = rw_xmalloc((size_t)n * sizeof(uint64_t));
  int * const bb = rw_xcalloc((size_t)m + 2, sizeof(int)); /* bucket sizes, then their starts */
  int npure = 0;
  for (int j = 0; j < n; j++) {
    const word * const lm = cm->mask + (size_t)j * kw;
    if (cp[j] == cp[j + 1]) { /* a zero column of H: a codeword of weight 1 if l_j != 0 */
      word any = 0;
      for (int w = 0; w < kw; w++) any |= lm[w];
      npure += (any != 0);
      h[j] = 0;
      continue;
    }
    uint64_t x = 0;
    for (int k = cp[j]; k < cp[j + 1]; k++) x ^= zr[ci[k]];
    for (int w = 0; w < kw; w++)
      for (word b = lm[w]; b; b &= b - 1) x ^= zl[w * 64 + __builtin_ctzll(b)];
    h[j] = x;
    bb[ci[cp[j]] + 1]++;
  }
  if (npure > 0) {
    t->pure = rw_xmalloc((size_t)npure * sizeof(int));
    for (int j = 0; j < n; j++) {
      if (cp[j] != cp[j + 1]) continue;
      const word * const lm = cm->mask + (size_t)j * kw;
      word any = 0;
      for (int w = 0; w < kw; w++) any |= lm[w];
      if (any) t->pure[t->npure++] = j;
    }
  }
  free(zr);
  free(zl);
  long long probes = 0;
  for (int r = 0; r < m; r++) probes += (long long)bb[r + 1] * (bb[r + 1] - 1) / 2;
  t->probes = probes;
  if (probes > budget) { /* skipped: too many pairs of columns with the same smallest row */
    t->complete = 0;
    free(cp);
    free(ci);
    free(h);
    free(bb);
    t->ms = 1e3 * (rw_wall() - t0);
    return t;
  }
  for (int r = 0; r <= m; r++) bb[r + 1] += bb[r];
  int * const bcol = rw_xmalloc((size_t)n * sizeof(int));
  int * const bpos = rw_xmalloc(((size_t)m + 1) * sizeof(int));
  memcpy(bpos, bb, ((size_t)m + 1) * sizeof(int));
  for (int j = 0; j < n; j++)
    if (cp[j] != cp[j + 1]) bcol[bpos[ci[cp[j]]]++] = j;
  free(bpos);
  /* the columns by their hash (chained) */
  uint32_t hsz = 16;
  while ((uint64_t)hsz < 2ULL * (uint64_t)n) hsz <<= 1;
  int * const head = rw_xmalloc((size_t)hsz * sizeof(int));
  int * const nxt = rw_xmalloc((size_t)n * sizeof(int));
  for (uint32_t s = 0; s < hsz; s++) head[s] = -1;
  for (int j = 0; j < n; j++) {
    if (cp[j] == cp[j + 1]) continue;
    const uint32_t s = (uint32_t)(h[j] >> 32) & (hsz - 1);
    nxt[j] = head[s];
    head[s] = j;
  }
  /* the triplets: for each pair {a, b} of columns with the same smallest row, the columns c with h_c = h_a + h_b */
  long long cap = 1024, nt = 0;
  int *tr = rw_xmalloc((size_t)cap * 3 * sizeof(int));
  for (int r = 0; r < m; r++) {
    for (int i1 = bb[r]; i1 < bb[r + 1]; i1++) {
      const int a = bcol[i1];
      for (int i2 = i1 + 1; i2 < bb[r + 1]; i2++) {
        const int b = bcol[i2];
        const uint64_t th = h[a] ^ h[b];
        if (!th) continue; /* equal columns a, b: c would be a zero column */
        for (int c = head[(uint32_t)(th >> 32) & (hsz - 1)]; c >= 0; c = nxt[c]) {
          if (h[c] != th || c == a || c == b || !t3_match(cp, ci, cm, a, b, c)) continue;
          if (nt == cap) {
            cap *= 2;
            tr = realloc(tr, (size_t)cap * 3 * sizeof(int));
            if (!tr) ERROR("memory allocation");
          }
          tr[3 * nt] = a;
          tr[3 * nt + 1] = b;
          tr[3 * nt + 2] = c;
          nt++;
        }
      }
    }
  }
  free(head);
  free(nxt);
  free(bcol);
  free(bb);
  free(h);
  free(cp);
  free(ci);
  if (nt == 0 || 6 * nt > ((long long)1 << 30)) { /* none (or far too many for the tables) */
    if (nt > 0) t->complete = 0;
    free(tr);
    t->ms = 1e3 * (rw_wall() - t0);
    return t;
  }
  t->ntrip = nt;
  /* the pair table (at most half full) */
  uint32_t ps = 16;
  while ((uint64_t)ps < 6ULL * (uint64_t)nt + 32) ps <<= 1;
  t->pmask = ps - 1;
  t->pkey = rw_xcalloc(ps, sizeof(uint64_t));
  t->pval = rw_xcalloc(ps, sizeof(int));
  for (long long k = 0; k < nt; k++) {
    const int a = tr[3 * k], b = tr[3 * k + 1], c = tr[3 * k + 2];
    t3_pair_insert(t, a, b, c);
    t3_pair_insert(t, a, c, b);
    t3_pair_insert(t, b, c, a);
  }
  /* the triplets of each column */
  t->tp = rw_xcalloc((size_t)n + 1, sizeof(int));
  for (long long k = 0; k < 3 * nt; k++) t->tp[tr[k] + 1]++;
  for (int u = 0; u < n; u++) t->tp[u + 1] += t->tp[u];
  t->tq = rw_xmalloc((size_t)nt * 6 * sizeof(int));
  int * const cur = rw_xmalloc(((size_t)n + 1) * sizeof(int));
  memcpy(cur, t->tp, ((size_t)n + 1) * sizeof(int));
  for (long long k = 0; k < nt; k++) {
    const int a = tr[3 * k], b = tr[3 * k + 1], c = tr[3 * k + 2];
    t->tq[2 * cur[a]] = b;
    t->tq[2 * cur[a] + 1] = c;
    cur[a]++;
    t->tq[2 * cur[b]] = a;
    t->tq[2 * cur[b] + 1] = c;
    cur[b]++;
    t->tq[2 * cur[c]] = a;
    t->tq[2 * cur[c] + 1] = b;
    cur[c]++;
  }
  free(tr);
  /* each list of triplets in a canonical order (the pairs of other columns, smaller first, sorted), so that the
   * stages do not depend on the order of the rows of H, which gives the order of the enumeration */
  for (int u = 0; u < n; u++) {
    for (int k = t->tp[u]; k < t->tp[u + 1]; k++) {
      if (t->tq[2 * k] > t->tq[2 * k + 1]) {
        const int x = t->tq[2 * k];
        t->tq[2 * k] = t->tq[2 * k + 1];
        t->tq[2 * k + 1] = x;
      }
    }
    if (t->tp[u + 1] - t->tp[u] > 1)
      qsort(t->tq + 2 * t->tp[u], (size_t)(t->tp[u + 1] - t->tp[u]), 2 * sizeof(int), rw_cmp_pair);
  }
  /* the triplet graph G_3: the distinct neighbors of each column */
  t->gp = rw_xmalloc(((size_t)n + 1) * sizeof(int));
  t->gi = rw_xmalloc((size_t)nt * 6 * sizeof(int));
  int ng = 0;
  for (int u = 0; u < n; u++) {
    t->gp[u] = ng;
    const int k0 = 2 * t->tp[u], k1 = 2 * t->tp[u + 1];
    memcpy(t->gi + ng, t->tq + k0, (size_t)(k1 - k0) * sizeof(int));
    rw_sort(t->gi + ng, k1 - k0);
    int q = 0;
    for (int k = 0; k < k1 - k0; k++)
      if (q == 0 || t->gi[ng + q - 1] != t->gi[ng + k]) t->gi[ng + q++] = t->gi[ng + k];
    ng += q;
  }
  t->gp[n] = ng;
  free(cur);
  t->ms = 1e3 * (rw_wall() - t0);
  return t;
}

rw_t3_t *rw_t3_free(rw_t3_t *t) {
  if (!t) return NULL;
  free(t->pure);
  free(t->pkey);
  free(t->pval);
  free(t->tp);
  free(t->tq);
  free(t->gp);
  free(t->gi);
  free(t);
  return NULL;
}

int rw_t3_contract(const rw_t3_t * const t, int * const supp, int cnt) {
  if (!t || t->ntrip <= 0) return cnt;
  int changed = 1;
  while (changed && cnt >= 2) {
    changed = 0;
    for (int i = 0; i < cnt && !changed; i++) {
      const int u = supp[i];
      for (int j = i + 1; j < cnt; j++) {
        const int w = rw_t3_pair(t, u, supp[j]);
        if (w < 0) continue;
        supp[j] = supp[--cnt]; /* (j > i: u stays at i) */
        int k = -1;
        for (int q = 0; q < cnt; q++) {
          if (q != i && supp[q] == w) {
            k = q;
            break;
          }
        }
        if (k >= 0) { /* 3 -> 0: remove u (at i) and w (at k), the larger position first */
          const int hi = (i > k) ? i : k, lo = (i > k) ? k : i;
          supp[hi] = supp[--cnt];
          supp[lo] = supp[--cnt];
        } else { /* 2 -> 1 */
          supp[i] = w;
        }
        changed = 1;
        break;
      }
    }
  }
  return cnt;
}

/* ---------------------------------------------------------------------------------------------------------------- *
 * Stage 2b: the known generators
 * ---------------------------------------------------------------------------------------------------------------- */

uint64_t *rw_zobrist_cols(const int n) {
  uint64_t * const z = rw_xmalloc(((size_t)n + 1) * sizeof(uint64_t));
  uint64_t st = 0x71e2b3c4d5e6f709ULL; /* fixed: the results do not depend on the hash */
  for (int j = 0; j < n; j++) z[j] = rw_splitmix64(&st);
  return z;
}

static inline uint32_t rw_hash_slot(const uint64_t key, const uint32_t mask) {
  return (uint32_t)(key ^ (key >> 32)) & mask;
}

rw_gdb_t *rw_gdb_build(const csr_t * const G, const uint64_t * const zcol) {
  if (!G || G->rows <= 0 || G->nz != -1) return NULL;
  const int n = G->cols;
  rw_gdb_t * const g = rw_xcalloc(1, sizeof(rw_gdb_t));
  g->n = n;
  g->rp = rw_xmalloc(((size_t)G->rows + 1) * sizeof(int));
  g->ri = rw_xmalloc(((size_t)G->p[G->rows] + 1) * sizeof(int));
  int nr = 0, nz = 0;
  for (int r = 0; r < G->rows; r++) {
    const int b = nz;
    for (int k = G->p[r]; k < G->p[r + 1]; k++) g->ri[nz++] = G->i[k];
    rw_sort(g->ri + b, nz - b);
    int q = b;
    for (int k = b; k < nz; k++) /* (a repeated column counts once, as in mzd_from_csr()) */
      if (q == b || g->ri[q - 1] != g->ri[k]) g->ri[q++] = g->ri[k];
    nz = q;
    if (nz == b) continue; /* a zero row */
    g->rp[nr++] = b;
    const int w = nz - b;
    if (w > g->wmax) g->wmax = w;
    g->even += !(w & 1);
  }
  g->rp[nr] = nz;
  g->nrows = nr;
  if (nr == 0) return rw_gdb_free(g);
  g->wpres = rw_xcalloc((size_t)g->wmax + 1, 1);
  g->cp = rw_xcalloc((size_t)n + 1, sizeof(int));
  for (int k = 0; k < nz; k++) g->cp[g->ri[k] + 1]++;
  for (int u = 0; u < n; u++) g->cp[u + 1] += g->cp[u];
  g->ci = rw_xmalloc(((size_t)nz + 1) * sizeof(int));
  int * const cur = rw_xmalloc(((size_t)n + 1) * sizeof(int));
  memcpy(cur, g->cp, ((size_t)n + 1) * sizeof(int));
  uint32_t hs = 16;
  while ((uint64_t)hs < 2ULL * (uint64_t)nr) hs <<= 1;
  g->hmask = hs - 1;
  g->hkey = rw_xcalloc(hs, sizeof(uint64_t));
  g->hval = rw_xcalloc(hs, sizeof(int));
  for (int r = 0; r < nr; r++) {
    const int w = g->rp[r + 1] - g->rp[r];
    g->wpres[w] = 1;
    for (int k = g->rp[r]; k < g->rp[r + 1]; k++) g->ci[cur[g->ri[k]]++] = r;
    const uint64_t key = rw_set_key(zcol, g->ri + g->rp[r], w, 0);
    uint32_t s = rw_hash_slot(key, g->hmask);
    while (g->hval[s]) s = (s + 1) & g->hmask;
    g->hkey[s] = key;
    g->hval[s] = r + 1;
  }
  free(cur);
  return g;
}

rw_gdb_t *rw_gdb_free(rw_gdb_t *g) {
  if (!g) return NULL;
  free(g->rp);
  free(g->ri);
  free(g->cp);
  free(g->ci);
  free(g->wpres);
  free(g->hkey);
  free(g->hval);
  free(g);
  return NULL;
}

rw_tws_t *rw_tws_init(const int n, const rw_gdb_t * const g, const int m) {
  rw_tws_t * const ws = rw_xcalloc(1, sizeof(rw_tws_t));
  ws->n = n;
  ws->mark = rw_xcalloc((size_t)n + 1, sizeof(uint32_t));
  ws->pos = rw_xmalloc(((size_t)n + 1) * sizeof(int));
  ws->tmp = rw_xmalloc(((size_t)n + 2) * sizeof(int));
  ws->sub = rw_xmalloc(((size_t)n + 2) * sizeof(int));
  ws->rest = rw_xmalloc(((size_t)n + 2) * sizeof(int));
  if (g) {
    ws->ov = rw_xcalloc((size_t)g->nrows + 1, sizeof(int));
    ws->rmark = rw_xcalloc((size_t)g->nrows + 1, sizeof(uint32_t));
    ws->rows = rw_xmalloc(((size_t)g->nrows + 1) * sizeof(int));
  }
  ws->irr = cols_ws_init(m);
  return ws;
}

rw_tws_t *rw_tws_free(rw_tws_t *ws) {
  if (!ws) return NULL;
  free(ws->mark);
  free(ws->pos);
  free(ws->tmp);
  free(ws->sub);
  free(ws->rest);
  free(ws->ov);
  free(ws->rmark);
  free(ws->rows);
  free(ws->ent);
  cols_ws_free(ws->irr);
  free(ws);
  return NULL;
}

/* A new generation of the column marks (cleared when the counter wraps around; 0 is never a generation) */
static inline uint32_t tws_gen(rw_tws_t * const ws) {
  if (++ws->gen == 0) {
    memset(ws->mark, 0, ((size_t)ws->n + 1) * sizeof(uint32_t));
    ws->gen = 1;
  }
  return ws->gen;
}

/* A new generation of the row marks of G (nrows rows) */
static inline uint32_t tws_rgen(rw_tws_t * const ws, const int nrows) {
  if (++ws->rgen == 0) {
    memset(ws->rmark, 0, ((size_t)nrows + 1) * sizeof(uint32_t));
    ws->rgen = 1;
  }
  ws->nrows_t = 0;
  return ws->rgen;
}

/* Add `d` to the overlaps of the rows of G with the column u (rows listed in ws->rows) */
static inline void tws_overlap(const rw_gdb_t * const g, rw_tws_t * const ws, const int u, const int d,
                               const uint32_t rg) {
  for (int k = g->cp[u]; k < g->cp[u + 1]; k++) {
    const int r = g->ci[k];
    if (ws->rmark[r] != rg) {
      ws->rmark[r] = rg;
      ws->ov[r] = 0;
      ws->rows[ws->nrows_t++] = r;
    }
    ws->ov[r] += d;
  }
}

int rw_gdb_shorten(const rw_gdb_t * const g, rw_tws_t * const ws, int * const supp, int cnt) {
  if (!g || cnt <= 0) return cnt;
  const uint32_t gen = tws_gen(ws);
  const uint32_t rg = tws_rgen(ws, g->nrows);
  for (int i = 0; i < cnt; i++) {
    const int u = supp[i];
    ws->mark[u] = gen;
    ws->pos[u] = i;
    tws_overlap(g, ws, u, 1, rg);
  }
  for (;;) {
    int best = -1, bgain = 0;
    for (int k = 0; k < ws->nrows_t; k++) {
      const int r = ws->rows[k];
      const int gain = 2 * ws->ov[r] - (g->rp[r + 1] - g->rp[r]);
      if (gain > 0 && (gain > bgain || (gain == bgain && r < best))) {
        best = r;
        bgain = gain;
      }
    }
    if (best < 0) break;
    for (int k = g->rp[best]; k < g->rp[best + 1]; k++) {
      const int u = g->ri[k];
      if (ws->mark[u] == gen) { /* remove u */
        const int i = ws->pos[u], last = supp[--cnt];
        supp[i] = last;
        ws->pos[last] = i;
        ws->mark[u] = 0;
        tws_overlap(g, ws, u, -1, rg);
      } else { /* add u */
        supp[cnt] = u;
        ws->pos[u] = cnt++;
        ws->mark[u] = gen;
        tws_overlap(g, ws, u, 1, rg);
      }
    }
  }
  return cnt;
}

int rw_gdb_is_row(const rw_gdb_t * const g, rw_tws_t * const ws, const int * const supp, const int cnt,
                  const uint64_t key) {
  if (!g || cnt <= 0 || cnt > g->wmax || !g->wpres[cnt]) return 0;
  int marked = 0;
  uint32_t gen = 0;
  for (uint32_t s = rw_hash_slot(key, g->hmask); g->hval[s]; s = (s + 1) & g->hmask) {
    if (g->hkey[s] != key) continue;
    const int r = g->hval[s] - 1;
    if (g->rp[r + 1] - g->rp[r] != cnt) continue;
    if (!marked) {
      gen = tws_gen(ws);
      for (int i = 0; i < cnt; i++) ws->mark[supp[i]] = gen;
      marked = 1;
    }
    int k = g->rp[r];
    while (k < g->rp[r + 1] && ws->mark[g->ri[k]] == gen) k++;
    if (k == g->rp[r + 1]) return 1;
  }
  return 0;
}

/* ---------------------------------------------------------------------------------------------------------------- *
 * Stage 5: the canonical form
 * ---------------------------------------------------------------------------------------------------------------- */

static int rw_cmp_ent(const void *a, const void *b) { /* triples (hub, column, partner) by hub, then column */
  const int * const x = a, * const y = b;
  if (x[0] != y[0]) return (x[0] > y[0]) - (x[0] < y[0]);
  if (x[1] != y[1]) return (x[1] > y[1]) - (x[1] < y[1]);
  return (x[2] > y[2]) - (x[2] < y[2]);
}

/* The 2 <-> 2 exchanges of the support c (marked with `gen`): the triples (z, u, v) for each column u of c and each
 * triplet {u, v, z} with v not in c (both choices of the hub z), sorted by z; two triples (z, u_a, v_a), (z, u_b, v_b)
 * with u_a != u_b and v_a != v_b give the exchange c + {u_a, v_a, u_b, v_b} of the same weight.  Returns the number of
 * triples (in ws->ent). */
static size_t t3_exchanges(const rw_t3_t * const t, rw_tws_t * const ws, const int * const c, const int cnt,
                           const uint32_t gen) {
  size_t need = 0;
  for (int a = 0; a < cnt; a++) need += 2 * (size_t)(t->tp[c[a] + 1] - t->tp[c[a]]);
  if (need > ws->entcap) {
    free(ws->ent);
    ws->entcap = 2 * need + 64;
    ws->ent = rw_xmalloc(ws->entcap * 3 * sizeof(int));
  }
  size_t ne = 0;
  int * const e = ws->ent;
  for (int a = 0; a < cnt; a++) {
    const int ua = c[a];
    for (int k = t->tp[ua]; k < t->tp[ua + 1]; k++) {
      const int x = t->tq[2 * k], y = t->tq[2 * k + 1];
      if (ws->mark[y] != gen) {
        e[3 * ne] = x;
        e[3 * ne + 1] = ua;
        e[3 * ne + 2] = y;
        ne++;
      }
      if (ws->mark[x] != gen) {
        e[3 * ne] = y;
        e[3 * ne + 1] = ua;
        e[3 * ne + 2] = x;
        ne++;
      }
    }
  }
  if (ne > 1) qsort(e, ne, 3 * sizeof(int), rw_cmp_ent);
  return ne;
}

/* Is the support c + Q1 lexicographically smaller than c + Q2 (Q1, Q2: sets of four columns; `inc`: the membership
 * marks of c with `gen`)?  The smallest column m in exactly one of Q1, Q2 decides: c + Q1 contains m iff m is added by
 * Q1 (not in c) or removed by Q2 (in c). */
static int q4_less(const rw_tws_t * const ws, const uint32_t gen, const int * const q1, const int * const q2) {
  int m = INT_MAX, in1 = 0;
  for (int i = 0; i < 4; i++) {
    int k = 0;
    while (k < 4 && q2[k] != q1[i]) k++;
    if (k == 4 && q1[i] < m) {
      m = q1[i];
      in1 = 1;
    }
    k = 0;
    while (k < 4 && q1[k] != q2[i]) k++;
    if (k == 4 && q2[i] < m) {
      m = q2[i];
      in1 = 0;
    }
  }
  if (m == INT_MAX) return 0;
  const int in_c = (ws->mark[m] == gen);
  return in1 ? !in_c : in_c;
}

/* canon_T3 (Proposition 4): the greedy lexicographic 2 <-> 2 canonicalization (c sorted, kept sorted) */
static void canon_t3(const rw_t3_t * const t, rw_tws_t * const ws, int * const c, const int cnt) {
  const int max_iter = 4 * cnt + 16;
  for (int it = 0; it < max_iter; it++) {
    const uint32_t gen = tws_gen(ws);
    for (int a = 0; a < cnt; a++) ws->mark[c[a]] = gen;
    const size_t ne = t3_exchanges(t, ws, c, cnt, gen);
    const int * const e = ws->ent;
    int best[4] = {0, 0, 0, 0}, found = 0;
    for (size_t g0 = 0; g0 < ne;) {
      size_t g1 = g0 + 1;
      while (g1 < ne && e[3 * g1] == e[3 * g0]) g1++;
      for (size_t i = g0; i < g1; i++) {
        for (size_t j = i + 1; j < g1; j++) {
          const int ua = e[3 * i + 1], va = e[3 * i + 2], ub = e[3 * j + 1], vb = e[3 * j + 2];
          if (ua == ub || va == vb) continue;
          /* smaller than c iff the smallest of the four columns is added */
          const int mu = (ua < ub) ? ua : ub, mv = (va < vb) ? va : vb;
          if (mv > mu) continue;
          const int q[4] = {ua, ub, va, vb};
          if (!found || q4_less(ws, gen, q, best)) {
            memcpy(best, q, sizeof(best));
            found = 1;
          }
        }
      }
      g0 = g1;
    }
    if (!found) break;
    for (int a = 0; a < cnt; a++) { /* replace u_a, u_b by v_a, v_b */
      if (c[a] == best[0]) c[a] = best[2];
      else if (c[a] == best[1]) c[a] = best[3];
    }
    rw_sort(c, cnt);
  }
}

/* Is the support c + r1 lexicographically smaller than c + r2 for two rows r1 != r2 of G (sorted)? */
static int row_less(const rw_gdb_t * const g, const rw_tws_t * const ws, const uint32_t gen, const int r1,
                    const int r2) {
  int i = g->rp[r1], j = g->rp[r2];
  const int ie = g->rp[r1 + 1], je = g->rp[r2 + 1];
  while (i < ie || j < je) {
    if (j >= je || (i < ie && g->ri[i] < g->ri[j])) {
      const int in_c = (ws->mark[g->ri[i]] == gen); /* in r1 only */
      return !in_c;
    }
    if (i >= ie || g->ri[j] < g->ri[i]) {
      const int in_c = (ws->mark[g->ri[j]] == gen); /* in r2 only */
      return in_c;
    }
    i++;
    j++;
  }
  return 0; /* equal rows */
}

/* canon_E (Proposition 6): the greedy lexicographic half-overlap exchanges with the even rows of G (c sorted) */
static void canon_e(const rw_gdb_t * const g, rw_tws_t * const ws, int * const c, const int cnt) {
  const int max_iter = 4 * cnt + 16;
  for (int it = 0; it < max_iter; it++) {
    const uint32_t gen = tws_gen(ws);
    const uint32_t rg = tws_rgen(ws, g->nrows);
    for (int a = 0; a < cnt; a++) {
      ws->mark[c[a]] = gen;
      tws_overlap(g, ws, c[a], 1, rg);
    }
    int best = -1;
    for (int k = 0; k < ws->nrows_t; k++) {
      const int r = ws->rows[k], w = g->rp[r + 1] - g->rp[r];
      if ((w & 1) || 2 * ws->ov[r] != w) continue;
      if (ws->mark[g->ri[g->rp[r]]] == gen) continue; /* smaller iff the smallest column of the row is added */
      if (best < 0 || row_less(g, ws, gen, r, best)) best = r;
    }
    if (best < 0) break;
    int q = 0;
    for (int a = 0; a < cnt; a++) { /* the columns of c not in the row (the row is sorted) */
      const int u = c[a];
      int lo = g->rp[best], hi = g->rp[best + 1];
      while (lo < hi) {
        const int mid = (lo + hi) >> 1;
        if (g->ri[mid] < u) lo = mid + 1;
        else hi = mid;
      }
      if (lo < g->rp[best + 1] && g->ri[lo] == u) continue;
      ws->tmp[q++] = u;
    }
    for (int k = g->rp[best]; k < g->rp[best + 1]; k++)
      if (ws->mark[g->ri[k]] != gen) ws->tmp[q++] = g->ri[k];
    if (q != cnt) break; /* (not reached: a half-overlap exchange keeps the weight) */
    memcpy(c, ws->tmp, (size_t)cnt * sizeof(int));
    rw_sort(c, cnt);
  }
}

void rw_canon(const rw_t3_t * const t, const rw_gdb_t * const g, rw_tws_t * const ws, int * const supp,
              const int cnt) {
  if (cnt < 2) return;
  if (t && t->ntrip > 0) canon_t3(t, ws, supp, cnt);
  else if (g && g->even > 0) canon_e(g, ws, supp, cnt);
}

int rw_cw_reduce(const csr_t * const HT, const colmask_t * const maskL, rw_tws_t * const ws, int * const arr, int w) {
  int k;
  while ((k = csr_cols_dependency(HT, w, arr, ws->irr, ws->sub)) > 0) {
    int nr = 0;
    for (int i = 0, j = 0; i < w; i++) {
      if (j < k && ws->sub[j] == arr[i]) j++;
      else ws->rest[nr++] = arr[i];
    }
    const int nt_sub = !maskL || colmask_syndrome_non_zero(maskL, k, ws->sub);
    const int nt_rest = !maskL || colmask_syndrome_non_zero(maskL, nr, ws->rest);
    const int keep_sub = nt_sub && (!nt_rest || k <= nr);
    memcpy(arr, keep_sub ? ws->sub : ws->rest, (size_t)(keep_sub ? k : nr) * sizeof(int));
    w = keep_sub ? k : nr;
  }
  return w;
}

/* ---------------------------------------------------------------------------------------------------------------- *
 * Stage 5: the table of the canonical classes
 * ---------------------------------------------------------------------------------------------------------------- */

rw_canon_tab_t *rw_canon_tab_new(void) {
  rw_canon_tab_t * const tab = rw_xcalloc(1, sizeof(rw_canon_tab_t));
  tab->cap = 1024;
  tab->key = rw_xcalloc(tab->cap, sizeof(uint64_t));
  tab->inc = rw_xcalloc(tab->cap, sizeof(uint32_t));
  tab->rcap = 1024;
  tab->rkey = rw_xcalloc(tab->rcap, sizeof(uint64_t));
  return tab;
}

rw_canon_tab_t *rw_canon_tab_free(rw_canon_tab_t *tab) {
  if (!tab) return NULL;
  free(tab->key);
  free(tab->inc);
  free(tab->rkey);
  free(tab);
  return NULL;
}

void rw_canon_tab_clear(rw_canon_tab_t * const tab) {
  memset(tab->key, 0, tab->cap * sizeof(uint64_t));
  memset(tab->inc, 0, tab->cap * sizeof(uint32_t));
  memset(tab->rkey, 0, tab->rcap * sizeof(uint64_t));
  tab->n = tab->rn = 0;
  tab->full = 0;
}

static inline size_t tab_slot(const uint64_t key, const size_t mask) { return (size_t)(key ^ (key >> 29)) & mask; }

uint32_t *rw_canon_tab_class(rw_canon_tab_t * const tab, const uint64_t key, int * const is_new) {
  *is_new = 0;
  size_t s = tab_slot(key, tab->cap - 1);
  while (tab->key[s]) {
    if (tab->key[s] == key) return &tab->inc[s];
    s = (s + 1) & (tab->cap - 1);
  }
  if (tab->n >= RW_CANON_MAX) {
    tab->full = 1;
    return NULL;
  }
  if (2 * (tab->n + 1) > tab->cap) { /* twice as many slots (at most half full) */
    const size_t cap = 2 * tab->cap;
    uint64_t * const k2 = rw_xcalloc(cap, sizeof(uint64_t));
    uint32_t * const i2 = rw_xcalloc(cap, sizeof(uint32_t));
    for (size_t q = 0; q < tab->cap; q++) {
      if (!tab->key[q]) continue;
      size_t t = tab_slot(tab->key[q], cap - 1);
      while (k2[t]) t = (t + 1) & (cap - 1);
      k2[t] = tab->key[q];
      i2[t] = tab->inc[q];
    }
    free(tab->key);
    free(tab->inc);
    tab->key = k2;
    tab->inc = i2;
    tab->cap = cap;
    s = tab_slot(key, cap - 1);
    while (tab->key[s]) s = (s + 1) & (cap - 1);
  }
  tab->key[s] = key;
  tab->inc[s] = 0;
  tab->n++;
  *is_new = 1;
  return &tab->inc[s];
}

int rw_canon_tab_raw(rw_canon_tab_t * const tab, const uint64_t key) {
  size_t s = tab_slot(key, tab->rcap - 1);
  while (tab->rkey[s]) {
    if (tab->rkey[s] == key) return 0;
    s = (s + 1) & (tab->rcap - 1);
  }
  if (tab->rn >= RW_CANON_MAX) {
    tab->full = 1;
    return 0;
  }
  if (2 * (tab->rn + 1) > tab->rcap) {
    const size_t cap = 2 * tab->rcap;
    uint64_t * const k2 = rw_xcalloc(cap, sizeof(uint64_t));
    for (size_t q = 0; q < tab->rcap; q++) {
      if (!tab->rkey[q]) continue;
      size_t t = tab_slot(tab->rkey[q], cap - 1);
      while (k2[t]) t = (t + 1) & (cap - 1);
      k2[t] = tab->rkey[q];
    }
    free(tab->rkey);
    tab->rkey = k2;
    tab->rcap = cap;
    s = tab_slot(key, cap - 1);
    while (tab->rkey[s]) s = (s + 1) & (cap - 1);
  }
  tab->rkey[s] = key;
  tab->rn++;
  return 1;
}

/* ---------------------------------------------------------------------------------------------------------------- *
 * rw_opt bit 32: the expansion of the collected codewords by T_3
 * ---------------------------------------------------------------------------------------------------------------- */

/* The queue of the breadth-first search: the supports of the codewords to expand */
typedef struct {
  int *pool;
  size_t len, cap;   /* ints used and allocated in pool */
  size_t head;       /* the next codeword: its weight at pool[head], then its columns */
} rw_bfs_t;

static void bfs_push(rw_bfs_t * const q, const int * const c, const int w) {
  if (q->len + (size_t)w + 1 > q->cap) {
    q->cap = 2 * (q->len + (size_t)w + 1) + 1024;
    q->pool = realloc(q->pool, q->cap * sizeof(int));
    if (!q->pool) ERROR("memory allocation");
  }
  q->pool[q->len++] = w;
  memcpy(q->pool + q->len, c, (size_t)w * sizeof(int));
  q->len += (size_t)w;
}

/* Add the codeword c (w columns, any order; sorted in place) to the hash if new and irreducible (irred=1); returns 1 if
 * added (and queued), -1 if maxC codewords are collected, 0 otherwise */
static int bfs_try(params_t * const p, const csr_t * const HT, rw_tws_t * const ws, rw_bfs_t * const q, int * const c,
                   const int w) {
  rw_sort(c, w);
  cw_vec_t *known = NULL;
  HASH_FIND(hh, p->codewords, c, (size_t)w * sizeof(int), known);
  if (known) return 0;
  if (codeword_maxc_reached(p)) return -1;
  if (p->irred && w >= 3 && csr_cols_rank(HT, w, c, ws->irr) != w - 1) return 0;
  p->codewords = codeword_add_maybe(p, c, w, CW_SRC_OTHER);
  HASH_FIND(hh, p->codewords, c, (size_t)w * sizeof(int), known);
  if (!known) return 0;
  bfs_push(q, c, w);
  return 1;
}

long long rw_t3_expand(params_t * const p, const rw_t3_t * const t, const csr_t * const HT, rw_tws_t * const ws,
                       const long long limit) {
  if (!t || t->ntrip <= 0 || p->min_w == INT_MAX || !p->codewords) return 0;
  const int wtop = p->min_w + (p->dW > 0 ? p->dW : 0);
  rw_bfs_t q = {NULL, 0, 0, 0};
  for (cw_vec_t *cw = p->codewords; cw != NULL; cw = (cw_vec_t *)(cw->hh.next))
    if (cw->weight <= wtop) bfs_push(&q, cw->arr, cw->weight);
  int * const c = rw_xmalloc(((size_t)wtop + 2) * sizeof(int));
  int * const nb = rw_xmalloc(((size_t)wtop + 2) * sizeof(int));
  long long added = 0;
  int stop = 0;
  while (q.head < q.len && !stop) {
    const int w = q.pool[q.head];
    memcpy(c, q.pool + q.head + 1, (size_t)w * sizeof(int));
    q.head += (size_t)w + 1;
    const uint32_t gen = tws_gen(ws);
    for (int a = 0; a < w; a++) ws->mark[c[a]] = gen;
    /* the equal-weight 2 <-> 2 exchanges (the triples are kept in ws->ent until the next generation) */
    const size_t ne = t3_exchanges(t, ws, c, w, gen);
    for (size_t g0 = 0; g0 < ne && !stop;) {
      size_t g1 = g0 + 1;
      while (g1 < ne && ws->ent[3 * g1] == ws->ent[3 * g0]) g1++;
      for (size_t i = g0; i < g1 && !stop; i++) {
        for (size_t j = i + 1; j < g1 && !stop; j++) {
          const int ua = ws->ent[3 * i + 1], va = ws->ent[3 * i + 2];
          const int ub = ws->ent[3 * j + 1], vb = ws->ent[3 * j + 2];
          if (ua == ub || va == vb) continue;
          for (int a = 0; a < w; a++) nb[a] = (c[a] == ua) ? va : ((c[a] == ub) ? vb : c[a]);
          const int r = bfs_try(p, HT, ws, &q, nb, w);
          if (r < 0) stop = 1;
          else if (r > 0 && ++added >= limit) stop = 2;
        }
      }
      g0 = g1;
    }
    /* the 1 -> 2 splits (weight w + 1); the marks of c are still valid (bfs_try() does not change them) */
    if (w + 1 > wtop) continue;
    for (int a = 0; a < w && !stop; a++) {
      const int z = c[a];
      for (int k = t->tp[z]; k < t->tp[z + 1] && !stop; k++) {
        const int x = t->tq[2 * k], y = t->tq[2 * k + 1];
        if (ws->mark[x] == gen || ws->mark[y] == gen) continue;
        for (int b = 0; b < w; b++) nb[b] = c[b];
        nb[a] = x;
        nb[w] = y;
        const int r = bfs_try(p, HT, ws, &q, nb, w + 1);
        if (r < 0) stop = 1;
        else if (r > 0 && ++added >= limit) stop = 2;
      }
    }
  }
  free(c);
  free(nb);
  free(q.pool);
  return (stop == 2) ? -1 - added : added;
}
