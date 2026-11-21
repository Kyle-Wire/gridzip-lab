// GridZip graft for cmix: grid location, spatial map, symmetry, topic,
// page-structure and multi-scale contexts (see github.com/Kyle-Wire/gridzip-lab).
// Enabled when the environment variable GRIDMAP names a map file written by
// gridzip-lab/src/write_map.py (with topics). GPL, as cmix.
#ifndef GRIDZIP_H
#define GRIDZIP_H
#include <vector>
#include <deque>
#include <cstdio>
#include <cstdlib>
#include "context.h"

class GridState {
 public:
  static GridState& Get() { static GridState g; return g; }
  bool enabled = false;
  unsigned long long mix_ctx = 0;
  enum { NCTX = 19 };
  unsigned long long h[NCTX];
  // GZFLAGS picks the context groups, numbered as lpaq1g's GFLAGS: 1 grid, 2 regions, 4 symmetry,
  // 8 topics, 16 page structure, 64 multi-scale; 128 adds the mixer selected by page structure;
  // 256 adds a mixer selected by fine map region (round 2's best lasting mixer change in lpaq1g);
  // 512 adds two field contexts: template and parameter name, or link target (lpaq1g round 3, GFLAGS 128).
  // 1024 adds a mixer selected by the cell's leaf in the adaptive region tree (lpaq1g round 6, MXSETS 128);
  // it needs a map with the tree appended (work/exp/make_tree_map.py), else the leaf is always 0.
  // 2048 adds two contexts on longer runs of cells: the last four cells + last byte, the last two cells + two bytes
  // (lpaq1g round 15, GFLAGS 8388608: level 9 on 10 MB -0.30% alone, -0.59% with the cell immune memories);
  // 4096 with 2048 keys the first on the last three cells instead (lpaq1g CTXMODE=2, rounds 16-17: better as a plain context).
  // 8192 adds a mixer selected by the current cell x the previous cell's coarse region (lpaq1g round 25, MXSETS 8192:
  // level 9 on 10 MB -0.45% on the full stack, -0.30% beyond a same-size byte-selected weight set).
  // Unset means 255, the graft as first written.
  unsigned flags = 255;
  unsigned long long mix_ctx2 = 0, mix_ctx3 = 0, mix_ctx4 = 0;
  static int Group(int i) { static const int g[NCTX] = {1, 1, 1, 2, 2, 2, 4, 4, 8, 8, 16, 16, 64, 64, 64, 512, 512, 2048, 2048}; return g[i]; }
  bool Uses(int i) const { return enabled && (flags & Group(i)); }

  void Update(unsigned int byte) {
    if (!enabled) return;
    c8_ = c8_ << 8 | (c4_ >> 24); c4_ = c4_ << 8 | byte;
    // grid automaton
    st_ = go_[st_ * 256 + byte]; ++pos_; hist_[pos_ & 65535] = out_[st_];
    int g1 = hist_[pos_ & 65535]; long j2 = pos_ - len_[g1]; int g2 = j2 >= 0 ? hist_[j2 & 65535] : 0;
    long j3 = j2 - len_[g2]; int g3 = j3 >= 0 ? hist_[j3 & 65535] : 0;
    // words
    bool letter = ((byte | 32) >= 'a' && (byte | 32) <= 'z');
    unsigned int done = 0;
    if (letter) word_ = word_ * 773 + (byte | 32); else if (word_) { done = word_; prevword_ = word_; word_ = 0; }
    Structure(byte); MultiScale(byte, done);
    if (len_[g1] >= 3 && tlab_[g1] >= 0) {
      for (int t = 0; t < nt_; ++t) tv_[t] *= 0.95;
      tv_[tlab_[g1]] += 1.0; int b = 0; for (int t = 1; t < nt_; ++t) if (tv_[t] > tv_[b]) b = t; topic_ = b;
    }
    unsigned int o3 = Orb(c4_ & 255) | Orb(c4_ >> 8 & 255) << 8 | Orb(c4_ >> 16 & 255) << 16;
    unsigned int o6 = Orb(c4_ >> 24 & 255) | Orb(c8_ & 255) << 8 | Orb(c8_ >> 8 & 255) << 16;
    unsigned int w = letter ? word_ : 0;
    h[0] = M(1, g1); h[1] = M(2, g1 * 8192ULL + g2); h[2] = M(3, M(g1, c4_ & 0xffffff));
    h[3] = M(4, reg_[g1] * 65536ULL + reg_[g2] * 256 + reg_[g3]); h[4] = M(5, M(g1, reg_[g2] * 256 + reg_[g3]));
    h[5] = M(6, M(regB_[g1] * 1024ULL + regB_[g2], reg_[g3] * 256 + (c4_ & 255)));
    h[6] = M(7, o3); h[7] = M(8, M(o3, o6));
    h[8] = M(9, topic_ * 65536ULL + (c4_ & 0xffff)); h[9] = M(10, M(topic_, w));
    h[10] = M(11, sstate_ * 65536ULL + (c4_ & 0xffff)); h[11] = M(12, M(sstate_ * 7 + link_ * 3 + (tmpl_ > 0), w));
    h[12] = M(13, M(M(sentw_, wc_ > 7 ? 7 : wc_), c4_ & 0xffff)); h[13] = M(14, M(sec_, w)); h[14] = M(15, M(sentw_, w));
    if (flags & 512) { Fields(byte); h[15] = M(16, M(fkey_, fh_)); h[16] = M(17, M(fkey_, c4_ & 0xffff)); }
    if (flags & 2048) { long j4 = j3 - len_[g3]; int g4 = j4 >= 0 ? hist_[j4 & 65535] : 0;
      if (flags & 4096) g4 = 0;
      h[17] = M(18, M(M(g1 * 8192ULL + g2, g3 * 8192ULL + g4), c4_ & 255)); h[18] = M(19, M(g1 * 8192ULL + g2, c4_ & 0xffff)); }
    mix_ctx = sstate_;
    mix_ctx2 = regB_[g1] & 127;
    mix_ctx3 = t256_.empty() ? 0 : t256_[g1];
    mix_ctx4 = g1 * 16ULL + (reg_[g2] & 15);
  }

  bool Load() {
    const char* f = getenv("GRIDMAP"); if (!f) return false;
    FILE* in = fopen(f, "rb"); if (!in) { fprintf(stderr, "cannot open GRIDMAP\n"); exit(1); }
    int m; if (fread(&m, 4, 1, in) != 1) exit(1);
    std::vector<std::vector<unsigned char>> cells(m);
    for (int i = 0; i < m; ++i) { int l = getc(in); cells[i].resize(l); if ((int)fread(&cells[i][0], 1, l, in) != l) exit(1); }
    reg_.resize(m); regB_.resize(m); tlab_.assign(m, -1);
    if ((int)fread(&reg_[0], 4, m, in) != m || (int)fread(&regB_[0], 4, m, in) != m) exit(1);
    if (fread(&nt_, 4, 1, in) == 1) { if ((int)fread(&tlab_[0], 4, m, in) != m) exit(1); } else nt_ = 1;
    int nl;   // optional region tree: three int32[m] cuts (16, 64, 256 leaves); only the 256-leaf one is used
    if (fread(&nl, 4, 1, in) == 1 && nl == 3) { std::vector<int> skip(m); t256_.resize(m);
      if ((int)fread(&skip[0], 4, m, in) != m || (int)fread(&skip[0], 4, m, in) != m || (int)fread(&t256_[0], 4, m, in) != m) exit(1); }
    fclose(in);
    std::vector<std::vector<int>> go(1, std::vector<int>(256, -1)); std::vector<int> term(1, -1);
    for (int c = 0; c < m; ++c) { int s = 0; for (unsigned char ch : cells[c]) {
        if (go[s][ch] < 0) { go[s][ch] = go.size(); go.push_back(std::vector<int>(256, -1)); term.push_back(-1); } s = go[s][ch]; }
      term[s] = c; }
    int S = go.size(); go_.assign(S * 256, 0); out_.assign(S, 0); std::vector<int> fail(S, 0); std::deque<int> q;
    for (int ch = 0; ch < 256; ++ch) { int nx = go[0][ch]; if (nx < 0) nx = 0; go_[ch] = nx; if (nx) q.push_back(nx); }
    while (!q.empty()) { int s = q.front(); q.pop_front(); out_[s] = term[s] >= 0 ? term[s] : out_[fail[s]];
      for (int ch = 0; ch < 256; ++ch) { int nx = go[s][ch]; if (nx < 0) go_[s * 256 + ch] = go_[fail[s] * 256 + ch];
        else { fail[nx] = go_[fail[s] * 256 + ch]; go_[s * 256 + ch] = nx; q.push_back(nx); } } }
    len_.resize(m); for (int c = 0; c < m; ++c) len_[c] = cells[c].size();
    hist_.assign(65536, 0); tv_.assign(64, 0.0);
    for (int i = 0; i < NCTX; ++i) h[i] = 0;
    if (getenv("GZFLAGS")) flags = (unsigned)atoi(getenv("GZFLAGS"));
    enabled = true; return true;
  }

 private:
  std::vector<int> go_, out_, len_, reg_, regB_, tlab_, hist_, t256_; std::vector<double> tv_;
  int st_ = 0, nt_ = 1, topic_ = 0; long pos_ = 0; unsigned int c4_ = 0, c8_ = 0, word_ = 0, prevword_ = 0;
  int sstate_ = 0, tag_ = 0, ent_ = 0, url_ = 0, urlk_ = 0, link_ = 0, tmpl_ = 0, line_ = 0, prev_ = 0, bol_ = 1;
  unsigned int sentw_ = 0, sec_ = 0, head_ = 0; int wc_ = 0, inhead_ = 0, mbol_ = 1;
  // fields (as wiki2_update in lpaq1g.cpp): template/parameter names inside {{ }}, link target inside [[ ]]
  unsigned long long fkey_ = 0, fh_ = 0, tname_[8] = {0}, tpar_[8] = {0}, ltgt_ = 0; int tst_[8] = {0}, fdepth_ = 0, fprev_ = 0, lpart_ = 0;
  void Fields(int c) {
    int p = fprev_; fprev_ = c;
    if (c == '{' && p == '{') { if (fdepth_ < 8) { tname_[fdepth_] = 0; tpar_[fdepth_] = 0; tst_[fdepth_] = 1; } ++fdepth_; fh_ = 0; }
    else if (c == '}' && p == '}' && fdepth_ > 0) { --fdepth_; fh_ = 0; }
    else if (c == '[' && p == '[') { lpart_ = 1; ltgt_ = 0; fh_ = 0; }
    else if (c == ']' && p == ']') { lpart_ = 0; fh_ = 0; }
    else if (c == '|' && lpart_ == 1) { lpart_ = 2; fh_ = 0; }
    else if (fdepth_ > 0 && fdepth_ <= 8 && c == '|' && !lpart_) { tst_[fdepth_ - 1] = 2; tpar_[fdepth_ - 1] = 0; fh_ = 0; }
    else if (fdepth_ > 0 && fdepth_ <= 8 && c == '=' && tst_[fdepth_ - 1] == 2) { tst_[fdepth_ - 1] = 3; fh_ = 0; }
    else {
      if (fdepth_ > 0 && fdepth_ <= 8) { int d = fdepth_ - 1; if (tst_[d] == 1 && c != '{') tname_[d] = M(tname_[d], Orb(c));
        else if (tst_[d] == 2 && c != ' ') tpar_[d] = M(tpar_[d], Orb(c)); }
      if (lpart_ == 1 && c != '[') ltgt_ = M(ltgt_, Orb(c));
      fh_ = M(fh_, c);
    }
    if (fdepth_ > 0 && fdepth_ <= 8) { int d = fdepth_ - 1; fkey_ = M(tname_[d] + 1, tst_[d] == 3 ? tpar_[d] : (unsigned long long)tst_[d]); }
    else if (fdepth_ > 8) fkey_ = 9;
    else fkey_ = lpart_ == 2 ? M(77, ltgt_) : lpart_ == 1 ? 78 : 0;
  }
  static unsigned long long M(unsigned long long a, unsigned long long b) {
    unsigned long long x = a * 0x9E3779B97F4A7C15ULL ^ (b + 0x7F4A7C15ULL) * 0xC2B2AE3D27D4EB4FULL;
    x ^= x >> 29; x *= 0xBF58476D1CE4E5B9ULL; return x ^ (x >> 32);
  }
  static unsigned int Orb(unsigned int c) { if (c >= 65 && c <= 90) return c + 32; if (c >= 48 && c <= 57) return 48; return c; }
  void Structure(int c) {
    if (bol_) { line_ = c; bol_ = 0; } if (c == '\n') { bol_ = 1; line_ = 0; }
    if (c == '<') tag_ = 1; else if (c == '>') tag_ = 0;
    if (c == '&') ent_ = 1; else if (c == ';' || c == ' ') ent_ = 0;
    if (c == '[' && prev_ == '[') link_++; else if (c == ']' && prev_ == ']' && link_ > 0) link_--;
    if (c == '{' && prev_ == '{') tmpl_++; else if (c == '}' && prev_ == '}' && tmpl_ > 0) tmpl_--;
    if (c == 'h') urlk_ = 1; else if (urlk_ == 1 && c == 't') urlk_ = 2; else if (urlk_ == 2 && c == 't') urlk_ = 3;
    else if (urlk_ == 3 && c == 'p') { url_ = 1; urlk_ = 0; } else urlk_ = 0;
    if (c == ' ' || c == ']' || c == '|' || c == '\n') url_ = 0;
    prev_ = c;
    int s = 0;
    if (tag_) s = 1; else if (ent_) s = 7; else if (url_) s = 8; else if (tmpl_) s = 3; else if (link_) s = 2;
    else if (line_ == '|' || line_ == '!' || line_ == '{') s = 4; else if (line_ == '=') s = 5;
    else if (line_ == '*' || line_ == '#' || line_ == ':' || line_ == ';') s = 6; else if (line_ == ' ') s = 9;
    sstate_ = s;
  }
  void MultiScale(int c, unsigned int done) {
    if (done) { if (wc_ == 0) sentw_ = done; wc_++; }
    if (c == '.' || c == '!' || c == '?' || c == '\n') { wc_ = 0; sentw_ = 0; }
    if (mbol_ && c == '=') { inhead_ = 1; head_ = 0; }
    if (inhead_ && c != '=') head_ = head_ * 31 + c;
    if (c == '\n') { if (inhead_) { sec_ = head_; inhead_ = 0; } mbol_ = 1; } else mbol_ = 0;
  }
};

class GridZipContext : public Context {
 public:
  explicit GridZipContext(int which) : which_(which) { context_ = 0; size_ = 1ULL << 32; }
  void Update() override { context_ = GridState::Get().h[which_]; }
 private:
  int which_;
};
#endif
