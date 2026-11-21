/* lpaq1.cpp file compressor, July 24, 2007.
(C) 2007, Matt Mahoney, matmahoney@yahoo.com

    LICENSE

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of
    the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details at
    Visit <http://www.gnu.org/copyleft/gpl.html>.

To compress:   N input output  (N=0..9, uses 3 + 3*2^N MB memory)
To decompress: d input output  (requires same memory)

For example:

  lpaq1 9 foo foo.lpq

compresses foo to foo.lpq using 1.5 GB memory.

  lpaq1 d foo.lpq foo

decompresses to foo, also using 1.5 GB memory.  Option 0 uses
6 MB memory.  The maximum file size is 2 GB.

lpaq1 is a "lite" version of PAQ, about 35 times faster than paq8l
at the cost of some compression (but similar to high end BWT and
PPM compressors).  It is a context-mixing compressor combining 7 
contexts: orders 1, 2, 3, 4, 6, a lowercase unigram word context 
(for ASCII text), and a "match" order, which predicts the next
bit in the last matching context.  The independent bit predictions of
the 7 models are combined using one of 80 neural networks (selected by
a small context), then adjusted using 2 SSE stages (order 0 and 1)
and arithmetic coded.

Prediction is bitwise.  This means that an order-n context consists
of the last n whole bytes plus any of the 0 to 7 previously coded
bits of the current byte starting with the most significant bit.
The unigram word context consists of a hash of the last (at most) 11
consecutive letters (A-Z, a-z) folded to lower case.  The context
does not include any nonalphabetic characters nor any characters
preceding the last nonalphabetic character.

The first 6 contexts (orders 1..4, 6, word) are used to index a
hash table to look up a bit-history represented by an 8-bit state.
The representable states are the same as in paq8l.  A state can
either represent all histories up to 4 bits long, or a pair of
0,1 counts plus a flag to indicate the most recent bit.  The counts
are bounded by (41,0), (40,1), (12,2), (5,3), (4,4) and likewise
for 1,0.  When a count is exceeded, the opposite count is reduced to 
approximately preserve the count ratio.  The last bit flag is present
only for states whose total count is less than 16.  There are 253
possible states.

A bit history is mapped to a probability using an adaptive table
(StateMap).  This differs from paq8l in that each table entry includes
a count so that adaptation is rapid at first.  Each table entry
has a 22-bit probability (initially p = 0.5) and 10-bit count (initially 
n = 0) packed into 32 bits.  After bit y is predicted, n is incremented
up to the limit (1023) and the probability is adjusted by 
p := p + (y - p)/(n + 0.5).  This model is stationary: 
p = (n1 + 0.5)/(n + 1), where n1 is the number of times y = 1 out of n.

The "match" model (MatchModel) looks up the current context in a
hash table, first using a longer context, then a shorter one.  If
a match is found, then the following bits are predicted until there is
a misprediction.  The prediction is computed by mapping the predicted
bit, the length of the match (1..15 or quantized by 4 in 16..62, max 62),
and the last whole byte as context into a StateMap.  If no match is found,
then the order 0 context (last 0..7 bits of the current byte) is used
as context to the StateMap.

The 7 predictions are combined using a neural network (Mixer) as in
paq8l, except it is a single level network without MMX code.  The
inputs p_i, i=0..6 are first stretched: t_i = log(p_i/(1 - p_i)), 
then the output is computed: p = squash(SUM_i t_i * w_i), where
squash(x) = 1/(1 + exp(-x)) is the inverse of stretch().  The weights
are adjusted to reduce the error: w_i := w_i + L * t_i * (y - p) where
(y - p) is the prediction error and L ~ 0.002 is the learning rate.
This is a standard single layer backpropagation network modified to
minimize coding cost rather than RMS prediction error (thus dropping
the factors p * (1 - p) from learning).

One of 80 neural networks are selected by a context that depends on
the 3 high order bits of the last whole byte plus the context order
(quantized to 0, 1, 2, 3, 4, 6, 8, 12, 16, 32).  The order is
determined by the number of nonzero bit histories and the length of
the match from MatchModel.

The Mixer output is adjusted by 2 SSE stages (called APM for adaptive
probability map).  An APM is a StateMap that accepts both a discrte
context and an input probability, pr.  pr is stetched and quantized
to 24 levels.  The output is interpolated between the 2 nearest
table entries, and then only the nearest entry is updated.  The entries
are initialized to p = pr and n = 6 (to slow down initial adaptation)
with a limit n <= 255.  The APM differs from paq8l in that it uses
the new StateMap rapid initial adaptation, does not update both
adjacent table entries, and uses 24 levels instead of 33.  The two
stages use a discrete order 0 context (last 0..7 bits) and a hashed order-1
context (14 bits).  Each output is averaged with its input weighted
by 1/4.

The output is arithmetic coded.  The code for a string s with probability
p(s) is a number between Q and Q+p(x) where Q is the total probability
of all strings lexicographically preceding s.  The number is coded as
a big-endian base-256 fraction.  A header is prepended as follows:

- "pQ" 2 bytes must be present or decompression gives an error.
- 1 (0x01) version number (other values give an error).
- memory option N as one byte '0'..'9' (0x30..0x39).
- file size as a 4 byte big-endian number.
- arithmetic coded data.

Two thirds of the memory (2 * 2^N MB) is used for a hash table mapping
the 6 regular contexts (orders 1-4, 6, word) to bit histories.  A lookup
occurs every 4 bits.  The input is a byte-oriented context plus possibly
the first nibble of the next byte.  The output is an array of 15 bit
histories (1 byte each) for all possible contexts formed by appending
0..3 more bits.  The table entries have the format:

 {checksum, "", 0, 1, 00, 10, 01, 11, 000, 100, 010, 110, 001, 101, 011, 111}

The second byte is the bit history for the context ending on a nibble
boundary.  It also serves as a priority for replacement.  The states
are ordered by increasing total count, where state 0 represents the
initial state (no history).  When a context is looked up, the 8 bit
checksum (part of the hash) is compared with 3 adjacent entries, and
if there is no match, the entry with the lowest priority is cleared
and the new checksum is stored.

The hash table is aligned on 64 byte cache lines.  A table lookup never
crosses a 64 byte address boundary.  Given a 32-bit hash h of the context,
8 bits are used for the checksum and 17 + N bits are used for the
index i.  Then the entries i, i XOR 1, and i XOR 2 are tried.  The hash h
is actually a collision-free permuation, consisting of multiplying the
input by a large odd number mod 2^32, a 16-bit rotate, and another multiply.

The order-1 context is mapped to a bit history using a 64K direct
lookup table, not a hash table.

One third of memory is used by MatchModel, divided equally between 
a rotating input buffer of 2^(N+19) bytes and an index (hash table)
with 2^(N+17) entries.  Two context hashes are maintained, a long one,
h1, of length ceil((N+17)/3) bytes and a shorter one, h2, of length 
ceil((N+17)/5) bytes, where ceil() is the ceiling function.  The index
does not use collision detection.  At each byte boundary, if there is 
not currently a match, then the bytes before the current byte are
compared with the location indexed by h1.  If less than 2 bytes match,
then h2 is tried.  If a match of length 1 or more is found, the
match is maintained until the next bit mismatches the predicted bit.
The table is updated at h1 and h2 after every byte.

To compile (g++ 3.4.5, upx 3.00w):
  g++ -Wall lpaq1.cpp -O2 -Os -march=pentiumpro -fomit-frame-pointer 
      -s -o lpaq1.exe
  upx -qqq lpaq1.exe

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <ctype.h>
#define NDEBUG  // remove for debugging
#include <assert.h>
#include <vector>
#include <deque>
#include <algorithm>

// 8, 16, 32 bit unsigned types (adjust as appropriate)
typedef unsigned char  U8;
typedef unsigned short U16;
typedef unsigned int   U32;

// Error handler: print message if any, and exit
void quit(const char* message=0) {
  if (message) printf("%s\n", message);
  exit(1);
}

// Create an array p of n elements of type T
template <class T> void alloc(T*&p, int n) {
  p=(T*)calloc(n, sizeof(T));
  if (!p) quit("out of memory");
}

///////////////////////////// Squash //////////////////////////////

// return p = 1/(1 + exp(-d)), d scaled by 8 bits, p scaled by 12 bits
int squash(int d) {
  static const int t[33]={
    1,2,3,6,10,16,27,45,73,120,194,310,488,747,1101,
    1546,2047,2549,2994,3348,3607,3785,3901,3975,4022,
    4050,4068,4079,4085,4089,4092,4093,4094};
  if (d>2047) return 4095;
  if (d<-2047) return 0;
  int w=d&127;
  d=(d>>7)+16;
  return (t[d]*(128-w)+t[(d+1)]*w+64) >> 7;
}

//////////////////////////// Stretch ///////////////////////////////

// Inverse of squash. stretch(d) returns ln(p/(1-p)), d scaled by 8 bits,
// p by 12 bits.  d has range -2047 to 2047 representing -8 to 8.  
// p has range 0 to 4095 representing 0 to 1.

class Stretch {
  short t[4096];
public:
  Stretch();
  int operator()(int p) const {
    assert(p>=0 && p<4096);
    return t[p];
  }
} stretch;

Stretch::Stretch() {
  int pi=0;
  for (int x=-2047; x<=2047; ++x) {  // invert squash()
    int i=squash(x);
    for (int j=pi; j<=i; ++j)
      t[j]=x;
    pi=i+1;
  }
  t[4095]=2047;
}

///////////////////////// state table ////////////////////////

// State table:
//   nex(state, 0) = next state if bit y is 0, 0 <= state < 256
//   nex(state, 1) = next state if bit y is 1
//
// States represent a bit history within some context.
// State 0 is the starting state (no bits seen).
// States 1-30 represent all possible sequences of 1-4 bits.
// States 31-252 represent a pair of counts, (n0,n1), the number
//   of 0 and 1 bits respectively.  If n0+n1 < 16 then there are
//   two states for each pair, depending on if a 0 or 1 was the last
//   bit seen.
// If n0 and n1 are too large, then there is no state to represent this
// pair, so another state with about the same ratio of n0/n1 is substituted.
// Also, when a bit is observed and the count of the opposite bit is large,
// then part of this count is discarded to favor newer data over old.

static const U8 State_table[256][2]={
{  1,  2},{  3,  5},{  4,  6},{  7, 10},{  8, 12},{  9, 13},{ 11, 14}, // 0
{ 15, 19},{ 16, 23},{ 17, 24},{ 18, 25},{ 20, 27},{ 21, 28},{ 22, 29}, // 7
{ 26, 30},{ 31, 33},{ 32, 35},{ 32, 35},{ 32, 35},{ 32, 35},{ 34, 37}, // 14
{ 34, 37},{ 34, 37},{ 34, 37},{ 34, 37},{ 34, 37},{ 36, 39},{ 36, 39}, // 21
{ 36, 39},{ 36, 39},{ 38, 40},{ 41, 43},{ 42, 45},{ 42, 45},{ 44, 47}, // 28
{ 44, 47},{ 46, 49},{ 46, 49},{ 48, 51},{ 48, 51},{ 50, 52},{ 53, 43}, // 35
{ 54, 57},{ 54, 57},{ 56, 59},{ 56, 59},{ 58, 61},{ 58, 61},{ 60, 63}, // 42
{ 60, 63},{ 62, 65},{ 62, 65},{ 50, 66},{ 67, 55},{ 68, 57},{ 68, 57}, // 49
{ 70, 73},{ 70, 73},{ 72, 75},{ 72, 75},{ 74, 77},{ 74, 77},{ 76, 79}, // 56
{ 76, 79},{ 62, 81},{ 62, 81},{ 64, 82},{ 83, 69},{ 84, 71},{ 84, 71}, // 63
{ 86, 73},{ 86, 73},{ 44, 59},{ 44, 59},{ 58, 61},{ 58, 61},{ 60, 49}, // 70
{ 60, 49},{ 76, 89},{ 76, 89},{ 78, 91},{ 78, 91},{ 80, 92},{ 93, 69}, // 77
{ 94, 87},{ 94, 87},{ 96, 45},{ 96, 45},{ 48, 99},{ 48, 99},{ 88,101}, // 84
{ 88,101},{ 80,102},{103, 69},{104, 87},{104, 87},{106, 57},{106, 57}, // 91
{ 62,109},{ 62,109},{ 88,111},{ 88,111},{ 80,112},{113, 85},{114, 87}, // 98
{114, 87},{116, 57},{116, 57},{ 62,119},{ 62,119},{ 88,121},{ 88,121}, // 105
{ 90,122},{123, 85},{124, 97},{124, 97},{126, 57},{126, 57},{ 62,129}, // 112
{ 62,129},{ 98,131},{ 98,131},{ 90,132},{133, 85},{134, 97},{134, 97}, // 119
{136, 57},{136, 57},{ 62,139},{ 62,139},{ 98,141},{ 98,141},{ 90,142}, // 126
{143, 95},{144, 97},{144, 97},{ 68, 57},{ 68, 57},{ 62, 81},{ 62, 81}, // 133
{ 98,147},{ 98,147},{100,148},{149, 95},{150,107},{150,107},{108,151}, // 140
{108,151},{100,152},{153, 95},{154,107},{108,155},{100,156},{157, 95}, // 147
{158,107},{108,159},{100,160},{161,105},{162,107},{108,163},{110,164}, // 154
{165,105},{166,117},{118,167},{110,168},{169,105},{170,117},{118,171}, // 161
{110,172},{173,105},{174,117},{118,175},{110,176},{177,105},{178,117}, // 168
{118,179},{110,180},{181,115},{182,117},{118,183},{120,184},{185,115}, // 175
{186,127},{128,187},{120,188},{189,115},{190,127},{128,191},{120,192}, // 182
{193,115},{194,127},{128,195},{120,196},{197,115},{198,127},{128,199}, // 189
{120,200},{201,115},{202,127},{128,203},{120,204},{205,115},{206,127}, // 196
{128,207},{120,208},{209,125},{210,127},{128,211},{130,212},{213,125}, // 203
{214,137},{138,215},{130,216},{217,125},{218,137},{138,219},{130,220}, // 210
{221,125},{222,137},{138,223},{130,224},{225,125},{226,137},{138,227}, // 217
{130,228},{229,125},{230,137},{138,231},{130,232},{233,125},{234,137}, // 224
{138,235},{130,236},{237,125},{238,137},{138,239},{130,240},{241,125}, // 231
{242,137},{138,243},{130,244},{245,135},{246,137},{138,247},{140,248}, // 238
{249,135},{250, 69},{ 80,251},{140,252},{249,135},{250, 69},{ 80,251}, // 245
{140,252},{  0,  0},{  0,  0},{  0,  0}};  // 252
#define nex(state,sel) State_table[state][sel]

//////////////////////////// StateMap, APM //////////////////////////

// A StateMap maps a context to a probability.  Methods:
//
// Statemap sm(n) creates a StateMap with n contexts using 4*n bytes memory.
// sm.p(y, cx, limit) converts state cx (0..n-1) to a probability (0..4095).
//     that the next y=1, updating the previous prediction with y (0..1).
//     limit (1..1023, default 1023) is the maximum count for computing a
//     prediction.  Larger values are better for stationary sources.

class StateMap {
protected:
  const int N;  // Number of contexts
  int cxt;      // Context of last prediction
  U32 *t;       // cxt -> prediction in high 22 bits, count in low 10 bits
  static int dt[1024];  // i -> 16K/(i+3)
  void update(int y, int limit) {
    assert(cxt>=0 && cxt<N);
    int n=t[cxt]&1023, p=t[cxt]>>10;  // count, prediction
    if (n<limit) ++t[cxt];
    else t[cxt]=t[cxt]&0xfffffc00|limit;
    t[cxt]+=(((y<<22)-p)>>3)*dt[n]&0xfffffc00;
  }
public:
  StateMap(int n=256);

  // update bit y (0..1), predict next bit in context cx
  int p(int y, int cx, int limit=1023) {
    assert(y>>1==0);
    assert(cx>=0 && cx<N);
    assert(limit>0 && limit<1024);
    update(y, limit);
    return t[cxt=cx]>>20;
  }
};

int StateMap::dt[1024]={0};

StateMap::StateMap(int n): N(n), cxt(0) {
  alloc(t, N);
  for (int i=0; i<N; ++i)
    t[i]=1<<31;
  if (dt[0]==0)
    for (int i=0; i<1024; ++i)
      dt[i]=16384/(i+i+3);
}

// An APM maps a probability and a context to a new probability.  Methods:
//
// APM a(n) creates with n contexts using 96*n bytes memory.
// a.pp(y, pr, cx, limit) updates and returns a new probability (0..4095)
//     like with StateMap.  pr (0..4095) is considered part of the context.
//     The output is computed by interpolating pr into 24 ranges nonlinearly
//     with smaller ranges near the ends.  The initial output is pr.
//     y=(0..1) is the last bit.  cx=(0..n-1) is the other context.
//     limit=(0..1023) defaults to 255.

class APM: public StateMap {
public:
  APM(int n);
  int pp(int y, int pr, int cx, int limit=255) {
    assert(y>>1==0);
    assert(pr>=0 && pr<4096);
    assert(cx>=0 && cx<N/24);
    assert(limit>0 && limit<1024);
    update(y, limit);
    pr=(stretch(pr)+2048)*23;
    int wt=pr&0xfff;  // interpolation weight of next element
    cx=cx*24+(pr>>12);
    assert(cx>=0 && cx<N-1);
    pr=(t[cx]>>13)*(0x1000-wt)+(t[cx+1]>>13)*wt>>19;
    cxt=cx+(wt>>11);
    return pr;
  }
};

APM::APM(int n): StateMap(n*24) {
  for (int i=0; i<N; ++i) {
    int p=((i%24*2+1)*4096)/48-2048;
    t[i]=(U32(squash(p))<<20)+6;
  }
}

//////////////////////////// Mixer /////////////////////////////

// Mixer m(N, M) combines models using M neural networks with
//     N inputs each using 4*M*N bytes of memory.  It is used as follows:
// m.update(y) trains the network where the expected output is the
//     last bit, y.
// m.add(stretch(p)) inputs prediction from one of N models.  The
//     prediction should be positive to predict a 1 bit, negative for 0,
//     nominally -2K to 2K.
// m.set(cxt) selects cxt (0..M-1) as one of M neural networks to use.
// m.p() returns the output prediction that the next bit is 1 as a
//     12 bit number (0 to 4095).  The normal sequence per prediction is:
//
// - m.add(x) called N times with input x=(-2047..2047)
// - m.set(cxt) called once with cxt=(0..M-1)
// - m.p() called once to predict the next bit, returns 0..4095
// - m.update(y) called once for actual bit y=(0..1).

inline void train(int *t, int *w, int n, int err) {
  for (int i=0; i<n; ++i)
    w[i]+=t[i]*err+0x8000>>16;
}

inline int dot_product(int *t, int *w, int n) {
  int sum=0;
  for (int i=0; i<n; ++i)
    sum+=t[i]*w[i];
  return sum>>8;
}

// ---- mixer knobs (round 2 experiments). Every default reproduces the champion exactly.
// MXLR1/MXLR2: layer-1/layer-2 learning rates (7, 2). MXBOOST/MXHALF: layer-1 rate starts at
// MXLR1+MXBOOST and falls back toward MXLR1 with a half-life of MXHALF KB. MXSETS: extra layer-1
// weight sets (bit 1 region x structure, 2 fine region, 4 topic, 8 sentence word count x structure,
// 16 byte two back x bit position). MXFINAL: final-layer selector (0 c0 x structure, 1 c0, 2 c0 x region).
// MXBIAS=1: a constant input. MXAPM: extra final APM (1 by structure, 2 by region), averaged in.
// Round 6 A: difficulty = moving average of the coding cost of the last few bytes, in half-bit steps (0..15).
// MXSETS bit 32 picks a weight set by it (with page structure); MXAPM=3 adds an APM on it.
static int D_LT[4096], D_cost=0, D_ema=0, D_q=0, D_on=-1;
static int MX_LR1=7, MX_LR2=2, MX_BOOST=0, MX_HALF=256, MX_SETS=0, MX_FINAL=0, MX_BIAS=0, MX_APM=0, MX_lr=7;
static long MX_nbytes=0;
static void load_mx() {
  const char* v;
  if ((v=getenv("MXLR1"))) MX_LR1=atoi(v); if ((v=getenv("MXLR2"))) MX_LR2=atoi(v);
  if ((v=getenv("MXBOOST"))) MX_BOOST=atoi(v); if ((v=getenv("MXHALF"))) MX_HALF=atoi(v);
  if ((v=getenv("MXSETS"))) MX_SETS=atoi(v); if ((v=getenv("MXFINAL"))) MX_FINAL=atoi(v);
  if ((v=getenv("MXBIAS"))) MX_BIAS=atoi(v); if ((v=getenv("MXAPM"))) MX_APM=atoi(v);
  MX_lr=MX_LR1+MX_BOOST;
}

class Mixer {
  // Layer 1: K weight sets over the same N inputs, each selected by its own
  // context (set k -> offset base[k]). Layer 2 (only when K>1): a small mixer
  // over the K layer-1 outputs, selected by its own context.
  // With K==1 this is exactly lpaq1's mixer.
  const int N, M;  // max inputs, total layer-1 weight sets
  int* tx;         // N inputs
  int* wx;         // N*M weights
  int cxt;         // context of set 0 (lpaq API)
  int nx;          // Number of inputs in tx, 0 to N
  int pr;          // final result (scaled 12 bits)
  enum {KMAX=8};
  int K, base[KMAX], cx[KMAX], prk[KMAX];
  int* fw; int fx[KMAX+1]; int fcx, FM;
  int bk, bcx2, bw;  // blended set (MXSETS 512): set bk mixes weights cx[bk] and bcx2, bw/16 and (16-bw)/16
public:
  Mixer(int n, int m);
  void layers(int k, const int* sizes, int final_sets) {  // call once, before use
    K=k; int off=0; for (int i=0;i<k;++i){base[i]=off; off+=sizes[i]; cx[i]=0;}
    if (off>M) quit("mixer sets");
    FM=final_sets; alloc(fw, (K+1)*FM);
    for (int i=0;i<FM;++i) for (int j=0;j<K;++j) fw[i*(K+1)+j]=65536/K;
  }
  void update(int y) {
    if (K<=1) {
      int err=((y<<12)-pr)*MX_lr;
      train(&tx[0], &wx[cxt*N], N, err);
    } else {
      for (int k=0;k<K;++k) { int err=((y<<12)-prk[k])*MX_lr;
        if (k==bk) { train(&tx[0], &wx[(base[k]+cx[k])*N], N, err*bw>>4); train(&tx[0], &wx[(base[k]+bcx2)*N], N, err*(16-bw)>>4); }
        else train(&tx[0], &wx[(base[k]+cx[k])*N], N, err); }
      int err=((y<<12)-pr)*MX_LR2; train(fx, &fw[fcx*(K+1)], K+1, err);
    }
    nx=0;
  }
  void add(int x) { assert(nx<N); tx[nx++]=x; }
  void set(int c) { assert(c>=0 && c<M); cxt=c; cx[0]=c; }
  void setk(int k, int c) { cx[k]=c; }
  void setfinal(int c) { fcx=c; }
  void setblend(int k, int c1, int c2, int w) { bk=k; cx[k]=c1; bcx2=c2; bw=w; }
  int p() {
    if (K<=1) return pr=squash(dot_product(&tx[0], &wx[cxt*N], N)>>8);
    for (int k=0;k<K;++k) { int d=dot_product(&tx[0], &wx[(base[k]+cx[k])*N], N);
      if (k==bk) d=(int)(((long long)d*bw + (long long)dot_product(&tx[0], &wx[(base[k]+bcx2)*N], N)*(16-bw))>>4);
      prk[k]=squash(d>>8); fx[k]=stretch(prk[k]); }
    fx[K]=256;
    return pr=squash(dot_product(fx, &fw[fcx*(K+1)], K+1)>>8);
  }
};

Mixer::Mixer(int n, int m):
    N(n), M(m), tx(0), wx(0), cxt(0), nx(0), pr(2048), K(1), fw(0), fcx(0), FM(0), bk(-1), bcx2(0), bw(16) {
  assert(n>0 && N>0 && M>0);
  alloc(tx, N);
  alloc(wx, N*M);
}

//////////////////////////// HashTable /////////////////////////

// A HashTable maps a 32-bit index to an array of B bytes.
// The first byte is a checksum using the upper 8 bits of the
// index.  The second byte is a priority (0 = empty) for hash
// replacement.  The index need not be a hash.

// HashTable<B> h(n) - create using n bytes  n and B must be 
//     powers of 2 with n >= B*4, and B >= 2.
// h[i] returns array [1..B-1] of bytes indexed by i, creating and
//     replacing another element if needed.  Element 0 is the
//     checksum and should not be modified.

template <int B>
class HashTable {
  U8* t;  // table: 1 element = B bytes: checksum priority data data
  const int N;  // size in bytes
public:
  HashTable(int n);
  U8* operator[](U32 i);
};

template <int B>
HashTable<B>::HashTable(int n): t(0), N(n) {
  assert(B>=2 && (B&B-1)==0);
  assert(N>=B*4 && (N&N-1)==0);
  alloc(t, N+B*4+64);
  t+=64-int(((long)t)&63);  // align on cache line boundary
}

template <int B>
inline U8* HashTable<B>::operator[](U32 i) {
  i*=123456791;
  i=i<<16|i>>16;
  i*=234567891;
  int chk=i>>24;
  i=i*B&N-B;
  if (t[i]==chk) return t+i;
  if (t[i^B]==chk) return t+(i^B);
  if (t[i^B*2]==chk) return t+(i^B*2);
  if (t[i+1]>t[i+1^B] || t[i+1]>t[i+1^B*2]) i^=B;
  if (t[i+1]>t[i+1^B^B*2]) i^=B^B*2;
  memset(t+i, 0, B);
  t[i]=chk;
  return t+i;
}

//////////////////////////// MatchModel ////////////////////////

// MatchModel(n) predicts next bit using most recent context match.
//     using n bytes of memory.  n must be a power of 2 at least 8.
// MatchModel::p(y, m) updates the model with bit y (0..1) and writes
//     a prediction of the next bit to Mixer m.  It returns the length of
//     context matched (0..62).

class MatchModel {
  const int N;  // last buffer index, n/2-1
  const int HN; // last hash table index, n/8-1
  enum {MAXLEN=62};   // maximum match length, at most 62
  U8* buf;    // input buffer
  int* ht;    // context hash -> next byte in buf
  int pos;    // number of bytes in buf
  int match;  // pointer to current byte in matched context in buf
  int len;    // length of match
  U32 h1, h2; // context hashes
  int c0;     // last 0-7 bits of y
  int bcount; // number of bits in c0 (0..7)
  StateMap sm;  // len, bit, last byte -> prediction
  const U32 M1, M2;  // hash multipliers: each byte shifts the hash by log2 of their power of two, so smaller = longer context
public:
  MatchModel(int n, U32 m1=3<<3, U32 m2=5<<5);  // n must be a power of 2 at least 8.
  int p(int y, Mixer& m);  // update bit y (0..1), predict next bit to m
};

MatchModel::MatchModel(int n, U32 m1, U32 m2): N(n/2-1), HN(n/8-1), buf(0), ht(0), pos(0), 
    match(0), len(0), h1(0), h2(0), c0(1), bcount(0), sm(56<<8), M1(m1), M2(m2) {
  assert(n>=8 && (n&n-1)==0);
  alloc(buf, N+1);
  alloc(ht, HN+1);
}

int MatchModel::p(int y, Mixer& m) {

  // update context
  c0+=c0+y;
  ++bcount;
  if (bcount==8) {
    bcount=0;
    h1=h1*M1+c0&HN;
    h2=h2*M2+c0&HN;
    buf[pos++]=c0;
    c0=1;
    pos&=N;

    // find or extend match
    if (len>0) {
      ++match;
      match&=N;
      if (len<MAXLEN) ++len;
    }
    else {
      match=ht[h1];
      if (match!=pos) {
        int i;
        while (len<MAXLEN && (i=match-len-1&N)!=pos
               && buf[i]==buf[pos-len-1&N])
          ++len;
      }
    }
    if (len<2) {
      len=0;
      match=ht[h2];
      if (match!=pos) {
        int i;
        while (len<MAXLEN && (i=match-len-1&N)!=pos
               && buf[i]==buf[pos-len-1&N])
          ++len;
      }
    }
  }

  // predict
  int cxt=c0;
  if (len>0 && (buf[match]+256>>8-bcount)==c0) {
    int b=buf[match]>>7-bcount&1;  // next bit
    if (len<16) cxt=len*2+b;
    else cxt=(len>>2)*2+b+24;
    cxt=cxt*256+buf[pos-1&N];
  }
  else
    len=0;
  m.add(stretch(sm.p(y, cxt)));

  // update index
  if (bcount==0) {
    ht[h1]=pos;
    ht[h2]=pos;
  }
  return len;
}


// ---- Reference model (genomics-style): match only inside the PREVIOUS article,
// which in chain order is the closest relative. No reference index is stored.
class RefModel {
  std::vector<U8> cur, ref; std::vector<int> ht4, ht8;
  int match, len, c0, bcount; U32 c4, c8; StateMap sm;
  enum {HB=20, MAXLEN=60};
public:
  RefModel(): match(-1), len(0), c0(1), bcount(0), c4(0), c8(0), sm(64<<8) { ht4.assign(1<<HB,-1); ht8.assign(1<<HB,-1); }
  int p(int y, Mixer& m) {
    c0+=c0+y; ++bcount;
    if (bcount==8) {
      int c=c0&255; bcount=0; c0=1;
      c8=c8<<8|c4>>24; c4=c4<<8|c;
      cur.push_back(c);
      if (c4==0x6167653Eu && c8==0x20203C70u) {   // "  <page>": new article begins
        ref.swap(cur); cur.clear(); std::fill(ht4.begin(),ht4.end(),-1); std::fill(ht8.begin(),ht8.end(),-1);
        U32 a=0,b=0; for (int i=0;i<(int)ref.size();++i){ b=b<<8|a>>24; a=a<<8|ref[i];
          if (i>=3) ht4[(a*2654435761u)>>(32-HB)]=i+1; if (i>=7) ht8[((a*2654435761u)^(b*0x85EBCA77u))>>(32-HB)]=i+1; }
        match=-1; len=0;
      } else if (!ref.empty()) {
        if (len>0 && match<(int)ref.size() && ref[match]==c) { ++match; if (len<MAXLEN) ++len; }
        else { len=0; match=-1; }
        if (len==0 && cur.size()>=8) { int q=ht8[((c4*2654435761u)^(c8*0x85EBCA77u))>>(32-HB)]; if (q>=0) { match=q; len=8; } }
        if (len==0 && cur.size()>=4) { int q=ht4[(c4*2654435761u)>>(32-HB)]; if (q>=0) { match=q; len=1; } }
        if (match>=(int)ref.size()) { len=0; match=-1; }
      }
    }
    int cxt=0;
    if (len>0 && match>=0 && match<(int)ref.size() && ((ref[match]+256)>>(8-bcount))==c0) {
      int b=ref[match]>>(7-bcount)&1;
      cxt=((len<16?len:15+len/16)*2+b)*256 + (c4&255) + 1;
    } else if (len>0 && bcount>0) len=0;
    m.add(stretch(sm.p(y, cxt)));
    return len;
  }
};


// ---- Immune memory: after a SURPRISE (the actual byte was poorly predicted),
// store "this context -> that byte" with a strength. When the context recurs,
// predict it strongly; correct predictions reinforce, wrong ones decay.
// Round 10 variants (env): IMMCAP strength cap (default 15, up to 63: clonal expansion), IMMTHR surprise threshold in
// bits (default 6). A second memory (IMMUNE2=mode) keys on the map instead of bytes: 1 = the last two cells
// (recognises a token pattern), 2 = the last three cells, 3 = the fine regions of the last three cells + last byte
// (fuzzy recognition: near-variants trigger the same memory).
// Modes 4/6/7/8: the last 4/5/6/7 cells.
static U32 imm_mapkey(int mode);
class ImmuneModel {
  enum {HB=22}; std::vector<U8> mb, ms; std::vector<U32> mk;
  U32 c4, key; int c0, bcount, hit, exp_byte, strength; double cost; StateMap sm; int mode, cap; double thr;
public:
  int pb=-1, ps=0;   // this bit's vote (-1 none, else the bit it predicts) and its strength, for co-stimulation
  ImmuneModel(int md=0): mb(1<<HB), ms(1<<HB), mk(1<<HB), c4(0), key(0), c0(1), bcount(0), hit(0), exp_byte(0), strength(0), cost(0), sm(64*2*256), mode(md),
    cap(getenv("IMMCAP") ? atoi(getenv("IMMCAP")) : 15), thr(getenv("IMMTHR") ? atof(getenv("IMMTHR")) : 6.0) { if (cap>63) cap=63; }
  void bitcost(int y, int pr) { int p = y ? pr : 4096-pr; if (p<1) p=1; cost += -log2(p/4096.0); }
  int p(int y, Mixer& m, int pr_last) {
    bitcost(y, pr_last);
    c0+=c0+y; ++bcount;
    if (bcount==8) {
      int c=c0&255; bcount=0; c0=1;
      U32 idx=key>>(32-HB);
      if (hit) { if (exp_byte==c) { if (ms[idx]<cap) ms[idx]++; } else { if (ms[idx]>0) ms[idx]--; if (ms[idx]==0) mb[idx]=c; } }
      if (cost > thr && !(hit && exp_byte==c)) { mk[idx]=key; mb[idx]=c; ms[idx]=1; }   // surprise: remember
      cost=0;
      c4=c4<<8|c; key = mode ? imm_mapkey(mode)*2654435761u ^ (c4&255)*0x9E3779B1u ^ 0x5bd1e995u : (c4&0xffffff)*2654435761u ^ 0x5bd1e995u;
      U32 j=key>>(32-HB); hit = mk[j]==key && ms[j]>0; exp_byte=mb[j]; strength=ms[j];
    }
    int cxt=0; pb=-1; ps=strength;
    if (hit && ((exp_byte+256)>>(8-bcount))==c0) { int b=exp_byte>>(7-bcount)&1; pb=b; cxt=(strength*2+b)*256+(c4&255)+1; }
    m.add(stretch(sm.p(y, cxt)));
    return hit;
  }
};

// ---- Learned regimes (finance-style regime switching): no hand-written rules.
// Recent text is summarised as a 8-dim profile of byte classes (decaying);
// an online k-means assigns it to one of K learned regimes, which then selects
// the mixer weights in place of the hand-written page-structure tracker.
static int REG_K=0, REG_cur=0; static double REG_prof[8], REG_cent[16][8]; static int REG_n[16];
static inline int bclass(int c){ if(c>=97&&c<=122) return 0; if(c>=65&&c<=90) return 1; if(c>=48&&c<=57) return 2; if(c==32) return 3;
  if(c==10) return 4; if(c=='['||c==']'||c=='{'||c=='}'||c=='|') return 5; if(c=='<'||c=='>'||c=='&'||c==';'||c=='='||c=='/') return 6; return 7; }
static void regime_update(int c) {
  if (!REG_K) return;
  for (int i=0;i<8;++i) REG_prof[i]*=0.9; REG_prof[bclass(c)]+=0.1;
  int best=0; double bd=1e9;
  for (int k=0;k<REG_K;++k){ double d=0; for(int i=0;i<8;++i){double e=REG_prof[i]-REG_cent[k][i]; d+=e*e;} if(d<bd){bd=d;best=k;} }
  REG_n[best]++; double lr = REG_n[best]<1000 ? 1.0/REG_n[best] : 0.001;
  for (int i=0;i<8;++i) REG_cent[best][i]+= lr*(REG_prof[i]-REG_cent[best][i]);
  REG_cur=best;
}

// ---- Multi-scale (renormalisation) state: sentence and section level.
static U32 MS_sent=0, MS_sentw=0, MS_sec=0, MS_head=0; static int MS_inhead=0, MS_wc=0, MS_bol=1;
static void multiscale_update(int c, U32 word_done) {
  if (word_done) { if (MS_wc==0) MS_sentw=word_done; MS_wc++; }
  if (c=='.' || c=='!' || c=='?' || c==10) { MS_sent=0; MS_wc=0; MS_sentw=0; }
  if (MS_bol && c=='=') { MS_inhead=1; MS_head=0; }
  if (MS_inhead && c!='=') MS_head=MS_head*31+c;
  if (c==10) { if (MS_inhead) { MS_sec=MS_head; MS_inhead=0; } MS_bol=1; } else MS_bol=0;
}

//////////////////////////// Predictor /////////////////////////

// A Predictor estimates the probability that the next bit of
// uncompressed data is 1.  Methods:
// Predictor(n) creates with 3*n bytes of memory.
// p() returns P(1) as a 12 bit number (0-4095).
// update(y) trains the predictor with the actual bit (0 or 1).

int MEM=0;  // Global memory usage = 3*MEM bytes (1<<20 .. 1<<29)

class Predictor {
  int pr;  // next prediction
public:
  Predictor();
  int p() const {assert(pr>=0 && pr<4096); return pr;}
  void update(int y);
};

Predictor::Predictor(): pr(2048) {}

// ---- GridZip graft: grid-location, spatial-map and symmetry contexts ----
// Map file (env GRIDMAP): int32 m; m x (uint8 len, bytes); int32 reg[m]; int32 regB[m].
// GFLAGS bitmask: 1 = grid contexts, 2 = spatial map contexts, 4 = symmetry (orbit) contexts.
#include <vector>
#include <deque>
static int GF=0, NCX=6;
static std::vector<int> G_go, G_out, G_len, G_reg, G_regB, G_tlab;
// Round 6 B (GFLAGS 262144): a nested region tree from the map file (after the topics): each cell's node when the
// tree had 16, 64 and 256 leaves. The tree splits the heaviest, most spread leaf first (adaptive refinement).
static std::vector<int> G_t16, G_t64, G_t256;
// Round 7 (GFLAGS 524288): the grid automaton's state (the longest suffix of the text that begins some cell: a
// variable-length context cut at the map's own vocabulary), alone and with the cell before the current one.
// GFLAGS 1048576: next-cell indirect: the cell that followed the current cell the last time, with the automaton state.
static std::vector<int> G_next;
static int G_p2=0, G_p3=0;
static std::vector<int> G_lastb;   // round 36: each cell's last byte (-1 for an empty cell)   // round 32: the cell trigram's guess (GFLAGS 67108864), read by MXSETS 131072
// Round 7b (GFLAGS 2097152; MXSETS 256): the correlogram's two signals. Each cell's position in the map plane
// (two 4-bit coordinates, an 'XY' section after the topics); the direction of the last step between cells,
// and the recent place: an exponential average of positions (~50 bytes) on an 8 x 8 grid.
// MXSETS 512 (the user's geographically weighted mixing): a weight set blended between the cell's nearest and
// second-nearest fine-region centroids, by inverse squared distance (G_w2 = weight of the nearest, in 16ths).
static std::vector<int> G_x, G_y, G_rb2, G_w2; static int XY_px=0, XY_py=0, XY_dir=8, XY_place=0;
// Fluid map (round 9, GFLAGS 4194304 contexts, MXSETS 1024 weight set): every cell is a particle on a 2-D torus
// (16-bit coordinates). Mass = how often the cell has ended a byte. Each byte, the cell that just ended (b) and the
// cell before it (a) pull toward each other, each moving in proportion to the other's mass, so heavy cells barely move
// (inertia); b is then pushed down the gradient of a 32x32 grid of particle counts (pressure, so the fluid doesn't
// collapse). Start: hashed positions (FLINIT=0, a "big bang", nothing sent) or the map plane (FLINIT=1, XY section).
// FLKA attraction (/256 of the gap, default 4), FLKP pressure step per count of difference (default 16, of 65536),
// FLQ place grid bits per axis (default 4: 256 places). Encoder and decoder run the same integer physics.
static std::vector<unsigned short> FL_x, FL_y; static std::vector<int> FL_m, FL_d;
static int FL_on=0, FL_ka=4, FL_kp=16, FL_init=0, FL_q=4, FL_p1=0, FL_p2=0;
static inline int fl_sq(int c) { return (FL_x[c]>>11)*32 + (FL_y[c]>>11); }
static inline int fl_place(int c) { return ((FL_x[c]>>(16-FL_q)) << FL_q) + (FL_y[c]>>(16-FL_q)); }
static void fl_step(int a, int b) {
  if (FL_m[b] < (1<<30)) ++FL_m[b];
  int s0 = fl_sq(b);
  if (a != b && FL_ka) {
    int dx = (short)(FL_x[a]-FL_x[b]), dy = (short)(FL_y[a]-FL_y[b]);      // the short way round the torus
    long long ma = FL_m[a], mb = FL_m[b], t = (ma+mb) << 8;
    int sa = fl_sq(a);
    FL_x[b] += (int)(dx*FL_ka*ma/t); FL_y[b] += (int)(dy*FL_ka*ma/t);
    FL_x[a] -= (int)(dx*FL_ka*mb/t); FL_y[a] -= (int)(dy*FL_ka*mb/t);
    int sa2 = fl_sq(a); if (sa2 != sa) { --FL_d[sa]; ++FL_d[sa2]; } }
  if (FL_kp) { int s = fl_sq(b), x = s>>5, y = s&31;
    int gx = FL_d[((x+31)&31)*32+y] - FL_d[((x+1)&31)*32+y], gy = FL_d[x*32+((y+31)&31)] - FL_d[x*32+((y+1)&31)];
    FL_x[b] += gx*FL_kp; FL_y[b] += gy*FL_kp; }
  int s1 = fl_sq(b); if (s1 != s0) { --FL_d[s0]; ++FL_d[s1]; }
}
static int G_NT=1, G_topic=0; static double G_tv[64]; static double G_decay=0.95;
static std::vector<int> G_hist;     // longest cell ending at each position
static int G_state=0;
static long G_pos=0;
static U32 imm_mapkey(int mode) {
  if (G_len.empty()) return 0;
  int g1=G_hist[G_pos&65535]; long j2=G_pos-G_len[g1]; int g2= j2>=0? G_hist[j2&65535]:0;
  long j3=j2-G_len[g2]; int g3= j3>=0? G_hist[j3&65535]:0;
  if (mode==1) return (U32)(g1*8192+g2)*0x2545F491u;
  if (mode==10) return (U32)(G_p3*8192+g1)*0x2545F491u ^ (U32)g2*0x85EBCA77u;   // round 35: the 4-gram's guessed next cell x the last two cells
  if (mode==11) return (U32)(G_p3*8192+G_p2)*0x2545F491u ^ (U32)g1*0x85EBCA77u;  // both guesses x current cell
  if (mode==9) return (U32)(g1*16+(G_reg[g2]&15))*0x2545F491u;   // round 27: current cell x previous cell's coarse region
  if (mode>=5 && mode<=8) { int n = mode-1; U32 hh = (U32)(g1*8192+g2)*0x2545F491u ^ (U32)g3*0x85EBCA77u; long j=j3; int g=g3;   // 5..8: the last 4..7 cells
    for (int i=3;i<n && j>=0;++i) { j-=G_len[g]; g = j>=0 ? G_hist[j&65535] : 0; hh = (hh + (U32)g*0xC2B2AE3Du)*0x27d4eb2fu; } return hh; }
  if (mode==4) { long j4=j3-G_len[g3]; int g4= j4>=0? G_hist[j4&65535]:0; return (((U32)(g1*8192+g2)*0x2545F491u) ^ (U32)g3*0x85EBCA77u) + (U32)g4*0xC2B2AE3Du; }
  if (mode==2) return ((U32)(g1*8192+g2)*0x2545F491u) ^ (U32)g3*0x85EBCA77u;
  return (U32)((G_regB[g1]&127)*16384+(G_regB[g2]&127)*128+(G_regB[g3]&127))*0x27D4EB2Fu;
}
// ---- evolvable contexts: GENOME="fmask:bmask,..." (hex). Each gene is one context:
// bmask bit i = include byte c(i+1) (i=0..7); fmask bits: 0 g1, 1 g2, 2 r1, 3 r2, 4 r3,
// 5 f1, 6 f2, 7 word, 8 previous word, 9 page structure, 10 topic, 11 fold bytes (symmetry), 12 link/template depth
static int NGENE=0; static int GENE_F[24], GENE_B[24];
static U32 W_cur=0, W_prev=0;
static void load_genome() {
  const char* g=getenv("GENOME"); if (!g) return;
  const char* p=g;
  while (*p && NGENE<17) { unsigned f=0,b=0; if (sscanf(p,"%x:%x",&f,&b)!=2) break; GENE_F[NGENE]=f; GENE_B[NGENE]=b; ++NGENE;
    const char* c=strchr(p,','); if (!c) break; p=c+1; }
}
static void load_grid() {
  const char* f=getenv("GRIDMAP"); const char* fl=getenv("GFLAGS");
  if (!f || !fl) return;
  GF=atoi(fl);
  FILE* in=fopen(f,"rb"); if (!in) quit("cannot open GRIDMAP");
  int m; if (fread(&m,4,1,in)!=1) quit("bad map");
  std::vector<std::vector<unsigned char> > cells(m);
  for (int i=0;i<m;++i){int l=getc(in); cells[i].resize(l); if ((int)fread(&cells[i][0],1,l,in)!=l) quit("bad map");}
  G_lastb.resize(m); for (int i=0;i<m;++i) G_lastb[i] = cells[i].empty() ? -1 : cells[i].back();
  G_reg.resize(m); G_regB.resize(m);
  if ((int)fread(&G_reg[0],4,m,in)!=m || (int)fread(&G_regB[0],4,m,in)!=m) quit("bad map");
  if (GF&8) { if (fread(&G_NT,4,1,in)!=1) quit("map has no topics"); G_tlab.resize(m);
    if ((int)fread(&G_tlab[0],4,m,in)!=m) quit("bad topics"); const char* d=getenv("GDECAY"); if (d) G_decay=atof(d); }
  // Optional sections after the topics, in this order: XY (tag 0x5859; per cell x<<4|y, second-nearest fine region,
  // weight of the nearest in 16ths) and the region tree (tag 3; 16-, 64- and 256-leaf labels).
  const char* mxs=getenv("MXSETS"); int ms = mxs ? atoi(mxs) : 0;
  { const char* v; if ((v=getenv("FLKA"))) FL_ka=atoi(v); if ((v=getenv("FLKP"))) FL_kp=atoi(v);
    if ((v=getenv("FLINIT"))) FL_init=atoi(v); if ((v=getenv("FLQ"))) FL_q=atoi(v); if (FL_q<1 || FL_q>6) quit("FLQ 1..6"); }
  FL_on = (GF&4194304) || (ms&1024);
  int want_xy = (GF&2097152) || (ms&(256|512)) || (FL_on && FL_init==1), want_tree = (GF&262144) || (ms&(64|128));
  if ((want_xy || want_tree) && !(GF&8)) { int nt; if (fread(&nt,4,1,in)!=1) quit("map has no topics"); fseek(in, (long)m*4, SEEK_CUR); }
  int tag=0, have=0; if (want_xy || want_tree) have = fread(&tag,4,1,in)==1;
  if (have && tag==0x5859) {
    if (want_xy) { G_x.resize(m); G_y.resize(m); G_rb2.resize(m); G_w2.resize(m);
      for (int i=0;i<m;++i) { int b=getc(in), r2=getc(in), w=getc(in); if (w<0) quit("bad XY"); G_x[i]=b>>4; G_y[i]=b&15; G_rb2[i]=r2&127; G_w2[i]=w; } }
    else fseek(in, (long)m*3, SEEK_CUR);
    have = fread(&tag,4,1,in)==1; }
  else if (want_xy) quit("map has no XY section");
  if (want_tree) { if (!have || tag!=3) quit("map has no tree");
    G_t16.resize(m); G_t64.resize(m); G_t256.resize(m);
    if ((int)fread(&G_t16[0],4,m,in)!=m || (int)fread(&G_t64[0],4,m,in)!=m || (int)fread(&G_t256[0],4,m,in)!=m) quit("bad tree"); }
  fclose(in);
  // Aho-Corasick automaton: longest cell that ends at the current byte
  std::vector<std::vector<int> > go(1, std::vector<int>(256,-1));
  std::vector<int> term(1,-1);
  for (int c=0;c<m;++c){int s=0; for (size_t j=0;j<cells[c].size();++j){int ch=cells[c][j];
      if (go[s][ch]<0){go[s][ch]=go.size(); go.push_back(std::vector<int>(256,-1)); term.push_back(-1);} s=go[s][ch];}
    term[s]=c;}
  int S=go.size(); G_go.assign(S*256,0); G_out.assign(S,0); std::vector<int> fail(S,0);
  std::deque<int> q;
  for (int ch=0;ch<256;++ch){int nx=go[0][ch]; if(nx<0) nx=0; G_go[ch]=nx; if(nx){fail[nx]=0;q.push_back(nx);}}
  while(!q.empty()){int s=q.front();q.pop_front(); G_out[s]= term[s]>=0? term[s] : G_out[fail[s]];
    for(int ch=0;ch<256;++ch){int nx=go[s][ch]; if(nx<0) G_go[s*256+ch]=G_go[fail[s]*256+ch];
      else {fail[nx]=G_go[fail[s]*256+ch]; G_go[s*256+ch]=nx; q.push_back(nx);}}}
  G_len.resize(m); for(int c=0;c<m;++c) G_len[c]=cells[c].size();
  G_hist.assign(1<<16,0); G_next.assign(m,0);
  if (FL_on) { FL_x.resize(m); FL_y.resize(m); FL_m.assign(m,1); FL_d.assign(1024,0);
    for (int c=0;c<m;++c) {
      if (FL_init==1) { FL_x[c]=(G_x[c]<<12)+2048; FL_y[c]=(G_y[c]<<12)+2048; }
      else { U32 h=(U32)(c+1)*0x9E3779B1u; h^=h>>15; h*=0x85EBCA77u; h^=h>>13; h*=0xC2B2AE3Du; h^=h>>16; FL_x[c]=h>>16; FL_y[c]=h&0xffff; }
      ++FL_d[fl_sq(c)]; } }
  if (GF&1) NCX+=3;
  if (GF&2) NCX+=3;
  if (GF&4) NCX+=2;
  if (GF&8) NCX+=2;
  if (GF&16) NCX+=2;
  if (GF&64) NCX+=3;
  if (GF&128) NCX+=2;
  if (GF&256) NCX+=2;
  if (GF&512) NCX+=1;
  if (GF&1024) NCX+=1;
  if (GF&2048) NCX+=1;
  if (GF&4096) NCX+=1;
  if (GF&8192) NCX+=2;
  if (GF&16384) NCX+=2;
  if (GF&32768) NCX+=3;
  if (GF&65536) NCX+=1;
  if (GF&131072) NCX+=1;
  if (GF&262144) NCX+=3;
  if (GF&524288) NCX+=2;
  if (GF&1048576) NCX+=1;
  if (GF&2097152) NCX+=2;
  if (GF&4194304) NCX+=2;
  if (GF&8388608) NCX+=2;
  if (GF&33554432) NCX+=2;
  if (GF&67108864) NCX+=2;
  if (GF&134217728) NCX+=2;
  if (GF&268435456) NCX+=2;
  if (GF&536870912) NCX+=2;
  if (GF&1073741824) NCX+=2;  // round 37: next-cell guess keyed on the fine regions of the last three cells  // round 36: the byte the cell guesses expect next, with its hit history  // round 33: confidence in the cell guesses (were the last ones right, do they agree)  // round 31: cell 4-gram: the cell that followed the last three cells the last time  // round 30: cell trigram: the cell that followed the last two cells the last time  // round 29: cell recency (bytes since the current / previous cell last ended here)
  if (GF&16777216) NCX+=2;  // round 18: the two cells before the current one + 2 bytes (skip), three cells + 2 bytes   // round 15: plain contexts on the last four cells (+ last byte) and the last two cells (+ 2 bytes)
  if (NCX>40) quit("too many contexts for NMAX");
  load_genome(); if (NGENE) NCX=1+NGENE;
}
static inline U32 hmix(U32 a, U32 b) {
  U32 h=a*0x9E3779B1u ^ (b+0x7F4A7C15u)*0x85EBCA77u; h^=h>>15; h*=0xC2B2AE3Du; return h^(h>>13);
}
// wiki-page structure: 0 prose, 1 xml tag, 2 [[link]], 3 {{template}}, 4 table line,
// 5 heading line, 6 list line, 7 &entity;, 8 url, 9 xml/text header lines
static int S_state=0, S_tag=0, S_link=0, S_tmpl=0, S_ent=0, S_url=0, S_line=0, S_prev=0, S_bol=1;
static int S_urlk=0;
static void struct_update(int c) {
  if (S_bol) { S_line = c; S_bol = 0; }
  if (c=='\n') { S_bol=1; S_line=0; }
  if (c=='<') S_tag=1; else if (c=='>') S_tag=0;
  if (c=='&') S_ent=1; else if (c==';' || c==' ') S_ent=0;
  if (c=='[' && S_prev=='[') S_link++; else if (c==']' && S_prev==']' && S_link>0) S_link--;
  if (c=='{' && S_prev=='{') S_tmpl++; else if (c=='}' && S_prev=='}' && S_tmpl>0) S_tmpl--;
  // url: "http" then until space or ] or |
  if (c=='h') S_urlk=1; else if (S_urlk==1 && c=='t') S_urlk=2; else if (S_urlk==2 && c=='t') S_urlk=3;
  else if (S_urlk==3 && c=='p') { S_url=1; S_urlk=0; } else S_urlk=0;
  if (c==' ' || c==']' || c=='|' || c=='\n') S_url=0;
  S_prev=c;
  int st=0;
  if (S_tag) st=1; else if (S_ent) st=7; else if (S_url) st=8; else if (S_tmpl) st=3; else if (S_link) st=2;
  else if (S_line=='|' || S_line=='!' || S_line=='{') st=4; else if (S_line=='=') st=5;
  else if (S_line=='*' || S_line=='#' || S_line==':' || S_line==';') st=6;
  else if (S_line==' ' ) st=9;
  S_state=st;
}
static inline int orb(int c){ if(c>=65&&c<=90) return c+32; if(c>=48&&c<=57) return 48; return c; }

// Round 3 contexts (GFLAGS 128 fields, 256 table columns, 512 previous word).
// Fields: inside {{template|param=value}} the template and parameter names; inside [[target|label]]
// the target while reading the label. F_h hashes the bytes since the field began.
// Table columns: in a {| ... |} table, the cell's column (cells since the last |- row line) and the
// bytes since the cell began.
static U32 F_key=0, F_h=0, B_prev=0, B_cur=0, C_cellh=0;
// Round 4 (GFLAGS 1024 two words back, 2048 the word before last skipping the previous one,
// 4096 previous word with the last two bytes, 8192 indirect (what followed these bytes before), 16384 sparse bytes)
static U32 B_prev2=0, I_t1[256]; static U16 I_t2[65536];
// Round 5 (GFLAGS 32768 more indirect masks, 65536 indirect after three bytes, 131072 word indirect:
// the word that followed the previous word the last time it appeared)
static U16 I_t3[65536]; static U32 I_w[65536];
static int F_depth=0, F_prev=0, L_part=0, C_intab=0, C_col=-1, C_bol=1, C_l1=0, C_pos=0;
static U32 T_name[8], T_par[8], L_tgt=0; static int T_st[8];   // T_st: 1 template name, 2 param name, 3 param value
static void wiki2_update(int c) {
  int p=F_prev; F_prev=c;
  if (c=='{' && p=='{') { if (F_depth<8) { T_name[F_depth]=0; T_par[F_depth]=0; T_st[F_depth]=1; } ++F_depth; F_h=0; }
  else if (c=='}' && p=='}' && F_depth>0) { --F_depth; F_h=0; }
  else if (c=='[' && p=='[') { L_part=1; L_tgt=0; F_h=0; }
  else if (c==']' && p==']') { L_part=0; F_h=0; }
  else if (c=='|' && L_part==1) { L_part=2; F_h=0; }
  else if (F_depth>0 && F_depth<=8 && c=='|' && !L_part) { int d=F_depth-1; T_st[d]=2; T_par[d]=0; F_h=0; }
  else if (F_depth>0 && F_depth<=8 && c=='=' && T_st[F_depth-1]==2) { T_st[F_depth-1]=3; F_h=0; }
  else {
    if (F_depth>0 && F_depth<=8) { int d=F_depth-1; if (T_st[d]==1 && c!='{') T_name[d]=hmix(T_name[d], orb(c));
      else if (T_st[d]==2 && c!=' ') T_par[d]=hmix(T_par[d], orb(c)); }
    if (L_part==1 && c!='[') L_tgt=hmix(L_tgt, orb(c));
    F_h=hmix(F_h, c);
  }
  if (F_depth>0 && F_depth<=8) { int d=F_depth-1; F_key=hmix(T_name[d]+1, T_st[d]==3 ? T_par[d] : (U32)T_st[d]); }
  else if (F_depth>8) F_key=9;
  else F_key = L_part==2 ? hmix(77, L_tgt) : L_part==1 ? 78 : 0;
  // tables
  if (C_bol) { C_l1=c; C_pos=0; if (C_intab && (c=='|' || c=='!')) { ++C_col; C_cellh=0; } }
  else if (C_pos==1) {
    if (C_l1=='{' && c=='|') { C_intab=1; C_col=-1; }
    else if (C_l1=='|' && c=='}') { C_intab=0; C_col=-1; }
    else if (C_l1=='|' && c=='-') C_col=-1;
    else C_cellh=hmix(C_cellh, c);
  }
  else if (C_intab && (c=='|' || c=='!') && p==c && F_depth==0 && !L_part) { ++C_col; C_cellh=0; }
  else C_cellh=hmix(C_cellh, c);
  ++C_pos; C_bol = c=='\n';
  // previous word
  if ((c|32)>='a' && (c|32)<='z') B_cur=hmix(B_cur, c|32); else if (B_cur) { B_prev2=B_prev; B_prev=B_cur; B_cur=0; I_w[B_prev2>>16]=B_prev; }
}

void Predictor::update(int y) {
  enum {NMAX=40};
  static U8 t0[0x10000];  // order 1 cxt -> state
  static HashTable<16> t(MEM*2);  // cxt -> state
  static int c0=1;  // last 0-7 bits with leading 1
  static int c4=0;  // last 4 bytes
  static int c8=0;  // the 4 bytes before those
  static U8 *cp[NMAX];
  static int init=0;
  static int bcount=0;  // bit count
  static StateMap sm[NMAX];
  static APM a1(0x100), a2(0x4000), a3(0x1000);
  static U32 h[NMAX];
  static Mixer* mp=0;
  static MatchModel mm(MEM);  // predicts next bit by matching context
  static MatchModel* mm2 = getenv("MATCH2") ? new MatchModel(MEM, 3<<1, 3<<2) : 0;  // round 6 C: longer contexts (~23 and ~12 bytes)
  static RefModel* rm = getenv("REFMODEL") ? new RefModel() : 0;
  static ImmuneModel* im = getenv("IMMUNE") ? new ImmuneModel() : 0;
  static ImmuneModel* im2 = getenv("IMMUNE2") ? new ImmuneModel(atoi(getenv("IMMUNE2"))) : 0;
  static ImmuneModel* im3 = getenv("IMMUNE3") ? new ImmuneModel(atoi(getenv("IMMUNE3"))) : 0;   // a third memory, another key
  static int regk = getenv("REGIME") ? atoi(getenv("REGIME")) : 0;
  // Co-stimulation (IMMCONS=1): one input from the memories' votes together: for each memory none / 0 / 1, the
  // strongest voter's strength, the bit position. Agreement between a byte memory and a place memory is evidence
  // neither has alone. IMMCONS=2 also adds the previous byte.
  static int immcons = getenv("IMMCONS") ? atoi(getenv("IMMCONS")) : 0; static StateMap* csm = immcons ? new StateMap(27*16*8*256) : 0;
  static U32 msw_cur=0;
  if (!init) { load_grid(); load_mx(); for (int i=0;i<NMAX;++i) cp[i]=t0, h[i]=0; { int s0=80*(GF&8 ? 64 : 1)*(GF&16 ? 10 : 1);
      int nin=NCX+1+(getenv("REFMODEL")?1:0)+(getenv("IMMUNE")?1:0)+(getenv("IMMUNE2")?1:0)+(getenv("IMMUNE3")?1:0)+(immcons?1:0)+(MX_BIAS?1:0)+(getenv("MATCH2")?1:0);
      if (getenv("MIX2")) { int sz[8]={s0, 256, 256*8, 32*10}; int K=4;
        if ((MX_SETS&1) && K<8) sz[K++]=16*10; if ((MX_SETS&2) && K<8) sz[K++]=128; if ((MX_SETS&4) && K<8) sz[K++]=64;
        if ((MX_SETS&8) && K<8) sz[K++]=8*10; if ((MX_SETS&16) && K<8) sz[K++]=256*8;
        if ((MX_SETS&32) && K<8) sz[K++]=16*10; if ((MX_SETS&64) && K<8) sz[K++]=64; if ((MX_SETS&128) && K<8) sz[K++]=256;
        if ((MX_SETS&256) && K<8) sz[K++]=64; if ((MX_SETS&512) && K<8) sz[K++]=128;
        if ((MX_SETS&1024) && K<8) sz[K++]=1<<(2*FL_q);
        if ((MX_SETS&2048) && K<8) sz[K++]=4096;   // round 23: a weight set per map cell (cell id & 4095)
        if ((MX_SETS&4096) && K<8) sz[K++]=4096;   // round 24 control: a weight set by the last two bytes, hashed to 4096
        if ((MX_SETS&8192) && K<8) sz[K++]=20480;  // round 24: by current cell x previous cell's coarse region
        if ((MX_SETS&16384) && K<8) sz[K++]=20480; // round 25 control: by the last two bytes, hashed to 20,480 (same size)
        if ((MX_SETS&32768) && K<8) sz[K++]=1280*128; // round 26: current cell x previous cell's fine region
        if ((MX_SETS&65536) && K<8) sz[K++]=65536;    // round 26: current cell x previous cell, hashed to 65,536
        if ((MX_SETS&131072) && K<8) sz[K++]=2048;    // round 32: guessed next cell's fine region x current cell's coarse region
        int tot=0; for (int i=0;i<K;++i) tot+=sz[i];
        mp=new Mixer(nin, tot); mp->layers(K, sz, MX_FINAL==2 ? 256*16 : MX_FINAL>=3 ? 256*256*(MX_FINAL==4 ? 10 : 1) : 256*10); }
      else mp=new Mixer(nin, s0); } init=1; }
  Mixer& m=*mp;
  const int N=NCX;
  if (D_on<0) { D_on = (MX_SETS&32) || MX_APM==3; for (int i=0;i<4096;++i) D_LT[i]=(int)(-log((i+0.5)/4096.0)/log(2.0)*256+0.5); }
  if (D_on) { int q=pr<0?0:pr>4095?4095:pr; D_cost += y ? D_LT[q] : D_LT[4095-q]; }

  for (int i=0;i<N;++i) *cp[i]=nex(*cp[i], y);
  m.update(y);

  ++bcount;
  c0+=c0+y;
  if (c0>=256) {
    c0-=256;
    c8=c8<<8|(c4>>24&255);
    c4=c4<<8|c0;
    int byte=c0;
    ++MX_nbytes;
    if (D_on) { D_ema += (D_cost*16 - D_ema) >> 3; D_cost=0; D_q = D_ema>>11; if (D_q>15) D_q=15; }
    if (MX_BOOST) MX_lr=MX_LR1+(int)((long long)MX_BOOST*MX_HALF*1024/((long long)MX_HALF*1024+MX_nbytes));
    h[0]=c0<<8;  // order 1
    h[1]=(c4&0xffff)<<5|0x57000000;  // order 2
    h[2]=(c4<<8)*3;  // order 3
    h[3]=c4*5;  // order 4
    h[4]=h[4]*(11<<5)+c0*13&0x3fffffff;  // order 6
    if (c0>=65 && c0<=90) c0+=32;  // lowercase unigram word order
    if (c0>=97 && c0<=122) h[5]=(h[5]+c0)*(7<<3);
    else h[5]=0;
    int k=6;
    if (GF) {
      if (regk) { if (!REG_K) { REG_K=regk; for(int a=0;a<regk;++a) for(int b=0;b<8;++b) REG_cent[a][b]=(a*7+b*3)%11/11.0; } regime_update(byte); }
      G_state=G_go[G_state*256+byte]; ++G_pos; G_hist[G_pos&65535]=G_out[G_state];
      int g1=G_hist[G_pos&65535]; long j2=G_pos-G_len[g1]; int g2= j2>=0? G_hist[j2&65535]:0;
      long j3=j2-G_len[g2]; int g3= j3>=0? G_hist[j3&65535]:0;
      if (GF&1) { h[k++]=hmix(101,g1); h[k++]=hmix(102,g1*8192+g2); h[k++]=hmix(103,hmix(g1,c4&0xffffff)); }
      if (GF&2) { int r1=G_reg[g1],r2=G_reg[g2],r3=G_reg[g3],f1=G_regB[g1],f2=G_regB[g2],f3=G_regB[g3];
        h[k++]=hmix(104,r1*65536+r2*256+r3); h[k++]=hmix(105,hmix(g1,r2*256+r3)); h[k++]=hmix(106,hmix(f1*1024+f2,f3*256+(c4&255))); }
      if (GF&4) { int o=orb(c4&255)|orb(c4>>8&255)<<8|orb(c4>>16&255)<<16;
        int o2=orb(c4>>24&255)|orb(c8&255)<<8|orb(c8>>8&255)<<16;
        h[k++]=hmix(107,o); h[k++]=hmix(108,hmix(o,o2)); }
      if (GF&8) {
        int gn=G_hist[G_pos&65535];
        if (G_len[gn]>=3 && G_tlab[gn]>=0) {   // content cells vote for the topic
          for (int t=0;t<G_NT;++t) G_tv[t]*=G_decay;
          G_tv[G_tlab[gn]]+=1.0; int b=0; for (int t=1;t<G_NT;++t) if (G_tv[t]>G_tv[b]) b=t; G_topic=b; }
        h[k++]=hmix(109,G_topic*65536+(c4&0xffff)); h[k++]=hmix(110,hmix(G_topic,h[5])); }
      if (GF&64) { U32 wd=0; if ((byte|32)>=97 && (byte|32)<=122) msw_cur=h[5]; else if (msw_cur) { wd=msw_cur; msw_cur=0; }
        multiscale_update(byte, wd);
        h[k++]=hmix(113,hmix(MS_sentw, MS_wc>7?7:MS_wc)*0x10000+(c4&0xffff));
        h[k++]=hmix(114,hmix(MS_sec, h[5]));
        h[k++]=hmix(115,hmix(MS_sentw, h[5])); }
      if (GF&16) { struct_update(byte);
        h[k++]=hmix(111,S_state*65536+(c4&0xffff)); h[k++]=hmix(112,hmix(S_state*7+S_link*3+(S_tmpl>0),h[5])); }
      if (GF&(128|256|512|1024|2048|4096|131072)) wiki2_update(byte);
      if (GF&128) { h[k++]=hmix(201,hmix(F_key,F_h)); h[k++]=hmix(202,hmix(F_key,c4&0xffff)); }
      if (GF&256) { int col = C_intab ? (C_col<0 ? 63 : C_col>62 ? 62 : C_col) : 64;
        h[k++]=hmix(203,hmix(col,C_cellh)); h[k++]=hmix(204,hmix(col,c4&0xffffff)); }
      if (GF&512) h[k++]=hmix(205,hmix(B_prev,h[5]));
      if (GF&1024) h[k++]=hmix(206,hmix(B_prev2,hmix(B_prev,h[5])));
      if (GF&2048) h[k++]=hmix(207,hmix(B_prev2,h[5]));
      if (GF&4096) h[k++]=hmix(208,hmix(B_prev,c4&0xffff));
      if (GF&(8192|32768)) { U32 d=c4&0xffff, cc=d&255; I_t1[d>>8]=I_t1[d>>8]<<8|cc; I_t2[c4>>8&0xffff]=(U16)(I_t2[c4>>8&0xffff]<<8|cc);
        if (GF&8192) { h[k++]=hmix(209,(cc|I_t1[cc]<<8)&0xffffff); h[k++]=hmix(210,d|(U32)I_t2[d]<<16); }
        if (GF&32768) { h[k++]=hmix(213,(cc|I_t1[cc]<<8)&0xffff); h[k++]=hmix(214,cc|I_t1[cc]<<8); h[k++]=hmix(215,(d|(U32)I_t2[d]<<16)&0xffffff); } }
      if (GF&65536) { U32 i0=(c4>>8&0xffffff)*2654435761u>>16, i1=(c4&0xffffff)*2654435761u>>16;
        I_t3[i0]=(U16)(I_t3[i0]<<8|(c4&255)); h[k++]=hmix(216,hmix(c4&0xffffff,I_t3[i1])); }
      if (GF&131072) h[k++]=hmix(217,hmix(I_w[B_prev>>16],h[5]));
      if (!G_x.empty()) { int dx=G_x[g1]-G_x[g2], dy=G_y[g1]-G_y[g2];
        XY_dir = (dx==0 && dy==0) ? 8 : (int)((atan2((double)dy,(double)dx)+3.14159265358979)*(8/6.28318530717959))%8;
        XY_px += ((G_x[g1]<<12) - XY_px) >> 6; XY_py += ((G_y[g1]<<12) - XY_py) >> 6; XY_place = (XY_px>>13)*8 + (XY_py>>13); }
      if (GF&2097152) { h[k++]=hmix(224,hmix(XY_dir*64+XY_place, c4&0xffff)); h[k++]=hmix(225,hmix(XY_place, h[5])); }
      if (FL_on) { fl_step(g2, g1); FL_p1 = fl_place(g1); FL_p2 = fl_place(g2); }
      if (GF&4194304) { h[k++]=hmix(226,hmix(FL_p1*4096+FL_p2, c4&0xffff)); h[k++]=hmix(227,hmix(FL_p1, h[5])); }
      static int ctxmode = getenv("CTXMODE") ? atoi(getenv("CTXMODE")) : 4;   // which cell key (imm_mapkey modes)
      if (GF&8388608) { h[k++]=hmix(228,imm_mapkey(ctxmode)^(c4&255)*0x9E3779B1u); h[k++]=hmix(229,hmix(g1*8192+g2, c4&0xffff)); }
      if (GF&16777216) { long jj=G_pos-G_len[g1]; int b2= jj>=0? G_hist[jj&65535]:0; long jk=jj-G_len[b2]; int b3= jk>=0? G_hist[jk&65535]:0;
        h[k++]=hmix(230,hmix(b2*8192+b3, c4&0xffff)); h[k++]=hmix(231,hmix(imm_mapkey(2), c4&0xffff)); }
      if (GF&33554432) { static std::vector<long> last; if (last.empty()) last.assign(G_len.size(),-1);
        auto lb=[&](int g){ long d= last[g]<0 ? 0 : G_pos-last[g]; int b=0; while (d>1 && b<24) { d>>=1; ++b; } return last[g]<0 ? 25 : b; };
        int r1=lb(g1), r2=lb(g2); last[g1]=G_pos;
        h[k++]=hmix(232,hmix(g1*32+r1, c4&0xffff)); h[k++]=hmix(233,hmix(r1*32+r2, imm_mapkey(2)^(c4&255)*0x9E3779B1u)); }
      if (GF&67108864) { static std::vector<int> nx2(1<<22,0); static int pg1=0,pg2=0;
        static int nxv = getenv("NXVOTE") ? atoi(getenv("NXVOTE")) : 0;   // round 34: majority vote (low 13 bits cell, high bits count)
        if (g1!=pg1 || g2!=pg2) { int& e=nx2[((U32)(g2*8192+g3)*0x9E3779B1u)>>10];
          if (!nxv) e=g1; else { int c=e>>13; if ((e&8191)==g1) { if (c<nxv) e+=8192; } else if (c>0) e-=8192; else e=g1; } pg1=g1; pg2=g2; }
        int p2=nx2[((U32)(g1*8192+g2)*0x9E3779B1u)>>10]&8191; G_p2=p2;
        h[k++]=hmix(234,hmix(p2*8192+g1,G_state)); h[k++]=hmix(235,hmix(p2, c4&0xffff)); }
      if (GF&134217728) { static std::vector<int> nx3(1<<22,0); static int qg1=0,qg2=0,qg3=0;
        long j4=j3-G_len[g3]; int g4= j4>=0? G_hist[j4&65535]:0;
        static int nxv3 = getenv("NXVOTE") ? atoi(getenv("NXVOTE")) : 0;
        if (g1!=qg1 || g2!=qg2 || g3!=qg3) { int& e=nx3[(((U32)(g2*8192+g3)*0x9E3779B1u)^(U32)g4*0x85EBCA77u)>>10];
          if (!nxv3) e=g1; else { int c=e>>13; if ((e&8191)==g1) { if (c<nxv3) e+=8192; } else if (c>0) e-=8192; else e=g1; } qg1=g1; qg2=g2; qg3=g3; }
        int p3=nx3[(((U32)(g1*8192+g2)*0x9E3779B1u)^(U32)g3*0x85EBCA77u)>>10]&8191; G_p3=p3;
        h[k++]=hmix(236,hmix(p3*8192+g1,G_state)); h[k++]=hmix(237,hmix(p3, c4&0xffff)); }
      if (GF&268435456) { static int o2=0,o3=0,hr=0; int hit=(o2==g1)+2*(o3==g1); hr=(hr*4+hit)&255;   // last four bytes' hits
        int ag= G_p2==G_p3; o2=G_p2; o3=G_p3;
        h[k++]=hmix(238,hmix(hr*2+ag, hmix(G_p3*8192+g1, c4&255))); h[k++]=hmix(239,hmix(hr*2+ag, c4&0xffff)); }
      if (GF&536870912) { static int e3=-1,hr=0; hr=(hr*2+(e3==(int)(c4&255)))&63;   // was the last expected byte right (6 bytes)
        e3=G_lastb[G_p3]; int e2=G_lastb[G_p2];
        h[k++]=hmix(240,(e3+1)*64+hr); h[k++]=hmix(241,hmix((e3+1)*512+(e2+1), c4&255)); }
      if (GF&1073741824) { static std::vector<int> nxr(1<<20,0); static U32 pk=0; long j4=j3-G_len[g3]; int g4= j4>=0? G_hist[j4&65535]:0;
        U32 wk=(((U32)(G_regB[g2]&127)*16384+(G_regB[g3]&127)*128+(G_regB[g4]&127))*0x9E3779B1u)>>12;   // the last three cells' regions, one step back
        if (wk!=pk) { nxr[wk]=g1; pk=wk; }
        int pr=nxr[(((U32)(G_regB[g1]&127)*16384+(G_regB[g2]&127)*128+(G_regB[g3]&127))*0x9E3779B1u)>>12];
        h[k++]=hmix(242,hmix(pr*8192+g1,G_state)); h[k++]=hmix(243,hmix(pr, c4&0xffff)); }
      if (GF&524288) { h[k++]=hmix(221,G_state); h[k++]=hmix(222,hmix(G_state,g2)); }
      if (GF&1048576) { G_next[g2]=g1; h[k++]=hmix(223,hmix(G_next[g1],G_state)); }
      if (GF&262144) { h[k++]=hmix(218,hmix(G_t16[g1]*256+G_t16[g2]*16+G_t16[g3], c4&255));
        h[k++]=hmix(219,hmix(G_t64[g1]*64+G_t64[g2], c4&0xffff)); h[k++]=hmix(220,hmix(G_t256[g1]*256+G_t256[g2], c4&255)); }
      if (GF&16384) { h[k++]=hmix(211,c4&0xffffff00); h[k++]=hmix(212,c4&0x00ff00ff); }
    }
    if (NGENE) {
      if ((byte|32)>=97 && (byte|32)<=122) W_cur=h[5]; else if (W_cur) { W_prev=W_cur; W_cur=0; }
      if (!(GF&16)) struct_update(byte);
      if (!(GF&64)) { U32 wd=0; if ((byte|32)>=97 && (byte|32)<=122) msw_cur=h[5]; else if (msw_cur) { wd=msw_cur; msw_cur=0; } multiscale_update(byte, wd); }
      int g1=G_hist[G_pos&65535]; long j2=G_pos-G_len[g1]; int g2= j2>=0? G_hist[j2&65535]:0;
      long j3=j2-G_len[g2]; int g3= j3>=0? G_hist[j3&65535]:0;
      for (int q=0;q<NGENE;++q) {
        int f=GENE_F[q], bm=GENE_B[q]; U32 hh=0x9E3779B1u*(q+1);
        for (int bi=0;bi<8;++bi) if (bm>>bi&1) { int cb = bi<4 ? (c4>>(8*bi))&255 : (c8>>(8*(bi-4)))&255; if (f&2048) cb=orb(cb); hh=hmix(hh, cb+256*bi); }
        if (f&1) hh=hmix(hh, 1000+g1); if (f&2) hh=hmix(hh, 20000+g2);
        if (f&4) hh=hmix(hh, 40000+G_reg[g1]); if (f&8) hh=hmix(hh, 41000+G_reg[g2]); if (f&16) hh=hmix(hh, 42000+G_reg[g3]);
        if (f&32) hh=hmix(hh, 43000+G_regB[g1]); if (f&64) hh=hmix(hh, 45000+G_regB[g2]);
        if (f&128) hh=hmix(hh, W_cur); if (f&256) hh=hmix(hh, W_prev^0x55555555u);
        if (f&512) hh=hmix(hh, 50000+S_state); if (f&1024) hh=hmix(hh, 51000+G_topic);
        if (f&4096) hh=hmix(hh, 52000+S_link*4+(S_tmpl>3?3:S_tmpl));
        if (f&8192) hh=hmix(hh, MS_sentw ^ 0x3C6EF372u); if (f&16384) hh=hmix(hh, MS_sec ^ 0xA54FF53Au);
        h[1+q]=hh;
      }
    }
    for (int i=1;i<N;++i) cp[i]=t[h[i]]+1;
    c0=1;
    bcount=0;
  }
  if (bcount==4) {
    for (int i=1;i<N;++i) cp[i]=t[h[i]+c0]+1;
  }
  else if (bcount>0) {
    int j=y+1<<(bcount&3)-1;
    for (int i=1;i<N;++i) cp[i]+=j;
  }
  cp[0]=t0+h[0]+c0;

  // predict
  int len=mm.p(y, m);
  if (mm2) mm2->p(y, m);
  if (rm) rm->p(y, m);
  if (im) im->p(y, m, pr);
  if (im2) im2->p(y, m, pr);
  if (im3) im3->p(y, m, pr);
  if (csm) { ImmuneModel* v[3]={im,im2,im3}; int combo=0, st=0;
    for (int i=0;i<3;++i) { int s3 = v[i] ? v[i]->pb+1 : 0; combo=combo*3+s3; if (s3 && v[i]->ps>st) st=v[i]->ps; }
    int cx = combo ? ((combo*16+(st>15?15:st))*8+bcount)*(immcons==2?256:1)+(immcons==2?(int)(c4&255):0) : 0;
    m.add(stretch(csm->p(y, cx))); }
  int order=0;
  if (len==0) {
    if (*cp[4]) ++order;
    if (*cp[3]) ++order;
    if (*cp[2]) ++order;
    if (*cp[1]) ++order;
  }
  else order=5+(len>=8)+(len>=12)+(len>=16)+(len>=32);
  for (int i=0;i<N;++i) m.add(stretch(sm[i].p(y, *cp[i])));
  if (MX_BIAS) m.add(256);
  int mx_r1 = GF ? G_reg[G_hist[G_pos&65535]] : 0, mx_s = (GF&16) ? S_state : 0;
  { static int msel = getenv("MIXSEL") ? atoi(getenv("MIXSEL")) : 0;   // 0 topic x struct, 1 struct only, 2 topic only
    int tp = (GF&8) && msel!=1 ? G_topic : 0, sp = (GF&16) && msel!=2 ? (regk ? REG_cur%10 : S_state) : 0;
    m.set(order+10*(h[0]>>13) + 80*(tp + ((GF&8) ? 64 : 1)*sp));
    static int mix2 = getenv("MIX2") ? 1 : 0;
    if (mix2) { m.setk(1, c0); m.setk(2, (c4&255)*8 + bcount); m.setk(3, ((len>15?15:len)*2 + (order>0))*10 + ((GF&16) ? S_state : 0));
      // round 28: MXFINAL=3 final layer by the current x previous cell's coarse regions, 4 = that x page structure
      int fr = 0; if (MX_FINAL>=3 && !G_len.empty()) { int f1=G_hist[G_pos&65535]; long fj=G_pos-G_len[f1]; int f2= fj>=0 ? G_hist[fj&65535] : 0; fr=(G_reg[f1]&15)*16+(G_reg[f2]&15); }
      m.setfinal(MX_FINAL==1 ? c0 : MX_FINAL==2 ? c0 + 256*mx_r1 : MX_FINAL==3 ? c0 + 256*fr : MX_FINAL==4 ? c0 + 256*(fr*10 + ((GF&16) ? S_state : 0)) : c0 + 256*((GF&16) ? S_state : 0));
      if (MX_SETS) { int kk=4;
        if ((MX_SETS&1) && kk<8) m.setk(kk++, mx_r1*10 + mx_s);
        if ((MX_SETS&2) && kk<8) m.setk(kk++, GF ? G_regB[G_hist[G_pos&65535]]&127 : 0);
        if ((MX_SETS&4) && kk<8) m.setk(kk++, G_topic&63);
        if ((MX_SETS&8) && kk<8) m.setk(kk++, (MS_wc>7?7:MS_wc)*10 + mx_s);
        if ((MX_SETS&16) && kk<8) m.setk(kk++, ((c4>>8)&255)*8 + bcount);
        if ((MX_SETS&32) && kk<8) m.setk(kk++, D_q*10 + mx_s);
        if ((MX_SETS&64) && kk<8) m.setk(kk++, G_t64.empty() ? 0 : G_t64[G_hist[G_pos&65535]]);
        if ((MX_SETS&128) && kk<8) m.setk(kk++, G_t256.empty() ? 0 : G_t256[G_hist[G_pos&65535]]);
        if ((MX_SETS&256) && kk<8) m.setk(kk++, XY_place);
        if ((MX_SETS&512) && kk<8) { int c=G_hist[G_pos&65535]; m.setblend(kk++, G_regB[c]&127, G_rb2[c], G_w2[c]); }
        if ((MX_SETS&1024) && kk<8) m.setk(kk++, FL_p1);
        if ((MX_SETS&2048) && kk<8) m.setk(kk++, G_len.empty() ? 0 : G_hist[G_pos&65535]&4095);
        if ((MX_SETS&4096) && kk<8) m.setk(kk++, (int)(((c4&0xffff)*2654435761u)>>20));
        if ((MX_SETS&8192) && kk<8) { int q1 = G_len.empty() ? 0 : G_hist[G_pos&65535]; long qj = G_pos-(G_len.empty()?0:G_len[q1]);
          int q2 = (qj>=0 && !G_len.empty()) ? G_hist[qj&65535] : 0; m.setk(kk++, q1*16 + (G_reg.empty() ? 0 : G_reg[q2]&15)); }
        if ((MX_SETS&16384) && kk<8) m.setk(kk++, (int)((((c4&0xffff)*2654435761u)>>8)%20480));
        if ((MX_SETS&131072) && kk<8) m.setk(kk++, G_len.empty() ? 0 : (G_regB[G_p2]&127)*16 + (G_reg[G_hist[G_pos&65535]]&15));
        if (MX_SETS&(32768|65536)) { int q1 = G_len.empty() ? 0 : G_hist[G_pos&65535]; long qj = G_pos-(G_len.empty()?0:G_len[q1]);
          int q2 = (qj>=0 && !G_len.empty()) ? G_hist[qj&65535] : 0;
          if ((MX_SETS&32768) && kk<8) m.setk(kk++, (q1%1280)*128 + (G_regB.empty() ? 0 : G_regB[q2]&127));
          if ((MX_SETS&65536) && kk<8) m.setk(kk++, (int)(((U32)(q1*8192+q2)*2654435761u)>>16)); } } } }
  pr=m.p();
  pr=pr+3*a1.pp(y, pr, c0)>>2;
  pr=pr+3*a2.pp(y, pr, c0^h[0]>>2)>>2;
  if (MX_APM>=4) {   // round 19: a final APM keyed on map cells (4: hash of the last three cells, 5: the current cell), averaged in
    // round 20: 6 = current cell x the previous cell's coarse region, 7 = the current cell's fine region, 8 = cell APM and
    // fine-region APM both, 9 = cell APM given 3/4 of the weight
    // round 21: 10 = control, the last two bytes hashed into 20,480 keys (same size as 6, no map); 11 = current cell x
    // previous cell, hashed into 32,768; 12 = current cell x previous cell's fine region, hashed into 32,768
    static APM* a4 = new APM(MX_APM==13 || MX_APM==14 ? 1 : MX_APM>=10 || MX_APM==6 ? 1<<23 : 1<<20); static APM* a5 = MX_APM==8 ? new APM(128*256) : 0;
    int cg1 = G_len.empty() ? 0 : G_hist[G_pos&65535]; long cj = G_pos-(G_len.empty()?0:G_len[cg1]); int cg2 = (cj>=0 && !G_len.empty()) ? G_hist[cj&65535] : 0;
    U32 kc = MX_APM==10 ? ((c4&0xffff)*2654435761u>>8)%20480 : MX_APM==11 ? ((U32)(cg1*8192+cg2)*2654435761u)>>17 :
      MX_APM==12 ? ((U32)(cg1*128+(G_regB[cg2]&127))*2654435761u)>>17 : MX_APM==4 ? imm_mapkey(2)>>20 : MX_APM==6 ? (U32)(cg1*16 + (G_reg[cg2]&15)) : MX_APM==7 ? (U32)(G_regB[cg1]&127) : (U32)(cg1&4095);
    // round 22: 13 = byte APM (as 10) + cell x previous region APM (as 6); 14 = control, two byte APMs (2 and 3 bytes)
    if (MX_APM==13 || MX_APM==14) { static APM* b4 = new APM(20480*256); static APM* b5 = new APM(20480*256);
      U32 kb = ((c4&0xffff)*2654435761u>>8)%20480, k2 = MX_APM==13 ? (U32)(cg1*16 + (G_reg[cg2]&15)) : ((c4&0xffffff)*2654435761u>>8)%20480;
      int p1=b4->pp(y, pr, kb*256 + c0), p2=b5->pp(y, pr, k2*256 + c0); pr=(2*pr+p1+p2+2)>>2; }
    else if (MX_APM==9) pr=pr+3*a4->pp(y, pr, kc*256 + c0)+2>>2;
    else if (MX_APM==8) { int p1=a4->pp(y, pr, kc*256 + c0), p2=a5->pp(y, pr, (G_regB[cg1]&127)*256 + c0); pr=(2*pr+p1+p2+2)>>2; }
    else pr=pr+a4->pp(y, pr, kc*256 + c0)+1>>1; }
  else if (MX_APM) pr=pr+a3.pp(y, pr, (MX_APM==1 ? mx_s : MX_APM==3 ? D_q : mx_r1)*256 + c0)+1>>1;
}

//////////////////////////// Encoder ////////////////////////////

// An Encoder does arithmetic encoding.  Methods:
// Encoder(COMPRESS, f) creates encoder for compression to archive f, which
//     must be open past any header for writing in binary mode.
// Encoder(DECOMPRESS, f) creates encoder for decompression from archive f,
//     which must be open past any header for reading in binary mode.
// code(i) in COMPRESS mode compresses bit i (0 or 1) to file f.
// code() in DECOMPRESS mode returns the next decompressed bit from file f.
// compress(c) in COMPRESS mode compresses one byte.
// decompress() in DECOMPRESS mode decompresses and returns one byte.
// flush() should be called exactly once after compression is done and
//     before closing f.  It does nothing in DECOMPRESS mode.

typedef enum {COMPRESS, DECOMPRESS} Mode;
class Encoder {
private:
  Predictor predictor;
  const Mode mode;       // Compress or decompress?
  FILE* archive;         // Compressed data file
  U32 x1, x2;            // Range, initially [0, 1), scaled by 2^32
  U32 x;                 // Decompress mode: last 4 input bytes of archive

  // Compress bit y or return decompressed bit
  FILE* dump; double dcost;  // env DUMPCOST (compress only): per byte, uint16 cost in 1/256 bit and uint16 map cell
  int code(int y=0) {
    int p=predictor.p();
    assert(p>=0 && p<4096);
    p+=p<2048;
    if (dump) dcost -= log2(y ? p/4096.0 : 1-p/4096.0);
    U32 xmid=x1 + (x2-x1>>12)*p + ((x2-x1&0xfff)*p>>12);
    assert(xmid>=x1 && xmid<x2);
    if (mode==DECOMPRESS) y=x<=xmid;
    y ? (x2=xmid) : (x1=xmid+1);
    predictor.update(y);
    while (((x1^x2)&0xff000000)==0) {  // pass equal leading bytes of range
      if (mode==COMPRESS) putc(x2>>24, archive);
      x1<<=8;
      x2=(x2<<8)+255;
      if (mode==DECOMPRESS) x=(x<<8)+(getc(archive)&255);  // EOF is OK
    }
    return y;
  }

public:
  Encoder(Mode m, FILE* f);
  void flush();  // call this when compression is finished

  // Pretraining: run bits through the model only (no coding); same on both sides
  void prime(int c) { for (int i=7; i>=0; --i) predictor.update((c>>i)&1); }

  // Compress one byte
  void compress(int c) {
    assert(mode==COMPRESS);
    for (int i=7; i>=0; --i)
      code((c>>i)&1);
    if (dump) { int q=(int)(dcost*256+0.5); if (q>65535) q=65535; int cell = GF ? G_hist[G_pos&65535] : 0;
      putc(q&255,dump); putc(q>>8,dump); putc(cell&255,dump); putc(cell>>8,dump); dcost=0; }
  }

  // Decompress and return one byte
  int decompress() {
    int c=0;
    for (int i=0; i<8; ++i)
      c+=c+code();
    return c;
  }
};

Encoder::Encoder(Mode m, FILE* f):
    mode(m), archive(f), x1(0), x2(0xffffffff), x(0), dump(0), dcost(0) {
  if (mode==COMPRESS && getenv("DUMPCOST")) { dump=fopen(getenv("DUMPCOST"),"wb"); if (!dump) quit("DUMPCOST"); }
  if (mode==DECOMPRESS) {  // x = first 4 bytes of archive
    for (int i=0; i<4; ++i)
      x=(x<<8)+(getc(archive)&255);
  }
}

void Encoder::flush() {
  if (mode==COMPRESS)
    putc(x1>>24, archive);  // Flush first unequal byte of range
}


//////////////////////////// User Interface ////////////////////////////

int main(int argc, char **argv) {

  // Check arguments
  if (argc!=4 || !isdigit(argv[1][0]) && argv[1][0]!='d') {
    printf(
      "lpaq1 file compressor (C) 2007, Matt Mahoney\n"
      "lpaq1g: lpaq1 + GridZip grid/spatial/symmetry contexts (GRIDMAP, GFLAGS env)\n"
      "Licensed under GPL, http://www.gnu.org/copyleft/gpl.html\n"
      "\n"
      "To compress:   lpaq1 N input output  (N=0..9, uses 3+3*2^N MB)\n"
      "To decompress: lpaq1 d input output  (needs same memory)\n");
    return 1;
  }

  // Get start time
  clock_t start=clock();

  // Open input file
  FILE *in=fopen(argv[2], "rb");
  if (!in) perror(argv[2]), exit(1);
  FILE *out=0;

  // Compress
  if (isdigit(argv[1][0])) {
    MEM=1<<(argv[1][0]-'0'+20);

    // Encode header: version 1, memory option, file size
    fseek(in, 0, SEEK_END);
    long size=ftell(in);
    if (size<0 || size>=0x7FFFFFFF) quit("input file too big");
    fseek(in, 0, SEEK_SET);
    out=fopen(argv[3], "wb");
    if (!out) perror(argv[3]), exit(1);
    fprintf(out, "pQ%c%c%c%c%c%c", 1, argv[1][0], 
      size>>24, size>>16, size>>8, size);

    // Compress
    Encoder e(COMPRESS, out);
    int c;
    if (getenv("PRIMER")) { FILE* pf=fopen(getenv("PRIMER"),"rb"); if (!pf) quit("no primer"); while ((c=getc(pf))!=EOF) e.prime(c); fclose(pf); }
    while ((c=getc(in))!=EOF)
      e.compress(c);
    e.flush();
  }

  // Decompress
  else {

    // Check header version, get memory option, file size
    if (getc(in)!='p' || getc(in)!='Q' || getc(in)!=1)
      quit("Not a lpaq1 file");
    MEM=getc(in);
    if (MEM<'0' || MEM>'9') quit("Bad memory option (not 0..9)");
    MEM=1<<(MEM-'0'+20);
    long size=getc(in)<<24;
    size|=getc(in)<<16;
    size|=getc(in)<<8;
    size|=getc(in);
    if (size<0) quit("Bad file size");

    // Decompress
    out=fopen(argv[3], "wb");
    if (!out) perror(argv[3]), exit(1);
    Encoder e(DECOMPRESS, in);
    if (getenv("PRIMER")) { FILE* pf=fopen(getenv("PRIMER"),"rb"); if (!pf) quit("no primer"); int c; while ((c=getc(pf))!=EOF) e.prime(c); fclose(pf); }
    while (size-->0)
      putc(e.decompress(), out);
  }

  // Report result
  assert(in);
  assert(out);
  printf("%ld -> %ld in %1.2f sec. using %d MB memory.\n", 
    ftell(in), ftell(out), double(clock()-start)/CLOCKS_PER_SEC, 
    3+(MEM>>20)*3);

  return 0;
}


