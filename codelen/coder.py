"""Exact integer arithmetic coder (the classic low/high/pending-bits scheme), written fresh for codelen.

Python integers, 62-bit state, frequency totals up to 2**44. Encoder and decoder must see the same (cum_lo, cum_hi,
total) for every symbol; codelen builds those from integer tables, so the round trip is exact.
"""
PREC = 62
TOP = (1 << PREC) - 1
HALF = 1 << (PREC - 1)
QTR = 1 << (PREC - 2)
MAX_TOTAL = 1 << 44


class Encoder:
    def __init__(self):
        self.low, self.high, self.pending, self.bits = 0, TOP, 0, []

    def _emit(self, b):
        self.bits.append(b)
        if self.pending:
            self.bits.extend([1 - b] * self.pending); self.pending = 0

    def encode(self, cum_lo, cum_hi, total):
        assert 0 <= cum_lo < cum_hi <= total <= MAX_TOTAL
        r = self.high - self.low + 1
        self.high = self.low + r * cum_hi // total - 1
        self.low = self.low + r * cum_lo // total
        while True:
            if self.high < HALF: self._emit(0)
            elif self.low >= HALF: self._emit(1); self.low -= HALF; self.high -= HALF
            elif self.low >= QTR and self.high < 3 * QTR: self.pending += 1; self.low -= QTR; self.high -= QTR
            else: break
            self.low <<= 1; self.high = (self.high << 1) | 1

    def finish(self):
        """Returns the bit list (list of 0/1). Call once."""
        self.pending += 1
        self._emit(0 if self.low < QTR else 1)
        return self.bits


class Decoder:
    def __init__(self, bits):
        self.bits, self.pos = bits, 0
        self.low, self.high, self.value = 0, TOP, 0
        for _ in range(PREC): self.value = (self.value << 1) | self._next()

    def _next(self):
        b = self.bits[self.pos] if self.pos < len(self.bits) else 0
        self.pos += 1; return b

    def target(self, total):
        """The cumulative count the next symbol's interval contains."""
        r = self.high - self.low + 1
        return ((self.value - self.low + 1) * total - 1) // r

    def consume(self, cum_lo, cum_hi, total):
        r = self.high - self.low + 1
        self.high = self.low + r * cum_hi // total - 1
        self.low = self.low + r * cum_lo // total
        while True:
            if self.high < HALF: pass
            elif self.low >= HALF: self.low -= HALF; self.high -= HALF; self.value -= HALF
            elif self.low >= QTR and self.high < 3 * QTR: self.low -= QTR; self.high -= QTR; self.value -= QTR
            else: break
            self.low <<= 1; self.high = (self.high << 1) | 1; self.value = (self.value << 1) | self._next()


def pack(bits):
    """Bits to bytes (zero padded); len(result) * 8 >= len(bits)."""
    out = bytearray((len(bits) + 7) // 8)
    for i, b in enumerate(bits):
        if b: out[i >> 3] |= 0x80 >> (i & 7)
    return bytes(out)


def unpack(data):
    return [(data[i >> 3] >> (7 - (i & 7))) & 1 for i in range(len(data) * 8)]
