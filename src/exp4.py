"""Transformer on GRID TOKENS vs on BYTES.

Same model shape, same steps, batch and context length (in tokens). The only
difference: the byte model reads raw bytes (256 symbols); the grid model
reads grid coordinates (256 + CELLS symbols), so each window covers ~3x more
text.

Train on enwik8[0:TRAIN], test on a held-out slice far later in the file.
Score: total bits to encode the test slice / test BYTES (bits per byte), so
both models are judged on the same thing. Grid map cost is added.
Model weights are NOT counted (a pre-trained predictor, like a pre-shared map).
"""
import json, math, os, sys, time
import numpy as np
import torch, torch.nn as nn, torch.nn.functional as F
from gridzip import learn_grid, apply_grid, map_bytes

CTX, BATCH, DM, LAYERS, HEADS = 128, 32, 128, 4, 4
STEPS = int(os.environ.get("STEPS", 4000))
torch.set_num_threads(int(os.environ.get("THREADS", 1)))
TRAIN, TEST_AT, TEST_N = 9_000_000, 95_000_000, 500_000


class Block(nn.Module):
    def __init__(s):
        super().__init__()
        s.ln1, s.ln2 = nn.LayerNorm(DM), nn.LayerNorm(DM)
        s.qkv, s.proj = nn.Linear(DM, 3 * DM), nn.Linear(DM, DM)
        s.mlp = nn.Sequential(nn.Linear(DM, 4 * DM), nn.GELU(), nn.Linear(4 * DM, DM))
    def forward(s, x):
        B, T, C = x.shape
        q, k, v = s.qkv(s.ln1(x)).view(B, T, 3, HEADS, C // HEADS).permute(2, 0, 3, 1, 4)
        a = F.scaled_dot_product_attention(q, k, v, is_causal=True)
        x = x + s.proj(a.transpose(1, 2).reshape(B, T, C))
        return x + s.mlp(s.ln2(x))


class GPT(nn.Module):
    def __init__(s, V):
        super().__init__()
        s.tok, s.pos = nn.Embedding(V, DM), nn.Embedding(CTX, DM)
        s.blocks = nn.Sequential(*[Block() for _ in range(LAYERS)])
        s.ln, s.head = nn.LayerNorm(DM), nn.Linear(DM, V)
    def forward(s, x):
        return s.head(s.ln(s.blocks(s.tok(x) + s.pos(torch.arange(x.shape[1])))))


@torch.no_grad()
def total_bits(model, data, V):
    model.eval()
    n = (len(data) - 1) // CTX
    tot = 0.0
    for i in range(0, n, 64):
        idx = range(i, min(n, i + 64))
        x = torch.stack([data[j * CTX:(j + 1) * CTX] for j in idx])
        y = torch.stack([data[j * CTX + 1:(j + 1) * CTX + 1] for j in idx])
        tot += F.cross_entropy(model(x).reshape(-1, V), y.reshape(-1), reduction="sum").item()
    model.train()
    covered = n * CTX
    return tot / math.log(2), covered


def run(kind, cells, seed, raw):
    train_b, test_b = raw[:TRAIN], raw[TEST_AT:TEST_AT + TEST_N]
    if kind == "grid":
        merges, _ = learn_grid(train_b[:2_000_000], cells)
        tr = apply_grid(train_b, merges); te = apply_grid(test_b, merges)
        V, mcost = 256 + len(merges), map_bytes(merges)
    else:
        tr = np.frombuffer(train_b, np.uint8).astype(np.int64)
        te = np.frombuffer(test_b, np.uint8).astype(np.int64)
        V, mcost = 256, 0.0
    tr, te = torch.tensor(tr), torch.tensor(te)
    torch.manual_seed(seed); g = torch.Generator().manual_seed(seed)
    model = GPT(V)
    opt = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=0.1)
    sched = torch.optim.lr_scheduler.LambdaLR(
        opt, lambda s: min(1, s / 100) * 0.5 * (1 + math.cos(math.pi * s / STEPS)))
    t0 = time.time(); curve = []
    tok_per_byte = len(te) / TEST_N
    for step in range(STEPS + 1):
        if step % (STEPS // 8) == 0:
            bits, cov = total_bits(model, te, V)
            bpb = bits / (cov / tok_per_byte)  # bits per test byte
            curve.append(dict(step=step, bpb=bpb))
            print(kind, cells, seed, curve[-1], f"{time.time()-t0:.0f}s", flush=True)
        if step == STEPS:
            break
        ix = torch.randint(0, len(tr) - CTX - 1, (BATCH,), generator=g)
        x = torch.stack([tr[i:i + CTX] for i in ix]); y = torch.stack([tr[i + 1:i + CTX + 1] for i in ix])
        loss = F.cross_entropy(model(x).reshape(-1, V), y.reshape(-1))
        opt.zero_grad(); loss.backward(); nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        opt.step(); sched.step()
    return dict(kind=kind, cells=cells, seed=seed, bpb=curve[-1]["bpb"], map_bytes=mcost,
                tokens_per_byte=tok_per_byte, params=sum(p.numel() for p in model.parameters()),
                curve=curve, secs=round(time.time() - t0))


if __name__ == "__main__":
    raw = open(sys.argv[1], "rb").read()
    out = sys.argv[2]
    for spec in sys.argv[3].split(","):        # e.g. bytes:0,grid:1024
        kind, cells = spec.split(":")
        for seed in [int(s) for s in sys.argv[4].split(",")]:
            r = run(kind, int(cells), seed, raw)
            open(out, "a").write(json.dumps(r) + "\n")
