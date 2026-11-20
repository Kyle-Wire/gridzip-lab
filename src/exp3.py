"""Kid -> PhD: does the order of training data change how well a small
transformer compresses specialist text?

Every condition gets the same model, the same number of training steps and
the same learning-rate schedule; only WHICH data each step sees changes.
Score = bits per byte on held-out text (= size an arithmetic coder driven
by the model would produce, not counting the model itself).
"""
import json, math, os, sys, time
import numpy as np
import torch, torch.nn as nn, torch.nn.functional as F

D = os.path.join(os.path.dirname(__file__), "..", "data")
CTX, BATCH, DM, LAYERS, HEADS = 128, 32, 128, 4, 4
STEPS = int(os.environ.get("STEPS", 4000))
NARROW_SHARE = float(os.environ.get("SHARE", 0.2))
torch.set_num_threads(int(os.environ.get("THREADS", 3)))

def load(name, n=None):
    b = open(f"{D}/{name}", "rb").read()
    return torch.tensor(np.frombuffer(b[:n] if n else b, dtype=np.uint8).astype(np.int64))

broad, narrow = load("train.txt"), load("docs_train.txt")
eval_sets = {"docs_test": load("docs_test.txt"), "wiki_valid": load("valid.txt", 200_000)}

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
    def __init__(s):
        super().__init__()
        s.tok, s.pos = nn.Embedding(256, DM), nn.Embedding(CTX, DM)
        s.blocks = nn.Sequential(*[Block() for _ in range(LAYERS)])
        s.ln, s.head = nn.LayerNorm(DM), nn.Linear(DM, 256)
    def forward(s, x):
        h = s.tok(x) + s.pos(torch.arange(x.shape[1]))
        return s.head(s.ln(s.blocks(h)))

def batch(src, g):
    ix = torch.randint(0, len(src) - CTX - 1, (BATCH,), generator=g)
    x = torch.stack([src[i:i + CTX] for i in ix])
    y = torch.stack([src[i + 1:i + CTX + 1] for i in ix])
    return x, y

@torch.no_grad()
def bpb(model, data):
    model.eval()
    n = (len(data) - 1) // CTX
    tot = 0.0
    for i in range(0, n, 64):
        idx = range(i, min(n, i + 64))
        x = torch.stack([data[j * CTX:(j + 1) * CTX] for j in idx])
        y = torch.stack([data[j * CTX + 1:(j + 1) * CTX + 1] for j in idx])
        tot += F.cross_entropy(model(x).reshape(-1, 256), y.reshape(-1), reduction="sum").item()
    model.train()
    return tot / (n * CTX) / math.log(2)

def source_for(cond, step, g):
    cut = int(STEPS * (1 - NARROW_SHARE))
    if cond == "broad_only": return broad
    if cond == "narrow_only": return narrow
    if cond == "broad_then_narrow": return broad if step < cut else narrow
    if cond == "narrow_then_broad": return narrow if step < STEPS - cut else broad
    if cond == "mixed": return narrow if torch.rand(1, generator=g).item() < NARROW_SHARE else broad
    raise ValueError(cond)

def run(cond, seed):
    torch.manual_seed(seed); g = torch.Generator().manual_seed(seed)
    model = GPT()
    opt = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=0.1)
    sched = torch.optim.lr_scheduler.LambdaLR(opt, lambda s: min(1, s / 100) * 0.5 * (1 + math.cos(math.pi * s / STEPS)))
    curve, t0 = [], time.time()
    for step in range(STEPS + 1):
        if step % (STEPS // 10) == 0:
            curve.append(dict(step=step, **{k: bpb(model, v) for k, v in eval_sets.items()}))
            print(cond, seed, curve[-1], f"{time.time()-t0:.0f}s", flush=True)
        if step == STEPS: break
        x, y = batch(source_for(cond, step, g), g)
        loss = F.cross_entropy(model(x).reshape(-1, 256), y.reshape(-1))
        opt.zero_grad(); loss.backward(); nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        opt.step(); sched.step()
    torch.save(model.state_dict(), os.path.join(os.path.dirname(__file__), f"model_{cond}_{seed}.pt"))
    return dict(cond=cond, seed=seed, final=curve[-1],
                best_docs=min(c["docs_test"] for c in curve), curve=curve,
                params=sum(p.numel() for p in model.parameters()))

if __name__ == "__main__":
    conds = sys.argv[1].split(",")
    seeds = [int(s) for s in sys.argv[2].split(",")]
    out = os.path.join(os.path.dirname(__file__), sys.argv[3] if len(sys.argv) > 3 else "exp3_results.jsonl")
    for seed in seeds:
        for c in conds:
            r = run(c, seed)
            open(out, "a").write(json.dumps(r) + "\n")
