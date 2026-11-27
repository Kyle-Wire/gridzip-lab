import sys, json, numpy as np
sys.path.insert(0,'.'); sys.path.insert(0,'physics')
from physlib import guards as Q, terms as T
import flood_walls as FW
res={}
for n in (25,50,100):
  hits=0; tries=3
  for seed in range(tries):
    rng=np.random.default_rng(100+seed)
    base=[{"site":"S","hwm_id":str(i),"kind":"level","quality":"Good","type":"Seed line","waterbody":f"w{i%3}","at_building":str(int(rng.random()<0.6)),"model_dry":"0","h_model_m":f"{rng.uniform(0.2,2.5):.3f}","hand_m":f"{rng.uniform(0,4):.3f}","d_channel_m":f"{rng.uniform(5,600):.1f}","d_edge_m":f"{rng.uniform(500,5000):.1f}","bld_frac_30":f"{rng.uniform(0,0.3):.4f}","slope":f"{rng.uniform(0,0.05):.4f}","v_max_ms":"","spill_rise_m":f"{rng.uniform(0,0.4):.3f}"} for i in range(n)]
    rr=[dict(r,err_m=f"{-0.8*FW._f(r['d_channel_m'])/(FW._f(r['d_channel_m'])+150)+rng.normal(0,0.15):.3f}") for r in base]
    ctx,y,s,wb,_=FW.build(rr)
    folds=[(f"w{k}",np.where(wb!=f"S:w{k}")[0],np.where(wb==f"S:w{k}")[0]) for k in range(3)]
    terms=[t for t in T.FLOOD_WALL if all(k in ctx for k in t.needs)]
    r=Q.analyse(terms,ctx,y,FW.GRID,folds=folds,bootstrap_B=0,null_N=19,name='p')
    w=r['winner']['terms']; p=r['guards']['null_bits_saved']['p']; ok=('storage_fill' in w or any('d_channel' in t or 'd_ch' in t for t in w)) and p<=0.05
    hits+=ok; print(n,seed,w,round(r['bits_saved_per_reading'],3),p,flush=True)
  res[n]=hits/tries
print(res); json.dump(res,open('results/physics/flood_walls_a13_power.json','w'))
