import json, csv, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sp={o['hwm_id']:o for o in json.load(open('results/physics/flood_walls_a13_lumber_spill.json'))}
ctx={x['hwm_id']:x for x in csv.DictReader(open('/home/user/climate-home-inspection/runs-summary/walls/lumber-terms.csv'))}
d=json.load(open('/home/user/climate-home-inspection/docs/validation/walls-tuning/lumber-sc-2016-edge-thinfilm-A13.json'))['depth']
fig,ax=plt.subplots(1,3,figsize=(16,5.2))
wet=[m for m in d if m.get('err_m') is not None]
xs=[float(ctx[str(m['hwm_id'])]['d_channel_m']) for m in wet]; ys=[m['err_m'] for m in wet]
ax[0].scatter(xs,ys,s=60,color='#1f6fb5',zorder=3); ax[0].axhspan(-0.3,0.3,color='#3a9a5b',alpha=.12,label='within 0.30 m')
ax[0].axhline(np.mean(ys),color='#d9822b',ls='--',label=f'whole-town shortfall {np.mean(ys):+.2f} m'); ax[0].axhline(0,color='k',lw=.8)
ax[0].set_xlabel('distance from the river channel (m)'); ax[0].set_ylabel('model minus mark (m)'); ax[0].legend(fontsize=8)
ax[0].set_title('Lumber, fix on: miss vs distance from the river\n(farther marks miss more; 10 marks cannot prove it)',fontsize=10,loc='left')
lv=[m for m in d if m.get('kind')=='level']
ids=[str(m['hwm_id']) for m in lv]; X=np.arange(len(ids))
ax[1].bar(X,[sp[i]['spill_from_channel']-15 for i in ids],0.5,bottom=15,color='#c9ced6',label='ridge the water must cross from the river')
ax[1].scatter(X,[m.get('model_m') or np.nan for m in lv],marker='s',color='#1f6fb5',s=50,zorder=3,label='model water level')
ax[1].scatter(X,[next(c['level_m'] for c in json.load(open('/tmp/claude-0/lg/crop.json'))['marks'] if c['hwm_id']==i) for i in ids],marker='*',color='k',s=120,zorder=3,label='high-water mark')
ax[1].axhline(16.62,color='#d9822b',ls='--',label='held downstream edge, 16.62 m'); ax[1].axhline(16.89,color='#7a3fb0',ls=':',label='Nichols gauge crest, 16.89 m')
for x,m in zip(X,lv):
    if m.get('model_m') is None: ax[1].text(x,15.1,'dry',ha='center',fontsize=8,color='#1f6fb5')
ax[1].set_xticks(X); ax[1].set_xticklabels([f'mark {i}' for i in ids],rotation=45,fontsize=8); ax[1].set_ylim(15,18.3); ax[1].set_ylabel('m above sea level (NAVD88)'); ax[1].legend(fontsize=7,loc='upper left',ncol=2)
ax[1].set_title('Lumber: no ridge cuts the town off (grey tops < edge),\nyet the model sits BELOW its own downstream edge',fontsize=10,loc='left')
ax[2].bar(['25','50','100'],[0,3,3],color=['#c44','#3a9a5b','#3a9a5b'])
ax[2].set_ylabel('planted tests where the cause was named (of 3)'); ax[2].set_xlabel('wet high-water marks in the search')
ax[2].text(0,0.2,'we have 15',ha='center',fontsize=9); ax[2].set_ylim(0,3.4)
ax[2].set_title('How many marks the compressor needs\nto name a 0.8 m distance effect',fontsize=10,loc='left')
fig.suptitle('Flood compressor, protocol 15 (exploratory): Lumber is a whole-town shortfall, not a blocked path',x=0.01,ha='left',fontsize=12)
fig.tight_layout(); fig.savefig('/mnt/project-files/flood-walls/p15-compressor.png',dpi=130)
