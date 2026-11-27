import json, numpy as np
from pyproj import Transformer
from skimage.morphology import reconstruction
D='/tmp/claude-0/lg/'; c=json.load(open(D+'crop.json')); nx,ny=c['cells']; x0,y0=c['box_utm'][0],c['box_utm'][1]
rd=lambda n: np.fromfile(D+n+'.f32','<f4').reshape(ny,nx)[::-1]  # row0 = south after flip? ESRI: first row is north
g=np.fromfile(D+'ground.f32','<f4').reshape(ny,nx); h=np.fromfile(D+'max_depth.f32','<f4').reshape(ny,nx); w=np.fromfile(D+'max_wse.f32','<f4').reshape(ny,nx)
g=np.where(g<-1000,np.nan,g)
print('ground range',np.nanmin(g),np.nanmax(g),'nan',np.isnan(g).mean())
# channel seeds: deep water in the free-edge run
seed_mask=(h>1.5)&(g<14.3)
print('seed share',seed_mask.mean())
G=np.where(np.isnan(g),np.nanmax(g)+10,g).astype(np.float64)
seed=np.where(seed_mask,G,G.max()); fill=reconstruction(seed,G,method='erosion')
T=Transformer.from_crs(4326,c['epsg'],always_xy=True)
out=[]
for m in c['marks']:
    X,Y=T.transform(m['lon'],m['lat']); i=int(X-x0); j=int(c['box_utm'][3]-Y)
    win=(slice(j-2,j+3),slice(i-2,i+3))
    gm=np.nanmin(g[win]); bar=np.min(fill[win])
    out.append(dict(hwm_id=m['hwm_id'],ground=round(float(gm),2),spill_from_channel=round(float(bar),2),obs=m['level_m'],obs_from_hag=round(float(gm+m['hag_m']),2),wse_free=round(float(np.nanmax(np.where(h[win]>0.01,w[win],np.nan))) if (h[win]>0.01).any() else None,2) if (h[win]>0.01).any() else None))
for o in out: print(o)
np.save('/tmp/claude-0/fw/lumber_fill.npy',fill.astype(np.float32)); json.dump(out,open('/tmp/claude-0/fw/lumber_barrier.json','w'),indent=1)
