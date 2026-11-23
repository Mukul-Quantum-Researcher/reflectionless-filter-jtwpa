"""Check the supplied surrogate network functions on a sampled frequency grid.

This checks the external models only, not harmonic-balance convergence or gain.
"""
from pathlib import Path
import json
import numpy as np

def circulator(f):
    x=(f-np.sqrt(4e9*10e9))/3e9
    window=np.where(abs(x)<=1,1,np.exp(-((abs(x)-1)*2.2)**2))
    return 10**(-(20*window+1*(1-window))/20)*np.exp(1j*2*np.pi*f*.35e-9)

def filter_s(f):
    envelope=np.sqrt(1/(1+(3.2e9/np.maximum(f,1))**8)/(1+(f/11.8e9)**8))
    s21=10**(-.75/20)*envelope*np.exp(-1j*2*np.pi*f*80e-12)
    s11=10**(-28/20)*np.exp(1j*(2*np.pi*f*35e-12+.3))
    return np.array([[s11,s21],[s21,s11*np.exp(.2j)]])

frequencies=np.linspace(1e9,14e9,400)
min_real=[];max_norm=0
for f in frequencies:
    sf=filter_s(f);gc=circulator(f)
    gfc=sf[0,0]+sf[0,1]*sf[1,0]*gc/(1-sf[1,1]*gc)
    z=50*(1+np.array([gc,sf[0,0],gfc]))/(1-np.array([gc,sf[0,0],gfc]))
    min_real.append(float(z.real.min()))
    max_norm=max(max_norm,float(np.linalg.svd(sf,compute_uv=False).max()))
assert min(min_real)>0 and max_norm<=1+1e-12
result={'sample_points':400,'frequency_GHz':[1,14],'minimum_real_impedance_ohm':min(min_real),'maximum_filter_singular_value':max_norm,'scope':'Sampled surrogate-network check only; no fresh Julia harmonic-balance run.'}
root=Path(__file__).resolve().parents[1];(root/'results').mkdir(exist_ok=True)
(root/'results/termination-check.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps(result,indent=2))
