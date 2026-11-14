"""E14: deterministic Boussinesq column, mass fractions and constant sensible cp.
Instantaneous conservative convective adjustment; no VLE, wind or resolved velocity.
Inputs inherited from reproduced E13 reference.csv; all forcing scenarios assumed.
"""
import csv,json,platform,itertools
from pathlib import Path
import numpy as np
D=Path(__file__).resolve().parent
M=np.array([.016043,.030070,.0280134]); x0=np.array(json.loads((D/'prior_E13/results.json').read_text())['reference_composition']); w0=x0*M/(x0@M)
CP=3000.; RHO=560.; DAY=86400.; END=11.4*DAY
MODELS={r['model']:r for r in csv.DictReader((D/'prior_E13/reference.csv').open()) if r['model'] in ('constant_shift','component_shift')}
FORCINGS={'quiet':(0,0,0),'eclipse':(0,5,0),'rain5':(.005,0,0),'rain25':(.025,0,0),'rain25_eclipse':(.025,5,0),'rain25_cooling':(.025,0,1)}
# Rain replaces equal mass of mixed surface liquid (constant layer mass). 24h pulse. Rain xN=0.12 and CH4/(CH4+C2H6)=0.8; 92K. Assumptions.
xrain=np.array([.8*.88,.2*.88,.12]); wrain=xrain*M/(xrain@M); train=92.; rainrho=450.
def mole(w):
    n=w/M;return n/n.sum(axis=-1,keepdims=True)
def slopes(model):
    row=MODELS[model]; a=float(row['dM_plus_5pt'])/.05; b=float(row['dN_minus_2pt'])/(-.02); t=float(row['dT_minus_0p3'])/(-.3)
    # Linearise in mass-fraction coordinates, with ethane as remainder; local secant E13 slopes.
    def sol(w):
        x=mole(w);return a*(x[0]/(1-x[2])-x0[0]/(1-x0[2]))+b*(x[2]-x0[2])
    eps=1e-6; d=[]
    for i in (0,2):
        e=np.zeros(3);e[i]=eps;e[1]=-eps; d.append((sol(w0+e)-sol(w0-e))/(2*eps))
    return np.array([d[0],0,d[1],t])
def adjust(s,b):
    # pooling adjacent unstable blocks, upper density must not exceed lower
    r=s@b
    if np.all(r[:-1]<=r[1:]+1e-10):return 0
    blocks=[];count=0
    for i in range(len(s)):
        blocks.append([i,i+1,s[i].copy(),r[i]])
        while len(blocks)>1 and blocks[-2][3]>blocks[-1][3]+1e-10:
            l,h=blocks[-2:]; n1=l[1]-l[0]; n2=h[1]-h[0]; v=(l[2]*n1+h[2]*n2)/(n1+n2);blocks[-2:]=[[l[0],h[1],v,float(v@b)]];count+=1
    for lo,hi,v,r in blocks:s[lo:hi]=v
    return count

def simulate(depth,K,model,forcing,n=20,dt=360,record=False,conv=True):
    b=slopes(model); dz=depth/n; m=RHO*dz; lam=K*dt/dz**2; assert lam<=.5
    # state: three mass fractions, T-94. Heat ledger in J/m2.
    s=np.tile(np.r_[w0,0.],(n,1)); init=s.sum(axis=0)*m; ledger=np.zeros(4); rain,eclipse,cooling=FORCINGS[forcing]; mixsteps=0; maxblocks=0; maxcontrast=0.; hist=[]; prev=0
    def snapshot(t):
        x=mole(s[:,:3]);r=s@b; return dict(day=t/DAY,top_T=94+s[0,3],bottom_T=94+s[-1,3],top_CH4=x[0,0],bottom_CH4=x[-1,0],top_N2=x[0,2],bottom_N2=x[-1,2],density_contrast_bottom_minus_top=r[-1]-r[0],min_adjacent_density_step=float(np.diff(r).min()),top_minus_bottom_CH4=x[0,0]-x[-1,0])
    times={0,DAY,DAY+15120,DAY+30240,3*DAY,7*DAY,END}; hist.append(snapshot(0)); at_eclipse=None; peak=0.; cumulative_q=0.; cumulative_rain=0.
    for j in range(round(END/dt)):
        t=j*dt
        # Rain replaces mass via conservative discrete exchange with pre-exchange surface composition.
        if t<DAY and rain:
            dm=rain*rainrho*dt/DAY; assert dm<=m
            delta=dm*(np.r_[wrain,train-94]-s[0]);s[0]+=delta/m;ledger+=delta;cumulative_rain+=dm
        # eclipse and equal-energy recovery start day 1; optional extended heat loss 3 days after rain.
        q=(-eclipse if DAY<=t<DAY+15120 else eclipse if DAY+15120<=t<DAY+30240 else 0.)-(cooling if DAY<=t<4*DAY else 0.)
        delta=q*dt/CP;s[0,3]+=delta/m;ledger[3]+=delta;cumulative_q+=q*dt
        f=lam*np.diff(s,axis=0);s[:-1]+=f;s[1:]-=f
        k=adjust(s,b) if conv else 0;mixsteps+=bool(k);maxblocks=max(maxblocks,k)
        peak=max(peak,float((s[-1]-s[0])@b))
        if t+dt in times:hist.append(snapshot(t+dt))
        if t+dt==DAY+15120:at_eclipse=snapshot(t+dt)
    err=s.sum(axis=0)*m-init-ledger
    assert np.min(s[:,:3])>=-1e-12;assert np.max(abs(s[:,:3].sum(axis=1)-1))<1e-10
    assert np.max(abs(err[:3]))<1e-7;assert abs(err[3]*CP)<1e-4
    if conv:assert np.min(np.diff(s@b))>=-1.01e-10
    last=snapshot(END); post=at_eclipse or last
    result=dict(depth_m=depth,K_m2_s=K,model=model,forcing=forcing,n=n,dt_s=dt,steps=round(END/dt),mixsteps=mixsteps,max_merges_per_step=maxblocks,peak_contrast_kg_m3=peak,**last,eclipse_end_contrast=post['density_contrast_bottom_minus_top'],eclipse_end_top_T=post['top_T'],mass_error_kg_m2=float(max(abs(err[:3]))),heat_error_J_m2=float(abs(err[3]*CP)),net_heat_J_m2=cumulative_q,rain_mass_kg_m2=cumulative_rain)
    return result,hist,s

def write(name,rows):
    with (D/name).open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def main():
    results=[];series=[];profiles=[]
    for depth,K,model,forcing in itertools.product([1.,10.,100.],[1e-8,1e-6],MODELS,FORCINGS):
        r,h,s=simulate(depth,K,model,forcing);results.append(r);keys={k:r[k] for k in ('depth_m','K_m2_s','model','forcing')};series.extend([dict(**keys,**v) for v in h]);x=mole(s[:,:3]);profiles.extend([dict(**keys,z_m=(i+.5)*depth/len(s),T_K=94+v[3],xCH4=x[i,0],xC2H6=x[i,1],xN2=x[i,2],density_anomaly=float((v-np.r_[w0,0])@slopes(model))) for i,v in enumerate(s)])
    convergence=[]
    for dep,model,force in itertools.product([1.,10.],MODELS,['eclipse','rain25_eclipse']):
        r=next(r for r in results if r['depth_m']==dep and r['K_m2_s']==1e-6 and r['model']==model and r['forcing']==force)
        # timestep refinement at fixed grid; not spatial convergence
        f,_,_=simulate(dep,1e-6,model,force,n=20,dt=180)
        convergence.append(dict(depth_m=dep,model=model,forcing=force,coarse_contrast=r['density_contrast_bottom_minus_top'],fine_contrast=f['density_contrast_bottom_minus_top'],abs_difference=abs(r['density_contrast_bottom_minus_top']-f['density_contrast_bottom_minus_top'])))
    spatial=[]
    for dep,model in itertools.product([1.,10.,100.],MODELS):
        r=next(r for r in results if r["depth_m"]==dep and r["K_m2_s"]==1e-6 and r["model"]==model and r["forcing"]=="rain25_eclipse")
        f,_,_=simulate(dep,1e-6,model,"rain25_eclipse",n=40,dt=90)
        spatial.append(dict(depth_m=dep,model=model,coarse_contrast=r["density_contrast_bottom_minus_top"],fine_contrast=f["density_contrast_bottom_minus_top"],abs_difference=abs(r["density_contrast_bottom_minus_top"]-f["density_contrast_bottom_minus_top"])))
    # analytical and conservation controls
    b=slopes('constant_shift');s=np.tile(np.r_[w0,0.],(4,1));s[0,3]=-1;old=s.sum(axis=0);assert adjust(s,b)>0;assert np.allclose(s,old/4);assert np.allclose(s.sum(axis=0),old)
    s=np.tile(np.r_[w0,0.],(4,1));s[0,3]=1;old=s.copy();assert adjust(s,b)==0;assert np.array_equal(s,old)
    assert all(r['peak_contrast_kg_m3']<1e-10 and r['mixsteps']==0 for r in results if r['forcing']=='quiet')
    write('sweep.csv',results);write('timeseries.csv',series);write('profiles.csv',profiles);write('timestep_check.csv',convergence);write('spatial_check.csv',spatial)
    summary=dict(python=platform.python_version(),numpy=np.__version__,seed=None,main_cases=len(results),time_series_rows=len(series),profile_rows=len(profiles),timestep_cases=len(convergence),spatial_cases=len(spatial),max_spatial_contrast_difference=max(r["abs_difference"] for r in spatial),max_mass_residual=max(r['mass_error_kg_m2'] for r in results),max_heat_residual=max(r['heat_error_J_m2'] for r in results),max_timestep_contrast_difference=max(r['abs_difference'] for r in convergence),slopes={m:slopes(m).tolist() for m in MODELS},controls=['quiet columns unchanged','unstable cooled top pooled conservatively','stable warm top unchanged','nonnegative species and conserved budgets','stable final density ordering'],assumed_RHO=RHO,cp=CP,rain_x=xrain.tolist(),rain_T=train,rain_rho=rainrho,forcings=FORCINGS)
    (D/'results.json').write_text(json.dumps(summary,indent=2)+'\n'); print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
