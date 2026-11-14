"""E15: deterministic capillary-gravity waves and two-layer mixing energy.
All EEB-1 parameters are fictional reference or sensitivity inputs.
No wind climate, entrainment rate, wave growth or sediment erosion is predicted.
"""
import csv, itertools, json, math, platform
from pathlib import Path

D = Path(__file__).resolve().parent
G, RHO, AIR = 1.33, 560., 6.
WINDS = [0., .25, .5, 1., 2., 4.]
HOURS = [1., 4.2, 24.]
TOP = [.05, .1, .5]
BOTTOM = [.5, 2., 10.]
CONTRAST = [.5, 5., 15.]
DRAG = [.001, .0015, .003]
EFF = [.1, .4, 1.]
SIGMA = [.0178, .023, .0318]

def write_csv(name, rows):
    with (D/name).open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

def barrier(h1,h2,dr):
    return .5*G*h1*h2*dr

def wind_power(u,cd,m):
    tau = AIR*cd*u*u
    us = math.sqrt(tau/RHO)
    return tau,us,m*tau*us

def wave(k,sigma,depth=None):
    return math.sqrt((G*k+sigma/RHO*k**3)*(1. if depth is None else math.tanh(k*depth)))

def main():
    rows=[]
    for u,hours,h1,h2,dr,cd,m in itertools.product(WINDS,HOURS,TOP,BOTTOM,CONTRAST,DRAG,EFF):
        tau,us,p=wind_power(u,cd,m)
        need=barrier(h1,h2,dr)
        energy=p*hours*3600
        # Invert the energy barrier for maximum energetically accessible underlying depth;
        # not a prediction of actual entrainment.
        accessible=min(h2,2*energy/(G*h1*dr))
        threshold=(need/(hours*3600*m*(AIR*cd)**1.5/math.sqrt(RHO)))**(1/3)
        rows.append(dict(wind_m_s=u,hours=hours,top_m=h1,bottom_m=h2,delta_rho=dr,drag=cd,mixing_coefficient=m,stress_Pa=tau,liquid_friction_velocity=us,power_W_m2=p,required_J_m2=need,available_J_m2=energy,energy_ratio=energy/need,energy_sufficient=energy>=need,accessible_bottom_m=accessible,threshold_wind_m_s=threshold))
    waves=[]
    for sigma,lam,depth in itertools.product(SIGMA,[.005,.01,.03,.1,.5,1.,5.,10.],[.1,1.,10.,None]):
        k=2*math.pi/lam
        omega=wave(k,sigma,depth)
        waves.append(dict(sigma_N_m=sigma,wavelength_m=lam,depth_m=depth,phase_speed_m_s=omega/k,period_s=2*math.pi/omega))
    minima=[]
    for s in SIGMA:
        k=math.sqrt(G*RHO/s)
        c=wave(k,s)/k
        minima.append(dict(sigma_N_m=s,wavelength_m=2*math.pi/k,phase_speed_m_s=c,period_s=2*math.pi/wave(k,s)))
    # Validate two-layer PE via explicit vertical integration, using bottom as z=0.
    for h1,h2,dr in itertools.product(TOP,BOTTOM,CONTRAST):
        a,b=RHO-dr,RHO
        H=h1+h2
        mean=(a*h1+b*h2)/H
        before=G*(b*h2*h2/2+a*(H*H-h2*h2)/2)
        after=G*mean*H*H/2
        assert math.isclose(after-before,barrier(h1,h2,dr),rel_tol=1e-9,abs_tol=1e-8)
    assert wind_power(0,.0015,.4)==(0.,0.,0.)
    assert math.isclose(wind_power(2,.0015,.4)[2],8*wind_power(1,.0015,.4)[2])
    assert barrier(.1,2,0)==0
    for s in SIGMA:
        k=math.sqrt(G*RHO/s)
        assert wave(k,s)/k < wave(.99*k,s)/(.99*k)
        assert wave(k,s)/k < wave(1.01*k,s)/(1.01*k)
    for r in rows:
        assert 0<=r["accessible_bottom_m"]<=r["bottom_m"]
        assert r["power_W_m2"]<=r["stress_Pa"]*r["wind_m_s"]
        assert math.isclose(wind_power(r["threshold_wind_m_s"],r["drag"],r["mixing_coefficient"])[2]*r["hours"]*3600,r["required_J_m2"],rel_tol=1e-12)
    write_csv("mixing.csv",rows)
    write_csv("waves.csv",waves)
    ref=[r for r in rows if r["hours"]==24 and r["top_m"]==.1 and r["delta_rho"]==5 and r["drag"]==.0015 and r["mixing_coefficient"]==.4]
    reference=[dict(wind=r["wind_m_s"],bottom_m=r["bottom_m"],energy_J_m2=r["available_J_m2"],barrier_J_m2=r["required_J_m2"],threshold_wind=r["threshold_wind_m_s"],accessible_bottom_m=r["accessible_bottom_m"]) for r in ref]
    summary=dict(python=platform.python_version(),mixing_cases=len(rows),wave_cases=len(waves),energy_sufficient_cases=sum(r["energy_sufficient"] for r in rows),reference=reference,deep_water_minima=minima,controls=["zero wind","zero density contrast","cubic wind scaling","independent PE integration","energy threshold inversion","accessible depth bounds","deep-water phase-speed minimum"],assumptions=dict(g=G,liquid_density=RHO,gas_density=AIR,winds=WINDS,hours=HOURS,top=TOP,bottom=BOTTOM,contrast=CONTRAST,drag=DRAG,mixing_coefficient=EFF,surface_tension=SIGMA))
    (D/"results.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2))
if __name__=="__main__":
    main()
