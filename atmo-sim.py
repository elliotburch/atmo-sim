"""
Terminal Burn unified planetary evolution proof-of-concept.

Couples:
  bulk/formation inventory
  atmosphere + photochemistry
  climate / surface temperature
  condensation and volatile reservoirs
  representative vertical column / interfaces
  silicate crust state (solid / partial melt / molten)
  Jeans + XUV/diffusion-limited H + solar-wind/ion escape
  clathrate storage/release
  interface-limited crust/volatile exchange

This is a gameplay-oriented reduced-order model, not a precision climate/EOS code.
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional
import math, random, copy

# --- constants ---
G=6.67430e-11; ME=5.972e24; RE=6.371e6; BAR=1e5
SIGMA=5.670374419e-8; KB=1.380649e-23; NA=6.02214076e23
RGAS=8.314462618; EARTH_FLUX=1361.; SEC_GYR=3.15576e16
EARTH_VESC=11186.
SPECIES=("N2","CO2","CO","CH4","Ar","H2","He","H2O")
MM={"N2":28.0134e-3,"CO2":44.01e-3,"CO":28.010e-3,"CH4":16.043e-3,
    "Ar":39.948e-3,"H2":2.01588e-3,"He":4.0026e-3,"H2O":18.01528e-3}

def clamp(x,a=0.,b=1.): return max(a,min(b,x))
def smooth(a,b,x):
    if x<=a:return 0.
    if x>=b:return 1.
    q=(x-a)/(b-a); return q*q*(3-2*q)

class Material(str,Enum):
    ATMOSPHERE="atmosphere"; ORGANICS="organic deposits"; FROST="volatile frost"
    METHANE_LIQUID="liquid methane"; ICE_I="water ice Ih"; WATER="liquid water"
    HP_ICE="high-pressure water ice"; SILICATE="solid silicate crust"
    PARTIAL="partially molten silicate"; MAGMA="molten silicate"
    CARBON_CRUST="carbon-rich silicate crust"

@dataclass
class Crust:
    thickness_km:float=35.
    mafic:float=.65; felsic:float=.20; ultramafic:float=.15
    carbon_mass_fraction:float=.001; sulfide_fraction:float=.002
    oxidation:float=.2       # -1 reducing, +1 oxidized
    porosity:float=.05
    geothermal_gradient_k_km:float=15.

@dataclass
class Planet:
    name:str
    mass_e:float; radius_e:float
    flux:float; formation_flux:float
    seed:int=1
    water_fraction:Optional[float]=None
    carbon_fraction:Optional[float]=None
    nitrogen_fraction:Optional[float]=None
    sulfur_fraction:Optional[float]=None
    carbon_richness:float=0.0
    formation_redox:Optional[float]=None
    albedo_rock:float=.20
    tidally_locked:bool=False
    crust:Crust=field(default_factory=Crust)

@dataclass
class Layer:
    material:Material; thickness_m:float
    p_top_bar:float; p_bottom_bar:float
    t_top_k:float; t_bottom_k:float
    coupling_to_next:float=1.

@dataclass
class ClimateColumn:
    name:str
    area_fraction:float
    insolation_factor:float
    temperature_k:float=250.
    albedo:float=.25
    cloud_fraction:float=0.
    ice_fraction:float=0.
    ocean_fraction:float=0.
    desert_fraction:float=1.
    local_water_kg_m2:float=0.
    local_ice_kg_m2:float=0.
    local_ocean_kg_m2:float=0.
    absorbed_w_m2:float=0.
    surface_relative_humidity:float=0.
    precipitable_water_kg_m2:float=0.
    cloud_base_m:float=0.
    cloud_top_m:float=0.
    cloud_water_kg_m2:float=0.
    methane_cloud_kg_m2:float=0.
    methane_cloud_base_m:float=0.
    co2_cloud_kg_m2:float=0.
    co2_cloud_base_m:float=0.
    cold_trap_h2o_x:float=0.
    vertical_profile:List[dict]=field(default_factory=list)

@dataclass
class State:
    # Lower atmosphere is the climate/surface reservoir; upper atmosphere is the
    # photochemistry/escape reservoir. Both are kg/m2 and their sum sets surface pressure.
    atm:Dict[str,float]=field(default_factory=lambda:{s:0. for s in SPECIES})
    upper_atm:Dict[str,float]=field(default_factory=lambda:{s:0. for s in SPECIES})
    frost:Dict[str,float]=field(default_factory=lambda:{s:0. for s in SPECIES})
    clath:Dict[str,float]=field(default_factory=lambda:{s:0. for s in SPECIES})
    crust_vol:Dict[str,float]=field(default_factory=lambda:{s:0. for s in SPECIES})
    water_ice:float=0.; ocean:float=0.; hp_ice:float=0.; supercrit:float=0.
    methane_liquid:float=0.; organics:float=0.; refractory_carbon:float=0.
    T:float=250.; albedo:float=.3; runaway:bool=False
    dipole:float=0.; wind_coupling:float=1.
    tropopause_temperature_k:float=200.; tropopause_pressure_bar:float=.1
    upper_temperature_k:float=500.; cold_trap_h2o_mole_fraction:float=0.
    formation_candidate_pressure_bar:float=0.
    initial_target_pressure_bar:float=0.
    formation_retained_fraction:float=1.
    formation_impact_loss_kg_m2:float=0.
    bulk_water_fraction:float=0.
    bulk_carbon_fraction:float=0.
    bulk_nitrogen_fraction:float=0.
    bulk_sulfur_fraction:float=0.
    layers:List[Layer]=field(default_factory=list)
    silicate_state:Material=Material.SILICATE; melt_fraction:float=0.
    ocean_crust_coupling:float=0.; atmosphere_crust_coupling:float=0.
    escaped_kg_m2:Dict[str,float]=field(default_factory=lambda:{s:0. for s in SPECIES})
    climate_columns:List[ClimateColumn]=field(default_factory=list)

# --- thermodynamics / phase proxies ---
def psat_water(T):
    if T>=647.096:return 220.64
    if T<=120:return 0.
    # Wagner saturation curve
    tc=647.096; pc=220.64
    tau=1-T/tc
    a=(-7.85951783,1.84408259,-11.7866497,22.6807411,-15.9618719,1.80122502)
    ex=(tc/T)*(a[0]*tau+a[1]*tau**1.5+a[2]*tau**3+a[3]*tau**3.5+a[4]*tau**4+a[5]*tau**7.5)
    return pc*math.exp(ex)

COND={"CO2":(194.7,25200.,304.1),"CH4":(111.7,8200.,190.6),
      "CO":(81.6,6000.,132.9),"N2":(77.4,5600.,126.2),"Ar":(87.3,6400.,150.7)}
def psat_simple(T,Tb,L):
    return math.exp(clamp(-L/RGAS*(1/T-1/Tb),-40,40)) # bar, anchored ~1 bar at Tb

def luminosity(age): return 1/(1+.4*(1-age/4.5))

def derive_redox(p):
    # Formation composition drives initial crust oxidation. Carbon-rich systems consume available O.
    if p.formation_redox is not None:return clamp(p.formation_redox,-1,1)
    flux_term=.15*clamp(math.log10(max(p.formation_flux,1)/100)/2,-1,1)
    return clamp(p.crust.oxidation - 1.15*p.carbon_richness + flux_term,-1,1)

def draw_bulk_inventories(p,rng):
    """Draw independent bulk H2O/C/N/S mass fractions.

    H2O varies by orders of magnitude with formation temperature; C/N/S are
    independent trace-to-minor inventories rather than fixed percentages of a
    single generic volatile pool. Explicit Planet values override the draws.
    """
    f=max(p.formation_flux,1e-9)
    cold=f<190
    very_cold=f<15

    if cold:
        water=10**rng.uniform(-2.3,-0.55)      # ~0.5--28% bulk H2O
        carbon=10**rng.uniform(-4.0,-2.0)      # 0.01--1% C
        nitrogen=10**rng.uniform(-5.0,-3.0)    # 0.001--0.1% N
        sulfur=10**rng.uniform(-4.0,-2.0)
        if very_cold:
            carbon*=10**rng.uniform(.0,.35)
            nitrogen*=10**rng.uniform(.0,.35)
    else:
        water=10**rng.uniform(-5.2,-2.7)       # ~6 ppm--0.2% H2O
        carbon=10**rng.uniform(-5.0,-2.7)      # ~10 ppm--0.2% C
        nitrogen=10**rng.uniform(-6.0,-3.7)    # ~1 ppm--0.02% N
        sulfur=10**rng.uniform(-4.5,-2.3)

    # Carbon-rich systems genuinely contain more bulk carbon and tend to be drier.
    carbon*=10**(.9*max(0.,p.carbon_richness))
    water*=10**(-.55*max(0.,p.carbon_richness))

    if p.water_fraction is not None: water=p.water_fraction
    if p.carbon_fraction is not None: carbon=p.carbon_fraction
    if p.nitrogen_fraction is not None: nitrogen=p.nitrogen_fraction
    if p.sulfur_fraction is not None: sulfur=p.sulfur_fraction
    return dict(H2O=clamp(water,0.,.50), C=clamp(carbon,0.,.10),
                N=clamp(nitrogen,0.,.03), S=clamp(sulfur,0.,.10))

def initial_species_inventory(p,redox,bulk,masscol):
    """Reduced-order conversion of bulk elemental inventories to accessible species.

    Most bulk material remains in deep/inaccessible reservoirs. The returned
    dictionaries are species columns (kg/m2) for accessible and deep reservoirs.
    """
    oxid=(redox+1)/2
    cold=p.formation_flux<190

    # Fraction of each bulk reservoir initially participating in surface evolution.
    # Water is far more surface-accessible on ice-rich worlds; C/N remain mostly deep.
    aw=.72 if cold else .18
    ac=.18 if cold else .10
    an=.28 if cold else .18

    water=bulk["H2O"]*masscol
    carbon=bulk["C"]*masscol
    nitrogen=bulk["N"]*masscol

    # Convert elemental C/N masses to molecular-carrier masses.
    c_access=carbon*ac
    n_access=nitrogen*an
    h2o_access=water*aw

    # Carbon carrier split by redox. Molecular masses account for attached O/H.
    fco=.04+.22*(1-oxid)
    fch4=.005+.055*(1-oxid)**2
    fco2=max(.05,1-fco-fch4)
    # carbon mass fractions in molecules: CO2 12/44, CO 12/28, CH4 12/16.
    carbon_species={
        "CO2": c_access*fco2/(12/44),
        "CO": c_access*fco/(12/28),
        "CH4": c_access*fch4/(12/16),
    }
    # Nitrogen begins mainly as N2 in this version; NH3 can later use the same N budget.
    n2_access=n_access  # N2 is entirely nitrogen by elemental mass.

    accessible={s:0. for s in SPECIES}
    accessible["H2O"]=h2o_access
    accessible["N2"]=n2_access
    accessible.update(carbon_species)

    # Small redox-derived H2 inventory, decoupled from total water abundance.
    accessible["H2"]=(.00002+.00035*(1-oxid)**2)*masscol
    accessible["Ar"]=2e-6*masscol

    # Deep bookkeeping uses existing species-shaped crust reservoir.
    deep={s:0. for s in SPECIES}
    deep["H2O"]=water-h2o_access
    deep["N2"]=nitrogen-n_access
    deep_c=carbon-c_access
    deep["CO2"]=deep_c*fco2/(12/44)
    deep["CO"]=deep_c*fco/(12/28)
    deep["CH4"]=deep_c*fch4/(12/16)
    return accessible,deep

def silicate_phase(T,c):
    solidus=1250+120*c.ultramafic-80*c.felsic-60*clamp(c.carbon_mass_fraction/.05)
    liquidus=1850+100*c.ultramafic-60*c.felsic
    if T<solidus:
        return (Material.CARBON_CRUST if c.carbon_mass_fraction>.03 else Material.SILICATE),0.,solidus,liquidus
    if T>=liquidus:return Material.MAGMA,1.,solidus,liquidus
    return Material.PARTIAL,clamp((T-solidus)/(liquidus-solidus)),solidus,liquidus

# --- initialization ---
def initialize(p):
    rng=random.Random(p.seed)
    M=p.mass_e*ME; R=p.radius_e*RE; area=4*math.pi*R*R
    g=G*M/R**2; masscol=M/area
    redox=derive_redox(p)
    bulk=draw_bulk_inventories(p,rng)
    accessible,deep=initial_species_inventory(p,redox,bulk,masscol)
    st=State(T=max(45,(p.formation_flux*(1-.3)/(4*SIGMA))**.25))
    st.bulk_water_fraction=bulk["H2O"]
    st.bulk_carbon_fraction=bulk["C"]
    st.bulk_nitrogen_fraction=bulk["N"]
    st.bulk_sulfur_fraction=bulk["S"]
    for s in SPECIES:
        st.atm[s]=accessible[s]
        st.crust_vol[s]=deep[s]

    # Cold formation: accessible volatiles do not begin as a warm well-mixed gas.
    # Water freezes; CO2/CO/CH4 are predominantly condensed/trapped in ice/clathrates.
    if p.formation_flux < 190:
        coldT=max(35.,(p.formation_flux*(1-.35)/(4*SIGMA))**.25)
        # Water is essentially a condensed structural reservoir.
        st.water_ice += st.atm["H2O"]; st.atm["H2O"]=0.
        # Fractions remaining initially gaseous are deliberately temperature dependent.
        gasfrac={
            "CO2": .001 + .03*smooth(80,150,coldT),
            "CO":  .01  + .18*smooth(45,100,coldT),
            "CH4": .005 + .10*smooth(55,110,coldT),
        }
        for s in ("CO2","CO","CH4"):
            q=st.atm[s]
            gas=q*gasfrac[s]
            trapped=q-gas
            st.atm[s]=gas
            # Clathrate capacity is limited by available water ice.
            cap=max(0.,st.water_ice*.12-sum(st.clath.values()))
            cl=min(trapped,cap)
            st.clath[s]+=cl
            st.frost[s]+=trapped-cl
        # H2/He are poor condensables and remain gas initially; escape handles them.
    # Carbon-richness mostly appears as refractory crustal C, not atmospheric methane.
    st.refractory_carbon=bulk["C"]*masscol*.25*clamp(p.carbon_richness,0.,1.)
    # Rare nebular envelope.
    pcap=min(.15,.006+.08*(p.mass_e/5)**1.6)/(1+p.formation_flux/4000)
    if rng.random()<pcap:
        fenv=min(.002,1.2e-5*max(p.mass_e,.1)**2.2*math.exp(.8*rng.gauss(0,1)))
        env=fenv*masscol
        st.atm["H2"]+=.75*env; st.atm["He"]+=.25*env

    # Primordial erosion is deliberately deferred until after the first complete
    # phase/layer-cake solve in evolve(). At this point reservoirs are only chemically
    # initialized; some nominal atmospheric material may still condense.
    return st,bulk,redox


def draw_initial_atmosphere_target_bar(p,rng):
    """Post-accretion atmosphere retained after giant impacts/formation violence.

    Main branch is lognormal in pressure and clipped to 0.1--100 bar. Small/hot
    bodies have an additional chance to emerge nearly airless. This is only an
    initialization prior; subsequent evolution may exceed 100 bar or lose all gas.
    """
    M=p.mass_e*ME; R=p.radius_e*RE
    vesc=math.sqrt(2*G*M/R)
    teq=max(35.,(max(p.flux,1e-9)*(1-.25)/(4*SIGMA))**.25)

    # Weak physical bias: larger/cooler planets retain somewhat more atmosphere,
    # but stochastic formation history remains dominant.
    mu=math.log10(1.5) + .32*math.log10(max(p.mass_e,.01)) - .16*math.log10(max(teq,80.)/255.)
    sigma=.78
    target=10**rng.gauss(mu,sigma)
    target=clamp(target,.1,100.)

    # Airless/tenuous branch for low-gravity and strongly irradiated bodies.
    escape_score=clamp((5000.-vesc)/3500.,0.,1.) + .45*clamp((teq-350.)/500.,0.,1.)
    if rng.random() < .22*clamp(escape_score,0.,1.):
        target=10**rng.uniform(-6.,-1.)
    return target

def formation_impact_erosion(st,p,rng):
    """Uniformly erode the *gaseous* candidate atmosphere to a retained-pressure target.

    The atmosphere's mole/mass ratios are preserved to first order: this represents
    impact/accretion stripping rather than equilibrium condensation. Condensed and
    crustal reservoirs are untouched. Lost gas is recorded in escaped_kg_m2.
    """
    M=p.mass_e*ME; R=p.radius_e*RE
    g=G*M/R**2
    candidate=g*total_atmosphere(st)/BAR
    target=draw_initial_atmosphere_target_bar(p,rng)
    st.formation_candidate_pressure_bar=candidate
    st.initial_target_pressure_bar=target

    if candidate <= 0:
        st.formation_retained_fraction=1.
        return
    if candidate <= target:
        # Do not manufacture atmosphere if the candidate inventory is already smaller.
        st.formation_retained_fraction=1.
        st.initial_target_pressure_bar=candidate
        return

    keep=target/candidate
    lost_total=0.
    for box in (st.atm,st.upper_atm):
        for s in SPECIES:
            lost=box[s]*(1-keep)
            box[s]*=keep
            st.escaped_kg_m2[s]+=lost
            lost_total+=lost
    st.formation_retained_fraction=keep
    st.formation_impact_loss_kg_m2=lost_total

# --- climate / repartition ---
def atmospheric_column(st,s):
    return st.atm[s]+st.upper_atm[s]

def total_atmosphere(st):
    return sum(st.atm.values())+sum(st.upper_atm.values())

def atmosphere_structure(st,p,T,dt_gyr=0.0):
    """
    Diagnose tropopause and exchange lower <-> upper atmosphere.
    The upper box is a small, continuously replenished photochemical reservoir.
    Water transfer is cold-trap limited; noncondensables mix much more freely.
    """
    M=p.mass_e*ME; R=p.radius_e*RE; g=G*M/R**2
    lower=max(sum(st.atm.values()),1e-30)
    upper=max(sum(st.upper_atm.values()),0.)
    psurf=g*(lower+upper)/BAR

    # Convective lapse to a stratospheric skin-temperature floor.
    te=max(35.,(max(p.flux,1e-9)*(1-clamp(st.albedo,.03,.8))/(4*SIGMA))**.25)
    skin=max(45.,te/(2**.25))
    # Thick CO2 atmospheres have cold middle atmospheres; moist/hot worlds lift the cold trap.
    mol_lower=sum(st.atm[s]/MM[s] for s in SPECIES)
    xco2=(st.atm["CO2"]/MM["CO2"])/mol_lower if mol_lower else 0.
    moist_lift=70.*smooth(315,390,T)
    co2_cool=25.*clamp(xco2/.5)
    ttrop=clamp(skin+moist_lift-co2_cool,55.,320.)

    # Tropopause pressure declines for deep convective columns, but never becomes an absurdly
    # tiny fraction in this two-box approximation.
    ratio=clamp((ttrop/max(T,ttrop+1e-6))**3.5,.015,.35)
    ptrop=max(1e-7,psurf*ratio)
    st.tropopause_temperature_k=ttrop
    st.tropopause_pressure_bar=ptrop

    # Saturation-limited H2O mole fraction above cold trap.
    xsat=clamp(psat_water(ttrop)/max(ptrop,1e-12),0.,1.)
    st.cold_trap_h2o_mole_fraction=xsat

    # Upper reservoir target is a small fraction of total atmosphere. Exchange times are
    # far shorter than a geological epoch, so relax strongly each epoch but conserve mass.
    target_upper_frac=clamp(.004+.018*math.sqrt(max(T,50.)/300.),.003,.035)
    target_mass=(lower+upper)*target_upper_frac

    if dt_gyr<=0: mix=1.
    else: mix=1-math.exp(-dt_gyr/0.003)  # ~3 Myr bulk exchange proxy

    # Noncondensables: target upper composition follows lower atmosphere.
    lower_moles={s:st.atm[s]/MM[s] for s in SPECIES}
    lm=sum(lower_moles.values())
    xlower={s:(lower_moles[s]/lm if lm else 0.) for s in SPECIES}
    mean_mu=sum(xlower[s]*MM[s] for s in SPECIES) if lm else .028

    for s in SPECIES:
        if s=="H2O": continue
        # diffusive enrichment of H2/He in upper atmosphere.
        enrich=1.
        if s=="H2": enrich=4.
        elif s=="He": enrich=2.5
        desired=target_mass*xlower[s]*MM[s]/max(mean_mu,1e-9)*enrich
        desired=min(desired, atmospheric_column(st,s))
        delta=(desired-st.upper_atm[s])*mix
        if delta>0:
            d=min(delta,st.atm[s]); st.atm[s]-=d; st.upper_atm[s]+=d
        else:
            d=min(-delta,st.upper_atm[s]); st.upper_atm[s]-=d; st.atm[s]+=d

    # Water: cold trap sets upper mole fraction, not lower surface humidity.
    nonwater_upper_moles=sum(st.upper_atm[s]/MM[s] for s in SPECIES if s!="H2O")
    desired_h2o_moles=(xsat/(max(1-xsat,1e-12)))*nonwater_upper_moles if xsat<.999 else lm*.02
    desired=min(desired_h2o_moles*MM["H2O"], atmospheric_column(st,"H2O"))
    delta=(desired-st.upper_atm["H2O"])*mix
    if delta>0:
        d=min(delta,st.atm["H2O"]); st.atm["H2O"]-=d; st.upper_atm["H2O"]+=d
    else:
        d=min(-delta,st.upper_atm["H2O"]); st.upper_atm["H2O"]-=d; st.atm["H2O"]+=d

def repartition(st,p,T):
    """Authoritative phase partition. Every accessible molecule is assigned exactly once."""
    M=p.mass_e*ME; R=p.radius_e*RE; g=G*M/R**2

    # Water total: collect all surface-accessible phases exactly once.
    W=max(0.,st.atm["H2O"]+st.upper_atm["H2O"]+st.water_ice+st.ocean+st.hp_ice+st.supercrit)
    st.atm["H2O"]=0.; st.upper_atm["H2O"]=0.; st.water_ice=0.; st.ocean=0.; st.hp_ice=0.; st.supercrit=0.

    if T>=647.096:
        active=min(W,300*BAR/g)
        st.atm["H2O"]=active; st.supercrit=W-active
    else:
        vapor=min(W,psat_water(T)*BAR/g)
        st.atm["H2O"]=vapor
        rem=W-vapor

        if rem>0:
            # First determine surface condensed state.
            if T<273.15:
                # Conductive Ice-I lid until ~melting; remaining water may be subglacial liquid.
                grad=max(p.crust.geothermal_gradient_k_km,1e-3)
                ice_depth=max(0.,(273.15-T)/grad)*1000.
                ice_cap=ice_depth*920.
                ice=min(rem,ice_cap)
                st.water_ice=ice; rem-=ice
            # Liquid layer extends only until the HP-ice stability pressure.
            if rem>0:
                overburden_bar=(st.water_ice*g/BAR)
                hp_threshold=max(4500.,6200.+5.*(max(T,260.)-273.))
                liquid_cap=max(0.,(hp_threshold-overburden_bar)*BAR/g)
                liq=min(rem,liquid_cap)
                st.ocean=liq; rem-=liq
            if rem>0:
                st.hp_ice=rem

    # Other condensables. Collect methane liquid exactly once.
    for s,(tb,L,crit) in COND.items():
        total=max(0.,st.atm[s]+st.upper_atm[s]+st.frost[s]+(st.methane_liquid if s=="CH4" else 0.))
        st.atm[s]=0.; st.upper_atm[s]=0.; st.frost[s]=0.
        if s=="CH4": st.methane_liquid=0.
        gas=total if T>=crit else min(total,psat_simple(T,tb,L)*BAR/g)
        if s=="CO" and T<120 and p.formation_flux<190:
            # CO is efficiently trapped/processed in cold H2O/N2-rich surface reservoirs.
            n2_col=max(st.atm["N2"],1e-30)
            co_ratio_cap=.02 + .08*smooth(60,110,T)
            gas=min(gas,n2_col*co_ratio_cap)
        st.atm[s]=gas; condensed=total-gas
        if s=="CH4" and 90.7<T<190.6: st.methane_liquid=condensed
        else: st.frost[s]=condensed

def build_column(st,p,T):
    """Build a representative vertical stack from current reservoirs WITHOUT moving mass."""
    M=p.mass_e*ME; R=p.radius_e*RE; g=G*M/R**2
    P=g*total_atmosphere(st)/BAR
    layers=[Layer(Material.ATMOSPHERE,max(1e3,4*8000*T/288*9.81/g),0,P,.75*T,T)]

    def add(mat,mass,rho,t0,t1,coup):
        nonlocal P
        if mass<=0:return
        p0=P; P+=mass*g/BAR
        layers.append(Layer(mat,mass/rho,p0,P,t0,t1,coup))

    if st.organics>0:add(Material.ORGANICS,st.organics,850,T,T,.2)
    frost=sum(st.frost.values())
    if frost>0:add(Material.FROST,frost,800,T,T,.25)
    if st.methane_liquid>0:add(Material.METHANE_LIQUID,st.methane_liquid,450,T,T,.7)

    grad=max(p.crust.geothermal_gradient_k_km,0.)
    currentT=T
    if st.water_ice>0:
        t1=min(273.15,currentT+grad*(st.water_ice/920.)/1000.)
        add(Material.ICE_I,st.water_ice,920,currentT,t1,.08); currentT=t1
    if st.ocean>0:
        t1=currentT+grad*(st.ocean/1000.)/1000.
        add(Material.WATER,st.ocean,1000,currentT,t1,.85); currentT=t1
    if st.hp_ice>0:
        t1=currentT+grad*(st.hp_ice/1300.)/1000.
        add(Material.HP_ICE,st.hp_ice,1300,currentT,t1,.015); currentT=t1
    if st.supercrit>0:
        # Dense supercritical water is represented structurally as a fluid layer.
        t1=currentT+grad*(st.supercrit/700.)/1000.
        add(Material.WATER,st.supercrit,700,currentT,t1,.5); currentT=t1

    over=sum(x.thickness_m for x in layers[1:])
    # Volatile shells convect; a linear geothermal gradient across 100+ km is unphysical.
    # Use a capped basal warming proxy. Exposed rock still equals surface T.
    volatile_over=sum(x.thickness_m for x in layers[1:] if x.material not in
                      (Material.SILICATE,Material.CARBON_CRUST,Material.PARTIAL,Material.MAGMA))
    basal_delta=min(450., grad*volatile_over/1000.)
    Trock=T+basal_delta
    mat,melt,_,_=silicate_phase(Trock,p.crust)
    crust_rho=2700+600*p.crust.mafic+400*p.crust.ultramafic
    add(mat,p.crust.thickness_km*1000*crust_rho,crust_rho,Trock,Trock+grad*p.crust.thickness_km,.8)
    st.layers=layers; st.silicate_state=mat; st.melt_fraction=melt

    atmcr=1.; oceancr=0.; seen_hp=False
    for L in layers[1:]:
        if L.material==Material.HP_ICE: seen_hp=True
        atmcr*=L.coupling_to_next
        if L.material==Material.WATER and not seen_hp: oceancr=.85
        if L.material==Material.HP_ICE: oceancr=0.
        if L.material in (Material.SILICATE,Material.CARBON_CRUST,Material.PARTIAL,Material.MAGMA):break
    if mat in (Material.PARTIAL,Material.MAGMA): atmcr=max(atmcr,.2+.7*melt)
    st.atmosphere_crust_coupling=clamp(atmcr)
    st.ocean_crust_coupling=clamp(oceancr)


LATENT={"H2O":2.50e6,"CH4":5.1e5,"CO2":5.7e5}
CP_MOLAR={"N2":29.1,"CO2":37.1,"CO":29.1,"CH4":35.7,"Ar":20.8,"H2":28.8,"He":20.8,"H2O":33.6}

def species_psat_bar(s,T):
    if s=="H2O": return psat_water(T)
    tb,L,crit=COND[s]
    if T>=crit:return 1e12
    return psat_simple(max(T,20.),tb,L)

def vertical_condensable_profile(st,p,c,nlev=32):
    """Diagnostic convective column in pressure coordinates.

    Unsaturated air follows a dry adiabat. Once a condensable reaches saturation,
    latent heating reduces the lapse rate and excess vapor is diagnosed as cloud
    condensate. This does not create/destroy global volatile mass; it diagnoses the
    atmospheric fraction that should remain vapor.
    """
    M=p.mass_e*ME; R=p.radius_e*RE; g=G*M/R**2
    ps=max(1e-5,g*total_atmosphere(st)/BAR)
    pt=max(1e-6,min(st.tropopause_pressure_bar,ps*.35))
    if pt>=ps:pt=max(1e-6,ps*.08)

    # Background composition and heat capacity.
    moles={s:atmospheric_column(st,s)/MM[s] for s in SPECIES}
    mt=max(sum(moles.values()),1e-30)
    x={s:moles[s]/mt for s in SPECIES}
    mu=max(.002,sum(x[s]*MM[s] for s in SPECIES))
    cp_molar=sum(x[s]*CP_MOLAR.get(s,30.) for s in SPECIES)
    cp=max(500.,cp_molar/mu)
    rd=RGAS/mu
    kappa=clamp(rd/cp,.08,.40)

    # Local surface RH comes from local surface state, not global water existence.
    rh=clamp(.10*c.desert_fraction+.38*c.ice_fraction+.78*c.ocean_fraction,.03,.95)
    c.surface_relative_humidity=rh

    # Pressure interfaces, surface -> tropopause.
    lnps=math.log(ps); lnpt=math.log(pt)
    pint=[math.exp(lnps+(lnpt-lnps)*i/nlev) for i in range(nlev+1)]
    T=c.temperature_k
    z=0.
    prof=[]
    vapor_x={}
    cloud_mass={"H2O":0.,"CH4":0.,"CO2":0.}
    cloud_base={"H2O":None,"CH4":None,"CO2":None}
    cloud_top={"H2O":None,"CH4":None,"CO2":None}

    # Surface vapor abundance. Water is humidity-limited; CH4/CO2 begin from global gas
    # mixing ratios and only condense if the local adiabat intersects their saturation curve.
    vapor_x["H2O"]=min(.95,rh*species_psat_bar("H2O",T)/ps)
    vapor_x["CH4"]=x.get("CH4",0.)
    vapor_x["CO2"]=x.get("CO2",0.)

    water_vapor_col=0.
    for i in range(nlev):
        p1,p2=pint[i],pint[i+1]
        pm=math.sqrt(p1*p2)
        dp=(p1-p2)*BAR
        # Saturation at current level; diagnose which species condense.
        active=[]
        satx={}
        for s in ("H2O","CH4","CO2"):
            satx[s]=clamp(species_psat_bar(s,T)/max(pm,1e-12),0.,.999)
            if vapor_x[s] > satx[s]:
                active.append(s)

        # Latent heating proxy: exact pseudoadiabat for one dominant condensable,
        # smoothly blended for multiple simultaneous condensables.
        moist_factor=1.
        for s in active:
            eps=MM[s]/mu
            qs=max(0.,eps*satx[s]/max(1-satx[s],1e-9))
            L=LATENT[s]
            gm=(1+L*qs/(rd*T))/max(1+L*L*qs/(cp*max(RGAS/MM[s],1.)*T*T),1e-9)
            moist_factor=min(moist_factor,clamp(gm,.22,1.))
        exponent=kappa*moist_factor
        T2=max(35.,T*(p2/p1)**exponent)

        # Hydrostatic layer thickness using mean T.
        tmid=.5*(T+T2)
        dz=max(0.,rd*tmid/g*math.log(p1/p2))
        layer_air=dp/g

        row={"z_m":z,"pressure_bar":p1,"temperature_k":T}
        for s in ("H2O","CH4","CO2"):
            xsat=clamp(species_psat_bar(s,T)/max(pm,1e-12),0.,.999)
            xv=min(vapor_x[s],xsat)
            # approximate mass mixing ratio from mole fraction
            q=clamp(xv*MM[s]/mu,0.,.95)
            vm=q*layer_air
            if s=="H2O": water_vapor_col += vm
            excess=max(0.,vapor_x[s]-xsat)
            if excess>0:
                cond=clamp(excess*MM[s]/mu,0.,.95)*layer_air
                cloud_mass[s]+=cond
                if cloud_base[s] is None:cloud_base[s]=z
                cloud_top[s]=z+dz
                vapor_x[s]=xsat
            row[s+"_x"]=xv
        prof.append(row)
        z+=dz; T=T2

    # top/cold-trap state
    c.precipitable_water_kg_m2=water_vapor_col
    c.cloud_water_kg_m2=cloud_mass["H2O"]
    c.methane_cloud_kg_m2=cloud_mass["CH4"]
    c.co2_cloud_kg_m2=cloud_mass["CO2"]
    c.cloud_base_m=cloud_base["H2O"] or 0.
    c.cloud_top_m=cloud_top["H2O"] or c.cloud_base_m
    c.methane_cloud_base_m=cloud_base["CH4"] or 0.
    c.co2_cloud_base_m=cloud_base["CO2"] or 0.
    c.cold_trap_h2o_x=min(vapor_x["H2O"],species_psat_bar("H2O",T)/max(pt,1e-12))
    c.vertical_profile=prof
    return c

def diagnose_vertical_columns(st,p):
    for c in st.climate_columns:
        vertical_condensable_profile(st,p,c)
    if st.climate_columns:
        # Area-weighted atmospheric water target. Keep this diagnostic conservative:
        # repartition owns total water mass; this only caps vapor to what the convective
        # profiles can support and returns excess to condensed reservoirs.
        target=sum(c.area_fraction*c.precipitable_water_kg_m2 for c in st.climate_columns)
        current=st.atm["H2O"]
        if current>target:
            excess=current-target
            st.atm["H2O"]=target
            # Deposit globally; the next three-column solve redistributes toward cold traps.
            if st.T<273.15: st.water_ice+=excess
            else: st.ocean+=excess
        st.cold_trap_h2o_mole_fraction=sum(c.area_fraction*c.cold_trap_h2o_x for c in st.climate_columns)

def column_geometry(p):
    if p.tidally_locked:
        # Area-weighted representative day/terminator/night regions. Factors have weighted mean 1.
        return [("substellar",.35,2.20),("terminator",.35,.80),("antistellar",.30,0.0)]
    # Latitude bands: equator |lat|<30, mid 30-60, poles >60.
    return [("equator",.50,1.24),("midlatitude",.366,0.88),("pole",.134,.30)]

def solve_three_columns(st,p):
    """
    Three horizontally coupled climate columns sharing one atmosphere.
    Condensed H2O is geographically redistributed between columns; atmospheric H2O remains global.
    """
    M=p.mass_e*ME; R=p.radius_e*RE; g=G*M/R**2
    geom=column_geometry(p)
    total_atm=max(total_atmosphere(st),1e-30)
    pbar=g*total_atm/BAR
    mol=sum(atmospheric_column(st,s)/MM[s] for s in SPECIES)
    xco2=(atmospheric_column(st,"CO2")/MM["CO2"])/mol if mol else 0.
    xch4=(atmospheric_column(st,"CH4")/MM["CH4"])/mol if mol else 0.
    ph2o=g*atmospheric_column(st,"H2O")/BAR

    # Shared dry-greenhouse amplification.
    pco2=g*atmospheric_column(st,"CO2")/BAR
    pch4=g*atmospheric_column(st,"CH4")/BAR
    ph2=g*atmospheric_column(st,"H2")/BAR
    pn2=g*atmospheric_column(st,"N2")/BAR
    tau=.42*math.log1p(pco2/.02)+.16*math.log1p(pch4/.003)+.10*math.sqrt(max(ph2,0))
    tau*=1+.08*math.log1p(max(pn2+pco2,0))
    greenhouse=(1+.75*max(tau,0))**.25

    # Horizontal transport efficiency grows strongly with pressure. Locked planets need
    # somewhat more atmosphere to erase day/night contrast.
    # Pressure-efficient horizontal heat transport.  The old logarithmic law left
    # Earth-pressure worlds with unrealistically isolated latitude bands.  This proxy
    # rises rapidly through the 0.1--3 bar regime and asymptotes toward near-isothermal
    # redistribution for dense atmospheres.
    transport=clamp(.10 + .72*(1-math.exp(-max(pbar,0)/0.85)), .10, .94)
    if p.tidally_locked:
        # Day/night circulation is less geometrically efficient than zonal/meridional
        # transport, but multi-bar atmospheres still approach strong redistribution.
        transport=clamp(.07 + .78*(1-math.exp(-max(pbar,0)/2.0)), .07, .90)

    condensed=max(0.,st.water_ice+st.ocean+st.hp_ice)
    if not st.climate_columns or [c.name for c in st.climate_columns]!=[x[0] for x in geom]:
        st.climate_columns=[ClimateColumn(n,a,f,temperature_k=st.T) for n,a,f in geom]

    # Initial local condensed-water allocation from previous step, otherwise uniform.
    if sum(c.local_water_kg_m2*c.area_fraction for c in st.climate_columns)<=0 and condensed>0:
        for c in st.climate_columns: c.local_water_kg_m2=condensed

    # Iterate radiation, clouds, ice-albedo, and horizontal sensible heat transport.
    temps=[clamp(c.temperature_k,45,2000) for c in st.climate_columns]
    for _ in range(45):
        raw=[]; albs=[]; clouds=[]
        tmean=sum(c.area_fraction*t for c,t in zip(st.climate_columns,temps))
        for j,(c,T) in enumerate(zip(st.climate_columns,temps)):
            # Surface fractions derive from local water and T.
            has_water=c.local_water_kg_m2>5
            cover=1-math.exp(-c.local_water_kg_m2/2500.) if has_water else 0.
            ice=cover*(1-smooth(245,278,T)) if has_water else 0.
            liquid=cover*(1-ice/max(cover,1e-12))*smooth(268,285,T) if has_water and T<647 else 0.
            # Approx relative humidity: wet surfaces humid, deserts dry; warm oceans convect.
            rh=(.72*liquid+.35*ice+.10*(1-liquid-ice))
            conv=smooth(285,330,T)*rh
            cloud=clamp(.06+.58*conv+.12*math.sqrt(max(pbar,0)/max(1+pbar,1e-9)),.03,.82)
            # Thick substellar cloud decks on wet locked worlds.
            if p.tidally_locked and c.name=="substellar" and liquid>.2:
                cloud=clamp(cloud+.18*smooth(290,330,T),.03,.9)
            surfaceA=(1-liquid-ice)*p.albedo_rock + liquid*.07 + ice*.62
            # Clouds mask the underlying surface; haze is global.
            haze=.10*clamp(xch4/.03)
            A=clamp(surfaceA*(1-cloud)+.52*cloud+haze,.03,.82)
            local_abs=p.flux*c.insolation_factor*(1-A)/4.
            te=max(30.,(max(local_abs,1e-9)/SIGMA)**.25)
            humid=clamp(math.log1p(ph2o/.01)/math.log(1+30/.01))
            waterwarm=(12+35*smooth(285,340,T))*humid*rh
            trad=te*greenhouse+waterwarm
            # Relax toward area mean according to atmospheric transport.
            tmix=(1-transport)*trad+transport*tmean
            raw.append(tmix); albs.append(A); clouds.append(cloud)

        newtemps=[.55*t+.45*n for t,n in zip(temps,raw)]
        if max(abs(a-b) for a,b in zip(newtemps,temps))<.03:
            temps=newtemps; break
        temps=newtemps

    # Redistribute condensed surface water toward cold traps, without changing global mass.
    # Warm habitable columns retain some water; very hot desert columns lose most to colder regions.
    if condensed>0:
        weights=[]
        for c,T in zip(st.climate_columns,temps):
            cold=math.exp(clamp((300-T)/32,-8,8))
            liquid_bonus=1.0 if 273<T<330 else .25
            weights.append(c.area_fraction*(cold+liquid_bonus))
        norm=sum(weights)
        for c,w in zip(st.climate_columns,weights):
            # local column mass such that area-weighted sum equals global condensed inventory
            c.local_water_kg_m2=condensed*(w/max(c.area_fraction,1e-9))/norm

    for c,T,A,cloud in zip(st.climate_columns,temps,albs,clouds):
        c.temperature_k=T; c.albedo=A; c.cloud_fraction=cloud
        cover=1-math.exp(-c.local_water_kg_m2/2500.) if c.local_water_kg_m2>5 else 0.
        c.ice_fraction=cover*(1-smooth(245,278,T)) if cover>0 else 0.
        c.ocean_fraction=(cover-c.ice_fraction)*smooth(268,285,T) if cover>0 and T<647 else 0.
        c.desert_fraction=max(0.,1-c.ice_fraction-c.ocean_fraction)
        c.local_ice_kg_m2=c.local_water_kg_m2*c.ice_fraction
        c.local_ocean_kg_m2=c.local_water_kg_m2*c.ocean_fraction
        c.absorbed_w_m2=p.flux*c.insolation_factor*(1-c.albedo)/4.

    # Area-weighted climate feeds the global chemistry/column model.
    st.T=sum(c.area_fraction*c.temperature_k for c in st.climate_columns)
    st.albedo=sum(c.area_fraction*c.albedo for c in st.climate_columns)
    diagnose_vertical_columns(st,p)
    return st.T,st.albedo

def climate(st,p,maxiter=50):
    """
    Reduced-order radiative-convective climate.

    Key behavior:
      * dry/noncondensing greenhouse gases warm via an emitting-level proxy;
      * H2O vapor follows saturation, but cannot create arbitrary grey feedback;
      * a wet-atmosphere OLR ceiling determines runaway rather than temperature alone;
      * below the ceiling, ocean-bearing states converge to a stable moist equilibrium.
    """
    M=p.mass_e*ME; R=p.radius_e*RE; g=G*M/R**2
    direct_floor=max(30.,(p.flux*(1-clamp(p.albedo_rock,.03,.78))/(4*SIGMA))**.25)
    T=clamp(st.T,45,2500)
    runaway=False

    for _ in range(maxiter):
        repartition(st,p,T)
        total_atm=max(total_atmosphere(st),1e-30)
        pbar=g*total_atm/BAR

        icecov=1-math.exp(-st.water_ice/2000) if st.water_ice else 0
        oceancov=(1-math.exp(-st.ocean/5000))*max(0,1-icecov)
        frost=sum(st.frost.values())
        frostcov=(1-math.exp(-frost/300))*max(0,1-icecov-oceancov)
        orgcov=(1-math.exp(-st.organics/500))*max(0,1-icecov-oceancov-frostcov)
        rock=max(0,1-icecov-oceancov-frostcov-orgcov)

        # Clouds rise on warm wet worlds and provide a stabilizing albedo feedback.
        wet_inventory=atmospheric_column(st,"H2O")+st.ocean+st.water_ice+st.hp_ice+st.supercrit
        water_mix=atmospheric_column(st,"H2O")/total_atm
        cloud=.16*smooth(270,330,T)*clamp(water_mix/.08) if wet_inventory>100 else 0.
        haze=.10*clamp(st.atm["CH4"]/total_atm/.03)
        A=clamp(rock*p.albedo_rock+icecov*.58+oceancov*.07+frostcov*.55+orgcov*.12+cloud+haze,.03,.78)
        if icecov+oceancov+frostcov+orgcov < 1e-6:
            A=clamp(p.albedo_rock+cloud+haze,.03,.78)

        absorbed=p.flux*(1-A)/4.
        Te=max(30.,(absorbed/SIGMA)**.25)

        # Non-condensing greenhouse: pressure broadening + composition-dependent opacity.
        pco2=g*atmospheric_column(st,"CO2")/BAR
        pch4=g*atmospheric_column(st,"CH4")/BAR
        ph2=g*atmospheric_column(st,"H2")/BAR
        pn2=g*atmospheric_column(st,"N2")/BAR
        tau_dry=0.
        tau_dry += .42*math.log1p(pco2/.02)
        tau_dry += .16*math.log1p(pch4/.003)
        tau_dry += .10*math.sqrt(max(ph2,0.))
        broad=1+.08*math.log1p(max(pn2+pco2,0.))
        tau_dry*=broad
        Tdry=Te*(1+.75*max(tau_dry,0))**.25

        # Water feedback below runaway: bounded radiative-convective warming.
        # Humidity effect approaches a finite increment rather than diverging with column.
        ph2o=g*atmospheric_column(st,"H2O")/BAR
        humid=clamp(math.log1p(ph2o/.01)/math.log(1+30/.01))
        water_warm=(18.+42.*smooth(285,340,T))*humid
        stable_T=Tdry+water_warm

        # Approximate Simpson-Nakajima OLR ceiling. Dry/background gases alter it modestly.
        # This is deliberately a broad physical proxy rather than a fitted Earth-only constant.
        olr_limit=282.
        olr_limit += 8.*math.log1p(max(pn2,0.))
        olr_limit -= 7.*math.log1p(max(pco2,0.))
        olr_limit=clamp(olr_limit,245.,315.)

        wet_surface=(st.ocean+st.water_ice+st.hp_ice+st.supercrit)>100.
        # Runaway is a radiative condition. Once a wet lower atmosphere has enough vapor
        # to approach the OLR ceiling, absorbed flux above that ceiling has no cool solution.
        moist_candidate=wet_surface and stable_T>315 and ph2o>.01
        runaway=bool(moist_candidate and absorbed>olr_limit)

        if runaway:
            excess=max(0.,absorbed-olr_limit)
            # Deterministic steam branch; avoid branch-flipping between climate iterations.
            Tnew=min(2500.,max(647.1,Tdry+220.+1.25*excess))
        else:
            # Stable branch. Only genuinely wet worlds are restricted to the subcritical
            # water branch. Dry/bare worlds may become arbitrarily hot and melt silicate.
            Tnew=max(Te,stable_T,direct_floor if rock>.98 else 0.)
            if wet_surface:
                Tnew=min(Tnew,646.)
            else:
                Tnew=min(Tnew,2500.)

        if abs(Tnew-T)<.05:
            T=Tnew; break
        T=.60*T+.40*Tnew

    st.T=clamp(T,45,2500); st.albedo=A; st.runaway=runaway
    repartition(st,p,st.T)
    # Resolve geographic climate after the global radiative branch. The area-weighted result
    # becomes the representative temperature for chemistry and vertical structure.
    solve_three_columns(st,p)
    atmosphere_structure(st,p,st.T,0.0)
    build_column(st,p,st.T)

# --- photochemistry and exchange ---
def chemistry_epoch(st,p,age,dt):
    uv=(p.flux/EARTH_FLUX)*(1+8*math.exp(-age/.6))
    red=derive_redox(p); oxid=(red+1)/2

    # Replenish upper atmosphere through mixing before UV chemistry.
    atmosphere_structure(st,p,st.T,dt)
    heavy=max(total_atmosphere(st),1e-30)

    # Methane photolysis occurs in the upper reservoir.
    shielding=1/(1+.35*math.log1p(max(heavy,0)))
    k_ch4=24*uv*shielding
    destroyed=st.upper_atm["CH4"]*(1-math.exp(-k_ch4*dt))
    st.upper_atm["CH4"]-=destroyed
    to_co2=.18+.52*oxid
    to_org=.12+.38*(1-oxid)
    to_co=max(0.,1-to_co2-to_org)
    # Heavy photochemical products settle/mix back into lower atmosphere.
    st.atm["CO2"]+=destroyed*to_co2
    st.atm["CO"]+=destroyed*to_co
    st.organics+=destroyed*to_org
    st.upper_atm["H2"]+=destroyed*.10

    # CO has both photochemical and surface/ocean sinks. Wet oxidizing surfaces are fast;
    # even dry atmospheres have a slow photochemical conversion pathway.
    wet_surface=smooth(235,320,st.T)*max(st.ocean_crust_coupling,st.atmosphere_crust_coupling*.35)
    photo_sink=.08*uv/(1+.25*math.log1p(max(heavy,0)))
    surface_sink=1.2*oxid*wet_surface
    k_co=photo_sink+surface_sink
    conv=st.atm["CO"]*(1-math.exp(-k_co*dt))
    st.atm["CO"]-=conv
    st.atm["CO2"]+=conv

    # On cryogenic ice worlds, surface irradiation and trapping remove CO over Gyr.
    # Route it to condensed/processed carbon rather than allowing primordial CO to dominate forever.
    if st.T<120 and p.formation_flux<190:
        cold_k=.35*(1+clamp((120-st.T)/70))
        cold_loss=st.atm["CO"]*(1-math.exp(-cold_k*dt))
        st.atm["CO"]-=cold_loss
        st.organics+=.35*cold_loss
        st.frost["CO2"]+=.65*cold_loss

    # Interior/surface-shell volatile release. Rocky worlds use volcanic speciation;
    # cold ice-rich bodies use slow cryovolcanic/leakage chemistry instead.
    coupling=max(.002,st.atmosphere_crust_coupling)
    ice_rich=(p.formation_flux<190 and st.bulk_water_fraction>.005)
    if ice_rich and st.T<180:
        # Cold icy shells strongly throttle deep exchange. Released C is mainly CH4/CO2;
        # CO and especially H2 are not treated as persistent bulk cryovolcanic products.
        activity=(.006+.025*clamp(p.mass_e/.1))*math.exp(-age/6.)
        frac=1-math.exp(-activity*dt)
        released={s:0. for s in SPECIES}
        for s in SPECIES:
            d=st.crust_vol[s]*frac
            st.crust_vol[s]-=d; released[s]=d
        carbon_release=released["CO2"]+released["CO"]+released["CH4"]
        fch4=.15+.45*(1-oxid)
        fco2=1-fch4
        st.atm["CO2"]+=carbon_release*fco2
        st.atm["CH4"]+=carbon_release*fch4
        # N2/Ar can leak from an icy interior; water largely refreezes locally.
        st.atm["N2"]+=released["N2"]; st.atm["Ar"]+=released["Ar"]
        st.frost["H2O"]+=0.0  # H2O handled below as condensed water inventory
        st.water_ice+=released["H2O"]
        # H2/He escape through porous/cryovolcanic pathways rather than accumulate.
        st.escaped_kg_m2["H2"]+=released["H2"]
        st.escaped_kg_m2["He"]+=released["He"]
    else:
        activity=(.08+.92*st.melt_fraction)*math.exp(-age/3.5)
        frac=1-math.exp(-.25*coupling*activity*dt)
        carbon_release=0.
        for s in SPECIES:
            d=st.crust_vol[s]*frac
            st.crust_vol[s]-=d
            if s in ("CO2","CO","CH4"):
                carbon_release+=d
            else:
                st.atm[s]+=d
        if carbon_release>0:
            fco=.025+.20*(1-oxid)
            fch4=.002+.035*(1-oxid)**2
            fco2=max(0.,1-fco-fch4)
            st.atm["CO2"]+=carbon_release*fco2
            st.atm["CO"]+=carbon_release*fco
            st.atm["CH4"]+=carbon_release*fch4

# --- escape ---
def magnetic_dipole(m,r,t,phase):
    life=1+3.5*math.sqrt(m); size=(m/r**2)**.7
    secular=.75+.25*math.sin(2*math.pi*(t/.7+phase))
    reversal=.15 if math.sin(2*math.pi*(t/.23+phase))>.92 else 1
    return max(0,size*math.exp(-max(0,t-life)/(0.8+life))*secular*reversal)

def escape_epoch(st,p,age,dt,phase):
    M=p.mass_e*ME; R=p.radius_e*RE; g=G*M/R**2
    vesc=math.sqrt(2*G*M/R)
    young=1+50*math.exp(-age/.3)
    xuv=4.6e-3*(p.flux/EARTH_FLUX)*young
    sw=(p.flux/EARTH_FLUX)*(1+35*math.exp(-age/.45))
    dip=magnetic_dipole(p.mass_e,p.radius_e,age,phase); st.dipole=dip

    atmosphere_structure(st,p,st.T,dt)

    rmp=max(1,10*dip**(1/3)*max(sw,.02)**(-1/6)/max(p.radius_e,.05))
    mol=sum(st.upper_atm[s]/MM[s] for s in SPECIES)
    mu=sum(st.upper_atm.values())/mol if mol else .028
    xco2=(st.upper_atm["CO2"]/MM["CO2"])/mol if mol else 0
    tex=(1.3*max(40,(p.flux*(1-st.albedo)/(4*SIGMA))**.25)+500*math.sqrt(max(p.flux/EARTH_FLUX,0)))*(1+2*math.exp(-age/.3))
    tex*=1-.65*xco2/(xco2+.05)
    st.upper_temperature_k=max(45.,tex)
    H=RGAS*max(tex,50)/(max(mu,.002)*g)
    exobase=1+min(.8,8*H/R)
    clearance=max(0,rmp-exobase)
    coupling=min(1,max(.10+.08*min(1,dip),1/(1+(clearance/2.5)**2)))
    st.wind_coupling=coupling

    # Thermal escape acts only on material actually present aloft.
    for s in ("H2","He","CH4","N2","CO","Ar","CO2"):
        lam=G*M*(MM[s]/NA)/(KB*max(tex,50)*R)
        rate=.8*math.exp(-max(lam-20,0)/3)
        frac=1-math.exp(-rate*dt)
        lost=st.upper_atm[s]*frac
        st.upper_atm[s]-=lost; st.escaped_kg_m2[s]+=lost

    if age < .8 and p.mass_e < .08:
        smallness=clamp((.08-p.mass_e)/.08)
        early=math.exp(-age/.18)
        f_blow=1-math.exp(-80*smallness*early*dt)
        for s in ("H2","He"):
            available=st.atm[s]+st.upper_atm[s]
            requested=available*f_blow
            du=min(st.upper_atm[s],requested)
            st.upper_atm[s]-=du
            rem=requested-du
            dl=min(st.atm[s],rem)
            st.atm[s]-=dl
            st.escaped_kg_m2[s]+=du+dl

    for s,boost in (("H2",2.5),("He",1.8)):
        lam=G*M*(MM[s]/NA)/(KB*max(tex,50)*R)
        hydro_rate=boost*young*(p.flux/EARTH_FLUX+0.015)*math.exp(-clamp(lam,0,80)/7)
        f=1-math.exp(-min(hydro_rate*dt,30))
        lost=st.upper_atm[s]*f
        st.upper_atm[s]-=lost; st.escaped_kg_m2[s]+=lost
        if p.mass_e < .01 and s in ("H2","He"):
            # Tiny icy bodies cannot hide a long-lived lower H2 reservoir from an escaping exosphere.
            supply=1-math.exp(-4*dt)
            dl=st.atm[s]*supply
            st.atm[s]-=dl; st.escaped_kg_m2[s]+=dl

    # H escape is limited by hydrogen actually delivered through the cold trap.
    heavy_moles=sum(st.upper_atm[s]/MM[s] for s in SPECIES if s not in ("H2","H2O","He"))
    h_moles=st.upper_atm["H2"]/MM["H2"]+2*st.upper_atm["H2O"]/MM["H2O"]
    fH=clamp(h_moles/(h_moles+heavy_moles+1e-30))
    phi_energy=.15*xuv*R/(4*G*M)
    phi_diff=2.5e13*fH*1e4*1.6735575e-27
    hcap=min(phi_energy,phi_diff)*dt*SEC_GYR

    d=min(st.upper_atm["H2"],hcap)
    st.upper_atm["H2"]-=d; st.escaped_kg_m2["H2"]+=d; hcap-=d
    if hcap>0:
        water_h_fraction=2*1.00794e-3/MM["H2O"]
        dw=min(st.upper_atm["H2O"],hcap/water_h_fraction)
        st.upper_atm["H2O"]-=dw
        # Still tracked as H2O-equivalent ocean loss until explicit oxygen accounting is added.
        st.escaped_kg_m2["H2O"]+=dw

    hazard=.035*dt*sw*coupling*exobase**2*(EARTH_VESC/max(vesc,300))**1.8
    frac=1-math.exp(-min(hazard,20))
    for s in ("N2","CO2","CO","CH4","Ar","He"):
        lost=st.upper_atm[s]*frac
        st.upper_atm[s]-=lost; st.escaped_kg_m2[s]+=lost


def species_total(st,s):
    total=st.atm[s]+st.upper_atm[s]+st.frost[s]+st.clath[s]+st.crust_vol[s]+st.escaped_kg_m2[s]
    if s=="H2O": total+=st.water_ice+st.ocean+st.hp_ice+st.supercrit
    if s=="CH4": total+=st.methane_liquid
    return total

def inventory_snapshot(st):
    return {s:species_total(st,s) for s in SPECIES}

def assert_finite_state(st,label=""):
    vals=[]
    for d in (st.atm,st.upper_atm,st.frost,st.clath,st.crust_vol,st.escaped_kg_m2): vals.extend(d.values())
    vals += [st.water_ice,st.ocean,st.hp_ice,st.supercrit,st.methane_liquid,st.organics,
             st.refractory_carbon,st.T,st.albedo]
    if not all(math.isfinite(x) and x>=0 for x in vals):
        raise AssertionError(f"Non-finite/negative state at {label}")

# --- evolution ---
def evolve(p,epochs=128):
    st,bulk,redox=initialize(p)
    rng=random.Random(p.seed+991); phase=rng.random()
    dt=4.5/epochs
    history=[]
    print('Evolving',p.name,'for',epochs,'epochs of',dt,'Gyr each')

    # Establish the primordial surface BEFORE impact erosion. This condenses water and
    # other volatiles into ocean/ice/frost/clathrate reservoirs and builds the layer cake.
    # Only gas remaining above that surface is eligible for formation-era stripping.
    climate(st,p)
    build_column(st,p,st.T)

    # Violent accretion / giant impacts remove a compositionally neutral fraction of the
    # actual post-condensation atmosphere. Condensed and crustal reservoirs are untouched.
    formation_impact_erosion(st,p,rng)

    # Re-equilibrate after the pressure drop before beginning secular evolution.
    climate(st,p)
    build_column(st,p,st.T)

    initial_mobile=sum(st.atm.values())+sum(st.upper_atm.values())+sum(st.frost.values())+sum(st.clath.values())+sum(st.crust_vol.values())+st.water_ice+st.ocean+st.hp_ice+st.supercrit+st.methane_liquid+sum(st.escaped_kg_m2.values())
    assert_finite_state(st,"initial")
    for e in range(epochs):
        print(f"Age {e*dt:.3f} Gyr: T={st.T:.1f} K, P={total_atmosphere(st)*G*(p.mass_e*ME)/(p.radius_e*RE)**2/BAR:.3f} bar, albedo={st.albedo:.3f}, runaway={st.runaway}")
        age=(e+.5)*dt
        u=smooth(.15,1.2,age)
        flux=math.exp(math.log(max(p.formation_flux,1e-9))*(1-u)+math.log(max(p.flux,1e-9))*u)*luminosity(age)
        pe=copy.copy(p); pe.flux=flux
        climate(st,pe)
        chemistry_epoch(st,pe,age,dt)
        escape_epoch(st,pe,age,dt,phase)
        climate(st,pe)
        assert_finite_state(st,f"epoch {e}")
        current_mobile=sum(st.atm.values())+sum(st.upper_atm.values())+sum(st.frost.values())+sum(st.clath.values())+sum(st.crust_vol.values())+st.water_ice+st.ocean+st.hp_ice+st.supercrit+st.methane_liquid+sum(st.escaped_kg_m2.values())
        # Photochemical bookkeeping is approximate but may never create orders of magnitude of mass.
        if current_mobile > initial_mobile*1.35 + 1e-6:
            raise AssertionError(f"Mass growth >35% at epoch {e}: {current_mobile/initial_mobile:.3g}x")
        if e in (0,epochs//4,epochs//2,3*epochs//4,epochs-1):
            history.append((age,st.T,total_atmosphere(st)*G*(p.mass_e*ME)/(p.radius_e*RE)**2/BAR,
                            st.silicate_state.value,st.melt_fraction))
    # final current-star solve
    climate(st,p)
    return st,bulk,redox,history

def result_row(p,st,fv,redox):
    g=G*(p.mass_e*ME)/(p.radius_e*RE)**2
    P=g*total_atmosphere(st)/BAR
    mol={s:st.atm[s]/MM[s] for s in SPECIES}; mt=sum(mol.values())
    comp={s:(mol[s]/mt if mt else 0) for s in SPECIES}
    water_depth=(st.water_ice/920+st.ocean/1000+st.hp_ice/1300)
    return dict(name=p.name,T_K=st.T,P_bar=P,albedo=st.albedo,runaway=st.runaway,
                silicate_state=st.silicate_state.value,melt_fraction=st.melt_fraction,
                surface_stack=" > ".join(L.material.value for L in st.layers),
                volatile_shell_km=sum(L.thickness_m for L in st.layers[1:-1])/1000,
                water_shell_km=water_depth/1000,
                water_total_kg_m2=atmospheric_column(st,"H2O")+st.water_ice+st.ocean+st.hp_ice+st.supercrit,
                upper_atm_kg_m2=sum(st.upper_atm.values()),
                tropopause_T_K=st.tropopause_temperature_k,
                tropopause_P_bar=st.tropopause_pressure_bar,
                cold_trap_H2O=st.cold_trap_h2o_mole_fraction,
                upper_T_K=st.upper_temperature_k,
                hp_ice_km=st.hp_ice/1300/1000,
                ocean_crust_coupling=st.ocean_crust_coupling,
                atmosphere_crust_coupling=st.atmosphere_crust_coupling,
                organics_m=st.organics/850,
                redox=redox,water_fraction=st.bulk_water_fraction,
                carbon_fraction=st.bulk_carbon_fraction,nitrogen_fraction=st.bulk_nitrogen_fraction,
                sulfur_fraction=st.bulk_sulfur_fraction,
                escaped_total_kg_m2=sum(st.escaped_kg_m2.values()),
                **{s:comp[s] for s in SPECIES})


def column_rows(p,st):
    return [dict(planet=p.name,geometry=("tidally_locked" if p.tidally_locked else "latitude"),
                 column=c.name,area_fraction=c.area_fraction,insolation_factor=c.insolation_factor,
                 T_K=c.temperature_k,albedo=c.albedo,cloud_fraction=c.cloud_fraction,
                 ice_fraction=c.ice_fraction,ocean_fraction=c.ocean_fraction,
                 desert_fraction=c.desert_fraction,local_water_kg_m2=c.local_water_kg_m2,
                 local_ocean_kg_m2=c.local_ocean_kg_m2,local_ice_kg_m2=c.local_ice_kg_m2,
                 absorbed_W_m2=c.absorbed_w_m2,surface_RH=c.surface_relative_humidity,
                 precipitable_water_kg_m2=c.precipitable_water_kg_m2,
                 cloud_base_m=c.cloud_base_m,cloud_top_m=c.cloud_top_m,
                 cloud_water_kg_m2=c.cloud_water_kg_m2,
                 methane_cloud_kg_m2=c.methane_cloud_kg_m2,
                 methane_cloud_base_m=c.methane_cloud_base_m,
                 co2_cloud_kg_m2=c.co2_cloud_kg_m2,co2_cloud_base_m=c.co2_cloud_base_m,
                 cold_trap_h2o_x=c.cold_trap_h2o_x)
            for c in st.climate_columns]
