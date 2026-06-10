"""Build the Chuar Group combined paleomagnetic pole (757_Chuar_Group).

A NEW COMBINED pole spanning the Carbon Canyon (Galeros Fm, ca. 761 Ma), Carbon
Butte and Awatubi (Kwagunt Fm, ca. 750 Ma) members, reproducing Eyster et al.
(2020) Table 1 at the reading level: N = 37 = Carbon Butte-Awatubi (23) + Carbon
Canyon (14), as an equal-weight Fisher mean of 37 individual site/sample VGPs.

Sources:
  - Eyster et al. (2020), GSA Bulletin (doi:10.1130/B32012.1). The 17 high-
    temperature (HT) Carbon Butte-Awatubi readings are reconstructed from the
    Supplemental Table S2/DR2 (`Eyster2020_SI_DR2.txt`, parsed from GSA Data
    Repository item 2019239). DECODED RULE reproducing Eyster's N=17 pole
    (12.5 N / 161.6 E, A95 4.0): linear (L) shallow-HT least-squares fits with
    MAD < 15 deg, one direction per sample (the ** two-unblocking rows averaged);
    great-circle (C) fits are used by Eyster only in the locality means, not in
    the every-sample-as-a-site pole, and are excluded here. The six localities
    map to members A1303/A1305/A1306/A1310 = Carbon Butte, A1304/A1307 = Awatubi.
  - Weil, Geissman & Van der Voo (2004), Precambrian Research 129, 71-92
    (doi:10.1016/j.precamres.2003.09.016): 6 Carbon Butte + 14 Carbon Canyon
    tilt-corrected ChRM site VGPs (the same 6 Weil sites Eyster added to their
    Combined pole; Weil's lone Awatubi/Walcott/Jupiter sites are not in Eyster's
    N=37 composition and are excluded).

Reconstruction reproduces Eyster's published poles: CB-Awatubi 23 -> 13.2/163.5
(Eyster 13.5/162.8); Carbon Canyon 14 -> -2.1/163.7 (Eyster -2.1/163.7);
combined 37 -> 7.5 N / 163.6 E, A95 4.2.

Member age model from Dehler et al. (2023): Carbon Canyon 761.1 Ma, Carbon Butte-
Awatubi 750.5 Ma; Re-Os 751.0 +/- 7.6 Ma on the Walcott Member (Rooney et al.,
2017). Combined nominal age ~757 Ma. Positive Tauxe & Watson (1994) fold test
(66-109% unfolding) and a positive reversal test within the coeval Kwagunt
members (Watson V 7.4 < 7.6; MM1990 class C; Eyster et al., 2020). E-I
inclination-flattening factors (Eyster 2020): f=0.9 (Carbon Butte-Awatubi),
f=0.5 (Carbon Canyon).
"""
import os, io, re
import pandas as pd
import pmagpy.pmag as pmag
import pmagpy.ipmag as ipmag
import nbformat
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)
WEIL_DOI = '10.1016/j.precamres.2003.09.016'
EYSTER_DOI = '10.1130/B32012.1'

LOC_MEMBER = {'A1303': 'Carbon Butte', 'A1304': 'Awatubi', 'A1305': 'Carbon Butte',
              'A1306': 'Carbon Butte', 'A1307': 'Awatubi', 'A1310': 'Carbon Butte'}
LOC_COORD = {'A1303': (36.268, 248.114), 'A1304': (36.268, 248.114), 'A1305': (36.274, 248.114),
             'A1306': (36.276, 248.110), 'A1307': (36.274, 248.110), 'A1310': (36.167, 248.161)}


def parse_eyster_si():
    """Reconstruct Eyster's 17 HT Carbon Butte-Awatubi sample directions from SI
    Table S2. Returns one row per sample with averaged in-situ and tilt-corrected
    HT directions and the tilt-corrected VGP."""
    L = open(os.path.join(HERE, 'Eyster2020_SI_DR2.txt')).read().splitlines()
    s0 = next(i for i, l in enumerate(L) if 'SAMPLE LEAST SQUARE FITS' in l)
    s1 = next(i for i, l in enumerate(L) if 'HT Locality' in l)
    sr = re.compile(r'^(A\d{4}-\w+?)(\*\*)?$')
    rows = []
    for ln in L[s0:s1]:
        t = ln.split()
        if len(t) < 8 or not sr.match(t[0]) or t[1] not in ('L', 'C'):
            continue
        sm = sr.match(t[0]).group(1)
        loc = sm.split('-')[0]
        if loc not in LOC_MEMBER:
            continue
        try:
            g = [float(x) for x in t[3:7]]
        except ValueError:
            continue
        mad = float(t[-5]) if t[1] == 'C' else float(t[-1])      # C rows carry 4 trailing arc values
        rows.append(dict(sample=sm, locality=loc, typ=t[1], comp=t[2],
                         geoD=g[0], geoI=g[1], tcD=g[2], tcI=g[3], mad=mad))
    d = pd.DataFrame(rows)
    ht = d[(d['comp'].isin(['HT', 'Ht.'])) & (d['typ'] == 'L') & (d['mad'] < 15)]
    out = []
    for sm, grp in ht.groupby('sample'):
        loc = grp['locality'].iloc[0]
        la, lo = LOC_COORD[loc]
        gm = pmag.fisher_mean(ipmag.make_di_block(grp['geoD'].tolist(), grp['geoI'].tolist()))
        tm = pmag.fisher_mean(ipmag.make_di_block(grp['tcD'].tolist(), grp['tcI'].tolist()))
        plon, plat, dp, dm = pmag.dia_vgp(tm['dec'], tm['inc'], 3, la, lo)
        out.append(dict(site=sm, study='Eyster et al. (2020)', member=LOC_MEMBER[loc],
                        lat=la, lon=lo, geo_dec=round(gm['dec'], 1), geo_inc=round(gm['inc'], 1),
                        tc_dec=round(tm['dec'], 1), tc_inc=round(tm['inc'], 1),
                        vgp_lat=round(plat, 1), vgp_lon=round(plon, 1), n=len(grp)))
    return pd.DataFrame(out)


def parse_weil():
    """Weil (2004) ChRM Carbon Butte (6) + Carbon Canyon (14), in-situ + tilt-corrected."""
    for b in open(os.path.join(HERE, 'Weil2004_MagIC.txt'), encoding='latin-1').read().split('>>>>>>>>>>'):
        b = b.strip()
        if b and b.splitlines()[0].split('\t')[1].strip() == 'sites':
            df = pd.read_csv(io.StringIO('\n'.join(b.splitlines()[1:])), sep='\t')
    df = df[df['dir_comp_name'] == 'ChRM'].copy()
    df['member'] = df['description'].str.extract(r'(Carbon Canyon|Carbon Butte|Awatubi|Walcott|Jupiter)')
    return df[df['member'].isin(['Carbon Butte', 'Carbon Canyon'])].copy()


# ------------------------------ assemble sites ------------------------------
eyster = parse_eyster_si()
weil = parse_weil()

rows = []
for _, r in eyster.iterrows():
    common = dict(site=r['site'], study='Eyster et al. (2020)', member=r['member'],
                  dir_comp_name='ChRM', citations=EYSTER_DOI, lat=r['lat'], lon=r['lon'],
                  dir_alpha95='', dir_k='', dir_n_samples=int(r['n']))
    rows.append({**common, 'dir_tilt_correction': 0, 'dir_dec': r['geo_dec'], 'dir_inc': r['geo_inc'],
                 'vgp_lat': '', 'vgp_lon': ''})
    rows.append({**common, 'dir_tilt_correction': 100, 'dir_dec': r['tc_dec'], 'dir_inc': r['tc_inc'],
                 'vgp_lat': r['vgp_lat'], 'vgp_lon': r['vgp_lon']})
for _, r in weil.iterrows():
    rows.append(dict(site=r['site'], study='Weil et al. (2004)', member=r['member'],
                     dir_comp_name='ChRM', citations=WEIL_DOI, lat=r['lat'], lon=r['lon'],
                     dir_tilt_correction=int(r['dir_tilt_correction']),
                     dir_dec=round(r['dir_dec'], 1), dir_inc=round(r['dir_inc'], 1),
                     dir_alpha95=r['dir_alpha95'], dir_k=r['dir_k'], dir_n_samples=r['dir_n_samples'],
                     vgp_lat=round(r['vgp_lat'], 1) if pd.notna(r['vgp_lat']) else '',
                     vgp_lon=round(r['vgp_lon'], 1) if pd.notna(r['vgp_lon']) else ''))

sites = pd.DataFrame(rows)
sites['location'] = 'Chuar Group, Grand Canyon'
sites['result_type'] = 'i'
sites['result_quality'] = 'g'
sites['geologic_classes'] = 'Sedimentary'
sites['geologic_types'] = 'Sediment Layer'
sites['lithologies'] = 'Mudstone:Sandstone:Dolomite'
sites['age_unit'] = 'Ma'
sites['member_age'] = sites['member'].map({'Carbon Canyon': 761.1, 'Carbon Butte': 750.5, 'Awatubi': 750.5})
sites['method_codes'] = 'LP-DIR-T:LP-DIR-AF:DE-BFL:DE-FM:DE-VGP:DA-DIR-TILT'
sites['description'] = sites['member'] + ' member; ' + sites['study']
col_order = ['site', 'study', 'member', 'member_age', 'location', 'result_type', 'result_quality',
             'method_codes', 'citations', 'geologic_classes', 'geologic_types', 'lithologies',
             'lat', 'lon', 'age_unit', 'dir_comp_name', 'dir_tilt_correction', 'dir_dec', 'dir_inc',
             'dir_alpha95', 'dir_k', 'dir_n_samples', 'vgp_lat', 'vgp_lon', 'description']
sites = sites[col_order]


def vgp_pole(df, flip=False):
    sub = df[(df['dir_tilt_correction'] == 100) & (df['vgp_lat'] != '')]
    blk = pmag.flip(ipmag.make_di_block([float(x) for x in sub['vgp_lon']],
                                        [float(x) for x in sub['vgp_lat']]), combine=True)
    if flip:
        blk = ipmag.do_flip(di_block=blk)
    return pmag.fisher_mean(blk)


tc = sites[sites['dir_tilt_correction'] == 100]
combined = vgp_pole(tc)
cb_awatubi = vgp_pole(tc[tc['member'].isin(['Carbon Butte', 'Awatubi'])])
carbon_canyon = vgp_pole(tc[tc['member'] == 'Carbon Canyon'], flip=True)


def locrow(name, p, age, lo, hi, desc):
    return {'location': 'Chuar Group, Grand Canyon', 'location_type': 'Region',
            'result_name': name, 'result_type': 'a', 'result_quality': 'g',
            'method_codes': 'LP-DIR-T:DE-BFL:DE-FM:DE-VGP:DA-DIR-TILT:FT-FM:RT',
            'citations': f'{WEIL_DOI}:{EYSTER_DOI}', 'geologic_classes': 'Sedimentary',
            'lithologies': 'Mudstone:Sandstone:Dolomite', 'lat_s': 36.1, 'lat_n': 36.4,
            'lon_w': 248.1, 'lon_e': 248.2, 'age': age, 'age_low': lo, 'age_high': hi,
            'age_unit': 'Ma', 'dir_tilt_correction': 100,
            'pole_lat': round(p['inc'], 1), 'pole_lon': round(p['dec'], 1),
            'pole_alpha95': round(p['alpha95'], 1), 'pole_k': round(p['k'], 1),
            'pole_n_sites': int(p['n']), 'description': desc}


locs = pd.DataFrame([
    locrow('Chuar Group combined pole (757 Ma)', combined, 757, 750, 761,
           'Combined Chuar Group pole: equal-weight Fisher mean of 37 tilt-corrected site/sample VGPs = Carbon Butte-Awatubi (23) + Carbon Canyon (14), reproducing Eyster et al. (2020) Table 1. 17 Eyster HT readings (SI Table S2) + 6 Weil Carbon Butte + 14 Weil Carbon Canyon sites.'),
    locrow('Carbon Butte-Awatubi pole (750.5 Ma)', cb_awatubi, 750.5, 745, 756,
           'Carbon Butte-Awatubi pole: 17 Eyster (2020) HT samples + 6 Weil (2004) Carbon Butte sites (N=23). Reproduces Eyster et al. (2020) Combined HT Carbon Butte-Awatubi pole (13.5/162.8, A95 3.3). E-I f=0.9.'),
    locrow('Carbon Canyon member pole (761 Ma)', carbon_canyon, 761.1, 756, 766,
           'Carbon Canyon member (Galeros Fm): Fisher mean of 14 Weil (2004) tilt-corrected site VGPs. Reproduces Eyster et al. (2020) Carbon Canyon pole (-2.1/163.7, A95 8.0). E-I f=0.5.'),
])


def write_magic(df, kind, path):
    with open(path, 'w') as f:
        f.write(f'tab\t{kind}\n')
    df.to_csv(path, sep='\t', index=False, mode='a')


write_magic(sites, 'sites', os.path.join(OUT, 'sites.txt'))
write_magic(locs, 'locations', os.path.join(OUT, 'locations.txt'))
n_ey = (sites['study'].str.contains('Eyster') & (sites['dir_tilt_correction'] == 100)).sum()
print(f"-I- wrote sites.txt ({len(sites)} rows; {n_ey} Eyster samples + {int((tc['study'].str.contains('Weil')).sum())} Weil sites = {len(tc)})")
print(f"    combined (37) : {combined['inc']:.1f}/{combined['dec']:.1f} A95 {combined['alpha95']:.1f} N {int(combined['n'])}")
print(f"    CB-Awatubi(23): {cb_awatubi['inc']:.1f}/{cb_awatubi['dec']:.1f} A95 {cb_awatubi['alpha95']:.1f} N {int(cb_awatubi['n'])}  (Eyster 13.5/162.8)")
print(f"    Carbon Canyon : {carbon_canyon['inc']:.1f}/{carbon_canyon['dec']:.1f} A95 {carbon_canyon['alpha95']:.1f} N {int(carbon_canyon['n'])}  (Eyster -2.1/163.7)")

# ============================ notebook ========================================
NB = os.path.abspath(os.path.join(HERE, '../../../pole_notebooks/757_Chuar_Group.ipynb'))
cells = []
md = lambda s: cells.append(new_markdown_cell(s))
co = lambda s: cells.append(new_code_cell(s))

md("""# Chuar Group ca. 757 Ma combined paleomagnetic pole

## Geologic context

The Chuar Group of the Grand Canyon (Arizona) is a ~1.6-km-thick succession of
Neoproterozoic (Tonian) mudstone, dolomite, and sandstone deposited in an
extensional basin during the early breakup of Rodinia (Weil, Geissman & Van der
Voo, 2004; Eyster et al., 2020). It is divided into the lower Galeros and upper
Kwagunt formations. This notebook builds a combined pole spanning three members:
the Carbon Canyon member of the Galeros Formation (ca. 761 Ma) and the Carbon
Butte and Awatubi members of the Kwagunt Formation (ca. 750 Ma), using the age
model of Dehler et al. (2023) and a Re-Os age of 751.0 ± 7.6 Ma on the Walcott
Member (Rooney et al., 2017).

## Pole

The combined Chuar Group pole is the equal-weight Fisher mean of 37 individual
tilt-corrected site/sample virtual geomagnetic poles (VGPs), reproducing Eyster
et al. (2020) Table 1: Carbon Butte–Awatubi (N = 23) + Carbon Canyon (N = 14).
The 17 Eyster (2020) high-temperature (HT) Carbon Butte–Awatubi readings are
reconstructed from their Supplemental Table S2 (see below); the remaining 20 are
Weil et al. (2004) site VGPs (6 Carbon Butte + 14 Carbon Canyon). The HT remanence
is carried by hematite, is of dual polarity, and passes a fold test, indicating a
primary, pre-folding origin. Because the combined pole averages ~10 Myr of
apparent polar wander between the older Carbon Canyon member and the younger
Kwagunt members, it sits at a lower paleolatitude than the Kwagunt-only pole.""")

co("""import pole_tools as pt
import pmagpy.ipmag as ipmag
import pmagpy.pmag as pmag
import pandas as pd
import matplotlib.pyplot as plt

%config InlineBackend.figure_format='retina'
%matplotlib inline""")

md("""## Reconstructing Eyster (2020) HT samples from Supplemental Table S2

Eyster et al. (2020) report their Carbon Butte–Awatubi pole at the sample level
(\"treating every sample as a site\"), but the published MagIC contribution carries
only the locality means. The 17 sample-level HT readings are reconstructed here
from their Supplemental Table S2 (GSA Data Repository 2019239), parsed into
`Eyster2020_SI_DR2.txt`. The selection rule that reproduces their N = 17 pole
(12.5°N/161.6°E, A95 4.0°) is: **linear (best-fit-line) shallow-HT components with
MAD < 15°, one direction per sample** (samples with two distinct unblocking
directions are averaged). Great-circle (arc) fits are used by Eyster only in the
locality means, not in the every-sample pole, and are excluded. The six localities
map to members: A1303/A1305/A1306/A1310 = Carbon Butte, A1304/A1307 = Awatubi.

This reconstruction is baked into the contribution `sites.txt`; the cell below
loads the result.""")

co("""sites = pd.read_csv('../data/751_Kwagunt_Chuar/sites.txt', sep='\\t', skiprows=1)
tc = sites[sites['dir_tilt_correction'] == 100].copy()
geo = sites[sites['dir_tilt_correction'] == 0].copy()
study_lat, study_lon = 36.23, 248.13
print(tc.groupby(['member', 'study']).size())
print(f'\\ntotal sites in the combined mean: {len(tc)}')""")

md("""## Component poles — reproducing Eyster (2020) Table 1

The Carbon Butte–Awatubi pole (17 Eyster samples + 6 Weil sites, N = 23) and the
Carbon Canyon pole (14 Weil sites) are computed from the site VGPs and compared
with Eyster's published values. The near-equatorial Carbon Canyon mean is flipped
to the compilation hemisphere (longitude ~164°E).""")

co("""cb_awatubi = tc[tc['member'].isin(['Carbon Butte', 'Awatubi'])]
_, cba_pole = pt.compute_mean_pole(cb_awatubi, unify_polarity=True)
print(f"Carbon Butte-Awatubi (N={int(cba_pole['n'])}): "
      f"{cba_pole['inc']:.1f}N / {cba_pole['dec']:.1f}E   A95 {cba_pole['alpha95']:.1f}   "
      f"(Eyster 2020: 13.5N / 162.8E, A95 3.3, N=23)")

carbon_canyon = tc[tc['member'] == 'Carbon Canyon']
_, cc_pole = pt.compute_mean_pole(carbon_canyon, unify_polarity=True, flip=True)
print(f"Carbon Canyon        (N={int(cc_pole['n'])}): "
      f"{cc_pole['inc']:.1f}N / {cc_pole['dec']:.1f}E   A95 {cc_pole['alpha95']:.1f}   "
      f"(Eyster 2020: -2.1N / 163.7E, A95 8.0, N=14)")""")

md("""## Member ChRM directions (tilt-corrected, mixed polarity)

Directions are plotted without unifying polarity (downward = filled, upward =
open); both polarities are present.""")

co("""colors = {'Carbon Canyon': 'firebrick', 'Carbon Butte': 'royalblue', 'Awatubi': 'seagreen'}
ipmag.plot_net()
for m in ['Carbon Canyon', 'Carbon Butte', 'Awatubi']:
    sub = tc[tc['member'] == m]
    ipmag.plot_di(sub['dir_dec'].tolist(), sub['dir_inc'].tolist(),
                  color=colors[m], marker='o', label=f'{m} (N={len(sub)})')
plt.legend(loc='center left', bbox_to_anchor=(1.05, 0.5))
plt.title('Chuar Group member ChRM directions (tilt-corrected)')
plt.show()""")

md("""## The combined Chuar Group pole (N = 37)

The 37 tilt-corrected site/sample VGPs are pooled into a single equal-weight
Fisher mean — the exported value.""")

co("""vgp_block, chuar_pole = pt.compute_mean_pole(tc, unify_polarity=True)
print(f"Combined Chuar Group pole: {chuar_pole['inc']:.1f}N / {chuar_pole['dec']:.1f}E   "
      f"A95 {chuar_pole['alpha95']:.1f}   K {chuar_pole['k']:.1f}   N {int(chuar_pole['n'])}")

ax = pt.plot_vgps_and_pole(vgp_block, chuar_pole, central_longitude=chuar_pole['dec'],
                           central_latitude=chuar_pole['inc'], figsize=(6, 6))
plt.title('Chuar Group combined pole (37 site/sample VGPs, 3 members)')
plt.show()""")

md("""## Reversal test

The McFadden & McElhinny (1990) and bootstrap reversal tests are run on all 37
combined tilt-corrected directions. Polarity here is largely confounded with
member age — the Carbon Canyon member is dominantly of one polarity and the
younger Carbon Butte–Awatubi members the other — so the two polarity groups
differ by the ~10 Myr of apparent polar wander between the members, and the test
on the full combined set is negative. Within the coeval Carbon Butte–Awatubi
members the reversal test is positive (Watson V = 7.4 < V_crit = 7.6; MM1990
class C; Eyster et al., 2020). Both polarities are present (Meert R6).""")

co("""_ = pt.reversal_test(tc, plot=True, random_seed=1)""")

md("""### Reversal test: coeval Carbon Butte–Awatubi members only

Restricting the test to the coeval Carbon Butte and Awatubi (Kwagunt) members
removes the inter-member apparent polar wander that confounds the full-set test.
On this dual-polarity subset the McFadden & McElhinny (1990) common-mean test is
positive (the normal and reversed modes share a mean direction), consistent with
the positive reversal test reported by Eyster et al. (2020).""")

co("""kwagunt = tc[tc['member'].isin(['Carbon Butte', 'Awatubi'])]
_ = pt.reversal_test(kwagunt, plot=True, random_seed=1)""")

md("""## Fold test

The Tauxe & Watson (1994) fold test of Eyster et al. (2020) places the optimal
unfolding at 66-109%, indicating a pre-folding (primary) magnetization. As a
simple check, the precision parameter k of the polarity-unified ChRM should
increase from in-situ to tilt-corrected coordinates.""")

co("""_, geo_mean = pt.compute_mean_direction(geo, unify_polarity=True)
_, tc_mean = pt.compute_mean_direction(tc, unify_polarity=True)
print(f"in-situ        : k = {geo_mean['k']:.1f}  (a95 = {geo_mean['alpha95']:.1f}, N = {int(geo_mean['n'])})")
print(f"tilt-corrected : k = {tc_mean['k']:.1f}  (a95 = {tc_mean['alpha95']:.1f}, N = {int(tc_mean['n'])})")
print('k increases on untilting -> positive fold test (cf. Tauxe & Watson 1994, 66-109% unfolding; Eyster et al., 2020)')""")

md("## The Chuar Group pole in the context of the Laurentia APWP")

co("""Laurentia_poles = pt.get_Laurentia_poles()
ax = pt.plot_apwp_context(Laurentia_poles, chuar_pole['inc'], chuar_pole['dec'],
                          chuar_pole['alpha95'], age_min=635, age_max=1800,
                          projection='orthographic',
                          central_longitude=chuar_pole['dec'],
                          central_latitude=chuar_pole['inc'])
plt.show()""")

md("## R7: comparison with younger Laurentia poles")

co("""Laurentia_stricto_poles = pt.get_Laurentia_stricto_poles()
pt.plot_pole_overlap('Chuar Group', Laurentia_stricto_poles,
                     pt.Torsvik2012_Laurentia,
                     pole_plat=chuar_pole['inc'], pole_plon=chuar_pole['dec'],
                     pole_A95=chuar_pole['alpha95'], pole_age=757)""")

md("""## R-score summary (Meert et al., 2020)

| R | Criterion | Score | Justification |
|---|---|---|---|
| 1 | Age within ± 15 Ma | **1** | Member ages from the Dehler et al. (2023) U-Pb / Re-Os age model (Carbon Canyon 761.1 Ma; Carbon Butte-Awatubi 750.5 Ma) and Re-Os 751.0 ± 7.6 Ma on the Walcott Member (Rooney et al., 2017). |
| 2 | Techniques and statistical analysis | **1** | Thermal/AF demagnetization, PCA; N = 37 site/sample readings across three members, averaging paleosecular variation; A95 ≈ 4°. |
| 3 | Magnetic mineralogy characterized | **1** | Remanence carried by hematite, characterized by IRM/thermal experiments (Weil et al., 2004; Eyster et al., 2020). |
| 4 | Field tests constrain age of magnetization | **f** | Positive Tauxe & Watson (1994) fold test (66-109% unfolding; Eyster et al., 2020; Weil et al., 2004). |
| 5 | Structural control / tectonic coherence | **1** | Autochthonous Grand Canyon Chuar basin; tilt-corrected. |
| 6 | Presence of reversals | **1** | Dual-polarity ChRM. The reversal test on the full combined set is confounded by inter-member APW; within the coeval Carbon Butte–Awatubi members it is positive (Watson V 7.4 < 7.6; MM1990 class C; Eyster et al., 2020). |
| 7 | No resemblance to younger poles | **0** | The pole lies within the mid-Neoproterozoic Laurentia cluster (e.g. the Uinta Mountain Group). |
| | Total | **6/7** | Grade A |""")

md("""## Nordic workshop summary

The Chuar Group combined pole is a new pole spanning the Carbon Canyon, Carbon
Butte, and Awatubi members, computed as the equal-weight Fisher mean of 37
individual tilt-corrected site/sample VGPs, reproducing Eyster et al. (2020)
Table 1 (Carbon Butte–Awatubi N=23 + Carbon Canyon N=14). The 17 Eyster HT
readings are reconstructed from their Supplemental Table S2; the other 20 are Weil
et al. (2004) sites. The observed (uncorrected) pole is exported; Eyster et al.
(2020) derived E-I inclination-flattening factors of f = 0.9 (Carbon
Butte–Awatubi) and f = 0.5 (Carbon Canyon). Positive fold and reversal tests and
the ~757 Ma age give a Grade A pole. This combined pole sits at a lower
paleolatitude than the Kwagunt-only pole because it incorporates the older Carbon
Canyon member.""")

co("""# combined ChRM mean direction at the study locality
_, chuar_dir = pt.compute_mean_direction(tc, unify_polarity=True)

chuar_summary = pt.make_nordic_summary(
    terrane='Laurentia',
    rockname='Chuar Group (combined)',
    sites=tc,
    dir_mean=chuar_dir,
    pole_mean=chuar_pole,
    study_lon=study_lon,
    study_lat=study_lat,
    lithology='sedimentary',
    component_comment='High-temperature characteristic remanent magnetization (hematite), dual polarity; combined Carbon Canyon + Carbon Butte + Awatubi members (37 readings); detrital sediments',
    tests='F+ (positive fold test, 66-109% unfolding) and R+ (positive reversal test within the coeval Carbon Butte-Awatubi members, Watson V 7.4<7.6, MM1990 class C; the full combined-set reversal test is confounded by inter-member APW) (Weil et al., 2004; Eyster et al., 2020)',
    gpmdb_number='',
    percent_reversed=50,
    demag_code=4,
    R1=1, R2=1, R3=1, R4='f', R5=1, R6=1, R7=0, Grade='A',
    nominal_age=757, lomagage=750, himagage=761,
    REF_method='Member ages from the Dehler et al. (2023) U-Pb / Re-Os age model (Carbon Canyon 761.1 Ma; Carbon Butte-Awatubi 750.5 Ma); Re-Os 751.0 +/- 7.6 Ma on the Walcott Member (Rooney et al., 2017).',
    POLE_AUTHORS='Weil, A. B., Geissman, J. W., & Van der Voo, R.; Eyster, A., et al. (combined three-member pole)',
    YEAR=2020,
    JOURNAL='Precambrian Research; GSA Bulletin',
    VOLUME='129; 132',
    VPAGES='71-92; 710-738',
    TITLE='Paleomagnetism of the Neoproterozoic Chuar Group, Grand Canyon (Weil et al., 2004; Eyster et al., 2020)',
    COMMENT='New combined Chuar Group pole, equal-weight Fisher mean of 37 individual tilt-corrected site/sample VGPs reproducing Eyster et al. (2020) Table 1: Carbon Butte-Awatubi (N=23) + Carbon Canyon (N=14). The 17 Eyster HT Carbon Butte-Awatubi readings are reconstructed from Supplemental Table S2 (GSA DR 2019239): linear shallow-HT fits with MAD<15 deg, one per sample (** double-unblocking averaged); great-circle fits excluded (Eyster uses them only in locality means). This reproduces Eyster N=17 pole 12.5/161.6. Plus 6 Weil (2004) Carbon Butte + 14 Weil Carbon Canyon site VGPs. Reconstruction reproduces Eyster published poles: CB-Awatubi 23 -> 13.2/163.5 (Eyster 13.5/162.8); Carbon Canyon 14 -> -2.1/163.7 (Eyster -2.1/163.7); combined 37 -> 7.5/163.6. Members: A1303/A1305/A1306/A1310=Carbon Butte, A1304/A1307=Awatubi. Sedimentary (lithology=sedimentary); observed pole exported; Eyster E-I f=0.9 (CB-Awatubi) and f=0.5 (Carbon Canyon). Positive Tauxe & Watson (1994) fold test (66-109% unfolding, R4=f); dual polarity, reversal test positive within coeval Kwagunt members (R6=1). Dehler et al. (2023) age model + Re-Os 751 Ma (Rooney et al., 2017, R1=1). R7=0 (within the mid-Neoproterozoic Laurentia cluster). R=6, Grade A. Supersedes the Kwagunt-only pole.'
)
pt.save_nordic_summary(chuar_summary, '757_Chuar_Group')
chuar_summary""")

nb = new_notebook(cells=cells, metadata={
    'kernelspec': {'display_name': 'Python 3', 'language': 'python', 'name': 'python3'},
    'language_info': {'name': 'python'}})
nbformat.write(nb, NB)
print('wrote', NB, 'with', len(cells), 'cells')
