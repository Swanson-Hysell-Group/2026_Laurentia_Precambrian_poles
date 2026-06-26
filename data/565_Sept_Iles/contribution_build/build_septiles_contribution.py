"""Build a MagIC contribution for the Sept-Iles layered mafic intrusion from
Tanczyk, Lapointe, Morris & Schmidt (1987), CJES 24, 1431-1438
(doi:10.1139/e87-135). Site-mean directions are transcribed from Table 1.

Two remanences (the intrusion is layered subhorizontally and treated in situ):
  Remanence A -- PRIMARY thermochemical remanence of the gabbro+anorthosite,
    acquired on initial cooling of the intrusion (ca. 565 Ma): D=333/I=-29,
    magnetite, unblocks 550-580 C. 10 sites. Pole 20 N / 141 E (= -20 / 321).
  Remanence B -- SECONDARY, carried by the cross-cutting diabase dykes (4 sites)
    and the remagnetized gabbro/anorthosite host adjacent to the dyke contacts
    (12 sites): D=188/I=-85, pole 59 S / 116 E. The inverse baked-contact test
    (host retains A away from the dykes, remagnetized to B at the contacts) shows
    A predates the dykes and is primary.
"""
import os
import pandas as pd
import pmagpy.pmag as pmag

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.dirname(HERE)
LAT, LON = 50.2, 293.5
N = None  # NaN marker

# Table 1 Remanence A (gabbro+anorthosite, primary): site, D, I, n, k, a95
A = [('1',320,-34,13,61,5),('3',336,-40,9,61,7),('4a',331,-41,3,16,31),
     ('8a',334,-21,5,21,17),('10a',342,-24,3,57,16),('11a',326,-39,6,100,8),
     ('13a',323,-30,4,28,18),('14',330,-10,5,173,6),('15a',344,-38,2,None,None),
     ('17a',343,-10,1,None,None)]
# Table 1 Remanence B in the host gabbro/anorthosite (secondary, near dyke contacts)
Bhost = [('4b',256,-79,4,22,20),('5',174,-88,2,None,None),('7',347,65,6,126,6),
         ('8b',31,83,2,None,None),('9',38,71,6,10,22),('10b',254,-78,4,30,17),
         ('11b',254,-88,2,None,None),('13b',263,-88,2,None,None),('15b',25,85,6,23,14),
         ('16',349,76,4,94,10),('17b',21,79,4,27,18),('19',73,-73,4,80,10)]
# Table 1 Remanence B in the diabase dykes (sites 2,6,12,18)
Bdyke = [('2',8,80,8,32,10),('6',262,89,8,18,13),('12',28,-83,8,32,10),('18',19,77,6,25,14)]

rows = []
def add(site, comp, dec, inc, n, k, a95, geol):
    plon, plat, dp, dm = pmag.dia_vgp(dec, inc, a95 if a95 else 5.0, LAT, LON)
    rows.append({'site': f'SI{site}', 'location': 'Sept-Iles intrusion', 'dir_comp_name': comp,
                 'result_type': 'i', 'result_quality': 'g',
                 'method_codes': 'LP-DIR-T:LP-DIR-AF:DE-BFL:DE-FM:DE-VGP',
                 'citations': '10.1139/e87-135', 'geologic_classes': 'Intrusive',
                 'geologic_types': geol, 'lithologies': 'Gabbro' if 'dyke' not in comp else 'Diabase',
                 'lat': LAT, 'lon': LON, 'dir_tilt_correction': 0, 'dir_dec': dec, 'dir_inc': inc,
                 'dir_k': k, 'dir_alpha95': a95, 'dir_n_samples': n,
                 'vgp_lat': round(plat, 1), 'vgp_lon': round(plon, 1),
                 'age': 565, 'age_low': 561, 'age_high': 569, 'age_unit': 'Ma'})

for s, d, i, n, k, a in A:      add(s, 'A', d, i, n, k, a, 'Layered Intrusion')
for s, d, i, n, k, a in Bhost:  add(s, 'B-host', d, i, n, k, a, 'Layered Intrusion')
for s, d, i, n, k, a in Bdyke:  add(s, 'B-dyke', d, i, n, k, a, 'Volcanic Dike')
sites = pd.DataFrame(rows)

# pole from the 10 Remanence-A site VGPs
Apole = pmag.fisher_mean(__import__('pmagpy.ipmag', fromlist=['ipmag']).make_di_block(
    sites[sites.dir_comp_name=='A']['vgp_lon'].tolist(), sites[sites.dir_comp_name=='A']['vgp_lat'].tolist()))

locs = pd.DataFrame([
    {'location':'Sept-Iles intrusion','location_type':'Region',
     'result_name':'Sept-Iles A component (primary) ca. 565 Ma pole','result_type':'a',
     'sites':':'.join(sites[sites.dir_comp_name=='A']['site']),
     'method_codes':'LP-DIR-T:DE-BFL:DE-FM:DE-VGP:ST-C-I','citations':'10.1139/e87-135:10.1086/516033:10.1130/G36247.1',
     'geologic_classes':'Intrusive','lithologies':'Gabbro','lat_s':LAT,'lat_n':LAT,'lon_w':LON,'lon_e':LON,
     'age':565,'age_low':561,'age_high':569,'age_unit':'Ma','dir_tilt_correction':0,
     'pole_lat':round(Apole['inc'],1),'pole_lon':round(Apole['dec'],1),
     'pole_alpha95':round(Apole['alpha95'],1),'pole_k':round(Apole['k'],1),'pole_n_sites':int(Apole['n']),
     'description':'Primary thermochemical remanence A of the gabbro/anorthosite, acquired on cooling ca. 565 Ma. Tanczyk et al. (1987) D=333/I=-29, pole 20 N/141 E. Positive inverse baked-contact test vs. the cross-cutting dykes. Bono & Tarduno (2015) single-crystal work confirms only this shallow component is primary (single-domain).'},
    {'location':'Sept-Iles intrusion','location_type':'Region',
     'result_name':'Sept-Iles B component (dyke / contact secondary) pole','result_type':'a',
     'sites':':'.join(sites[sites.dir_comp_name=='B-dyke']['site']),
     'method_codes':'LP-DIR-T:DE-BFL:DE-FM:DE-VGP','citations':'10.1139/e87-135',
     'geologic_classes':'Intrusive','lithologies':'Diabase','lat_s':LAT,'lat_n':LAT,'lon_w':LON,'lon_e':LON,
     'age':540,'age_low':500,'age_high':565,'age_unit':'Ma','dir_tilt_correction':0,
     'pole_lat':-59.0,'pole_lon':116.0,'pole_alpha95':10.0,'pole_k':81.0,'pole_n_sites':4,
     'description':'Steep secondary remanence B carried by the diabase dykes and the remagnetized host at dyke contacts: D=188/I=-85, pole 59 S/116 E. Tanczyk et al. (1987) related it to the dyke event; Bono & Tarduno (2015) reinterpret the steep direction as a young, soft/multidomain overprint that fails a reversal test. Not the Sept-Iles pole.'},
])

def write_magic(df, kind, path):
    with open(path,'w') as f: f.write(f'tab\t{kind}\n')
    df.to_csv(path, sep='\t', index=False, mode='a')
write_magic(sites,'sites',os.path.join(OUT,'sites.txt'))
write_magic(locs,'locations',os.path.join(OUT,'locations.txt'))
print(f"-I- wrote sites.txt ({len(sites)}: A={sum(sites.dir_comp_name=='A')}, "
      f"B-host={sum(sites.dir_comp_name=='B-host')}, B-dyke={sum(sites.dir_comp_name=='B-dyke')})")
print(f"   A pole: {Apole['inc']:.1f}/{Apole['dec']:.1f} A95 {Apole['alpha95']:.1f} K {Apole['k']:.1f} N {int(Apole['n'])}")
print("   compilation A pole: -20/321 (= 20 N/141 E), dp5/dm9, A95 6.7, B=10")
