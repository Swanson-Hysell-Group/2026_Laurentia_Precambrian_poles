"""Build the Franklin LIP MagIC contribution (data/719_Franklin_LIP/) from the
audited Cortopassi / Denyszyn et al. (2009) student contribution.

Source: Denyszyn, Halls, Davis & Evans (2009), "Paleomagnetism and U-Pb
geochronology of Franklin dykes, High Arctic Canada and Greenland", CJES 46,
689-705 (doi:10.1139/E09-042). The student also compiled the Franklin grand mean
(adc9430 = Pu et al., 2022 geochronology).

Audited fixes (instructor review): NU1 site vgp_lat -6.6 -> -16.6 (Denyszyn 2009
Table 1, reverse-polarity site); grand-mean citation DOI sciadv.adc9431 ->
adc9430. The repo's 719_Franklin_LIP notebook already reports the prior-compilation
Franklin event grand mean (6.7 N / 162.1 E, A95 3.0, B=56); this contribution adds
the 27 site-level Franklin dykes (15 Arctic Canada + 12 Greenland) for the
site-level workflow.
"""
import os, io
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)
SRC = os.path.join(HERE, 'Denyszyn2009_Franklin_magic_source.txt')

data = open(SRC, encoding='latin-1').read()
sites = locs = None
for b in data.split('>>>>>>>>>>'):
    b = b.strip()
    if not b:
        continue
    kind = b.splitlines()[0].split('\t')[1].strip()
    df = pd.read_csv(io.StringIO('\n'.join(b.splitlines()[1:])), sep='\t')
    if kind == 'sites':
        sites = df
    elif kind == 'locations':
        locs = df

# NU1 reverse-polarity site vgp_lat transcription fix
if 'vgp_lat' in sites.columns:
    sites.loc[sites['site'] == 'NU1', 'vgp_lat'] = -16.6

def write_magic(df, kind, path):
    with open(path, 'w') as f:
        f.write(f'tab\t{kind}\n')
    df.to_csv(path, sep='\t', index=False, mode='a')

write_magic(sites, 'sites', os.path.join(OUT, 'sites.txt'))
write_magic(locs, 'locations', os.path.join(OUT, 'locations.txt'))
print(f'-I- wrote sites.txt ({len(sites)}), locations.txt ({len(locs)} poles)')
