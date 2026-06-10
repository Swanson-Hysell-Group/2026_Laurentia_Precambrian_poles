"""Build the Sudbury Dikes MagIC contribution (sites.txt + locations.txt).

Source: a student MagIC contribution (id 20644) for Palmer, Merz & Hayatsu
(1977), CJES 14, 1867-1887 (doi:10.1139/e77-158), audited against the paper.
The contribution is essentially clean (instructor review found only minor items);
this script re-emits it in the repository's two-file MagIC layout, with the
student identity removed and small provenance improvements.

Two location-level results, both from Palmer et al. (1977):
- Sudbury Dike Swarm primary pole (the ca. 1238 Ma westerly, shallow "Sudbury
  dike direction"): mean dir D=265 deg, I=+2 deg; pole 168 deg W / 2.5 deg S
  (= -2.5 N, 192 E), A95 2.5, N=38; baked-contact test positive (site 90 +
  Sopher 1963; Schwarz 1977). Single (primary) polarity.
- Grenville-Front overprint pole (ESE-and-down, ca. 1000 Ma Grenvillian
  remagnetization): -3.7 N / 343.4 E, A95 6.3, N=23.

Age 1238 +/- 4 Ma is the U-Pb baddeleyite age of the Sudbury swarm
(Krogh et al., 1987; Fahrig & West, 1986), superseding the paper's K-Ar estimate.
"""
import os
import io
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)
SRC = os.path.join(HERE, 'Palmer1977_Sudbury_magic_20644_source.txt')

# --- read the multi-table MagIC source -------------------------------------
data = open(SRC, encoding='latin-1').read()
tables = {}
for b in data.split('>>>>>>>>>>'):
    b = b.strip()
    if not b:
        continue
    lines = b.splitlines()
    kind = lines[0].split('\t')[1].strip()
    tables[kind] = pd.read_csv(io.StringIO('\n'.join(lines[1:])), sep='\t')

sites = tables['sites'].copy()
locs = tables['locations'].copy()

# --- minor audited improvements (issues.md) --------------------------------
# Grenville-Front overprint pole uses Palmer's N=23, but 26 sites carry the
# "Grenville Front direction" tag; 3 ESE-down sites were not in Palmer's Table 4
# pole. We cannot identify the exact three from the paper, so leave the tag and
# document the N=23 vs 26 in the location description rather than guess.
gf_desc_add = (' Palmer et al. (1977) Table 4 pole uses N=23 of the 26 sites '
               'carrying the Grenville-Front overprint direction.')
mask = locs['result_name'].str.contains('Grenville', na=False)
locs.loc[mask, 'description'] = locs.loc[mask, 'description'].astype(str) + gf_desc_add

# --- write repository-format sites.txt / locations.txt ---------------------
def write_magic(df, kind, path):
    with open(path, 'w') as f:
        f.write(f'tab\t{kind}\n')
    df.to_csv(path, sep='\t', index=False, mode='a')

write_magic(sites, 'sites', os.path.join(OUT, 'sites.txt'))
write_magic(locs, 'locations', os.path.join(OUT, 'locations.txt'))
print(f'-I- wrote sites.txt ({len(sites)} sites) and locations.txt ({len(locs)} poles)')
print('    Sudbury Dike direction sites:',
      int((sites['dir_comp_name'] == 'Sudbury Dike direction').sum()))
