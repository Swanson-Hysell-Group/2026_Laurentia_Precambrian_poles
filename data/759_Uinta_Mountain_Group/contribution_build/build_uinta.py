"""Build the Uinta Mountain Group MagIC contribution (sites.txt + locations.txt).

The pole notebook (``pole_notebooks/759_Uinta.ipynb``) is hand-maintained and
edited directly; there is no notebook builder.

Source: a student MagIC contribution (id 20680) for Weil, Geissman & Ashby
(2006), Precambrian Research 147, 234-259 (doi:10.1016/j.precamres.2006.01.017),
audited against the paper.

Audited fixes (instructor review): site result_type a -> i; the eastern
overprint localities (Irish Canyon, Cross Mountain, Lone Mountain, Juniper
Mountain), which carry only a north-directed recent VRM (poles at 80-88 N), are
relabeled as "UMG present-field overprint" rather than the headline UMG pole;
primary-locality result_name updated to the ca. 759 Ma UMG pole.

Weil et al. (2006) UMG ChRM pole: 0.8 N, 161.3 E, a95 4.6, N=9 sampling
localities (79 sites). Hematite-cemented sandstone/quartzite; dual-polarity ChRM
(west = normal, east = reverse), shallow. Age: late Tonian (~730-766 Ma; CA-ID-TIMS
maximum depositional age 766.3 +/- 0.5 Ma, Dehler et al., 2023; upper UMG
correlated with the upper Chuar Group, ~730 Ma).
"""
import os, io
import pandas as pd
import pmagpy.pmag as pmag, pmagpy.ipmag as ipmag

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.dirname(HERE)
SRC = os.path.join(HERE, 'Weil2006_UMG_magic_20680_source.txt')

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

OVERPRINT = ['Irish Canyon Section', 'Cross Mountain', 'Lone Mountain', 'Juniper Mountain']
# relabel overprint vs primary
locs.loc[locs['location'].isin(OVERPRINT), 'result_name'] = 'UMG present-field overprint'
prim_mask = (~locs['location'].isin(OVERPRINT)) & (~locs['location'].str.contains('Uinta Mountain Group ca', na=False))
locs.loc[prim_mask | locs['location'].str.contains('Uinta Mountain Group ca', na=False),
         'result_name'] = 'Uinta Mountain Group ca. 759 Ma pole'
# site result_type a -> i
if 'result_type' in sites.columns:
    sites['result_type'] = 'i'

def write_magic(df, kind, path):
    with open(path, 'w') as f:
        f.write(f'tab\t{kind}\n')
    df.to_csv(path, sep='\t', index=False, mode='a')

write_magic(sites, 'sites', os.path.join(OUT, 'sites.txt'))
write_magic(locs, 'locations', os.path.join(OUT, 'locations.txt'))

# the 9 primary locality poles (each one VGP for the grand mean)
prim = locs[prim_mask].copy()
blk = ipmag.make_di_block(prim['pole_lon'].tolist(), prim['pole_lat'].tolist())
p = pmag.fisher_mean(blk)
print(f'-I- wrote sites.txt ({len(sites)}), locations.txt; UMG locality-mean pole '
      f'{p["inc"]:.1f}/{p["dec"]:.1f} A95 {p["alpha95"]:.1f} N {int(p["n"])} ({len(prim)} localities). '
      f'The pole notebook recomputes the preferred site-level pole.')
