"""Reusable CV + formatting audit for a MagIC combined upload file."""
import json, glob, re, sys, io
import pandas as pd

REPO = '/Users/unimos/0000_Github/2026_Laurentia_Precambrian_poles'
dm = json.load(open(glob.glob(f'{REPO}/resources/MagIC Data Model*.json')[0], encoding='utf-8-sig'))['tables']
cv = json.load(open(f'{REPO}/resources/EarthRef Controlled Vocabularies.json', encoding='utf-8-sig'))
mcj = json.load(open(f'{REPO}/resources/MagIC Method Codes.json', encoding='utf-8-sig'))
mcodes = set()
for k, v in mcj.items():
    if isinstance(v, dict) and 'codes' in v:
        for c in v['codes']:
            mcodes.add(c.get('code'))
def cvitems(n): return {x['item'] for x in cv[n]['items']}

def audit(path):
    raw = open(path, encoding='utf-8', errors='replace').read()
    print('=== ', path.split('/')[-1], ' ===')
    print('mojibake Â:', 'Â' in raw, '| double-quote:', '"' in raw)
    markers = [l for l in raw.split('\n') if l.startswith('tab')]
    print('markers:', markers)
    parts = raw.split('>>>>>>>>>>')
    def load(n):
        for p in parts:
            L = p.strip().split('\n')
            if L[0].split('\t')[1].strip() == n:
                return pd.read_csv(io.StringIO('\n'.join(L[1:])), sep='\t', dtype=str).fillna('')
    problems = []
    for tname in ['locations', 'sites', 'ages']:
        df = load(tname)
        if df is None: continue
        tcols = dm[tname]['columns']
        for col in df.columns:
            if col not in tcols:
                problems.append(f'{tname}.{col}: NOT IN DATA MODEL'); continue
            vals = tcols[col].get('validations', []) or []
            cvref = [re.search(r'cv\("([^"]+)"\)', v).group(1) for v in vals if v.startswith('cv(')]
            is_list = tcols[col].get('type') == 'List'
            if col == 'method_codes': allowed = mcodes
            elif cvref: allowed = cvitems(cvref[0])
            else: continue
            bad = set()
            for cell in df[col]:
                if cell == '': continue
                for t in (cell.split(':') if is_list else [cell]):
                    if t.strip() not in allowed: bad.add(t.strip())
            if bad: problems.append(f'{tname}.{col} [{cvref[0] if cvref else "methodcodes"}]: INVALID {sorted(bad)}')
        # required cols
        req = [c for c, d in tcols.items() if 'required()' in (d.get('validations') or [])]
        miss = [c for c in req if c not in df.columns or (df[c] == '').any()]
        if miss: problems.append(f'{tname}: MISSING/blank required {miss}')
    print('CV/required PROBLEMS:', problems or 'NONE — data-model clean')
    return not problems

if __name__ == '__main__':
    audit(sys.argv[1])
