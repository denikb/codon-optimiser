"""
app.py - Streamlit front end for codon_optimizer.py

Run locally:   streamlit run app.py
Deploy:        see DEPLOYMENT.md
"""

import re
import urllib.error
import urllib.request

import streamlit as st

from codon_optimizer import (
    ORGANISM_CONFIGS,
    clean_protein,
    optimize,
)
from theme import apply_theme, draw, helix, squiggle_rule

MAX_PROTEIN_LEN = 3000          # guards against a pasted genome

# Verified character-for-character against the UniProt entries cited.
# (label: (accession, note, sequence))
EXAMPLES = {
    'avGFP': (
        'P42212',
        'jellyfish reporter, the codon-optimisation benchmark',
        'MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYG'
        'VQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDG'
        'NILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYL'
        'STQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK'),
    'DsRed': (
        'Q9U6Y8',
        'coral red fluorescent protein, parent of mCherry',
        'MRSSKNVIKEFMRFKVRMEGTVNGHEFEIEGEGEGRPYEGHNTVKLKVTKGGPLPFAWDILSPQFQY'
        'GSKVYVKHPADIPDYKKLSFPEGFKWERVMNFEDGGVVTVTQDSSLQDGCFIYKVKFIGVNFPSDGP'
        'VMQKKTMGWEASTERLYPRDGVLKGEIHKALKLKDGGHYLVEFKSIYMAKKPVQLPGYYYVDSKLDI'
        'TSHNEDYTIVEQYERTEGRHHLFL'),
    'TDH3': (
        'P00359',
        'native S. cerevisiae GAPDH - with yeast selected, compare the '
        'design against the real gene',
        'MVRVAINGFGRIGRLVMRIALSRPNVEVVALNDPFITNDYAAYMFKYDSTHGRYAGEVSHDDKHIIV'
        'DGKKIATYQERDPANLPWGSSNVDIAIDSTGVFKELDTAQKHIDAGAKKVVITAPSSTAPMFVMGVN'
        'EEKYTSDLKIVSNASCTTNCLAPLAKVINDAFGIEEGLMTTVHSLTATQKTVDGPSHKDWRGGRTAS'
        'GNIIPSSTGAAKAVGKVLPELQGKLTGMAFRVPTVDVSVVDLTVKLNKETTYDEIKKVVKAAAEGKL'
        'KGVLGYTEDAVVSSDFLGDSHSSIFDASAGIQLSPKFVKLVSWYDNEYGYSTRVVDLVEHVAKA'),
    'Insulin': (
        'P01308',
        'human preproinsulin, a classic yeast-expressed therapeutic',
        'MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVE'
        'LGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN'),
}

UNIPROT_URL = 'https://rest.uniprot.org/uniprotkb/{}.fasta'


@st.cache_data(show_spinner=False, ttl=3600)
def fetch_uniprot(accession: str) -> tuple:
    """Look up one UniProt accession and return (description, sequence).

    Cached for an hour so repeated lookups of the same id do not hammer the
    API, and wrapped so a network failure surfaces as a readable message
    rather than a stack trace in the browser.
    """
    acc = accession.strip().upper()
    if not re.fullmatch(r'[A-Z0-9]{6,10}', acc):
        raise ValueError('That does not look like a UniProt accession. '
                         'Expected something like P42212.')
    req = urllib.request.Request(
        UNIPROT_URL.format(acc),
        headers={'User-Agent': 'codon-optimiser-teaching-app'})
    try:
        with urllib.request.urlopen(req, timeout=15) as resp:
            text = resp.read().decode('utf-8')
    except urllib.error.HTTPError as exc:
        # UniProt answers 400 (not 404) for a well-formed accession that does
        # not exist, so both codes mean the same thing to the user.
        if exc.code in (400, 404):
            raise ValueError(f'No UniProt entry found for {acc}.') from exc
        raise ValueError(f'UniProt returned HTTP {exc.code}. '
                         'The service may be temporarily unavailable.') from exc
    except urllib.error.URLError as exc:
        raise ValueError(f'Could not reach UniProt: {exc.reason}') from exc

    lines = text.splitlines()
    if not lines or not lines[0].startswith('>'):
        raise ValueError(f'Unexpected response for {acc}.')
    head = lines[0]
    desc = head.split('|')[-1].split(' OS=')[0] if '|' in head else acc
    return desc, ''.join(lines[1:])

st.set_page_config(page_title='Codon Optimiser', page_icon='DNA', layout='wide')
apply_theme()


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
def wrap_dna(seq: str, width: int = 60) -> str:
    return '\n'.join(f'{i + 1:>5}  {seq[i:i + width]}'
                     for i in range(0, len(seq), width))


def to_fasta(seq: str, header: str, width: int = 60) -> str:
    body = '\n'.join(seq[i:i + width] for i in range(0, len(seq), width))
    return f'>{header}\n{body}\n'


def badge(ok: bool, good: str, bad: str) -> None:
    (st.success if ok else st.warning)(good if ok else bad, icon=None)


# --------------------------------------------------------------------------
# Sidebar - the design parameters
# --------------------------------------------------------------------------
with st.sidebar:
    draw(helix(width=300, height=44, period=38, opacity=0.9))

st.sidebar.header('Design parameters')

if len(ORGANISM_CONFIGS) == 1:
    # A one-item dropdown is pointless; show the target instead. The moment a
    # second organism is added to ORGANISM_CONFIGS this becomes a selector
    # again, with no other change needed.
    organism = next(iter(ORGANISM_CONFIGS))
    st.sidebar.markdown(f"**Target organism**  \n{ORGANISM_CONFIGS[organism]['label']}")
else:
    organism = st.sidebar.selectbox(
        'Target organism',
        list(ORGANISM_CONFIGS),
        format_func=lambda k: ORGANISM_CONFIGS[k]['label'],
        help='Codon usage table, restriction panel, motif set and GC window '
             'all switch with this choice.',
    )

expression = st.sidebar.radio(
    'Desired expression level',
    ['high', 'medium', 'low'],
    help='high = always the most-used codon (maximum CAI). '
         'medium = sample from real usage frequencies. '
         'low = deliberately de-optimised, for tuning expression down.',
)

cfg = ORGANISM_CONFIGS[organism]
lo, hi = cfg['gc_range']
st.sidebar.caption(
    f"**Target GC window:** {lo:.0%}-{hi:.0%}  \n"
    f"**Restriction panel:** {', '.join(cfg['restriction_sites'])}  \n"
    f"**Motifs avoided:** {len(cfg['avoid'])} patterns"
)

with st.sidebar.expander('Advanced'):
    extra_raw = st.text_input(
        'Additional sites to avoid',
        placeholder='GGTCTC, CGTCTC',
        help='Comma-separated DNA sequences, e.g. BsaI/BsmBI for Golden Gate. '
             'Both strands are blocked automatically.',
    )
    upstream = st.text_input(
        "5' UTR before the ATG",
        placeholder='AAAAAA',
        help="The last bases of your vector's 5' UTR, immediately before the "
             'start codon. Used only to score the Kozak initiation context - '
             'this tool does not design the UTR. Leave blank if unknown.',
    )
    add_stop = st.checkbox('Append stop codon', value=True)
    five_prime_check = st.checkbox(
        "5' structure sanity check", value=True,
        help="If the CAI-optimal design has a stable 5' hairpin, re-roll the "
             "first 20 codons and keep the weakest-folding variant that costs "
             "<= 0.05 CAI. Uncheck for a pure CAI baseline.",
    )
    seed = st.number_input('Random seed', value=0, step=1,
                           help='Same seed = same sequence, every time.')

extra_sites = [s.strip().upper() for s in extra_raw.split(',') if s.strip()]


# --------------------------------------------------------------------------
# Main - input
# --------------------------------------------------------------------------
st.title('Codon Optimiser')
st.caption('Protein sequence in, expression-optimised DNA out.')
draw(helix(width=1200, height=70, period=110, opacity=0.8))

if 'protein_text' not in st.session_state:
    st.session_state.protein_text = ''

cols = st.columns(len(EXAMPLES) + 1)
for col, (label, (acc, note, seq)) in zip(cols, EXAMPLES.items()):
    if col.button(f'{label}', use_container_width=True,
                  help=f'UniProt {acc} - {note} ({len(seq)} aa)'):
        st.session_state.protein_text = seq
if cols[-1].button('Clear', use_container_width=True):
    st.session_state.protein_text = ''

# This block must stay ABOVE the text_area: Streamlit forbids writing to a
# widget's session_state key after that widget has been created.
with st.expander('Fetch any protein from UniProt'):
    fc = st.columns([3, 1])
    accession = fc[0].text_input('Accession', placeholder='P42212',
                                 label_visibility='collapsed')
    if fc[1].button('Fetch', use_container_width=True):
        try:
            desc, seq = fetch_uniprot(accession)
            st.session_state.protein_text = seq
            st.success(f'{desc} - {len(seq)} aa')
        except ValueError as exc:
            st.error(str(exc))
    st.caption('Engineered variants (EGFP, mCherry, mKate) have no UniProt '
               'entry - they are not natural gene products. Use the parent '
               'protein (DsRed = Q9U6Y8) or paste the sequence directly.')

protein_text = st.text_area(
    'Protein sequence',
    key='protein_text',
    height=160,
    placeholder='Paste a protein sequence (one-letter code). '
                'FASTA headers, line breaks and spaces are handled.',
)

run = st.button('Optimise', type='primary')

# --------------------------------------------------------------------------
# Main - run and report
# --------------------------------------------------------------------------
if run:
    st.session_state.pop('result', None)

    if not protein_text.strip():
        st.info('Paste a protein sequence, or load one of the examples above.')
        st.stop()

    try:
        cleaned = clean_protein(protein_text)
    except ValueError as exc:
        st.error(str(exc))
        st.stop()

    if len(cleaned) > MAX_PROTEIN_LEN:
        st.error(f'Sequence is {len(cleaned)} aa; this app is capped at '
                 f'{MAX_PROTEIN_LEN} aa.')
        st.stop()

    with st.spinner('Optimising...'):
        st.session_state.result = optimize(
            protein_text, organism, expression,
            extra_sites=extra_sites, add_stop=add_stop,
            seed=int(seed), five_prime_check=five_prime_check,
            upstream=upstream)

# Render from session state, not straight from the button. A button is True
# for one rerun only, so rendering inside `if run:` makes the whole report
# vanish the moment anything else on the page reruns the script.
r = st.session_state.get('result')

if r:
    # Use the config the RESULT was built with, not whatever the sidebar
    # currently shows - the two differ as soon as the user changes organism
    # without re-running.
    rcfg = ORGANISM_CONFIGS[r['organism']]
    lo, hi = r['gc_range']

    draw(squiggle_rule())
    st.subheader(f"{r['organism_label']} - {r['expression']} expression")

    m = st.columns(4)
    m[0].metric('CAI', f"{r['cai']:.3f}")
    m[1].metric('GC content', f"{r['gc']:.1%}",
                delta='in window' if r['gc_ok'] else 'out of window',
                delta_color='normal' if r['gc_ok'] else 'inverse')
    m[2].metric('Length', f"{r['dna_length']} nt",
                delta=f"{r['protein_length']} aa")
    m[3].metric('Back-translation',
                'PASS' if r['translates_back'] else 'FAIL')

    if not r['translates_back']:
        st.error('The designed DNA does not translate back to the input '
                 'protein. Do not use this sequence; please report it.')

    tab_seq, tab_qc, tab_gc, tab_props = st.tabs(
        ['Sequence', 'Constraints', 'GC profile', 'Other properties'])

    # ---- sequence -------------------------------------------------------
    with tab_seq:
        st.code(wrap_dna(r['dna']), language=None)
        header = f"{organism}_{expression}_CAI{r['cai']}_GC{r['gc']}"
        c1, c2 = st.columns(2)
        c1.download_button('Download FASTA',
                           to_fasta(r['dna'], header),
                           file_name=f'{organism}_{expression}.fasta',
                           mime='text/plain', use_container_width=True)
        c2.download_button('Download plain sequence', r['dna'],
                           file_name=f'{organism}_{expression}.txt',
                           mime='text/plain', use_container_width=True)

    # ---- constraints ----------------------------------------------------
    with tab_qc:
        sites = r['restriction_sites_found']
        badge(not sites,
              f"No restriction sites from the {rcfg['label']} panel.",
              f'Restriction sites still present: {sites}')

        left = r['motifs_avoided_but_present']
        badge(not left,
              'All avoided motifs eliminated.',
              f'Motifs that could not be removed: '
              f'{ {k: len(v) for k, v in left.items()} }')

        hp = r['homopolymers']
        long_runs = [h for h in hp if h['length'] >= 8]
        if long_runs:
            st.warning(
                'Long single-base runs (DNA synthesis risk): '
                + ', '.join(f"{h['base']}x{h['length']} at nt {h['start'] + 1}"
                            for h in long_runs))
        elif hp:
            st.info('Short runs (>=6 nt), usually tolerated: '
                    + ', '.join(f"{h['base']}x{h['length']} at nt {h['start'] + 1}"
                                for h in hp))
        else:
            st.success('No homopolymer runs of 6 nt or more.')

        rep = r['motifs_reported']
        if rep:
            st.caption('Reported but not designed against (too common to '
                       'eliminate without distorting the sequence):')
            for k, v in rep.items():
                st.caption(f'  - {k}: {len(v)} occurrence(s)')

    # ---- GC profile -----------------------------------------------------
    with tab_gc:
        prof = r['gc_profile']
        if prof:
            try:
                import pandas as pd
                from theme import PALETTE
                df = pd.DataFrame(prof, columns=['position', 'GC'])
                df['target min'] = lo
                df['target max'] = hi
                st.line_chart(
                    df.set_index('position'),
                    color=[PALETTE['primary'], PALETTE['soft'], PALETTE['soft']],
                )
            except ImportError:
                st.line_chart([p[1] for p in prof])
            st.caption(f'30 nt sliding window, 3 nt step. Target window '
                       f'{lo:.0%}-{hi:.0%}. A mean inside the window can still '
                       f'hide local GC-rich or GC-poor stretches, which is '
                       f'what this plot is for.')
        else:
            st.caption('Sequence too short to profile.')

    # ---- other properties ----------------------------------------------
    with tab_props:
        for name, c in r['extra_checks'].items():
            label = name.replace('_', ' ')
            detail = {k: v for k, v in c.items() if k not in ('flag', 'note')}
            if c.get('flag'):
                st.warning(f"**{label}** - {c['note']}\n\n`{detail}`")
            elif c.get('status') in ('unknown', 'not_applicable'):
                # Neither "not assessed" nor "does not apply" means "fine" -
                # never show either as a green pass.
                word = ('not assessed' if c['status'] == 'unknown'
                        else 'not applicable to this host')
                st.info(f"**{label}** - {word}\n\n{c['note']}")
            else:
                st.success(f"**{label}** - OK\n\n`{detail}`")
            if name == 'five_prime_structure':
                st.caption('Heuristic nearest-neighbour stem score, not a '
                           'thermodynamic folding prediction. Treat it as a '
                           'relative indicator.')

        fp = r['five_prime_optimization']
        st.divider()
        st.markdown(f"**5' tuning** - {fp.get('reason', 'n/a')}")
        if fp.get('applied'):
            st.caption(f"dG {fp['baseline_dG']} -> {fp['final_dG']} kcal/mol, "
                       f"CAI {fp['baseline_cai']} -> {fp['final_cai']} "
                       f"(cost {fp['cai_cost']})")

st.divider()
st.caption(
    'Codon usage: Kazusa Codon Usage Database (Nakamura, Gojobori & Ikemura, '
    '1999 Nucleic Acids Res 28:292). CAI after Sharp & Li (1987). '
    'This is a teaching tool - verify any construct before synthesis.'
)
