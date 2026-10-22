"""
codon_optimizer.py - reverse-translation with codon-usage optimisation.

Primary target organism : Saccharomyces cerevisiae
Comparison organism     : Escherichia coli K-12

Pure standard library (math / random / re) so it runs anywhere, including
a Streamlit Cloud container, with no compiled dependencies.

Codon frequencies are within-family fractions (each amino-acid family sums
to 1.0) transcribed from the Kazusa Codon Usage Database:
  Nakamura Y., Gojobori T., Ikemura T. (2000) Nucleic Acids Res. 28:292.
See LICENSING NOTES at the bottom of this file before publishing.
"""

from __future__ import annotations

import math
import random
import re

# ==========================================================================
# 1. Genetic code
# ==========================================================================
CODON_TABLE = {
    'A': ['GCT', 'GCC', 'GCA', 'GCG'],
    'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'N': ['AAT', 'AAC'], 'D': ['GAT', 'GAC'], 'C': ['TGT', 'TGC'],
    'Q': ['CAA', 'CAG'], 'E': ['GAA', 'GAG'],
    'G': ['GGT', 'GGC', 'GGA', 'GGG'], 'H': ['CAT', 'CAC'],
    'I': ['ATT', 'ATC', 'ATA'],
    'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
    'K': ['AAA', 'AAG'], 'M': ['ATG'], 'F': ['TTT', 'TTC'],
    'P': ['CCT', 'CCC', 'CCA', 'CCG'],
    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
    'T': ['ACT', 'ACC', 'ACA', 'ACG'], 'W': ['TGG'], 'Y': ['TAT', 'TAC'],
    'V': ['GTT', 'GTC', 'GTA', 'GTG'], '*': ['TAA', 'TAG', 'TGA'],
}

AA_BY_CODON = {c: aa for aa, cods in CODON_TABLE.items() for c in cods}

# Families with a single codon carry no usage information and are excluded
# from CAI by the original Sharp & Li (1987) definition, together with stops.
CAI_EXCLUDED_AA = {'M', 'W', '*'}

VALID_AA = set(CODON_TABLE)

_COMPLEMENT = str.maketrans('ACGT', 'TGCA')


def revcomp(dna: str) -> str:
    """Reverse complement of an uppercase A/C/G/T string."""
    return dna.translate(_COMPLEMENT)[::-1]


# ==========================================================================
# 2. Codon usage frequencies (fraction within each synonymous family)
# ==========================================================================
FREQS = {
    'S_cerevisiae': {
        'A': {'GCT': 0.38, 'GCC': 0.22, 'GCA': 0.29, 'GCG': 0.11},
        'R': {'CGT': 0.14, 'CGC': 0.06, 'CGA': 0.07, 'CGG': 0.04,
              'AGA': 0.48, 'AGG': 0.21},
        'N': {'AAT': 0.59, 'AAC': 0.41}, 'D': {'GAT': 0.65, 'GAC': 0.35},
        'C': {'TGT': 0.63, 'TGC': 0.37}, 'Q': {'CAA': 0.69, 'CAG': 0.31},
        'E': {'GAA': 0.70, 'GAG': 0.30},
        'G': {'GGT': 0.47, 'GGC': 0.19, 'GGA': 0.22, 'GGG': 0.12},
        'H': {'CAT': 0.64, 'CAC': 0.36},
        'I': {'ATT': 0.46, 'ATC': 0.26, 'ATA': 0.27},
        'L': {'TTA': 0.28, 'TTG': 0.29, 'CTT': 0.13, 'CTC': 0.06,
              'CTA': 0.14, 'CTG': 0.11},
        'K': {'AAA': 0.58, 'AAG': 0.42}, 'M': {'ATG': 1.00},
        'F': {'TTT': 0.59, 'TTC': 0.41},
        'P': {'CCT': 0.31, 'CCC': 0.15, 'CCA': 0.42, 'CCG': 0.12},
        'S': {'TCT': 0.26, 'TCC': 0.16, 'TCA': 0.21, 'TCG': 0.10,
              'AGT': 0.16, 'AGC': 0.11},
        'T': {'ACT': 0.35, 'ACC': 0.22, 'ACA': 0.30, 'ACG': 0.14},
        'W': {'TGG': 1.00}, 'Y': {'TAT': 0.56, 'TAC': 0.44},
        'V': {'GTT': 0.39, 'GTC': 0.21, 'GTA': 0.21, 'GTG': 0.19},
        '*': {'TAA': 0.47, 'TAG': 0.23, 'TGA': 0.30},
    },
    'E_coli_K12': {
        'A': {'GCT': 0.16, 'GCC': 0.27, 'GCA': 0.21, 'GCG': 0.36},
        'R': {'CGT': 0.38, 'CGC': 0.40, 'CGA': 0.06, 'CGG': 0.10,
              'AGA': 0.04, 'AGG': 0.02},
        'N': {'AAT': 0.45, 'AAC': 0.55}, 'D': {'GAT': 0.63, 'GAC': 0.37},
        'C': {'TGT': 0.44, 'TGC': 0.56}, 'Q': {'CAA': 0.35, 'CAG': 0.65},
        'E': {'GAA': 0.69, 'GAG': 0.31},
        'G': {'GGT': 0.34, 'GGC': 0.41, 'GGA': 0.11, 'GGG': 0.15},
        'H': {'CAT': 0.57, 'CAC': 0.43},
        'I': {'ATT': 0.51, 'ATC': 0.42, 'ATA': 0.07},
        'L': {'TTA': 0.13, 'TTG': 0.13, 'CTT': 0.10, 'CTC': 0.10,
              'CTA': 0.04, 'CTG': 0.50},
        'K': {'AAA': 0.76, 'AAG': 0.24}, 'M': {'ATG': 1.00},
        'F': {'TTT': 0.57, 'TTC': 0.43},
        'P': {'CCT': 0.16, 'CCC': 0.12, 'CCA': 0.19, 'CCG': 0.53},
        'S': {'TCT': 0.15, 'TCC': 0.15, 'TCA': 0.12, 'TCG': 0.15,
              'AGT': 0.15, 'AGC': 0.28},
        'T': {'ACT': 0.16, 'ACC': 0.44, 'ACA': 0.13, 'ACG': 0.27},
        'W': {'TGG': 1.00}, 'Y': {'TAT': 0.57, 'TAC': 0.43},
        'V': {'GTT': 0.26, 'GTC': 0.22, 'GTA': 0.15, 'GTG': 0.37},
        '*': {'TAA': 0.64, 'TAG': 0.07, 'TGA': 0.29},
    },
}


def _validate_freqs() -> None:
    """Fail loudly if a table is incomplete or badly rounded.

    Kazusa fractions are published to two decimals, so a family can sum to
    0.99 or 1.01. Anything further off means a transcription error.
    """
    for organism, table in FREQS.items():
        if set(table) != VALID_AA:
            missing = VALID_AA ^ set(table)
            raise ValueError(f'{organism}: amino-acid families off: {missing}')
        for aa, fam in table.items():
            if set(fam) != set(CODON_TABLE[aa]):
                raise ValueError(f'{organism}/{aa}: codon set does not match '
                                 'the genetic code')
            total = sum(fam.values())
            if not 0.95 <= total <= 1.05:
                raise ValueError(f'{organism}/{aa}: fractions sum to {total}')


_validate_freqs()


def build_weights(freqs: dict) -> dict:
    """Relative adaptiveness w = f / f_max within each family (best codon = 1.0).

    Fractions are renormalised per family first, so two-decimal rounding in
    the published tables cannot bias one amino acid against another.
    """
    weights = {}
    for aa, fam in freqs.items():
        total = sum(fam.values())
        norm = {c: f / total for c, f in fam.items()}
        fmax = max(norm.values())
        for c, f in norm.items():
            weights[c] = f / fmax
    return weights


# ==========================================================================
# 3. Per-organism design rules
# ==========================================================================
# Motifs are split into two groups:
#   'avoid'  -> actively designed against during back-translation
#   'report' -> only scanned afterwards (too common to eliminate entirely)
ORGANISM_CONFIGS = {
    'S_cerevisiae': {
        'label': 'Saccharomyces cerevisiae (baker\'s yeast)',
        'restriction_sites': {
            'EcoRI': 'GAATTC', 'BamHI': 'GGATCC', 'HindIII': 'AAGCTT',
            'XhoI': 'CTCGAG', 'SalI': 'GTCGAC', 'NotI': 'GCGGCCGC',
            'SpeI': 'ACTAGT', 'XbaI': 'TCTAGA',
        },
        'avoid': {
            'premature polyA signal': r'AATAAA',
            'AU-rich decay element': r'ATTTA',
            'cryptic 5\' splice donor': r'GTATGT',
            'rare-codon cluster': r'(?:CGG|CGC|CGA|CCG|TCG|ACG|GCG){3,}',
        },
        'report': {
            'AT-rich cryptic terminator': r'[AT]{6,}',
        },
        'gc_range': (0.35, 0.45),
    },
    'E_coli_K12': {
        'label': 'Escherichia coli K-12',
        'restriction_sites': {
            'EcoRI': 'GAATTC', 'BamHI': 'GGATCC', 'HindIII': 'AAGCTT',
            'XhoI': 'CTCGAG', 'SalI': 'GTCGAC', 'NotI': 'GCGGCCGC',
            'SpeI': 'ACTAGT', 'XbaI': 'TCTAGA',
            'NdeI': 'CATATG',   # contains ATG - classic re-cloning trap
            'NcoI': 'CCATGG',   # contains ATG - classic re-cloning trap
        },
        'avoid': {
            'internal Shine-Dalgarno': r'AGG[AG]GG|AGGAG|GAGGG',
            'terminator-like poly-T': r'T{6,}',
            'rare-codon cluster': r'(?:AGA|AGG|ATA|CCC|CTA|GGA){3,}',
            'T7-terminator-like': r'TAAAT[AT]{5,}',
        },
        'report': {
            'dam methylation site (GATC)': r'GATC',
            'dcm methylation site (CCWGG)': r'CC[AT]GG',
        },
        'gc_range': (0.40, 0.60),
    },
}


# ==========================================================================
# 4. Input handling
# ==========================================================================
def clean_protein(protein: str, allow_stop: bool = True) -> str:
    """Normalise user input into a bare uppercase one-letter sequence.

    Strips FASTA headers, whitespace, digits and a single trailing '*'.
    Raises ValueError on any character that is not a standard amino acid,
    which is what you want in a web app: a clear message, not a KeyError.
    """
    lines = [ln for ln in protein.splitlines() if not ln.startswith('>')]
    seq = re.sub(r'[\s\d]', '', ''.join(lines)).upper()
    if not seq:
        raise ValueError('Empty protein sequence.')
    if seq.endswith('*'):
        if not allow_stop:
            raise ValueError('Stop symbol not allowed here.')
        seq = seq[:-1]
    bad = sorted(set(seq) - VALID_AA | ({'*'} & set(seq)))
    if bad:
        raise ValueError(
            'Unsupported symbol(s): ' + ', '.join(bad) +
            '. Use the 20 standard amino acids in one-letter code '
            '(ambiguity codes B/Z/J/X and selenocysteine U are not supported).'
        )
    return seq


# ==========================================================================
# 5. Forbidden-pattern machinery
# ==========================================================================
class Forbidden:
    """A set of patterns the designed sequence should not contain.

    Holds compiled regexes (literal restriction sites are regex-escaped) and
    answers one question efficiently during back-translation: does appending
    this codon CREATE a new hit? Only the junction is re-examined, so a hit
    that already exists upstream never blocks further codon choices.
    """

    #  longest pattern we expect to match, in nucleotides
    PAD = 18

    def __init__(self, patterns: dict | None = None):
        self.patterns = dict(patterns or {})
        self._rx = [re.compile(p) for p in self.patterns.values()]

    def __bool__(self) -> bool:
        return bool(self._rx)

    def introduces(self, prefix: str, codon: str) -> bool:
        """True if `prefix + codon` ends in a hit that `prefix` did not have."""
        if not self._rx:
            return False
        tail = prefix[-self.PAD:]
        boundary = len(tail)
        window = tail + codon
        for rx in self._rx:
            for m in rx.finditer(window):
                if m.end() > boundary:     # the match uses the new codon
                    return True
        return False

    def count(self, dna: str) -> int:
        return sum(len(rx.findall(dna)) for rx in self._rx)

    def scan(self, dna: str) -> dict:
        """{label: [(start, end), ...]} for every pattern that hits."""
        out = {}
        for label, pat in self.patterns.items():
            spans = [m.span() for m in re.finditer(pat, dna)]
            if spans:
                out[label] = spans
        return out


def build_forbidden(cfg: dict, extra_sites=()) -> Forbidden:
    """Assemble the avoid-set: restriction sites (both strands) + hard motifs."""
    patterns = {}
    for name, site in cfg['restriction_sites'].items():
        patterns[f'{name} ({site})'] = re.escape(site)
        rc = revcomp(site)
        if rc != site:                      # non-palindromic: block both strands
            patterns[f'{name} rev ({rc})'] = re.escape(rc)
    for site in extra_sites:
        site = site.strip().upper()
        if site:
            patterns[f'user site ({site})'] = re.escape(site)
            rc = revcomp(site)
            if rc != site:
                patterns[f'user site rev ({rc})'] = re.escape(rc)
    patterns.update(cfg['avoid'])
    return Forbidden(patterns)


# ==========================================================================
# 6. Back-translation
# ==========================================================================
def reverse_translate(protein: str, weights: dict, forbidden: Forbidden | None = None,
                      expression: str = 'high', seed: int | None = None,
                      head_sample_codons: int = 0) -> str:
    """Back-translate a protein, codon by codon, avoiding forbidden patterns.

    expression:
      'high'   - always the highest-weight codon that stays clean (max CAI)
      'medium' - sample proportionally to real usage frequency
      'low'    - sample proportionally to rarity (deliberate de-optimisation)

    head_sample_codons:
      Number of codons at the 5' end to pick by weighted sampling instead of
      argmax, even in 'high' mode. This is the only knob the 5' structure
      search turns: it generates alternative starts for the SAME protein
      while the rest of the gene stays strictly CAI-maximal.
    """
    if expression not in ('high', 'medium', 'low'):
        raise ValueError("expression must be 'high', 'medium' or 'low'")
    forbidden = forbidden or Forbidden()
    rng = random.Random(seed)
    out = []
    dna = ''                                # running string: no O(n^2) rejoin

    for idx, aa in enumerate(protein):
        cods = CODON_TABLE[aa]

        if expression == 'high' and idx >= head_sample_codons:
            ranked = sorted(cods, key=lambda c: weights.get(c, 0.0), reverse=True)
        else:
            w = {c: max(weights.get(c, 0.0), 1e-4) for c in cods}
            if expression == 'low':
                w = {c: max((1.0 - v) ** 2, 1e-4) for c, v in w.items()}
            # Draw a random ORDER once, then walk it. Unlike resampling in a
            # loop this is guaranteed to try every synonym before giving up.
            ranked = []
            pool = dict(w)
            while pool:
                pick = rng.choices(list(pool), weights=list(pool.values()))[0]
                ranked.append(pick)
                del pool[pick]

        chosen = next((c for c in ranked if not forbidden.introduces(dna, c)), None)
        if chosen is None:
            # Nothing is clean here (can happen inside a Trp/Met stretch).
            # Fall back to the first-ranked codon and let the QC report say so.
            chosen = ranked[0]
        out.append(chosen)
        dna += chosen

    return dna


# ==========================================================================
# 7. Metrics
# ==========================================================================
def codons_of(dna: str) -> list:
    return [dna[i:i + 3] for i in range(0, len(dna) - len(dna) % 3, 3)]


def cai(dna: str, weights: dict) -> float:
    """Codon Adaptation Index: geometric mean of relative adaptiveness.

    Met, Trp and stop codons are excluded (Sharp & Li 1987): they have no
    synonymous choice, so including them would inflate the score toward 1.
    """
    ws = [weights[c] for c in codons_of(dna)
          if c in weights and AA_BY_CODON.get(c) not in CAI_EXCLUDED_AA]
    ws = [w for w in ws if w > 0]
    if not ws:
        return 0.0
    return math.exp(sum(math.log(w) for w in ws) / len(ws))


def gc_content(dna: str) -> float:
    return (dna.count('G') + dna.count('C')) / len(dna) if dna else 0.0


def gc_profile(dna: str, window: int = 30, step: int = 3) -> list:
    """Sliding-window GC, for plotting and for spotting local GC spikes."""
    if len(dna) < window:
        return [(0, gc_content(dna))] if dna else []
    return [(i, gc_content(dna[i:i + window]))
            for i in range(0, len(dna) - window + 1, step)]


def find_restriction_sites(dna: str, sites: dict) -> dict:
    """{'EcoRI (GAATTC)': [positions]} - searches both strands."""
    hits = {}
    for name, site in sites.items():
        for tag, seq in (('', site), (' [rev]', revcomp(site))):
            if tag and seq == site:
                continue
            pos = [m.start() for m in re.finditer(f'(?={re.escape(seq)})', dna)]
            if pos:
                hits[f'{name}{tag} ({seq})'] = pos
    return hits


def find_motifs(dna: str, motifs: dict) -> dict:
    out = {}
    for label, pat in motifs.items():
        spans = [m.span() for m in re.finditer(pat, dna)]
        if spans:
            out[label] = spans
    return out


def homopolymer_runs(dna: str, min_length: int = 6) -> list:
    """Single-nucleotide runs >= min_length (DNA-synthesis slippage risk)."""
    rx = r'(A|C|G|T)\1{%d,}' % (min_length - 1)
    return [{'base': m.group(1), 'length': len(m.group(0)), 'start': m.start()}
            for m in re.finditer(rx, dna)]


# ==========================================================================
# 8. GC repair
# ==========================================================================
def repair_gc(dna: str, weights: dict, gc_range: tuple,
              forbidden: Forbidden | None = None, max_passes: int = 8,
              prefer: str = 'common') -> str:
    """Walk GC content into the target window using synonymous swaps.

    This resolves the CAI-vs-GC conflict: a max-CAI yeast design comes out
    AT-rich because yeast prefers NNT/NNA codons, below the 35-45% window.

    Three rules keep it honest:
      * GC is re-evaluated after EVERY single swap, so the loop stops the
        moment it enters the window instead of overshooting and oscillating.
      * A swap is rejected if it increases the number of forbidden hits
        ANYWHERE. Comparing counts (rather than asking "is the sequence
        clean?") means a pre-existing hit cannot freeze the whole repair.
      * `prefer` keeps the repair aligned with the expression setting. With
        prefer='rare' the tie-break favours low-weight synonyms, so fixing
        GC on a 'low' design does not silently re-optimise it back up.
    """
    forbidden = forbidden or Forbidden()
    codons = codons_of(dna)
    if not codons:
        return dna

    def gc_of(cs):
        s = ''.join(cs)
        return (s.count('G') + s.count('C')) / len(s)

    lo, hi = gc_range
    base_bad = forbidden.count(''.join(codons))
    tie_rng = random.Random(0)              # deterministic tie-spreading
    EPS = 0.05                              # keeps zero-cost swaps finite

    for _ in range(max_passes):
        if lo <= gc_of(codons) <= hi:
            break
        want_more_gc = gc_of(codons) < lo

        # Score EVERY legal swap in the whole gene, then spend the cheapest
        # ones first. Ranking globally (rather than scanning left to right)
        # is what keeps the repair from dumping all its CAI damage into the
        # 5' end, which is both the costliest place to lose CAI and the one
        # region where raising GC hurts translation initiation.
        cands = []
        for i, c in enumerate(codons):
            cur_gc = c.count('G') + c.count('C')
            log_cur = math.log(max(weights.get(c, 0.0), 1e-4))
            for alt in CODON_TABLE[AA_BY_CODON[c]]:
                if alt == c:
                    continue
                alt_gc = alt.count('G') + alt.count('C')
                if alt_gc == cur_gc or (alt_gc > cur_gc) != want_more_gc:
                    continue
                # CAI is a geometric mean, so log-weight differences are the
                # quantity that actually adds up across the sequence.
                dlog = math.log(max(weights.get(alt, 0.0), 1e-4)) - log_cur
                penalty = max(0.0, -dlog) if prefer == 'common' else max(0.0, dlog)
                score = abs(alt_gc - cur_gc) / (penalty + EPS)
                cands.append((score, i, alt))

        if not cands:
            break

        tie_rng.shuffle(cands)                      # spread equal-cost swaps
        cands.sort(key=lambda t: -t[0])             # stable: ties stay spread

        used = set()
        applied = False
        for _score, i, alt in cands:
            if lo <= gc_of(codons) <= hi:
                break
            if i in used:                   # one swap per position per pass
                continue
            trial = ''.join(codons[:i]) + alt + ''.join(codons[i + 1:])
            if forbidden.count(trial) > base_bad:
                continue
            codons[i] = alt
            used.add(i)
            applied = True

        if not applied:
            break                           # no legal move left; stop early

    return ''.join(codons)


# ==========================================================================
# 9. Extra sequence-property checks
# ==========================================================================
# Nearest-neighbour stacking free energies, kcal/mol, 1 M NaCl 37 C
# (SantaLucia 1998 unified parameters). Used only as a coarse 5'-structure
# proxy - see LICENSING NOTES for the real-folding option.
_NN_DG = {
    'AA': -1.00, 'AT': -0.88, 'TA': -0.58, 'CA': -1.45, 'GT': -1.44,
    'CT': -1.28, 'GA': -1.30, 'CG': -2.17, 'GC': -2.24, 'GG': -1.84,
}


def _stack_dg(pair: str) -> float:
    return _NN_DG.get(pair, _NN_DG.get(revcomp(pair), -1.0))


def five_prime_structure_score(dna: str, window: int = 48) -> dict:
    """Coarse estimate of secondary structure in the ribosome-landing region.

    Scores the most stable local stem the 5' window can form with itself:
    stronger (more negative) means a hairpin more likely to block initiation.
    This is a heuristic, NOT a thermodynamic folding prediction.
    """
    head = dna[:window]
    best_dg, best_at = 0.0, None
    n = len(head)
    for i in range(n):
        for j in range(i + 7, n):           # >= 3 nt loop between the arms
            stem = 0
            dg = 0.0
            while (i - stem >= 0 and j + stem < n
                   and head[i - stem] == revcomp(head[j + stem])):
                if stem:
                    dg += _stack_dg(head[i - stem:i - stem + 2])
                stem += 1
            if stem >= 4 and dg < best_dg:
                best_dg, best_at = dg, (i - stem + 1, j + stem)
    return {
        'window_nt': min(window, len(dna)),
        'estimated_dG': round(best_dg, 2),
        'stem_span': best_at,
        'flag': best_dg < -8.0,
        'note': 'heuristic stem score, not a folding prediction',
    }


# Destabilising N-terminal residues. Yeast: Bachmair/Varshavsky N-end rule.
# E. coli: the ClpS/Leu-Phe-transferase branch.
N_END_DESTABILISING = {
    'S_cerevisiae': set('RKHFLWIYND'),
    'E_coli_K12': set('RKLFWY'),
}


def extra_checks(dna: str, protein: str, organism: str) -> dict:
    """Properties beyond codon usage that affect transcription/translation."""
    cfg = ORGANISM_CONFIGS[organism]
    checks = {}

    checks['five_prime_structure'] = five_prime_structure_score(dna)

    # 5'-terminal GC: a GC-poor ramp eases ribosome loading
    head = dna[:48]
    checks['five_prime_gc'] = {
        'value': round(gc_content(head), 3),
        'flag': gc_content(head) > cfg['gc_range'][1] + 0.05,
        'note': 'high 5\' GC can occlude the ribosome-binding region',
    }

    # Second codon: position 2 has a measurable effect on initiation.
    # The note is organism-specific so a report never mentions the host the
    # user did not select.
    second = dna[3:6]
    if organism == 'E_coli_K12':
        second_flag = second.startswith('A')
        second_note = ('A-starting second codons can promote internal '
                       're-initiation in E. coli')
    else:
        second_flag = False
        second_note = ('position 2 influences initiation efficiency; '
                       'no strong consensus rule in S. cerevisiae')
    checks['second_codon'] = {
        'codon': second,
        'aa': AA_BY_CODON.get(second),
        'flag': second_flag,
        'note': second_note,
    }

    # N-end rule: the residue exposed after Met excision sets protein half-life
    exposed = protein[1] if len(protein) > 1 and protein[0] == 'M' else protein[0]
    checks['n_end_rule'] = {
        'n_terminal_residue': protein[0],
        'residue_after_Met_excision': exposed,
        'flag': exposed in N_END_DESTABILISING[organism],
        'note': 'destabilising N-terminal residue - short half-life expected',
    }

    if organism == 'E_coli_K12':
        dam = [m.start() for m in re.finditer(r'(?=GATC)', dna)]
        dcm = [m.start() for m in re.finditer(r'(?=CC[AT]GG)', dna)]
        checks['dam_dcm_sites'] = {
            'dam_GATC': dam, 'dcm_CCWGG': dcm,
            'flag': bool(dam or dcm),
            'note': 'methylated in dam+/dcm+ strains; blocks some REs on cloning',
        }

    return checks


# ==========================================================================
# 10. 5' structure sanity check
# ==========================================================================
def five_prime_sanity_check(protein: str, weights: dict, cfg: dict,
                            forbidden: Forbidden, expression: str,
                            baseline_dna: str,
                            head_codons: int = 20, n_variants: int = 40,
                            max_cai_drop: float = 0.05,
                            dg_threshold: float = -8.0, seed: int = 0) -> dict:
    """Tie-break on 5' structure WITHOUT displacing CAI as the objective.

    CAI remains the thing being maximised. This step only runs if the
    CAI-maximal design has a 5' hairpin stable enough to impair ribosome
    loading, and it is only allowed to re-roll the first `head_codons`
    codons. Everything downstream of the head is still built strictly
    CAI-first by the same reverse_translate call.

    A variant is accepted only if its CAI stays within `max_cai_drop` of the
    baseline. If no variant clears that bar, the baseline is kept and the
    report says so - the check never silently trades away CAI.
    """
    base_struct = five_prime_structure_score(baseline_dna)
    base_cai = cai(baseline_dna, weights)

    report = {
        'applied': False,
        'head_codons': head_codons,
        'variants_tested': 0,
        'baseline_dG': base_struct['estimated_dG'],
        'baseline_cai': round(base_cai, 3),
        'final_dG': base_struct['estimated_dG'],
        'final_cai': round(base_cai, 3),
        'cai_cost': 0.0,
        'dg_threshold': dg_threshold,
        'max_cai_drop': max_cai_drop,
    }

    # No 5' problem -> change nothing. The common case costs one scan.
    if base_struct['estimated_dG'] >= dg_threshold:
        report['reason'] = ('baseline 5\' structure already acceptable '
                            '(dG >= threshold); sequence unchanged')
        return {'dna': baseline_dna, 'report': report}

    prefer = 'rare' if expression == 'low' else 'common'
    best = None                             # (dG, -cai_drop, dna)

    for v in range(n_variants):
        dna = reverse_translate(protein, weights, forbidden, expression,
                                seed=seed + 1 + v, head_sample_codons=head_codons)
        dna = repair_gc(dna, weights, cfg['gc_range'], forbidden, prefer=prefer)

        v_cai = cai(dna, weights)
        if base_cai - v_cai > max_cai_drop:     # too expensive in CAI terms
            continue

        v_dg = five_prime_structure_score(dna)['estimated_dG']
        # rank on structure first, then prefer the higher-CAI variant
        cand = (-v_dg, v_cai, dna)
        if best is None or cand > best:
            best = cand

    report['variants_tested'] = n_variants

    if best is None:
        report['reason'] = (f'no variant stayed within {max_cai_drop} CAI of '
                            'the baseline; CAI kept, 5\' structure not improved')
        return {'dna': baseline_dna, 'report': report}

    best_dg, best_cai, best_dna = -best[0], best[1], best[2]

    if best_dg <= base_struct['estimated_dG']:
        report['reason'] = ('no affordable variant folded more weakly than '
                            'the baseline; sequence unchanged')
        return {'dna': baseline_dna, 'report': report}

    report.update({
        'applied': True,
        'final_dG': round(best_dg, 2),
        'final_cai': round(best_cai, 3),
        'cai_cost': round(base_cai - best_cai, 3),
        'reason': (f'5\' structure relaxed from {base_struct["estimated_dG"]} '
                   f'to {round(best_dg, 2)} kcal/mol for '
                   f'{round(base_cai - best_cai, 3)} CAI'),
    })
    return {'dna': best_dna, 'report': report}


# ==========================================================================
# 11. One-call pipeline
# ==========================================================================
def optimize(protein: str, organism: str, expression: str = 'high',
             extra_sites=(), add_stop: bool = True, seed: int = 0,
             five_prime_check: bool = True) -> dict:
    """Protein sequence in -> optimised DNA plus a full QC report.

    five_prime_check: if the CAI-maximal design has a stable 5' hairpin,
    re-roll only the first 20 codons and keep the weakest-folding variant
    that costs <= 0.05 CAI. Set False for a pure CAI baseline.
    """
    if organism not in ORGANISM_CONFIGS:
        raise ValueError(f'Unknown organism: {organism}. '
                         f'Choose from {sorted(ORGANISM_CONFIGS)}.')

    cfg = ORGANISM_CONFIGS[organism]
    weights = build_weights(FREQS[organism])
    aa = clean_protein(protein)
    coding_aa = aa + '*' if add_stop else aa

    forbidden = build_forbidden(cfg, extra_sites)
    dna = reverse_translate(coding_aa, weights, forbidden, expression, seed)
    dna = repair_gc(dna, weights, cfg['gc_range'], forbidden,
                    prefer='rare' if expression == 'low' else 'common')

    if five_prime_check:
        fp = five_prime_sanity_check(coding_aa, weights, cfg, forbidden,
                                     expression, dna, seed=seed)
        dna, fp_report = fp['dna'], fp['report']
    else:
        fp_report = {'applied': False, 'reason': 'disabled (pure CAI baseline)'}

    gc = gc_content(dna)
    back = ''.join(AA_BY_CODON[c] for c in codons_of(dna))

    return {
        'organism': organism,
        'organism_label': cfg['label'],
        'expression': expression,
        'protein': aa,
        'protein_length': len(aa),
        'dna': dna,
        'dna_length': len(dna),
        'cai': round(cai(dna, weights), 3),
        'gc': round(gc, 3),
        'gc_range': cfg['gc_range'],
        'gc_ok': cfg['gc_range'][0] <= gc <= cfg['gc_range'][1],
        'gc_profile': gc_profile(dna),
        'restriction_sites_found': find_restriction_sites(dna, cfg['restriction_sites']),
        'motifs_avoided_but_present': find_motifs(dna, cfg['avoid']),
        'motifs_reported': find_motifs(dna, cfg['report']),
        'homopolymers': homopolymer_runs(dna),
        'extra_checks': extra_checks(dna, coding_aa, organism),
        'five_prime_optimization': fp_report,
        # Hard correctness gate: the DNA must translate back to the input.
        'translates_back': back == coding_aa,
    }


def compare_organisms(protein: str, expression: str = 'high', **kw) -> dict:
    """Design for both hosts and cross-score each design in the other host.

    The cross-CAI is the interesting number: a yeast-optimised gene scored
    with E. coli weights shows directly how host-specific the design is.
    """
    results = {org: optimize(protein, org, expression, **kw)
               for org in ORGANISM_CONFIGS}
    for org, res in results.items():
        res['cross_cai'] = {
            other: round(cai(res['dna'], build_weights(FREQS[other])), 3)
            for other in ORGANISM_CONFIGS if other != org
        }
    return results


# ==========================================================================
# 11. Demo
# ==========================================================================
if __name__ == '__main__':
    # 34-aa fragment of GFP, used as the standard test case
    DEMO = 'MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGE'

    for organism in ORGANISM_CONFIGS:
        for level in ('high', 'medium', 'low'):
            r = optimize(DEMO, organism, level)
            flags = (len(r['restriction_sites_found'])
                     + len(r['motifs_avoided_but_present']))
            fp = r['five_prime_optimization']
            print(f"{r['organism_label']:<40} {level:<7} "
                  f"CAI={r['cai']:.3f}  GC={r['gc']:.3f} "
                  f"{'OK ' if r['gc_ok'] else 'OUT'}  "
                  f"flags={flags}  5'dG={fp['final_dG']:>6}"
                  f"{' (tuned)' if fp['applied'] else '        '}  "
                  f"back-translation={'ok' if r['translates_back'] else 'FAIL'}")

    print()
    cmp = compare_organisms(DEMO)
    for org, r in cmp.items():
        print(f"{org}: native CAI {r['cai']:.3f}, "
              f"cross-scored {r['cross_cai']}")


# ==========================================================================
# LICENSING NOTES - read before publishing this as a web app
# ==========================================================================
# Safe to publish as-is:
#   * This file: standard library only (math, random, re) under the PSF
#     licence. No third-party code is vendored here.
#   * The genetic code and the restriction-enzyme recognition sequences:
#     scientific facts, not protectable expression.
#   * Streamlit: Apache-2.0. Biopython (if you add it): Biopython Licence,
#     BSD-style. Both fine for a public web app.
#
# Needs attention:
#   * Codon frequencies. The numbers come from the Kazusa database
#     (python_codon_tables, CC0, redistributes them). Individual facts are
#     not copyrightable, but the EU sui generis database right (Directive
#     96/9/EC, implemented in Germany as UrhG section 87a ff.) can protect a
#     substantial extraction from a database. Two complete species tables is
#     arguably substantial. Mitigation, in order of safety:
#       (a) cite Nakamura, Gojobori & Ikemura (2000) NAR 28:292 and state
#           the source in the app - done above and required either way;
#       (b) or compute your own tables from NCBI/Ensembl CDS downloads,
#           which makes the tables your own work product;
#       (c) or pull at runtime via python_codon_tables (CC0) rather than
#           hard-coding, which inherits its public-domain dedication.
#     For a university assignment (a) is proportionate. For anything
#     commercial, do (b).
#   * REBASE is the authoritative restriction-enzyme database and is free
#     for academic use only. Do NOT scrape it for an expanded enzyme panel
#     in a public app; use NEB's published recognition sequences or type
#     them in by hand, as here.
#   * ViennaRNA, if you later want real hairpin dG instead of the heuristic
#     in five_prime_structure_score: its licence permits research and
#     commercial use but forbids selling the package, requires credit to the
#     University of Vienna, and bundles files Debian classifies as non-free.
#     It is also a compiled dependency, awkward on Streamlit Cloud.
#     MIT-licensed `seqfold` is the cleaner drop-in if you want real folding.
#   * Do not copy UI text, figures or scoring formulas from commercial codon
#     optimisers (GenSmart, GeneOptimizer, OptimumGene). Those are
#     proprietary and some are patented.
