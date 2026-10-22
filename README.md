# Codon Optimiser

Reverse-translates a protein sequence into an expression-optimised DNA
sequence, with a QC report. Built as a Streamlit web app.

**Hosts:** *S. cerevisiae* · *E. coli* K-12 · *B. subtilis* 168 ·
*H. sapiens* / HEK293

```bash
pip install -r requirements.txt
python -m streamlit run app.py
```

## What it does

Protein in, DNA out, with the user choosing a target host, an expression
level (high / medium / low), optional extra restriction sites to avoid, and
optionally the vector's 5' UTR.

The pipeline is an explicit priority order:

| Rank | Property | Treatment |
|-----:|----------|-----------|
| 1 | Restriction sites, hard motifs | absolute - avoided during construction |
| 2 | Codon usage (CAI) | the objective, maximised within those limits |
| 3 | GC window | repaired afterwards; may spend CAI, never breaks rank 1 |
| 4 | 5' structure | tie-break only, capped at 0.05 CAI, only if there is a problem |
| 5 | Kozak, N-end rule, CpG/UpA, homopolymers | measured and reported, not optimised |

Every design is gated on back-translating to the exact input protein.

## Honest limitations

- **CAI is the objective, but the evidence for it is mixed.** Kudla et al.
  (2009) found 5' folding energy dominated expression in *E. coli* while CAI
  showed no significant effect. The counter-argument for yeast is Presnyak
  et al. (2015): codon optimality controls mRNA stability there, a mechanism
  absent from the *E. coli* work.
- **The 5' dG is a nearest-neighbour stem heuristic, not a folding
  prediction.** Treat it as a relative indicator. `seqfold` (MIT) would be
  the drop-in for real folding.
- **CAI, not tAI.** tAI weights by tRNA gene copy number and is more
  mechanistic.
- **Homopolymer runs >= 8 are reported, not avoided.** They are a real DNA
  synthesis risk.
- **Kozak +4 is flagged but not optimised**, even where residue 2 has a
  synonymous codon that would supply it.
- **In mammalian hosts, the GC window and the CpG check pull against each
  other.** GC-rich codons carry CpG. Both are reported; neither is
  reconciled.
- **Promoter, terminator, UTRs and copy number are out of scope**, and often
  matter more than codon choice.
- The *B. subtilis* N-end rule set is indicative; that pathway is less
  characterised than the *E. coli* one.

## Data sources

Codon usage: Kazusa Codon Usage Database, via
[python_codon_tables](https://github.com/Edinburgh-Genome-Foundry/python_codon_tables)
(CC0). Cite Nakamura, Gojobori & Ikemura (2000) *Nucleic Acids Res* 28:292.
CAI after Sharp & Li (1987). Example proteins are verified
character-for-character against their UniProt entries (P42212, Q9U6Y8,
P00359, P01308).

See the LICENSING NOTES at the foot of `codon_optimizer.py` before reusing
this outside a course - the EU database right is a live question once four
complete species tables are shipped together.

## Adding a host

Add one entry to `FREQS`, one to `ORGANISM_CONFIGS`, one to
`N_END_DESTABILISING`, and one branch to `second_codon_check`. Three
import-time checks (`_validate_freqs`, `_validate_configs`,
`_check_host_notes`) will tell you if you missed any of them. The UI picks
up the new host automatically.

## This is a teaching tool

Not validated for constructs going to synthesis. Verify anything before you
order it.
