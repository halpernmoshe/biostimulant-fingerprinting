"""Clean MANUSCRIPT_DRAFT_v2.md -> v3: remove em dashes and LLM-isms."""
import re, sys
from pathlib import Path
sys.stdout.reconfigure(encoding='utf-8')

src = Path('C:/Users/moshe/Dropbox/ISF 2025/agent_coordination/MANUSCRIPT_DRAFT_v2.md')
dst = Path('C:/Users/moshe/Dropbox/ISF 2025/agent_coordination/MANUSCRIPT_DRAFT_v3.md')

text = src.read_text(encoding='utf-8')
EM = '\u2014'

# ── EM DASH REPLACEMENTS (each handled by context) ────────────────────────────
replacements = [
    # Table cell
    ('PGPR volatile exposure \u2014 GMV bacterial mixture',
     'PGPR volatile exposure: GMV bacterial mixture'),
    # Methods 2.2 — long wrapped parenthetical
    ('values from this reference contrast \u2014 comprising all expressed genes, not restricted to statistically significant genes \u2014 was retained as the axis signature.',
     'values from this reference contrast, comprising all expressed genes (not restricted to statistically significant genes), was retained as the axis signature.'),
    # Methods 2.5
    ('were enumerated \u2014 C(*n*_total, *n*_treat) permutations \u2014 and',
     'were enumerated (C(*n*_total, *n*_treat) permutations), and'),
    # Results 3.3 — wrapped lines
    ('The direction of the N-axis effect \u2014 the\nprimary amino acid fingerprint after correction of the library-size artifact in GSE297649 \u2014 was\nconcordant',
     'The direction of the N-axis effect, which is the primary amino acid fingerprint after correction of the library-size artifact in GSE297649, was concordant'),
    # Results 3.4
    ('N +11.8%; GB03: P +13.7%, N +10.6%) \u2014 a dual-nutrient fingerprint mechanistically consistent',
     'N +11.8%; GB03: P +13.7%, N +10.6%), a dual-nutrient fingerprint mechanistically consistent'),
    ('that complexity predicts JA-mode \u2014 it is that the origin-based PGPR label',
     'that complexity predicts JA-mode; the finding is that the origin-based PGPR label'),
    ('produced by origin-based labeling \u2014 precisely\nthe heterogeneity the framework is designed to detect.',
     'produced by origin-based labeling, which is precisely the heterogeneity the framework is designed to detect.'),
    ('Fig. S3 \u2014 P-axis F = 52.8, N-axis F = 49.2 \u2014 are from a separate discriminant analysis',
     'Fig. S3 (P-axis F = 52.8, N-axis F = 49.2) are from a separate discriminant analysis'),
    ('Extended PGPR analyses \u2014 including nutritional context gating of ISR mode (GSE129395) and community composition effects on Bacillus fingerprints (GSE285599) \u2014 are presented in Supplementary Note S1.',
     'Extended PGPR analyses, including nutritional context gating of ISR mode (GSE129395) and community composition effects on Bacillus fingerprints (GSE285599), are presented in Supplementary Note S1.'),
    # Results 3.5
    ('PSK is a sulfated pentapeptide derived from tyrosine \u2014 an\namino acid-derived signaling molecule \u2014 and might be expected',
     'PSK is a sulfated pentapeptide derived from tyrosine, an amino acid-derived signaling molecule, and might be expected'),
    ('across all four comparisons \u2014 the largest\nABA signal of any biostimulant in the framework',
     'across all four comparisons, the largest ABA signal of any biostimulant in the framework'),
    ('fingerprint \u2014 high ABA, negative P, negative Fe, modest JA \u2014 places it in the ABA/stress class',
     'fingerprint (high ABA, negative P, negative Fe, modest JA) places it in the ABA/stress class'),
    # Results 3.6
    ('the evolutionarily conserved PSR gene core \u2014 genes whose phosphate-starvation function has been\nretained across the approximately 100 million year divergence between Arabidopsis and tomato.',
     'the evolutionarily conserved PSR gene core, comprising genes whose phosphate-starvation function has been retained across the approximately 100 million year divergence between Arabidopsis and tomato.'),
    # Results 3.7
    ('r = 0.059 (p = 5.5 \u00d7 10\u207b\u00b9\u2078) \u2014 honest but weak, as expected given the large differences in tissue, duration, and compound between the two experiments.',
     'r = 0.059 (p = 5.5 \u00d7 10\u207b\u00b9\u2078), a weak but real result expected given the large differences in tissue, duration, and compound between the two experiments.'),
    # Results 3.8
    ('Fingerprinting of a seventh biostimulant class \u2014 the beneficial endophytic fungus *Serendipita indica* \u2014 is presented in Supplementary Note S2.',
     'Fingerprinting of a seventh biostimulant class, the beneficial endophytic fungus *Serendipita indica*, is presented in Supplementary Note S2.'),
    # Discussion 4.1
    ('This approach has practical appeal \u2014 origin is verifiable and stable \u2014 but is poorly predictive of physiological outcome.',
     'This approach has practical appeal (origin is verifiable and stable) but is poorly predictive of physiological outcome.'),
    ('The axes we built \u2014 calibrated against known physiological states \u2014 serve as a common reference frame',
     'The axes, calibrated against known physiological states, serve as a common reference frame'),
    ('per treatment \u2014 sufficient to demonstrate directional consistency and within-class coherence, but insufficient to place formal confidence intervals on individual axis scores',
     'per treatment, sufficient to demonstrate directional consistency and within-class coherence but insufficient to place formal confidence intervals on individual axis scores'),
    ('The principal conclusions \u2014 that the PGPR origin class is transcriptomically incoherent, that PSK-origin misleads about mode of action, and that cross-species conservation ratios distinguish biostimulants with different translational potential \u2014 rest on large effect sizes',
     'The principal conclusions (that the PGPR origin class is transcriptomically incoherent, that PSK-origin misleads about mode of action, and that cross-species conservation ratios distinguish biostimulants with different translational potential) rest on large effect sizes'),
    # Discussion 4.2
    ('enabling simultaneous N and P acquisition \u2014 a dual-nutrient fingerprint unique to the SA-ISR subgroup.',
     'enabling simultaneous N and P acquisition, a dual-nutrient fingerprint unique to the SA-ISR subgroup.'),
    ("the plant's primary line of resistance \u2014 the opposite of the intended effect.",
     "the plant's primary line of resistance, the opposite of the intended effect."),
    ('Fingerprint-based product matching \u2014 selecting JA-ISR or SA-ISR products based on the anticipated pathogen pressure \u2014 is therefore not merely a precision-agriculture refinement but a basic agronomic precaution.',
     'Fingerprint-based product matching, selecting JA-ISR or SA-ISR products based on anticipated pathogen pressure, is therefore not merely a precision-agriculture refinement but a basic agronomic precaution.'),
    # Discussion 4.3
    ('The PHR1 validation result \u2014 94% of SynCom P-axis effect lost in phr1 phl1 \u2014 confirms PHR1 as the key effector node.',
     'The PHR1 validation result (94% of SynCom P-axis effect lost in phr1 phl1) confirms PHR1 as the key effector node.'),
    # Discussion 4.4
    ('while suppressing the P-axis (\u22127.8 to \u221211.3%) \u2014 the opposite direction from amino acid hydrolysates.',
     'while suppressing the P-axis (\u22127.8 to \u221211.3%), the opposite direction from amino acid hydrolysates.'),
    # Discussion 4.5
    ('suggests stress priming without ISR activation \u2014 a different mechanism from PGPR.',
     'suggests stress priming without ISR activation, a distinct mechanism from PGPR.'),
    # Discussion 4.6
    ('+0.79 (AMF root under salt stress \u2014 most stable) to \u22120.52 (silicon under pathogen infection \u2014 near-reversal)',
     '+0.79 (AMF root under salt stress, most stable) to \u22120.52 (silicon under pathogen infection, near-reversal)'),
    # Discussion 4.7
    ('in dicotyledonous crops \u2014 consistent with PHR1 orthologs performing the same function in all angiosperms',
     'in dicotyledonous crops, consistent with PHR1 orthologs performing the same function in all angiosperms'),
    ('and rice biostimulant datasets \u2014 a key next step.',
     'and rice biostimulant datasets, a key next step.'),
    # Discussion 4.8
    ('(PRJNA526343, PRJNA592062) \u2014 no new wet-lab work is needed.',
     '(PRJNA526343, PRJNA592062); no new wet-lab work is needed.'),
    ('The cosine similarity framework introduced here \u2014 quantifying fingerprint rotation angle between normal and stressed conditions \u2014 provides a ready metric for such experiments.',
     'The cosine similarity framework introduced here, which quantifies fingerprint rotation angle between normal and stressed conditions, provides a ready metric for such experiments.'),
    ('generates specific testable predictions \u2014 particularly for seaweed biostimulants, where species-specific fingerprints are predicted based on known active compounds.',
     'generates specific testable predictions, particularly for seaweed biostimulants, where species-specific fingerprints are predicted based on known active compounds.'),
    # Figure legends
    ('before declining toward baseline \u2014 consistent with a transient transcriptional burst',
     'before declining toward baseline, consistent with a transient transcriptional burst'),
    ('GMV and SynCom separate from diacetyl, GB03, and P. megaterium \u2014 demonstrating',
     'GMV and SynCom separate from diacetyl, GB03, and P. megaterium, demonstrating'),
    # Supplementary notes
    ('**Finding 1 \u2014 Third PGPR mode',
     '**Finding 1: Third PGPR mode'),
    # Word count footer (not in paper body, just note)
    ('~50 words) \u2014 marginally over limit',
     '~50 words) -- marginally over limit'),
]

v3 = text
not_found = []
for old, new in replacements:
    if old in v3:
        v3 = v3.replace(old, new)
    else:
        not_found.append(old[:70])

# ── LLM-ISMS ─────────────────────────────────────────────────────────────────
llm_fixes = [
    ('Importantly, ', ''),
    ('Importantly,\n', '\n'),
    ('Crucially, ', ''),
    ('Crucially,\n', '\n'),
    ('Furthermore, ', ''),
    ('Furthermore,\n', '\n'),
    ('It is worth noting that ', ''),
    ('It should be noted that ', ''),
    ('showcase', 'demonstrate'),
    ('pivotal', 'key'),
    ('delve into', 'examine'),
    ('intricate', 'complex'),
    (' tapestry', ''),
    ('In summary, ', ''),
    ('In conclusion, ', ''),
    ('This is consistent with', 'Consistent with'),
    # Header/footer
    ('# MANUSCRIPT DRAFT v1', '# MANUSCRIPT DRAFT v3'),
    ('## Plant Cell & Environment submission', '## Plant, Cell & Environment submission'),
    ('*MANUSCRIPT_DRAFT_v1.md compiled 2026-02-23*', '*MANUSCRIPT_DRAFT_v3.md compiled 2026-02-24*'),
]
for old, new in llm_fixes:
    v3 = v3.replace(old, new)

# ── REPORT ────────────────────────────────────────────────────────────────────
remaining = v3.count(EM)
print(f'Remaining em dashes: {remaining}')
if remaining:
    for m in re.finditer(EM, v3):
        ctx = v3[max(0, m.start()-60):m.end()+60]
        print(f'  ...{ctx}...')

if not_found:
    print(f'\nNot found ({len(not_found)}):')
    for s in not_found:
        print(f'  {repr(s)}')

dst.write_text(v3, encoding='utf-8')
print(f'\nSaved: {dst}')
print(f'Size: {len(v3):,} chars')
