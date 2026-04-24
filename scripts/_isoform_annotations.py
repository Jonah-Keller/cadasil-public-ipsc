# Shared FLT1 / FLT4 / FN1 isoform annotations (GENCODE v44 / Ensembl release 112).
# Single source of truth imported by 24_annotate_isoforms.py and
# 25_cross_cohort_isoform.py so annotations cannot drift between tables.
#
# Biological provenance:
#   - MANE Select for canonical transcripts
#   - Kendall & Thomas 1996 sFLT1-i13 (FLT1-204, intron-13 readthrough, poly-Ser)
#   - Thomas 2007 / Heydarian 2009 sFLT1-e15a (a.k.a. sFLT1-i15) is ~789 aa;
#     its closest Ensembl representation is FLT1-207 or FLT1-213 and neither
#     passes filterByExpr in this cohort, so it is NOT assayed here.
#   - VEGFR3 long vs short (FLT4-201 / FLT4-202) distinguished by 3' terminal
#     phosphotyrosine retention
#   - FN1 EDA / EDB splice state (cellular vs plasma fibronectin)

ANNOT = {
    # ══ FLT1 / VEGFR1 ══════════════════════════════════════════════════════
    "ENST00000282397": {
        "ensembl_name": "FLT1-201",
        "isoform_label": "VEGFR1 full-length (membrane)",
        "aa_length": 1338, "n_exons": 30, "biotype": "protein_coding",
        "notes": "MANE Select; signal peptide + 7 Ig-like domains + transmembrane + split kinase; receptor tyrosine kinase on endothelium"},
    "ENST00000541932": {
        "ensembl_name": "FLT1-204",
        "isoform_label": "sFLT1-i13 (soluble VEGFR1)",
        "aa_length": 733, "n_exons": 15, "biotype": "protein_coding",
        "notes": "Intron-13 readthrough. Diagnostic poly-serine C-terminus (...PSSSSSSSSSSS). Lacks TM + kinase; secreted. Classical Kendall & Thomas 1996 sFLT1; VEGF-A / PlGF decoy"},
    "ENST00000615840": {
        "ensembl_name": "FLT1-207",
        "isoform_label": "FLT1 truncated (KSTRNDCTTQSNVKH tail)",
        "aa_length": 687, "n_exons": 13, "biotype": "protein_coding",
        "notes": "Clean C-terminus ...KSTRNDCTTQSNVKH, no poly-Ser. Closest Ensembl representation of sFLT1-e15a / sFLT1-i15 (Thomas 2007) but filtered out at our depth"},
    "ENST00000639477": {
        "ensembl_name": "FLT1-209",
        "isoform_label": "alternate truncated soluble-like VEGFR1",
        "aa_length": 665, "n_exons": 14, "biotype": "protein_coding",
        "notes": "Clean C-terminus ...EEMLLPK (no poly-Ser). Ends after Ig-like domain region; distinct termination from sFLT1-i13. NOT sFLT1-e15a"},
    "ENST00000615611": {
        "ensembl_name": "FLT1-206",
        "isoform_label": "noncoding (CDS not defined)",
        "aa_length": None, "n_exons": 12, "biotype": "protein_coding_CDS_not_defined",
        "notes": "1476 nt, no annotated ORF; low abundance"},
    "ENST00000543394": {
        "ensembl_name": "FLT1-205",
        "isoform_label": "noncoding (CDS not defined)",
        "aa_length": None, "n_exons": 10, "biotype": "protein_coding_CDS_not_defined",
        "notes": "1275 nt, no annotated ORF; low abundance"},

    # ══ FLT4 / VEGFR3 ══════════════════════════════════════════════════════
    "ENST00000261937": {
        "ensembl_name": "FLT4-201",
        "isoform_label": "VEGFR3 long (canonical full-length)",
        "aa_length": 1363, "n_exons": 30, "biotype": "protein_coding",
        "notes": "MANE Select; retains 3 C-terminal phosphotyrosines (Y1337/Y1363...) required for full signaling; lymphangiogenic + angiogenic"},
    "ENST00000393347": {
        "ensembl_name": "FLT4-202",
        "isoform_label": "VEGFR3 short (alt C-terminus)",
        "aa_length": 1298, "n_exons": 30, "biotype": "protein_coding",
        "notes": "Alternative 3' splicing in terminal exon; lacks most distal C-terminal tyrosines found in the long form; attenuated signaling"},
    "ENST00000619105": {
        "ensembl_name": "FLT4-212",
        "isoform_label": "truncated VEGFR3 (premature stop)",
        "aa_length": 295, "n_exons": 28, "biotype": "protein_coding",
        "notes": "Short ORF ending early; likely translationally non-productive or regulatory"},
    "ENST00000424276": {
        "ensembl_name": "FLT4-203",
        "isoform_label": "noncoding (CDS not defined)",
        "aa_length": None, "n_exons": 14, "biotype": "protein_coding_CDS_not_defined",
        "notes": "2298 nt, no annotated ORF"},
    "ENST00000513527": {
        "ensembl_name": "FLT4-210",
        "isoform_label": "noncoding (CDS not defined)",
        "aa_length": None, "n_exons": 6, "biotype": "protein_coding_CDS_not_defined",
        "notes": "684 nt, no annotated ORF"},
    "ENST00000502603": {
        "ensembl_name": "FLT4-205",
        "isoform_label": "retained intron",
        "aa_length": None, "n_exons": 4, "biotype": "retained_intron",
        "notes": "Unspliced intron retained; likely nuclear / NMD substrate"},
    "ENST00000507059": {
        "ensembl_name": "FLT4-207",
        "isoform_label": "retained intron",
        "aa_length": None, "n_exons": 21, "biotype": "retained_intron",
        "notes": "Partial splicing intermediate"},
    "ENST00000510000": {
        "ensembl_name": "FLT4-208",
        "isoform_label": "retained intron",
        "aa_length": None, "n_exons": 2, "biotype": "retained_intron",
        "notes": "Short retained-intron variant"},
    "ENST00000514810": {
        "ensembl_name": "FLT4-211",
        "isoform_label": "retained intron",
        "aa_length": None, "n_exons": 4, "biotype": "retained_intron",
        "notes": "Short retained-intron variant"},

    # ══ FN1 / Fibronectin 1 ════════════════════════════════════════════════
    "ENST00000354785": {
        "ensembl_name": "FN1-203",
        "isoform_label": "FN1 canonical (MANE Select)",
        "aa_length": 2477, "n_exons": 46, "biotype": "protein_coding",
        "notes": "MANE Select reference; longest canonical V region; defines FN1 contig for downstream cross-referencing"},
    "ENST00000323926": {
        "ensembl_name": "FN1-201",
        "isoform_label": "FN1 cellular (EDA+ EDB+)",
        "aa_length": 2446, "n_exons": 47, "biotype": "protein_coding",
        "notes": "Includes both extra type-III domains EDA (EIIIA) and EDB (EIIIB); marker of activated / remodeling endothelium"},
    "ENST00000357867": {
        "ensembl_name": "FN1-205",
        "isoform_label": "FN1 plasma-like (EDA- EDB-)",
        "aa_length": 2176, "n_exons": 44, "biotype": "protein_coding",
        "notes": "Lacks both EDA and EDB; closest to hepatic plasma fibronectin"},
    "ENST00000336916": {
        "ensembl_name": "FN1-202",
        "isoform_label": "FN1 46-exon isoform",
        "aa_length": 2355, "n_exons": 46, "biotype": "protein_coding",
        "notes": "Single extra domain included (EDA or EDB); intermediate cellular/plasma"},
    "ENST00000446046": {
        "ensembl_name": "FN1-212",
        "isoform_label": "FN1 46-exon isoform",
        "aa_length": 2330, "n_exons": 46, "biotype": "protein_coding",
        "notes": "Single extra domain included"},
    "ENST00000359671": {
        "ensembl_name": "FN1-206",
        "isoform_label": "FN1 45-exon isoform",
        "aa_length": 2386, "n_exons": 45, "biotype": "protein_coding",
        "notes": "45-exon protein-coding FN1; V-region / IIICS variant"},
    "ENST00000421182": {
        "ensembl_name": "FN1-207",
        "isoform_label": "FN1 45-exon isoform",
        "aa_length": 2240, "n_exons": 45, "biotype": "protein_coding",
        "notes": "45-exon protein-coding FN1; V-region variant"},
    "ENST00000443816": {
        "ensembl_name": "FN1-211",
        "isoform_label": "FN1 45-exon isoform",
        "aa_length": 2265, "n_exons": 45, "biotype": "protein_coding",
        "notes": "45-exon protein-coding FN1"},
    "ENST00000432072": {
        "ensembl_name": "FN1-209",
        "isoform_label": "FN1 45-exon isoform",
        "aa_length": 2267, "n_exons": 45, "biotype": "protein_coding",
        "notes": "45-exon protein-coding FN1"},
    "ENST00000356005": {
        "ensembl_name": "FN1-204",
        "isoform_label": "FN1 44-exon isoform",
        "aa_length": 2296, "n_exons": 44, "biotype": "protein_coding",
        "notes": "44-exon protein-coding FN1"},
    "ENST00000426059": {
        "ensembl_name": "FN1-208",
        "isoform_label": "FN1 truncated (N-terminal)",
        "aa_length": 657, "n_exons": 13, "biotype": "protein_coding",
        "notes": "Short N-terminal fragment; ends early; regulatory or secreted truncation"},
    "ENST00000438981": {
        "ensembl_name": "FN1-210",
        "isoform_label": "FN1 truncated (very short)",
        "aa_length": 241, "n_exons": 4, "biotype": "protein_coding",
        "notes": "241 aa fragment; unlikely to assemble into fibrillar FN"},
    "ENST00000494446": {"ensembl_name": "FN1-225", "isoform_label": "retained intron",
        "aa_length": None, "n_exons": 4, "biotype": "retained_intron", "notes": ""},
    "ENST00000461974": {"ensembl_name": "FN1-215", "isoform_label": "retained intron",
        "aa_length": None, "n_exons": 3, "biotype": "retained_intron", "notes": ""},
    "ENST00000492816": {"ensembl_name": "FN1-224", "isoform_label": "retained intron",
        "aa_length": None, "n_exons": 19, "biotype": "retained_intron", "notes": "Long retained-intron variant"},
    "ENST00000480024": {"ensembl_name": "FN1-220", "isoform_label": "retained intron",
        "aa_length": None, "n_exons": 2, "biotype": "retained_intron", "notes": ""},
    "ENST00000469569": {"ensembl_name": "FN1-216", "isoform_label": "retained intron",
        "aa_length": None, "n_exons": 2, "biotype": "retained_intron", "notes": ""},
    "ENST00000473614": {"ensembl_name": "FN1-218", "isoform_label": "retained intron",
        "aa_length": None, "n_exons": 2, "biotype": "retained_intron", "notes": ""},
    "ENST00000498719": {"ensembl_name": "FN1-227", "isoform_label": "retained intron",
        "aa_length": None, "n_exons": 2, "biotype": "retained_intron", "notes": ""},
    "ENST00000474036": {"ensembl_name": "FN1-219", "isoform_label": "retained intron",
        "aa_length": None, "n_exons": 5, "biotype": "retained_intron", "notes": ""},
    "ENST00000496542": {"ensembl_name": "FN1-226", "isoform_label": "retained intron",
        "aa_length": None, "n_exons": 2, "biotype": "retained_intron", "notes": ""},
    "ENST00000485567": {"ensembl_name": "FN1-222", "isoform_label": "retained intron",
        "aa_length": None, "n_exons": 2, "biotype": "retained_intron", "notes": ""},
    "ENST00000460217": {"ensembl_name": "FN1-214", "isoform_label": "retained intron",
        "aa_length": None, "n_exons": 4, "biotype": "retained_intron", "notes": ""},
    "ENST00000471193": {"ensembl_name": "FN1-217", "isoform_label": "retained intron",
        "aa_length": None, "n_exons": 2, "biotype": "retained_intron", "notes": ""},
}


def strip_version(tid: str) -> str:
    return tid.split(".")[0]
