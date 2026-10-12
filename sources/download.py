"""Downloads the external sources to data/raw/<source>/.

    py -3 -m sources.download

- Open Targets Platform (CC0): drugs, mechanisms of action with direction (from
  ChEMBL), clinical indications with phase, FAERS adverse effects and
  gene-disease associations by evidence type.
- HGNC: Ensembl / Entrez / UniProt / gene symbol mapping.
- MONDO (SSSOM): disease mappings (EFO, DOID, OMIM, MeSH…) -> MONDO.
- STRING (CC BY 4.0): human protein interactions with a confidence score.
- CTD (non-commercial use): chemical-gene interactions with direction and chemical-disease relations.
"""
import re
import sys
import urllib.request
from pathlib import Path

RAW = Path(__file__).parent.parent / "data" / "raw"
OT_VERSION = "26.06"  # in 26.09 multiple sclerosis is inside a measurement term (README)
OT = f"https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/{OT_VERSION}/output"
OT_TABLES = [
    "drug_molecule",
    "drug_mechanism_of_action",
    "drug_warning",
    "clinical_indication",
    "clinical_target",
    "clinical_report",
    "openfda_significant_adverse_drug_reactions",
    "disease",
    "target",
    "association_by_datatype_direct",
    # evidence with direction (expression in the disease, effect of the variant on the gene)
    "evidence_expression_atlas",
    "evidence_gwas_credible_sets",
    "evidence_gene_burden",
    "evidence_eva",
    "evidence_uniprot_variants",
    "evidence_gene2phenotype",
    "evidence_orphanet",
    "evidence_clingen",
    "evidence_genomics_england",
    "evidence_impc",  # knockout mice: loss of the gene -> phenotype similar to the disease
    # how to attack a target (sources/targets.py): prioritisation factors and essentiality (DepMap)
    "target_prioritisation",
    "target_essentiality",
]
HGNC = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
MONDO_SSSOM = "http://purl.obolibrary.org/obo/mondo/mappings/mondo.sssom.tsv"
STRING_VERSION = "12.5"
STRING = {
    "links": f"https://stringdb-downloads.org/download/protein.links.v{STRING_VERSION}/9606.protein.links.v{STRING_VERSION}.txt.gz",
    "info": f"https://stringdb-downloads.org/download/protein.info.v{STRING_VERSION}/9606.protein.info.v{STRING_VERSION}.txt.gz",
}
# CTD (Comparative Toxicogenomics Database): non-commercial use, which is the case for this project.
CTD = [f"https://ctdbase.org/reports/{f}" for f in
       ("CTD_chem_gene_ixns.tsv.gz", "CTD_chemicals_diseases.tsv.gz", "CTD_chemicals.tsv.gz")]
HEADERS = {"User-Agent": "orevid/0.1"}


def open_url(url):
    return urllib.request.urlopen(urllib.request.Request(url, headers=HEADERS))


def fetch(url, dest):
    if dest.exists():
        return 0
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_suffix(dest.suffix + ".part")
    size = 0
    with open_url(url) as r, open(part, "wb") as f:
        while block := r.read(1 << 20):
            f.write(block)
            size += len(block)
    part.rename(dest)
    return size


def ot_table(name):
    with open_url(f"{OT}/{name}/") as r:
        files = sorted(set(re.findall(r'href="([^"/?]+\.parquet)"', r.read().decode())))
    total = 0
    for f in files:
        total += fetch(f"{OT}/{name}/{f}", RAW / "opentargets" / name / f)
    print(f"  {name}: {len(files)} files ({total / 1e6:.0f} MB new)", flush=True)


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    print(f"Open Targets {OT_VERSION}")
    for t in OT_TABLES:
        ot_table(t)
    (RAW / "opentargets" / "VERSION").write_text(OT_VERSION)
    print(f"HGNC: {fetch(HGNC, RAW / 'hgnc' / 'hgnc_complete_set.txt') / 1e6:.0f} MB")
    print(f"MONDO: {fetch(MONDO_SSSOM, RAW / 'mondo' / 'mondo.sssom.tsv') / 1e6:.0f} MB")
    for name, url in STRING.items():
        print(f"STRING {name}: {fetch(url, RAW / 'string' / url.rsplit('/', 1)[1]) / 1e6:.0f} MB")
    for url in CTD:
        print(f"CTD {url.rsplit('/', 1)[1]}: {fetch(url, RAW / 'ctd' / url.rsplit('/', 1)[1]) / 1e6:.0f} MB")


if __name__ == "__main__":
    main()
