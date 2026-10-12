"""Twin nodes: the same drug or the same disease represented by more than one node (e.g. the
DrugBank/PrimeKG node and the ChEMBL/Open Targets one when they could not be linked by identifier).
Detected in the second independent review: it made approved drugs show up under «no previous use» and hid
the trials of a disease (docs/DECISIONS.md, P22).

They are linked by normalised name (lowercase, without extra whitespace).
"""
from collections import defaultdict

import numpy as np


def _key(name):
    return " ".join(str(name).lower().split())


def twins(g, node_type):
    """{node: set of nodes of the same type with the same name (itself included)}."""
    groups = defaultdict(set)
    nodes = np.flatnonzero(g.types == node_type)
    for i in nodes:
        groups[_key(g.names[i])].add(int(i))
    return {i: groups[_key(g.names[i])] for i in nodes}


def twin_neighbors(g, nodes, rel):
    """Union of the neighbours through relation `rel` of all the given nodes."""
    res = set()
    for n in nodes:
        res.update(int(x) for x in g.neighbors(n, rel))
    return res


def names_with_use(g, rel):
    """(Normalised) names of drugs that have some neighbour through `rel` in any of their twins."""
    degree = np.diff(g.M(rel).indptr)
    drugs = np.flatnonzero((g.types == "drug") & (degree[: len(g.types)] > 0))
    return {_key(g.names[i]) for i in drugs}


key = _key


# Salts and esters: the same active moiety (expert review: «prednisolone» showed up as a candidate for IPF
# where «prednisolone phosphate» is already approved).
_SALTS = ("sodium phosphate", "disodium phosphate", "phosphate", "hydrochloride", "hydrobromide", "acetate",
          "sodium succinate", "hemisuccinate", "succinate", "dipropionate", "propionate", "valerate", "butyrate",
          "mesylate", "maleate", "tartrate", "bitartrate", "citrate", "besylate", "fumarate", "sulfate", "bromide",
          "chloride", "sodium", "potassium", "calcium", "disodium", "tromethamine", "hyclate", "monohydrate",
          "dihydrate", "trihydrate", "anhydrous", "pivalate", "furoate", "benzoate", "dienanthate", "cypionate",
          "enanthate", "undecanoate", "palmitate", "decanoate", "lysine", "meglumine", "magnesium")


_GENERICS = {"dimethyl", "diethyl", "monomethyl", "methyl", "ethyl"}  # «dimethyl fumarate» is the drug


def active_moiety(name):
    """Name without trailing salts or esters (lowercase). If it would be left too short, the whole name."""
    s = _key(name)
    changed = True
    while changed:
        changed = False
        for x in _SALTS:
            if s.endswith(" " + x) and len(s) - len(x) - 1 >= 4 and s[: -len(x) - 1] not in _GENERICS:
                s = s[: -len(x) - 1].rstrip(" ,")
                changed = True
    return s


# Prodrugs, active metabolites and enantiomers of a drug (2nd expert review: teriflunomide proposed
# for RA where leflunomide is approved; dexibuprofen, fosphenytoin…). Hand-written map + enantiomer
# prefixes. It only serves to recognise «already approved here», never to propose anything.
_FAMILY = {
    "teriflunomide": "leflunomide", "fosphenytoin": "phenytoin", "mycophenolic acid": "mycophenolate mofetil",
    "valaciclovir": "aciclovir", "acyclovir": "aciclovir", "valganciclovir": "ganciclovir", "famciclovir": "penciclovir",
    "prednisone": "prednisolone", "paliperidone": "risperidone", "desvenlafaxine": "venlafaxine",
    "desloratadine": "loratadine", "fexofenadine": "terfenadine", "lisdexamfetamine": "amphetamine",
    "dextroamphetamine": "amphetamine", "dexmethylphenidate": "methylphenidate", "armodafinil": "modafinil",
    "eszopiclone": "zopiclone", "levocetirizine": "cetirizine", "escitalopram": "citalopram",
    "esomeprazole": "omeprazole", "dexlansoprazole": "lansoprazole", "levofloxacin": "ofloxacin",
    "levalbuterol": "salbutamol", "albuterol": "salbutamol", "levosalbutamol": "salbutamol",
    "dexketoprofen": "ketoprofen", "dexibuprofen": "ibuprofen", "levomilnacipran": "milnacipran",
    "levoleucovorin": "leucovorin", "folinic acid": "leucovorin", "dexmedetomidine": "medetomidine",
    "carbidopa": "levodopa", "oxcarbazepine": "carbamazepine", "eslicarbazepine": "carbamazepine",
    "eslicarbazepine acetate": "carbamazepine", "tenofovir alafenamide": "tenofovir", "tenofovir disoproxil": "tenofovir",
    "colestipol": "colestipol", "methylprednisolone": "methylprednisolone", "hydroxychloroquine": "hydroxychloroquine",
    "norethindrone": "norethisterone", "calcitriol": "calcitriol", "dihydroergotamine": "ergotamine",
}


def drug_family(name):
    """Keys of the same active moiety: without salts, with the prodrug map and without enantiomer prefix."""
    p = active_moiety(name)
    keys = {p, _FAMILY.get(p, p)}
    for pre in ("dex", "levo", "lev", "es", "ar", "r-", "s-"):
        if p.startswith(pre) and len(p) - len(pre) >= 5:
            keys.add(_FAMILY.get(p[len(pre):], p[len(pre):]))
    return keys
