# Open Targets: wrong «approved» indications (issues #4573 and #4579)

## What happened

- **#4573**: in Open Targets 26.09, multiple sclerosis (MONDO_0005301) was merged into a measurement
  term (EFO_0803536). OREVID uses 26.06 because of it.
- **#4579** (reported 12-10-2026): two kinds of false approvals in 26.06, lincomycin for multiple
  sclerosis (from a label with 3 molecules) and ADHD drugs approved for Alzheimer disease (the
  PMDA's «ad»).
- **Reply from David Ochoa (12-10-2026, 22:36)**, in short:
  1. The lincomycin record was curated by hand by ChEMBL, and the ChEMBL website also lists MS as an
     indication of lincomycin (CHEMBL1447). Open Targets respects curation over automation. A
     colleague (@vivienho) will follow up with ChEMBL. He asks: «if you have a script or a dataset
     that can identify more cases like this, it could be very valuable.»
  2. The acronym problem is the same bug as #4573: both are grounded with OnToma, which struggles
     with acronyms. They are discussing how to handle it; knowing it also affects the clinical
     corpus is useful.

## What I think

- It is a good reply: they confirm both errors, name the root causes (ChEMBL curation; OnToma
  grounding) and ask for exactly what we can provide. It is also a chance for OREVID to be useful
  to a major resource, and a natural way to become known by the Open Targets and ChEMBL teams.
- The lincomycin case is not a one-off. ChEMBL's curation attached every molecule mentioned in the
  FDA label 214962 (fingolimod lauryl sulfate) to multiple sclerosis, so mycophenolate sodium is
  «approved» for MS too. The same pattern (all approval reports for a pair list several molecules,
  none lists the drug alone) finds 729 pairs in 26.06. A few are clear errors of the same kind
  (berotralstat for Alzheimer disease from an EMA report that also contains flutemetamol F-18, an
  amyloid PET tracer). Most are components of fixed-dose combinations credited with the
  combination's indication: omeprazole/esomeprazole for arthritis (naproxen–esomeprazole),
  vitamins and minerals for anaemia (multivitamin labels), tetanus toxoid for hepatitis B
  (combined vaccine). Those are approved only as part of the combination: not wrong in the same
  sense, but misleading for anyone who reads «drug X approved for disease Y».
- The acronym problem is narrow and very clean: 14 pairs in 26.06 where every report gives the
  disease as ≤ 3 characters, all from the PMDA, 13 of them wrong: «AD/HD» split into «ad»
  (Alzheimer disease) and «hd» (Huntington disease) for four ADHD drugs; «sos» (sinusoidal
  obstruction syndrome) grounded as spondylo-ocular syndrome for defibrotide; «mac» (Mycobacterium
  avium complex) as microphthalmia with coloboma and mitral annular calcification for rifabutin;
  «cad» (probably congenital antithrombin deficiency) as cold agglutinin disease and coronary
  artery disorder for antithrombin gamma. Only «vwd» → von Willebrand disease is right. Not
  grounding ≤ 3-character strings without their expansion would remove all the wrong ones.

## What was done (12-10-2026, night)

- `sources/ot_suspect_approvals.py`: reads only the Open Targets 26.06 tables and writes
  `data/ext/ot_suspect_approvals.csv` (copy for sharing: `resources/opentargets/ot_suspect_approvals_26.06.csv`).
  Columns: pattern (`acronym` or `passenger`), ChEMBL and disease ids and names, the disease text in
  the report, the molecules in the report, which of them are approved alone for the disease, report
  ids, sources and URL, whether DrugBank (PrimeKG) has the indication, how many official reports
  carry the same set of molecules, and a review priority (`high`: acronyms, and passengers whose
  molecule set appears in a single report and that DrugBank does not have).
- Counts on 26.06: 11,223 approval-stage pairs; 14 acronym; 729 passenger (355 high priority);
  4,593 supported only by curated resources (TTD, ATC, …), not flagged here.
- Not checked by hand beyond the examples above: the passenger list needs a manual review.

## How to proceed (proposal)

1. Reply in #4579 with the draft below (the maintainer posts it; nothing has been posted).
2. Offer to rerun the script on each release and to adapt it to their schema; a pull request to
   their quality-control tooling could come later if they want it.
3. If they confirm the combination-component cases are unwanted, suggest a flag («approved as part
   of a fixed-dose combination») rather than deleting them.
4. OREVID itself needs no change: `sources/ot_approvals.py` already counts as confirmed only
   approvals backed by a single-drug report with the disease spelled out, so all the cases above
   are «unconfirmed» on the website.

## Draft reply for #4579 (to post from the maintainer's account)

> Thanks @d0choa, that helps. Here is a first pass on 26.06, in case it is useful.
>
> Script: `sources/ot_suspect_approvals.py` in https://github.com/orevid-org/orevid. Output for 26.06:
> `resources/opentargets/ot_suspect_approvals_26.06.csv` in the same repository. It only reads the Open
> Targets tables: for every approval-stage pair in `clinical_indication` it looks at the approval-level
> reports behind it and flags two patterns.
>
> **1. Acronyms** (14 pairs, all PMDA): every report gives the disease as ≤ 3 characters. 13 look wrong:
> «AD/HD» seems to be split into «ad» → Alzheimer disease and «hd» → Huntington disease for atomoxetine,
> guanfacine, lisdexamfetamine and methylphenidate; «sos» → spondylo-ocular syndrome for defibrotide
> (sinusoidal obstruction syndrome); «mac» → microphthalmia with coloboma and mitral annular
> calcification for rifabutin (Mycobacterium avium complex); «cad» → cold agglutinin disease and
> coronary artery disorder for antithrombin gamma (probably congenital antithrombin deficiency). Only
> «vwd» → von Willebrand disease is right.
>
> **2. Several molecules in every report** (729 pairs; 355 marked high priority, where the same set of
> molecules appears in only one report and DrugBank does not have the indication): all label or agency
> reports for the pair list several molecules and none lists the drug alone. This catches lincomycin and
> mycophenolate sodium for multiple sclerosis (FDA label 214962, fingolimod lauryl sulfate) and, for
> example, berotralstat for Alzheimer disease (an EMA Orladeyo report that also contains flutemetamol F-18).
> Most of the list, though, are components of fixed-dose combinations credited with the combination's
> indication (omeprazole/esomeprazole for arthritis from naproxen–esomeprazole, vitamins and minerals for
> anaemia, tetanus toxoid for hepatitis B from a combined vaccine). Whether those should count as
> approvals of the single component is your call; a flag might be enough.
>
> It is heuristic and I have only checked the examples above by hand, so the second list needs review.
> Happy to rerun it on new releases or adapt it if that is useful.
