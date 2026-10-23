# N4. Drug status audit against current labels — first version, not released

Protocol and release gate: `docs/IMPROVEMENT_PLAN.md`, N4 (registered 22-10-2026). Code: `sources/drug_status.py`
(`py -3 -m sources.drug_status build` → `data/ext/drug_status.parquet`, `data/ext/pair_label.parquet`), sample for the
manual check: `evaluation/drug_status_check.py`. Website code ready but switched off (`web/build.py` `N4_RELEASED`).

## Sources (downloaded 22-10-2026)

- **openFDA drug labels** (export of 22-10-2026, 14 files, 1.8 GB; US public domain): «indications and usage» of
  every current SPL label; 17,008 active moieties have a label with an indication text.
- **EMA medicines table** (already downloaded for I3): medicine status, withdrawal / refusal dates, therapeutic
  indication text (centrally authorised medicines only).
- **FDA Orange Book and Purple Book** (I3): marketed vs discontinued products.
- DrugCentral was not needed for this version (its licence, CC BY-SA 4.0, would have to be checked before publishing
  derived fields).

## Coverage (approved drugs in DrugBank, 2,596)

US status found for 79.5 % (a marketed product, a current FDA label, or only discontinued products), EU status for
27.3 %, a current FDA label with indications for 65.6 %, an EMA indication text for 23.3 %.

## Pairs the website calls «approved», by category (version checked)

| Source | Confirmed | Negated | Only in a combination product | Not found in a current label | No current label |
|---|---|---|---|---|---|
| PrimeKG indications | 2,967 | 26 | 285 | 3,918 | 2,192 |
| Open Targets approvals, confirmed | 3,627 | 30 | 101 | 1,219 | 413 |
| Open Targets approvals, «unconfirmed» | 606 | 6 | 259 | 1,581 | 2,373 |

(All pairs, including diseases without a website page.)

## The cases found on 22-10

| Pair | Website today | N4 |
|---|---|---|
| griseofulvin — candidiasis | approved (Open Targets) | **negated** («not effective in … candidiasis») |
| ublituximab — multiple sclerosis and its subtypes | unconfirmed | **confirmed** (BRIUMVI; subtypes through the broader term) |
| butoconazole, terconazole — candidiasis | unconfirmed | **confirmed** (Gynazole-1, terconazole labels) |
| zinc oxide — cutaneous candidiasis | approved (PrimeKG) | **only in a combination product** (VUSION, miconazole / zinc oxide) |
| benzoic acid — systemic sclerosis | approved (PrimeKG) | **not found in a current label** |
| olaratumab | approved | **no current label**; US discontinued; EU authorisation revoked |

All six come out right.

## Manual check (gate: ≤ 5 % wrong)

Two random samples of 50 pairs (`data/eval/drug_status_sample.jsonl`, seeds 1 and 2), each with the label sentence
behind the decision, were checked one by one by an independent reviewer (an AI agent acting as a clinical
pharmacologist, without web search; its verdicts below).

| Sample | Correct | Wrong | Unsure |
|---|---|---|---|
| 50 «approved» pairs | 40 | 9 | 1 |
| 50 pairs that change category | 31 | 14 | 5 |
| **Total** | **71** | **23** | **6** |

**Decision: not released** (23 % wrong against the 5 % gate). The website keeps the approvals as before; nothing of
N4 is shown.

**Deviation from the registered check (noted 23-10-2026, work PC).** The plan (`docs/IMPROVEMENT_PLAN.md`, N4,
«Checks before release») asked for a manual check «against the source pages». This one was done by an AI agent
without web search, from its own knowledge and the label sentence the sample carried. Its verdicts can be wrong in
both directions (some «wrong» may be right and the other way round). It does not change tonight's decision (nothing
was released, and 23 % is far above the gate), but **it does not count as the release check**: before anything of
N4 reaches the website, each sample is checked by a person with the source pages open (the current FDA label on
DailyMed, the EMA product page, the Orange / Purple Book entry), or at least by an agent that opens those pages and
quotes the sentence it relies on; the record of each verdict keeps the URL of the page used.

### Where it goes wrong

- **Synonyms missed** (largest group): hypertensive disorder / hypertension, ocular / intraocular hypertension,
  tinea corporis / ringworm, tinea cruris / jock itch, nasopharyngitis / colds, IgA glomerulonephritis / IgA
  nephropathy, urinary tract infection (plural in the label).
- **Narrower diseases covered by a broader label term** more than one ontology step away: schizophreniform disorder
  by «psychotic disorders», an Alzheimer subtype by «Alzheimer's disease», ovarian mucinous adenocarcinoma by
  «carcinoma of the ovary», Hodgkin lymphoma by «leukemias and lymphomas», arteriosclerosis by «coronary
  atherosclerosis».
- **Other forms of the same drug:** triamcinolone acetonide (fixed after the check), injectable betamethasone,
  the formulated form of bortezomib.
- **«Confirmed» on a sentence that is not an indication:** diltiazem — ventricular tachycardia (a clinical-study
  sentence), empagliflozin — type 1 diabetes (the label names type 2).
- A combination product not recognised (glecaprevir / Mavyret).

### What a second version needs (to register before it is checked again)

1. Disease matching on the label's own wording: MONDO / UMLS synonyms and plural forms, and the label term mapped to
   a disease with all its descendants (not one step), so that «psychotic disorders» covers schizophreniform disorder.
2. Only the indication sentences («indicated for …», «treatment of …»), not clinical-study or limitation text, and
   the diabetes type (or any qualifier) taken into account.
3. Drug matching by UNII / RxNorm ingredient (openFDA has both) instead of names.
4. Show first only the robust parts, each with its own check: «negated» (all of them reviewed by hand: 62 pairs),
   withdrawn / revoked / refused EU medicines, and «unconfirmed» Open Targets approvals that a label confirms.

## Verdicts of the manual check

| # | Drug | Disease | Category | Verdict | Reason |
|---|---|---|---|---|---|
| 1 | Testolactone | female breast carcinoma | no current label | correct | discontinued in the US around 2008; no EU authorisation |
| 2 | Aminolevulinic acid | skin disease | not found | correct | labels name actinic keratosis and glioma imaging; the pair is too broad |
| 3 | Glasdegib | unclassified acute myeloid leukemia | confirmed | correct | label covers newly diagnosed AML |
| 4 | Benzylpenicillin | inhalational botulism | not found | correct | botulism is treated with antitoxin |
| 5 | Triprolidine | familial dermatographia | not found | correct | only OTC cold / allergy combinations |
| 6 | Cetyl alcohol | newborn respiratory distress syndrome | not found | correct | excipient; Exosurf discontinued |
| 7 | Albiglutide | type 2 diabetes mellitus | no current label | correct | withdrawn in the US and EU (2017–2018) |
| 8 | Diltiazem | ventricular tachycardia | confirmed | **wrong** | clinical-study sentence; IV diltiazem is contraindicated in VT |
| 9 | Captopril | myocardial infarction | confirmed | correct | survival after MI with LV dysfunction |
| 10 | Amphotericin B | esophageal candidiasis | negated | correct | «should not be used» in non-neutropenic patients |
| 11 | Pemetrexed | malignant pleural mesothelioma | confirmed | correct | explicit indication |
| 12 | Sarecycline | acne | confirmed | correct | acne vulgaris |
| 13 | Stannous fluoride | dental caries | confirmed | correct | prevention of dental caries |
| 14 | Epitizide | hypertensive disorder | no current label | correct | not marketed in the US; no central EU authorisation |
| 15 | Naloxone | substance abuse / dependence | only in a combination | correct | buprenorphine / naloxone |
| 16 | Mometasone furoate | primary cutaneous T-cell lymphoma | not found | correct | off-label |
| 17 | Triamcinolone | acquired thrombocytopenia | no current label | **wrong** | Kenalog-40 (acetonide) lists it; ingredient match missed (fixed after the check) |
| 18 | Oxetorone | migraine disorder | no current label | correct | French national product only |
| 19 | Chlorpromazine | schizophreniform disorder | not found | **wrong** | «psychotic disorders» covers it |
| 20 | Grazoprevir | chronic hepatitis C | only in a combination | correct | Zepatier |
| 21 | Insulin human | type 1 diabetes mellitus | confirmed | correct | «diabetes mellitus» as broader term |
| 22 | Dexamethasone | systemic lupus erythematosus | confirmed | correct | named among collagen diseases |
| 23 | Telmisartan | stroke disorder | confirmed | correct | CV risk reduction (snippet misleading) |
| 24 | Almasilate | gastroesophageal reflux disease | no current label | correct | antacid without FDA / EMA label |
| 25 | Acetylsalicylic acid | nasopharyngitis | not found | **wrong** | «pain and fever of colds»: synonym missed |
| 26 | Bortezomib d-mannitol | mantle cell lymphoma | no current label | **wrong** | formulated bortezomib (Velcade) lists it |
| 27 | Lasofoxifene | osteoporosis | no current label | correct | EU authorisation lapsed; never US |
| 28 | Capecitabine | breast cancer | confirmed | correct | explicit |
| 29 | Betamethasone | Hodgkin's paragranuloma | not found | unsure | «leukemias and lymphomas» (palliative), generic |
| 30 | Diltiazem | hypertension | confirmed | correct | explicit |
| 31 | Carteolol | ocular hypertension | not found | **wrong** | «intraocular hypertension»: synonym |
| 32 | Cenobamate | partial epilepsy | confirmed | correct | focal-onset seizures |
| 33 | Tinidazole | amebiasis | confirmed | correct | explicit |
| 34 | Omeprazole | gastroesophageal reflux disease | confirmed | correct | explicit |
| 35 | Tolnaftate | tinea pedis | confirmed | correct | athlete's foot |
| 36 | Hydrocortisone | classic Hodgkin lymphoma | not found | **wrong** | «palliative management of leukemias and lymphomas» |
| 37 | Camostat | esophagitis | no current label | correct | Japan only |
| 38 | Pheniramine | papillary conjunctivitis | not found | correct | allergic itchy eyes only |
| 39 | Zoledronic acid | osteoporosis | confirmed | correct | explicit |
| 40 | Eplerenone | myocardial infarction | confirmed | correct | heart failure after acute MI |
| 41 | Clonidine | ADHD | confirmed | correct | extended-release clonidine |
| 42 | Vernakalant | atrial fibrillation | confirmed | correct | Brinavess |
| 43 | Guanfacine | hypertensive disorder | not found | **wrong** | «hypertension»: synonym |
| 44 | Norfloxacin | salmonellosis | no current label | correct | Noroxin discontinued 2014 |
| 45 | Naltrexone | substance abuse / dependence | confirmed | correct | alcohol / opioid dependence |
| 46 | Simvastatin | myocardial infarction | confirmed | correct | reduces risk of non-fatal MI |
| 47 | Cycrimine | Parkinson disease | no current label | correct | discontinued long ago |
| 48 | Butenafine | tinea corporis | not found | **wrong** | «ringworm»: synonym |
| 49 | Prednisolone | posterior uveitis | confirmed | correct | explicit |
| 50 | Dexamethasone sodium phosphate | epicondylitis | confirmed | correct | listed under rheumatic disorders |
| 51 | Influenza A virus, inactivated whole virion | influenza | no current label | unsure | US vaccines are split or subunit |
| 52 | Meclocycline | acute gonococcal cervicitis | no current label | correct | discontinued |
| 53 | Somatropin | endocrine system disease | not found | correct | labels name GH deficiency, not the generic term |
| 54 | Cefixime | urinary tract infection | not found | **wrong** | «Urinary Tract Infections» (plural) |
| 55 | Isavuconazonium | white piedra | not found | correct | aspergillosis and mucormycosis only |
| 56 | Glecaprevir | hepatitis C virus infection | not found | **wrong** | should be «only in a combination product» (Mavyret) |
| 57 | Romiplostim | aplastic anemia | not found | correct | approved in Japan only |
| 58 | Eptacog alfa (activated) | factor VII deficiency | not found | **wrong** | NovoSeven lists congenital FVII deficiency |
| 59 | Droxidopa | hemiparkinsonism-hemiatrophy syndrome | not found | correct | neurogenic orthostatic hypotension only |
| 60 | Norfloxacin | maxillary sinusitis | no current label | correct | discontinued |
| 61 | Tolcapone | hepatic failure | not found | correct | adverse effect, not an indication |
| 62 | Empagliflozin | type 1 diabetes mellitus | confirmed | **wrong** | labels name type 2 only |
| 63 | Iron | anemia of prematurity | not found | correct | only homeopathic hits |
| 64 | Cefapirin | bacterial arthritis | no current label | correct | veterinary only now |
| 65 | Tolazamide | type 2 diabetes mellitus | no current label | correct | discontinued in the US |
| 66 | Alipogene tiparvovec | familial apolipoprotein C-II deficiency | no current label | correct | Glybera withdrawn 2017 |
| 67 | Sulfisoxazole | ornithosis | no current label | correct | discontinued |
| 68 | Clonazepam | early myoclonic encephalopathy | not found | correct | not this syndrome |
| 69 | Aztreonam | gram-negative bacterial infections | confirmed | correct | snippet misleading, indication right |
| 70 | Betamethasone | adrenogenital syndrome | not found | **wrong** | Celestone Soluspan lists congenital adrenal hyperplasia |
| 71 | Methapyrilene | papillary conjunctivitis | no current label | correct | withdrawn 1979 |
| 72 | Benzoic acid | dermatomyositis | not found | correct | excipient / homeopathic |
| 73 | Fluocinonide | lupus erythematosus | not found | correct | generic dermatoses only |
| 74 | Liothyronine | familial thyroid carcinoma | not found | unsure | TSH suppression in differentiated thyroid cancer |
| 75 | Lindane | mite infestation | no current label | unsure | largely discontinued |
| 76 | Lovastatin | arteriosclerosis disorder | not found | **wrong** | «coronary atherosclerosis» |
| 77 | Ketorolac | papillary conjunctivitis | not found | correct | seasonal allergic conjunctivitis only |
| 78 | Butenafine | dermatophytosis of groin and perianal area | not found | **wrong** | «jock itch»: synonym |
| 79 | Cortisone acetate | fetal erythroblastosis | not found | correct | never an indication |
| 80 | Rivastigmine | Alzheimer disease without neurofibrillary tangles | not found | **wrong** | «Alzheimer's disease» covers the subtype |
| 81 | Chlordiazepoxide | peptic ulcer disease | only in a combination | correct | Librax |
| 82 | Amifampridine | nervous system disorder | not found | correct | LEMS only |
| 83 | Estramustine | prostate cancer | no current label | unsure | may have been discontinued recently |
| 84 | Iptacopan hydrochloride | IgA glomerulonephritis | not found | **wrong** | «IgA nephropathy»: synonym |
| 85 | Triamcinolone | adrenocortical insufficiency | no current label | **wrong** | Kenalog-40 lists it (fixed after the check) |
| 86 | Magnesium lactate | dyspepsia | no current label | correct | supplement salt |
| 87 | Paclitaxel | ovarian mucinous adenocarcinoma | not found | **wrong** | «carcinoma of the ovary» covers it |
| 88 | Methotrexate | female breast carcinoma | not found | **wrong** | injection labels list breast cancer |
| 89 | Etretinate | parapsoriasis | no current label | correct | withdrawn 1998 |
| 90 | Tacrine | Alzheimer disease | no current label | correct | withdrawn 2013 |
| 91 | Elacestrant | breast cancer | confirmed | correct | explicit |
| 92 | Olmutinib | non-small cell lung carcinoma | no current label | correct | Korea only, withdrawn |
| 93 | Atorvastatin | Prinzmetal angina | only in a combination | correct | literally true; the indication belongs to amlodipine |
| 94 | Dexamethasone sodium phosphate | herpes simplex infectious disease | not found | correct | contraindicated in ocular herpes simplex |
| 95 | Besifloxacin | acute contagious conjunctivitis | not found | unsure | bacterial vs viral «pink eye» |
| 96 | Amifampridine | autoimmune disease | not found | correct | LEMS only |
| 97 | Enviomycin | pulmonary tuberculosis | no current label | correct | Japan only |
| 98 | Methyldopa | hypertensive disorder | not found | **wrong** | «hypertension»: synonym |
| 99 | Simvastatin | hyperlipoproteinemia | not found | **wrong** | Fredrickson types in the label |
| 100 | Cefalotin | streptococcal pneumonia | no current label | correct | discontinued in the US |

## Second version (23-10-2026, home PC): release in parts, checked against the source pages

Changes (plan «Added 23-10-2026»): labels re-extracted with their SPL set id (DailyMed link) and the EMA product page;
a confirmation needs an indication sentence («indicated», «treatment of»…) that respects type qualifiers (type 1 /
type 2 diabetes); plurals; more ester / salt forms (triamcinolone acetonide). The six cases of 22-10 still come out
right. Each part was checked by an independent agent that **opened every source page** (DailyMed SPL, EMA product page,
DailyMed search) and kept the URL and the quoted sentence of each verdict (`data/eval/n4_part_{a,b,c}.jsonl`, verdicts
in the session record).

| Part | Pairs / drugs | Checked | Wrong | Gate (≤ 5 %) |
|---|---|---|---|---|
| (a) «negated» pairs | 61 pairs | all 61 | **31 (51 %)** | fails |
| (b) current status of drugs no longer on the market | 504 drugs | 50 | **3 (6 %)** | fails |
| (c) «unconfirmed» Open Targets approvals a label confirms | 545 pairs | 50 | **4 (8 %)** | fails |

- **(a)** The negation pattern is not reliable: a limitation about another use for the same disease («not indicated for
  induction of immune tolerance» in haemophilia B) or another population («not recommended for patients with end-stage
  renal disease») was read as a negation, and indications listed as bare items (cyclophosphamide: myeloma,
  retinoblastoma, breast carcinoma) were outweighed by a negated sentence elsewhere. **Not releasable as a rule.**
  The 25 pairs the check confirmed as clear negations of the drug for that disease (borderline ones and a product
  mismatch left out) are kept as a curated list with their quotes and URLs: `resources/n4_verified_negations.csv`
  (griseofulvin – candidiasis, insulins – diabetic ketoacidosis, procaine benzylpenicillin – syphilis and other
  treponematoses, hCG excluded, …).
- **(b)** The EU side was right in every checked case. The three errors are US: «no marketed product» for drugs with a
  current label on DailyMed that the openFDA export lacks (abametapir / Xeglyze; inotersen / Tegsedi, two nodes).
- **(c)** Four errors: a broader term qualified by the label («relapsing forms of multiple sclerosis» does not cover
  progressive MS: fingolimod), a negation not caught («has not been evaluated in the treatment of infections»:
  minocycline), and two mentions in an indication sentence that are not the indication (hydrochlorothiazide –
  diabetes, doxepin – anxiety).

Next: (b) only claim «no US product» when DailyMed's own search returns nothing; (c) the three fixes; then a new check
of 50 new random rows each, against the source pages, before release.
