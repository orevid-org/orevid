# Launch plan: from a working tool to one that researchers use (written 23-10-2026)

Model work is in `docs/ROADMAP.md` and `docs/IMPROVEMENT_PLAN.md`; this file is about publishing OREVID and getting
it used. State of each channel and contacts: `docs/STATUS.md` («Done on 22-10-2026», outreach).

**Audience:** researchers and clinicians working on repurposing, not patients (no general press; disease pages stay
non-indexable; the first-visit notice stays).

**Main rule: fix before announcing.** The website still shows known errors (griseofulvin «approved» for
candidiasis, insulin as a candidate for multiple sclerosis…; STATUS next steps 9). A field expert who finds them on
launch day forms an opinion that is hard to change, so the visible errors go first.

## Phase 1 (late October – early November): ready to show (home PC)

1. **N4, robust parts** (`docs/DRUG_STATUS.md`, plan «Added 23-10-2026»): negated pairs (griseofulvin), withdrawn /
   revoked medicines (olaratumab), Open Targets approvals that a current label confirms (ublituximab, butoconazole,
   terconazole). Checked against the source pages before release.
2. **G1 (GenCC)**, registered 23-10 (`docs/IMPROVEMENT_PLAN.md`).
3. **Preprint v3** on Preprints.org: R1b and the new hit rates, N1 positive (first evidence that the shared genes
   carry some mechanism), N2 / N4 negative results, the fixed figure 3, G1 if done. The v2 no longer matches the
   website.
4. **Submit v3 to PCI Mathematical and Computational Biology** (preparation of 20-10 holds: figures from
   `evaluation/figures.py`, PCI conflict-of-interest wording, line numbers). Cover letter and 8 recommenders
   drafted on 29-10 (STATUS, that day). Release v1.3.0 (matches v3) on Zenodo since 30-10,
   doi:10.5281/zenodo.23800955; its prediction zip is still to be attached from the home PC.
5. Afterwards: v3 in README, `CITATION.cff`, the website and ORCID (as for v2); a reply with the v3 on the Open
   Targets forum topic.

## Phase 2 (about two weeks later): launch day

- **r/bioinformatics** post (to draft; check the subreddit's self-promotion rules first) and the **Bluesky thread**
  (`OREVID_anunci_Bluesky_Mastodon.txt` on the work PC's desktop; numbers updated to the v3), the same day; pin the
  thread to the profile. Mastodon: the same thread.
- **LinkedIn** post (to draft): the personal network and industry.
- Answer every comment in the first days; that matters more than the post itself.
- Bluesky beforehand: follow at most 5–10 accounts a day (the account was suspended once for bulk following).

## Phase 3 (first and second month): reasons to come back

- **OREVID-bench as an open leaderboard** on Hugging Face (`docs/ROADMAP.md`, item 3): a second announcement,
  «submit your model and compare it honestly», aimed at machine-learning groups.
- **A plain-language article** on what OREVID is and what it taught: why the popularity baseline matters, why
  shared genes are not a mechanism. Shared far more than a preprint.
- **Expert feedback as validation labels** (STATUS next steps 7): GitHub issues and e-mails turned into labelled
  pairs.
- Follow up the outreach: TxGNN authors, Every Cure, ReDO, Cures Within Reach, CURE ID, EURORDIS / REMEDi4ALL,
  Bioteque (IRB), Open Targets (forum and David Ochoa), the three awesome-list pull requests.

## Phase 4 (third month on): long-term credibility

| When | What |
|---|---|
| when the review ends (months) | PCI recommendation, then Peer Community Journal |
| December 2026 | Open Targets 26.12: check the multiple sclerosis merge fix before upgrading (STATUS, data pitfalls) |
| around April 2027 | **BOSC 2027** abstract (Open Bioinformatics Foundation, within ISMB); **JOSS** submission (repository public for 6 months) |
| October 2027 | **12-month check of the prospective snapshot** (`docs/PROSPECTIVE_PLAN.md`): which candidates of 12-10-2026 moved forward in trials; the strongest evidence OREVID can have |

## What to measure

GoatCounter visits to disease pages; Hugging Face and Zenodo downloads (GitHub clone counts are mostly bots: 575
clones and 247 cloners against 6 page views in the first two weeks); GitHub stars, forks, issues and referring
sites; replies from the groups contacted; later, citations. Compare the launch week with the weeks before.

## Account and housekeeping

- **GitHub two-factor authentication on `quimqu` before 27-11-2026**, or the account is restricted from making
  changes (including OREVID). The author does it; keep the recovery codes.
- The `quimqu` forks of the three awesome lists can be deleted once their pull requests are merged or closed.
- Commits in this repository carry no Claude co-author line (author's decision of 21-10; the earlier history is not
  rewritten because Zenodo and the prospective plan cite its tags).
