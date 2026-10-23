"use strict";
// OREVID: explorer of drug-repurposing candidates (static website, data in data/).
// Every user-visible text lives in i18n/<lang>.json (UI strings, the Method and Legal pages, and a
// «terms» map that translates the English data values shown to users).

// --- language -----------------------------------------------------------------------------
const LANGS = ["en", "es", "ca"];
const LOCALES = {};  // lang -> Promise of the prepared locale (cached)
let lang = (() => {
  let stored = null;
  try { stored = localStorage.getItem("lang"); } catch (e) { /* no storage */ }
  // English by default for everyone; Spanish and Catalan only when the visitor picks them (remembered)
  return LANGS.includes(stored) ? stored : "en";
})();
let L = null;  // current locale

// Whole-word replacement, longest keys first («activates» is not «activate» + s)
function termsRegex(map) {
  const keys = Object.keys(map);
  if (!keys.length) return null;
  return new RegExp("(?<!\\p{L})(" + keys.sort((a, b) => b.length - a.length)
    .map((k) => k.replace(/[.*+?^${}()|[\]\\]/g, "\\$&")).join("|") + ")(?!\\p{L})", "gu");
}
function prepareLocale(j) {
  const terms = j.terms || {}, targetTerms = {...terms, ...(j.target_terms || {})};
  return {strings: j.strings || {}, terms, targetTerms, termsRe: termsRegex(terms), targetRe: termsRegex(targetTerms),
    method: j.method || "", legal: j.legal || "", about: j.about || ""};
}
function loadLocale(l) {
  return (LOCALES[l] ||= fetch(`i18n/${l}.json`, {cache: "no-cache"})
    .then((r) => { if (!r.ok) throw new Error(`HTTP ${r.status}: i18n/${l}.json`); return r.json(); })
    .then(prepareLocale)
    .catch((e) => { delete LOCALES[l]; throw e; }));
}

const t = (k) => (L && L.strings[k]) ?? k;
const replaceTerms = (s, map, re) => (s == null || !re ? s : String(s).replace(re, (m) => map[m]));
// Data values (status, phases, actions, reasons) shown in the UI language
const tr = (s) => (L ? replaceTerms(s, L.terms, L.termsRe) : s);
// Targets table values (some differ in gender in Catalan)
const trTarget = (s) => (L ? replaceTerms(s, L.targetTerms, L.targetRe) : s);
// Evidence source of a target (internal Open Targets and PrimeKG values)
const evTarget = (s) => String(s || "").split(", ").map((x) => (L && L.targetTerms[x]) || trTarget(x)).join(", ");
const esc = (s) => String(s ?? "").replace(/[&<>"']/g, (c) => ({"&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;"}[c]));

// --- data ---------------------------------------------------------------------------------
let INDEX = null, FRAG = 500, INDEX_PROMISE = null;
const cache = {};

async function loadIndex() {
  const r = await fetch("data/index.json", {cache: "no-cache"});
  const j = await r.json();
  FRAG = j.fragments;
  INDEX = j.diseases.map(([id, name, hasData]) => ({id, name, hasData, lower: name.toLowerCase()}));
}

async function diseasePage(id) {
  const f = id % FRAG;
  if (!cache[f]) cache[f] = fetch(`data/d/${f}.json`, {cache: "no-cache"}).then((r) => (r.ok ? r.json() : {}));
  const d = await cache[f];
  return d[String(id)];
}

// --- search -------------------------------------------------------------------------------
const q = document.getElementById("q"), sug = document.getElementById("suggestions");
let sel = -1, found = [];
function search(text) {
  const s = text.trim().toLowerCase();
  if (!INDEX || s.length < 2) return [];
  const starts = [], contains = [];
  for (const m of INDEX) {
    if (m.lower.startsWith(s)) starts.push(m);
    else if (m.lower.includes(s)) contains.push(m);
    if (starts.length > 40) break;
  }
  const order = (a, b) => (b.hasData - a.hasData) || (a.name.length - b.name.length);
  return starts.sort(order).concat(contains.sort(order)).slice(0, 12);
}
function drawSuggestions() {
  sug.innerHTML = found.map((m, i) =>
    `<li role="option" data-id="${m.id}" aria-selected="${i === sel}"><span>${esc(m.name)}</span>` +
    (m.hasData ? "" : `<span class="no-data">${t("no_data_short")}</span>`) + "</li>").join("");
  sug.style.display = found.length ? "block" : "none";
}
q.addEventListener("input", () => { found = search(q.value); sel = found.length ? 0 : -1; drawSuggestions(); });
q.addEventListener("keydown", (e) => {
  if (e.key === "ArrowDown") { sel = Math.min(sel + 1, found.length - 1); drawSuggestions(); e.preventDefault(); }
  else if (e.key === "ArrowUp") { sel = Math.max(sel - 1, 0); drawSuggestions(); e.preventDefault(); }
  else if (e.key === "Enter" && found[sel]) { location.hash = `#/d/${found[sel].id}`; }
  else if (e.key === "Escape") { sug.style.display = "none"; }
});
sug.addEventListener("click", (e) => { const li = e.target.closest("li"); if (li) location.hash = `#/d/${li.dataset.id}`; });
document.addEventListener("click", (e) => { if (!e.target.closest(".search-box")) sug.style.display = "none"; });

// --- disease view -------------------------------------------------------------------------
const list = (xs) => (xs && xs.length ? xs.map(esc).join(", ") : t("none"));
const ext = (url, text, title) => `<a href="${url}" target="_blank" rel="noopener"${title ? ` title="${esc(title)}"` : ""}>${text}</a>`;

// One explanation line (in a candidate's detail). The context ones are not validated evidence.
function explanation(e) {
  switch (e.t) {
    case "target": return `${t("e_target")} ${e.g.map(([g, a]) => `<b>${esc(g)}</b>` + (a ? ` (${esc(tr(a))})` : "")).join(", ")}`;
    case "interaction": return `${t("e_interaction")} ${e.g.map(([a, b]) => `${esc(a)}↔${esc(b)}`).join(", ")}${e.n > e.g.length ? ` (+${e.n - e.g.length})` : ""}`;
    case "similar": return `${t("e_similar")}<ul>${e.m.map((m) => `<li${validSimilar(m) ? "" : ' class="muted"'}>${esc(m[0])} — <b>${m[1]}</b> ${t("shared_genes")}` +
      (m[3] != null ? ` <span class="pval">(p = ${pval(m[3])}${m[4] != null ? `; ${m[4]} ${t("specific")}` : ""})</span>` : "") +
      (m[2] && m[2].length ? `: ${m[2].map(esc).join(", ")}${m[1] > m[2].length ? ", …" : ""}` : "") + "</li>").join("")}</ul>`;
    case "neighbour": return `${t("e_neighbour")} <b>${e.m.map(esc).join(", ")}</b>`;
    case "neighbour2": return `${t("e_neighbour2")} ${e.m.map((m) => Array.isArray(m) ? `${esc(m[0])} (${t("via")} ${esc(m[1])})` : esc(m)).join(", ")}`;
    case "phenotype": return `${t("e_phenotype")} ${e.m.map(([m, c]) => `${esc(m)} (${c})`).join("; ")}`;
    case "similar_drug": return `${t("e_similar_drug")} ${e.g.map(([f2, m2, d]) => `${esc(f2)} → ${esc(m2)} (${t("via")} ${esc(d)})`).join("; ")}`;
    case "adverse_effect": return `${t("e_adverse_effect")} ${e.f.map(esc).join(", ")}`;
    case "clinical": return `${t("e_clinical")} ${esc(tr(e.phase))}` + (e.n ? `, ${e.n} ${t("trials_n")}` : "") +
      (e.year ? `, ${t("since")} ${e.year}` : "") + (e.stopped ? ` — <span class="warning">${e.stopped} ${t("stopped_n")}</span>` : "");
    case "stopped": return `<span class="warning">${t("e_stopped")} ${esc(tr(e.reason))}</span>`;
    case "literature": return `${t("e_literature")} ${e.treats} ${t("treats")}, ${e.causes} ${t("causes")}`;
    default: return "";
  }
}
const EVIDENCE_T = new Set(["neighbour", "similar", "target"]);
const CLINICAL_T = new Set(["clinical", "stopped", "literature"]);
// Shared genes: ≥ 3 and significant (hypergeometric over the genes with a disease, p < 1e-5; P28, P40).
// Must match MAX_P_SHARED in web/build.py.
const validSimilar = (m) => m[1] >= 3 && (m[3] == null || m[3] < 1e-5);
const isEvidence = (e) => (e.t === "similar" ? e.m.some(validSimilar) : EVIDENCE_T.has(e.t));
const pval = (p) => {
  if (p === 0) return "< 1e-300";
  const [m, ex] = p.toExponential(0).split("e");
  return Number(ex) > -3 ? p.toFixed(3) : `${m}×10<sup>${Number(ex)}</sup>`;
};

// Development status (DrugBank)
function statusTag(c) {
  const s = c.status || "unknown";
  const cls = s === "approved" ? "ok" : s === "investigational" ? "warn" : s === "unknown" ? "neutral" : "risk";
  return `<span class="tag ${cls}">${esc((L && L.terms[s]) || s)}</span>${nowTags(c.status_now)}`;
}

// Current regulatory status when the drug is no longer on the market (FDA Orange / Purple Book, EMA; N4)
function nowTags(v) {
  if (!v) return "";
  const out = [];
  if (v.us) out.push(`<span class="tag risk small" title="${esc(t("now_us_help"))}">${t("now_us")}</span>`);
  if (v.eu) out.push(`<span class="tag risk small" title="${esc(t("now_eu_help"))}">${t("now_eu")} ${esc(tr(v.eu))}</span>`);
  return " " + out.join(" ");
}

// «Approved» drugs that no current FDA / EMA label confirms for this disease (display only; docs/DRUG_STATUS.md)
function checkedList(xs, d) {
  if (!xs || !xs.length) return t("none");
  const lc = d.label_check || {}, now = d.status_now || {};
  return xs.map((x) => {
    const c = lc[x];
    const tag = c ? ` <span class="tag ${c[0] === "negated" ? "risk" : "warn"} small" title="${esc(t("lc_" + c[0].replace(/ /g, "_") + "_help") + (c[1] ? " (" + c[1] + ")" : ""))}">${t("lc_" + c[0].replace(/ /g, "_"))}</span>` : "";
    return esc(x) + tag + nowTags(now[x]);
  }).join(", ");
}

// Mechanism: targets with the action (↓ inhibits, ↑ activates)
function mechanism(c) {
  if (!c.mechanism || !c.mechanism.length) return `<span class="muted">${t("none")}</span>`;
  const sym = {"inhibits": "↓", "activates": "↑"};
  return c.mechanism.map(([g, a]) => g === "+" ? `<span class="muted">+${a}</span>` :
    `<span class="gen" title="${esc(a ? tr(a) : t("action_unknown"))}">${esc(g)}${sym[a] ? `<sup>${sym[a]}</sup>` : ""}</span>`).join(" ");
}

// Lines of evidence (facts, no hierarchy; P27)
function evidenceChips(c) {
  const lines = new Set(c.lines || []), xs = [];
  if (lines.has("neighbour")) xs.push(`<span class="chip ev" title="${esc(t("ev_neighbour_help"))}">${t("ev_neighbour")}</span>`);
  if (lines.has("shared_genes")) {
    const s = c.expl.find((e) => e.t === "similar");
    const v = s ? s.m.find(validSimilar) : null;
    const n = v ? v[1] : "";
    const strong = lines.has("specific_genes");
    xs.push(`<span class="chip ${strong ? "ev" : "ev2"}" title="${esc(t(strong ? "ev_specific_help" : "ev_generic_help"))}">${t("ev_shared")} (${n})${strong ? "" : " · " + t("ev_generic")}</span>`);
  }
  if (lines.has("target")) xs.push(`<span class="chip ev2" title="${esc(t("ev_target_help"))}">${t("ev_target")}</span>`);
  return xs.join(" ") || `<span class="muted">${t("none")}</span>`;
}

// Safety and clinical history in this disease
function safetyChips(c) {
  const xs = [], r = c.risk_reasons || {};
  if (c.contra) xs.push(`<span class="chip contra">${t("contraindicated")}</span>`);
  if (c.risk) {
    const reasons = [];
    if (r.contraindicated) reasons.push(t("r_contra").replace("{c}", r.contraindicated).replace("{i}", r.indicated || 0));
    if (r.checkpoint) reasons.push(t("r_checkpoint"));
    xs.push(`<span class="chip risk" title="${esc(reasons.join("; "))}">⚠ ${t("possible_risk")}</span>`);
  }
  const ef = c.expl.find((e) => e.t === "adverse_effect");
  if (ef) xs.push(`<span class="chip neutral" title="${esc(t("e_adverse_effect") + " " + ef.f.join(", "))}">${t("ae_match")}</span>`);
  const cl = c.expl.find((e) => e.t === "clinical");
  if (cl) xs.push(`<span class="chip info" title="${esc(explanation(cl).replace(/<[^>]+>/g, ""))}">${t("tested")}: ${esc(tr(cl.phase))}</span>`);
  if (c.expl.some((e) => e.t === "stopped") || (cl && cl.stopped)) xs.push(`<span class="chip risk">${t("trial_stopped")}</span>`);
  if ((c.pgx || []).length) xs.push(`<span class="chip neutral" title="${esc(t("pgx_help"))}">${t("pgx")}: ${c.pgx.map((x) => esc(x.g)).join(", ")}</span>`);
  if ((c.genetics || []).length) {
    const agree = c.genetics.some((x) => x.agree), oppose = c.genetics.some((x) => !x.agree);
    xs.push(`<span class="chip ${agree && !oppose ? "ev2" : oppose && !agree ? "risk" : "neutral"}" title="${esc(t("gen_help"))}">${t(agree && !oppose ? "gen_agrees" : oppose && !agree ? "gen_opposes" : "gen_mixed")}</span>`);
  }
  return xs.join(" ");
}

function drugLinks(c, disease) {
  const term = encodeURIComponent(`"${c.name}" AND "${searchName(disease)}"`);
  return (c.id && c.id.startsWith("DB") ? ext(`https://go.drugbank.com/drugs/${esc(c.id)}`, "DrugBank") : "") +
    ext(`https://pubmed.ncbi.nlm.nih.gov/?term=${term}`, "PubMed") +
    ext(`https://clinicaltrials.gov/search?cond=${encodeURIComponent(searchName(disease))}&intr=${encodeURIComponent(c.name)}`, "ClinicalTrials.gov") +
    ext(`https://platform.opentargets.org/search?q=${encodeURIComponent(c.name)}`, "Open Targets");
}

function candidateDetail(c, disease) {
  const ev = c.expl.filter(isEvidence), cli = c.expl.filter((e) => CLINICAL_T.has(e.t));
  const ctx = c.expl.filter((e) => !isEvidence(e) && !CLINICAL_T.has(e.t));
  const ul = (xs) => `<ul>${xs.map((e) => `<li>${explanation(e)}</li>`).join("")}</ul>`;
  const r = c.risk_reasons || {};
  return `<div class="detail">
    ${c.neighbour_diseases ? `<p><b>${t("already_indicated")}:</b> ${c.neighbour_diseases.map(esc).join(", ")}</p>` : ""}
    ${ev.length ? `<h4>${t("h_evidence")}</h4>${ul(ev)}` : ""}
    ${cli.length ? `<h4>${t("h_clinical")}</h4>${ul(cli)}` : ""}
    ${c.risk ? `<h4>${t("h_risk")}</h4><ul>${r.contraindicated ? `<li>${t("r_contra").replace("{c}", r.contraindicated).replace("{i}", r.indicated || 0)}</li>` : ""}${r.checkpoint ? `<li>${t("r_checkpoint")}</li>` : ""}</ul>` : ""}
    ${ctx.length ? `<h4>${t("h_context")}</h4>${ul(ctx)}` : ""}
    ${(c.genetics || []).length ? `<h4>${t("h_genetics")}</h4><ul>${c.genetics.map((x) => `<li>${t("gen_line")
      .replace("{g}", `<span class="gen">${esc(x.g)}</span>`).replace("{risk}", t(x.s > 0 ? "gen_higher" : "gen_lower"))
      .replace("{act}", t(x.a < 0 ? "gen_lowers" : "gen_raises"))} → <b>${t(x.agree ? "gen_agrees" : "gen_opposes")}</b>
      <span class="muted">(${esc(x.q)}; GWAS: ${x.st ? ext(`https://platform.opentargets.org/study/${esc(x.st)}`, esc(x.tr || x.st)) : ""}; ${ext(`https://platform.opentargets.org/target/${esc(x.e)}`, "Open Targets")})</span></li>`).join("")}</ul>
      <p class="muted small">${t("gen_note")}</p>` : ""}
    ${(c.pgx || []).length ? `<h4>${t("h_pgx")}</h4><ul>${c.pgx.map((x) => `<li><span class="gen">${esc(x.g)}</span> · CPIC ${esc(x.l)}${x.u ? " · " + ext(x.u, t("guideline")) : ""}</li>`).join("")}</ul>` : ""}
    <p class="links">${drugLinks(c, disease)}</p>
    <p class="feedback small">${t("fb_ask")} ${feedbackLinks(c, disease)}</p></div>`;
}

// Expert feedback: a pre-filled GitHub issue (form .github/ISSUE_TEMPLATE/candidate_feedback.yml) or an email.
// No server involved; the feedback becomes new validation data.
const REPO = "https://github.com/orevid-org/orevid";
function feedbackLinks(c, disease) {
  const page = location.origin + location.pathname + location.hash;
  const q = (o) => Object.entries(o).map(([k, v]) => `${k}=${encodeURIComponent(v)}`).join("&");
  const issue = `${REPO}/issues/new?` + q({template: "candidate_feedback.yml", title: `[Candidate] ${c.name} for ${disease}`,
    disease, drug: c.name, page});
  const mail = `mailto:${CONTACT}?` + q({subject: `OREVID feedback: ${c.name} for ${disease}`, body: page + "\n\n"});
  return ext(issue, t("fb_github")) + ` <span class="muted">${t("fb_or")}</span> ` +
    `<a href="${mail}">${t("fb_email")}</a>`;
}

// Candidate table (the detail row expands)
const TABLES = {};
function candidateTable(id, xs, disease, opts = {}) {
  if (!xs || !xs.length) return `<p class="muted">${opts.empty || t("no_candidates")}</p>`;
  TABLES[id] = {xs, disease};
  const rows = xs.map((c, i) => `<tr class="row" data-i="${i}" tabindex="0" aria-expanded="false"
      data-under="${c.underexplored ? 1 : 0}" data-generic="${hasGeneric(c) ? 1 : 0}">
      <td class="n">${i + 1}</td>
      <td><span class="drug-name">${esc(c.name)}</span>${c.also ? `<div class="muted small">${t("also")}: ${c.also.map(esc).join(", ")}</div>` : ""}
        <div class="muted small" title="${esc(c.atc || "")}">${esc(c.atc ? shortAtc(c.atc) : "")}</div></td>
      <td>${statusTag(c)}</td>
      <td class="small">${genericCell(c)}</td>
      <td>${mechanism(c)}</td>
      <td class="small">${c.indications && c.indications.length ? indications(c.indications) : `<span class="muted">${t("none")}</span>`}</td>
      <td>${opts.extensions ? `<span class="small">${(c.neighbour_diseases || []).map(esc).join("; ")}</span>${c.equivalent ? `<div class="tag equivalent" title="${esc(c.equivalent)}">${t("equivalent")}</div>` : ""}` : evidenceChips(c)}</td>
      <td class="small">${studiedCell(c)}</td>
      <td>${safetyChips(c)}</td>
      <td class="toggle" aria-hidden="true">▸</td></tr>
    <tr class="detail-row" hidden><td></td><td colspan="9">${candidateDetail(c, disease)}</td></tr>`).join("");
  return `<div class="tools"><button class="button" data-csv="${id}">${t("download_csv")}</button>
      <label class="filter"><input type="checkbox" data-filter="under" data-table="${id}"> ${t("f_underexplored")}</label>
      <label class="filter"><input type="checkbox" data-filter="generic" data-table="${id}"> ${t("f_generic")}</label></div>
    <div class="table-scroll"><table class="table" id="${id}"><thead><tr>
      <th>#</th><th>${t("c_drug")}</th><th>${t("c_status")}</th><th title="${esc(t("generic_help"))}">${t("c_generic")}</th><th>${t("c_mechanism")}</th><th>${t("c_indications")}</th>
      <th>${opts.extensions ? t("c_neighbour") : t("c_evidence")}</th><th title="${esc(t("studied_help"))}">${t("c_studied")}</th><th>${t("c_safety")}</th><th></th></tr></thead>
      <tbody>${rows}</tbody></table></div>`;
}
// Practical information (generic / biosimilar status, how much the pair has been studied; display only)
const hasGeneric = (c) => { const g = c.generic || {}; return !!(g.generic_us || g.biosimilar_us || g.generic_eu || g.biosimilar_eu); };
function genericCell(c) {
  const g = c.generic || {};
  const tags = [["generic_us", "g_generic_us"], ["biosimilar_us", "g_biosimilar_us"], ["generic_eu", "g_generic_eu"],
    ["biosimilar_eu", "g_biosimilar_eu"]].filter(([k]) => g[k]).map(([, l]) => `<span class="tag generic">${t(l)}</span>`);
  const routes = g.routes ? `<div class="muted">${esc(g.routes)}</div>` : "";
  if (tags.length) return tags.join(" ") + routes;
  if (g.exclusivity_until) return `<span class="muted">${t("protected_until")} ${esc(g.exclusivity_until.slice(0, 4))}</span>` + routes;
  return `<span class="muted">${t("none")}</span>` + routes;
}
function studiedCell(c) {
  if (c.n_trials == null) return "";
  const s = `${t("trials_short")}: ${c.n_trials} · ${t("papers_short")}: ${c.n_papers}`;
  return c.underexplored ? `${s}<div class="tag underexplored">${t("underexplored")}</div>` : s;
}
// Current indications: the first 3 and how many more
function indications(xs) {
  const more = xs.filter((x) => x.startsWith("+")).map((x) => Number(x.slice(1)))[0] || 0;
  const names = xs.filter((x) => !x.startsWith("+"));
  const extra = Math.max(names.length - 3, 0) + more;
  return names.slice(0, 3).map(esc).join("; ") + (extra ? ` <span class="muted">+${extra}</span>` : "");
}
// Some DrugBank ATC classes are many chained classes: the first one is shown.
function shortAtc(s) {
  return s.length > 48 ? s.slice(0, s.lastIndexOf(" ", 46)) + " …" : s;
}

function csv(id) {
  // Always in English (to share and process), whatever the website language; the data values are English
  const {xs, disease} = TABLES[id];
  const q = (v) => `"${String(v ?? "").replace(/"/g, '""')}"`;
  const LINES = {neighbour: "neighbour_disease", shared_genes: "shared_genes", specific_genes: "specific_shared_genes",
    target: "target_on_disease_gene"};
  const header = ["rank", "drug", "drugbank_id", "status", "atc_class", "targets", "current_indications", "evidence_lines",
    "shared_genes_disease", "shared_genes_n", "shared_genes_p", "shared_genes_specific_n", "shared_genes_examples",
    "neighbour_disease", "risk_flag", "risk_reason", "tested_in_this_disease", "generic_us", "biosimilar_us",
    "generic_eu_centralised", "biosimilar_eu_centralised", "us_protection_until", "routes", "trials_for_pair", "papers_for_pair",
    "underexplored"];
  const lines = xs.map((c, i) => {
    const s = c.expl.find((e) => e.t === "similar"), v = c.expl.find((e) => e.t === "neighbour"), cl = c.expl.find((e) => e.t === "clinical");
    const m = s ? (s.m.find(validSimilar) || s.m[0]) : null;
    const r = c.risk_reasons || {};
    return [i + 1, c.name, c.id, c.status || "", c.atc, (c.mechanism || []).filter(([g]) => g !== "+").map(([g, a]) => g + (a ? `(${a})` : "")).join(" "),
      (c.indications || []).join("; "), (c.lines || []).map((x) => LINES[x] || x).join(" "),
      m ? m[0] : "", m ? m[1] : "", m && m[3] != null ? m[3] : "", m && m[4] != null ? m[4] : "", m ? (m[2] || []).join(" ") : "",
      v ? v.m.join("; ") : (c.neighbour_diseases || []).join("; "),
      c.risk ? "yes" : "", [r.contraindicated ? `contraindicated in ${r.contraindicated} related diseases, indicated in ${r.indicated || 0}` : "",
        r.checkpoint ? "immune checkpoint blocker in immune-mediated disease" : ""].filter(Boolean).join("; "),
      cl ? cl.phase || "" : "", ...["generic_us", "biosimilar_us", "generic_eu", "biosimilar_eu"].map((k) => (c.generic || {})[k] ? "yes" : ""),
      (c.generic || {}).exclusivity_until || "", (c.generic || {}).routes || "", c.n_trials ?? "", c.n_papers ?? "",
      c.underexplored ? "yes" : ""].map(q).join(",");
  });
  const blob = new Blob(["﻿" + [header.join(","), ...lines].join("\n")], {type: "text/csv;charset=utf-8"});
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob);
  a.download = `orevid_${disease.replace(/[^a-z0-9]+/gi, "_")}_${id}.csv`;
  a.click();
  setTimeout(() => URL.revokeObjectURL(a.href), 1000);
}

document.addEventListener("click", (e) => {
  const b = e.target.closest("[data-csv]");
  if (b) { csv(b.dataset.csv); return; }
  const row = e.target.closest("tr.row");
  if (row && !e.target.closest("a")) toggleRow(row);
});
document.addEventListener("change", (e) => {
  const f = e.target.closest("[data-filter]");
  if (!f) return;
  const table = document.getElementById(f.dataset.table);
  const on = [...document.querySelectorAll(`[data-filter][data-table="${f.dataset.table}"]`)].filter((x) => x.checked).map((x) => x.dataset.filter);
  for (const row of table.querySelectorAll("tr.row")) {
    const show = on.every((k) => row.dataset[k] === "1");
    row.hidden = !show;
    if (!show) row.nextElementSibling.hidden = true;
  }
});
document.addEventListener("keydown", (e) => {
  const row = e.target.closest && e.target.closest("tr.row");
  if (row && (e.key === "Enter" || e.key === " ")) { toggleRow(row); e.preventDefault(); }
});
function toggleRow(row) {
  const detail = row.nextElementSibling, opening = detail.hidden;
  detail.hidden = !opening;
  row.setAttribute("aria-expanded", String(opening));
  row.querySelector(".toggle").textContent = opening ? "▾" : "▸";
}

function targets(d) {
  if (!d || !d.new.length) return `<p class="muted">${t("none")}</p>`;
  const sorted = [...d.new].sort((a, b) => b.s - a.s);
  const direct = sorted.filter((x) => x.direct), others = sorted.filter((x) => !x.direct);
  const header = `<thead><tr><th>${t("d_gene")}</th><th>${t("evidence")}</th><th>${t("direction")}</th>
      <th>${t("modalities")}</th><th>${t("safety")}</th><th>${t("drugs_acting")}</th><th></th></tr></thead>`;
  const table = (xs) => `<div class="table-scroll"><table class="table targets">${header}<tbody>${targetRows(xs)}</tbody></table></div>`;
  return (d.exploited.length ? `<p class="muted">${t("exploited")} ${d.exploited.map(esc).join(", ")}</p>` : "") +
    (direct.length ? table(direct) : `<p class="muted">${t("no_direct_target")}</p>`) +
    (others.length ? `<details class="others"><summary>${t("indirect_targets")} (${others.length})</summary>${table(others)}</details>` : "");
}
function targetRows(xs) {
  return xs.map((x) => {
    const flags = [];
    if (x.ess === 1 || (x.dep != null && x.dep > 0.5)) flags.push(`<span class="chip risk" title="${esc(t("essential"))}">DepMap</span>`);
    if (x.loeuf != null && x.loeuf < 0.35) flags.push(`<span class="chip risk" title="${esc(t("loeuf_help"))}">LOEUF ${x.loeuf.toFixed(2)}</span>`);
    if (x.tissue != null && x.tissue > 0.5) flags.push(`<span class="chip neutral" title="${esc(t("tissue"))}">${t("tissue_short")}</span>`);
    return `<tr><td><b>${esc(x.gene)}</b></td>
      <td class="small">${esc(evTarget(x.ev))}</td><td class="small">${esc(trTarget(x.dir))}</td><td class="small">${esc(trTarget(x.mod))}</td>
      <td>${flags.join(" ")}${x.safety ? `<div class="muted small">${esc(x.safety)}</div>` : ""}</td>
      <td class="small">${x.drugs ? esc(trTarget(x.drugs)) : `<span class="muted">${t("none")}</span>`}</td>
      <td class="links">${x.uniprot ? ext(`https://www.uniprot.org/uniprotkb/${esc(x.uniprot)}`, "UniProt") +
        ext(`https://www.ebi.ac.uk/chembl/g/#search_results/targets/query=${esc(x.uniprot)}`, "ChEMBL") : ""}${ext(`https://platform.opentargets.org/search?q=${encodeURIComponent(x.gene)}`, "Open Targets")}</td></tr>`;
  }).join("");
}

// MONDO: PrimeKG groups some terms (ids joined with «_»); the first one is shown.
function mondo(id, name) {
  const ids = String(id || "").split("_").filter(Boolean);
  if (!ids.length || !/^\d+$/.test(ids[0])) return "";
  if (ids.length > 1) {  // PrimeKG groups terms: the first id is not necessarily the one of the name
    return ext(`https://monarchinitiative.org/search?q=${encodeURIComponent(name)}`, "MONDO") +
      ` <span class="muted">(${ids.length} ${t("grouped")})</span>`;
  }
  const m = `MONDO:${ids[0].padStart(7, "0")}`;
  return ext(`https://monarchinitiative.org/${m}`, m);
}
// Name to search PubMed and ClinicalTrials.gov: without qualifiers such as «(disease)»
const searchName = (s) => String(s).replace(/\s*\((disease|disorder)\)/gi, "").trim();

// Loading in steps: the user sees what is happening (done steps in grey, the current one with animated dots)
function steps(done, current) {
  document.getElementById("view").innerHTML = `<div class="loading" role="status" aria-live="polite">` +
    done.map((x) => `<div class="step done">✓ ${x}</div>`).join("") +
    (current ? `<div class="step current">${current}<span class="dots"></span></div>` : "") + `</div>`;
}
// let the step paint before building the tables (setTimeout: requestAnimationFrame does not run in a hidden tab)
const yieldFrame = () => new Promise((r) => setTimeout(r, 20));

async function showDisease(id) {
  document.querySelector(".search").classList.add("compact");
  const done = [];
  if (!INDEX) {
    steps(done, t("step_index"));
    await (INDEX_PROMISE ||= loadIndex());
    done.push(t("step_index_done").replace("{n}", INDEX.length.toLocaleString(lang)));
  }
  const view = document.getElementById("view");
  const info = INDEX && INDEX.find((m) => m.id === id);
  steps(done, t("step_page").replace("{m}", esc(info ? info.name : id)));
  const d = await diseasePage(id);
  if (d) {
    done.push(t("step_page_done").replace("{c}", d.candidates.length).replace("{e}", d.extensions.length));
    steps(done, t("step_tables"));
    await yieldFrame();
  }
  if (!d) {
    view.innerHTML = `<div class="disease-page"><h2>${esc(info ? info.name : id)}</h2><p class="muted">${t("no_page")}</p></div>`;
    return;
  }
  document.title = `${d.name} — OREVID`;
  const name = encodeURIComponent(searchName(d.name));
  const trials = Object.entries(d.trials || {}).map(([phase, a]) => {
    const xs = Array.isArray(a) ? a : a.drugs, n = Array.isArray(a) ? a.length : a.n;
    return `<b>${esc(tr(phase))}</b> (${n}): ${list(xs.slice(0, 15))}${n > 15 ? ` … <span class="muted">+${n - 15}</span>` : ""}`;
  }).join("<br>");
  const gs = d.gene_sources || {};
  const sources = ["genetic", "primekg", "animal", "expression"].filter((k) => gs[k]).map((k) => `${t("gs_" + k)} ${gs[k]}`).join(", ");
  const genes = d.genes || [];
  const nav = (xs) => xs.map(([id, n, page]) => page ? `<a href="#/d/${id}">${esc(n)}</a>` : `<span class="muted">${esc(n)}</span>`).join(", ");
  const broader = d.broader || [], narrower = d.narrower || [];
  const hierarchy = broader.length || narrower.length ? `<div class="hierarchy">
      ${broader.length ? `<p><b>${t("broader")}:</b> ${nav(broader)}</p>` : ""}
      ${narrower.length ? `<p><b>${t("narrower")}:</b> ${nav(narrower.slice(0, 12))}${narrower.length > 12 ? ` <details class="more"><summary>${t("more_forms").replace("{n}", narrower.length - 12)}</summary>${nav(narrower.slice(12))}</details>` : ""}</p>` : ""}
    </div>` : "";
  const broad = (d.n_descendants || 0) >= 5 || d.n_genes >= 500;
  view.innerHTML = `<div class="disease-page">
    <h2>${esc(d.name)}</h2>
    <div class="meta">${mondo(d.mondo, d.name)} · ${d.n_genes} ${t("genes")}${sources ? ` <span class="muted">(${sources})</span>` : ""} ·
      ${ext(`https://platform.opentargets.org/search?q=${name}`, "Open Targets")} ·
      ${ext(`https://clinicaltrials.gov/search?cond=${name}`, "ClinicalTrials.gov")}${d.prevalence ? ` · <span title="${esc(t("prev_help"))}">${t("prevalence")}: ${esc(d.prevalence.c)}</span> (${ext(`https://www.orpha.net/en/disease/detail/${esc(d.prevalence.o)}`, "Orphanet")})` : ""}</div>
    ${hierarchy}
    ${broad ? `<p class="warning broad">${t("broad_term")}</p>` : ""}
    ${(d.markers || []).length ? `<p class="markers"><b>${t("markers")}:</b> ${d.markers.map(([m, st]) => `${esc(m)} ${t("marker_" + st)}`).join(", ")}</p>` : ""}
    ${genes.length ? `<details class="genes"><summary>${t("genes_list")} (${d.n_genes})</summary><p>${genes.map((g) => `<span class="gen">${esc(g)}</span>`).join(" ")}${d.n_genes > genes.length ? " …" : ""}</p></details>` : ""}
    <div class="card known"><h3>${t("known")}</h3><dl>
      <dt>${t("approved")}</dt><dd>${checkedList(d.approved, d)}</dd>
      ${d.approved_ot && d.approved_ot.length ? `<dt>${t("approved_ot")}</dt><dd>${checkedList(d.approved_ot, d)}</dd>` : ""}
      ${d.approved_unconfirmed && d.approved_unconfirmed.length ? `<dt>${t("approved_unconfirmed")}</dt><dd class="muted">${list(d.approved_unconfirmed)} <div class="small">${t("approved_unconfirmed_help")}</div></dd>` : ""}
      ${d.offlabel.length ? `<dt>${t("offlabel")}</dt><dd>${list(d.offlabel)}</dd>` : ""}
      ${trials ? `<dt>${t("trials")}</dt><dd>${trials}</dd>` : ""}
      ${d.stopped.length ? `<dt>${t("stopped")}</dt><dd>${list(d.stopped)}</dd>` : ""}
      ${d.contra.length ? `<dt>${t("contra")}</dt><dd>${list(d.contra)}</dd>` : ""}
    </dl></div>
    ${d.extensions.length ? `<div class="card"><h3>${t("extensions_title")}</h3><p class="help">${t("extensions_help")}</p>
      ${candidateTable("extensions", d.extensions.map((c) => (typeof c === "object" && !Array.isArray(c) ? c : null)).filter(Boolean), d.name, {extensions: true})}</div>` : ""}
    <div class="card"><h3>${t("candidates")}</h3><p class="help">${t("candidates_help")}</p>
      ${d.little_info ? `<p class="warning">${t("little_info")}</p>` : ""}
      ${candidateTable("candidates", d.candidates, d.name)}
      ${d.hidden ? `<p class="muted">${t("hidden").replace("{n}", d.hidden)}</p>` : ""}
      ${(d.excluded_marker || []).length ? `<p class="muted">${t("excluded_marker")} ${d.excluded_marker.map(([n, m]) => `${esc(n)} (${esc(m)})`).join(", ")}</p>` : ""}</div>
    <div class="card"><h3>${t("targets")}</h3><p class="help">${t("targets_help")}</p>${targets(d.targets)}</div>
    ${d.no_use.length || d.excluded.length ? `<details class="card exploratory"><summary><h3>${t("no_use")}</h3></summary>
      <p class="help">${t("no_use_help")}</p>
      ${candidateTable("no_use", d.no_use, d.name, {empty: t("none")})}
      ${d.excluded.length ? `<p class="muted">${t("excluded")} ${d.excluded.map(esc).join(", ")}</p>` : ""}</details>` : ""}
  </div>`;
}

// --- method and legal notice --------------------------------------------------------------
// Owner and contact address shown in the legal notice
const OWNER = "Joaquim Quadrada Jané", CONTACT = "contact@orevid.org";

function showLegal() {
  document.querySelector(".search").classList.add("compact");
  document.getElementById("view").innerHTML = L.legal.replaceAll("{owner}", esc(OWNER))
    .replaceAll("{contact}", `<a href="mailto:${esc(CONTACT)}">${esc(CONTACT)}</a>`);
  document.title = `${t("nav_legal")} — OREVID`;
}

function showAbout() {
  document.querySelector(".search").classList.add("compact");
  document.getElementById("view").innerHTML = L.about;
  document.title = `${t("nav_about")} — OREVID`;
}

function showMethod() {
  document.querySelector(".search").classList.add("compact");
  document.getElementById("view").innerHTML = L.method;
  document.title = `${t("nav_method")} — OREVID`;
}

// --- navigation and language --------------------------------------------------------------
function applyLang() {
  document.documentElement.lang = lang;
  document.querySelectorAll("[data-t]").forEach((el) => { el.innerHTML = t(el.dataset.t); });
  q.setAttribute("aria-label", t("search_label"));
  document.querySelectorAll("[data-lang]").forEach((b) => b.setAttribute("aria-pressed", String(b.dataset.lang === lang)));
}
async function switchLang(e) {
  const b = e.target.closest("[data-lang]");
  if (!b || b.dataset.lang === lang) return;
  lang = b.dataset.lang;
  try { localStorage.setItem("lang", lang); } catch (e) { /* no storage */ }
  try { L = await loadLocale(lang); } catch (e) { showError(e); return; }
  applyLang();
  route();
}
document.getElementById("langs").addEventListener("click", switchLang);
document.querySelector(".notice-langs").addEventListener("click", switchLang);

// --- first-visit notice -------------------------------------------------------------------
// A modal on the first visit, whatever page it starts on (search engines and shared links land on any page): the
// candidates are hypotheses, not treatments. Only «I understand» dismisses it and is remembered in this browser;
// Escape does not, and the legal-notice link closes it for this visit only. Without storage it shows on every visit.
const NOTICE_KEY = "notice_ok_v1";
function showNotice() {
  const dlg = document.getElementById("notice");
  let seen = false;
  try { seen = localStorage.getItem(NOTICE_KEY) === "1"; } catch (e) { /* no storage */ }
  if (seen || typeof dlg.showModal !== "function") return;
  dlg.addEventListener("cancel", (e) => e.preventDefault());
  document.getElementById("notice-ok").addEventListener("click", () => {
    try { localStorage.setItem(NOTICE_KEY, "1"); } catch (e) { /* no storage */ }
    dlg.close();
  });
  document.getElementById("notice-legal").addEventListener("click", () => dlg.close());
  dlg.showModal();
}

// Never stay stuck on «Loading…»: say what failed and how to try again
function showError(e) {
  console.error(e);
  document.getElementById("view").innerHTML = `<div class="card"><p class="warning">${t("load_error")}</p>
    <p class="muted small">${esc(e && e.message)}</p><button class="button" onclick="location.reload()">${t("retry")}</button></div>`;
}

async function route() {
  const h = location.hash;
  const m = h.match(/^#\/d\/(\d+)/);
  sug.style.display = "none";
  try {
    if (m) { await showDisease(Number(m[1])); window.scrollTo(0, 0); }
    else if (h === "#/method") showMethod();
    else if (h === "#/about") showAbout();
    else if (h === "#/legal") showLegal();
    else { document.querySelector(".search").classList.remove("compact"); document.getElementById("view").innerHTML = ""; document.title = "OREVID"; }
  } catch (e) {
    showError(e);
  }
}

async function start() {
  INDEX_PROMISE = loadIndex();
  INDEX_PROMISE.catch((e) => console.error(e));
  try { L = await loadLocale(lang); } catch (e) { showError(e); return; }
  applyLang();
  showNotice();
  window.addEventListener("hashchange", route);
  route();
}
start();
