"""Research report of a splicing study (Workflow step 13, `splicing.py report`).

Writes <study-dir>/report.md in Portuguese with the sections Estudo, Decisões, Coorte e rejeitados,
Threshold, Priorizadas, Skills utilizadas, Limitações and Citações. Everything comes from files of the
study directory: study.yaml, decisions.yaml, manifest.json, processed/* and, when present,
raw/atlas_links.md and raw/interpret_*.txt (optional: they are included only if the researcher produced
them). The skills listed are the ones recorded in manifest.json, so the report names what was actually
used. Nothing here is specific to one gene. Research use only: the report describes the predicted
molecular impact of variants and never a diagnosis.

Library: ReportError, render_report, write_report, summarize_interpretation.
"""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import polars as pl

from core import clinical as core_clinical
from core import locus as core_locus
from core import prioritize as core_prioritize
from core import threshold as core_threshold
from core.threshold import RocResult
from study import load_decisions, load_study

REPORT_FILE = "report.md"
MANIFEST_FILE = "manifest.json"
ATLAS_LINKS_FILE = "raw/atlas_links.md"
INTERPRET_GLOB = "interpret_*.txt"
SATURATION_FILE = "processed/saturation_scored.parquet"
DEFAULT_TOP_N = 20
BULK_REASON = "outside_region"
MAX_SUMMARY_LINES = 14
PRIMARY_CITATION = ("Avsec, Ž. et al. (2026). Advancing regulatory variant effect prediction with AlphaGenome. "
                    "Nature. (AlphaGenome model and Atlas predictions used throughout this study.)")
_DEEPMIND = "https://github.com/google-deepmind/science-skills"

# Citation of each dependency skill (key = `skill` field of manifest.json steps).
SKILL_CITATIONS: dict[str, str] = {
    "alphagenome-variant-impact-score": (
        f"DeepMind science skills, `alphagenome-variant-impact-score` (GENCODE gene structure, AlphaGenome Atlas "
        f"variant impact scores and splice scorers). {_DEEPMIND}. Underlying model: Avsec et al. (2026)."),
    "alphagenome-single-variant-analysis": (
        f"DeepMind science skills, `alphagenome-single-variant-analysis` (per-variant AlphaGenome predictions and "
        f"splicing quantification). {_DEEPMIND}. Underlying model: Avsec et al. (2026)."),
    "alphagenome-atlas-website-links": (
        f"DeepMind science skills, `alphagenome-atlas-website-links` (Atlas deep links). {_DEEPMIND}."),
    "clinvar-database": (
        f"DeepMind science skills, `clinvar-database`: NCBI ClinVar (Landrum et al., Nucleic Acids Res.). {_DEEPMIND}; "
        "https://www.ncbi.nlm.nih.gov/clinvar/."),
    "gnomad-database": (
        f"DeepMind science skills, `gnomad-database`: gnomAD (Karczewski et al., Nature 2020; Chen et al., Nature 2024). "
        f"{_DEEPMIND}; https://gnomad.broadinstitute.org."),
    "dbsnp-database": (
        f"DeepMind science skills, `dbsnp-database`: NCBI dbSNP (Sherry et al., Nucleic Acids Res. 2001). {_DEEPMIND}; "
        "https://www.ncbi.nlm.nih.gov/snp/."),
    "ensembl-database": (
        f"DeepMind science skills, `ensembl-database`: Ensembl REST API (reference sequence, Variant Effect Predictor). "
        f"{_DEEPMIND}; https://rest.ensembl.org."),
    "alphagenome": (
        "K-Dense scientific agent skills, `alphagenome` (AlphaGenome API: per-track scores, indels). "
        "https://github.com/K-Dense-AI/scientific-agent-skills. Underlying model: Avsec et al. (2026)."),
}
LOCAL_SKILL = "splicing-analysis"


class ReportError(Exception):
    """A study file the report needs is missing or unreadable; the CLI prints it and exits with 1."""


def _read_json(path: Path, step: str) -> Any:
    if not path.is_file():
        raise ReportError(f"{path.name} not found in {path.parent}; {step}.")
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ReportError(f"{path}: invalid JSON ({exc}).") from None


def _need(study_dir: Path, name: str, step: str) -> Path:
    path = study_dir / name
    if not path.is_file():
        raise ReportError(f"{name} not found in {study_dir}; {step}.")
    return path


def _cell(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.4g}"
    return str(value).replace("|", "\\|").replace("\n", " ")


def _table(headers: Sequence[str], rows: Sequence[Sequence[Any]]) -> list[str]:
    out = ["| " + " | ".join(headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
    out += ["| " + " | ".join(_cell(c) for c in row) + " |" for row in rows]
    return out


def _skills_used(manifest: Mapping[str, Any]) -> list[tuple[str, list[str]]]:
    """Skills of the manifest steps, in first-use order, with the step names that used them."""
    used: dict[str, list[str]] = {}
    for step in manifest.get("steps", []):
        used.setdefault(step["skill"], []).append(step["step"])
    return list(used.items())


def _fmt(value: float | None, digits: int = 4) -> str:
    return "n/d" if value is None else f"{value:.{digits}f}"


def summarize_interpretation(text: str) -> list[str]:
    """Key lines of an `alphagenome-single-variant-analysis` splicing text: the header, the splice-site usage
    lines, the junction flags (`->`) and the table of the most changed junctions. Falls back to the first lines."""
    lines = [ln.rstrip() for ln in text.splitlines()]
    out: list[str] = []
    for ln in lines:
        if ln.strip():
            out.append(ln)
            break
    out += [ln for ln in lines if ln.startswith("Position ") and "Delta=" in ln]
    out += [ln for ln in lines if ln.startswith("->")]
    for i, ln in enumerate(lines):
        if ln.startswith("Top ") and "Most Changed Junctions" in ln:
            block = [ln]
            for nxt in lines[i + 1:]:
                if not nxt.strip():
                    break
                block.append(nxt)
            out += block
            break
    if len(out) <= 1:
        out = [ln for ln in lines if ln.strip()][:MAX_SUMMARY_LINES]
    return out[:MAX_SUMMARY_LINES]


# --- sections -------------------------------------------------------------------------------


def _section_study(study: Any, locus: core_locus.SpliceLocus, manifest: Mapping[str, Any],
                   rules_name: str) -> list[str]:
    strand = "+" if locus.strand == 1 else "-"
    mane = ", MANE Select" if locus.is_mane_select else ""
    introns = ", ".join(f"{i.number} ({locus.chromosome}:{i.start:,}-{i.end:,}, {i.length:,} pb)"
                        for i in locus.introns_of_interest)
    sites = ", ".join(f"{s.label} ({s.position:,})" for s in locus.splice_sites)
    lines = [
        "## Estudo", "",
        f"- Gene: **{locus.gene_symbol}** (`{locus.gene_id}`), {study.assembly}, fita {strand}.",
        f"- Transcrito: `{locus.transcript_id}` ({locus.transcript_type}{mane}), {len(locus.exons)} éxons.",
        f"- Íntron(s) de interesse: {introns}.",
        f"- Sítios de splicing: {sites}.",
    ]
    sat = locus.saturation
    if sat is not None:
        rule = ("`saturation_region` do study.yaml" if sat.rule == core_locus.SATURATION_EXPLICIT else
                f"±{sat.flank_bp} pb em torno de cada doador e aceitador dos íntrons de interesse, janelas "
                "sobrepostas unidas")
        spans = ", ".join(f"{locus.chromosome}:{w.start:,}-{w.end:,}" for w in sat.windows)
        lines.append(f"- Varredura de saturação: {len(sat.windows)} janela(s), {sat.total_bp:,} pb (cerca de "
                     f"{3 * sat.total_bp:,} SNVs; regra: {rule}): {spans}. As variantes priorizadas vêm só dessas "
                     "janelas.")
    if locus.isoforms:
        seg = locus.isoform_segment
        lines.append(f"- Isoformas comparadas: referência `{locus.isoforms['reference']}` e alternativa "
                     f"`{locus.isoforms['alternative']}`"
                     + (f"; segmento específico {locus.chromosome}:{seg.start:,}-{seg.end:,}" if seg else "") + ".")
    for m in locus.manual_elements:
        lines.append(f"- Elemento manual **{m.label}** ({locus.chromosome}:{m.start:,}-{m.end:,}): {m.citation}")
    zones = ", ".join(locus.zone_labels())
    lines += [
        f"- Janela de visualização: {locus.chromosome}:{locus.view.start:,}-{locus.view.end:,}; zonas de distância ao "
        f"sítio mais próximo: {zones}.",
        f"- Regras de classe da coorte: `{rules_name}`; peso do termo de junção no composto: "
        f"{(study.data.get('composite') or {}).get('junction_weight', 0.2):g}.",
    ]
    origin = manifest.get("origin", {})
    lines += ["", f"Origem dos dados brutos: `{origin.get('kind', 'desconhecida')}`"
              + (f". {origin['note']}" if origin.get("note") else "."), ""]
    return lines


def _section_decisions(study_dir: Path) -> list[str]:
    decisions = load_decisions(study_dir)
    lines = ["## Decisões", ""]
    if not decisions:
        return lines + ["Nenhuma decisão registrada em `decisions.yaml`.", ""]
    lines += ["Respostas do pesquisador registradas em `decisions.yaml` (passo do workflow, chave, opções oferecidas, "
              "recomendada, valor escolhido); \"—\" quando a pergunta foi registrada sem opções:", ""]

    def show(value: Any) -> str:
        if isinstance(value, list):
            return ", ".join(map(str, value))
        if isinstance(value, dict):
            return ", ".join(f"{k}: {v}" for k, v in value.items())
        return f"`{value}`"

    rows = [(e.workflow_step, f"`{key}`", "; ".join(f"`{o}`" for o in e.options) if e.options else "—",
             f"`{e.recommended}`" if e.recommended is not None else "—", show(e.value), e.answered_at[:10])
            for key, e in sorted(decisions.items(), key=lambda kv: (kv[1].workflow_step, kv[0]))]
    return lines + _table(["Passo", "Decisão", "Opções", "Recomendada", "Valor", "Data"], rows) + [""]


def _section_cohort(study_dir: Path, rules: Sequence[core_clinical.ClassRule]) -> list[str]:
    cohort = pl.read_parquet(_need(study_dir, core_clinical.COHORT_FILE, "run splicing.py cohort (Workflow step 7)"))
    rejected = pl.read_csv(_need(study_dir, core_clinical.REJECTED_FILE, "run splicing.py cohort (Workflow step 7)"),
                           separator="\t", infer_schema_length=0)
    counts = dict(cohort.group_by("class_key").len().iter_rows())
    lines = ["## Coorte e rejeitados", "",
             f"A coorte tem **{cohort.height:,}** variantes (região do íntron de interesse e éxons vizinhos). "
             "Variantes por classe (papel `case` e `control` entram no ROC; `excluded` fica na coorte e nas "
             "figuras, fora do ROC):", ""]
    lines += _table(["Classe", "Papel", "Variantes"],
                    [(f"{r.label} (`{r.key}`)", r.role, counts.get(r.key, 0)) for r in rules])
    comp = dict(cohort.group_by("compartment").len().iter_rows())
    kinds = dict(cohort.group_by("variant_type").len().iter_rows())
    lines += ["", "Compartimentos: " + ", ".join(f"{k or 'fora do transcrito'} {n}" for k, n in sorted(
        comp.items(), key=lambda kv: str(kv[0]))) + ". Tipos: "
              + ", ".join(f"{k} {n}" for k, n in sorted(kinds.items(), key=lambda kv: str(kv[0]))) + ".", ""]
    if rejected.is_empty():
        return lines + ["Nenhum registro foi rejeitado (`processed/rejected.tsv` vazio).", ""]
    bulk = rejected.filter(pl.col("reason") == BULK_REASON)
    rest = rejected.filter(pl.col("reason") != BULK_REASON)
    lines += [f"Registros rejeitados (`processed/rejected.tsv`, {rejected.height:,} no total):", ""]
    if bulk.height:
        by_source = ", ".join(f"{s} {n:,}" for s, n in bulk.group_by("source").len().sort("source").iter_rows())
        lines += [f"- **{bulk.height:,}** registros fora da região da coorte (`{BULK_REASON}`: {by_source}), "
                  "agregados nesta linha por serem o volume esperado de uma consulta regional.", ""]
    if rest.height:
        grouped = rest.group_by("source", "reason").len().sort(["source", "reason"])
        lines += _table(["Fonte", "Motivo", "Registros"], list(grouped.iter_rows()))
        lines += [""]
        meaningful = rest.sort("source", "reason", "record_id").group_by(["source", "reason"], maintain_order=True).head(3)
        if meaningful.height:
            lines += ["Exemplos de rejeitados (até 3 por fonte e motivo):", ""]
            lines += _table(["Fonte", "Registro", "Título", "Motivo", "Detalhe"],
                            [(r["source"], r["record_id"], r["title"], r["reason"], r["detail"]) for r in meaningful.iter_rows(named=True)])
            lines += [""]
    return lines


def _section_threshold(roc: RocResult, tau_name: str) -> list[str]:
    y = roc.youden
    lines = ["## Threshold", ""]
    if roc.calibrated:
        lines += [f"O limiar foi **calibrado** nesta coorte (casos contra controles, composto `{roc.score_column}`; "
                  "uma chamada é `composto >= τ`). O τ vale para este estudo e este peso de junção; não é um "
                  "limiar universal.", ""]
        lines += _table(["Métrica", "Valor"], [
            ("Casos / controles com escore", f"{roc.n_case} / {roc.n_control}"),
            ("AUROC", _fmt(roc.auroc, 6)), ("AUPRC", _fmt(roc.auprc, 6)),
            (f"{tau_name} (τ)", _fmt(y.threshold, 6)),
            ("Sensibilidade no τ", _fmt(y.sensitivity, 4)), ("Especificidade no τ", _fmt(y.specificity, 4)),
        ])
    else:
        lines += [f"**O limiar NÃO foi calibrado.** A coorte tem {roc.n_case} caso(s) e {roc.n_control} controle(s) com "
                  f"escore, e o ROC precisa dos dois; o pesquisador escolheu o τ fixo {y.threshold:g} "
                  f"(decisão `{core_threshold.NO_CASES_KEY}`). AUROC e AUPRC não existem para este estudo.", ""]
        lines += _table(["Métrica", "Valor"], [
            ("Casos / controles com escore", f"{roc.n_case} / {roc.n_control}"), ("τ fixo", f"{y.threshold:g}")])
    lines += ["", "Limiares avaliados nos mesmos casos e controles:", ""]
    points = [(f"{y.label} (estudo)", y), *((p.label, p) for p in roc.fixed)]
    lines += _table(["Limiar", "τ", "Sensibilidade", "Especificidade", "TP", "FP", "FN", "TN"],
                    [(name, f"{p.threshold:.6g}", _fmt(p.sensitivity, 4), _fmt(p.specificity, 4), p.tp, p.fp, p.fn, p.tn)
                     for name, p in points])
    lines += [""]
    return lines


def _section_prioritized(study_dir: Path, roc: RocResult, top_n: int) -> list[str]:
    frame = core_prioritize.load_prioritized(study_dir)
    counts = core_prioritize.summary_counts(frame, roc)
    name = "Youden" if roc.calibrated else "τ fixo (não calibrado)"
    lines = ["## Priorizadas", "",
             f"{counts['n_candidates']:,} SNVs da varredura de saturação fora do ClinVar foram comparados ao limiar do "
             f"estudo: **{counts['n_hypothetical']:,} hipotéticas** (sem registro no ClinVar nem no gnomAD consultados, "
             "contando também os registros que a coorte rejeitou) e "
             f"**{counts['n_observed']:,} observadas** no gnomAD. Acima do τ ({name}, {roc.tau:.4f}): "
             f"{counts['hypothetical_ge_youden']} hipotéticas e {counts['observed_ge_youden']} observadas; "
             f"com composto >= 1,0: {counts['hypothetical_ge_1']} hipotéticas e {counts['observed_ge_1']} observadas. "
             "Observada não significa benigna e hipotética não significa patogênica: são previsões de impacto "
             "no splicing para escolher experimentos.", "",
             f"Principais candidatas (top {min(top_n, frame.height)} por composto; `processed/prioritized.tsv`):", ""]
    rows = [(i, f"`{r['variant_id']}`", r["hgvs_c"], r["score"], r["tier"], r["origin"], r["top_modality"])
            for i, r in enumerate(frame.head(top_n).iter_rows(named=True), start=1)]
    lines += _table(["#", "Variante", "HGVS c.", "Composto", "Faixa", "Origem", "Modalidade dominante"], rows)
    lines += ["", "Faixa: `above_tau` (composto >= τ do estudo), `above_fixed` (>= menor limiar fixo e < τ). A figura D "
              "(`figures/D_prioritized.png`) mostra as variantes acima do menor limiar.", ""]
    links = study_dir / ATLAS_LINKS_FILE
    if links.is_file():
        lines += ["### Links do Atlas", "", links.read_text(encoding="utf-8").strip(), ""]
    interprets = sorted((study_dir / "raw").glob(INTERPRET_GLOB))
    if interprets:
        lines += ["### Interpretação por variante (AlphaGenome)", "",
                  "Resumo de cada `raw/interpret_*.txt` (cabeçalho, uso do sítio de splicing, alertas de junção e as "
                  "junções mais alteradas); o texto completo está no arquivo.", ""]
        for path in interprets:
            summary = summarize_interpretation(path.read_text(encoding="utf-8", errors="replace"))
            lines += [f"`{path.name}`", "", "```", *summary, "```", ""]
    return lines


def _section_skills(manifest: Mapping[str, Any]) -> list[str]:
    lines = ["## Skills utilizadas", "",
             "Lidas de `manifest.json` (cada passo registra a skill e o comando que gerou os arquivos de `raw/`):", ""]
    used = _skills_used(manifest)
    lines += _table(["Skill", "Passos do workflow"],
                    [(f"`{skill}`", ", ".join(f"`{s}`" + (f" (x{steps.count(s)})" if steps.count(s) > 1 else "")
                                              for s in dict.fromkeys(steps))) for skill, steps in used])
    lines += ["", f"Versão da skill `{manifest.get('skill', LOCAL_SKILL)}`: {manifest.get('skill_version', 'desconhecida')}.", ""]
    return lines


def _section_limits(study_dir: Path, study: Any, manifest: Mapping[str, Any], roc: RocResult,
                    rules: Sequence[core_clinical.ClassRule]) -> list[str]:
    cohort = pl.read_parquet(study_dir / core_clinical.COHORT_FILE)
    scored = pl.read_parquet(study_dir / core_threshold.COHORT_SCORED_FILE)
    lines = [
        "## Limitações", "",
        "- **Uso em pesquisa apenas.** São previsões de impacto molecular no splicing (AlphaGenome Atlas) para priorizar "
        "experimentos; não são medições nem achados clínicos e não devem ser usadas como diagnóstico.",
        "- **O composto é um ranking, não uma probabilidade.** O máximo sobre as trilhas mistura tecidos; leia a tabela de "
        "trilhas de uma variante antes de propor um mecanismo (ver `references/interpretation.md`).",
        f"- **Casos e controles são poucos:** {roc.n_case} casos e {roc.n_control} controles com escore; o τ e as métricas "
        "do ROC têm incerteza grande quando há poucos casos.",
        "- **Controles populacionais** (gnomAD) são tratados como controles, mas podem conter variantes patogênicas de baixa "
        "penetrância ou de início tardio; ClinVar enviesa para variantes já estudadas.",
    ]
    if not roc.calibrated:
        lines.append("- **O τ não foi calibrado** (sem casos ou sem controles); ele é o valor fixo escolhido pelo pesquisador.")
    origin = (manifest.get("origin") or {}).get("kind")
    gnomad_rows = cohort.filter(pl.col("sources").fill_null("").str.contains(core_clinical.GNOMAD))
    no_filter = gnomad_rows.filter(pl.col("gnomad_filters").is_null()).height
    wanted = (study.data.get("cohort") or {}).get("gnomad_filter")
    if origin != "snapshot" and no_filter:
        lines.append(f"- **Controles gnomAD ao vivo sem filtro PASS:** as saídas da skill `gnomad-database` não trazem o campo "
                     f"de filtros, então {no_filter:,} variantes gnomAD entraram sem verificar `{wanted or 'PASS'}`; "
                     "variantes de baixa qualidade podem estar entre os controles.")
    unscored_cohort = int(scored["composite"].is_null().sum())
    sat_path = study_dir / SATURATION_FILE
    unscored_sat = int(pl.read_parquet(sat_path)["composite"].is_null().sum()) if sat_path.is_file() else 0
    detail = ", ".join(f"{role} {n}" for role, n in sorted(roc.unscored.items())) or "nenhuma"
    lines.append(f"- **Sem composto (fora do ROC e da priorização, nunca lidas como 0):** {unscored_cohort} variantes da coorte "
                 f"(por papel: {detail}) e {unscored_sat} da varredura de saturação.")
    if origin == "snapshot":
        lines.append("- **Dados de referência (snapshot):** a coorte e a saturação foram copiadas de uma análise anterior; "
                     "a reprodução ao vivo pode diferir (ver `manifest.json`).")
    lines.append("- **Escopo:** só o(s) íntron(s) de interesse e os éxons vizinhos foram analisados; efeitos distantes ou em "
                 "outros transcritos ficam de fora, e CNVs e eventos em escala gênica não entram na coorte (`rejected.tsv`).")
    return lines + [""]


def _section_citations(manifest: Mapping[str, Any]) -> list[str]:
    lines = ["## Citações", "", f"- {PRIMARY_CITATION}"]
    for skill, _ in _skills_used(manifest):
        if skill in SKILL_CITATIONS:
            lines.append(f"- {SKILL_CITATIONS[skill]}")
        elif skill.split("/")[0] == LOCAL_SKILL:
            lines.append(f"- Script local da skill `{skill}` (splicing-analysis, versão {manifest.get('skill_version', '?')}).")
        else:
            lines.append(f"- Fonte `{skill}` (ver o `manifest.json` do estudo).")
    return lines + [""]


# --- entry point ----------------------------------------------------------------------------


def render_report(study_dir: Path) -> str:
    """Report text of a study that completed the prioritize step.

    Raises:
        ReportError: manifest.json or a processed file is missing.
        study.StudyError, core.* errors: a required study file is invalid.
    """
    study_dir = Path(study_dir)
    study = load_study(study_dir)
    locus = core_locus.load_locus(study_dir)
    manifest = _read_json(study_dir / MANIFEST_FILE, "record the workflow steps with manifest.py add")
    _need(study_dir, core_prioritize.PRIORITIZED_TSV, "run splicing.py prioritize (Workflow step 12)")
    roc = core_threshold.load_threshold(study_dir)
    rules_name, rules = core_clinical.resolve_rules(study, load_decisions(study_dir))
    top_n = int((study.data.get("report") or {}).get("top_n", DEFAULT_TOP_N))
    tau_name = "τ de Youden" if roc.calibrated else "τ fixo (não calibrado)"
    lines = [f"# Relatório: {locus.gene_symbol} ({study.study_id})", "",
             "> Uso em pesquisa apenas. Predições de impacto molecular no splicing, não diagnóstico.", ""]
    lines += _section_study(study, locus, manifest, rules_name)
    lines += _section_decisions(study_dir)
    lines += _section_cohort(study_dir, rules)
    lines += _section_threshold(roc, tau_name)
    lines += _section_prioritized(study_dir, roc, top_n)
    lines += _section_skills(manifest)
    lines += _section_limits(study_dir, study, manifest, roc, rules)
    lines += _section_citations(manifest)
    return "\n".join(lines)


def write_report(study_dir: Path) -> Path:
    """Writes <study_dir>/report.md and returns its path."""
    study_dir = Path(study_dir)
    path = study_dir / REPORT_FILE
    path.write_text(render_report(study_dir), encoding="utf-8")
    return path
