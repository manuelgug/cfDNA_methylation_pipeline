/*
 * Module: VALIDATE_SAMPLESHEET
 * Validates input CSV structure and file paths before pipeline starts.
 */

process VALIDATE_SAMPLESHEET {
    tag "samplesheet"
    label 'process_low'

    publishDir "${params.outdir}/pipeline_info", mode: params.publish_mode

    input:
    path samplesheet

    output:
    path "validated_samplesheet.csv", emit: validated_csv
    path "samplesheet_report.txt",    emit: report

    script:
    """
    python3 - << 'PYEOF'
import csv, sys, os, re
from pathlib import Path

required_cols = {"sample_id","patient_id","condition","sample_type",
                 "library_type","read1"}
valid_conditions   = {"CRC","control","adenoma","polyp"}
valid_sample_types = {"cfDNA","tissue","cell_line"}
valid_libraries    = {"WGBS","targeted","hybrid_capture","RRBS"}

errors   = []
warnings = []
rows_out = []

with open("${samplesheet}") as fh:
    reader = csv.DictReader(fh)
    cols   = set(reader.fieldnames or [])

    missing = required_cols - cols
    if missing:
        sys.exit(f"FATAL: missing columns: {missing}")

    seen_ids = set()
    for i, row in enumerate(reader, start=2):
        sid = row["sample_id"].strip()
        if not sid:
            errors.append(f"Row {i}: empty sample_id")
            continue
        if sid in seen_ids:
            errors.append(f"Row {i}: duplicate sample_id '{sid}'")
        seen_ids.add(sid)

        # Condition check
        cond = row["condition"].strip()
        if cond not in valid_conditions:
            warnings.append(f"Row {i} [{sid}]: unexpected condition '{cond}'")

        # Library check
        lib = row["library_type"].strip()
        if lib not in valid_libraries:
            warnings.append(f"Row {i} [{sid}]: unexpected library_type '{lib}'")

        # File existence
        for col in ("read1","read2"):
            fp = row.get(col,"").strip()
            if fp and not Path(fp).exists():
                errors.append(f"Row {i} [{sid}]: {col} not found: {fp}")

        rows_out.append(row)

report_lines = ["=== Samplesheet Validation Report ===",
                f"Samples    : {len(rows_out)}",
                f"Errors     : {len(errors)}",
                f"Warnings   : {len(warnings)}",""]
report_lines += [f"ERROR: {e}" for e in errors]
report_lines += [f"WARN:  {w}" for w in warnings]

with open("samplesheet_report.txt","w") as fh:
    fh.write("\\n".join(report_lines))

print("\\n".join(report_lines))

if errors:
    sys.exit(f"Samplesheet validation FAILED with {len(errors)} error(s).")

# Write validated copy
with open("validated_samplesheet.csv","w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=reader.fieldnames)
    writer.writeheader()
    writer.writerows(rows_out)

print(f"\\n✅  Samplesheet OK: {len(rows_out)} samples validated.")
PYEOF
    """
}
