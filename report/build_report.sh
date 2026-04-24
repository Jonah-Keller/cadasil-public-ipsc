#!/bin/bash
# Copy CSVs into report directory and render the Rmd
# Run from login node with R available

REPORT_DIR="/sc/arion/work/kellej10/cadasil_public_ipsc/report"
DATA_DIR="/sc/arion/work/kellej10/cadasil_public_ipsc"

# Copy all CSVs into report directory
cp "${DATA_DIR}/gene_expression_matrix.csv" "$REPORT_DIR/"
cp "${DATA_DIR}/sample_metadata.csv" "$REPORT_DIR/"
cp "${DATA_DIR}/transcript_variant_counts.csv" "$REPORT_DIR/"
cp "${DATA_DIR}/notch3_variant_summary.csv" "$REPORT_DIR/"
cp "${DATA_DIR}/genes_of_interest_counts.csv" "$REPORT_DIR/"

echo "CSVs copied to ${REPORT_DIR}/"
ls -lh "${REPORT_DIR}/"

echo ""
echo "To render the report, run:"
echo "  module load R"
echo "  cd ${REPORT_DIR}"
echo "  Rscript -e \"rmarkdown::render('E-MTAB-16303_reanalysis_report.Rmd')\""
echo ""
echo "To create a zip for distribution:"
echo "  cd ${DATA_DIR}"
echo "  zip -r E-MTAB-16303_reanalysis.zip report/"
