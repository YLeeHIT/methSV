#!/bin/bash
# ===============================
# Script: download_giab_pod5.sh
# Description:
#   Download GIAB HG002–HG007 pod5 files from ONT S3 bucket
#
# Usage:
#   ./download_giab_pod5.sh
#
# Notes:
#   - Requires AWS CLI with --no-sign-request support
#   - Outputs are saved in directories named HG002 to HG007
# ===============================

# ----------- Main Download Loop ----------
for i in {2..7}; do
    sample_id="HG00${i}"
    dest_dir="./${sample_id}"

    echo "[INFO] Downloading pod5 files for ${sample_id} into ${dest_dir}..."

    # Create destination directory if it doesn't exist
    mkdir -p "${dest_dir}"

    # Download using AWS S3 CLI (public access, no auth required)
    aws s3 cp "s3://ont-open-data/giab_2025.01/flowcells/${sample_id}/" "${dest_dir}/" --recursive --no-sign-request

    # Check if the download was successful
    if [ $? -eq 0  ]; then
        echo "[INFO] Successfully downloaded data for ${sample_id}"
    else
        echo "[ERROR] Failed to download data for ${sample_id}" >&2
    fi
done

echo "[INFO] All downloads completed."
