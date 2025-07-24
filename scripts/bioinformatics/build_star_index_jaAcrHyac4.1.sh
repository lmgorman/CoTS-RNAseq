build_star_index_jaAcrHyac4.1.sh
if [ ! -f "$GTF_FILE" ]; then
    echo "Converting GFF3 to GTF..."
    gffread "$GFF3_FILE" -T -o "$GTF_FILE"
    if [ $? -ne 0 ]; then
        echo "ERROR: gffread failed to convert GFF3 to GTF"
        exit 1
    fi
else
    echo "GTF file already exists, skipping conversion."
fi

# Create output directory
mkdir -p "$INDEX_DIR"

# Run STAR genomeGenerate with annotation
STAR --runMode genomeGenerate \
    --genomeDir "$INDEX_DIR" \
    --genomeFastaFiles "$GENOME_FA" \
    --sjdbGTFfile "$GTF_FILE" \
    --runThreadN 8 \
    --sjdbOverhang 99

echo "[$(date)] Genome index build completed"
