# Esercitazione di Metagenomica — Pipeline QIIME2

Questo repository contiene i materiali per l'esercitazione pratica di analisi metagenomica del microbioma (16S rRNA, regione V3-V4) tramite la pipeline **QIIME2**, usando dati paired-end Illumina.

---

## Campioni

| Sample ID | Gruppo | File R1 | File R2 |
|---|---|---|---|
| H79 | CTRL | `NG-A1943_V3V4a_H79_240117a_libLAF6019_1.fastq.gz` | `NG-A1943_V3V4a_H79_240117a_libLAF6019_2.fastq.gz` |
| H81 | TRT  | `NG-A1943_V3V4a_H81_240117a_libLAF6015_1.fastq.gz` | `NG-A1943_V3V4a_H81_240117a_libLAF6015_2.fastq.gz` |

---

## Requisiti

- [Git](https://git-scm.com/downloads) installato
- [Docker](https://docs.docker.com/get-docker/) installato e funzionante
- ~10 GB di spazio su disco

---

## 0. Setup iniziale

### 0a. Scarica il repository

Clona il repository con tutti i file necessari (dati, manifest, metadata):

```bash
git clone https://github.com/BCGL-Unime/esercitazione_meta.git
cd esercitazione_meta
```

> Se non hai Git installato, puoi scaricare il repository come archivio ZIP dalla pagina GitHub:
> **Code → Download ZIP**, poi estrai la cartella e aprila nel terminale.

### 0b. Scarica l'immagine Docker

Tutta la pipeline è preconfigurata in un'immagine Docker. Scaricala con:

```bash
docker pull anbonomo/esercitazione_qiime2
```

Avvia il container montando la cartella di lavoro corrente:

```bash
docker run -it --rm \
  -v $(pwd):/data \
  -w /data \
  anbonomo/esercitazione_qiime2 bash
```

> Tutti i comandi seguenti vanno eseguiti **all'interno del container**.

---

## Struttura della directory di lavoro

```
.
├── manifest/
│   └── manifest.tsv          # Manifest paired-end per l'import in QIIME2
├── metadata.txt              # Metadati: sample-id e gruppo (CTRL / TRT)
├── NG-A1943_V3V4a_H79_240117a_libLAF6019_1.fastq.gz   # CTRL — R1
├── NG-A1943_V3V4a_H79_240117a_libLAF6019_2.fastq.gz   # CTRL — R2
├── NG-A1943_V3V4a_H81_240117a_libLAF6015_1.fastq.gz   # TRT  — R1
└── NG-A1943_V3V4a_H81_240117a_libLAF6015_2.fastq.gz   # TRT  — R2
```

Gli output generati dalla pipeline vengono salvati nelle cartelle:
- `classification/` — tassonomia, tabelle collassate, analisi differenziale
- `diversity/` — albero filogenetico, metriche alfa/beta, PCoA plots

---

## Step 1 — Import dei dati paired-end

Il manifest (`manifest/manifest.tsv`) elenca i percorsi assoluti di forward (R1) e reverse (R2) per ciascun campione:

```tsv
sample-id	forward-absolute-filepath	reverse-absolute-filepath
H79	$PWD/NG-A1943_V3V4a_H79_240117a_libLAF6019_1.fastq.gz	$PWD/NG-A1943_V3V4a_H79_240117a_libLAF6019_2.fastq.gz
H81	$PWD/NG-A1943_V3V4a_H81_240117a_libLAF6015_1.fastq.gz	$PWD/NG-A1943_V3V4a_H81_240117a_libLAF6015_2.fastq.gz
```

Esegui l'import e genera il sommario:

```bash
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest/manifest.tsv \
  --output-path import.qza \
  --input-format PairedEndFastqManifestPhred33V2

qiime demux summarize \
  --i-data import.qza \
  --o-visualization summary_import.qzv
```

**Output:** `import.qza`, `summary_import.qzv`

> Apri `summary_import.qzv` su [view.qiime2.org](https://view.qiime2.org) per visualizzare la qualità delle read e scegliere i valori di `--p-trunc-len-f` e `--p-trunc-len-r` per lo step successivo.

---

## Step 2 — Denoising con DADA2 (paired-end)

DADA2 corregge gli errori di sequenziamento, rimuove le chimere e produce gli **ASV**. La modalità paired-end unisce le read forward e reverse per ottenere sequenze più accurate.

```bash
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs import.qza \
  --p-n-threads 0 \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-q 2 \
  --p-trunc-len-f 248 \
  --p-trunc-len-r 220 \
  --p-min-fold-parent-over-abundance 4 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza
```

> **`--p-trunc-len-f 248` / `--p-trunc-len-r 220`**: adatta questi valori in base alla qualità visualizzata in `summary_import.qzv`. Le read reverse tendono a calare prima in qualità, quindi il valore di R2 è spesso inferiore a quello di R1. Assicurati che le regioni troncate si sovrappongano di almeno 12 bp per permettere il merging.

Visualizzazioni:

```bash
qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file metadata.txt

qiime metadata tabulate \
  --m-input-file denoising-stats.qza \
  --o-visualization denoising-stats.qzv
```

**Output:** `table.qza`, `rep-seqs.qza`, `denoising-stats.qza`, `table.qzv`, `denoising-stats.qzv`

---

## Step 3 — Classificazione tassonomica e filtraggio

### 3a. Classificazione con SILVA

```bash
mkdir -p classification

qiime feature-classifier classify-sklearn \
  --i-reads rep-seqs.qza \
  --i-classifier silva-138-99-nb-classifier.qza \
  --o-classification classification/taxonomy_denoised.qza \
  --p-n-jobs 30
```

### 3b. Filtraggio (rimozione eucarioti e non assegnati)

```bash
qiime taxa filter-table \
  --i-table table.qza \
  --i-taxonomy classification/taxonomy_denoised.qza \
  --p-exclude eukaryota,unassigned \
  --o-filtered-table table-no-eukaryota-no-unassigned.qza
```

### 3c. Visualizzazioni

```bash
qiime feature-table summarize \
  --i-table table-no-eukaryota-no-unassigned.qza \
  --o-visualization table-no-eukaryota-no-unassigned.qzv \
  --m-sample-metadata-file metadata.txt

qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --i-taxonomy classification/taxonomy_denoised.qza \
  --o-visualization rep-seqs.qzv

qiime metadata tabulate \
  --m-input-file classification/taxonomy_denoised.qza \
  --o-visualization classification/taxonomy_denoised_metadata.qzv

qiime taxa barplot \
  --i-table table-no-eukaryota-no-unassigned.qza \
  --i-taxonomy classification/taxonomy_denoised.qza \
  --m-metadata-file metadata.txt \
  --output-dir classification/otu_final \
  --o-visualization classification/taxa-bar-plots.qzv
```

### 3d. Collasso per livello tassonomico e frequenze relative

I livelli SILVA corrispondono a: 1=Domain, 2=Phylum, 5=Family, 6=Genus.

```bash
for LEVEL in 2 5 6; do
  mkdir -p classification/level${LEVEL}/DA

  qiime taxa collapse \
    --i-table table-no-eukaryota-no-unassigned.qza \
    --i-taxonomy classification/taxonomy_denoised.qza \
    --p-level ${LEVEL} \
    --output-dir classification/level${LEVEL}

  qiime feature-table relative-frequency \
    --i-table classification/level${LEVEL}/collapsed_table.qza \
    --o-relative-frequency-table classification/level${LEVEL}/level${LEVEL}_relative.qza

  qiime metadata tabulate \
    --m-input-file classification/level${LEVEL}/level${LEVEL}_relative.qza \
    --o-visualization classification/level${LEVEL}/level${LEVEL}_relative_tabulated.qzv

  qiime tools export \
    --input-path classification/level${LEVEL}/collapsed_table.qza \
    --output-path classification/level${LEVEL}/

  biom convert \
    -i classification/level${LEVEL}/feature-table.biom \
    -o classification/level${LEVEL}/level${LEVEL}.tsv \
    --to-tsv
done
```

**Output:** tabelle collassate in formato `.qza`, `.qzv` e `.tsv` per i livelli 2, 5 e 6.

---

## Step 4 — Diversità alfa e beta

Usiamo `qiime diversity core-metrics` (versione non-filogenetica), che non richiede un albero evolutivo ed è sufficiente per le metriche composizionali standard.

### 4a. Metriche core (alfa + beta)

Prima di eseguire questo step, verifica il numero minimo di read per campione aprendo `table-no-eukaryota-no-unassigned.qzv` su [view.qiime2.org](https://view.qiime2.org) e sostituisci `INSERISCI_PROFONDITA` con quel valore.

```bash
mkdir -p diversity

qiime diversity core-metrics \
  --i-table table-no-eukaryota-no-unassigned.qza \
  --p-sampling-depth INSERISCI_PROFONDITA \
  --m-metadata-file metadata.txt \
  --output-dir diversity/core-metrics
```

Questo comando calcola in un colpo solo:

| Tipo | Metriche |
|---|---|
| **Alfa** | `observed_features`, `shannon`, `evenness` |
| **Beta** | `bray_curtis`, `jaccard` |
| **Visualizzazioni** | Emperor PCoA plots per ciascuna metrica beta |

**Output:** tutti i file in `diversity/core-metrics/`

### 4b. Curva di rarefazione

```bash
qiime diversity alpha-rarefaction \
  --i-table table-no-eukaryota-no-unassigned.qza \
  --p-max-depth INSERISCI_PROFONDITA \
  --m-metadata-file metadata.txt \
  --o-visualization diversity/alpha-rarefaction.qzv
```

> La curva di rarefazione mostra se la profondità di sequenziamento è sufficiente a catturare la diversità del campione (plateau = saturazione).

### 4c. Significatività statistica — alfa diversità

```bash
for METRIC in evenness shannon observed_features; do
  qiime diversity alpha-group-significance \
    --i-alpha-diversity diversity/core-metrics/${METRIC}_vector.qza \
    --m-metadata-file metadata.txt \
    --o-visualization diversity/core-metrics/${METRIC}-group-significance.qzv
done
```

> Test di Kruskal-Wallis: confronta la distribuzione di ciascuna metrica alfa tra i gruppi. Con solo 2 campioni il test non raggiunge potenza statistica, ma il workflow è lo stesso su dataset più ampi.

### 4d. Significatività statistica — beta diversità

```bash
for METRIC in bray_curtis jaccard; do
  qiime diversity beta-group-significance \
    --i-distance-matrix diversity/core-metrics/${METRIC}_distance_matrix.qza \
    --m-metadata-file metadata.txt \
    --m-metadata-column group \
    --p-pairwise \
    --o-visualization diversity/core-metrics/${METRIC}-group-significance.qzv
done
```

> PERMANOVA: testa se la composizione microbica differisce significativamente tra i gruppi. `--p-pairwise` produce confronti per coppia di gruppi.

**Output:** `diversity/core-metrics/*-group-significance.qzv` (apribili su view.qiime2.org)

---

## Step 5 — Analisi differenziale (ANCOM-BC)

ANCOM-BC identifica i taxa statisticamente differenziati tra il gruppo **TRT** e il gruppo di riferimento **CTRL**.

```bash
for LEVEL in 2 5 6; do
  qiime composition ancombc \
    --i-table classification/level${LEVEL}/collapsed_table.qza \
    --m-metadata-file metadata.txt \
    --p-formula 'group' \
    --p-reference-levels group::CTRL \
    --o-differentials classification/level${LEVEL}/DA/group_CTRL_ref_differentials.qza

  qiime composition da-barplot \
    --i-data classification/level${LEVEL}/DA/group_CTRL_ref_differentials.qza \
    --p-significance-threshold 0.05 \
    --p-level-delimiter ';' \
    --o-visualization classification/level${LEVEL}/DA/group_CTRL_ref_differentials.qzv
done
```

**Output:** `classification/levelX/DA/group_CTRL_ref_differentials.qza/.qzv` per i livelli 2, 5 e 6.

> I taxa con valore positivo sono enrichiti nel gruppo **TRT** rispetto a **CTRL**; quelli con valore negativo sono depleti.

---

## Visualizzazione dei risultati

Tutti i file `.qzv` possono essere aperti direttamente su **[view.qiime2.org](https://view.qiime2.org)** trascinando il file nella pagina.

---

## Riepilogo degli output principali

| File | Descrizione |
|---|---|
| `import.qza` | Dati paired-end importati in formato QIIME2 |
| `summary_import.qzv` | Qualità delle read forward e reverse |
| `table.qza` | Feature table (ASV) |
| `rep-seqs.qza` | Sequenze rappresentative degli ASV |
| `denoising-stats.qzv` | Statistiche del denoising (merged reads, chimere rimosse) |
| `classification/taxonomy_denoised.qza` | Classificazione tassonomica SILVA |
| `classification/taxa-bar-plots.qzv` | Barplot della composizione tassonomica |
| `classification/levelX/levelX.tsv` | Tabelle collassate per livello tassonomico (TSV) |
| `classification/levelX/DA/*.qzv` | Risultati ANCOM-BC: taxa differenziali CTRL vs TRT |
| `diversity/core-metrics/` | Metriche alfa (Shannon, evenness, observed_features) e beta (Bray-Curtis, Jaccard) + Emperor PCoA |
| `diversity/alpha-rarefaction.qzv` | Curve di rarefazione per la saturazione della diversità |
| `diversity/core-metrics/*-group-significance.qzv` | Test statistici (Kruskal-Wallis per alfa, PERMANOVA per beta) |
