# Publishing checklist

Everything below is a one-time setup. The repository tree and the seven data tarballs
are already built and verified.

## 1. Create the GitHub repository

```bash
gh repo create krferrier/Methylation-EWAS-tutorial --public \
  --description "A reproducible eleven-chapter EWAS walkthrough: raw IDATs to functional annotation"
cd Methylation-EWAS-tutorial          # this directory
git remote add origin git@github.com:krferrier/Methylation-EWAS-tutorial.git
git push -u origin main
```

The commit is already made, so `git push` is all that is needed.

## 2. Turn on GitHub Pages

Settings → Pages → Source: **Deploy from a branch** → branch `main`, folder **`/docs`**.
The rendered site appears at `https://krferrier.github.io/Methylation-EWAS-tutorial/` within a minute.
`docs/.nojekyll` is committed so Jekyll does not strip the `_files` asset directories.

## 3. Deposit the data on Zenodo

The tarballs are in the `release/` directory of the build workspace:

```
ewas-tutorial-data-B_qc.tar.gz             587 MB
ewas-tutorial-data-C_normalized.tar.gz    1301 MB
ewas-tutorial-data-D_filtered.tar.gz      1102 MB
ewas-tutorial-data-E_model_inputs.tar.gz   511 MB
ewas-tutorial-data-F_ewas_results.tar.gz   117 MB
ewas-tutorial-data-G_pipeline_run.tar.gz   220 MB
ewas-tutorial-data-H_annotation.tar.gz     306 MB
                                          ------
                                          4.14 GB
```

Zenodo's per-record limit is 50 GB and the per-file limit is 50 GB, so all seven upload
to a single record with room to spare. Browser upload will struggle at 1.3 GB per file —
use the API instead:

```bash
# Create a draft deposition, then for each tarball:
BUCKET=<bucket url from the draft's links.bucket>
for f in release/ewas-tutorial-data-*.tar.gz; do
  curl -T "$f" "$BUCKET/$(basename $f)?access_token=$ZENODO_TOKEN"
done
```

`zenodo.json` in this repository holds the deposition metadata (title, description,
creators, keywords, license, related identifiers) — POST it to
`https://zenodo.org/api/deposit/depositions` to create the draft, then publish.

Consider **not** using the automatic GitHub–Zenodo integration here: it archives the
repository source, not these tarballs, which is the opposite of what you want. Create the
record manually and link it with the `isSupplementTo` relation already in `zenodo.json`.

## 4. Wire the DOI back into the repository — done

The record was published on 2026-08-28 as **22135216**. Zenodo issued two DOIs:

| DOI | resolves to | use for |
|---|---|---|
| `10.5281/zenodo.22135215` | newest version, always | prose, badges, citation |
| `10.5281/zenodo.22135216` | this exact set of files | reproducing *this* build |

Both are recorded in `CITATION.cff`; the concept DOI is the badge in `README.md`. The
numeric record id is the default in `get_data.sh`, so a student who clones the repo can
run `./get_data.sh` with no environment setup. Exporting `ZENODO_RECORD` still overrides
it, which is how you would pin an older version of the data after publishing a new one.

All seven tarball MD5s on the record match the local files in `release/`, so the upload
is verified byte-for-byte.

### One correction outstanding on the record

The published `related_identifiers` carry `isSupplementTo:
https://github.com/krferrier/EWAS-tutorial`, which 404s — the repository was renamed to
`Methylation-EWAS-tutorial` after that metadata was drafted. `zenodo.json` in this
repository has the correct URL. Fix it on the record via **Edit** on the record page
(metadata edits do not mint a new DOI or require re-uploading files).

### Publishing a new version later

Use **New version** on the record page, not a new record: the concept DOI keeps
resolving to the newest version, so the badge and `CITATION.cff` never need touching.
Only rebuild the tiers whose contents actually changed, then update `SHA256SUMS.txt`,
`MANIFEST.md`, and the record id in `get_data.sh` if you want clones to default to the
new version.

## 5. Verify the student experience

From a clean directory, as a student would:

```bash
git clone https://github.com/krferrier/Methylation-EWAS-tutorial.git
cd Methylation-EWAS-tutorial
./get_data.sh F_ewas_results H_annotation
bash qrender.sh render 08_annotation.qmd
```

That is the smallest download (422 MB) that exercises the fetch script, the checksum
verification, the extraction paths and a chapter render.
