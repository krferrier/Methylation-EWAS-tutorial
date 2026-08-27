# Publishing checklist

Everything below is a one-time setup. The repository tree and the seven data tarballs
are already built and verified.

## 1. Create the GitHub repository

```bash
gh repo create krferrier/EWAS-tutorial --public \
  --description "A reproducible eleven-chapter EWAS walkthrough: raw IDATs to functional annotation"
cd EWAS-tutorial          # this directory
git remote add origin git@github.com:krferrier/EWAS-tutorial.git
git push -u origin main
```

The commit is already made, so `git push` is all that is needed.

## 2. Turn on GitHub Pages

Settings → Pages → Source: **Deploy from a branch** → branch `main`, folder **`/docs`**.
The rendered site appears at `https://krferrier.github.io/EWAS-tutorial/` within a minute.
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
ewas-tutorial-data-H_annotation.tar.gz     305 MB
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

## 4. Wire the DOI back into the repository

After publishing, Zenodo issues two DOIs: a **concept DOI** that always resolves to the
newest version, and a version DOI. Use the concept DOI in prose.

```bash
# in get_data.sh, replace REPLACE_WITH_RECORD_ID with the numeric record id,
# or have users export ZENODO_RECORD themselves.
sed -i 's/REPLACE_WITH_RECORD_ID/1234567/' get_data.sh
```

Then add the badge to the top of `README.md`:

```markdown
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1234567.svg)](https://doi.org/10.5281/zenodo.1234567)
```

and commit both changes.

## 5. Verify the student experience

From a clean directory, as a student would:

```bash
git clone https://github.com/krferrier/EWAS-tutorial.git
cd EWAS-tutorial
export ZENODO_RECORD=1234567
./get_data.sh F_ewas_results H_annotation
bash qrender.sh render 08_annotation.qmd
```

That is the smallest download (422 MB) that exercises the fetch script, the checksum
verification, the extraction paths and a chapter render.
