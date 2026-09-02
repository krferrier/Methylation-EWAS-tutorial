# Maintainer notes

**These are working notes for the tutorial's author, not part of the tutorial.**

If you are here to learn how to run an EWAS, you want the rendered tutorial at
<https://krferrier.github.io/Methylation-EWAS-tutorial/> or the chapter sources
in [`../tutorial/`](../tutorial/). Nothing in this directory is needed to read,
render, or reuse the tutorial, and some of it will not make sense without the
context of how the tutorial was built.

They are kept in the open because the reasoning behind an analysis is usually
more instructive than the analysis itself — if you want to see *why* a
particular mask version, covariate set, or threshold was chosen, and what was
tried and rejected along the way, that record is here rather than hidden.

| file | what it is |
|---|---|
| `DECISIONS.md` | the analysis decisions — mask source and columns, probe scope, cell-composition handling, smoking, SVA and ComBat rulings — and the reasoning for each |
| `CHANGES.md` | a running log of every revision round: what was edited, what was recomputed, and which published numbers moved as a result |
| `PUBLISHING.md` | the operational checklist for releasing the site and depositing the checkpoint data on Zenodo |

A caution if you read `CHANGES.md`: it is written chronologically, so it
records numbers that were later superseded. The authoritative values are always
the ones in the rendered tutorial, never a figure quoted in an old changelog
entry.
