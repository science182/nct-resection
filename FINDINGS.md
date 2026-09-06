# Degree, weighting, and the stability of deletion-based connectome maps

**Standalone summary. For the full working record including failures, see [README.md](README.md).**

A short account of what was tested, what held, what did not, and what follows
for anyone computing resection risk from a structural connectome.

Sujaal Gelle. Code, data acquisition and verification:
https://github.com/science182/nct-resection

---

## Summary

Simulated resection analyses rank brain parcels by how much removing them
degrades a network measure. This tested whether network controllability adds
anything to that ranking beyond node degree, using the contiguous deletion
design of Lin et al. (2024) on 1065 HCP connectomes and a second, independent
70-subject dataset.

Three results.

**Deletion damage is largely node degree.** Average controllability damage
correlates with node strength at +0.896 within subject in HCP and +0.883 in an
unrelated pipeline and parcellation. PageRank correlates at +0.997. Global
efficiency is the least degree-dependent of the three at +0.539.

**Degree-corrected maps can be determined by the edge weighting convention.**
The same 1065 brains, weighted two standard ways, produce risk maps that are
uncorrelated (+0.021) while each is reproducible across independent halves at
+0.999. One flags left perisylvian language cortex (9 of the top 36 parcels
against 2.3 expected, spin test p = 0.005); the other flags ventromedial visual
cortex (25 of 36 against 6.0, p = 0.008). Each survives a spin test and
Bonferroni correction. At most one can be true.

**The cause is the weight distribution, not density.** On a grid separating tail
heaviness from density, stability tracks log tail ratio at Spearman -0.917 while
density contributes almost nothing once the tail is known (partial +0.267).
Compressing the tail at full density lifts stability from +0.319 to +0.843.
Thinning edges at an already-compressed tail does nothing (-0.070).

The practical consequence is a diagnostic and a remedy, both given below.

---

## What was done

**Data.** Rosen & Halgren (2021) HCP connectomes, HCP-MMP1 360 parcels, 1065
subjects, in two native weightings: fractional probability (Fpt) and raw
streamline counts. Independent replication used a 70-subject Lausanne-atlas
dataset (Zenodo 2872624) from a different site, acquisition, pipeline, and
parcellation, and at 11 percent density against HCP's fully dense matrices.

**Deletion design.** Following Lin et al., parcels are removed and the damage to
a network measure is scored over the surviving parcels. Contiguity uses true
HCP-MMP1 surface adjacency, derived from the fsaverage white-matter mesh and the
parcellation annotation files (1050 borders, mean parcel degree 5.83, consistent
with a cortical tessellation).

**Measures.** Average and modal controllability (Gu et al. 2015), weighted
global efficiency, PageRank, node strength, and minimum control energy.

**Degree correction.** Removing a parcel removes exactly its own strength, so
every raw damage measure inherits degree. The corrected score is the residual of
damage on within-subject strength rank. The result is invariant to how that
correction is done: linear, quadratic, cubic and rank-only give the same
enrichment at p ~ 0.01, and a regression-free strength-matched-pair test agrees.

**Nulls.** Network enrichment uses a spin test: 1000 rotations of the map on the
fsaverage sphere with Hungarian matching for genuine one-to-one parcel
permutations. This preserves spatial autocorrelation, which free label shuffling
does not; most secondary enrichments that looked significant under the naive
null did not survive it.

---

## What did not work

Recorded because the negative results are most of the value.

**Control energy.** With full control the network barely enters, and the effect
was under 0.01 percent of the intact energy with the sign wrong in 66 percent of
cases. Restricting the driver set makes the network matter but leaves the
Gramian near-singular, energies of order 1e12, and sign violations still near
half. It needs reachable target states from measured activation and a
regularized formulation, neither of which was available here.

**Modal controllability.** Less degree-dependent (-0.483) but at any
normalization where lesion deltas are usable it is close to a relabeling of
average controllability, so it is not independent evidence.

**Individual variation.** Nodal controllability does carry subject-specific
structure beyond strength: the residual is 2.17 times a parcel-shuffled null and
yields 159 distinct top parcels against strength's 23. That does not transfer to
resection damage, where the individual-specific gain is +0.039.

---

## What this means in practice

Four changes, in order of how cheap they are to make.

**1. Transform the weights before computing network measures.** This is one line
and loses nothing. Raw streamline counts have a tail ratio around 10,000, which
is the least stable regime measured. Raising weights to about the 0.35 power, or
taking log10(1 + w) when counts are large, brings the ratio under 30 and takes
degree-corrected map stability from +0.32 to +0.84.

Do not threshold instead. Thresholding to 11 percent density reaches only +0.64
and discards 89 percent of the edges to get there. It works at all only because
it compresses the tail as a side effect.

**2. Report node strength beside any network measure.** If a network finding
correlates with strength above about 0.9, summing a row of the connectivity
matrix would have produced it. PageRank at +0.997 is the clearest case: hubness
results on dense weighted connectomes are close to statements about degree, and
saying so costs nothing and prevents an overclaim.

**3. Where the risk actually lies, and where it does not.** Raw global
efficiency is reasonably stable across weighting conventions (+0.58), comparable
to raw controllability (+0.57). The instability is specific to *degree-corrected*
derivatives, which fall to +0.04 and +0.29. So an existing raw-measure pipeline
is not the thing at risk. The exposure is in the natural next step: correcting
for degree in order to surface regions that plain connectivity misses. That is
exactly where the map becomes convention-dependent.

**4. Report comparisons with a robustness grade, not a point estimate.** The
clinical question is usually whether corridor A is safer than corridor B, and
that can be robust even when the whole-brain map is not. Scoring each comparison
across weighting conventions and damage measures separates decisions that hold
under every choice from those that flip. Measured here, clearly different
corridors gave 5 of 6 robust comparisons with no coin flips, while corridors
matched on total connectivity removed gave 1 robust and 1 coin flip.

`audit.py` implements 1 and 2 as a per-dataset check; `resection_report.py`
implements 4.

---

## What is new, and what is not

Most of this replicates known results, and the repository says so.

Average controllability tracking node strength is documented: Gu et al. (2015)
report r = 0.91, and the +0.896 here reproduces it. That edge weighting affects
graph metrics is established, with its own literature; a 2024 Network
Neuroscience comparison of streamline count, FA and axon-diameter weightings
concludes that the choice affects interpretation and that no one weighting is
superior. Log-transforming skewed
connectome weights is already common practice. Multiverse analysis is a
recognised methodology, including in network neuroscience.

There is also a standing rebuttal. Parkes et al. argue the strength correlation
is spatial rather than between-subject, and that average controllability beats
strength at out-of-sample prediction. Testing both axes here, the
degree-corrected measure fails on both: spatial agreement +0.021 and
between-subject +0.059 against a shuffled null of -0.008, with 149 of 360
parcels negative. That does not refute their prediction result, which measured
something different, but it does mean the defence does not rescue this
application.

What appears not to have been done is the combination: the tail ratio as a
quantitative predictor of when weighting choice determines the answer, applied
to control measures, in the deletion and resection setting. The network control
theory protocol paper names edge-weight variability as a limitation without
measuring it, and the multiverse tooling that exists for network neuroscience
does not include controllability measures.

So the contribution is a mechanism and a diagnostic for a known problem, not a
new phenomenon. It is offered as a tool rather than a discovery.

---

## Limitations

Two datasets, one of them 70 subjects. The Lausanne release ships a single
native weighting, so its four conventions are derived rather than native, and
the HCP figure used for that comparison is the matching derived-weighting number.

The identification grid used 8 subjects per cell and 9 cells. The density sweep
used 10 to 12 per condition. Both are adequate to estimate a map correlation and
not to put a tight interval on it.

Global efficiency is this implementation, not the original authors', so the
comparative statements about it should be read as a flag rather than a result.

No outcome data. The central clinical question, whether any of these maps
predict post-operative deficits, cannot be answered from public connectomes.
That requires resection extent and domain-specific neuropsychological outcomes.

---

## Reproducing

```bash
pip install -r requirements.txt
python3 download_data.py --full     # public sources, ~1.3 GB
python3 reproduce.py                # checks every number quoted above
```

`reproduce.py` recomputes each headline figure and compares it against the
written value, exiting non-zero on any mismatch.

Nine separate results in this project were clean, significant, and wrong before
a control caught them: a normalization that made hub resection look protective,
a unit-induced collapse of two measures into one, a spectral offset that swamped
the signal and flipped a correlation's sign, an invalid fingerprinting design
that returned chance for a known-positive control, an anti-conservative
enrichment null, and a causal reading of density that an orthogonal manipulation
overturned. Each is documented in the repository with the test that now guards
it. That record is the main thing worth taking from this work.

---

## References

Gu S, et al. Controllability of structural brain networks. Nat Commun 6:8414 (2015).
Parkes L, et al. A network control theory pipeline for studying the dynamics of the structural connectome. Nat Protoc 19(12) (2024).
Lin YH, et al. Discernible interindividual patterns of global efficiency decline during theoretical brain surgery. Sci Rep 14:14573 (2024).
Yeung JT, et al. Unexpected hubness: a proof-of-concept study of the human connectome using PageRank centrality. J Neurooncol 151(2):249-256 (2021).
Rosen BQ, Halgren E. A whole-cortex probabilistic diffusion tractography connectome. eNeuro 8(1) (2021).
Alexander-Bloch AF, et al. On testing for spatial correspondence between maps of human brain structure and function. NeuroImage 178:540-551 (2018).
