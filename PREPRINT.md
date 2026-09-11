# Edge-weight distribution determines the conclusions of node-deletion analyses in structural connectomes

**Sujaal Gelle**

Independent researcher. Correspondence: gellesujaal@gmail.com
Code, data acquisition and verification: https://github.com/science182/nct-resection

---

## Abstract

Node-deletion analyses rank brain regions by how much removing them degrades a
network measure, and are used to reason about which tissue is safe to resect.
Such analyses require an edge-weighting convention, a choice usually reported in
a single clause of a methods section. We asked how much that choice determines
the result.

Using 1065 Human Connectome Project subjects parcellated into 360 regions, we
computed deletion damage under average controllability, global efficiency,
PageRank and node strength. Deletion damage was largely a restatement of node
degree: average controllability tracked node strength at rho = +0.896 within
subject, PageRank at +0.997, global efficiency at +0.539. This replicates a
known spatial correlation.

Correcting damage for degree produces a map that is highly reproducible within a
weighting convention and incompatible across conventions. Under fractional
probability weighting the corrected map concentrated on left perisylvian
language cortex (9 of the top 36 parcels against 2.3 expected, spin test
p = 0.005); under raw streamline counts, on the same brains, it concentrated on
ventromedial visual cortex (25 of 36 against 6.0, p = 0.008) with zero language
parcels. The two maps correlate at +0.021 while each reproduces across
independent halves of the cohort at +0.999.

Separating tail heaviness from density on a factorial grid, map stability
tracked log tail ratio at rho = -0.917 while density contributed little once the
tail was known. Compressing the weight distribution at full density raised
stability from +0.319 to +0.843; thinning edges at an already-compressed tail
did not help.

At the level of conclusions rather than correlations, the cost is substantial.
Comparing the two native weightings, the ten highest-risk parcels shared 43
percent of members for raw damage and none at all after degree correction, and
clear-cut pairwise risk judgements reversed in 29.6 and 47.6 percent of cases
respectively. This is not confined to degree-corrected derivatives: for raw
deletion damage across four conventions, the single most central parcel and the
most over-represented network each changed in 83 percent of comparisons, at a
map correlation of +0.58.

We conclude that reporting map correlation understates how contingent a
conclusion is, that the relevant diagnostic is the weight distribution's tail,
and that compressing weights is preferable to thresholding as a remedy. Two
open-source tools implement the diagnostic and a robustness-graded comparison of
candidate resections.

---

## 1. Introduction

Structural connectomes derived from diffusion tractography are increasingly used
to reason about surgical risk. A common design removes a region, or a contiguous
set of regions, recomputes a network measure, and ranks candidate resections by
the resulting damage [1]. Related work identifies hub regions by centrality and
argues that damage to them carries disproportionate cost [2].

Every such analysis requires a decision about what an edge weight means.
Tractography can be summarized as a raw streamline count, as a fraction of
streamlines seeded, as a density normalized by region size and streamline
length, or as a diffusion scalar. These are different measurements, not
rescalings of one another, and comparisons of weighting schemes have concluded
that the choice affects interpretation without one being superior [3].

Separately, network control theory has been applied to connectomes, and average
controllability is known to correlate strongly with weighted degree [4]. A
protocol paper for these methods names variability in edge-weight distributions
across preprocessing pipelines as a limitation, without quantifying it [5].

We combined these concerns. Our question was not whether weighting matters in
principle, which is established, but how much of a specific, clinically framed
conclusion it determines in practice, and whether the sensitivity can be
predicted from a property of the data before any analysis is run.

## 2. Methods

### 2.1 Data

Primary analyses used the connectomes of Rosen and Halgren [6]: 1065 Human
Connectome Project subjects, HCP-MMP1 parcellation, 360 cortical parcels,
released in two native weightings. Fpt is the fraction of streamlines seeded
from a parcel that reach a target; the raw streamline count file reports counts
directly. Subject ordering was verified identical between the two files before
any paired comparison.

Independent replication used 70 subjects from a Lausanne-atlas release [7], at
scale 219, differing in site, acquisition, tractography pipeline, parcellation,
edge-weight definition (fiber density normalized by streamline length and region
surface area), and density (0.111 against 1.000).

### 2.2 Deletion design

Following [1], parcels were removed and damage scored over the surviving
parcels. Anatomical contiguity for multi-parcel resections used true HCP-MMP1
adjacency derived from the fsaverage white-matter surface and the parcellation
annotation files, giving 1050 parcel borders with mean parcel degree 5.83,
consistent with a cortical tessellation. Adjacency was block diagonal by
hemisphere, since cortical tissue is not contiguous across the midline.

### 2.3 Network measures

Average and modal controllability follow [4], computed from the eigendecomposition
of the adjacency matrix normalized to a target spectral radius. Global efficiency
is the weighted form, with edge weights inverted to distances. PageRank and node
strength are standard.

Two normalization details materially affect lesion analyses. First, both networks
must be normalized by the *intact* network's spectral scale: normalizing each by
its own largest singular value means removing a hub lowers the denominator,
inflating the lesioned network's controllability and making hub resection appear
protective. Second, the target spectral radius must suit the weight units. The
conventional additive form assumes streamline counts, and applied to probability
weights it places the system in a regime where average and modal controllability
both linearize around the same quantity and become mutually redundant.

### 2.4 Degree correction

Removing a parcel removes exactly its own strength, so every raw damage measure
inherits degree by construction. The corrected score is the residual of damage on
within-subject strength rank. Results were invariant to the form of that
correction: linear, quadratic, cubic and rank-only fits gave identical network
enrichment at p ~ 0.01, and a regression-free comparison restricted to
strength-matched parcel pairs agreed.

### 2.5 Weighting conventions

Four conventions were derived from a single matrix, holding tractography fixed:
raw, log10(1 + w), row-normalized, and binarized at the top 30 percent of edges.
Because these are derived, the two *native* weightings of the same subjects, Fpt
against raw streamline counts, were compared separately as the stronger test.

To separate tail heaviness from density, base connectomes were built on a
factorial grid before any convention was applied: weights raised to a power alpha
in {1.0, 0.5, 0.25}, which compresses the distribution without removing an edge,
crossed with thinning to densities {1.0, 0.30, 0.11}, which removes edges while
leaving surviving weights untouched. Tail heaviness is summarized as the ratio of
the maximum to the median nonzero off-diagonal weight.

### 2.6 Conclusion types

Five statement types were scored for stability, chosen to match what is
routinely reported: the identity of the top-ranked parcel; the membership of the
top ten; the full rank order by Spearman correlation; pairwise judgements of the
form "removing A is worse than removing B"; and which network is most
over-represented among the highest-risk parcels.

Pairwise judgements were counted over "clear-cut" pairs only, defined as those
the first map separates by more than half its interquartile spread. Reversing a
near-tie is not a substantive disagreement, and counting all pairs inflates the
reported rate by roughly a third. Both figures are reported in the code output.

### 2.7 Spatial nulls

Network enrichment was tested against a spin null: 1000 random rotations of the
map on the fsaverage sphere, with parcels matched to rotated positions by the
Hungarian algorithm to produce genuine one-to-one permutations [8]. This
preserves the spatial autocorrelation of the map and the geometry of the
parcellation, and destroys only the alignment between map and anatomy. Free
label shuffling is anti-conservative here, and most secondary enrichments
significant under that null did not survive the spin test.

## 3. Results

### 3.1 Deletion damage is largely node degree

Across 1065 subjects, average controllability deletion damage correlated with
node strength at rho = +0.896 (sd 0.018) within subject. In the independent
Lausanne dataset, with a different pipeline and parcellation, the same
correlation was +0.883. PageRank correlated with strength at +0.997; global
efficiency, the least degree-dependent of the measures examined, at +0.539.

This held for single-parcel and for anatomically contiguous multi-parcel
resections, at group level and per subject. The spatial correlation between
average controllability and weighted degree is documented [4]; the contribution
here is that the deletion delta inherits it more strongly than the nodal measure
does (+0.896 against +0.654), which follows from removing a parcel removing
exactly its own strength.

### 3.2 Two native weightings give incompatible maps

Degree-corrected maps were highly reproducible within a convention. Splitting the
cohort into independent halves of 532 and 533 subjects, the mean corrected map
of one half correlated with the other at +0.999.

Across the two native conventions, on the same subjects, the corrected maps
correlated at +0.021. Under Fpt, 9 of the 36 highest-risk parcels belonged to the
language network against 2.3 expected (spin test p = 0.005), and none to the
visual network. Under raw streamline counts, 25 of 36 belonged to the visual
network against 6.0 expected (p = 0.008), and none to language. Both survive
Bonferroni correction across ten networks. At most one can reflect anatomy.

Testing both axes of comparison, the corrected measure failed on each: spatial
agreement between conventions +0.021, between-subject agreement +0.059 against a
subject-shuffled null of -0.008, with 149 of 360 parcels negative. Node strength,
by contrast, was stable on both axes (+0.870 and +0.729).

### 3.3 Tail heaviness, not density, determines stability

On the factorial grid, map stability across conventions tracked log tail ratio at
rho = -0.917, while density tracked it at -0.211 and contributed a partial
correlation of +0.267 once the tail was known. Compressing the tail at full
density raised stability from +0.319 to +0.843, exceeding what the sparser
Lausanne connectomes achieve. Thinning edges at an already-compressed tail
changed stability by -0.070, that is, not at all.

Cells with matched tail ratios but very different densities agreed closely: tail
193 at density 0.11 gave +0.638, tail 94 at density 1.00 gave +0.629. An earlier
analysis varying only density appeared to show density as the cause; that sweep
confounded the two, since thresholding discards the smallest weights and thereby
compresses the tail as a side effect.

### 3.4 What the choice costs

Correlations understate the instability. For raw deletion damage across four
conventions in HCP, the top-ranked parcel changed in 83 percent of comparisons,
top-ten sets overlapped 54 percent, clear-cut pairwise judgements reversed 21.1
percent of the time, and the most over-represented network changed in 83 percent
of comparisons. This occurred at a mean map correlation of +0.58.

After degree correction the same figures were 100 percent, 11 percent, 35.7
percent and 83 percent.

Between the two native weightings the contrast is sharpest. For raw damage,
top-ten overlap was 43 percent and 29.6 percent of clear-cut judgements
reversed. After degree correction, the ten highest-risk parcels shared no members
at all and 47.6 percent of clear-cut judgements reversed, which is
indistinguishable from chance.

### 3.5 Independent dataset

The degree result replicated (+0.883). The instability did not, in the same
magnitude: degree-corrected map stability across conventions was +0.756 in the
Lausanne connectomes against +0.289 in HCP. This is consistent with the
mechanism, since the Lausanne matrices have a tail ratio near 190 against roughly
10,000 for raw HCP streamline counts. Conclusion-level flip rates were
correspondingly lower but not negligible: top hub changed in 83 percent of
comparisons, top-ten overlap was 68 percent, and 7.6 percent of clear-cut
judgements reversed.

### 3.6 Measures that did not work

Minimum control energy was examined as a candidate that is not a spectral summary
of the matrix. With full control the network barely enters: the effect of a
deletion was under 0.01 percent of the intact energy, with the sign contradicting
the model in 66 percent of cases. Restricting the driver set makes the network
matter but leaves the Gramian near-singular, with energies of order 1e12 and sign
violations near half. A usable formulation requires reachable target states from
measured activation and a regularized objective.

Modal controllability is less degree-dependent (-0.483) but at any normalization
where lesion deltas are well conditioned it is close to a relabeling of average
controllability.

## 4. Discussion

### 4.1 What is and is not new

The degree correlation is documented [4]. That weighting affects graph metrics is
established [3]. Log-transforming skewed connectome weights is common practice,
and multiverse analysis is a recognised methodology.

There is also a standing counter-argument: that the strength correlation is
spatial rather than between-subject, and that average controllability
outperforms strength at out-of-sample prediction of clinical variables [5]. Our
analyses are spatial, which is the axis where the correlation is conceded. We do
not test prediction of an external variable and therefore do not contradict that
result. We do show that the defence does not extend to this application, since
the corrected measure is unstable on the between-subject axis as well.

What appears not to have been quantified is the combination reported here: the
tail ratio as a predictor of when weighting determines the answer, applied to
control measures, in the deletion setting, and expressed as flip rates over
conclusion types rather than map correlations.

### 4.2 Recommendations

**Transform the weights rather than thresholding.** Raw streamline counts have
tail ratios near 10,000. Raising weights to approximately the 0.35 power, or a
log transform where counts are large, brings the ratio below 30 and lifts
stability from +0.32 to +0.84 while retaining every edge. Thresholding to 11
percent density reaches only +0.64 and discards 89 percent of the graph.

**Report node strength alongside any network measure.** Where a finding
correlates with strength above roughly 0.9, the simpler quantity would have
produced it.

**Do not report map correlation as evidence of conclusion stability.** A
correlation of +0.58 coexisted here with the top hub and the dominant network
each changing in five comparisons out of six.

**Grade comparisons rather than reporting point estimates.** Scoring a candidate
comparison across conventions and measures separates judgements that hold under
every choice from those that do not. In our data, clearly separated resection
candidates gave 5 of 6 robust comparisons with no coin flips, while candidates
matched on total connectivity removed gave 1 robust and 1 coin flip. An unstable
whole-brain map does not by itself make every specific comparison unreliable, and
distinguishing the two is tractable.

### 4.3 Limitations

Two datasets, one of 70 subjects. The Lausanne release ships a single native
weighting, so its four conventions are derived; the native two-weighting
comparison exists only in HCP.

The factorial grid used 8 subjects per cell across 9 cells, and the density sweep
10 to 12 per condition. These are adequate to estimate a map correlation and not
to place a tight interval on it.

Global efficiency is our implementation rather than the original authors', so
comparative statements about it should be read as a flag rather than a result.

Most importantly, there is no outcome data. Whether any of these maps predicts
post-operative deficit is the question that matters clinically and cannot be
answered from public connectomes. It requires resection extent and
domain-specific neuropsychological outcomes. The present work establishes only
that the maps disagree with one another, not which, if any, is right.

## 5. Data and code availability

All inputs are public. `download_data.py` retrieves them from Zenodo and
TemplateFlow. `reproduce.py` recomputes every numerical claim in this manuscript
and compares it against the reported value, exiting non-zero on any mismatch; 22
of 22 currently pass. Two command-line tools, `audit.py` and
`resection_report.py`, implement the diagnostic and the graded comparison.

https://github.com/science182/nct-resection

## References

[1] Lin YH, Dadario NB, Tang SJ, et al. Discernible interindividual patterns of global efficiency decline during theoretical brain surgery. Sci Rep. 2024;14:14573.

[2] Yeung JT, Taylor HM, Young IM, Nicholas PJ, Doyen S, Sughrue ME. Unexpected hubness: a proof-of-concept study of the human connectome using PageRank centrality and implications for intracerebral neurosurgery. J Neurooncol. 2021;151(2):249-256.

[3] Weighting the structural connectome: exploring its impact on network properties and predicting cognitive performance in the human brain. Netw Neurosci. 2024;8(1):119-137.

[4] Gu S, Pasqualetti F, Cieslak M, et al. Controllability of structural brain networks. Nat Commun. 2015;6:8414.

[5] Parkes L, Kim JZ, Stiso J, et al. A network control theory pipeline for studying the dynamics of the structural connectome. Nat Protoc. 2024;19(12).

[6] Rosen BQ, Halgren E. A whole-cortex probabilistic diffusion tractography connectome. eNeuro. 2021;8(1):ENEURO.0416-20.2020.

[7] Structural and functional connectome of 70 healthy adults, Lausanne atlas. Zenodo 2872624.

[8] Alexander-Bloch AF, Shou H, Liu S, et al. On testing for spatial correspondence between maps of human brain structure and function. NeuroImage. 2018;178:540-551.
