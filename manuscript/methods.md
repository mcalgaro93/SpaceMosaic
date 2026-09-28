# Methods

## Method overview

SpaceMosaic partitions a prespecified population of target cells into spatially contiguous patches. Patch construction combines an elliptical spatial model with terms that preserve variation in the design variables and favor similarity of the surrounding cellular context. The current implementation consists of a multiscale neighborhood embedding followed by an iterative, two-pass assignment procedure and a graph-based contiguity filter.

## Input data and indices

Let the complete spatial dataset be represented by

$$
\mathcal D = (\mathbf C,\mathbf S,\mathbf t,\mathbf M),
$$

where cells are indexed by $i=1,\ldots,N$ and genes by $g=1,\ldots,G$. The matrix $\mathbf C$ contains cell-by-gene expression measurements, $\mathbf S$ contains spatial coordinates, $\mathbf t$ identifies independent spatial domains, and $\mathbf M$ contains cell-level metadata and annotations. Their dimensions are

$$
\mathbf C\in\mathbb R^{N\times G},
\qquad
\mathbf S\in\mathbb R^{N\times2},
\qquad
\mathbf t=(t_1,\ldots,t_N)^\top.
$$

The entries of $\mathbf C$ may be nonnegative counts or a specified transformed expression measure. For each cell $i$, the $i$-th row of $\mathbf S$ is

$$
\mathbf s_i=(x_i,y_i)^\top.
$$

The label $t_i$ denotes the tissue section, field of view, or other independent coordinate system containing cell $i$. Cell identifiers and row order must agree across $\mathbf C$, $\mathbf S$, $\mathbf t$, and $\mathbf M$.

## Target cells and design matrix

Let

$$
\mathcal R=\{r_1,\ldots,r_n\}\subseteq\{1,\ldots,N\}
$$

be the ordered set of $n$ target cells. The index $a=1,\ldots,n$ refers to the $a$-th target cell, whose complete-data index is $r_a$. Patches are constructed from the target-cell coordinates

$$
\mathbf S_{\mathcal R}
=
\begin{pmatrix}
\mathbf s_{r_1}^\top\\
\vdots\\
\mathbf s_{r_n}^\top
\end{pmatrix}
\in\mathbb R^{n\times2}.
$$

Let $K$ denote the number of design variables and let $\ell=1,\ldots,K$ index their columns. The target-cell design matrix is

$$
\mathbf X=[x_{a\ell}]\in\mathbb R^{n\times K},
$$

with row vector $\mathbf x_a^\top$. Design variables may include distance from a biological structure, local abundance of a cell population, or another spatially resolved covariate.

For gene $g$, the expression outcome among target cells may be written as

$$
\mathbf y_g
=
(C_{r_1g},\ldots,C_{r_ng})^\top,
$$

or as the corresponding normalized or transformed outcome.

## Cellular embedding and multiscale neighborhood representation

Let

$$
\mathbf E\in\mathbb R^{N\times Q},
\qquad
\mathbf e_i^\top=\mathbf E[i,\cdot],
$$

be a $Q$-dimensional embedding of the measured state of all $N$ cells.

The embedding provides a compact representation of cell state. When PCA scores are used, correlated gene-level measurements are summarized into orthogonal components, reducing dimensionality before spatial neighborhoods are characterized.

Let $h=1,\ldots,H$ index the selected neighborhood scales and let $k_h$ be the number of neighbors at scale $h$. For cell $i$, define

$$
\mathcal N_{k_h}(i;t_i)
$$

as its $k_h$ nearest cells in Euclidean coordinate space among cells with the same spatial-domain label $t_i$. The focal cell is excluded. The neighborhood embedding at scale $h$ is

$$
\mathbf z_i^{(h)}
=
\frac{1}{k_h}
\sum_{j\in\mathcal N_{k_h}(i;t_i)}\mathbf e_j.
$$

This average replaces the state of any single neighboring cell with the typical cell state observed around $i$, reducing sensitivity to cell-level variation. Excluding the focal cell prevents its own embedding from directly determining its neighborhood representation, while restricting neighbors by $t_i$ prevents information from being shared across independent spatial domains.

The $H$ scale-specific vectors are concatenated in their supplied order:

$$
\mathbf z_i
=
\mathbf z_i^{(1)}\oplus\cdots\oplus\mathbf z_i^{(H)}
\in\mathbb R^{QH}.
$$

Concatenation retains the information from each spatial scale as a separate block rather than averaging the scales together. Small neighborhoods capture local variation, whereas larger neighborhoods provide a more stable description of the broader tissue context.

Stacking these rows gives the all-cell context matrix

$$
\mathbf Z_{\mathrm{all}}
=
(\mathbf z_1,\ldots,\mathbf z_N)^\top
\in\mathbb R^{N\times QH}.
$$

After neighborhood construction, the target-cell rows are selected:

$$
\mathbf Z
=
\mathbf Z_{\mathrm{all}}[\mathcal R,\cdot]
\in\mathbb R^{n\times QH}.
$$

Neighborhoods are computed before this restriction so that the context of a target cell can be informed by all measured cells in its spatial domain, including non-target populations. Only the rows to be partitioned are retained for patch construction.

The default neighborhood sizes are $k_1=5$ and $k_2=50$, providing complementary local and broader representations of the same position. Each spatial domain must contain at least $1+\max_h k_h$ cells.

![Figure 1. Construction of the multiscale cellular-context representation. Cell-level embeddings are averaged separately over local and broader neighborhoods and concatenated to form $\mathbf Z$. For target cells, the concatenated vectors are subsequently averaged over the 10 nearest target cells.](figures/multiscale_embedding.png)

## Preparation of the patching inputs

### Scaling of the design matrix

Each column of $\mathbf X$ is divided by its sample standard deviation to make the design variables comparable in scale:

$$
x^*_{a\ell}
=
\frac{x_{a\ell}}
{\operatorname{sd}(x_{1\ell},\ldots,x_{n\ell})},
\qquad
a=1,\ldots,n,
\quad
\ell=1,\ldots,K.
$$

The resulting matrix is denoted by $\mathbf X^*$.

This rescaling makes distances and variances in $\mathbf X$ invariant to the measurement units of its columns, so that a variable does not receive more weight solely because it has a larger numerical range. The columns are not centered because the algorithm uses within-patch differences and variances, which are unchanged by subtracting a common column mean.

### Internal smoothing of the context matrix

When $\mathbf Z$ is supplied, it is smoothed over the 10 nearest target cells to get a local average of the neighborhood embedding of each target cell $a$. Let $\mathcal N^{\mathcal R}_{10}(a)$ denote the 10 nearest target cells to target cell $a$. Then

$$
\mathbf z_a^*
=
\frac{1}{10}
\sum_{b\in\mathcal N^{\mathcal R}_{10}(a)}\mathbf z_b.
$$

This second averaging step regularizes the context representation across nearby target cells. It reduces the influence of an atypical cell or an unstable neighborhood estimate and makes the context representation reflect differences that persist locally in space.

A context-distance scale is estimated from a random sample $\mathcal J$ of

$$
m=\min(n,10000)
$$

target cells:

$$
\sigma_Z
=
\frac{1}{10m}
\sum_{a\in\mathcal J}
\sum_{b\in\mathcal N^{\mathcal R}_{10}(a)}
\left\|\mathbf z_a^*-\mathbf z_b^*\right\|_2.
$$

$\sigma_Z$ estimates the typical context difference between spatially adjacent target cells. Expressing context distances relative to this empirical scale makes their contribution invariant to a uniform rescaling of $\mathbf Z$. The cap of 10,000 sampled cells limits computation in large datasets while retaining a broad estimate of local context variability.

If $\sigma_Z=0$, it is replaced by one to avoid division by zero. The internal neighbor search uses only the supplied coordinates and does not receive a spatial-domain label. Independent coordinate systems should therefore be processed separately. When $\mathbf Z$ is supplied, the current implementation also requires at least 11 target cells.

### Spatial limits and contiguity graph

When `max_radius` is not supplied, the algorithm calculates the area $A_{\mathrm{hull}}$ of the convex hull of the target-cell coordinates and sets

$$
r_{\max}
=
3\left(\frac{A_{\mathrm{hull}}}{P\pi}\right)^{1/2},
$$

where $P$ is the requested number of patches.

The quantity $A_{\mathrm{hull}}/P$ is the average area available to one patch. Approximating this area by a circle gives the characteristic radius $\{A_{\mathrm{hull}}/(P\pi)\}^{1/2}$; the factor of three provides tolerance around this scale while preventing assignments to spatially remote centroids.

Before the iterative procedure, a symmetric nearest-neighbor contiguity graph $G_c=(\mathcal V,\mathcal E_c)$ is constructed on the target cells, using

$$
k_c=\min(10,n-1).
$$

An undirected edge is included when either cell is among the $k_c$ nearest neighbors of the other:

$$
(a,b)\in\mathcal E_c
\iff
b\in\mathcal N_{k_c}(a)
\ \text{or}\
a\in\mathcal N_{k_c}(b).
$$

Symmetrization creates an adjacency whenever either cell identifies the other as locally close. This is less sensitive to asymmetric nearest-neighbor relations caused by spatially varying cell density and provides a sparse representation of local spatial continuity.

## Patch initialization

Under the default initialization, k-means is applied to $\mathbf S_{\mathcal R}$ using $P$ centers, five random starts, and at most 50 k-means iterations. This produces the initial patch labels

$$
c_a^{(0)}\in\{1,\ldots,P\},
\qquad a=1,\ldots,n.
$$

K-means provides a spatially compact starting partition without using the design or context variables. Multiple random starts reduce sensitivity to a poor initial set of centroids.

Because k-means is stochastic, a random seed is required for exact reproduction.

The optional `gradient_ellipse` initialization is available when $K=1$. It first obtains k-means labels and then estimates a local spatial gradient of the scaled design variable $X^*$ for each target cell. Let $\mathbf A_a$ contain the coordinate displacements from cell $a$ to its selected neighbors and let $\mathbf d_a$ contain the corresponding differences in $X^*$. The local gradient is

$$
\widehat{\mathbf g}_a
=
(\mathbf A_a^\top\mathbf A_a+10^{-8}\mathbf I)^{-1}
\mathbf A_a^\top\mathbf d_a.
$$

Within each initial patch, the major axis of the covariance ellipse is oriented along the mean estimated gradient. The default gradient neighborhood size is 30 and the target initial eigenvalue ratio is 4, capped by `max_elongation`.

This optional initialization allows an elongated starting patch to follow a smooth spatial gradient in the design variable, rather than forcing the first iteration to begin from approximately isotropic regions. The small ridge term $10^{-8}\mathbf I$ stabilizes the local gradient estimate when neighbor displacements are nearly collinear.

## Iterative patch construction

Let $u=1,\ldots,U$ index the outer iterations, with $U=15$ by default. Empty patch labels are removed from the active set. Each iteration performs the following steps.

### Ellipse estimation

For active patch $p$, let

$$
\mathcal I_p^{(u-1)}
=
\{a:c_a^{(u-1)}=p\},
\qquad
n_p^{(u-1)}
=
\left|\mathcal I_p^{(u-1)}\right|.
$$

To simplify notation within one iteration, write these as $\mathcal I_p$ and $n_p$. The spatial centroid is

$$
\boldsymbol\mu_p
=
\frac{1}{n_p}
\sum_{a\in\mathcal I_p}\mathbf s_{r_a}.
$$

For $n_p\geq3$, $\boldsymbol\Sigma_p$ is the sample covariance of the coordinates in the patch. For $n_p<3$, it is set to

$$
\boldsymbol\Sigma_p
=
v_{\mathrm{global}}\mathbf I_2,
\qquad
v_{\mathrm{global}}
=
\frac{1}{2}
\sum_{d=1}^2
\operatorname{Var}
\left\{
s_{r_a,d}:a=1,\ldots,n
\right\}.
$$

The centroid and covariance summarize the location, orientation, and spatial extent of each patch. The global isotropic fallback permits these quantities to remain defined when too few cells are available to estimate a two-dimensional covariance.

Let $\lambda_{p1}\geq\lambda_{p2}$ be the eigenvalues of $\boldsymbol\Sigma_p$. Both eigenvalues are truncated below at $10^{-10}$. If

$$
\frac{\lambda_{p1}}{\lambda_{p2}}>L_{\max},
$$

the smaller eigenvalue is replaced by

$$
\lambda_{p2}
\leftarrow
\frac{\lambda_{p1}}{L_{\max}},
$$

where $L_{\max}$ is `max_elongation`. The regularized covariance, its inverse, and its log determinant define the elliptical spatial model for patch $p$.

The eigenvalue floor prevents singular covariance matrices, whereas the elongation constraint prevents an ellipse from collapsing into an arbitrarily thin spatial structure. Both operations stabilize Mahalanobis distances and the inversion of $\boldsymbol\Sigma_p$.

### Joint assignment using spatial X Z and hunger terms

Patch-level means are calculated from the assignments entering the step:

$$
\overline{\mathbf x}_p^*
=
\frac{1}{n_p}
\sum_{a\in\mathcal I_p}\mathbf x_a^*,
\qquad
\overline{\mathbf z}_p^*
=
\frac{1}{n_p}
\sum_{a\in\mathcal I_p}\mathbf z_a^*.
$$

These means provide the current reference values of the design variables and cellular context within each patch. Candidate cells are assessed relative to these patch-specific references rather than to a global average.

The total within-patch variation of the design variables is

$$
V_p
=
n_p
\sum_{\ell=1}^K
\operatorname{Var}
\left\{
x^*_{a\ell}:a\in\mathcal I_p
\right\}.
$$

$V_p$ measures the total spread of the scaled design variables within patch $p$, with multiplication by $n_p$ accounting for patch size. Because the variables have been rescaled, their contributions can be summed without one column dominating through its original units.

$V_p$ is set to zero for a patch containing at most one cell. Let

$$
\overline V
=
\frac{1}{|\mathcal P|}
\sum_{p\in\mathcal P}V_p
$$

be the mean across the set $\mathcal P$ of active patches, and let $w$ be `hunger_weight`. The normalized hunger term is

$$
h'_p
=
\frac{1}{(1-w)\overline V+wV_p},
\qquad
h_p
=
\frac{h'_p}{\sum_{q\in\mathcal P}h'_q}.
$$

Patches with smaller current design variation receive larger hunger values and are therefore more competitive for additional cells. The parameter $w$ interpolates between equal hunger across patches and stronger adaptation to their current design variation.

For computational efficiency, target cell $a$ is evaluated only against the $C$ nearest active patch centroids, where

$$
C
=
\min\{\texttt{n\_candidates},|\mathcal P|\}.
$$

Restricting the comparison to nearby centroids reduces computation and avoids evaluating patches with implausibly distant spatial centers. The subsequent radius constraints still determine whether any candidate is spatially admissible.

Let $\mathcal C_a\subseteq\mathcal P$ denote this candidate set. For candidate patch $p$, the squared Mahalanobis distance is

$$
m_{ap}^2
=
(\mathbf s_{r_a}-\boldsymbol\mu_p)^\top
\boldsymbol\Sigma_p^{-1}
(\mathbf s_{r_a}-\boldsymbol\mu_p).
$$

The spatial log score is

$$
G_{ap}
=
-\frac{1}{2}m_{ap}^2
-\frac{1}{2}\log|\boldsymbol\Sigma_p|.
$$

The Mahalanobis term measures distance in the orientation and scale of the patch ellipse, so displacement along a long patch axis is treated differently from displacement across its narrow axis. The log-determinant term accounts for ellipse size and prevents diffuse patches from receiving high spatial scores solely because of their large covariance.

The score is set to $-\infty$ if

$$
\|\mathbf s_{r_a}-\boldsymbol\mu_p\|_2>r_{\max}
$$

or

$$
m_{ap}>r_M,
$$

where $r_M$ is `mahal_radius`. Otherwise, the complete assignment score is

$$
Q_{ap}
=
G_{ap}
+\beta\left\|\mathbf x_a^*-\overline{\mathbf x}_p^*\right\|_2^2
-\alpha
\frac{
\left\|\mathbf z_a^*-\overline{\mathbf z}_p^*\right\|_2^2
}{\sigma_Z^2}
+\log h_p.
$$

The positive $\mathbf X$ term favors assignments that preserve within-patch contrast in the design variables. The negative $\mathbf Z$ term favors coherence of the local cellular context, and $\sigma_Z^2$ expresses this mismatch relative to typical local context variation. The coefficients $\beta$ and $\alpha$ control the relative strengths of these two contributions, while $\log h_p$ discourages low-variation patches from becoming inactive.

![Figure 2. Geometry and scoring of cell-to-patch assignment. **A**, the patch covariance defines elliptical Mahalanobis contours, while a separate Euclidean radius imposes an absolute spatial cap. **B**, two cells at the same Euclidean distance can have different Mahalanobis distances because displacement along the major and minor axes is scaled differently. **C**, the final score combines signed contributions from spatial fit, design-variable contrast, cellular-context coherence, and patch hunger. Values are illustrative.](figures/patch_geometry.png)

When $\mathbf Z$ is absent, the context-penalty term is omitted. The provisional assignment is

$$
c_a^{\mathrm{joint}}
=
\underset{p\in\mathcal C_a}{\operatorname{arg\,max}}
\ Q_{ap}.
$$

A cell is provisionally unassigned when all candidate scores are $-\infty$.

### Ellipse refitting and spatial-only assignment

The patch ellipses are re-estimated from the provisional joint assignments. By default, all assigned cells receive equal weight.

If `x_weighted_ellipse_second_pass=TRUE` and $K=1$, cell $a$ in patch $p$ receives weight

$$
\omega_a
=
1+\gamma
\min\left\{
\frac{|x_a^*-\overline x_p^*|}
{\operatorname{MAD}_p(x^*)+10^{-8}},
w_{\max}
\right\}.
$$

A sample standard deviation replaces the median absolute deviation if the latter is not positive, and unit weights are used if neither scale is positive. Weighted centroids and covariances are then calculated before applying the same eigenvalue regularization.

This option gives greater geometric influence to cells that extend the range of $X$ within a patch. The standardized deviation and the cap $w_{\max}$ limit the influence of isolated extreme values.

Using the refitted ellipses, every target cell is reassigned with the spatial score alone:

$$
c_a^{\mathrm{spatial}}
=
\underset{p\in\mathcal C_a}{\operatorname{arg\,max}}
\ G_{ap}.
$$

Thus, $\mathbf X$, $\mathbf Z$, and hunger affect the geometry estimated in the first pass, whereas the labels passed to the contiguity step are produced by the second, spatial-only pass.

The second pass therefore converts the jointly informed patch geometry into a spatial partition without applying the $\mathbf X$, $\mathbf Z$, or hunger terms a second time.

### Contiguity filtering

For each spatial-only patch $p$, the algorithm extracts the subgraph of $G_c$ induced by the cells assigned to $p$. Let $\mathcal L_p$ denote its largest connected component. The iteration-level assignment is

$$
c_a^{(u)}
=
\begin{cases}
p, & a\in\mathcal L_p,\\
\mathrm{NA}, & a\notin\mathcal L_p.
\end{cases}
$$

Single-cell patches are retained. Cells outside the largest component are marked as unassigned and are not automatically restored in the same iteration.

Retaining only the largest connected component prevents a single patch label from describing multiple spatially separated islands. Marking disconnected cells as unassigned is preferable to forcing them into a component unsupported by the local adjacency graph.

![Figure 3. Initialization and iterative refinement of the patch partition. **A**, all target cells enter the procedure together, with no preliminary hotspot selection. **B**, spatial k-means provides a complete initial partition from which patch centroids and covariance ellipses are estimated. **C**, the joint score can reassign boundary cells, shown with black outlines, after which the ellipses are refitted. **D**, a second spatial-only assignment is followed by the contiguity filter; cells belonging to a disconnected island of a patch label are set to `NA`. The refinement cycle is repeated for the requested number of iterations. Synthetic data and the disconnected island are illustrative.](figures/patch_iteration.png)

## Iteration logs and final output

Ellipse estimation, joint assignment, ellipse refitting, spatial-only assignment, and contiguity filtering are repeated for the fixed number $U$ of outer iterations. The current implementation does not use an early-stopping criterion.

The final patch vector is

$$
\mathbf c
=
(c_1^{(U)},\ldots,c_n^{(U)})^\top,
$$

with one patch label or `NA` for each target cell.

When `log_iters=TRUE`, the algorithm also returns the $n\times U$ membership history and, for every active patch and iteration, the within-patch sum of squared centered design values:

$$
SS_p^{(u)}
=
\sum_{a:c_a^{(u)}=p}
\left\|
\mathbf x_a^*-\overline{\mathbf x}_p^{*(u)}
\right\|_2^2.
$$

## Patch diagnostics

The diagnostic function summarizes the final assignment without modifying it. For each nonempty patch, it reports the number of assigned cells and, for each column of the user-supplied $\mathbf X$, the number of finite observations and their sample standard deviation.

Spatial connectivity is evaluated from symmetric k-nearest-neighbor graphs. For patch $p$ and neighbor count $k$, let $\mathcal L_p(k)$ be the largest connected component of the patch-specific graph. The connectivity curve is

$$
F_p(k)
=
\frac{|\mathcal L_p(k)|}{n_p}.
$$

The strict component fraction is $F_p(k_{\mathrm{strict}})$, with

$$
k_{\mathrm{strict}}=\min(5,k)
$$

by default. The minimum connectivity $k$ is the smallest evaluated value for which $F_p(k)=1$. It is zero for a single-cell patch and `NA` if full connectivity is not reached.

If at least two iterations were logged, membership stability for final patch $p$ is the fraction of its final cells that already carried label $p$ in iteration $U-1$:

$$
S_p
=
\frac{
|\{a:c_a^{(U)}=p\ \land\ c_a^{(U-1)}=p\}|
}{
|\{a:c_a^{(U)}=p\}|
}.
$$

Global output summaries include the number and fraction of assigned and unassigned target cells and the number of nonempty patches.

## Default implementation settings

The following values are the defaults in the current implementation. They describe software behavior and should not be interpreted as universally optimal settings.

| Quantity | Implementation argument | Default |
|---|---|---:|
| Neighborhood scales | `ks` | `c(5, 50)` |
| Internal smoothing neighbors | internal | 10 |
| Context penalty | `alpha` | 0.5 |
| Design diversity weight | `beta` | 1 |
| Hunger interpolation | `hunger_weight` | 0.5 |
| Maximum covariance elongation | `max_elongation` | 4 |
| Euclidean radius | `max_radius` | Automatic from convex-hull area |
| Mahalanobis radius | `mahal_radius` | 3 |
| Candidate patches per cell | `n_candidates` | 20, capped at the number of active patches |
| Outer iterations | `n_iters` | 15 |
| Initialization | `init_method` | `"kmeans"` |
| Gradient neighborhood size | `init_gradient_k` | 30 |
| Initial gradient elongation | `init_gradient_elongation` | 4 |
| X-weighted second ellipse fit | `x_weighted_ellipse_second_pass` | `FALSE` |
| X-weighting strength | `x_ellipse_gamma` | 1 |
| X-weighting cap | `x_ellipse_wmax` | 3 |
| Iteration logging | `log_iters` | `TRUE` |

## Reproducibility and implementation scope

The random seed affects k-means initialization and, when $n>10000$, the sample used to estimate $\sigma_Z$. Reproducible analyses must therefore report the seed, all nondefault arguments, the target-cell definition, the construction of $\mathbf E$ and $\mathbf X$, and the spatial units processed together.

The functions described here assume that segmentation, cell-level quantification, expression preprocessing, cell annotation, and construction of $\mathbf E$ and $\mathbf X$ have already been completed. The patching algorithm itself begins with aligned target-cell coordinates, $\mathbf X$, and optionally $\mathbf Z$.

`embedCellNeighborhoods()` can restrict the initial neighborhood graph by tissue or field of view. In contrast, the internal smoothing and contiguity graphs constructed by `getPatches()` do not receive the spatial-domain vector $\mathbf t$. Each independent coordinate system should therefore be passed to `getPatches()` separately.
