latex input:    mmd-beamer-header-11pt
Title:          Reduced Order Modelling for Solar System Barycentring
Date:           12 July 2017
Author:         Matthew Pitkin
Affiliation:    University of Glasgow
LaTeX xslt:     beamer
latex mode:     beamer
Theme:          m
Event:          LIGO-G1701327
latex input:    mmd-beamer-begin-doc
latex footer:   mmd-beamer-footer

<!--
% A presentation to the LVC CW group on using reduced order modelling to speed up solar
% system barycentring time delay calculations
%
% Note: comments can be included in the LaTeX file by surrounding them with html style comment
% blocks and a % sign
-->

### Background ###

The phase evolution of a continuous gravitational wave signal arriving at a detector on
Earth can be written in the form:
\\[
\phi(t) = \phi_0 + 2\pi \left(f_0\left[t+\tau(t)\right] + \frac{\dot{f}}{2}\left[t + \tau(t)\right]^2 \ldots \right),
\\]
where $t$ is the time at the detector, and $\tau(t)$ is the time delay between the signal arrival
at the detector and the signal arrival at the solar system barycentre (SSB)[^fn].

[^fn]: in, e.g., the <!-- \href{https://en.wikipedia.org/wiki/Barycentric_Dynamical_Time}{TDB}--> coordinate time standard

### Background ###

The time delay term depends on the location of the source, and the position and velocity of the detector with
respect to the SSB, and comprises:
\\[
\tau(t) = \Delta T_{\rm R}(t) + \Delta T_{\rm E} - \Delta T_{\rm S}, 
\\]
where:

 * $\Delta T_{\rm R}(t)$ is the geometric Roemer delay
 * $\Delta T_{\rm E}$ is the relativistic Einstein delay
 * $\Delta T_{\rm S}$ is the Shapiro delay

### Background ###

Shapiro delay and Einstein delay over one year for source at: $\alpha=0^{\rm h}0^{\rm m}0.0^{\rm s}$, $\delta = 0^{\circ}0'0''.0$

![][shap_ein]

[shap_ein]: images/shap_ein_delay.pdf "Shapiro and Einstein delays" height="170px"

### Background ###

Roemer delay (over one month) between Earth and SSB, and H1 and geocentre

![][roemer]

[roemer]: images/roemer_delay.pdf "Roemer delays" height="170px"

### Background ###

For long duration signals the phase needs to be tracked precisely -- time delay needs to be known
to $\mathcal{O}(\mu{\rm s})$ accuracy.

<!-- %For an observation time of a year, and a signal at 100 Hz, sky locations
%separated by $15 \mu$ rad (near the ecliptic) produce 10% mismatches, so long wide-sky-area
%searches require the SSB calculation to be repeated for many sky positions. -->

The code for the SSB calculation, e.g.<!-- \href{http://software.ligo.org/docs/lalsuite/lalpulsar/\_l\_a\_l\_barycenter\_8c\_source.html\#l00078}{\texttt{XLALBarycenterEarth()}} --> and <!-- \href{http://software.ligo.org/docs/lalsuite/lalpulsar/\_l\_a\_l\_barycenter\_8c\_source.html\#l00828}{\texttt{XLALBarycenter()}}, -->
consists of many calls to `sin` and `cos`, so if needed many times this could be a computational bottleneck.

**Is there a way to speed up the calculation _and_ maintain accuracy?**

### Reduced Order Modelling ###

<!-- \href{https://en.wikipedia.org/wiki/Model_order_reduction}{Reduced-order modelling} --> (ROM) is basically a compression technique (similar to, e.g., <!-- \href{https://en.wikipedia.org/wiki/Principal_component_analysis}{Principal
Component Analysis} -->).

 * generate a ``training set'' of signal model vectors (each with length $M$) over a required parameter space
 * use modified <!-- \href{https://en.wikipedia.org/wiki/Gram\%E2\%80\%93Schmidt\_process}{Gram-Schmidt process} --> to form a
   minimal set of orthonormal bases from the ``training set'' that satisfy some constraint:
    * a small projection error of the _current_ bases onto the remaining training data
    * a small **resdiual** when generating an interpolant from the current bases and
     comparing the interpolated models to the training set models

### Reduced Order Modelling ###

The set of $N$ _reduced_ orthonormal bases can then be used to form an interpolant:

 * find $N$ best points (interpolation nodes) in the basis with which to form the interpolant
 * straightforward linear algebra to find a $N\times N$ matrix for interpolation

Linear algebra using the interpolant, the _reduced_ basis, and values of the full model function evaluated at only
the $N$ nodes (as opposed to $M$ points) gives an approximation of the full function at all $M$ points.

### Reduced Order Modelling ###

See, e.g., Appendix A & B of [Canizares _et al_, PRD, 124005 (2013)](http://ukads.nottingham.ac.uk/abs/2013PhRvD..87l4005C) for algorithms,
and [`greedycpp`](https://bitbucket.org/sfield83/greedycpp) code for more details.

### Example analysis ###

Form a reduced bases, and interpolant for any sky position, for the SSB time delay $\tau$ spanning
1 year for H1

 * generate 2000 random sky positions $[\alpha, \delta]$ drawn uniformly over the sky
 * for each sky point calculate[^fn1] $\tau(t)$ over one year in 60 s steps (a $2000 \times 525960$ array)
   using `XLALBarycenterEarth()` and `XLALBarycenter()` 
 * form a reduced basis with the constraint that the interpolant produces residual
   time delays of $< 0.1\mu{\rm s}$
 * validation and enrichment performed to fill in any sky gaps ($36\,000$ sky points tested
   in total)

[^fn1]: excluding Shapiro delay due to cuspy nature

### Example analysis ###

Time delay (without Shapiro delay) can be reduced to 5 bases

![][bases]

[bases]: images/reduced_bases.pdf "Reduced bases" height="200px"

### Example analysis ###

From the (now $525960 \times 5$) reduced basis, $\mathbf{B}$, calculate a $5\times 5$ interpolation matrix.

Then, for **any point in the sky**, we can produce an approximation of the time delay over a year by calculating
the delay using `XLALBarycenterEarth()` and `XLALBarycenter()` at only the 5 interpolation nodes.

See additional slides for linear algebra.

### Example analysis ###

Therefore, we have 5 calls to the barycentring routines[^fn3] (plus some linear algebra) versus 525960 calls.

Using SWIG wrapped versions of the functions a single evaluation of the barycentring routines
takes[^fn4] $\sim 1.2\!\times\!10^{-3}$ ms, whilst using `numpy` to perform the required linear algebra
for interpolation takes $\sim 3.1$ ms.

 * Full calculation: $525960 \times 1.2\!\times\!10^{-3} \approx 630$ ms
 * ROM calculation: $(5 \times 1.2\!\times\!10^{-3}) + 3.1 \approx 3.1$ ms

There's roughly a 200 times speed-up!

[^fn3]: `XLALBarycenterEarth()`+`XLALBarycenter()`
[^fn4]: After initialising everything required. Setting the SWIG-wrapped `LIGOTimeGPS` structure
actually takes roughly as long as each call to the barycentering routines!

### Example analysis ###

**It can be quick, but is it accurate?**

![][validation]

[validation]: images/validation.pdf "Validation points" height="190px"


### Example analysis ###

**But, what about Shapiro delay?**

Shapiro delay depends on sky position and detector location, but two main components of its calculation are sky
position independent, so can be computed once at an initialisation stage. 

For new sky positions these components then be used to calculate the full Shapiro delay over all times steps in the basis,
which (in `python`)[^fn5] takes $\sim 27$ ms.

So, when including Shapiro delay the ROM takes $\sim 30$ ms, which is an overall speed improvement of $\sim 20$.

[^fn5]: This could probably be reduced if coded in `C`

### Example analysis ###

Example residuals for a random sky location (including Shapiro delay)

![][residuals]

[residuals]: images/residuals.pdf "Example residuals" height="170px"

### Next steps ###

 * 5 basis vectors seem to be all that is required for day $\rightarrow$ month $\rightarrow$ year time-scales
 * create `LAL` functionality to use ROM for SSB calculation
    * calculate reduced basis and interpolant matrices and nodes for successive year long time stretches
    * store matrices, nodes and sky locations of training models that produced the reduced basis
    * create initialisation function to reproduce the required reduced basis, and Shapiro delay
      components
    * create function to return a vector of times delays at a set of given GPS times (may use linear
      interpolation between basis set time steps)

### Usage ###

Although potentially not relevant for current all-sky searches, this could be useful for:

 * **long coherent follow-ups** of candidates with broad sky error regions
 * searches that **stochasticly sample the sky** (e.g. using an MCMC) at many points

Reduced order modelling (and the related _Reduced Order Quadrature_ for likelihood calculations)
is already implemented (but, yet to be written-up!) in the standard Bayesian pipeline known pulsar search
for narrow parameter searches including frequency and binary parameters.

### Future ###

Even if this has limited applicability in current CW searches it may be relevant in other cases:

 * CBC signals may need days-long templates in 3G detector era that require referencing to SSB
 * Many signals in LISA require referencing to the SSB
 * PTAs currently have quite sparse time samples, but in the future (e.g. SKA) there may be larger
   numbers of samples

### Acknowlegdments ###

Rory Smith has provided invaluable information on Reduced Order Modelling.

Initial stages of this work were performed by two undergraduate students: Stuart Doolan
and Lisa McMenanmin.

I am grateful to Scott Field and the other developers of [`greedycpp`](https://bitbucket.org/sfield83/greedycpp).

The modified version of `greedycpp` used for this analysis can be found [here](https://github.com/mattpitkin/greedycpp/tree/redordbar/).

### Additional slides: linear algebra ###

We have an $M \times N$ reduced basis matrix $\mathbf{B}$, where $M$ is the length of each base and $N$ is the number of bases.
Due to their orthonomality we know that a linear combination of bases should be able to produce an approximation to any _true_
full model. If we have $N$ nodes that are optimal points for interpolation, then we can calculate the coefficients $\vec{C}$ of the
linear superposition of bases:
<!-- \begin{align}
\vec{t} &= \vec{C}~\mathbf{b}, \nonumber \\
\vec{C} &= \vec{t}~\mathbf{b}^{-1} \nonumber
\end{align} -->
where $\vec{t}$ is the vector of model values computed at the $N$ nodes, and $\mathbf{b}$ is an $N \times N$ matrix of the reduced
basis values at the $N$ nodes.

### Additional slides: linear algebra ###

We can then form the full model vector (at all $M$ points) with
\\[
\vec{T} = \vec{C}~\mathbf{B}.
\\]

Alternatively, we can reorder this, given that $\vec{C} = \vec{t}~\mathbf{b}^{-1}$, so that
\\[
\vec{T} = \vec{t}\,(\mathbf{b}^{-1}~\mathbf{B}),
\\]
where the $\mathbf{b}^{-1}~\mathbf{B}$ part can be precomputed.

