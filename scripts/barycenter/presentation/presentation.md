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
at the detector and the signal arrival at the solar system barycentre (SSB).

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
to $\mathcal{O}(\mu{\rm s})$ accuracy. For an observation time of a year, and a signal at 100 Hz, sky locations
separated by $15 \mu$ rad (near the ecliptic) produce 10% mismatches, so long wide-sky-area
searches require the SSB calculation to be repeated for many sky positions.

The code for the SSB calculation, e.g. <!-- \href{http://software.ligo.org/docs/lalsuite/lalpulsar/\_l\_a\_l\_barycenter\_8c\_source.html\#l00078}{\texttt{XLALBarycenter()}}, -->
consists of many calls to `sin` and `cos`, so if needed many times could be a computational bottleneck.

**Is there a way to speed up the calculation _and_ maintain accuracy?**

### Reduced Order Modelling ###

Reduced order modelling (ROM) is basically a compression technique (similar to, e.g., Principal
Component Analysis).

 * generate a ``training set'' of signal model vectors (each with length $M$) over a required parameter space
 * use modified <!-- \href{https://en.wikipedia.org/wiki/Gram\%E2\%80\%93Schmidt\_process}{Gram-Schmidt process} --> to form a
   minimal set of orthonormal bases from the ``training set'' that satisfy some constraint:
    * a small projection error of the _current_ bases onto the remaining training data
    * a small **resdiual** when generating an interpolant from the current bases and
     comparing the interpolated models to the training set models

### Reduced Order Modelling ###

The set of $N$ _reduced_ orthonormal bases can then be used to form an interpolant:

 * find $N$ best points (interpolation nodes) in the bases with which to form the interpolant
 * straightforward linear algebra to find a $N\times N$ matrix for interpolation

Using the interpolant, the _reduced_ bases, and evaluating the full model function at only
the $N$ nodes (as opposed to $M$ points) gives an approximation of the full function at all
$M$ points.

See, e.g., Appendix A & B of [Canizares _et al_, PRD, 124005 (2013)](http://ukads.nottingham.ac.uk/abs/2013PhRvD..87l4005C) for algorithms,
and [`greedycpp`](https://bitbucket.org/sfield83/greedycpp) code for more details.

### Example analysis ###

Form a reduced bases, and interpolant for any sky position, for the SSB time delay $\tau$ spanning
1 year for H1

 * generate 2000 sky positions drawn uniformly over the sky
 * for each sky point calculate $\tau(t)$ over one year in 60 s steps (a $2000 \times 525960$ array)
 * form a reduced basis with the constraint that the interpolant produces residual
   time delays of $< 0.1\mu{\rm s}$
 * validation and enrichment steps performed to check for any gaps ($30\,000$ sky points tested
   in total)

### Example analysis ###

The full sky the time delay can be reduced to 5 bases

![][bases]

[bases]: images/reduced_bases.pdf "Reduced bases" height="190px"

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

