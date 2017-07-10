latex input:    mmd-beamer-header-11pt
Title:          Reduced Order Modelling for Solar System Barycentring
Date:           10 July 2017
Author:         Matthew Pitkin
Affiliation:    University of Glasgow
LaTeX xslt:     beamer
latex mode:     beamer
Theme:          m
Event:          LVC CW group
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

For long duration signals, where the phase needs to be tracked precisely, this time delay needs to be known
to $\mathcal{O}(\mu{\rm s})$ accuracy. For an observation time of a year, and a signal at 100 Hz, sky locations
separated by $15 \mu$ rad (near the ecliptic) produce 10% mismatches between signals, so long wide sky area
searches require the SSB calculation to be repeated for many sky positions.

The code for the SSB calculated, see e.g. <!-- \href{http://software.ligo.org/docs/lalsuite/lalpulsar/\_l\_a\_l\_barycenter\_8c\_source.html\#l00078}{\texttt{XLALBarycenter()}}, -->
consists of many calls to `sin` and `cos`, so if needed many times could be a computational bottleneck. Is there a way to speed
up the calculation, but maintain accuracy?

### Reduced Order Modelling ###

Reduced order modelling (ROM) is basically a compression technique (similar to, e.g., Principal
Component Analysis). It essentially

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

