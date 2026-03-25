Comprehensive Digital Communications Study Guide (Midterm 1)
Section 1 — Equation Reference Sheet
Wireless Digital Communication Equations (Up-Down Converter Layer)

Passband/Baseband Signal Representation

    x(t)=xI​(t)+jxQ​(t)

        Explanation: The complex baseband signal x(t) decomposed into its real in-phase component xI​(t) and imaginary quadrature component xQ​(t).

        Application: Used to represent baseband information compactly before up-conversion. 

I/Q Up Converter (Time Domain)

    xc​(t)=I/Q↑(x(t))≜ℜ[x(t)2​ej2πfc​t]

    xc​(t)=xI​(t)2​cos(2πfc​t)−xQ​(t)2​sin(2πfc​t)

        Explanation: Converts a baseband signal x(t) to a real-valued passband signal xc​(t) centered at the carrier frequency fc​. The 2​ factor ensures energy is preserved. 

        Application: Hardware implementation of the transmitter's I/Q up converter to prepare a baseband signal for efficient RF antenna transmission. 

I/Q Up Converter (Frequency Domain)

    x~c​(f)=2​x~(f−fc​)+x~∗(−f−fc​)​

        Explanation: The Fourier transform of the passband signal, showing that the baseband spectrum x~(f) is shifted to both +fc​ and −fc​. 

        Application: Used for analyzing the spectral bandwidth and frequency content of the transmitted passband signal. 

I/Q Down Converter (Time Domain)

    x(t)=I/Q↓(xc​(t))≜(xc​(t)2​e−j2πfc​t)∗hL​(t)

        Explanation: Recovers the complex baseband signal x(t) by multiplying the received passband signal xc​(t) by a complex exponential at −fc​ and passing it through a lowpass filter hL​(t). 

        Application: Hardware implementation of the receiver's I/Q down converter. 

I/Q Down Converter (Frequency Domain)

    x~(f)=2​x~c​(f+fc​)h~L​(f)

        Explanation: The Fourier transform of the down-converted signal, showing the passband spectrum shifted left by fc​ and filtered by the ideal lowpass filter h~L​(f)=rectb​(f). 

        Application: Filtering out the high-frequency double-carrier component (at −2fc​) generated during down-conversion. 

Equivalent Baseband Filter (Time & Frequency Domains)

    h~(f)=h~c​(f+fc​)h~L​(f)=2​1​I/Q↓(h~c​(f))

    h(t)=(hc​(t)e−j2πfc​t)∗hL​(t)=2​1​I/Q↓(hc​(t))

        Explanation: Relates the physical passband channel filter hc​(t) to the equivalent baseband filter h(t), which fully captures the channel's effect on the baseband signal such that y(t)=x(t)∗h(t). 

        Application: Simplifying passband channel analysis by allowing all signal processing to be modeled strictly at baseband. 

AWGN Power Spectral Density (Passband & Baseband)

    r~nc​​(f)=2N0​​

    r~n​(f)=N0​rectb​(f)

        Explanation: The passband noise nc​(t) has a flat two-sided PSD of N0​/2. After passing through the I/Q down converter's lowpass filter of bandwidth b, the equivalent baseband noise n(t) has a PSD of N0​ restricted to the filter bandwidth. 

        Application: Calculating the total noise power or variance entering the digital demodulator. 

Time Varying Multi-Path Baseband Channel Model

    y(t)=∑i​ai​(t)e−j2πfc​τi​(t)x(t−τi​(t))

        Explanation: The baseband received signal y(t) resulting from a transmission over multiple paths, each with a time-varying complex gain ai​(t) and time-varying delay τi​(t). 

        Application: General modeling of fading environments, specifically frequency-selective and fast fading wireless channels. 

Doppler Frequency

    fd,i​(t)=−fc​τi′​(t)=−ri′​(t)cfc​​

        Explanation: The frequency shift induced by the relative motion of the transmitter, receiver, or scatterer, where ri′​(t) is the rate of change of the path length. 

        Application: Modeling the phase rotations in the baseband equivalent channel h(t) for mobile users. 

General Probability Theory and Statistics Equations

Gaussian PDF (Scalar, Real)

    fx​(x)=2πσx2​​1​e−2σx2​(x−μx​)2​

        Explanation: The probability density function of a real-valued Gaussian random variable with mean μx​ and variance σx2​. 

        Application: Modeling thermal noise variations in the real or imaginary branch of a receiver individually. 

Circularly Symmetric Complex Gaussian (CSCG) PDF

    fx​(x)=πσ21​e−σ2∣x−μx​∣2​

        Explanation: The PDF of a complex Gaussian random variable x=xI​+jxQ​ where the real and imaginary parts are independent and identically distributed, and σ2 is the total variance. 

        Application: Evaluating the distribution of equivalent baseband noise samples n(t) at the output of the I/Q down converter. 

Gaussian Random Vector PDF (Real)

    fx​(x)=det(2πΣx​)​1​e−21​(x−μ​x​)TΣx−1​(x−μ​x​)

        Explanation: The joint PDF of a real-valued Gaussian random vector x, defined entirely by its mean vector μ​x​ and covariance matrix Σx​. 

        Application: Joint probability analysis of multiple correlated noise samples or channel states. 

CSCG Random Vector PDF

    fx​(x)=det(πΣx​)1​e−(x−μ​x​)†Σx−1​(x−μ​x​)

        Explanation: The joint PDF of a complex-valued CSCG random vector, where † denotes the conjugate transpose (Hermitian) operator. 

        Application: Optimal multi-sample detection and MAP/ML decoding where the noise across samples is complex and potentially correlated. 

Covariance Matrix Definition

    Σx​≜E[(x−μ​x​)(x−μ​x​)†]

        Explanation: The matrix representing the auto-covariance and cross-covariance of all elements in the random vector x. 

        Application: Structuring the covariance matrix needed for multivariate Gaussian density calculations. 

Section 2 — Practice Exam: Full Worked Solutions
Problem 1

Problem Statement
Find xI​(t), xQ​(t) of the baseband signal x(t) that produces the following passband signal xc​(t). 

    (5/10 pts) xc​(t)=8​m(t)cos(2πfc​t+4π​) where m(t) is a real-valued baseband message signal. 

    (5/10 pts) xc​(t)=4cos(2π(fc​+fΔ​)t) where fΔ​ is a constant frequency offset (fΔ​≪fc​). 

Solution Strategy
The general representation of a passband signal is xc​(t)=ℜ[2​x(t)ej2πfc​t]. To find the in-phase and quadrature components, we must manipulate the given xc​(t) expression into the standard complex exponential form, isolate the complex baseband envelope x(t), and extract its real and imaginary parts, as x(t)=xI​(t)+jxQ​(t). 

Relevant Equations

    xc​(t)=ℜ[2​x(t)ej2πfc​t] 

    cos(θ)=ℜ[ejθ] 

    x(t)=xI​(t)+jxQ​(t) 

    Euler's Formula: ejθ=cos(θ)+jsin(θ) 

Step-by-Step Solution (Part 1)

    Rewrite the cosine function in complex exponential form: xc​(t)=8​m(t)ℜ[ej(2πfc​t+4π​)]. 

    Because 8​m(t) is real, we can bring it inside the real-part operator: xc​(t)=ℜ[8​m(t)ej4π​ej2πfc​t]. 

    Equate the argument to the standard form ℜ[2​x(t)ej2πfc​t] to extract x(t):

    2​x(t)=8​m(t)ej4π​⟹x(t)=2​8​​m(t)ej4π​=2m(t)ej4π​. 

    Apply Euler's formula to expand the exponential:

    x(t)=2m(t)(cos(4π​)+jsin(4π​))=2m(t)(2​1​+j2​1​)=2​m(t)+j2​m(t). 

    Identify the real and imaginary parts: xI​(t)=2​m(t) and xQ​(t)=2​m(t). 

Final Answer (Part 1)

xI​(t)=2​m(t), xQ​(t)=2​m(t) 

Step-by-Step Solution (Part 2)

    Rewrite the passband expression using the complex exponential: xc​(t)=ℜ[4ej2π(fc​+fΔ​)t]. 

    Separate the carrier frequency fc​ from the offset fΔ​: xc​(t)=ℜ[4ej2πfΔ​tej2πfc​t]. 

    Equate to the standard form:

    2​x(t)=4ej2πfΔ​t⟹x(t)=2​4​ej2πfΔ​t=22​ej2πfΔ​t. 

    Apply Euler's formula:

    x(t)=22​(cos(2πfΔ​t)+jsin(2πfΔ​t)). 

    Isolate the real and imaginary components:

    xI​(t)=22​cos(2πfΔ​t) and xQ​(t)=22​sin(2πfΔ​t). 

Final Answer (Part 2)

xI​(t)=22​cos(2πfΔ​t), xQ​(t)=22​sin(2πfΔ​t) 
Problem 2

Problem Statement
A zero-mean, white, real, stationary, and Gaussian noise process w(t) with rw​(τ)=2N0​​δ(τ) is the input to a passband filter h~c​(f) and the output is nc​(t). The passband noise nc​(t) is the input to an I/Q down converter, which produces output n(t)=nI​(t)+jnQ​(t). We know fc​=1000Hz and the filter has an asymmetric frequency response: 
h~c​(f)=⎩⎨⎧​2​4−500∣f∣​​0​500≤∣f∣≤10001000<∣f∣[cites​tart]≤2000elsewhere​ 

    (2/10 pts) Calculate r~nc​​(f), the power spectral density of nc​(t). Draw the picture. 

    (2/10 pts) Calculate r~n​(f), the complex baseband power spectral density of n(t). Draw the picture. 

Solution Strategy
The passband filter h~c​(f) shapes the flat AWGN spectrum into a specific PSD according to standard linear system filtering rules. Specifically, r~nc​​(f)=r~w​(f)∣h~c​(f)∣2. To find the baseband PSD r~n​(f), apply the frequency-domain I/Q down converter transformation n~(f)=2​n~c​(f+fc​)h~L​(f), which requires shifting the passband spectrum left by fc​, right by fc​ (due to the complex conjugate), summing them, and scaling appropriately. 

Relevant Equations

    r~w​(f)=2N0​​ 

    r~nc​​(f)=r~w​(f)∣h~c​(f)∣2 

    r~n​(f)=2r~nc​​(f+fc​)∣h~L​(f)∣2 

Step-by-Step Solution

    First, find the PSD of the input white noise w(t) by taking the Fourier transform of the autocorrelation rw​(τ)=2N0​​δ(τ), which gives r~w​(f)=2N0​​. 

    Apply the passband filter to find r~nc​​(f). Multiply r~w​(f) by ∣h~c​(f)∣2. 
    Square the piecewise terms of h~c​(f):
    For 500≤∣f∣[cites​tart]≤1000: ∣2​∣2=2. 
    For 1000<∣f∣[cites​tart]≤2000: ∣4−500∣f∣​​∣2=4−500∣f∣​. 

    Construct r~nc​​(f):
    r~nc​​(f)=2N0​​×⎩⎨⎧​24−500∣f∣​0​500≤∣f∣≤10001000<∣f∣[cites​tart]≤2000elsewhere​ 
    r~nc​​(f)=⎩⎨⎧​N0​N0​(2−1000∣f∣​)0​500≤∣f∣≤10001000<∣f∣[cites​tart]≤2000elsewhere​ 

    To find the baseband PSD r~n​(f), apply the down-conversion relation: r~n​(f)=2r~nc​​(f+fc​)∣h~L​(f)∣2. Since fc​=1000, we shift the positive frequency part of r~nc​​(f) to the left by 1000 Hz, and scale by 2.  (The negative frequency part shifts out of the baseband filter bounds). 
    Shift 500≤f≤1000 left by 1000: baseband range −500≤f≤0. The value is 2×N0​=2N0​. 
    Shift 1000<f≤2000 left by 1000: baseband range 0<f≤1000. Substitute fnew​=f−1000⟹f=fnew​+1000 into N0​(2−1000f​):
    2×N0​(2−1000f+1000​)=2N0​(2−1000f​−1)=2N0​(1−1000f​). 

    Assemble r~n​(f):

    r~n​(f)=⎩⎨⎧​2N0​2N0​(1−1000f​)0​−500≤f≤00<f≤1000elsewhere​ 

Final Answer
r~nc​​(f)=⎩⎨⎧​N0​N0​(2−1000∣f∣​)0​500≤∣f∣≤10001000<∣f∣[cites​tart]≤2000else​ r~n​(f)=⎩⎨⎧​2N0​2N0​(1−1000f​)0​−500≤f≤00<f≤1000else​ 
(The requested plots are graphical step functions/ramps mapping directly to these piecewise definitions).
Additional Topics Not Covered by Practice Exam

    Equivalent Baseband Filter Calculation: The practice exam primarily covered finding signals and PSDs, but it did not test solving for the time-domain equivalent baseband filter h(t). Expect a problem where you are given a passband filter impulse response hc​(t)=δ(t)+aδ(t−τ1​) and asked to compute h(t)=(hc​(t)e−j2πfc​t)∗hL​(t). 

    Time Varying Multipath Doppler Derivation: Expect questions involving geometric scattering. For example, given a receiver moving at velocity v, you may need to compute the specific Doppler frequency fd,i​=−fc​τi′​(t0​)=λv​cos(αi​) given the arrival angle αi​ of an electromagnetic wave from a scatterer. 

    MAP/ML Sequence Detection & Sufficient Statistics: The midterm spans Part 1 in entirety. You may be asked to prove whether a specific sample vector forms a sufficient statistic for detecting an underlying symbol, using likelihood functions and continuous/discrete Gaussian pdf integrations. 

Section 3 — Key Concepts & Definitions
General Communication Architecture

    Shannon's Information Theory: The theoretical foundation (1948) establishing that reliable communication is possible with finite SNR whenever the data rate is less than the channel capacity. 

    Layered Physical System Partition: Optimal partitioning of physical layer transmitters and receivers includes the Error Control Codec Layer, Modulator/Demodulator (Modem) Layer, and the I/Q Up/Down Converter Layer. 

    Why Passband?: Antennas efficiently radiate signals whose wavelengths are comparable to their physical size (λ=fc​). Baseband signals (near zero frequency) would require impractically large antennas. Up-converting to high carrier frequencies (fc​) minimizes antenna scale. 

I/Q Conversion Characteristics

    Bandwidth Usage: Real-valued baseband signals waste half their spectral bandwidth because their frequency domain representation x~(f) has complex conjugate symmetry (x~(f)=x~∗(−f)). Complex baseband signals x(t)∈C eliminate this symmetry, doubling spectral efficiency in the same physical bandwidth. 

    Linearity: The I/Q up converter is not a linear system with respect to complex inputs because it takes the real part of the modulated signal, effectively throwing away imaginary scaling components. Conversely, the I/Q down converter is a linear system, enabling standard LTI filter analysis. 

Stochastic Models and Noise

    AWGN: Additive White Gaussian Noise describes standard thermal circuit interference. "White" implies constant power density across all frequencies, theoretically yielding infinite total power until restricted by a receiver's lowpass/passband filters. 

    CSCG Random Vectors: "Circularly Symmetric" means the complex random vector x has a distribution invariant to phase rotations. Consequently, its real and imaginary parts possess identical distributions and their cross-covariance matrix ΣxI​xQ​​ is antisymmetric. If the components' correlation is zero, the real and imaginary components are strictly independent. 

Wireless Fading Classifications

    Multipath Propagation: Electromagnetic waves bounce off scatterers resulting in multiple transmission paths, each experiencing unique attenuation, phase shift, and delay. 

    Frequency Flat Fading: Occurs when multiple path delays are extremely similar (relative to the signal bandwidth). The entire frequency band is multiplied by a single time-varying complex coefficient h(t). 

    Frequency Selective Fading: Occurs when path delays vary significantly. Different frequency components of the signal experience different phase and amplitude distortions. Modeled via convolution with h(t)