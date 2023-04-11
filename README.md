# Aging_NIRS__NVC-FC
Code used in the publication: "Impaired neurovascular coupling and increased functional connectivity in the frontal cortex predict age-related cognitive dysfunction"

Software requirement:
- MATLAB 2018a or newer
- MATLAB Wavelet Toolbox 
- Brain Connectivity Toolbox (https://sites.google.com/site/bctnet/) 
- AnalyzIR Toolbox (https://github.com/huppertt/nirs-toolbox)

Measurements were carried out in the Translational Geroscience Laboratory of Oklahoma University (Center for Healthy Brain Aging), in a quiet and dimly lit room. All participants performed cognitive tests in an uninterrupted environment using the touchscreen 10.5” iOS tablet device running the CANTAB application. The study protocol was approved by the Institutional Review Board of the University of Oklahoma Health Sciences Center and was conducted in compliance with the Helsinki Declaration.

Functional NIRS measurements were performed using a NIRScout platform (NIRx Medical Technologies LLC, NY, USA). The system was equipped with 16 sources (F3, AF7, AF3, Fz, Fpz, AF4, F4, AF8, FC6, C4, FC2, CP2, FC1, CP1, C3, FC5) emitting light at two different wavelengths (760 and 850nm) and 16 photodetectors (F5, F1, Fp1, AFz, F2, Fp2, F6, AFF6h, C6, CC4, CP4, C2, C1, FC3, CP3, C5) defining 48 channels.
Detailed definition of channels: Channel 1: F3-F5, Channel 2: F3-F1, Channel 3: F3-FC3, Channel 4: AF7-F5, Channel 5: AF7-Fp1, Channel 6: AF3-F5, Channel 7: AF3-F1, Channel 8: AF3-Fp1, Channel 9: AF3-AFz, Channel 10: Fz-F1, Channel 11: Fz-AFz, Channel 12: Fz-F2, Channel 13: Fpz-Fp1, Channel 14: Fpz-AFz, Channel 15: Fpz-Fp2, Channel 16: AF4-AFz, Channel 17: AF4-F2, Channel 18: AF4-Fp2, Channel 19: AF4-F6, Channel 20: F4-F2, Channel 21: F4-F6, Channel 22: F4-CC4, Channel 23: AF8-Fp2, Channel 24: AF8-F6, Channel 25: FC6-F6, Channel 26: FC6-C6, Channel 27: FC6-CC4, Channel 28: C4-C6, Channel 29: C4-CC4, Channel 30: C4-CP4, Channel 31: C4-C2, Channel 32: FC2-F2, Channel 33: FC2-CC4, Channel 34: FC2-C2, Channel 35: CP2-CP4, Channel 36: CP2-C2, Channel 37: FC1-F1, Channel 38: FC1-C1, Channel 39: FC1-FC3, Channel 40: CP1-C1, Channel 41: CP1-CP3, Channel 42: C3-C1, Channel 43: C3-FC3, Channel 44: C3-CP3, Channel 45: C3-C5, Channel 46: FC5-F5, Channel 47: FC5-FC3, Channel 48: FC5-C5.
Functional NIRS examinations were performed using the NIRScout platform (NIRx Medical Technologies LLC, NY, USA) equipped with 16 light source and 16 detector optodes. A128-port Easycap headcap (Easycap GmbH, Woerthsee-Etterschlag, Germany) was positioned over the head to cover the area of the international 10-10 system. The line between Fpz and Iz ports on the headcap was aligned with the sagittal plane of the head, and the optode in the Fpz position of the cap was aligned with Fpz on the subject. The cap was set up with custom spacers that limit the variability of distance between optodes to average source-detector separation of 3 cm. The placement of optodes covered the prefrontal cortex, dorsolateral prefrontal cortex, and also included the medial motor cortex. Sufficient coverage of these regions was determined by the projection of channel position to the cortical surface within the Montreal Neurological Institute coordinate space [8, 9].
To evoke NVC responses during fNIRS recording, the cognitive n-back paradigm was used. The protocol consisted of n-back runs of various difficulties, which have been well-established in cognitive neuroscience research. Prior to an assessment, we had explained the n-back paradigm in detail to participants and confirmed if they understood the instructions. The n-back runs were presented by the ePrime 3 software (Psychology Software Tools, Sharpsburg, PA) in the following order: 0-back, 1-back, 0-back, 2-back. During the n-back task, participants had to continuously remember a series of rapidly flashing letters (60 letter / session, interleaved by a pause of 250 ms) that were presented by a custom script in a randomized sequence and for a random time interval (1950±100 ms). Participants were requested to provide response only to target stimuli. We denote a letter as a stimulus, which is the same as the one shown previously as a target. Definition of target or non-target stimuli depended on the actual session; in general, the n-back task requires participants to react when a stimulus is the same as the n-th letter before the current stimulus letter. In case of 0-back, the target stimulus was ‘W’, while non-target stimulus was any other letter.

We assessed the neurovascular coupling responses and functional connectivity in the frontal brain cortex, similarly to our previous studies as described in the following papers.

1.	Csipo T, Lipecz A, Owens C, Mukli P, Perry JW, Tarantini S, Balasubramanian P, Nyul-Toth A, Yabluchanska V, Sorond FA, Kellawan JM, Purebl G, Sonntag WE, Csiszar A, Ungvari Z, Yabluchanskiy A. Sleep deprivation impairs cognitive performance, alters task-associated cerebral blood flow and decreases cortical neurovascular coupling-related hemodynamic responses. Sci Rep. 2021;11(1):20994. Epub 20211025. doi: 10.1038/s41598-021-00188-8. PubMed PMID: 34697326; PMCID: PMC8546061.

2.	Mukli P, Csipo T, Lipecz A, Stylianou O, Racz FS, Owens CD, Perry JW, Tarantini S, Sorond FA, Kellawan JM, Purebl G, Yang Y, Sonntag WE, Csiszar A, Ungvari ZI, Yabluchanskiy A. Sleep deprivation alters task-related changes in functional connectivity of the frontal cortex: A near-infrared spectroscopy study. Brain Behav. 2021;11(8):e02135. Epub 20210622. doi: 10.1002/brb3.2135. PubMed PMID: 34156165; PMCID: PMC8413792.

For functional NIRS data analysis, we use a pipeline based on General Linear Model (GLM) approach created using the Brain AnalyzIR toolbox (commit 46c645d). Briefly, measured optical densities are converted to change of hemoglobin concentration using the Beer-Lambert law, pre-whitening of data is performed with an autoregressive model based algorithm, and a discrete cosine transform based high-pass filter (0.08 Hz) is used to remove slow drift. The design matrix included four boxcar regressors (one for each n-back session), which are convolved with a canonical hemodynamic response function to predict brain activation. Parameter estimates (beta-weights), scaling the predictors, are then used for group level statistics. Group level statistics are performed using a mixed effects model. The model is defined in a Wilkinson-Rogers formula of ‘beta ~ -1 + Age-group:n-back condition + (1|Subject)’. Output of the mixed effects model statistics were used for a t-test, where t-contrasts of hemodynamic responses were calculated: [2-back - 0-back] within both the young and elderly group, and [elderly:[2-back-0-back]-young:[2-back-0-back]]. Increased activation was considered significant where q < 0.05 (after false discovery rate correction).

Prior to functional connectivity analysis, artifacts are removed (a) by thresholding after discrete wavelet transformation and then in the frequency domain using a 5th order Butterworth filter with a low‐pass frequency of 0.4 Hz. Preprocessed optical densities are converted into concentration changes of oxy‐ (HbO) and deoxyhemoglobin (HbR) according to the modified Beer–Lambert law, and total hemoglobin (HbT) is calculated as the sum of HbO and HbR.  In the next step, correlation‐based signal improvement (CBSI) is applied to enhance the representation of signal components associated with task‐related brain activity. 

Functional brain networks can be determined using preprocessed HbT time series. We use the entire signal from each n-back session. The strength of functional connection between brain regions (channels) is characterized in terms of temporal correlation of the measured signal pairs. For each subject and each state, we obtain the connection matrices of Pearson's correlation coefficients representing all pairwise combinations of channels. In order to eliminate spurious connections, we apply surrogate thresholding which eliminates all negative and positive Pearson-coefficients with p>0.05.

In the next step, we define the graph theoretical parameters both for binary and weighted undirected networks. During binarization, all nonzero elements of the thresholded connection matrix is replaced by 1, while the latter representation considered these coefficient values as edge weights in the graph. In this study, functional brain networks are characterized using density (D). It can be calculated for each node (local, denoted by a subscript loc) and for the whole network as an average across the nodes. For definitions of these paerameters please check out the following paper:

Rubinov M, Sporns O. Complex network measures of brain connectivity: uses and interpretations. Neuroimage. 2010 Sep;52(3):1059-69. doi: 10.1016/j.neuroimage.2009.10.003. Epub 2009 Oct 9. PMID: 19819337.
