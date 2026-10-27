# Q4 occurrence response metrics

Computed from the existing baseline-corrected AQuA2 dF/F curves at lags -8..+8; no model was refit.

- Response magnitude: maximum baseline-corrected dF/F in the 17-frame window.
- Response latency: lag of that maximum, in frames, defined as Ca frame minus motion peak; negative means Ca peak precedes motion. Earliest lag wins ties.
- Response consistency: median and IQR of Pearson correlations for every pair of occurrence curves from the same pair, calculated over the 17 aligned lag samples. This measures shape similarity while being insensitive to amplitude scaling.
- Pair-level magnitude and latency are reported as median and 25th/75th percentiles across occurrences. A non-positive peak is retained as a signed maximum and flagged; its timing should be interpreted cautiously.
