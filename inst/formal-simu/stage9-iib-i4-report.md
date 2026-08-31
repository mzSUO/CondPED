# Stage 9 batch 12: II-B-I4 strict same-sign comparison (freeze §1.58 ext.)

- Sigma_P = I4; R1/R3 with strict ++++ rejection, R=500 each, 1000/1000 ok

| scenario | n | rep_exact | rep_card | rho MAE | thr acc | Irr | cand exact | direction | eta bias | eta RMSE | mean attempts |
|---|---|---|---|---|---|---|---|---|---|---|---|
| highly_representable | 500 | 0.740 | 0.740 | 0.0222 | 0.958 | 0.740 | 1.000 | 0.588 | 0.0000 | 0.0444 | 8.4 |
| strongly_nonredundant | 500 | 0.790 | 0.790 | 0.0607 | 0.986 | 0.790 | 1.000 | 1.000 | -0.0010 | 0.0457 | 7.9 |

- empirical acceptance rate: highly_representable 0.119; strongly_nonredundant 0.127 (probe expected ~0.13; material deviation must be flagged)
- wall 10078 s
