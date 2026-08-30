# Stage 9 batch 2: sample-size sensitivity (freeze §1.55)

- n in {500, 750, 1250}, R=200/cell, paired, 600/600 ok

| n | locus det | secondary power | overall recovery | cand exact | direction | count_exact | extra rate | mean rep runtime (s) |
|---|---|---|---|---|---|---|---|---|
| 500 | 0.285 | 0.200 | 0.000 | 0.439 | 0.439 | 0.351 | 0.068 | 16.1 |
| 750 | 0.605 | 0.535 | 0.000 | 0.607 | 0.620 | 0.521 | 0.132 | 25.5 |
| 1250 | 0.980 | 0.950 | 0.000 | 0.901 | 0.903 | 0.781 | 0.281 | 46.2 |

- reference: main n=1000 at R=500 (Stage 8.2): full/secondary 0.778
- 核对点：recovery 随 n 单调不降；runtime 随 n 增长可接受

- 严格匹配口径（Stage 8.2 约定）：full recovery 0.530, secondary 0.530

### 严格匹配口径按 n 分组

| n | full recovery | secondary recovery |
|---|---|---|---|
| 500 | 0.155 | 0.155 |
| 750 | 0.490 | 0.490 |
| 1250 | 0.945 | 0.945 |

### 严格匹配口径按 n 分组

| n | full recovery | secondary recovery |
|---|---|---|---|
| 500 | 0.155 | 0.155 |
| 750 | 0.490 | 0.490 |
| 1250 | 0.945 | 0.945 |
