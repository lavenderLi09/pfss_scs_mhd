# off_2270: anomalous high-speed component diagnosis

- Sign convention: `+r` is outward.
- Component definition: largest connected high-speed component in each frame.
- Focus: plasma beta, Alfvenic speed, momentum force budget, residual, acceleration and direction.

## Global stats

- Frames analyzed: `19`
- Mean `jxb_dom_ratio_p50`: `433.7` (>1 means |(JxB)_r| dominates |-(dp/dr)|+|rho g_r| on median)
- Mean `comp_a_res_r_pos_frac`: `0.2588`
- Mean `comp_divv_p50`: `94.21`

- Peak high-speed (`comp_vabs_p95`) at `off0090`: `12.8322`, beta_p50=`0.00732351`, vA_p50=`21.7943`
- Lowest beta (`comp_beta_p50`) at `off0010`: `0.000313352`, vabs_p95=`8.75001`
- Highest Alfven speed (`comp_va_p50`) at `off0010`: `108.27`, beta_p50=`0.000313352`

## Top 10 high-speed frames

| frame | time | vabs_p95 | beta_p50 | va_p50 | jxb_p50 | f_res_p50 | a_res_p50 | frac(a_res>0) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| off0090 | 9.300000e-01 | 12.8322 | 0.00732351 | 21.7943 | -0.767557 | -0.76782 | -1348.28 | 0.0023 |
| off0080 | 8.300000e-01 | 12.739 | 0.00683053 | 22.6396 | -0.834193 | -0.834025 | -1349.13 | 0.0025 |
| off0070 | 7.300000e-01 | 12.5998 | 0.00634123 | 23.338 | -0.975246 | -0.975573 | -1351.99 | 0.0027 |
| off0160 | 1.630000e+00 | 12.5478 | 0.273146 | 3.695 | -4.28484e-05 | -5.6914e-06 | -0.227553 | 0.4739 |
| off0170 | 1.730000e+00 | 12.54 | 0.286439 | 3.61892 | -5.38711e-05 | -2.42634e-05 | -0.512823 | 0.4539 |
| off0180 | 1.830000e+00 | 12.5227 | 0.241865 | 3.5804 | -4.86217e-05 | -4.77321e-05 | -0.619987 | 0.4316 |
| off0150 | 1.530000e+00 | 12.4924 | 0.261267 | 3.77736 | -3.49254e-05 | -6.43185e-07 | -0.0206959 | 0.4926 |
| off0140 | 1.430000e+00 | 12.4835 | 0.244713 | 3.9072 | -1.65071e-05 | 4.1253e-07 | 0.0342353 | 0.5133 |
| off0060 | 6.300000e-01 | 12.3811 | 0.0059853 | 24.5028 | -1.65336 | -1.65365 | -1766.61 | 0.0031 |
| off0130 | 1.330000e+00 | 12.3786 | 0.232328 | 4.02709 | -6.86328e-06 | 8.50424e-07 | 0.126079 | 0.5341 |
