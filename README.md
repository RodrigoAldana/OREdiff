# OREdiff
Optimal Robust Exact (ORE) numerical differentiation scripts

OREdiff provides scripts for optimal and robust exact numerical differentiation of noisy signals. The project includes implementations in both MATLAB and Python, with testbenches for each version.

If you use this repository in your work, please cite it as:
```
@ARTICLE{OREdiff2025,
  author={Aldana-López, Rodrigo and Seeber, Richard and Haimovich, Hernan and Gómez-Gutiérrez, David},
  journal={IEEE Transactions on Automatic Control}, 
  title={Optimal robust exact first-order differentiators with Lipschitz-continuous output}, 
  year={2025},
  volume={},
  number={},
  pages={1-8},
  keywords={Noise;Accuracy;Convergence;Robustness;Noise measurement;Tuning;Fault diagnosis;Doppler effect;Vibrations;Upper bound},
  doi={10.1109/TAC.2025.3555481}}

```

## Usage

### MATLAB
To test the MATLAB version, run the following script:
```
matlab matlab/diff_testbench.m
```
### Python
To test the Python version, run:
```
python python/diff.py
```

