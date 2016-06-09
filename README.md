# Gravastar
Crappy code to reproduce the graphs for an anisotropic gravastar (https://arxiv.org/abs/0706.1513)
Running "itp-multiple-run.sh" on a cluster with slurm calls the program "anicomp.py" 6000 times for varying masses M and inner radius r1, returning the thickness/M and compactness M/r2 for every configuration, which creates the plot "anisotropic compactness.png" (and a lot of errors)
