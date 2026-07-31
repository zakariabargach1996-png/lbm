# LBM Jupyter notebooks

Two complementary notebooks are available:

- `lbm_results_report.ipynb`: results-first executable report, with validation and
  scaling evidence before a compact methods section;
- `lbm_results_report.pdf`: print-ready version of the executed results report;
- `lbm_code_presentation.ipynb`: theory and implementation presentation with
  slideshow metadata;
- `lbm_code_presentation.slides.html`: exported Reveal.js browser presentation.

Launch the notebook from the repository root:

```bash
conda run -n meep jupyter notebook presentation/lbm_code_presentation.ipynb
```

In Jupyter, use the notebook normally as an executable report. To present it,
open the exported HTML file in a browser and navigate with Space or the arrow
keys.

To rebuild and execute the results report:

```bash
conda run -n meep python presentation/build_lbm_report.py
conda run -n meep jupyter nbconvert \
  --to notebook --execute --inplace \
  presentation/lbm_results_report.ipynb
```

To rebuild the notebook after editing `build_lbm_presentation.py`:

```bash
conda run -n meep python presentation/build_lbm_presentation.py
conda run -n meep jupyter nbconvert \
  --to notebook --execute --inplace \
  presentation/lbm_code_presentation.ipynb
conda run -n meep jupyter nbconvert \
  --to slides --no-input \
  presentation/lbm_code_presentation.ipynb
```

The compact density-wave asset was reduced from a fresh 128 x 128 solver run by
`prepare_density_profiles.py`. Other figures and tables are read directly from
the archived validation and performance results in the repository.
