# LBM Jupyter presentation

The presentation is available in two forms:

- `lbm_code_presentation.ipynb`: executable Jupyter notebook with slideshow metadata;
- `lbm_code_presentation.slides.html`: exported Reveal.js browser presentation.

Launch the notebook from the repository root:

```bash
conda run -n meep jupyter notebook presentation/lbm_code_presentation.ipynb
```

In Jupyter, use the notebook normally as an executable report. To present it,
open the exported HTML file in a browser and navigate with Space or the arrow
keys.

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
