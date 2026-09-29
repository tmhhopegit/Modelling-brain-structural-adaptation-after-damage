# Lesion growth: predicting Chilean patients from acute scans

This repository is built on code used in Loughnan et al., 2019: https://pmc.ncbi.nlm.nih.gov/articles/PMC6831940/

The question: a model trained on PLORAS patients predicts language scores from
demographics and the lesion. Applied to Chilean patients, its predictions are
over-optimistic. The Chilean patients were scanned within two weeks of the
stroke but assessed months later. The hypothesis is that their lesions were
still growing, so the acute scan understates the damage present at
assessment.

The analysis tests this in three steps:

1. **Train on PLORAS and predict Chile.** Train on PLORAS patients' first scans, predict the Chilean patients from their acute scans, and test whether their errors are more optimistic than PLORAS patients' at the same actual score.
2. **Learn how lesions change.** Use PLORAS patients scanned more than once.
3. **Project and predict again.** Project each Chilean lesion forward to the date of the assessment, predict again, and compare again.

## Layout

```
lesion_growth/
  README.md
  config_lesion_growth.m     every setting: file paths, column names, scores, predictors, model
  setup_lesion_growth.m      adds the folders below to the path
  example_lesion_growth.m    the whole analysis, section by section
  data/
    load_ploras.m            PLORAS spreadsheet + lesion images   (PLORAS_Raw's constructor)
    load_chile.m             Chilean lesion folder + demographics  (create_chile_mat.m)
    chile_summary.m          sample description (create_chile_mat's summary_stats)
    code_demographic.m       handedness / sex / language coded the same way for both datasets
    read_lesion_images.m, lesion_volumes.m
  analysis/
    lesion_growth_analysis.m the three steps; returns predictions, comparisons and growth models
    lesion_pca_fit.m / lesion_pca_project.m   lesion components
    predictor_matrix.m, fit_gp.m, cv_predict_gp.m
    compare_errors.m         error ~ actual score + dataset (the original's anovan), plus group summaries
    longitudinal_pairs.m     earlier -> later scan pairs from PLORAS
    learn_lesion_growth.m    one GP per later lesion measure, with cross-validated accuracy
    project_lesions_forward.m
    plot_lesion_growth.m     predicted vs actual: PLORAS, Chile acute, Chile projected
  tests/
    LesionGrowthTest.m       runtests('LesionGrowthTest')
    make_synthetic_lesion_data.m
```

## Setup

You need:

- MATLAB R2019b or later, with the Statistics and Machine Learning Toolbox (`fitrgp`, `anovan`, `cvpartition`);
- SPM12 on the path to read the lesion images.

```matlab
run('...\lesion_growth\setup_lesion_growth.m')
runtests('LesionGrowthTest')     % about a minute; the loader test needs SPM
edit config_lesion_growth        % put in your column names
edit example_lesion_growth
```

**Column names.** The data spreadsheet columns are named in the config. The
original read them by number (age = data column 4, time post stroke = 5, and so
on; the numbers are listed in the config), and I can't tell from the code what
the headers are called. Please fill them in.

## Output

`R = lesion_growth_analysis(P, C, cfg)` gives:

- **`R.summary`**: one row per score and condition (`acute scan` or `projected to assessment`). For PLORAS (A) and Chile (B) it gives N, mean error (predicted − actual; positive means over-optimistic), mean absolute error and the correlation of predicted with actual. It then gives the original's ANOVA of error on actual score plus dataset, with `AdjustedDifference` (Chile's extra error at the same actual score), `F_group` and `p_group`.
- **`R.scores(k)`**: every patient's actual, predicted and (Chile only) projected-then-predicted score.
- **`R.growth.cv`**: how well each later lesion measure is predicted for patients the growth model did not see, next to how well "no change" does. If lesion growth can't be predicted better than "no change", step 3 cannot work. The original never checked this.
- **`R.chile_summary`**: the sample description that `create_chile_mat` saved.

