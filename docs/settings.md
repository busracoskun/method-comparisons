# Settings

`config/default.json` contains the public defaults used by `run_benchmark.py`

## Current defaults

| Setting | Current value | Change when |
|---|---:|---|
| `adni_input_dir` | `data/adni` | ADNI inputs are stored elsewhere |
| `tcga_input_dir` | `data/tcga` | TCGA inputs are stored elsewhere |
| `cache_input_dir` | `data/cache` | cached private inputs are stored elsewhere |
| `output_dir` | `private_outputs` | outputs should use another project-local folder |
| `rscript` | `Rscript` | a specific R executable is required |
| `python` | `python` | a specific Python executable is required |
| `bootstrap_count` | `10` | a different number of observed-data bootstraps is needed |
| `bootstrap_start` | `1` | resuming at another bootstrap index |
| `cores` | `1` | more local compute is available |
| `save_models` | `false` | fitted model objects should be retained |
| `repetitions` | `100` | fixed manuscript design — do not change for reproduction |
| `r_home` | `null` | R cannot find its installation |
| `r_library` | `null` | packages use a private R library |
| `library_paths` | `[]` | additional shared libraries are required |

`observed_methods` and `semisynthetic_methods` list the methods run by default

Method order is preserved because it can affect the original random-number streams

## What users can change

For one run, use command-line options

```bash
python run_benchmark.py adni-observed --methods mofa,cmf --count 5 --start 1 --output private_outputs/test
```

Available run options

| Option | Purpose |
|---|---|
| `pipeline` | choose ADNI or TCGA and observed or semi-synthetic analysis |
| `--config` | load a JSON file that overrides public defaults |
| `--methods` | run a comma-separated method subset |
| `--engine` | run all methods, R methods, or Python methods |
| `--mode` | select a valid stage for the chosen pipeline |
| `--count` | override the observed-data bootstrap count |
| `--start` | override the first bootstrap index |
| `--task` | select semi-synthetic task `1..600` |
| `--output` | select a project-local output folder |
| `--check` | validate settings and inputs without fitting |

Use `config/host.local.json` for machine-specific paths and executables

The file is Git-ignored and overrides `config/default.json` automatically

```json
{
  "adni_input_dir": "/authorized/path/to/adni",
  "tcga_input_dir": "/authorized/path/to/tcga",
  "cores": 4,
  "python": "python3"
}
```

Do not edit constants inside the analysis scripts for an ordinary run

Those values preserve the manuscript analysis unless a nearby comment explicitly identifies them as user settings
