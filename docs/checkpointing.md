# Checkpointing in GWAT's MCMC

GWAT's MCMC runs on the BayesShip sampler, and checkpointing is done entirely by BayesShip. BayesShip's
`docs/checkpointing.md` describes the files it writes, when it writes them, how a resume works, and what each
option means. This page only covers how those options are set from GWAT.

## Where each option is set

GWAT's sampler setup functions in `src/mcmc_gw.cpp` set the options on the sampler before calling `sample()`:
- `PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW` (including the `GWModel` overload in `include/gwat/mcmc_gw.h`);
- `RJPTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW`.

**Function arguments:**

| Argument | Sampler member |
|---|---|
| `outputDir` | `outputDir` (needs a trailing `/`) |
| `outputFileMoniker` | `outputFileMoniker` |
| `batchSize` | `batchSize` |
| `ignoreExistingCheckpoint` | `ignoreExistingCheckpoint` |

**Fields of `MCMC_modification_struct`** (`include/gwat/mcmc_gw.h`). They have the same names and defaults as
the `bayesshipSampler` members, and `set_checkpoint_options` copies them onto the sampler:

| Field | Default |
|---|---|
| `hotChainHistoryLength` | 2000 |
| `checkpointGenerationsToKeep` | 2 |
| `resumeValueRelTol` | 1e-10 |
| `wallTimeLimitSeconds` | 0 (no limit of its own) |
| `wallTimeMarginSeconds` | 120 |
| `handleStopSignals` | true |

A struct that is left untouched gives BayesShip's own defaults. To change them:

```cpp
MCMC_modification_struct mod_struct;
mod_struct.hotChainHistoryLength = 500;
mod_struct.wallTimeLimitSeconds = 3600;
auto* sampler = PTMCMC_MH_dynamic_PT_alloc_uncorrelated_GW(/* ... */, &mod_struct, /* ... */);
if (sampler->stoppedEarly) {
  // A checkpoint is committed; call again with the same outputDir to resume.
}
```

The environment variables `BAYESSHIP_WALL_TIME_LIMIT`, `BAYESSHIP_WALL_TIME_MARGIN` and `SLURM_JOB_END_TIME`
work as described in BayesShip's notes.

## Notes

- **Defaults:** the defaults are copied from BayesShip by hand. If BayesShip changes one of its defaults,
  update the matching field here.
- **Resuming:** these options are read again on every start; none of them is restored from the checkpoint.
- **`ignoreExistingCheckpoint`:** prefer an empty `outputDir` over `true`. A fresh run in a directory that
  already holds a checkpoint writes over the old run's batch files.
