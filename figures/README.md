# Scripted figures

Regenerate with `python -m pipeline.cli figures` after held-out evaluation. PNG files support README viewing; PDF files support export.

## Figure 1. Humanization and interface audit

Substantial conformation mismatch and residual repaired clashes limit interpretation.

![Humanization and interface audit](figure_1.png)

## Figure 2. Discovery E2 interface conservation

Some interface positions vary in the retrieved discovery sequences; frequencies are sampling-dependent.

![Discovery E2 interface conservation](figure_2.png)

## Figure 3. Experimental mutation score benchmark

The interface scoring gate failed; candidate model outputs cannot support improved-binding claims.

![Experimental mutation score benchmark](figure_3.png)

## Figure 4. Single-state and multi-state discovery tradeoff

Single-state score and worst paired delta measure different model properties.

![Single-state and multi-state discovery tradeoff](figure_4.png)

## Figure 5. Frozen candidates across discovery E2 states

WT is the zero baseline on every state; positive deltas favor WT in the model.

![Frozen candidates across discovery E2 states](figure_5.png)

## Figure 6. Frozen held-out evaluation

The headline effect is a model-score comparison under a failed experimental validation gate.

![Frozen held-out evaluation](figure_6.png)

## Figure 7. Sequence and exposure tradeoff

Soluble-model differences are geometric proxies, not expression or solubility measurements.

![Sequence and exposure tradeoff](figure_7.png)

## Figure 8. Measured resource profile

Per-stage maxima across successful executions; cache skips do not replace full computation time. Original bootstrap telemetry remains unavailable.

![Measured resource profile](figure_8.png)
