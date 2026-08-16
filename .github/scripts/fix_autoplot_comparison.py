from pathlib import Path

p = Path("vignettes/an_introduction_to_markovchain_package.Rmd")
s = p.read_text(encoding="utf-8")
old = '''```{r mcPlotComparison, fig.cap="The same Markov chain represented with the traditional `plot` method (left) and the optional `ggplot2`-based `autoplot` method (right).", fig.width=12, fig.height=5}\nif (requireNamespace("ggplot2", quietly = TRUE)) {\n  old_par <- par(mfrow = c(1, 2), mar = c(1, 1, 3, 1))\n  plot(mcWeather, main = "Traditional plot()")\n  print(ggplot2::autoplot(mcWeather) +\n          ggplot2::labs(title = "ggplot2::autoplot()"))\n  par(old_par)\n}\n```'''
new = '''```{r mcPlotComparison, fig.cap="The same Markov chain represented with the traditional `plot` method (left) and the optional `ggplot2`-based `autoplot` method (right).", fig.width=6, fig.height=5, fig.show="hold", out.width="49%"}\nplot(mcWeather, main = "Traditional plot()")\nif (requireNamespace("ggplot2", quietly = TRUE)) {\n  print(ggplot2::autoplot(mcWeather) +\n          ggplot2::labs(title = "ggplot2::autoplot()"))\n}\n```'''
if old not in s:
    raise SystemExit("comparison chunk not found")
p.write_text(s.replace(old, new, 1), encoding="utf-8")
Path(".github/workflows/fix-autoplot-comparison.yml").unlink()
Path(".github/scripts/fix_autoplot_comparison.py").unlink()
