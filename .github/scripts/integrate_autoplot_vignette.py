from pathlib import Path

path = Path('vignettes/an_introduction_to_markovchain_package.Rmd')
text = path.read_text(encoding='utf-8')
needle = "The transition matrix of a `markovchain` object can be displayed using `print` or `show` methods (the latter being less verbose). Similarly, the underlying transition probability diagram can be plotted by the use of `plot` method (as shown in Figure \\@ref(fig:mcPlot)) which is based on \\pkg{igraph} package [@pkg:igraph]."
if needle not in text:
    raise SystemExit('Expected plotting paragraph was not found')

addition = r'''

### ggplot2-based plotting

As an alternative to the standard `plot` method, the package also provides an `autoplot` method based on `ggplot2`. The method returns a regular `ggplot` object, so the resulting graph can be further customised using the grammar of graphics. States are represented as nodes, transition probabilities as directed edges, and communicating classes are shown with different node fills. By default, only transitions with positive probability are displayed and their probabilities are labelled on the edges.

The method is available when the optional `ggplot2` package is installed:

```{r mcAutoplot, eval=FALSE}
ggplot2::autoplot(mcWeather)
```

The returned object can be customised like any other `ggplot` object. For example, transition labels can be hidden and very small probabilities can be omitted:

```{r mcAutoplotCustom, eval=FALSE}
ggplot2::autoplot(mcWeather,
                  show_probabilities = FALSE,
                  threshold = 0.2)
```
'''
text = text.replace(needle, needle + addition, 1)
path.write_text(text, encoding='utf-8')

Path('.github/workflows/integrate-autoplot-vignette.yml').unlink()
Path('.github/scripts/integrate_autoplot_vignette.py').unlink()
