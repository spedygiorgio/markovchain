from pathlib import Path
p = Path('vignettes/an_introduction_to_markovchain_package.Rmd')
s = p.read_text(encoding='utf-8')
needle = '\n `plot` method for `markovchain` objects is a wrapper of `plot.igraph`'
if needle not in s:
    raise SystemExit('expected leading space before plot paragraph not found')
s = s.replace(needle, '\n`plot` method for `markovchain` objects is a wrapper of `plot.igraph`', 1)
p.write_text(s, encoding='utf-8')
Path('.github/workflows/fix-autoplot-vignette.yml').unlink()
Path('.github/scripts/fix_vignette.py').unlink()
