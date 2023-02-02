# PWN: Prioritization with a Warped Network

![Graphical overview](preview.png)

[![Paper on BioMed Central](https://img.shields.io/badge/BMC%20Bioinformatics-10.1186%2Fs12859--023--05227--x-04caa8?labelColor=1b3051&logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQBAMAAADt3eJSAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAAD1BMVEUPbY8NZYkn3+kAKFr////4rfcQAAAAAnRSTlPAwCYN2MYAAAABYktHRASPaNlRAAAACXBIWXMAAABlAAAAZQGq1hevAAAAB3RJTUUH5AEODjgScgZNDgAAADRJREFUCNc1yEkNACAQwECCBCTgoFn/3tiD9jXp2ndaiDgC8VcCMauA6NVA1Bogcn0g4ggeZjsUofG3VpgAAAAldEVYdGRhdGU6Y3JlYXRlADIwMjAtMDEtMTRUMTQ6NTY6MTgrMDE6MDDu7x8XAAAAJXRFWHRkYXRlOm1vZGlmeQAyMDIwLTAxLTE0VDE0OjU2OjE4KzAxOjAwn7KnqwAAABl0RVh0U29mdHdhcmUAd3d3Lmlua3NjYXBlLm9yZ5vuPBoAAABXelRYdFJhdyBwcm9maWxlIHR5cGUgaXB0YwAAeJzj8gwIcVYoKMpPy8xJ5VIAAyMLLmMLEyMTS5MUAxMgRIA0w2QDI7NUIMvY1MjEzMQcxAfLgEigSi4A6hcRdPJCNZUAAAAASUVORK5CYII=)](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-023-05227-x)
[![Licensed under BSD-3-Clause-Clear](https://img.shields.io/github/license/Standigm/PWN)](LICENSE)

## Installation

```bash
> pip install git+https://github.com/Standigm/PWN
```

## Experiments

Notice that all experiment scripts and results are located under [`./experiments`](experiments).
Following section explains how to reproduce the results.

### Collect Experiments Results

```bash
> pip install "pwn[experiment] @ git+https://github.com/Standigm/PWN"  # additional requirements
> python run.py
```

### Visualize

To generate figures, install [R](https://www.r-project.org) and run the following commands:

```bash
> Rscript -e 'chooseCRANmirror(ind=1); install.packages("tidyverse")'  # requirements
> Rscript figure.R
```

## License

The source code is distributed under [BSD 3-Clause Clear License](LICENSE), while each data used in experiments is
distributed under different licenses.

Please use the following command for license information of other software used in this project:

```bash
> poetry show --only main | grep -E -o '^[^ ]+' | xargs pip-licenses -o license -p
```
