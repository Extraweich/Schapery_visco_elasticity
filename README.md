# Schapery_visco_elasticity
UMAT and VUMAT of a non-linear visco-elastic model after Schapery of 3 Kelvin-Voigt elements:

<br/><br/>

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://github.com/Extraweich/Schapery_visco_elasticity/blob/main/Figures/Kelvin-Voigt3_light.svg?raw=true">
    <img alt="Schematic of implemented homogenization procedures" src="https://github.com/Extraweich/Schapery_visco_elasticity/blob/main/Figures/Kelvin-Voigt3.svg?raw=true", width="50%">
  </picture>
</p>

<br/><br/>

The model is part of my PhD thesis, where you can also find a full derivation of the 3D formulation. In case you would like to use the model within your own research, please cite the resource correctly. Here is a [link](https://publikationen.bibliothek.kit.edu/1000194396) to my thesis. Alternatively, you can use the following BibTex:

```bibtex
@phdthesis{Christ.2026,
    author       = {Christ, Nicolas},
    year         = {2026},
    title        = {On the influence of temperature and humidity on interfaces in carbon fiber reinforced polyamide 6},
    doi          = {10.5445/IR/1000194396},
    publisher    = {{Karlsruher Institut für Technologie (KIT)}},
    keywords     = {nonlinear viscoelasticity, Schapery, single fiber pull-out test, climbing drum peel test},
    pagetotal    = {261},
    school       = {Karlsruher Institut für Technologie (KIT)},
    language     = {english}
}
```


Remark: The current version has some errors, which affect suddenly applied load cases (e.g. creep experiments). These errors have been resolved and will be uploaded within the next weeks.
