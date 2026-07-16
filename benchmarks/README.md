# Molfig benchmarks

This suite measures complete Typst compiler processes against the molecular structure data already published in `package/examples/data`. It covers PDB, text mmCIF, and BinaryCIF inputs as well as cartoon, spacefill, and molecular surface mesh generation.

The default `export` mode measures parsing, Model/Structure/Unit construction, mesh generation, and OBJ serialization without maquette rendering. The optional `render` mode additionally measures OBJ material handling, maquette rasterization, and PDF assembly.

Run every case with Typst 0.15.0 and the Molfig version declared in `package/typst.toml`:

```sh
nix develop ./benchmarks --command benchmarks/run.sh
```

Molfig's supported compiler versions are available as locked development shells. The default shell is Typst 0.15.0:

| Typst | Nix development shell |
| --- | --- |
| 0.14.0 | `./benchmarks#typst-0_14_0` |
| 0.14.1 | `./benchmarks#typst-0_14_1` |
| 0.14.2 | `./benchmarks#typst-0_14_2` |
| 0.15.0 | `./benchmarks` |

For example, run the same case with Typst 0.14.0:

```sh
nix develop './benchmarks#typst-0_14_0' --command \
  benchmarks/run.sh 1crn-bcif-spacefill
```

Select a released version and one or more cases:

```sh
nix develop ./benchmarks --command \
  benchmarks/run.sh --version "$VERSION" \
  1crn-bcif-spacefill 9r1o-pdb-cartoon
```

Measure the end-to-end render path:

```sh
nix develop ./benchmarks --command \
  benchmarks/run.sh --mode render --runs 10 9r1o-pdb-cartoon
```

For a locally installed package, use `--namespace local`. The package version can also be supplied through `MOLFIG_VERSION`; command-line options take precedence.

The runner requires Nix with flakes enabled and must be invoked through `nix develop` as shown above. `flake.lock` fixes the Nixpkgs sources and content hashes for all four Typst versions and hyperfine 1.20.0. The runner verifies the selected compiler version instead of falling back to the host `PATH`, and prints the complete environment with every result. The first invocation may download the pinned Nix store paths; downloads happen before hyperfine starts measuring.

Benchmark inputs are PDB archive data from RCSB PDB / wwPDB and are available under CC0 1.0. Per-entry provenance and attribution are recorded in `package/examples/data/README.md` and `package/examples/data/ATTRIBUTION.tsv`.
