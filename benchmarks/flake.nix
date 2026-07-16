{
  description = "Reproducible Molfig benchmark environment";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/35d3407a3816f3b341d8cf1d60abaf2b7b8166ac";
  inputs.nixpkgs-typst-0-14-0.url = "github:NixOS/nixpkgs/1d0bb7b61b251a261b0963aacf4b141e770a4f1d";
  inputs.nixpkgs-typst-0-14-1.url = "github:NixOS/nixpkgs/533f8cf6e2396e335b5fa9041dc095658566d0a9";
  inputs.nixpkgs-typst-0-14-2.url = "github:NixOS/nixpkgs/b6e4a72e837aee46a59a76e1a3edcf10fa6eebb1";

  outputs = inputs@{ nixpkgs, ... }:
    let
      systems = [
        "aarch64-darwin"
        "aarch64-linux"
        "x86_64-linux"
      ];
      forAllSystems = function:
        nixpkgs.lib.genAttrs systems (system: function {
          inherit system;
          pkgs = nixpkgs.legacyPackages.${system};
        });
      mkBenchmarkShell = {
        pkgs,
        typst,
        typstVersion,
        nixpkgsRevision,
      }:
        pkgs.mkShellNoCC {
          packages = [
            typst
            pkgs.hyperfine
          ];

          MOLFIG_BENCH_ENV = "1";
          MOLFIG_BENCH_TYPST_VERSION = typstVersion;
          MOLFIG_BENCH_NIXPKGS_REV = nixpkgsRevision;
        };
    in
    {
      devShells = forAllSystems ({ system, pkgs }: {
        default = mkBenchmarkShell {
          inherit pkgs;
          typst = pkgs.typst;
          typstVersion = "0.15.0";
          nixpkgsRevision = "35d3407a3816f3b341d8cf1d60abaf2b7b8166ac";
        };
        typst-0_14_0 = mkBenchmarkShell {
          inherit pkgs;
          typst = inputs.nixpkgs-typst-0-14-0.legacyPackages.${system}.typst;
          typstVersion = "0.14.0";
          nixpkgsRevision = "1d0bb7b61b251a261b0963aacf4b141e770a4f1d";
        };
        typst-0_14_1 = mkBenchmarkShell {
          inherit pkgs;
          typst = inputs.nixpkgs-typst-0-14-1.legacyPackages.${system}.typst;
          typstVersion = "0.14.1";
          nixpkgsRevision = "533f8cf6e2396e335b5fa9041dc095658566d0a9";
        };
        typst-0_14_2 = mkBenchmarkShell {
          inherit pkgs;
          typst = inputs.nixpkgs-typst-0-14-2.legacyPackages.${system}.typst;
          typstVersion = "0.14.2";
          nixpkgsRevision = "b6e4a72e837aee46a59a76e1a3edcf10fa6eebb1";
        };
      });
    };
}
