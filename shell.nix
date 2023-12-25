{ pkgs ? import <nixpkgs> {} }:
let
  my-python-packages = ps: with ps; [
    pandas
    rdkit
    numpy
    scikit-learn
    seaborn
    matplotlib
    (
    buildPythonPackage rec {
      pname = "padelpy";
      version = "0.1.14";
      src = fetchPypi {
        inherit pname version;
	sha256 = "8ba7f2478f6e7732b96ee385b98d00d2af3d06f69225a1b6b39a8ad9dfb40a3d";
	};
    })
  ];
  my-python = pkgs.python311.withPackages my-python-packages;
in
pkgs.mkShell {
  buildInputs = [
    my-python
    pkgs.jre8
    pkgs.openbabel2
    pkgs.unrar
    ];
}
