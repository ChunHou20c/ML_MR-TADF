{
  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixos-unstable";
  };
  
  outputs = { self, nixpkgs }:
    let
      pkgs = import nixpkgs { system = "x86_64-linux"; };
      # A list of shell names and their Python versions
      pythonVersions = {
        python311 = pkgs.python311;
        default = pkgs.python311;
      };
      my-python-packages = ps: with ps; [
        pandas
        rdkit
        numpy
        scikit-learn
        seaborn
        matplotlib
        python-dotenv
        (
          buildPythonPackage rec {
            pname = "padelpy";
            version = "0.1.14";
            src = fetchPypi {
              inherit pname version;
              sha256 = "8ba7f2478f6e7732b96ee385b98d00d2af3d06f69225a1b6b39a8ad9dfb40a3d";
            };
          })

        (
          buildPythonPackage rec {
            pname = "Boruta";
            version = "0.4.3";
            src = fetchPypi {
              inherit pname version;
              sha256 = "abd67119bb88c7d1675b0dfc78072e50b5ba1d32b7cc7664258a640ad3a495a4";
            };
          })
        ];
        my-python = pkgs.python311.withPackages my-python-packages;
      # A function to make a shell with a python version
      makePythonShell = shellName: pythonPackage: pkgs.mkShell {
        # You could add extra packages you need here too

        packages = [
          pkgs.jre8
          my-python
          pkgs.openbabel2
        ]; 
        # You can also add commands that run on shell startup with shellHook
        shellHook = ''
          echo "Now entering ${shellName} environment."
        '';
      };
    in
    {
      # mapAttrs runs the given function (makePythonShell) against every value
      # in the attribute set (pythonVersions) and returns a new set
      devShells.x86_64-linux = builtins.mapAttrs makePythonShell pythonVersions;
    };
}

