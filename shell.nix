with import <nixpkgs> {};

stdenvNoCC.mkDerivation rec {
  name = "boltzbub-env";
  env = buildEnv { name = name; paths = buildInputs; };

  buildInputs = [
    git
    gcc8
    cmake
  ];

  shellHook = ''
    export NIX_SHELL_NAME="${name}"
  '';
}
