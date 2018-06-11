with import <nixpkgs> {};

stdenv.mkDerivation rec {
  name = "boltzbub-env";
  env = buildEnv { name = name; paths = buildInputs; };

  buildInputs = [
    git
    gcc
    cmake
  ];

  shellHook = ''
    export NIX_SHELL_NAME="${name}"
  '';
}
