{ pkgs ? import <nixpkgs> {} }:

pkgs.mkShell {
  packages = with pkgs; [
    clang-tools
    gmp
  ] ++ (if system == "aarch64-darwin" then [ ] else [ gdb ]);
}
