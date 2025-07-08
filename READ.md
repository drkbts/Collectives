Search for "C/C++: Edit configurations (UI)".

In the configuration settings, add the generated compile_commands.json file to the compileCommands field: (without quotes)

"compileCommands": "${workspaceFolder}/compile_commands.json"

compile_commands generate by bcc (compiled from source; copied from bazel-bin/bcc to /usr/local/bin/bcc)

https://github.com/kiron1/bazel-compile-commands?tab=readme-ov-file

# comment