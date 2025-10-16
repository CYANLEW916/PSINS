# Repository Agent Instructions

## Scope

These instructions apply to the entire repository unless superseded by a more specific `AGENTS.md` placed in a subdirectory.

## Coding guidelines

- This codebase primarily uses MATLAB scripts. Follow existing naming conventions and prefer functions over scripts when practical.
- Keep line endings Unix-style (LF) and limit lines to 100 characters when feasible.
- Document any new functions with a brief comment header describing inputs and outputs.
- Maintain existing indentation (typically 4 spaces) and avoid tabs.

## Testing and validation

- When modifying algorithms, add or update example usage scripts under `demos/` if applicable.

## Documentation
- Update `README.md` or the appropriate documentation under `doc/` for user-facing changes.
- Include clear inline comments for complex mathematical steps or transformations.

## Pull requests
- Provide a concise summary of changes and list any tests executed.
- Reference affected modules or demo scripts so reviewers can reproduce results easily.
