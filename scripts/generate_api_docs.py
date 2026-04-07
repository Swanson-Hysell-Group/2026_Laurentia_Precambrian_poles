"""Generate markdown API documentation for pole_tools.py using griffe.

Parses the source file directly (no imports needed) and writes a
MyST-compatible markdown file to resources/pole_tools_api.md.

Usage:
    python scripts/generate_api_docs.py
"""
from pathlib import Path
from griffe import load, parse_google

root = Path(__file__).resolve().parent.parent

# Load the module from source (no import required)
module = load("pole_tools", search_paths=[root / "pole_notebooks"])

lines = [
    "# `pole_tools` API Reference\n",
]

if module.docstring:
    lines.append(f"{module.docstring.value}\n")

for name, func in sorted(module.functions.items()):
    # Build signature
    params = []
    for param in func.parameters:
        if param.name == "self":
            continue
        if param.default is not None:
            params.append(f"{param.name}={param.default}")
        else:
            params.append(param.name)
    sig = ", ".join(params)

    lines.append(f"## `{name}`\n")
    lines.append(f"```python\n{name}({sig})\n```\n")

    if not func.docstring:
        lines.append("---\n")
        continue

    # Parse the Google-style docstring into structured sections
    sections = parse_google(func.docstring)

    for section in sections:
        kind = section.kind.value

        if kind == "text":
            lines.append(f"{section.value}\n")

        elif kind == "parameters":
            lines.append("**Parameters**\n")
            for param in section.value:
                annotation = f" (`{param.annotation}`)" if param.annotation else ""
                desc = param.description.replace("\n", " ")
                lines.append(f"- **{param.name}**{annotation} — {desc}")
            lines.append("")

        elif kind == "returns":
            lines.append("**Returns**\n")
            # Join all return entries into a single description since
            # griffe sometimes splits multi-line returns into separate items
            parts = []
            for ret in section.value:
                annotation = f"`{ret.annotation}` — " if ret.annotation else ""
                desc = ret.description.replace("\n", " ")
                parts.append(f"{annotation}{desc}")
            lines.append(f"- {' '.join(parts)}")
            lines.append("")

    lines.append("---\n")

# Remove trailing separator
if lines and lines[-1] == "---\n":
    lines.pop()

out_path = root / "resources" / "pole_tools_api.md"
out_path.write_text("\n".join(lines))
print(f"API docs written to {out_path}")
