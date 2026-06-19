import re

import json

import argparse

from collections import OrderedDict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

# Read matched JSON data and build YAML configs for each group of measurements


def loadJSON(path: Path) -> Dict[str, Dict[str, Dict[str, Any]]]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def formatFilterToken(value: Optional[float]) -> str:
    if value is None:
        return "NONE"
    text = f"{value}"
    return text.replace(".", "_")


def parseLayer(name: str) -> str:
    match = re.search(r"layer(\d+)", name, re.IGNORECASE)
    if match:
        # Swap layer 1 and layer 2 in the name (connector issue) 
        layer_num = int(match.group(1))
        if layer_num == 1:
            return "layer2"
        elif layer_num == 2:
            return "layer1"
        else:
            return f"layer{layer_num}"
    return "layer_unknown"


def makeFilename(
    measurement_type: Optional[str],
    source: Optional[str],
    filter_value: Optional[float],
    lv_setting: Optional[int],
    mixture: Optional[str],
    used: Dict[str, int],
) -> str:

    type_token = (measurement_type or "unknown").upper().replace(" ", "_")
    source_token = "SOURCE" + (source or "UNKNOWN")
    filter_token = "" if source == "OFF" else formatFilterToken(filter_value)
    lv_token = f"LVSET{lv_setting}" if lv_setting is not None else "LVSETNONE"

    if mixture is None:
        mix_suffix = "MIXNONE"
    else:
        mix_token = "STD" if str(mixture).lower() == "std" else str(mixture)
        mix_suffix = f"MIX{mix_token}"

    name = f"{type_token}_{source_token}_{filter_token}_{lv_token}_{mix_suffix}"

    used[name] = used.get(name, 0) + 1
    return f"{name}.yaml"


def groupEntries(data: Dict[str, Dict[str, Dict[str, Any]]]) -> Dict[Tuple, List[Dict[str, Any]]]:
    groups: Dict[Tuple, List[Dict[str, Any]]] = {}
    for _, measurements in data.items():
        for measurement_name, info in measurements.items():
            measurement_type = info.get("type")
            mixture = info.get("mixture")
            lv_setting = info.get("lv_setting")
            source = info.get("source_status")
            filter_value = info.get("filter") if source == "ON" else None
            group_key = (measurement_type, mixture, filter_value, lv_setting, source)

            entry = OrderedDict()
            entry["name"] = measurement_name
            if measurement_type != "noise scan" and measurement_type != "source scan":
                entry["layer"] = parseLayer(measurement_name)
            entry["scanned_hv"] = info.get("scanned_hv")
            if "other_hv" in info:
                entry["other_hv"] = info.get("other_hv")
            entry["path_plot"] = info.get("path_plot")
            entry["path_raw"] = info.get("path_raw")

            groups.setdefault(group_key, []).append(entry)
    return groups


def yamlScalar(value: Any) -> str:
    if value is None:
        return "null"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)):
        return str(value)
    text = str(value)
    if text == "" or re.search(r"[:#\n]", text) or text.strip() != text:
        escaped = text.replace("\\", "\\\\").replace('"', '\\"')
        return f"\"{escaped}\""
    return text


def yamlDump(obj: Any, indent: int = 0) -> List[str]:
    lines: List[str] = []
    prefix = " " * indent
    if isinstance(obj, dict):
        for key, value in obj.items():
            if indent == 0 and lines and key in {"mixture", "data", "measurements"}:
                lines.append("")
            if isinstance(value, (dict, list)):
                lines.append(f"{prefix}{key}:")
                lines.extend(yamlDump(value, indent + 2))
            else:
                lines.append(f"{prefix}{key}: {yamlScalar(value)}")
        return lines
    if isinstance(obj, list):
        for item in obj:
            if isinstance(item, dict) and item:
                first_key = next(iter(item))
                first_value = item[first_key]
                lines.append(f"{prefix}- {first_key}: {yamlScalar(first_value)}")
                remaining = OrderedDict((k, v) for k, v in item.items() if k != first_key)
                if remaining:
                    lines.extend(yamlDump(remaining, indent + 2))
            elif isinstance(item, (dict, list)):
                lines.append(f"{prefix}-")
                lines.extend(yamlDump(item, indent + 2))
            else:
                lines.append(f"{prefix}- {yamlScalar(item)}")
        return lines
    lines.append(f"{prefix}{yamlScalar(obj)}")
    return lines


def writeYAML(path: Path, data: OrderedDict) -> None:
    lines = yamlDump(data)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def buildConfigs(data: Dict[str, Dict[str, Dict[str, Any]]], output_dir: Path) -> None:
    groups = groupEntries(data)
    used_names: Dict[str, int] = {}
    output_dir.mkdir(parents=True, exist_ok=True)

    for (measurement_type, mixture, filter_value, lv_setting, source), entries in groups.items():
        config = OrderedDict()
        config["measurement_type"] = measurement_type
        config["mixture"] = mixture
        config["source"] = source
        config["filter"] = "OFF" if source == "OFF" else filter_value
        config["lv_setting"] = lv_setting
        if measurement_type == "efficiency scan":
            data_map: Dict[str, List[OrderedDict]] = {}
            for entry in entries:
                layer = entry.get("layer") or "layer_unknown"
                item = OrderedDict()
                item["hv"] = entry.get("scanned_hv")
                item["name"] = entry.get("name")
                item["raw files"] = entry.get("path_raw")
                if entry.get("path_plot"):
                    item["plot"] = entry.get("path_plot")
                if "other_hv" in entry:
                    item["other_hv"] = entry.get("other_hv")
                data_map.setdefault(layer, []).append(item)

            ordered_data = OrderedDict()
            for layer in sorted(data_map.keys(), key=lambda k: ("layer" not in k, k)):
                ordered_data[layer] = sorted(
                    data_map[layer],
                    key=lambda item: (item.get("hv") is None, item.get("hv")),
                )
            config["data"] = ordered_data
        else:
            if measurement_type in {"noise scan", "source scan"}:
                config["data"] = sorted(
                    entries,
                    key=lambda item: (item.get("scanned_hv") is None, item.get("scanned_hv")),
                )
            else:
                config["data"] = entries

        filename = makeFilename(
            measurement_type,
            source,
            filter_value,
            lv_setting,
            mixture,
            used_names,
        )
        writeYAML(output_dir / filename, config)


# Parser for building YAML configs from matched JSON data
def parseArgs() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build YAML configs from matched JSON data."
    )
    parser.add_argument(
        "json_path",
        help="Path to input JSON file."
    )
    parser.add_argument(
        "output_dir",
        help="Directory to write YAML configs."
    )
    return parser.parse_args()


def main() -> int:
    args = parseArgs()
    json_path = Path(args.json_path).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()

    data = loadJSON(json_path)
    buildConfigs(data, output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
