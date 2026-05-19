import argparse
import re, os
import sys

import json

from collections import OrderedDict
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# Match timestamps in formats like "YYYY-MM-DD_HH-MM-SS", "YYYY-MM-DD-HH-MM-SS", 
# "YYYY-MM-DD_HH_MM_SS", "YYYY-MM-DD-HH_MM_SS", and with optional seconds
TIMESTAMP_RE = re.compile(
    r"(?P<date>\d{4}-\d{2}-\d{2})[_-](?P<hh>\d{2})[-_](?P<mm>\d{2})(?:[-_](?P<ss>\d{2}))?"
)

def normalizeTimestamp(name: str) -> Optional[str]:
    match = TIMESTAMP_RE.search(name)
    if not match: return None

    date = match.group("date")
    hh = match.group("hh")
    mm = match.group("mm")
    ss = match.group("ss") or "00"
    return f"{date}_{hh}_{mm}_{ss}"

# Match measurement types based on keywords in the name
typesMap = {"longrun": ["longrun", "long_run"],
            "noise scan": ["noise"],
            "efficiency scan": ["layer"],
            "debug": ["debug", "test", "prova"]}
def parseMeasurementType(lowerName: str) -> Optional[str]:
    for type_name, keywords in typesMap.items():
        if any(keyword in lowerName for keyword in keywords):
            return type_name
    return None

# Match mixture types based on keywords like "std", "mix1", "mix2", etc.
def parseMixtureType(lowerName: str) -> Optional[str]:
    match = re.search(r"(?:^|[_-])(std|mix(\d+))(?:$|[_-])", lowerName)
    if not match:
        return None
    if match.group(1) == "std":
        return "std"
    return match.group(2)

# Match source status based on keywords like "sourceon", "source_off", etc.
sourceMap = {"ON": ["sourceon", "source_on", "source-on"],
             "OFF": ["sourceoff", "source_off", "source-off"]}
def parseSource(lowerName: str) -> Optional[str]:
    for status, keywords in sourceMap.items():
        if any(keyword in lowerName for keyword in keywords):
            return status
    return None

# Match filter settings based on patterns like "absX", "abs_X", "abs-X_Y", etc.
def parseFilterSetting(lowerName: str) -> Optional[float]:
    match = re.search(r"(?:^|[_-])abs[_-]?(\d+(?:[_-]\d+)*)(?:$|[_-])", lowerName)
    if not match:
        return None
    token = match.group(1)
    if "_" in token or "-" in token:
        parts = re.split(r"[_-]", token)
        if len(parts) >= 2:
            return float(f"{parts[0]}.{parts[1]}")
    return float(token)

# Match high voltage settings based on patterns like "HVXXX" and "lyXY_HVXXX"
def parseVoltages(name: str) -> Tuple[Optional[int], Optional[int]]:
    hv_matches = list(re.finditer(r"HV_?(\d+)", name))
    scanned = int(hv_matches[0].group(1)) if hv_matches else None
    other_match = re.search(r"ly\d+_HV_?(\d+)", name, re.IGNORECASE)
    other = int(other_match.group(1)) if other_match else None
    all_match = re.search(r"all[^H]*HV_?(\d+)", name, re.IGNORECASE)
    if all_match:
        hv_value = int(all_match.group(1))
        if scanned is None:
            scanned = hv_value
        other = scanned
    return scanned, other

# Match LV settings based on given patterns
lvPatternMap = {"LVset": r"(?:^|[_-])lvset(\d+)(?:$|[_-])",
                "LV": r"(?:^|[_-])lv(\d+)(?:$|[_-])",
                "LVsetting": r"(?:^|[_-])lvsetting(\d+)(?:$|[_-])"}
def parseLVSetting(lowerName: str) -> Optional[int]:
    match = re.search(r"(?:^|[_-])lvset(\d+)(?:$|[_-])", lowerName)
    if match: return int(match.group(1))
    match = re.search(r"(?:^|[_-])lv(\d+)(?:$|[_-])", lowerName)
    return int(match.group(1)) if match else None

# Index raw data directories
def indexRawData(raw_dir: Path) -> Dict[str, List[str]]:
    index: Dict[str, List[str]] = {}
    for path in raw_dir.rglob("*"):
        if not path.is_dir(): continue

        ts = normalizeTimestamp(path.name)
        if not ts: continue

        files = [str(p) for p in path.iterdir() if p.is_file()]
        if files: index[ts] = files
    return index

# Load existing JSON data if necessary
def loadExisting(output_path: Path) -> Dict:
    if not output_path.exists():
        return {}
    try:
        with output_path.open("r", encoding="utf-8") as handle:
            content = handle.read().strip()
            return json.loads(content) if content else {}
    except json.JSONDecodeError:
        return {}

# Build the info dictionary for a measurement file, extracting all relevant metadata.
def buildInfo(
    measure_path: Path,
    raw_paths: List[str],
    context_name: str,
) -> "OrderedDict[str, object]":

    # Pass lowercase name to some parsers for consistent matching
    name = measure_path.name
    lower_context_name = context_name.lower()

    scan_type = parseMeasurementType(lower_context_name)
    mixture = parseMixtureType(lower_context_name)
    lv_setting = parseLVSetting(lower_context_name)
    source = parseSource(lower_context_name)
    filter_setting = parseFilterSetting(lower_context_name) if source == "ON" else None
    if scan_type == "noise scan" and source == "ON": scan_type = "source scan"
    scanned_voltage, other_voltage = parseVoltages(name)

    info = OrderedDict()
    info["type"] = scan_type
    info["mixture"] = mixture
    info["lv_setting"] = lv_setting
    info["source_status"] = source
    info["filter"] = filter_setting
    info["scanned_hv"] = scanned_voltage
    if scan_type != "noise scan": info["other_hv"] = other_voltage
    info["path_plot"] = str(measure_path)
    info["path_raw"] = raw_paths
    return info


def warnIfConflict(
    field_name: str,
    file_value: Optional[object],
    dir_value: Optional[object],
    file_name: str,
    dir_name: str,
) -> None:
    if file_value is None or dir_value is None:
        return
    if file_value != dir_value:
        print(
            f"Warning: {field_name} conflict between file '{file_name}' and dir '{dir_name}'",
            file=sys.stderr,
        )

# Arguments parserer
def buildParser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Match measurement files to raw data by timestamp.")
    parser.add_argument(
        "--measure-dir",
        "-m",
        required=True, help="Directory with measurement files.")
    parser.add_argument(
        "--raw-dir",
        "-r",
        required=True, help="Directory with raw data (timestamp dirs).")
    parser.add_argument(
        "--output",
        "-o",
        default=str(Path(__file__).resolve().parents[1] / "rawData.json"),
        help="Output JSON path (default: rpc-daq/viktorg/rawData.json).",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing JSON instead of merging.",
    )
    return parser


def main() -> int:

    args = buildParser().parse_args()

    measure_dir = Path(args.measure_dir).expanduser().resolve()
    raw_dir = Path(args.raw_dir).expanduser().resolve()
    output_path = Path(args.output).expanduser().resolve()

    raw_index = indexRawData(raw_dir)
    if args.overwrite:
        merged = {}
    else:
        existing = loadExisting(output_path)
        merged = dict(existing)

    for measure_path in sorted(measure_dir.rglob("*")):
        if not measure_path.is_file(): continue

        relative_parts = measure_path.relative_to(measure_dir).parts
        relative_name = "_".join(relative_parts)
        dir_name = "_".join(relative_parts[:-1])
        measurement_name = measure_path.stem

        ts = normalizeTimestamp(measure_path.name)
        if not ts: continue

        file_lower = measure_path.name.lower()
        dir_lower = dir_name.lower() if dir_name else ""
        warnIfConflict("type", parseMeasurementType(file_lower), parseMeasurementType(dir_lower), measure_path.name, dir_name)
        warnIfConflict("mixture", parseMixtureType(file_lower), parseMixtureType(dir_lower), measure_path.name, dir_name)
        warnIfConflict("lv_setting", parseLVSetting(file_lower), parseLVSetting(dir_lower), measure_path.name, dir_name)
        warnIfConflict("source_status", parseSource(file_lower), parseSource(dir_lower), measure_path.name, dir_name)
        warnIfConflict("filter", parseFilterSetting(file_lower), parseFilterSetting(dir_lower), measure_path.name, dir_name)

        raw_paths = raw_index.get(ts)
        if not raw_paths: continue

        if ts not in merged: merged[ts] = OrderedDict()

        info = buildInfo(measure_path, sorted(raw_paths), relative_name)
        merged[ts][measurement_name] = info

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        sorted_merged = OrderedDict(sorted(merged.items(), key=lambda item: item[0]))
        json.dump(sorted_merged, handle, indent=4)
        handle.write("\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())