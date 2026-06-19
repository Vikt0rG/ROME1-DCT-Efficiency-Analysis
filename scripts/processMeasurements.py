import os
import sys
import subprocess
import shutil

import tarfile
import yaml

import argparse


# Safely extract tar files to prevent path traversal vulnerabilities
def safe_extract(tar, path):
    for member in tar.getmembers():
        member_path = os.path.join(path, member.name)
        abs_path = os.path.abspath(path)
        abs_member = os.path.abspath(member_path)
        if not abs_member.startswith(abs_path + os.sep):
            raise RuntimeError(f"Unsafe path in tar: {member.name}")
    tar.extractall(path)


# Find the bin_data_files.tar.gz in the list of raw files
def find_bin_tar(raw_files):
    for item in raw_files:
        if "bin_data_files.tar.gz" in item:
            return item
    return None


# Check if there is a compressed tarball containing text files
def find_text_tar(raw_files):
    for item in raw_files:
        if "text_data_files.tar.gz" in item:
            return item
    return None


# Check if there is directly an uncompressed txt file
def find_raw_txt(raw_files, name):
    txt_name = f"{name}.txt"
    for file in raw_files:
        if os.path.basename(file) == txt_name:
            return file
    for file in raw_files:
        if file.endswith(".txt"):
            return file
    return None


# Check if a path exists, or if it can be found in the fallback directory
def resolve_path(path, fallback_dir):
    if not path:
        return None
    if os.path.exists(path):
        return path
    candidate = os.path.join(fallback_dir, os.path.basename(path))
    return candidate if os.path.exists(candidate) else None


# Collect all entries from the config's data section
def collect_entries(config):
    data = config.get("data", {})
    if isinstance(data, list):
        return data
    if not isinstance(data, dict):
        return []
    entries = []
    for _, runs in data.items():
        if isinstance(runs, list):
            entries.extend(runs)
    return entries


def normalize_raw_files(raw_files):
    if isinstance(raw_files, str):
        if os.path.isdir(raw_files):
            return [os.path.join(raw_files, f) for f in os.listdir(raw_files)]
        return [raw_files]
    if isinstance(raw_files, list):
        expanded = []
        for item in raw_files:
            if isinstance(item, str) and os.path.isdir(item):
                expanded.extend(
                    os.path.join(item, f) for f in os.listdir(item)
                )
            else:
                expanded.append(item)
        return expanded
    return []


# Parse command-line arguments
def parse_args():
    parser = argparse.ArgumentParser(
        description="Process measurement configs into ROOT files."
    )
    parser.add_argument("configs", nargs="+", help="Config YAML file(s)")
    parser.add_argument("--dt-max", type=int, default=-100,
                        help="Max time window (default: -100)")
    parser.add_argument("--dt-min", type=int, default=-180,
                        help="Min time window (default: -180)")
    parser.add_argument("--no-external", action="store_true",
                        help="Disable external trigger usage in analysis")
    parser.add_argument("--force", action="store_true",
                        help="Rebuild txt/root even if they exist")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Enable verbose logging")
    return parser.parse_args()


# Main processing function
def main():
    args = parse_args()

    def vprint(message):
        if args.verbose:
            print(message)

    # Step 1: Resolve paths and tools
    root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    target_bin_dir = os.path.join(root_dir, "data", "raw", "bin")
    target_txt_dir = os.path.join(root_dir, "data", "raw", "txt")
    target_root_dir = os.path.join(root_dir, "data", "root")
    bin_to_txt = os.path.join(root_dir, "scripts", "bin_to_txt.sh")
    analysis_bin = os.path.join(root_dir, "bin", "analysis")

    # Step 2: Ensure output directories exist
    os.makedirs(target_bin_dir, exist_ok=True)
    os.makedirs(target_txt_dir, exist_ok=True)
    os.makedirs(target_root_dir, exist_ok=True)

    # Recompile the binary once before processing configs
    vprint("Recompiling analysis binary...")
    subprocess.run(["make", "build"], check=True, cwd=root_dir)

    # Step 3: Process each config file
    for cfg_path in args.configs:
        if not os.path.isfile(cfg_path):
            print(f"Config not found: {cfg_path}", file=sys.stderr)
            continue

        vprint(f"\nConfig: {cfg_path}")

        with open(cfg_path, "r", encoding="utf-8") as handle:
            config = yaml.safe_load(handle) or {}

        # Step 4: Process each run entry
        config_changed = False
        for entry in collect_entries(config):
            name = entry.get("name")
            raw_files = normalize_raw_files(
                entry.get("raw files", entry.get("path_raw", []))
            )
            if not name or not raw_files:
                continue

            # If a finalized txt already exists, skip re-processing unless
            # forced
            expected_txt = os.path.join(target_txt_dir, f"{name}.txt")
            if os.path.exists(expected_txt) and not args.force:
                vprint(
                    f"Skipping raw processing; txt already exists: {expected_txt}"
                )
                txt_paths = [expected_txt]
            else:
                # Step 5: Determine the txt source
                # (raw txt -> text tar -> bin tar)
                raw_txt_path = find_raw_txt(raw_files, name)
                text_tar_path = find_text_tar(raw_files)
                bin_tar_path = find_bin_tar(raw_files)

                txt_paths = []

                # Check if the source raw txt file exists or can be found in
                # the fallback directory
                if raw_txt_path:
                    source_txt = resolve_path(raw_txt_path, target_txt_dir)
                    if not source_txt:
                        print(f"Missing raw txt for {name}: {raw_txt_path}")
                        continue

                    # Step 6a: Use raw txt directly
                    target_txt = os.path.join(target_txt_dir, f"{name}.txt")
                    if os.path.abspath(source_txt) != os.path.abspath(target_txt):

                        # If the target already exists and --force is not set,
                        # skip moving
                        if os.path.exists(target_txt) and not args.force:
                            pass
                        else:
                            # Remove existing target for --force if it exists
                            # to avoid issues with shutil.copy2
                            if os.path.exists(target_txt):
                                os.remove(target_txt)
                            shutil.copy2(source_txt, target_txt)
                            vprint(f"Copied raw txt to {target_txt}")
                    else:
                        vprint(f"Using existing raw txt {target_txt}")

                    txt_paths = [target_txt]

                    # Cleanup any tar files that might have been used in this
                    # process to avoid confusion later
                    for cleanup_path in (text_tar_path, bin_tar_path):
                        resolved = resolve_path(cleanup_path, target_bin_dir)
                        if resolved and os.path.abspath(resolved).startswith(root_dir):
                            os.remove(resolved)

                # If no raw txt file, check for a text tarball and extract it
                elif text_tar_path:
                    # Step 6b: Extract text tar into txt dir
                    source_tar = resolve_path(text_tar_path, target_bin_dir)
                    if not source_tar:
                        print(f"Missing text tar for {name}: {text_tar_path}")
                        continue

                    with tarfile.open(source_tar, "r:gz") as tar:
                        safe_extract(tar, target_txt_dir)
                        members = [m.name for m in tar.getmembers() if m.name.endswith(".txt")]
                    vprint(f"Extracted text tar {source_tar}")
                    filedump_txt = os.path.join(target_txt_dir, "filedump_0.txt")
                    if os.path.exists(filedump_txt):
                        target_txt = os.path.join(target_txt_dir, f"{name}.txt")
                        if os.path.exists(target_txt) and not args.force:
                            pass
                        else:
                            if os.path.exists(target_txt):
                                os.remove(target_txt)
                            shutil.move(filedump_txt, target_txt)
                            vprint(f"Renamed filedump_0.txt to {target_txt}")
                        txt_paths = [target_txt]
                    else:
                        txt_paths = [os.path.join(target_txt_dir, os.path.basename(m)) for m in members]
                        if len(txt_paths) == 1:
                            target_txt = os.path.join(target_txt_dir, f"{name}.txt")
                            if os.path.basename(txt_paths[0]) != f"{name}.txt":
                                if os.path.exists(target_txt) and not args.force:
                                    pass
                                else:
                                    if os.path.exists(target_txt):
                                        os.remove(target_txt)
                                    shutil.move(txt_paths[0], target_txt)
                                    vprint(f"Renamed {txt_paths[0]} to {target_txt}")
                                txt_paths = [target_txt]

                    if os.path.exists(source_tar):
                        os.remove(source_tar)

                else:
                    # Step 6c: Extract bin tar and convert pcap to txt
                    if not bin_tar_path:
                        print(f"No bin_data_files.tar.gz for {name}")
                        continue

                    source_tar = resolve_path(bin_tar_path, target_bin_dir)
                    if not source_tar:
                        print(f"Missing binary tar for {name}: {bin_tar_path}")
                        continue

                    target_tar = os.path.join(target_bin_dir, f"{name}.tar.gz")
                    if os.path.abspath(source_tar) != os.path.abspath(target_tar):
                        if os.path.exists(target_tar) and not args.force:
                            pass
                        else:
                            shutil.copy2(source_tar, target_tar)

                    if not os.path.exists(target_tar):
                        print(f"Missing target tar for {name}: {target_tar}")
                        continue

                    with tarfile.open(target_tar, "r:gz") as tar:
                        safe_extract(tar, target_bin_dir)
                        members = [m.name for m in tar.getmembers() if m.name.endswith(".pcapng")]
                    vprint(f"Extracted bin tar {target_tar}")

                    pcap_paths = [os.path.join(target_bin_dir, os.path.basename(m)) for m in members]
                    if not pcap_paths:
                        print(f"No .pcapng files in {target_tar}")
                        continue

                    for pcap in pcap_paths:
                        txt_path = os.path.join(
                            target_txt_dir,
                            os.path.basename(pcap).replace(".pcapng", ".txt")
                        )
                        txt_paths.append(txt_path)
                        if os.path.exists(txt_path) and not args.force:
                            continue

                        vprint(f"Converting {pcap} to txt")
                        subprocess.run(
                            ["bash", bin_to_txt, "--input", pcap, "--force" if args.force else ""],
                            check=True,
                        )

                    if len(txt_paths) == 1 and os.path.exists(txt_paths[0]):
                        if os.path.basename(txt_paths[0]) != f"{name}.txt":
                            target_txt = os.path.join(target_txt_dir, f"{name}.txt")
                            if os.path.exists(target_txt) and not args.force:
                                pass
                            else:
                                if os.path.exists(target_txt):
                                    os.remove(target_txt)
                                shutil.move(txt_paths[0], target_txt)
                                vprint(f"Renamed {txt_paths[0]} to {target_txt}")
                            txt_paths = [target_txt]

                    for pcap in pcap_paths:
                        if os.path.exists(pcap) and os.path.abspath(pcap).startswith(root_dir):
                            os.remove(pcap)
                    if os.path.exists(target_tar) and os.path.abspath(target_tar).startswith(root_dir):
                        os.remove(target_tar)

                if os.path.exists(expected_txt) and expected_txt not in txt_paths:
                    txt_paths.append(expected_txt)

            # Step 7: Run analysis and write ROOT output
            root_path = os.path.join(target_root_dir, f"{name}.root")
            if os.path.exists(root_path) and not args.force:
                continue

            existing_txt = [p for p in txt_paths if os.path.exists(p)]
            if len(existing_txt) != 1:
                raise ValueError(
                    f"Expected 1 txt for {name}, found {len(existing_txt)}.",
                    "Skipping analysis."
                )

                continue

            command = "process"
            analysis_cmd = [analysis_bin, command,
                            existing_txt[0],
                            str(args.dt_max), str(args.dt_min)]
            if args.no_external:
                analysis_cmd.append("--no-external")

            subprocess.run(analysis_cmd, check=True, cwd=root_dir)

            output_root = os.path.join(root_dir, "output.root")
            if not os.path.exists(output_root):
                print(f"Missing output.root for {name}")
                continue

            shutil.move(output_root, root_path)
            print(f"Wrote {root_path}")

            # Step 8: Update config with root file path if not already set
            if entry.get("root file") != root_path:
                entry["root file"] = root_path
                config_changed = True

        if config_changed:
            with open(cfg_path, "w", encoding="utf-8") as handle:
                yaml.safe_dump(config, handle, sort_keys=False)
            vprint(f"Updated config file: {cfg_path}")


if __name__ == "__main__":
    main()
